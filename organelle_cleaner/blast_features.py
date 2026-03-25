#!/usr/bin/env python3
"""Parse BLAST tabular outputs into contig-level organelle support summaries."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable


BLAST_SUPPORT_STRONG = "strong_blast_support"
BLAST_SUPPORT_MODERATE = "moderate_blast_support"
BLAST_SUPPORT_NONE = "no_blast_support"
BLAST_SOURCE_PLASTID = "plastid"
BLAST_SOURCE_MITOCHONDRIAL = "mitochondrial"


@dataclass(frozen=True)
class BlastThresholdConfig:
    strong_identity: float = 98.0
    strong_aligned_fraction: float = 0.8
    # The minimum default BLAST support cutoff is 90% identity over 80% of the
    # subject contig. Strong and exceptional tiers remain stricter overlays.
    moderate_identity: float = 90.0
    moderate_aligned_fraction: float = 0.8
    exceptional_identity: float = 99.0
    exceptional_aligned_fraction: float = 0.95


@dataclass(frozen=True)
class BlastHit:
    query_id: str
    contig_id: str
    identity: float
    alignment_length: int
    subject_length: int
    subject_start: int
    subject_end: int
    source: str

    @property
    def aligned_fraction(self) -> float:
        return self.alignment_length / self.subject_length


def _to_finite_float(value: str) -> float | None:
    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return None

    if not math.isfinite(numeric_value):
        return None
    return numeric_value


def _to_positive_int(value: str) -> int | None:
    try:
        numeric_value = int(value)
    except (TypeError, ValueError):
        return None

    if numeric_value <= 0:
        return None
    return numeric_value


def _normalize_interval(start: int, end: int) -> tuple[int, int]:
    if start <= end:
        return start, end
    return end, start


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    normalized = sorted(_normalize_interval(start, end) for start, end in intervals)
    if not normalized:
        return []

    merged: list[tuple[int, int]] = [normalized[0]]
    for start, end in normalized[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def merged_interval_length(intervals: Iterable[tuple[int, int]]) -> int:
    return sum((end - start) + 1 for start, end in merge_intervals(intervals))


def _support_rank(level: str) -> int:
    if level == BLAST_SUPPORT_STRONG:
        return 2
    if level == BLAST_SUPPORT_MODERATE:
        return 1
    return 0


def _source_summary(
    source_hits: list[BlastHit],
    *,
    contig_length: int,
    config: BlastThresholdConfig,
) -> dict[str, object]:
    if not source_hits:
        return {
            "hit_count": 0,
            "best_identity": None,
            "best_hsp_aligned_fraction": None,
            "merged_aligned_bp": 0,
            "merged_coverage_fraction": 0.0,
            "strong_merged_aligned_bp": 0,
            "strong_merged_coverage_fraction": 0.0,
            "moderate_merged_aligned_bp": 0,
            "moderate_merged_coverage_fraction": 0.0,
            "support_level": BLAST_SUPPORT_NONE,
            "has_exceptional_support": False,
            "support_description": "",
        }

    best_identity = max(hit.identity for hit in source_hits)
    best_hsp_aligned_fraction = max(hit.aligned_fraction for hit in source_hits)
    merged_aligned_bp = merged_interval_length((hit.subject_start, hit.subject_end) for hit in source_hits)
    coverage_fraction = merged_aligned_bp / contig_length if contig_length > 0 else 0.0

    strong_hits = [hit for hit in source_hits if hit.identity >= config.strong_identity]
    moderate_hits = [hit for hit in source_hits if hit.identity >= config.moderate_identity]
    exceptional_hits = [hit for hit in source_hits if hit.identity >= config.exceptional_identity]

    strong_merged_aligned_bp = merged_interval_length(
        (hit.subject_start, hit.subject_end) for hit in strong_hits
    )
    moderate_merged_aligned_bp = merged_interval_length(
        (hit.subject_start, hit.subject_end) for hit in moderate_hits
    )
    exceptional_merged_aligned_bp = merged_interval_length(
        (hit.subject_start, hit.subject_end) for hit in exceptional_hits
    )

    strong_coverage_fraction = (
        strong_merged_aligned_bp / contig_length if contig_length > 0 else 0.0
    )
    moderate_coverage_fraction = (
        moderate_merged_aligned_bp / contig_length if contig_length > 0 else 0.0
    )
    exceptional_coverage_fraction = (
        exceptional_merged_aligned_bp / contig_length if contig_length > 0 else 0.0
    )

    if strong_coverage_fraction >= config.strong_aligned_fraction:
        support_level = BLAST_SUPPORT_STRONG
        support_coverage_fraction = strong_coverage_fraction
        support_identity_floor = config.strong_identity
    elif moderate_coverage_fraction >= config.moderate_aligned_fraction:
        support_level = BLAST_SUPPORT_MODERATE
        support_coverage_fraction = moderate_coverage_fraction
        support_identity_floor = config.moderate_identity
    else:
        support_level = BLAST_SUPPORT_NONE
        support_coverage_fraction = coverage_fraction
        support_identity_floor = None

    has_exceptional_support = (
        exceptional_coverage_fraction >= config.exceptional_aligned_fraction
    )
    support_description = ""
    if support_level != BLAST_SUPPORT_NONE:
        support_description = (
            f"{source_hits[0].source} BLAST merged {merged_aligned_bp} bp "
            f"overall and {support_coverage_fraction:.3f} contig coverage from hits "
            f"at >= {support_identity_floor:.2f}% identity; best HSP identity {best_identity:.2f}%"
        )

    return {
        "hit_count": len(source_hits),
        "best_identity": round(best_identity, 4),
        "best_hsp_aligned_fraction": round(best_hsp_aligned_fraction, 4),
        "merged_aligned_bp": merged_aligned_bp,
        "merged_coverage_fraction": round(coverage_fraction, 4),
        "strong_merged_aligned_bp": strong_merged_aligned_bp,
        "strong_merged_coverage_fraction": round(strong_coverage_fraction, 4),
        "moderate_merged_aligned_bp": moderate_merged_aligned_bp,
        "moderate_merged_coverage_fraction": round(moderate_coverage_fraction, 4),
        "support_level": support_level,
        "has_exceptional_support": has_exceptional_support,
        "support_description": support_description,
    }


def parse_blast_hits(
    path: Path,
    *,
    source: str,
) -> list[BlastHit]:
    """Parse BLAST tabular output with organelle references as query and assembly contigs as subject."""
    hits: list[BlastHit] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, start=1):
            if not fields:
                continue
            if fields[0].startswith("#"):
                continue
            if fields[0] == "qseqid":
                continue
            if len(fields) < 11:
                raise ValueError(
                    f"{path}: line {line_number} has {len(fields)} columns; expected at least 11"
                )

            identity = _to_finite_float(fields[2])
            alignment_length = _to_positive_int(fields[3])
            subject_length = _to_positive_int(fields[4])
            subject_start = _to_positive_int(fields[9])
            subject_end = _to_positive_int(fields[10])
            if identity is None:
                raise ValueError(f"{path}: line {line_number} has invalid pident value {fields[2]!r}")
            if alignment_length is None:
                raise ValueError(f"{path}: line {line_number} has invalid alignment length {fields[3]!r}")
            if subject_length is None:
                raise ValueError(f"{path}: line {line_number} has invalid subject length {fields[4]!r}")
            if subject_start is None:
                raise ValueError(f"{path}: line {line_number} has invalid subject start {fields[9]!r}")
            if subject_end is None:
                raise ValueError(f"{path}: line {line_number} has invalid subject end {fields[10]!r}")

            hits.append(
                BlastHit(
                    query_id=fields[0],
                    # Assembly support is tracked on the subject contig ID and subject span.
                    contig_id=fields[1],
                    identity=identity,
                    alignment_length=alignment_length,
                    subject_length=subject_length,
                    subject_start=subject_start,
                    subject_end=subject_end,
                    source=source,
                )
            )
    return hits


def aggregate_blast_hits(
    hits: Iterable[BlastHit],
    *,
    config: BlastThresholdConfig | None = None,
) -> dict[str, dict[str, object]]:
    if config is None:
        config = BlastThresholdConfig()

    hits_by_contig: dict[str, list[BlastHit]] = {}
    for hit in hits:
        hits_by_contig.setdefault(hit.contig_id, []).append(hit)

    rows: dict[str, dict[str, object]] = {}
    for contig_id, contig_hits in hits_by_contig.items():
        contig_length = max(hit.subject_length for hit in contig_hits)
        source_hits = {
            BLAST_SOURCE_PLASTID: [
                hit for hit in contig_hits if hit.source == BLAST_SOURCE_PLASTID
            ],
            BLAST_SOURCE_MITOCHONDRIAL: [
                hit for hit in contig_hits if hit.source == BLAST_SOURCE_MITOCHONDRIAL
            ],
        }
        source_summaries = {
            source: _source_summary(source_hits_for_source, contig_length=contig_length, config=config)
            for source, source_hits_for_source in source_hits.items()
        }

        strongest_source = max(
            source_summaries,
            key=lambda source: (
                _support_rank(str(source_summaries[source]["support_level"])),
                float(source_summaries[source]["merged_coverage_fraction"]),
                float(source_summaries[source]["best_identity"] or 0.0),
            ),
        )
        strongest_summary = source_summaries[strongest_source]
        strongest_support = str(strongest_summary["support_level"])
        has_exceptional_support = any(
            bool(summary["has_exceptional_support"]) for summary in source_summaries.values()
        )
        support_sources: set[str] = {
            source
            for source, summary in source_summaries.items()
            if str(summary["support_level"]) != BLAST_SUPPORT_NONE
        }
        support_descriptions: list[str] = []
        for source in (BLAST_SOURCE_PLASTID, BLAST_SOURCE_MITOCHONDRIAL):
            summary = source_summaries[source]
            support_description = str(summary["support_description"])
            if support_description:
                support_descriptions.append(support_description)

        combined_merged_aligned_bp = merged_interval_length(
            (hit.subject_start, hit.subject_end) for hit in contig_hits
        )
        combined_coverage_fraction = (
            combined_merged_aligned_bp / contig_length if contig_length > 0 else 0.0
        )
        best_identity = max(hit.identity for hit in contig_hits)
        best_aligned_fraction = max(hit.aligned_fraction for hit in contig_hits)

        if strongest_support == BLAST_SUPPORT_STRONG:
            blast_score = 4
        elif strongest_support == BLAST_SUPPORT_MODERATE:
            blast_score = 2
        else:
            blast_score = 0

        rows[contig_id] = {
            "contig_id": contig_id,
            "blast_contig_length": contig_length,
            "blast_support_level": strongest_support,
            "blast_best_identity": round(best_identity, 4),
            "blast_best_aligned_fraction": round(best_aligned_fraction, 4),
            "blast_merged_aligned_bp": combined_merged_aligned_bp,
            "blast_merged_coverage_fraction": round(combined_coverage_fraction, 4),
            "blast_score": blast_score,
            "blast_hit_count": len(contig_hits),
            "blast_support_sources": ",".join(sorted(support_sources)),
            "blast_selected_source": strongest_source,
            "blast_supporting_evidence": "; ".join(support_descriptions),
            "has_strong_blast_support": strongest_support == BLAST_SUPPORT_STRONG,
            "has_moderate_blast_support": strongest_support in {
                BLAST_SUPPORT_STRONG,
                BLAST_SUPPORT_MODERATE,
            },
            "has_exceptional_blast_support": has_exceptional_support,
            "plastid_support_level": source_summaries[BLAST_SOURCE_PLASTID]["support_level"],
            "plastid_best_identity": source_summaries[BLAST_SOURCE_PLASTID]["best_identity"],
            "plastid_hit_count": source_summaries[BLAST_SOURCE_PLASTID]["hit_count"],
            "plastid_merged_aligned_bp": source_summaries[BLAST_SOURCE_PLASTID]["merged_aligned_bp"],
            "plastid_merged_coverage_fraction": source_summaries[BLAST_SOURCE_PLASTID]["merged_coverage_fraction"],
            "mit_support_level": source_summaries[BLAST_SOURCE_MITOCHONDRIAL]["support_level"],
            "mit_best_identity": source_summaries[BLAST_SOURCE_MITOCHONDRIAL]["best_identity"],
            "mit_hit_count": source_summaries[BLAST_SOURCE_MITOCHONDRIAL]["hit_count"],
            "mit_merged_aligned_bp": source_summaries[BLAST_SOURCE_MITOCHONDRIAL]["merged_aligned_bp"],
            "mit_merged_coverage_fraction": source_summaries[BLAST_SOURCE_MITOCHONDRIAL]["merged_coverage_fraction"],
        }

    return rows


def load_blast_support_by_contig(
    *,
    plastid_blast_tsv: Path | None = None,
    chl_blast_tsv: Path | None = None,
    mit_blast_tsv: Path | None = None,
    config: BlastThresholdConfig | None = None,
) -> dict[str, dict[str, object]]:
    hits: list[BlastHit] = []
    plastid_paths: list[Path] = []
    if plastid_blast_tsv is not None:
        plastid_paths.append(plastid_blast_tsv)
    if chl_blast_tsv is not None:
        plastid_paths.append(chl_blast_tsv)
    seen_plastid_paths: set[Path] = set()
    for plastid_path in plastid_paths:
        normalized_path = plastid_path.resolve()
        if normalized_path in seen_plastid_paths:
            continue
        seen_plastid_paths.add(normalized_path)
        hits.extend(parse_blast_hits(plastid_path, source=BLAST_SOURCE_PLASTID))
    if mit_blast_tsv is not None:
        hits.extend(parse_blast_hits(mit_blast_tsv, source=BLAST_SOURCE_MITOCHONDRIAL))
    return aggregate_blast_hits(hits, config=config)
