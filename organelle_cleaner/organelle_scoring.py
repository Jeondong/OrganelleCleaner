#!/usr/bin/env python3
"""Score contigs for graph-only, BLAST-only, and hybrid organelle detection."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import math
import statistics
import sys
from pathlib import Path
from typing import Iterable

import networkx as nx

from .blast_features import (
    BLAST_SUPPORT_MODERATE,
    BLAST_SUPPORT_NONE,
    BLAST_SUPPORT_STRONG,
    BlastThresholdConfig,
    load_blast_support_by_contig,
)
from .graph_analysis import summarize_contigs
from .parse_gfa import parse_gfa
from .sequence_features import calculate_contig_features, estimate_nuclear_coverage


ScoreRow = dict[str, object]

MODE_GRAPH_ONLY = "graph-only"
MODE_BLAST_ONLY = "blast-only"
MODE_HYBRID = "hybrid"
VALID_MODES = (MODE_GRAPH_ONLY, MODE_BLAST_ONLY, MODE_HYBRID)

HIGH_CONFIDENCE = "high_confidence"
MEDIUM_CONFIDENCE = "medium_confidence"
LOW_CONFIDENCE = "low_confidence"
NOT_FLAGGED = "not_flagged"

SIZE_BIN_VERY_SMALL = "very_small"
SIZE_BIN_SMALL_FRAGMENT = "small_fragment"
SIZE_BIN_PREFERRED_RANGE = "preferred_range"
SIZE_BIN_INTERMEDIATE = "intermediate"
SIZE_BIN_LARGE = "large"
SIZE_BIN_VERY_LARGE = "very_large"

DEFAULT_HIGH_CONFIDENCE_THRESHOLD = 5
DEFAULT_MEDIUM_CONFIDENCE_THRESHOLD = 3
DEFAULT_LOW_CONFIDENCE_THRESHOLD = 1
DEFAULT_MIN_ORGANELLE_LENGTH = 50_000
DEFAULT_PREFERRED_MAX_LENGTH = 150_000
DEFAULT_LARGE_CONTIG_LENGTH = 300_000
DEFAULT_VERY_LARGE_CONTIG_LENGTH = 500_000
DEFAULT_ISOLATED_LENGTH_CAP = 150_000
DEFAULT_LARGE_CONTIG_PENALTY = 2
DEFAULT_VERY_LARGE_CONTIG_PENALTY = 2
DEFAULT_LARGE_CONTIG_HIGH_MIN_SIGNALS = 3
DEFAULT_LARGE_CONTIG_MEDIUM_MIN_SIGNALS = 2
DEFAULT_COMPACT_COMPONENT_SIZE = 6
DEFAULT_GC_DEVIATION_THRESHOLD = 5.0
DEFAULT_COVERAGE_MULTIPLIER_THRESHOLD = 5.0
DEFAULT_VERY_SMALL_MAX = 10_000
DEFAULT_SMALL_FRAGMENT_MAX = 49_999
DEFAULT_INTERMEDIATE_MAX = 300_000
GC_TRIM_FRACTION = 0.10


@dataclass(frozen=True)
class ScoringConfig:
    high_confidence_threshold: int = DEFAULT_HIGH_CONFIDENCE_THRESHOLD
    medium_confidence_threshold: int = DEFAULT_MEDIUM_CONFIDENCE_THRESHOLD
    low_confidence_threshold: int = DEFAULT_LOW_CONFIDENCE_THRESHOLD
    min_organelle_length: int = DEFAULT_MIN_ORGANELLE_LENGTH
    preferred_max_length: int = DEFAULT_PREFERRED_MAX_LENGTH
    large_contig_length: int = DEFAULT_LARGE_CONTIG_LENGTH
    very_large_contig_length: int = DEFAULT_VERY_LARGE_CONTIG_LENGTH
    isolated_length_cap: int = DEFAULT_ISOLATED_LENGTH_CAP
    large_contig_penalty: int = DEFAULT_LARGE_CONTIG_PENALTY
    very_large_contig_penalty: int = DEFAULT_VERY_LARGE_CONTIG_PENALTY
    large_contig_high_min_signals: int = DEFAULT_LARGE_CONTIG_HIGH_MIN_SIGNALS
    large_contig_medium_min_signals: int = DEFAULT_LARGE_CONTIG_MEDIUM_MIN_SIGNALS
    compact_component_size: int = DEFAULT_COMPACT_COMPONENT_SIZE
    gc_deviation_threshold: float = DEFAULT_GC_DEVIATION_THRESHOLD
    coverage_multiplier_threshold: float = DEFAULT_COVERAGE_MULTIPLIER_THRESHOLD
    very_small_max_length: int = DEFAULT_VERY_SMALL_MAX
    small_fragment_max_length: int = DEFAULT_SMALL_FRAGMENT_MAX
    intermediate_max_length: int = DEFAULT_INTERMEDIATE_MAX


def _to_finite_float(value: object) -> float | None:
    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric_value):
        return None
    return numeric_value


def _to_nonnegative_float(value: object) -> float | None:
    numeric_value = _to_finite_float(value)
    if numeric_value is None or numeric_value < 0:
        return None
    return numeric_value


def _to_nonnegative_int(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        numeric_value = int(value)
    except (TypeError, ValueError):
        return None
    if numeric_value < 0:
        return None
    return numeric_value


def _normalize_nuclear_coverage_baseline(value: object) -> float | None:
    numeric_value = _to_nonnegative_float(value)
    if numeric_value is None or numeric_value <= 0:
        return None
    return numeric_value


def estimate_gc_baseline(values: Iterable[object]) -> float | None:
    gc_values: list[float] = []
    for value in values:
        gc_value = _to_finite_float(value)
        if gc_value is None or gc_value < 0 or gc_value > 100:
            continue
        gc_values.append(gc_value)

    if not gc_values:
        return None

    gc_values.sort()
    trim_count = math.floor(len(gc_values) * GC_TRIM_FRACTION)
    if trim_count > 0 and len(gc_values) > (trim_count * 2):
        gc_values = gc_values[trim_count:-trim_count]

    if not gc_values:
        return None
    return float(statistics.median(gc_values))


def _tier_rank(confidence_tier: str) -> int:
    if confidence_tier == HIGH_CONFIDENCE:
        return 3
    if confidence_tier == MEDIUM_CONFIDENCE:
        return 2
    if confidence_tier == LOW_CONFIDENCE:
        return 1
    return 0


def _cap_tier(confidence_tier: str, maximum: str) -> str:
    if _tier_rank(confidence_tier) > _tier_rank(maximum):
        return maximum
    return confidence_tier


def _graph_positive_tier(confidence_tier: str) -> str | None:
    if confidence_tier in {HIGH_CONFIDENCE, MEDIUM_CONFIDENCE, LOW_CONFIDENCE}:
        return confidence_tier
    return None


def _format_value(value: float | int) -> str:
    if isinstance(value, int) or float(value).is_integer():
        return str(int(value))
    return f"{value:.2f}"


def classify_contig(confidence_tier: str) -> str:
    if confidence_tier == HIGH_CONFIDENCE:
        return "organelle"
    return "nuclear"


def determine_size_bin(length_value: int | None, config: ScoringConfig) -> str:
    if length_value is None:
        return SIZE_BIN_SMALL_FRAGMENT
    if length_value < config.very_small_max_length:
        return SIZE_BIN_VERY_SMALL
    if length_value <= config.small_fragment_max_length:
        return SIZE_BIN_SMALL_FRAGMENT
    if length_value <= config.preferred_max_length:
        return SIZE_BIN_PREFERRED_RANGE
    if length_value <= config.intermediate_max_length:
        return SIZE_BIN_INTERMEDIATE
    if length_value < config.large_contig_length:
        return SIZE_BIN_INTERMEDIATE
    if length_value <= config.very_large_contig_length:
        return SIZE_BIN_LARGE
    return SIZE_BIN_VERY_LARGE


def _is_compact_component(
    *,
    component_size: int | None,
    component_edge_count: int | None,
    component_branching_nodes: int | None,
    component_total_length: int | None,
    config: ScoringConfig,
) -> bool:
    if component_size is None or component_edge_count is None or component_branching_nodes is None:
        return False
    if component_size > config.compact_component_size:
        return False
    if component_branching_nodes != 0:
        return False
    if component_edge_count > component_size + 1:
        return False
    if component_total_length is not None and component_total_length >= config.large_contig_length:
        return False
    return True


def _has_preferred_range_correlated_pattern(
    *,
    size_bin: str,
    is_compact_component: bool,
    is_isolated: bool,
    has_gc_signal: bool,
) -> bool:
    return (
        size_bin == SIZE_BIN_PREFERRED_RANGE
        and is_compact_component
        and is_isolated
        and has_gc_signal
    )


def _build_decision_explanation(
    *,
    mode: str,
    confidence_tier: str,
    support_reasons: list[str],
    penalty_reasons: list[str],
    decision_notes: list[str],
) -> str:
    if confidence_tier == HIGH_CONFIDENCE:
        prefix = f"High-confidence organelle candidate in {mode} mode."
    elif confidence_tier in {MEDIUM_CONFIDENCE, LOW_CONFIDENCE}:
        prefix = f"Flagged for review in {mode} mode."
    else:
        prefix = f"Not flagged as an organelle candidate in {mode} mode."

    parts = [prefix]
    if support_reasons:
        parts.append(f"Evidence: {'; '.join(support_reasons)}.")
    if penalty_reasons:
        parts.append(f"Penalties: {'; '.join(penalty_reasons)}.")
    if decision_notes:
        parts.append(f"Decision: {'; '.join(decision_notes)}.")
    return " ".join(parts)


def evaluate_graph_features(
    *,
    feature_row: dict[str, object],
    component_row: dict[str, object],
    nuclear_coverage: float | None,
    gc_baseline: float | None,
    config: ScoringConfig,
) -> dict[str, object]:
    score = 0
    total_signal_count = 0
    strong_signal_count = 0
    topology_signal_count = 0
    topology_strong_signal_count = 0
    support_reasons: list[str] = []
    penalty_reasons: list[str] = []

    length_value = _to_nonnegative_int(feature_row.get("length"))
    coverage_value = _to_nonnegative_float(feature_row.get("coverage"))
    gc_content_value = _to_finite_float(feature_row.get("gc_content"))

    component_size = _to_nonnegative_int(component_row.get("component_size"))
    component_edge_count = _to_nonnegative_int(component_row.get("component_edge_count"))
    component_branching_nodes = _to_nonnegative_int(component_row.get("component_branching_nodes"))
    component_total_length = _to_nonnegative_int(component_row.get("component_total_length"))

    is_circular = bool(component_row.get("is_circular"))
    is_isolated = bool(component_row.get("is_isolated"))
    has_self_loop = bool(component_row.get("has_self_loop"))
    is_branching = bool(component_row.get("is_branching"))
    size_bin = determine_size_bin(length_value, config)
    is_compact_component = _is_compact_component(
        component_size=component_size,
        component_edge_count=component_edge_count,
        component_branching_nodes=component_branching_nodes,
        component_total_length=component_total_length,
        config=config,
    )

    if has_self_loop:
        score += 3
        total_signal_count += 1
        strong_signal_count += 1
        topology_signal_count += 1
        topology_strong_signal_count += 1
        support_reasons.append("graph self-loop supports circular topology")
    elif is_circular:
        score += 3
        total_signal_count += 1
        strong_signal_count += 1
        topology_signal_count += 1
        topology_strong_signal_count += 1
        support_reasons.append("graph component forms a simple directed cycle")

    if is_compact_component:
        score += 1
        total_signal_count += 1
        topology_signal_count += 1
        support_reasons.append("graph component is compact and non-branching")

    if is_isolated:
        if length_value is None or length_value <= config.isolated_length_cap:
            score += 1
            total_signal_count += 1
            topology_signal_count += 1
            support_reasons.append("graph component is an isolated singleton")
        else:
            penalty_reasons.append(
                "isolated status alone was not credited because the contig is unusually large"
            )

    has_coverage_signal = False
    if (
        coverage_value is not None
        and nuclear_coverage is not None
        and coverage_value > config.coverage_multiplier_threshold * nuclear_coverage
    ):
        score += 2
        total_signal_count += 1
        strong_signal_count += 1
        has_coverage_signal = True
        support_reasons.append(
            "coverage exceeds the nuclear baseline by the configured multiplier"
        )

    has_gc_signal = False
    if (
        gc_content_value is not None
        and gc_baseline is not None
        and 0 <= gc_content_value <= 100
        and 0 <= gc_baseline <= 100
        and abs(gc_content_value - gc_baseline) > config.gc_deviation_threshold
    ):
        score += 1
        total_signal_count += 1
        has_gc_signal = True
        support_reasons.append("GC content deviates from the assembly baseline")

    has_preferred_length_signal = False
    if length_value is not None:
        if size_bin == SIZE_BIN_PREFERRED_RANGE:
            score += 1
            total_signal_count += 1
            has_preferred_length_signal = True
            support_reasons.append("length is within the preferred candidate range")
        elif size_bin == SIZE_BIN_VERY_SMALL:
            penalty_reasons.append("very small contigs are weak graph-only candidates")
        elif size_bin == SIZE_BIN_SMALL_FRAGMENT:
            penalty_reasons.append("small fragments need BLAST support for confident recovery")

        if length_value >= config.large_contig_length:
            score -= config.large_contig_penalty
            penalty_reasons.append(
                f"length {_format_value(length_value)} exceeds the large-contig threshold"
            )
        if length_value > config.very_large_contig_length:
            score -= config.very_large_contig_penalty
            penalty_reasons.append(
                f"length {_format_value(length_value)} exceeds the very-large-contig threshold"
            )

    if is_branching or (component_branching_nodes is not None and component_branching_nodes > 0):
        score -= 1
        penalty_reasons.append("branching topology reduces organelle specificity")

    return {
        "graph_score": score,
        "support_signal_count": total_signal_count,
        "strong_signal_count": strong_signal_count,
        "topology_signal_count": topology_signal_count,
        "topology_strong_signal_count": topology_strong_signal_count,
        "graph_support_reasons": support_reasons,
        "graph_penalty_reasons": penalty_reasons,
        "is_compact_component": is_compact_component,
        "size_bin": size_bin,
        "has_gc_signal": has_gc_signal,
        "has_coverage_signal": has_coverage_signal,
        "has_preferred_length_signal": has_preferred_length_signal,
        "has_preferred_range_promotion_pattern": _has_preferred_range_correlated_pattern(
            size_bin=size_bin,
            is_compact_component=is_compact_component,
            is_isolated=is_isolated,
            has_gc_signal=has_gc_signal,
        ),
        "is_large": size_bin in {SIZE_BIN_LARGE, SIZE_BIN_VERY_LARGE},
    }


def evaluate_blast_features(
    *,
    blast_row: dict[str, object] | None,
) -> dict[str, object]:
    if blast_row is None:
        return {
            "blast_support_level": BLAST_SUPPORT_NONE,
            "blast_best_identity": None,
            "blast_best_aligned_fraction": None,
            "blast_merged_aligned_bp": 0,
            "blast_merged_coverage_fraction": 0.0,
            "blast_support_sources": "",
            "blast_selected_source": None,
            "plastid_support_level": BLAST_SUPPORT_NONE,
            "plastid_best_identity": None,
            "plastid_merged_aligned_bp": 0,
            "plastid_merged_coverage_fraction": 0.0,
            "mit_support_level": BLAST_SUPPORT_NONE,
            "mit_best_identity": None,
            "mit_merged_aligned_bp": 0,
            "mit_merged_coverage_fraction": 0.0,
            "blast_score": 0,
            "has_strong_blast_support": False,
            "has_moderate_blast_support": False,
            "has_exceptional_blast_support": False,
            "blast_support_reasons": [],
        }

    reasons: list[str] = []
    support_level = str(blast_row.get("blast_support_level", BLAST_SUPPORT_NONE))
    if support_level == BLAST_SUPPORT_STRONG:
        reasons.append("strong BLAST support")
    elif support_level == BLAST_SUPPORT_MODERATE:
        reasons.append("moderate BLAST support")
    support_text = blast_row.get("blast_supporting_evidence")
    if isinstance(support_text, str) and support_text:
        reasons.append(support_text)

    return {
        "blast_support_level": support_level,
        "blast_best_identity": blast_row.get("blast_best_identity"),
        "blast_best_aligned_fraction": blast_row.get("blast_best_aligned_fraction"),
        "blast_merged_aligned_bp": int(blast_row.get("blast_merged_aligned_bp", 0)),
        "blast_merged_coverage_fraction": blast_row.get("blast_merged_coverage_fraction"),
        "blast_support_sources": blast_row.get("blast_support_sources", ""),
        "blast_selected_source": blast_row.get("blast_selected_source"),
        "plastid_support_level": blast_row.get("plastid_support_level", BLAST_SUPPORT_NONE),
        "plastid_best_identity": blast_row.get("plastid_best_identity"),
        "plastid_merged_aligned_bp": int(blast_row.get("plastid_merged_aligned_bp", 0)),
        "plastid_merged_coverage_fraction": blast_row.get("plastid_merged_coverage_fraction"),
        "mit_support_level": blast_row.get("mit_support_level", BLAST_SUPPORT_NONE),
        "mit_best_identity": blast_row.get("mit_best_identity"),
        "mit_merged_aligned_bp": int(blast_row.get("mit_merged_aligned_bp", 0)),
        "mit_merged_coverage_fraction": blast_row.get("mit_merged_coverage_fraction"),
        "blast_score": int(blast_row.get("blast_score", 0)),
        "has_strong_blast_support": bool(blast_row.get("has_strong_blast_support")),
        "has_moderate_blast_support": bool(blast_row.get("has_moderate_blast_support")),
        "has_exceptional_blast_support": bool(blast_row.get("has_exceptional_blast_support")),
        "blast_support_reasons": reasons,
    }


def determine_graph_only_tier(
    graph_features: dict[str, object],
    config: ScoringConfig,
) -> tuple[str, list[str]]:
    score = int(graph_features["graph_score"])
    support_signal_count = int(graph_features["support_signal_count"])
    strong_signal_count = int(graph_features["strong_signal_count"])
    topology_signal_count = int(graph_features["topology_signal_count"])
    topology_strong_signal_count = int(graph_features["topology_strong_signal_count"])
    size_bin = str(graph_features["size_bin"])
    is_large = bool(graph_features["is_large"])
    has_preferred_range_promotion_pattern = bool(
        graph_features["has_preferred_range_promotion_pattern"]
    )
    notes: list[str] = []

    if topology_signal_count == 0:
        notes.append("graph-only mode found no graph-topology support")
        return NOT_FLAGGED, notes

    if size_bin == SIZE_BIN_VERY_SMALL:
        notes.append("very small contigs are not promoted by graph evidence alone")
        return NOT_FLAGGED, notes

    if size_bin == SIZE_BIN_SMALL_FRAGMENT and topology_strong_signal_count == 0:
        if score >= config.low_confidence_threshold:
            notes.append("small fragments are capped at low confidence without strong topology")
            return LOW_CONFIDENCE, notes
        notes.append("small fragments need stronger graph structure or BLAST support")
        return NOT_FLAGGED, notes

    if (
        size_bin == SIZE_BIN_PREFERRED_RANGE
        and not is_large
        and has_preferred_range_promotion_pattern
        and support_signal_count >= 4
    ):
        notes.append(
            "preferred-range compact singleton pattern with GC deviation promoted the graph-only call"
        )
        return HIGH_CONFIDENCE, notes

    if score >= config.high_confidence_threshold and support_signal_count >= 2:
        if is_large:
            if (
                support_signal_count >= config.large_contig_high_min_signals
                and topology_strong_signal_count >= 1
                and strong_signal_count >= 2
            ):
                return HIGH_CONFIDENCE, notes
            notes.append(
                "large contigs need multiple strong signals and a strong structural signal for high confidence"
            )
        elif topology_strong_signal_count >= 1:
            return HIGH_CONFIDENCE, notes
        else:
            notes.append(
                "high confidence requires a strong structural signal such as a self-loop or simple cycle"
            )

    if score >= config.medium_confidence_threshold and support_signal_count >= 2:
        if size_bin == SIZE_BIN_SMALL_FRAGMENT:
            notes.append("small fragments are capped at low confidence in graph-only mode")
            return LOW_CONFIDENCE, notes
        if is_large:
            if (
                support_signal_count >= config.large_contig_medium_min_signals
                and topology_strong_signal_count >= 1
            ):
                return MEDIUM_CONFIDENCE, notes
            notes.append(
                "large contigs were downgraded because weak topology alone is not sufficient"
            )
        else:
            return MEDIUM_CONFIDENCE, notes

    if score >= config.low_confidence_threshold:
        if size_bin == SIZE_BIN_VERY_LARGE and topology_strong_signal_count == 0:
            notes.append("very large contigs without strong structure stay suppressed")
            return NOT_FLAGGED, notes
        return LOW_CONFIDENCE, notes

    notes.append("combined graph evidence did not reach the review threshold")
    return NOT_FLAGGED, notes


def determine_blast_only_tier(
    *,
    size_bin: str,
    blast_features: dict[str, object],
) -> tuple[str, list[str]]:
    notes: list[str] = []
    has_exceptional = bool(blast_features["has_exceptional_blast_support"])
    has_strong = bool(blast_features["has_strong_blast_support"])
    has_moderate = bool(blast_features["has_moderate_blast_support"])

    if not has_moderate:
        notes.append("no BLAST support was detected")
        return NOT_FLAGGED, notes

    if size_bin == SIZE_BIN_VERY_SMALL:
        if has_strong:
            notes.append("strong BLAST support rescued a very small candidate")
            return MEDIUM_CONFIDENCE, notes
        return LOW_CONFIDENCE, notes

    if size_bin == SIZE_BIN_SMALL_FRAGMENT:
        if has_strong:
            return MEDIUM_CONFIDENCE, notes
        return LOW_CONFIDENCE, notes

    if size_bin == SIZE_BIN_PREFERRED_RANGE:
        if has_strong:
            return HIGH_CONFIDENCE, notes
        return MEDIUM_CONFIDENCE, notes

    if size_bin == SIZE_BIN_INTERMEDIATE:
        if has_strong:
            return MEDIUM_CONFIDENCE, notes
        return LOW_CONFIDENCE, notes

    if size_bin == SIZE_BIN_LARGE:
        if has_exceptional:
            notes.append("exceptional BLAST support kept a large contig under review")
            return MEDIUM_CONFIDENCE, notes
        if has_strong:
            notes.append("large contigs with BLAST support are capped below high confidence")
            return LOW_CONFIDENCE, notes
        return NOT_FLAGGED, notes

    if has_exceptional:
        notes.append("very large contigs remain capped even with exceptional BLAST support")
        return LOW_CONFIDENCE, notes
    notes.append("very large contigs require more than BLAST-only evidence by default")
    return NOT_FLAGGED, notes


def determine_hybrid_tier(
    *,
    size_bin: str,
    graph_features: dict[str, object],
    blast_features: dict[str, object],
    graph_only_tier: str,
    config: ScoringConfig,
) -> tuple[str, list[str]]:
    notes: list[str] = []
    has_exceptional_blast = bool(blast_features["has_exceptional_blast_support"])
    has_strong_blast = bool(blast_features["has_strong_blast_support"])
    has_moderate_blast = bool(blast_features["has_moderate_blast_support"])
    graph_has_support = int(graph_features["support_signal_count"]) > 0
    graph_has_strong_topology = int(graph_features["topology_strong_signal_count"]) > 0
    graph_multi_signal = int(graph_features["support_signal_count"]) >= 2
    preferred_range_promotion = bool(graph_features["has_preferred_range_promotion_pattern"])

    if size_bin == SIZE_BIN_VERY_SMALL:
        if has_strong_blast and graph_has_strong_topology:
            notes.append("hybrid mode used both strong BLAST support and strong graph structure")
            return HIGH_CONFIDENCE, notes
        if has_strong_blast:
            notes.append("hybrid mode recovered a very small fragment using BLAST support")
            return MEDIUM_CONFIDENCE, notes
        if has_moderate_blast:
            return LOW_CONFIDENCE, notes
        notes.append("very small contigs need BLAST support in hybrid mode")
        return NOT_FLAGGED, notes

    if size_bin == SIZE_BIN_SMALL_FRAGMENT:
        if has_strong_blast and graph_has_support:
            notes.append("hybrid mode combined BLAST support with graph corroboration for a fragment")
            return HIGH_CONFIDENCE, notes
        if has_strong_blast:
            notes.append("hybrid mode recovered a small fragment primarily from BLAST support")
            return MEDIUM_CONFIDENCE, notes
        if has_moderate_blast:
            return LOW_CONFIDENCE, notes
        notes.append("small fragments are not promoted without BLAST support")
        return NOT_FLAGGED, notes

    if size_bin == SIZE_BIN_PREFERRED_RANGE:
        if has_strong_blast:
            if graph_has_support:
                notes.append("hybrid mode promoted this preferred-range contig using BLAST and graph evidence")
            else:
                notes.append("hybrid mode promoted this preferred-range contig primarily from strong BLAST support")
            return HIGH_CONFIDENCE, notes
        if has_moderate_blast:
            if graph_has_support or preferred_range_promotion:
                notes.append("hybrid mode promoted a preferred-range contig through BLAST plus graph corroboration")
                return HIGH_CONFIDENCE, notes
            return MEDIUM_CONFIDENCE, notes
        if graph_only_tier == HIGH_CONFIDENCE:
            if graph_has_strong_topology:
                notes.append(
                    "hybrid mode allowed a graph-only rescue in the preferred range because strong topology was present"
                )
                return HIGH_CONFIDENCE, notes
            notes.append(
                "hybrid mode capped this BLAST-negative preferred-range candidate below high confidence because singleton-style graph evidence lacked strong topology"
            )
            return MEDIUM_CONFIDENCE, notes
        return graph_only_tier, notes

    if size_bin == SIZE_BIN_INTERMEDIATE:
        if has_strong_blast:
            if graph_has_strong_topology or graph_only_tier == HIGH_CONFIDENCE:
                notes.append("hybrid mode required BLAST plus strong graph corroboration for intermediate size")
                return HIGH_CONFIDENCE, notes
            return MEDIUM_CONFIDENCE, notes
        if has_moderate_blast:
            if graph_has_support:
                return MEDIUM_CONFIDENCE, notes
            return LOW_CONFIDENCE, notes
        return _cap_tier(graph_only_tier, MEDIUM_CONFIDENCE), notes

    if size_bin == SIZE_BIN_LARGE:
        if has_strong_blast:
            if (
                graph_has_strong_topology
                and int(graph_features["strong_signal_count"]) >= 2
                and int(graph_features["support_signal_count"]) >= config.large_contig_high_min_signals
                and has_exceptional_blast
            ):
                notes.append("hybrid mode accepted a large contig only with exceptional BLAST and strong graph support")
                return HIGH_CONFIDENCE, notes
            if graph_has_support:
                notes.append("hybrid mode kept a large BLAST-supported contig under review")
                return MEDIUM_CONFIDENCE, notes
            return LOW_CONFIDENCE, notes
        if has_moderate_blast:
            if graph_has_support:
                return LOW_CONFIDENCE, notes
            return NOT_FLAGGED, notes
        if graph_only_tier in {HIGH_CONFIDENCE, MEDIUM_CONFIDENCE}:
            notes.append("large graph-only candidates are capped below high confidence in hybrid mode")
            return MEDIUM_CONFIDENCE, notes
        return graph_only_tier, notes

    if not has_moderate_blast:
        notes.append("very large contigs without BLAST support remain suppressed")
        return NOT_FLAGGED, notes
    if has_exceptional_blast and graph_has_strong_topology and graph_multi_signal:
        notes.append("very large contigs remain capped even with exceptional combined evidence")
        return MEDIUM_CONFIDENCE, notes
    if has_strong_blast and graph_has_strong_topology:
        notes.append("very large contigs with mixed evidence stay in review")
        return LOW_CONFIDENCE, notes
    notes.append("very large contigs require exceptional corroboration in hybrid mode")
    return NOT_FLAGGED, notes


def _blank_blast_row(contig_id: str) -> dict[str, object]:
    return {
        "contig_id": contig_id,
        "blast_support_level": BLAST_SUPPORT_NONE,
        "blast_best_identity": None,
        "blast_best_aligned_fraction": None,
        "blast_score": 0,
        "has_strong_blast_support": False,
        "has_moderate_blast_support": False,
        "has_exceptional_blast_support": False,
        "blast_supporting_evidence": "",
    }


def evaluate_contig(
    *,
    contig_id: str,
    mode: str,
    feature_row: dict[str, object],
    component_row: dict[str, object],
    blast_row: dict[str, object] | None,
    nuclear_coverage: float | None,
    gc_baseline: float | None,
    config: ScoringConfig,
) -> ScoreRow:
    graph_features = evaluate_graph_features(
        feature_row=feature_row,
        component_row=component_row,
        nuclear_coverage=nuclear_coverage,
        gc_baseline=gc_baseline,
        config=config,
    )
    blast_features = evaluate_blast_features(blast_row=blast_row)

    graph_only_tier, graph_only_notes = determine_graph_only_tier(graph_features, config)

    if mode == MODE_GRAPH_ONLY:
        confidence_tier = graph_only_tier
        decision_notes = graph_only_notes
    elif mode == MODE_BLAST_ONLY:
        confidence_tier, decision_notes = determine_blast_only_tier(
            size_bin=str(graph_features["size_bin"]),
            blast_features=blast_features,
        )
    else:
        confidence_tier, decision_notes = determine_hybrid_tier(
            size_bin=str(graph_features["size_bin"]),
            graph_features=graph_features,
            blast_features=blast_features,
            graph_only_tier=graph_only_tier,
            config=config,
        )

    final_score = int(graph_features["graph_score"]) + int(blast_features["blast_score"])
    support_reasons = list(graph_features["graph_support_reasons"])
    support_reasons.extend(blast_features["blast_support_reasons"])
    penalty_reasons = list(graph_features["graph_penalty_reasons"])

    classification = classify_contig(confidence_tier)
    return {
        "contig_id": contig_id,
        "mode": mode,
        "size_bin": graph_features["size_bin"],
        "graph_score": int(graph_features["graph_score"]),
        "blast_score": int(blast_features["blast_score"]),
        "final_score": final_score,
        "confidence_tier": confidence_tier,
        "classification": classification,
        "support_signal_count": int(graph_features["support_signal_count"]),
        "strong_signal_count": int(graph_features["strong_signal_count"]),
        "topology_signal_count": int(graph_features["topology_signal_count"]),
        "topology_strong_signal_count": int(graph_features["topology_strong_signal_count"]),
        "is_compact_component": bool(graph_features["is_compact_component"]),
        "blast_support_level": blast_features["blast_support_level"],
        "blast_best_identity": blast_features["blast_best_identity"],
        "blast_best_aligned_fraction": blast_features["blast_best_aligned_fraction"],
        "blast_merged_aligned_bp": blast_features["blast_merged_aligned_bp"],
        "blast_merged_coverage_fraction": blast_features["blast_merged_coverage_fraction"],
        "blast_support_sources": blast_features["blast_support_sources"],
        "blast_selected_source": blast_features["blast_selected_source"],
        "plastid_support_level": blast_features["plastid_support_level"],
        "plastid_best_identity": blast_features["plastid_best_identity"],
        "plastid_merged_aligned_bp": blast_features["plastid_merged_aligned_bp"],
        "plastid_merged_coverage_fraction": blast_features["plastid_merged_coverage_fraction"],
        "mit_support_level": blast_features["mit_support_level"],
        "mit_best_identity": blast_features["mit_best_identity"],
        "mit_merged_aligned_bp": blast_features["mit_merged_aligned_bp"],
        "mit_merged_coverage_fraction": blast_features["mit_merged_coverage_fraction"],
        "supporting_evidence": "; ".join(support_reasons),
        "penalty_reasons": "; ".join(penalty_reasons),
        "decision_explanation": _build_decision_explanation(
            mode=mode,
            confidence_tier=confidence_tier,
            support_reasons=support_reasons,
            penalty_reasons=penalty_reasons,
            decision_notes=decision_notes,
        ),
    }


def score_contigs(
    contigs: dict[str, dict[str, object]],
    graph: nx.MultiDiGraph | None,
    *,
    mode: str = MODE_GRAPH_ONLY,
    config: ScoringConfig | None = None,
    blast_support_by_contig: dict[str, dict[str, object]] | None = None,
    nuclear_coverage: float | None = None,
    gc_baseline: float | None = None,
) -> list[ScoreRow]:
    if config is None:
        config = ScoringConfig()
    if mode not in VALID_MODES:
        raise ValueError(f"Unsupported mode {mode!r}; expected one of {', '.join(VALID_MODES)}")
    if mode in {MODE_HYBRID, MODE_BLAST_ONLY} and blast_support_by_contig is None:
        raise ValueError(
            f"{mode} mode requires BLAST support input; provide blast_support_by_contig, even if it is empty after BLAST filtering."
        )
    if graph is None:
        graph = nx.MultiDiGraph()
        for contig_id, record in contigs.items():
            graph.add_node(
                contig_id,
                name=contig_id,
                length=record.get("length"),
                coverage=record.get("coverage"),
            )
    if blast_support_by_contig is None:
        blast_support_by_contig = {}

    feature_rows = calculate_contig_features(contigs)
    feature_by_contig = {str(row["contig_id"]): row for row in feature_rows}

    if nuclear_coverage is None:
        nuclear_coverage = estimate_nuclear_coverage(contigs)
    nuclear_coverage = _normalize_nuclear_coverage_baseline(nuclear_coverage)

    if gc_baseline is None:
        gc_baseline = estimate_gc_baseline(row["gc_content"] for row in feature_rows)

    if nuclear_coverage is None and gc_baseline is None and mode != MODE_BLAST_ONLY:
        print(
            "Warning: nuclear coverage and GC baselines are unavailable; "
            "graph-based classification reliability is reduced.",
            file=sys.stderr,
        )

    component_by_contig = {str(row["contig_id"]): row for row in summarize_contigs(graph)}
    rows: list[ScoreRow] = []
    for contig_id in sorted(contigs):
        rows.append(
            evaluate_contig(
                contig_id=contig_id,
                mode=mode,
                feature_row=feature_by_contig.get(contig_id, {}),
                component_row=component_by_contig.get(contig_id, {}),
                blast_row=blast_support_by_contig.get(contig_id),
                nuclear_coverage=nuclear_coverage,
                gc_baseline=gc_baseline,
                config=config,
            )
        )
    return rows


def write_score_table(rows: Iterable[ScoreRow]) -> None:
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "contig_id",
            "mode",
            "size_bin",
            "graph_score",
            "blast_score",
            "final_score",
            "blast_support_level",
            "blast_merged_coverage_fraction",
            "confidence_tier",
            "classification",
            "supporting_evidence",
            "penalty_reasons",
            "decision_explanation",
        ]
    )
    for row in rows:
        writer.writerow(
            [
                row["contig_id"],
                row["mode"],
                row["size_bin"],
                row["graph_score"],
                row["blast_score"],
                row["final_score"],
                row["blast_support_level"],
                row.get("blast_merged_coverage_fraction"),
                row["confidence_tier"],
                row["classification"],
                row["supporting_evidence"],
                row["penalty_reasons"],
                row["decision_explanation"],
            ]
        )


def add_scoring_arguments(
    parser: argparse.ArgumentParser,
    *,
    mode_group: argparse._ArgumentGroup | None = None,
    blast_group: argparse._ArgumentGroup | None = None,
    advanced_group: argparse._ArgumentGroup | None = None,
) -> None:
    mode_arguments = mode_group or parser
    blast_arguments = blast_group or parser
    advanced_arguments = advanced_group or parser

    mode_arguments.add_argument(
        "--mode",
        choices=VALID_MODES,
        default=MODE_GRAPH_ONLY,
        help=(
            "Analysis mode: graph-only is the default; hybrid is usually the best practical choice when organelle FASTA inputs are available; blast-only uses only BLAST-derived scoring."
        ),
    )
    blast_arguments.add_argument(
        "--plastid-fasta",
        dest="plastid_fasta",
        type=Path,
        default=None,
        help="Plastid or chloroplast reference FASTA for internal BLAST.",
    )
    blast_arguments.add_argument(
        "--mit-fasta",
        type=Path,
        default=None,
        help="Mitochondrial reference FASTA for internal BLAST.",
    )
    advanced_arguments.add_argument(
        "--nuclear-coverage",
        type=float,
        default=None,
        help="Override the inferred nuclear coverage baseline.",
    )
    advanced_arguments.add_argument(
        "--gc-baseline",
        type=float,
        default=None,
        help="Override the inferred GC baseline percentage.",
    )
    advanced_arguments.add_argument(
        "--high-confidence-threshold",
        "--organelle-threshold",
        dest="high_confidence_threshold",
        type=int,
        default=DEFAULT_HIGH_CONFIDENCE_THRESHOLD,
        help="Minimum graph score for a high-confidence graph call.",
    )
    advanced_arguments.add_argument(
        "--medium-confidence-threshold",
        type=int,
        default=DEFAULT_MEDIUM_CONFIDENCE_THRESHOLD,
        help="Minimum graph score for a medium-confidence graph call.",
    )
    advanced_arguments.add_argument(
        "--low-confidence-threshold",
        type=int,
        default=DEFAULT_LOW_CONFIDENCE_THRESHOLD,
        help="Minimum graph score for a low-confidence graph call.",
    )
    advanced_arguments.add_argument(
        "--very-small-max-length",
        type=int,
        default=DEFAULT_VERY_SMALL_MAX,
        help="Upper bound for the very-small contig bin.",
    )
    advanced_arguments.add_argument(
        "--small-fragment-max-length",
        type=int,
        default=DEFAULT_SMALL_FRAGMENT_MAX,
        help="Upper bound for the small-fragment contig bin.",
    )
    advanced_arguments.add_argument(
        "--min-organelle-length",
        type=int,
        default=DEFAULT_MIN_ORGANELLE_LENGTH,
        help="Lower bound of the preferred-size range.",
    )
    advanced_arguments.add_argument(
        "--preferred-max-length",
        type=int,
        default=DEFAULT_PREFERRED_MAX_LENGTH,
        help="Upper bound of the preferred-size range.",
    )
    advanced_arguments.add_argument(
        "--intermediate-max-length",
        type=int,
        default=DEFAULT_INTERMEDIATE_MAX,
        help="Upper bound of the intermediate-size bin.",
    )
    advanced_arguments.add_argument(
        "--large-contig-length",
        type=int,
        default=DEFAULT_LARGE_CONTIG_LENGTH,
        help="Lower bound of the large-contig bin.",
    )
    advanced_arguments.add_argument(
        "--very-large-contig-length",
        type=int,
        default=DEFAULT_VERY_LARGE_CONTIG_LENGTH,
        help="Lower bound of the very-large-contig bin.",
    )
    advanced_arguments.add_argument(
        "--isolated-length-cap",
        type=int,
        default=DEFAULT_ISOLATED_LENGTH_CAP,
        help="Maximum length for isolated singleton status to count as positive graph evidence.",
    )
    advanced_arguments.add_argument(
        "--large-contig-penalty",
        type=int,
        default=DEFAULT_LARGE_CONTIG_PENALTY,
        help="Penalty applied once a contig exceeds --large-contig-length.",
    )
    advanced_arguments.add_argument(
        "--very-large-contig-penalty",
        type=int,
        default=DEFAULT_VERY_LARGE_CONTIG_PENALTY,
        help="Additional penalty applied once a contig exceeds --very-large-contig-length.",
    )
    advanced_arguments.add_argument(
        "--large-contig-high-min-signals",
        type=int,
        default=DEFAULT_LARGE_CONTIG_HIGH_MIN_SIGNALS,
        help="Minimum graph signals for large contigs to reach high confidence.",
    )
    advanced_arguments.add_argument(
        "--large-contig-medium-min-signals",
        type=int,
        default=DEFAULT_LARGE_CONTIG_MEDIUM_MIN_SIGNALS,
        help="Minimum graph signals for large contigs to reach medium confidence.",
    )
    advanced_arguments.add_argument(
        "--gc-deviation-threshold",
        type=float,
        default=DEFAULT_GC_DEVIATION_THRESHOLD,
        help="Minimum GC deviation from baseline required for positive graph evidence.",
    )
    advanced_arguments.add_argument(
        "--coverage-multiplier-threshold",
        type=float,
        default=DEFAULT_COVERAGE_MULTIPLIER_THRESHOLD,
        help="Coverage multiplier over the nuclear baseline required for positive graph evidence.",
    )


def config_from_args(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser | None = None,
) -> tuple[ScoringConfig, BlastThresholdConfig, float | None, float | None]:
    def fail(message: str) -> None:
        if parser is not None:
            parser.error(message)
        raise ValueError(message)

    if args.high_confidence_threshold < 0:
        fail("High confidence threshold must be non-negative")
    if args.medium_confidence_threshold < 0:
        fail("Medium confidence threshold must be non-negative")
    if args.low_confidence_threshold < 0:
        fail("Low confidence threshold must be non-negative")
    if args.high_confidence_threshold < args.medium_confidence_threshold:
        fail("High confidence threshold must be greater than or equal to medium confidence threshold")
    if args.medium_confidence_threshold < args.low_confidence_threshold:
        fail("Medium confidence threshold must be greater than or equal to low confidence threshold")

    for field_name in (
        "very_small_max_length",
        "small_fragment_max_length",
        "min_organelle_length",
        "preferred_max_length",
        "intermediate_max_length",
        "large_contig_length",
        "very_large_contig_length",
        "isolated_length_cap",
        "large_contig_penalty",
        "very_large_contig_penalty",
        "large_contig_high_min_signals",
        "large_contig_medium_min_signals",
    ):
        if getattr(args, field_name) < 0:
            fail(f"{field_name.replace('_', '-')} must be non-negative")

    if args.very_small_max_length >= args.small_fragment_max_length:
        fail("--very-small-max-length must be less than --small-fragment-max-length")
    if args.small_fragment_max_length >= args.min_organelle_length:
        fail("--small-fragment-max-length must be less than --min-organelle-length")
    if args.preferred_max_length < args.min_organelle_length:
        fail("--preferred-max-length must be greater than or equal to --min-organelle-length")
    if args.intermediate_max_length < args.preferred_max_length:
        fail("--intermediate-max-length must be greater than or equal to --preferred-max-length")
    if args.large_contig_length < args.intermediate_max_length:
        fail("--large-contig-length must be greater than or equal to --intermediate-max-length")
    if args.very_large_contig_length < args.large_contig_length:
        fail("--very-large-contig-length must be greater than or equal to --large-contig-length")

    nuclear_coverage = None
    if args.nuclear_coverage is not None:
        nuclear_coverage = _normalize_nuclear_coverage_baseline(args.nuclear_coverage)
        if nuclear_coverage is None:
            fail("Nuclear coverage must be a finite positive value")

    gc_baseline = None
    if args.gc_baseline is not None:
        gc_baseline = _to_finite_float(args.gc_baseline)
        if gc_baseline is None or gc_baseline < 0 or gc_baseline > 100:
            fail("GC baseline must be a finite percentage between 0 and 100")

    if _to_nonnegative_float(args.gc_deviation_threshold) is None:
        fail("GC deviation threshold must be a finite non-negative value")
    if _normalize_nuclear_coverage_baseline(args.coverage_multiplier_threshold) is None:
        fail("Coverage multiplier threshold must be a finite positive value")
    has_organelle_fasta = args.plastid_fasta is not None or args.mit_fasta is not None
    internal_blast_requested = args.mode in {MODE_HYBRID, MODE_BLAST_ONLY}
    if internal_blast_requested and not has_organelle_fasta:
        fail(
            f"{args.mode} mode requires at least one organelle FASTA input via --plastid-fasta or --mit-fasta"
        )
    if internal_blast_requested and args.assembly_fasta is None:
        fail("--assembly-fasta is required whenever --plastid-fasta or --mit-fasta is used")

    config = ScoringConfig(
        high_confidence_threshold=args.high_confidence_threshold,
        medium_confidence_threshold=args.medium_confidence_threshold,
        low_confidence_threshold=args.low_confidence_threshold,
        min_organelle_length=args.min_organelle_length,
        preferred_max_length=args.preferred_max_length,
        large_contig_length=args.large_contig_length,
        very_large_contig_length=args.very_large_contig_length,
        isolated_length_cap=args.isolated_length_cap,
        large_contig_penalty=args.large_contig_penalty,
        very_large_contig_penalty=args.very_large_contig_penalty,
        large_contig_high_min_signals=args.large_contig_high_min_signals,
        large_contig_medium_min_signals=args.large_contig_medium_min_signals,
        gc_deviation_threshold=args.gc_deviation_threshold,
        coverage_multiplier_threshold=args.coverage_multiplier_threshold,
        very_small_max_length=args.very_small_max_length,
        small_fragment_max_length=args.small_fragment_max_length,
        intermediate_max_length=args.intermediate_max_length,
    )
    return config, BlastThresholdConfig(), nuclear_coverage, gc_baseline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Score GFA contigs for graph-only, BLAST-only, or hybrid organelle-like evidence. graph-only remains the default CLI mode for backward compatibility."
    )
    parser.add_argument("input_gfa", type=Path, help="Path to the input GFA file")
    parser.add_argument(
        "--assembly-fasta",
        type=Path,
        default=None,
        help="Assembly FASTA required when internal BLAST needs sequence access.",
    )
    add_scoring_arguments(parser)
    return parser


def _run_self_checks() -> None:
    assert determine_size_bin(5_000, ScoringConfig()) == SIZE_BIN_VERY_SMALL
    assert determine_size_bin(80_000, ScoringConfig()) == SIZE_BIN_PREFERRED_RANGE
    assert determine_size_bin(700_000, ScoringConfig()) == SIZE_BIN_VERY_LARGE


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        config, blast_config, nuclear_coverage, gc_baseline = config_from_args(args, parser)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    try:
        contigs, graph = parse_gfa(args.input_gfa)
        blast_support_by_contig = load_blast_support_by_contig(
            plastid_blast_tsv=args.plastid_blast_tsv,
            mit_blast_tsv=args.mit_blast_tsv,
            config=blast_config,
        )
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (ImportError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1

    rows = score_contigs(
        contigs,
        graph,
        mode=args.mode,
        config=config,
        blast_support_by_contig=blast_support_by_contig,
        nuclear_coverage=nuclear_coverage,
        gc_baseline=gc_baseline,
    )
    write_score_table(rows)
    return 0


if __name__ == "__main__":
    _run_self_checks()
    raise SystemExit(main())
