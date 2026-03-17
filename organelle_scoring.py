#!/usr/bin/env python3
"""Score contigs for organelle-like properties."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path
from typing import Iterable

import networkx as nx

from graph_analysis import summarize_contigs
from parse_gfa import parse_gfa
from sequence_features import calculate_contig_features, estimate_nuclear_coverage


ScoreRow = dict[str, object]

DEFAULT_ORGANELLE_SCORE_THRESHOLD = 4
MIN_ORGANELLE_LENGTH = 50_000
MAX_ORGANELLE_LENGTH = 1_000_000
GC_DEVIATION_THRESHOLD = 5.0
COVERAGE_MULTIPLIER_THRESHOLD = 5.0
GC_TRIM_FRACTION = 0.10


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


def organelle_score(
    *,
    length: object,
    coverage: object,
    nuclear_coverage: object,
    gc_content: object,
    gc_baseline: object,
    is_circular: object,
    is_isolated: object,
) -> int:
    score = 0

    length_value = _to_finite_float(length)
    coverage_value = _to_nonnegative_float(coverage)
    nuclear_coverage_value = _normalize_nuclear_coverage_baseline(nuclear_coverage)
    gc_content_value = _to_finite_float(gc_content) if gc_content is not None else None
    gc_baseline_value = _to_finite_float(gc_baseline) if gc_baseline is not None else None

    # Coverage is assumed to be depth-normalized. If it comes from raw read-count
    # support such as GFA RC, coverage-based scoring may be systematically biased.
    if (
        coverage_value is not None
        and nuclear_coverage_value is not None
        and coverage_value > COVERAGE_MULTIPLIER_THRESHOLD * nuclear_coverage_value
    ):
        score += 2

    if bool(is_circular) or bool(is_isolated):
        score += 2

    if (
        length_value is not None
        and MIN_ORGANELLE_LENGTH < length_value < MAX_ORGANELLE_LENGTH
    ):
        score += 1

    if (
        gc_content_value is not None
        and gc_baseline_value is not None
        and 0 <= gc_content_value <= 100
        and 0 <= gc_baseline_value <= 100
        and abs(gc_content_value - gc_baseline_value) > GC_DEVIATION_THRESHOLD
    ):
        score += 1

    return score


def classify_contig(score: int, organelle_threshold: int) -> str:
    if score >= organelle_threshold:
        return "organelle"
    return "nuclear"


def score_contigs(
    contigs: dict[str, dict[str, object]],
    graph: nx.MultiDiGraph,
    *,
    organelle_threshold: int = DEFAULT_ORGANELLE_SCORE_THRESHOLD,
    nuclear_coverage: float | None = None,
    gc_baseline: float | None = None,
) -> list[ScoreRow]:
    feature_rows = calculate_contig_features(contigs)
    feature_by_contig = {
        str(row["contig_id"]): row
        for row in feature_rows
    }

    if nuclear_coverage is None:
        nuclear_coverage = estimate_nuclear_coverage(contigs)
    nuclear_coverage = _normalize_nuclear_coverage_baseline(nuclear_coverage)

    if gc_baseline is None:
        gc_baseline = estimate_gc_baseline(row["gc_content"] for row in feature_rows)

    if nuclear_coverage is None and gc_baseline is None:
        print(
            "Warning: nuclear coverage and GC baselines are unavailable; "
            "classification reliability is reduced.",
            file=sys.stderr,
        )

    component_by_contig = {str(row["contig_id"]): row for row in summarize_contigs(graph)}

    rows: list[ScoreRow] = []
    for contig_id in sorted(contigs):
        feature_row = feature_by_contig.get(contig_id, {})
        component_row = component_by_contig.get(contig_id, {})
        score = organelle_score(
            length=feature_row.get("length"),
            coverage=feature_row.get("coverage"),
            nuclear_coverage=nuclear_coverage,
            gc_content=feature_row.get("gc_content"),
            gc_baseline=gc_baseline,
            is_circular=component_row.get("is_circular", False),
            is_isolated=component_row.get("is_isolated", False),
        )
        rows.append(
            {
                "contig_id": contig_id,
                "organelle_score": score,
                "classification": classify_contig(score, organelle_threshold),
            }
        )

    return rows


def write_score_table(rows: Iterable[ScoreRow]) -> None:
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(["contig_id", "organelle_score", "classification"])
    for row in rows:
        writer.writerow(
            [row["contig_id"], row["organelle_score"], row["classification"]]
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Score GFA contigs for organelle-like features."
    )
    parser.add_argument("input_gfa", type=Path, help="Path to the input GFA file")
    parser.add_argument(
        "--nuclear-coverage",
        type=float,
        default=None,
        help="Override the nuclear coverage baseline",
    )
    parser.add_argument(
        "--gc-baseline",
        type=float,
        default=None,
        help="Override the GC baseline percentage",
    )
    parser.add_argument(
        "--organelle-threshold",
        type=int,
        default=DEFAULT_ORGANELLE_SCORE_THRESHOLD,
        help="Minimum score required to classify a contig as organelle",
    )
    return parser


def _run_self_checks() -> None:
    assert organelle_score(
        length=100_000,
        coverage=60,
        nuclear_coverage=10,
        gc_content=55,
        gc_baseline=45,
        is_circular=True,
        is_isolated=False,
    ) == 6
    assert organelle_score(
        length=40_000,
        coverage=None,
        nuclear_coverage=10,
        gc_content=50,
        gc_baseline=47,
        is_circular=False,
        is_isolated=True,
    ) == 2
    assert organelle_score(
        length=100_000,
        coverage=60,
        nuclear_coverage=0,
        gc_content=55,
        gc_baseline=45,
        is_circular=True,
        is_isolated=True,
    ) == 4
    assert estimate_gc_baseline([40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 80, 81]) == 45.5
    assert classify_contig(4, DEFAULT_ORGANELLE_SCORE_THRESHOLD) == "organelle"


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.organelle_threshold < 0:
        print("Organelle threshold must be non-negative", file=sys.stderr)
        return 1
    if args.nuclear_coverage is not None and _to_nonnegative_float(args.nuclear_coverage) is None:
        print("Nuclear coverage must be a finite non-negative value", file=sys.stderr)
        return 1
    if args.gc_baseline is not None:
        gc_baseline_value = _to_finite_float(args.gc_baseline)
        if gc_baseline_value is None or gc_baseline_value < 0 or gc_baseline_value > 100:
            print("GC baseline must be a finite percentage between 0 and 100", file=sys.stderr)
            return 1

    try:
        contigs, graph = parse_gfa(args.input_gfa)
    except FileNotFoundError:
        print(f"Input file not found: {args.input_gfa}", file=sys.stderr)
        return 1
    except (ImportError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1

    rows = score_contigs(
        contigs,
        graph,
        organelle_threshold=args.organelle_threshold,
        nuclear_coverage=args.nuclear_coverage,
        gc_baseline=args.gc_baseline,
    )

    write_score_table(rows)
    return 0


if __name__ == "__main__":
    _run_self_checks()
    raise SystemExit(main())
