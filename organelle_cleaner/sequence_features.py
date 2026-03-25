#!/usr/bin/env python3
"""Calculate simple sequence features for contigs."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path
from typing import Iterable


FeatureRow = dict[str, object]
ContigRecord = dict[str, object]


def gc_content(sequence: str | None) -> float | None:
    """Return GC fraction as a percentage."""
    if not sequence:
        return None

    normalized = sequence.upper()
    if not normalized:
        return None

    gc_count = sum(1 for base in normalized if base in {"G", "C"})
    return (gc_count / len(normalized)) * 100.0


def calculate_contig_features(contigs: dict[str, ContigRecord]) -> list[FeatureRow]:
    """Build a feature table from contig sequence records."""
    rows: list[FeatureRow] = []

    for contig_id in sorted(contigs):
        record = contigs[contig_id]
        sequence_value = record.get("sequence", "")
        sequence = sequence_value if isinstance(sequence_value, str) else ""
        length = record.get("length")
        if not isinstance(length, int):
            length = len(sequence)

        coverage = record.get("coverage")
        if not isinstance(coverage, (int, float)):
            coverage = None

        gc_value = gc_content(sequence)

        rows.append(
            {
                "contig_id": contig_id,
                "length": length,
                "gc_content": round(gc_value, 4) if gc_value is not None else None,
                "coverage": coverage,
            }
        )

    return rows


def estimate_nuclear_coverage_from_values(values: Iterable[float]) -> float | None:
    """Estimate nuclear coverage after removing the top 5% highest values, rounded up.

    Processing order:
    1. Convert each value to ``float``.
    2. Remove non-finite values such as NaN and inf.
    3. Sort the remaining values.
    4. Trim the top 5% highest coverage values, rounded up.
    5. Compute the median of the remaining values.

    For small datasets, at least one highest value may be removed. This
    intentionally removes high-coverage outliers such as organelle contigs.

    Examples:
        >>> estimate_nuclear_coverage_from_values([10, 11, 12, 13, 1000])
        11.5
        >>> estimate_nuclear_coverage_from_values(["10", 11.0, "bad", None, 1000])
        10.5
        >>> estimate_nuclear_coverage_from_values(["10", "nan", "20"])
        10.0
        >>> estimate_nuclear_coverage_from_values(["inf", "100", "200"])
        100.0
        >>> estimate_nuclear_coverage_from_values([None, "bad"])
        None
    """
    coverages: list[float] = []

    for value in values:
        try:
            coverage = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(coverage):
            coverages.append(coverage)

    if not coverages:
        return None

    coverages.sort()
    trim_count = math.ceil(len(coverages) * 0.05)
    if len(coverages) >= 2:
        trim_count = max(1, trim_count)

    trimmed_coverages = coverages[:-trim_count] if trim_count > 0 else coverages
    if not trimmed_coverages:
        return None

    if len(trimmed_coverages) == 1:
        return float(trimmed_coverages[0])

    return float(statistics.median(trimmed_coverages))


def estimate_nuclear_coverage(contigs: dict[str, ContigRecord]) -> float | None:
    """Estimate nuclear coverage from contig coverages, accepting numeric strings.

    The top 5% highest coverage values are removed using round-up trimming before
    taking the median. Invalid or non-numeric coverage values are skipped.
    """
    return estimate_nuclear_coverage_from_values(
        record.get("coverage") for record in contigs.values()
    )


def _run_self_checks() -> None:
    assert estimate_nuclear_coverage_from_values(["10", "nan", "20"]) == 10.0
    assert estimate_nuclear_coverage_from_values(["inf", "100", "200"]) == 100.0
    assert estimate_nuclear_coverage_from_values(["nan", "inf"]) is None
    assert estimate_nuclear_coverage_from_values([10, 20]) == 10.0


def read_fasta(fasta_path: Path) -> dict[str, ContigRecord]:
    """Read contig sequences from a FASTA file."""
    contigs: dict[str, ContigRecord] = {}
    current_id: str | None = None
    sequence_chunks: list[str] = []

    with fasta_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    sequence = "".join(sequence_chunks)
                    contigs[current_id] = {
                        "sequence": sequence,
                        "length": len(sequence),
                        "coverage": None,
                    }

                current_id = line[1:].split()[0]
                if current_id in contigs:
                    raise ValueError(f"Duplicate FASTA ID found: {current_id}")
                sequence_chunks = []
            else:
                if current_id is None:
                    raise ValueError("Invalid FASTA: sequence found before header")
                sequence_chunks.append(line)

    if current_id is not None:
        sequence = "".join(sequence_chunks)
        contigs[current_id] = {
            "sequence": sequence,
            "length": len(sequence),
            "coverage": None,
        }

    return contigs


def write_feature_table(rows: Iterable[FeatureRow]) -> None:
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(["contig_id", "length", "gc_content", "coverage"])
    for row in rows:
        writer.writerow(
            [row["contig_id"], row["length"], row["gc_content"], row["coverage"]]
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Calculate contig length, GC content, and coverage."
    )
    parser.add_argument(
        "input_fasta",
        type=Path,
        help="Path to a FASTA file containing contig sequences",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        contigs = read_fasta(args.input_fasta)
        rows = calculate_contig_features(contigs)
    except FileNotFoundError:
        print(f"Input file not found: {args.input_fasta}", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    write_feature_table(rows)
    return 0


if __name__ == "__main__":
    _run_self_checks()
    raise SystemExit(main())
