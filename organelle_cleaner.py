#!/usr/bin/env python3
"""Run the organelle-cleaning pipeline on a GFA assembly."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from graph_analysis import summarize_contigs
from organelle_scoring import score_contigs
from parse_gfa import parse_gfa
from sequence_features import calculate_contig_features


REPORT_NAME = "report.tsv"
ORGANELLE_CONTIGS_NAME = "organelle_contigs.txt"
NUCLEAR_CONTIGS_NAME = "nuclear_contigs.txt"
CLEANED_ASSEMBLY_NAME = "cleaned_assembly.fa"
FASTA_LINE_WIDTH = 80


def write_contig_id_list(path: Path, contig_ids: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for contig_id in contig_ids:
            handle.write(f"{contig_id}\n")


def _sequence_lines(sequence: str) -> list[str]:
    if not sequence:
        return [""]
    return [
        sequence[index : index + FASTA_LINE_WIDTH]
        for index in range(0, len(sequence), FASTA_LINE_WIDTH)
    ]


def write_fasta(path: Path, contigs: dict[str, dict[str, object]], contig_ids: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for contig_id in contig_ids:
            sequence_value = contigs[contig_id].get("sequence", "")
            sequence = sequence_value if isinstance(sequence_value, str) else ""
            if sequence == "":
                print(
                    f"Contig {contig_id} has no sequence; skipped from FASTA output",
                    file=sys.stderr,
                )
                continue
            handle.write(f">{contig_id}\n")
            for line in _sequence_lines(sequence):
                handle.write(f"{line}\n")


def write_report(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    fieldnames = [
        "contig_id",
        "length",
        "coverage",
        "gc_content",
        "component_id",
        "component_size",
        "is_circular",
        "is_isolated",
        "has_sequence",
        "organelle_score",
        "classification",
    ]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def build_pipeline_rows(
    contigs: dict[str, dict[str, object]],
    graph,
) -> list[dict[str, object]]:
    feature_rows = calculate_contig_features(contigs)
    feature_by_contig = {str(row["contig_id"]): row for row in feature_rows}
    topology_rows = summarize_contigs(graph)
    topology_by_contig = {str(row["contig_id"]): row for row in topology_rows}
    scoring_rows = score_contigs(contigs, graph)
    scoring_by_contig = {str(row["contig_id"]): row for row in scoring_rows}

    rows: list[dict[str, object]] = []
    for contig_id in sorted(contigs):
        feature_row = feature_by_contig.get(contig_id, {})
        topology_row = topology_by_contig.get(contig_id, {})
        scoring_row = scoring_by_contig.get(contig_id, {})
        rows.append(
            {
                "contig_id": contig_id,
                "length": feature_row.get("length"),
                "coverage": feature_row.get("coverage"),
                "gc_content": feature_row.get("gc_content"),
                "component_id": topology_row.get("component_id"),
                "component_size": topology_row.get("component_size"),
                "is_circular": topology_row.get("is_circular"),
                "is_isolated": topology_row.get("is_isolated"),
                "has_sequence": bool(contigs[contig_id].get("sequence")),
                "organelle_score": scoring_row.get("organelle_score"),
                "classification": scoring_row.get("classification"),
            }
        )

    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Parse, score, classify, and clean a GFA assembly."
    )
    parser.add_argument("input_gfa", type=Path, help="Path to the input assembly GFA")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for organelle_contigs.txt, nuclear_contigs.txt, cleaned_assembly.fa, and report.tsv",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        contigs, graph = parse_gfa(args.input_gfa)
    except FileNotFoundError:
        print(f"Input file not found: {args.input_gfa}", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    report_rows = build_pipeline_rows(contigs, graph)
    organelle_ids = [
        str(row["contig_id"])
        for row in report_rows
        if row.get("classification") == "organelle"
    ]
    nuclear_ids: list[str] = []
    for row in report_rows:
        if row.get("classification") != "nuclear":
            continue
        contig_id = str(row["contig_id"])
        if bool(row.get("has_sequence")):
            nuclear_ids.append(contig_id)
            continue
        print(
            f"Contig {contig_id} is classified as nuclear but has no sequence; excluded from nuclear outputs",
            file=sys.stderr,
        )

    write_contig_id_list(output_dir / ORGANELLE_CONTIGS_NAME, organelle_ids)
    write_contig_id_list(output_dir / NUCLEAR_CONTIGS_NAME, nuclear_ids)
    write_fasta(output_dir / CLEANED_ASSEMBLY_NAME, contigs, nuclear_ids)
    write_report(output_dir / REPORT_NAME, report_rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
