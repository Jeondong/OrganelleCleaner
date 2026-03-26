#!/usr/bin/env python3
"""Run the organelle-cleaning pipeline on a GFA assembly."""

from __future__ import annotations

import argparse
import csv
from multiprocessing import Pool
import sys
from pathlib import Path

from .blast_features import load_blast_support_by_contig
from .graph_analysis import summarize_contigs
from .internal_blast import run_internal_blast
from .organelle_scoring import (
    HIGH_CONFIDENCE,
    add_scoring_arguments,
    config_from_args,
    score_contigs,
)
from .parse_gfa import parse_gfa
from .sequence_features import calculate_contig_features, read_fasta


REPORT_NAME = "report.tsv"
ORGANELLE_CONTIGS_NAME = "organelle_contigs.txt"
NUCLEAR_CONTIGS_NAME = "nuclear_contigs.txt"
CLEANED_ASSEMBLY_NAME = "cleaned_assembly.fa"
DEFAULT_CANDIDATE_REPORT_NAME = "organelle_candidates.tsv"
FASTA_LINE_WIDTH = 80
_ROW_BUILD_CONTEXT: tuple[
    dict[str, dict[str, object]],
    dict[str, dict[str, object]],
    dict[str, dict[str, object]],
    dict[str, dict[str, object]],
] | None = None


class CliHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    """Formatter that keeps examples readable and shows defaults."""


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
        "mode",
        "size_bin",
        "contig_length",
        "coverage",
        "gc_content",
        "component_id",
        "component_size",
        "component_edge_count",
        "component_branching_nodes",
        "component_total_length",
        "is_circular",
        "is_isolated",
        "in_degree",
        "out_degree",
        "total_degree",
        "self_loop_count",
        "has_self_loop",
        "is_branching",
        "is_compact_component",
        "has_sequence",
        "blast_support_level",
        "blast_support_sources",
        "blast_selected_source",
        "blast_best_identity",
        "blast_best_aligned_fraction",
        "blast_merged_aligned_bp",
        "blast_merged_coverage_fraction",
        "plastid_support_level",
        "plastid_best_identity",
        "plastid_merged_aligned_bp",
        "plastid_merged_coverage_fraction",
        "mit_support_level",
        "mit_best_identity",
        "mit_merged_aligned_bp",
        "mit_merged_coverage_fraction",
        "graph_score",
        "blast_score",
        "final_score",
        "support_signal_count",
        "strong_signal_count",
        "topology_signal_count",
        "topology_strong_signal_count",
        "confidence_tier",
        "classification",
        "supporting_evidence",
        "penalty_reasons",
        "decision_explanation",
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


def _init_row_build_context(
    contigs: dict[str, dict[str, object]],
    feature_by_contig: dict[str, dict[str, object]],
    topology_by_contig: dict[str, dict[str, object]],
    scoring_by_contig: dict[str, dict[str, object]],
) -> None:
    global _ROW_BUILD_CONTEXT
    _ROW_BUILD_CONTEXT = (
        contigs,
        feature_by_contig,
        topology_by_contig,
        scoring_by_contig,
    )


def _build_pipeline_row(
    contig_id: str,
    contigs: dict[str, dict[str, object]] | None = None,
    feature_by_contig: dict[str, dict[str, object]] | None = None,
    topology_by_contig: dict[str, dict[str, object]] | None = None,
    scoring_by_contig: dict[str, dict[str, object]] | None = None,
) -> dict[str, object]:
    if contigs is None or feature_by_contig is None or topology_by_contig is None or scoring_by_contig is None:
        if _ROW_BUILD_CONTEXT is None:
            raise RuntimeError("Row build context was not initialized")
        contigs, feature_by_contig, topology_by_contig, scoring_by_contig = _ROW_BUILD_CONTEXT

    feature_row = feature_by_contig.get(contig_id, {})
    topology_row = topology_by_contig.get(contig_id, {})
    scoring_row = scoring_by_contig.get(contig_id, {})
    return {
        "contig_id": contig_id,
        "mode": scoring_row.get("mode"),
        "size_bin": scoring_row.get("size_bin"),
        "contig_length": feature_row.get("length"),
        "coverage": feature_row.get("coverage"),
        "gc_content": feature_row.get("gc_content"),
        "component_id": topology_row.get("component_id"),
        "component_size": topology_row.get("component_size"),
        "component_edge_count": topology_row.get("component_edge_count"),
        "component_branching_nodes": topology_row.get("component_branching_nodes"),
        "component_total_length": topology_row.get("component_total_length"),
        "is_circular": topology_row.get("is_circular"),
        "is_isolated": topology_row.get("is_isolated"),
        "in_degree": topology_row.get("in_degree"),
        "out_degree": topology_row.get("out_degree"),
        "total_degree": topology_row.get("total_degree"),
        "self_loop_count": topology_row.get("self_loop_count"),
        "has_self_loop": topology_row.get("has_self_loop"),
        "is_branching": topology_row.get("is_branching"),
        "is_compact_component": scoring_row.get("is_compact_component"),
        "has_sequence": bool(contigs[contig_id].get("sequence")),
        "blast_support_level": scoring_row.get("blast_support_level"),
        "blast_support_sources": scoring_row.get("blast_support_sources"),
        "blast_selected_source": scoring_row.get("blast_selected_source"),
        "blast_best_identity": scoring_row.get("blast_best_identity"),
        "blast_best_aligned_fraction": scoring_row.get("blast_best_aligned_fraction"),
        "blast_merged_aligned_bp": scoring_row.get("blast_merged_aligned_bp"),
        "blast_merged_coverage_fraction": scoring_row.get("blast_merged_coverage_fraction"),
        "plastid_support_level": scoring_row.get("plastid_support_level"),
        "plastid_best_identity": scoring_row.get("plastid_best_identity"),
        "plastid_merged_aligned_bp": scoring_row.get("plastid_merged_aligned_bp"),
        "plastid_merged_coverage_fraction": scoring_row.get("plastid_merged_coverage_fraction"),
        "mit_support_level": scoring_row.get("mit_support_level"),
        "mit_best_identity": scoring_row.get("mit_best_identity"),
        "mit_merged_aligned_bp": scoring_row.get("mit_merged_aligned_bp"),
        "mit_merged_coverage_fraction": scoring_row.get("mit_merged_coverage_fraction"),
        "graph_score": scoring_row.get("graph_score"),
        "blast_score": scoring_row.get("blast_score"),
        "final_score": scoring_row.get("final_score"),
        "support_signal_count": scoring_row.get("support_signal_count"),
        "strong_signal_count": scoring_row.get("strong_signal_count"),
        "topology_signal_count": scoring_row.get("topology_signal_count"),
        "topology_strong_signal_count": scoring_row.get("topology_strong_signal_count"),
        "confidence_tier": scoring_row.get("confidence_tier"),
        "classification": scoring_row.get("classification"),
        "supporting_evidence": scoring_row.get("supporting_evidence"),
        "penalty_reasons": scoring_row.get("penalty_reasons"),
        "decision_explanation": scoring_row.get("decision_explanation"),
    }


def build_pipeline_rows(
    contigs: dict[str, dict[str, object]],
    graph,
    *,
    mode: str,
    scoring_config,
    blast_support_by_contig: dict[str, dict[str, object]] | None = None,
    nuclear_coverage: float | None = None,
    gc_baseline: float | None = None,
    threads: int = 1,
) -> list[dict[str, object]]:
    feature_rows = calculate_contig_features(contigs)
    feature_by_contig = {str(row["contig_id"]): row for row in feature_rows}
    if graph is None:
        topology_rows = []
    else:
        topology_rows = summarize_contigs(graph)
    topology_by_contig = {str(row["contig_id"]): row for row in topology_rows}
    scoring_rows = score_contigs(
        contigs,
        graph,
        mode=mode,
        config=scoring_config,
        blast_support_by_contig=blast_support_by_contig,
        nuclear_coverage=nuclear_coverage,
        gc_baseline=gc_baseline,
    )
    scoring_by_contig = {str(row["contig_id"]): row for row in scoring_rows}

    contig_ids = sorted(contigs)
    if threads == 1:
        return [
            _build_pipeline_row(
                contig_id,
                contigs,
                feature_by_contig,
                topology_by_contig,
                scoring_by_contig,
            )
            for contig_id in contig_ids
        ]

    with Pool(
        processes=threads,
        initializer=_init_row_build_context,
        initargs=(contigs, feature_by_contig, topology_by_contig, scoring_by_contig),
    ) as pool:
        return pool.map(_build_pipeline_row, contig_ids)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=CliHelpFormatter,
        usage=(
            "%(prog)s [input_gfa] [-o OUTPUT_DIR] "
            "[--mode {graph-only,hybrid,blast-only}] [--plastid-fasta PLASTID_FASTA] "
            "[--mit-fasta MIT_FASTA] [--threads THREADS]"
        ),
        description=(
            "Detect organelle-derived contigs from a GFA assembly.\n\n"
            "graph-only uses graph evidence only and is the default mode.\n"
            "hybrid combines graph evidence with organelle FASTA-based internal BLAST evidence.\n"
            "blast-only uses organelle FASTA-based internal BLAST evidence only.\n"
            "Most users only need the basic input, output, mode, and optional organelle FASTA inputs.\n"
            "Advanced scoring and threshold parameters are optional and usually do not need adjustment."
        ),
        epilog=(
            "Examples:\n"
            "  organelle-cleaner assembly.gfa --output-dir results\n"
            "  organelle-cleaner assembly.gfa --mode hybrid --plastid-fasta plastid.fa --mit-fasta mit.fa --output-dir results\n"
            "  organelle-cleaner assembly.gfa --mode blast-only --plastid-fasta plastid.fa --output-dir results"
        ),
    )

    required_group = parser.add_argument_group("Required Input")
    basic_io_group = parser.add_argument_group("Basic I/O")
    mode_group = parser.add_argument_group("Run Mode")
    blast_group = parser.add_argument_group("BLAST Inputs")
    performance_group = parser.add_argument_group("Performance")
    advanced_group = parser.add_argument_group(
        "Advanced Scoring / Tuning Parameters (Optional Expert Settings)"
    )

    required_group.add_argument(
        "input_gfa",
        nargs="?",
        type=Path,
        help="Input assembly GFA. Required for graph-only and hybrid; normally provided for blast-only as well.",
    )
    basic_io_group.add_argument(
        "--assembly-fasta",
        type=Path,
        default=None,
        help=argparse.SUPPRESS,
    )
    basic_io_group.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        type=Path,
        default=Path("."),
        help="Output directory for reports and FASTA files.",
    )
    basic_io_group.add_argument(
        "--all-candidates-name",
        type=str,
        default=None,
        help=argparse.SUPPRESS,
    )
    performance_group.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of worker processes and the exact thread count used for internal blastn.",
    )
    add_scoring_arguments(
        parser,
        mode_group=mode_group,
        blast_group=blast_group,
        advanced_group=advanced_group,
    )
    return parser


def _validate_cli_inputs(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
) -> None:
    has_organelle_fasta = args.plastid_fasta is not None or args.mit_fasta is not None
    has_contig_source = args.input_gfa is not None or args.assembly_fasta is not None

    if args.mode == "graph-only" and args.input_gfa is None:
        parser.error("graph-only mode requires input_gfa")
    if args.mode == "hybrid" and args.input_gfa is None:
        parser.error("hybrid mode requires input_gfa")
    if args.mode in {"hybrid", "blast-only"} and not has_organelle_fasta:
        parser.error(
            f"{args.mode} mode requires at least one organelle FASTA input via --plastid-fasta or --mit-fasta"
        )
    if args.mode == "blast-only" and not has_contig_source:
        parser.error("blast-only mode requires either input_gfa or --assembly-fasta")


def merge_fasta_sequences(
    contigs: dict[str, dict[str, object]],
    fasta_contigs: dict[str, dict[str, object]],
    *,
    add_missing: bool = False,
) -> None:
    for contig_id, fasta_record in fasta_contigs.items():
        if contig_id not in contigs:
            if add_missing:
                contigs[contig_id] = dict(fasta_record)
            continue
        sequence_value = fasta_record.get("sequence", "")
        sequence = sequence_value if isinstance(sequence_value, str) else ""
        if sequence:
            contigs[contig_id]["sequence"] = sequence
            if not isinstance(contigs[contig_id].get("length"), int):
                contigs[contig_id]["length"] = len(sequence)


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    _validate_cli_inputs(args, parser)
    try:
        scoring_config, blast_config, nuclear_coverage, gc_baseline = config_from_args(args, parser)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        if args.input_gfa is not None:
            contigs, graph = parse_gfa(args.input_gfa)
        else:
            contigs = {}
            graph = None
        if args.assembly_fasta is not None:
            merge_fasta_sequences(
                contigs,
                read_fasta(args.assembly_fasta),
                add_missing=graph is None,
            )
        plastid_blast_tsv = None
        mit_blast_tsv = None
        if args.mode in {"hybrid", "blast-only"} and (
            args.plastid_fasta is not None or args.mit_fasta is not None
        ):
            # Internal BLAST preserves the existing TSV-driven support code path.
            internal_blast_outputs = run_internal_blast(
                assembly_fasta=args.assembly_fasta,
                contigs=contigs,
                plastid_fasta=args.plastid_fasta,
                mit_fasta=args.mit_fasta,
                output_dir=output_dir,
                threads=args.threads,
            )
            plastid_blast_tsv = internal_blast_outputs.plastid_tsv
            mit_blast_tsv = internal_blast_outputs.mit_tsv
        blast_support_by_contig = load_blast_support_by_contig(
            plastid_blast_tsv=plastid_blast_tsv,
            mit_blast_tsv=mit_blast_tsv,
            config=blast_config,
        )
    except FileNotFoundError:
        print("One or more input files were not found", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    try:
        report_rows = build_pipeline_rows(
            contigs,
            graph,
            mode=args.mode,
            scoring_config=scoring_config,
            blast_support_by_contig=blast_support_by_contig,
            nuclear_coverage=nuclear_coverage,
            gc_baseline=gc_baseline,
            threads=args.threads,
        )
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    organelle_ids = [
        str(row["contig_id"])
        for row in report_rows
        if row.get("confidence_tier") == HIGH_CONFIDENCE
    ]
    candidate_rows = [
        row for row in report_rows if row.get("confidence_tier") != "not_flagged"
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
    if args.all_candidates_name:
        write_report(output_dir / args.all_candidates_name, candidate_rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
