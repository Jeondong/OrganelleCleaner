#!/usr/bin/env python3
"""Parse a GFA assembly graph and print contigs plus graph edges."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import networkx as nx


def format_parse_error(line_number: int, line: str, message: str) -> ValueError:
    return ValueError(f"Line {line_number}: {message}. Content: {line}")


def parse_tags(fields: list[str], line_number: int, line: str) -> dict[str, object]:
    tags: dict[str, object] = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) != 3:
            continue

        key, value_type, value = parts
        if key in {"LN", "RC"} and value_type != "i":
            raise format_parse_error(
                line_number,
                line,
                f"Tag {key} must use integer type i, found {value_type}:{value}",
            )

        if value_type == "i":
            try:
                tags[key] = int(value)
            except ValueError:
                if key in {"LN", "RC"}:
                    raise format_parse_error(
                        line_number,
                        line,
                        f"Tag {key} has invalid integer value {value!r}",
                    )
                tags[key] = value
        elif value_type == "f":
            try:
                tags[key] = float(value)
            except ValueError:
                tags[key] = value
        else:
            tags[key] = value
    return tags


def parse_gfa(gfa_path: Path) -> tuple[dict[str, dict[str, object]], nx.MultiDiGraph]:
    contigs: dict[str, dict[str, object]] = {}
    graph = nx.MultiDiGraph()

    with gfa_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            record_type = fields[0]

            if record_type == "S":
                if len(fields) < 3:
                    raise format_parse_error(line_number, line, "Invalid S line")

                name = fields[1]
                if name in contigs:
                    raise format_parse_error(
                        line_number,
                        line,
                        f"Duplicate segment ID {name!r}",
                    )

                sequence = fields[2]
                tags = parse_tags(fields[3:], line_number, line)

                if sequence == "*":
                    sequence = ""

                length_tag = tags.get("LN")
                length = length_tag
                if not isinstance(length, int):
                    length = len(sequence)

                read_count = tags.get("RC")
                # Coverage is approximated as RC / length and may differ from true depth.
                if isinstance(length_tag, int) and length_tag > 0 and isinstance(read_count, int):
                    coverage = read_count / length_tag
                else:
                    coverage = None

                contig_data = {
                    "name": name,
                    "sequence": sequence,
                    "length": length,
                    "coverage": coverage,
                }
                contigs[name] = contig_data
                graph.add_node(name, name=name, length=length, coverage=coverage)

            elif record_type == "L":
                if len(fields) < 6:
                    raise format_parse_error(line_number, line, "Invalid L line")

                from_name = fields[1]
                from_orient = fields[2]
                to_name = fields[3]
                to_orient = fields[4]
                overlap = fields[5]
                tags = parse_tags(fields[6:], line_number, line)

                if from_orient not in {"+", "-"}:
                    raise format_parse_error(
                        line_number,
                        line,
                        f"Invalid from orientation {from_orient!r}",
                    )
                if to_orient not in {"+", "-"}:
                    raise format_parse_error(
                        line_number,
                        line,
                        f"Invalid to orientation {to_orient!r}",
                    )
                if from_name not in contigs:
                    raise format_parse_error(
                        line_number,
                        line,
                        f"Link references undefined segment {from_name!r}",
                    )
                if to_name not in contigs:
                    raise format_parse_error(
                        line_number,
                        line,
                        f"Link references undefined segment {to_name!r}",
                    )

                graph.add_edge(
                    from_name,
                    to_name,
                    from_orient=from_orient,
                    to_orient=to_orient,
                    overlap=overlap,
                    **tags,
                )

    return contigs, graph


def write_contig_table(contigs: dict[str, dict[str, object]]) -> None:
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(["name", "length", "coverage"])
    for name in sorted(contigs):
        contig = contigs[name]
        writer.writerow([contig["name"], contig["length"], contig["coverage"]])


def write_graph_edges(graph: nx.MultiDiGraph) -> None:
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(["from", "from_orient", "to", "to_orient", "overlap"])
    edge_rows = sorted(
        graph.edges(data=True),
        key=lambda item: (
            item[0],
            item[2].get("from_orient", ""),
            item[1],
            item[2].get("to_orient", ""),
            item[2].get("overlap", ""),
        ),
    )
    for source, target, edge_data in edge_rows:
        writer.writerow(
            [
                source,
                edge_data.get("from_orient"),
                target,
                edge_data.get("to_orient"),
                edge_data.get("overlap"),
            ]
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Parse a GFA assembly graph into contig and edge tables."
    )
    parser.add_argument("input_gfa", type=Path, help="Path to the input GFA file")
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

    print("Contig table")
    write_contig_table(contigs)
    print()
    print("Graph edges")
    write_graph_edges(graph)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
