#!/usr/bin/env python3
"""Analyze connected components and circularity in a GFA-derived graph."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import networkx as nx

from parse_gfa import parse_gfa


def _sorted_components(graph: nx.MultiDiGraph) -> list[set[str]]:
    undirected_graph = graph.to_undirected(as_view=True)
    components = [set(component) for component in nx.connected_components(undirected_graph)]
    components.sort(key=lambda component: min(str(node) for node in component))
    return components


def connected_components(graph: nx.MultiDiGraph) -> list[set[str]]:
    """Return connected components using the undirected view of the graph."""
    return _sorted_components(graph)


def component_has_cycle(graph: nx.MultiDiGraph, component_nodes: set[str]) -> bool:
    """Return True if the component forms a simple directed cycle (no branching)."""
    component_graph = graph.subgraph(component_nodes)

    if not component_nodes:
        return False

    if not nx.is_strongly_connected(component_graph):
        return False

    for node in component_nodes:
        if component_graph.in_degree(node) != 1 or component_graph.out_degree(node) != 1:
            return False

    return True


def _component_summaries(graph: nx.MultiDiGraph) -> list[dict[str, object]]:
    summaries: list[dict[str, object]] = []

    for index, component_nodes in enumerate(_sorted_components(graph), start=1):
        component_id = f"component_{index}"
        component_size = len(component_nodes)
        is_circular = component_has_cycle(graph, component_nodes)
        is_isolated = (
            component_size == 1 and graph.degree(next(iter(component_nodes))) == 0
        )
        summaries.append(
            {
                "component_id": component_id,
                "component_nodes": component_nodes,
                "component_size": component_size,
                "is_circular": is_circular,
                "is_isolated": is_isolated,
            }
        )

    return summaries


def circular_components(graph: nx.MultiDiGraph) -> dict[str, bool]:
    """Map component IDs to whether the component contains a cycle."""
    return {
        str(summary["component_id"]): bool(summary["is_circular"])
        for summary in _component_summaries(graph)
    }


def isolated_components(graph: nx.MultiDiGraph) -> dict[str, bool]:
    """Map component IDs to whether the component is isolated."""
    return {
        str(summary["component_id"]): bool(summary["is_isolated"])
        for summary in _component_summaries(graph)
    }


def summarize_contigs(graph: nx.MultiDiGraph) -> list[dict[str, object]]:
    """Build per-contig annotations for component membership and graph topology."""
    rows: list[dict[str, object]] = []

    for summary in _component_summaries(graph):
        component_id = str(summary["component_id"])
        component_nodes = summary["component_nodes"]
        component_size = int(summary["component_size"])
        is_circular = bool(summary["is_circular"])
        is_isolated = bool(summary["is_isolated"])

        for contig_id in sorted(component_nodes):
            rows.append(
                {
                    "contig_id": contig_id,
                    "component_id": component_id,
                    "component_size": component_size,
                    "is_circular": is_circular,
                    "is_isolated": is_isolated,
                }
            )

    rows.sort(key=lambda row: str(row["contig_id"]))
    return rows


def write_component_table(rows: list[dict[str, object]]) -> None:
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "contig_id",
            "component_id",
            "component_size",
            "is_circular",
            "is_isolated",
        ]
    )
    for row in rows:
        writer.writerow(
            [
                row["contig_id"],
                row["component_id"],
                row["component_size"],
                row["is_circular"],
                row["is_isolated"],
            ]
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Report connected components and circularity for a GFA graph."
    )
    parser.add_argument("input_gfa", type=Path, help="Path to the input GFA file")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        _, graph = parse_gfa(args.input_gfa)
        rows = summarize_contigs(graph)
    except FileNotFoundError:
        print(f"Input file not found: {args.input_gfa}", file=sys.stderr)
        return 1
    except (ImportError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1

    write_component_table(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
