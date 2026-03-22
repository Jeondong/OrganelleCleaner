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


def _component_total_length(
    graph: nx.MultiDiGraph,
    component_nodes: set[str],
) -> int | None:
    total_length = 0
    has_length = False

    for node in component_nodes:
        length = graph.nodes[node].get("length")
        if not isinstance(length, int) or length < 0:
            continue
        total_length += length
        has_length = True

    if not has_length:
        return None

    return total_length


def _component_summaries(graph: nx.MultiDiGraph) -> list[dict[str, object]]:
    summaries: list[dict[str, object]] = []

    for index, component_nodes in enumerate(_sorted_components(graph), start=1):
        component_graph = graph.subgraph(component_nodes)
        component_id = f"component_{index}"
        component_size = len(component_nodes)
        is_circular = component_has_cycle(graph, component_nodes)
        is_isolated = (
            component_size == 1 and graph.degree(next(iter(component_nodes))) == 0
        )
        component_edge_count = component_graph.number_of_edges()
        component_branching_nodes = sum(
            1
            for node in component_nodes
            if component_graph.in_degree(node) > 1 or component_graph.out_degree(node) > 1
        )
        summaries.append(
            {
                "component_id": component_id,
                "component_nodes": component_nodes,
                "component_size": component_size,
                "component_edge_count": component_edge_count,
                "component_branching_nodes": component_branching_nodes,
                "component_total_length": _component_total_length(graph, component_nodes),
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
        component_edge_count = int(summary["component_edge_count"])
        component_branching_nodes = int(summary["component_branching_nodes"])
        is_circular = bool(summary["is_circular"])
        is_isolated = bool(summary["is_isolated"])
        component_total_length = summary["component_total_length"]

        for contig_id in sorted(component_nodes):
            in_degree = int(graph.in_degree(contig_id))
            out_degree = int(graph.out_degree(contig_id))
            self_loop_count = int(graph.number_of_edges(contig_id, contig_id))
            rows.append(
                {
                    "contig_id": contig_id,
                    "component_id": component_id,
                    "component_size": component_size,
                    "component_edge_count": component_edge_count,
                    "component_branching_nodes": component_branching_nodes,
                    "component_total_length": component_total_length,
                    "is_circular": is_circular,
                    "is_isolated": is_isolated,
                    "in_degree": in_degree,
                    "out_degree": out_degree,
                    "total_degree": in_degree + out_degree,
                    "self_loop_count": self_loop_count,
                    "has_self_loop": self_loop_count > 0,
                    "is_branching": in_degree > 1 or out_degree > 1,
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
        ]
    )
    for row in rows:
        writer.writerow(
            [
                row["contig_id"],
                row["component_id"],
                row["component_size"],
                row["component_edge_count"],
                row["component_branching_nodes"],
                row["component_total_length"],
                row["is_circular"],
                row["is_isolated"],
                row["in_degree"],
                row["out_degree"],
                row["total_degree"],
                row["self_loop_count"],
                row["has_self_loop"],
                row["is_branching"],
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
