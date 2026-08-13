#!/usr/bin/env python3
"""Create auditable rectangular and radial program-aware regulon views."""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import networkx as nx

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from bulk_rna_frame.config import load_project  # noqa: E402
from bulk_rna_frame.figures import normalized_gene_panels  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def program_map(project) -> tuple[dict[str, str], dict[str, str]]:
    mapping: dict[str, str] = {}
    colors: dict[str, str] = {}
    palette = ["#D97706", "#0F9D78", "#C76C9E", "#167BB5", "#7A5195", "#6B8E23", "#8C8C8C"]
    if project.panel_config:
        for panel_index, (panel_id, panel) in enumerate(normalized_gene_panels(project.panel_config).items()):
            for group, genes in panel["groups"].items():
                colors.setdefault(group, panel.get("color") or palette[panel_index % len(palette)])
                for gene in genes:
                    mapping.setdefault(gene.upper(), group)
    colors["Unassigned"] = "#A0A8AE"
    return mapping, colors


def numeric(value: object, default: float) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else default


def program_modularity(graph: nx.Graph, labels: dict[str, str]) -> float:
    groups = defaultdict(set)
    for node in graph:
        groups[labels.get(node, "Unassigned")].add(node)
    communities = [nodes for nodes in groups.values() if nodes]
    return nx.community.modularity(graph, communities, weight="weight") if graph.number_of_edges() and len(communities) > 1 else 0.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--edges", required=True)
    parser.add_argument("--regulators", required=True)
    parser.add_argument("--de", required=True)
    parser.add_argument("--contrast-id", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    project = load_project(args.project_config)
    outdir = Path(args.outdir).resolve(); tables = outdir / "tables"; figures = outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True); figures.mkdir(parents=True, exist_ok=True)
    edges = [row for row in read_tsv(Path(args.edges)) if row.get("measured", "").lower() in {"true", "1", "t"}]
    regulator_rows = read_tsv(Path(args.regulators))
    de_rows = read_tsv(Path(args.de))
    de = {row["gene_symbol"].upper(): row for row in de_rows if row.get("gene_symbol")}
    top_regulators = [row["regulator"] for row in sorted(regulator_rows, key=lambda row: numeric(row.get("adj.P.Val"), math.inf))[:15]]
    regulator_activity = {row["regulator"]: numeric(row.get("logFC"), 0) for row in regulator_rows}
    candidate_edges = [row for row in edges if row["source"] in top_regulators]
    selected_edges: list[dict[str, Any]] = []
    for regulator in top_regulators:
        current = [row for row in candidate_edges if row["source"] == regulator]
        current.sort(key=lambda row: (
            numeric(de.get(row["target"].upper(), {}).get("adjusted_p_value"), math.inf),
            -abs(numeric(de.get(row["target"].upper(), {}).get("log2_fold_change"), 0)),
            row["target"],
        ))
        selected_edges.extend(current[:6])
    mapping, colors = program_map(project)
    graph = nx.DiGraph()
    for row in selected_edges:
        source, target = row["source"], row["target"]
        weight = abs(float(row.get("mor", 1))) * abs(float(row.get("likelihood", 1)))
        graph.add_edge(source, target, weight=max(weight, 0.05), regulation="activating" if float(row.get("mor", 1)) >= 0 else "repressing")
    target_program = {node: mapping.get(node.upper(), "Unassigned") for node in graph if node not in top_regulators}
    labels = dict(target_program)
    for regulator in top_regulators:
        programs = [target_program.get(target, "Unassigned") for target in graph.successors(regulator)] if regulator in graph else []
        labels[regulator] = Counter(programs).most_common(1)[0][0] if programs else "Unassigned"
    undirected = graph.to_undirected()
    communities = nx.community.louvain_communities(undirected, weight="weight", seed=project.config["analysis"].get("random_seed", 1)) if undirected.number_of_edges() else []
    community_map = {node: index + 1 for index, group in enumerate(communities) for node in group}
    node_rows = []
    for node in graph:
        is_regulator = node in top_regulators
        de_row = de.get(node.upper(), {})
        node_rows.append({
            "node": node, "node_type": "regulator" if is_regulator else "target",
            "program": labels.get(node, "Unassigned"), "community": community_map.get(node, 0),
            "log2_fold_change": de_row.get("log2_fold_change", ""),
            "adjusted_p_value": de_row.get("adjusted_p_value", ""),
            "activity_logfc": regulator_activity.get(node, ""), "degree": undirected.degree(node),
        })
    edge_rows = [{"source": source, "target": target, **data} for source, target, data in graph.edges(data=True)]
    write_tsv(tables / "grn_nodes_displayed.tsv", node_rows, ["node", "node_type", "program", "community", "log2_fold_change", "adjusted_p_value", "activity_logfc", "degree"])
    write_tsv(tables / "grn_edges_displayed.tsv", edge_rows, ["source", "target", "weight", "regulation"])

    observed = program_modularity(undirected, labels)
    rng = random.Random(project.config["analysis"].get("random_seed", 1))
    nodes = list(undirected); values = [labels[node] for node in nodes]
    permutations = []
    for _ in range(1000):
        shuffled = values[:]; rng.shuffle(shuffled)
        permutations.append(program_modularity(undirected, dict(zip(nodes, shuffled))))
    p_value = (1 + sum(value >= observed for value in permutations)) / (len(permutations) + 1)
    separation = [{"observed_program_modularity": observed, "permutation_p_value": p_value, "permutations": len(permutations), "communities": len(communities)}]
    write_tsv(tables / "grn_program_separation_test.tsv", separation, list(separation[0]))
    sector_rows = [{"program": program, "color": colors.get(program, colors["Unassigned"]), "nodes": sum(row["program"] == program for row in node_rows), "regulators": sum(row["program"] == program and row["node_type"] == "regulator" for row in node_rows), "targets": sum(row["program"] == program and row["node_type"] == "target" for row in node_rows)} for program in sorted(set(labels.values()))]
    write_tsv(tables / "grn_sector_summary.tsv", sector_rows, ["program", "color", "nodes", "regulators", "targets"])

    # Rectangular breadth-preserving view.
    regulators = [node for node in graph if node in top_regulators]
    targets = [node for node in graph if node not in top_regulators]
    positions = {node: (0, index) for index, node in enumerate(regulators)}
    positions.update({node: (1, index * max(len(regulators), 1) / max(len(targets), 1)) for index, node in enumerate(targets)})
    figure, axis = plt.subplots(figsize=(10.5, max(6.2, 0.22 * max(len(targets), len(regulators)) + 2)))
    nx.draw_networkx_edges(graph, positions, edge_color=["#D66B5D" if graph.edges[edge]["regulation"] == "activating" else "#4E88A8" for edge in graph.edges], width=[0.6 + graph.edges[edge]["weight"] for edge in graph.edges], alpha=0.38, arrows=True, arrowsize=8, ax=axis)
    nx.draw_networkx_nodes(graph, positions, nodelist=regulators, node_shape="s", node_color=[colors.get(labels[node], colors["Unassigned"]) for node in regulators], node_size=140, edgecolors="white", ax=axis)
    nx.draw_networkx_nodes(graph, positions, nodelist=targets, node_shape="o", node_color=[colors.get(labels[node], colors["Unassigned"]) for node in targets], node_size=75, edgecolors="white", ax=axis)
    nx.draw_networkx_labels(graph, positions, font_size=6.2, font_color="#183B56", ax=axis)
    axis.set_title("Regulator–target network", loc="left", weight="bold", color="#183B56"); axis.axis("off")
    figure.savefig(figures / "grn_rectangular.pdf", bbox_inches="tight"); figure.savefig(figures / "grn_rectangular.png", dpi=300, bbox_inches="tight"); plt.close(figure)

    # Compact radial view with program sectors and breadth retained in tables.
    ordered_programs = [row["program"] for row in sector_rows]
    sector_angles = {program: 2 * math.pi * index / max(len(ordered_programs), 1) for index, program in enumerate(ordered_programs)}
    radial_positions = {}
    for program in ordered_programs:
        program_regulators = [node for node in regulators if labels[node] == program]
        program_targets = [node for node in targets if labels[node] == program]
        center = sector_angles[program]
        for index, node in enumerate(program_regulators):
            angle = center + (index - (len(program_regulators) - 1) / 2) * 0.10
            radial_positions[node] = (0.48 * math.cos(angle), 0.48 * math.sin(angle))
        for index, node in enumerate(program_targets):
            angle = center + (index - (len(program_targets) - 1) / 2) * min(0.08, 0.65 / max(len(program_targets), 1))
            radial_positions[node] = (math.cos(angle), math.sin(angle))
    figure, axis = plt.subplots(figsize=(8.4, 8.4))
    nx.draw_networkx_edges(graph, radial_positions, edge_color=["#D66B5D" if graph.edges[edge]["regulation"] == "activating" else "#4E88A8" for edge in graph.edges], width=0.65, alpha=0.25, arrows=True, arrowsize=7, connectionstyle="arc3,rad=0.08", ax=axis)
    nx.draw_networkx_nodes(graph, radial_positions, nodelist=regulators, node_shape="s", node_color=[colors.get(labels[node], colors["Unassigned"]) for node in regulators], node_size=130, edgecolors="white", ax=axis)
    nx.draw_networkx_nodes(graph, radial_positions, nodelist=targets, node_shape="o", node_color=[colors.get(labels[node], colors["Unassigned"]) for node in targets], node_size=65, edgecolors="white", ax=axis)
    label_nodes = regulators + sorted(targets, key=lambda node: -undirected.degree(node))[:25]
    nx.draw_networkx_labels(graph, radial_positions, labels={node: node for node in label_nodes}, font_size=6, font_color="#183B56", ax=axis)
    for program, angle in sector_angles.items():
        axis.text(1.22 * math.cos(angle), 1.22 * math.sin(angle), program.replace("_", " "), ha="center", va="center", color=colors.get(program), weight="bold", fontsize=8)
    axis.set_title("Program-aware radial regulator network", loc="left", weight="bold", color="#183B56"); axis.axis("off"); axis.set_aspect("equal")
    figure.savefig(figures / "grn_radial.pdf", bbox_inches="tight"); figure.savefig(figures / "grn_radial.png", dpi=300, bbox_inches="tight"); plt.close(figure)
    (outdir / "grn_summary.json").write_text(json.dumps({
        "schema_version": 1, "contrast_id": args.contrast_id, "displayed_regulators": len(regulators),
        "displayed_targets": len(targets), "displayed_edges": len(edge_rows),
        "selection": "top 15 differential regulators and top 6 DE-ranked measured targets per regulator",
        "full_edges": str(Path(args.edges).resolve()), "program_modularity": observed,
        "program_separation_permutation_p": p_value,
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
