#!/usr/bin/env python3
"""Ollivier-Ricci curvature of the WGCNA co-expression graph (opt-in, exploratory).

For every edge (x, y) of the sparsified co-expression graph that wgcna.R exports,
compute the Ollivier-Ricci curvature

    kappa(x, y) = 1 - W1(m_x, m_y) / d(x, y)

where m_x is a lazy-random-walk probability measure (mass ``alpha`` retained at x,
the rest spread over x's neighbours in proportion to the signed co-expression
adjacency), d(x, y) is the graph hop distance (1 for an edge), and W1 is the exact
1-Wasserstein (earth-mover) distance solved as a small transportation linear
program with scipy.optimize.linprog. Positive curvature marks locally redundant,
robust neighbourhoods (many alternative short paths); negative curvature marks
bottleneck edges bridging otherwise separate regions.

The ground metric is hop distance (the classic graph Ollivier-Ricci convention);
the co-expression weight shapes only the walk's transition probabilities. Outputs:
per-edge curvature, gene-level mean curvature, a per-module robustness scalar (mean
intra-module curvature), a negative-curvature inter-module bridge-gene table, a
curvature histogram, and a curvature-coloured network view. Fully deterministic.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import linprog  # noqa: E402

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from tifzoret.config import load_project  # noqa: E402

# Deterministic layout seed and a cap on how many nodes the network figure draws
# (all TABLES cover the whole graph; only the drawn view is subset to stay legible).
LAYOUT_SEED = 42
DRAW_MAX_NODES = 400
# Diverging colour map for signed curvature (blue = negative/bottleneck, red =
# positive/redundant, near-white midpoint) -- polarity, so a two-hue diverging ramp.
CURVATURE_CMAP = "coolwarm"
NAVY = "#183B56"


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a tab-delimited file into a list of column-keyed dictionaries."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    """Write ``rows`` as a tab-delimited file with a ``fields`` header, creating parent directories."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_graph(edges: list[dict[str, str]]) -> nx.Graph:
    """Assemble the undirected weighted co-expression graph from the edge table."""
    graph = nx.Graph()
    for row in edges:
        weight = float(row["weight"])
        graph.add_edge(row["source"], row["target"], weight=weight)
    return graph


def measure(graph: nx.Graph, node: str, alpha: float) -> dict[str, float]:
    """Lazy-random-walk probability measure at ``node``.

    Mass ``alpha`` stays at ``node``; the remaining ``1 - alpha`` is spread over its
    neighbours in proportion to edge weight (uniformly if the weights sum to zero).
    An isolated node keeps all its mass.
    """
    neighbours = list(graph.neighbors(node))
    if not neighbours:
        return {node: 1.0}
    weights = np.array([graph[node][other]["weight"] for other in neighbours], dtype=float)
    total = float(weights.sum())
    distribution = {node: alpha}
    if total > 0:
        for other, weight in zip(neighbours, weights):
            distribution[other] = (1.0 - alpha) * (weight / total)
    else:
        share = (1.0 - alpha) / len(neighbours)
        for other in neighbours:
            distribution[other] = share
    return distribution


def wasserstein_distance(graph: nx.Graph, source_measure: dict[str, float], target_measure: dict[str, float]) -> float:
    """Exact 1-Wasserstein distance between two measures under graph hop distance.

    Costs are shortest-path hop distances computed within the subgraph induced by
    the two measures' supports (which always contains the connecting edge, so any
    x-neighbour to y-neighbour path is present); pairs unreachable inside that
    local subgraph are capped at 3 hops. The balanced transportation problem is
    solved as a dense linear program (supports are small, ~degree-sized).
    """
    support_x = list(source_measure)
    support_y = list(target_measure)
    induced = graph.subgraph(set(support_x) | set(support_y))
    lengths = {node: nx.single_source_shortest_path_length(induced, node) for node in support_x}
    cost = np.empty((len(support_x), len(support_y)), dtype=float)
    for i, a in enumerate(support_x):
        row = lengths[a]
        for j, b in enumerate(support_y):
            cost[i, j] = 0.0 if a == b else float(row.get(b, 3.0))

    n_x, n_y = len(support_x), len(support_y)
    equality = np.zeros((n_x + n_y, n_x * n_y), dtype=float)
    for i in range(n_x):
        equality[i, i * n_y:(i + 1) * n_y] = 1.0
    for j in range(n_y):
        equality[n_x + j, j::n_y] = 1.0
    demand = np.array([source_measure[a] for a in support_x] + [target_measure[b] for b in support_y], dtype=float)
    result = linprog(cost.reshape(-1), A_eq=equality, b_eq=demand, bounds=(0.0, None), method="highs")
    if not result.success:
        return float("nan")
    return float(result.fun)


def compute_edge_curvature(graph: nx.Graph, alpha: float) -> dict[tuple[str, str], float]:
    """Ollivier-Ricci curvature for every edge; d(x, y) = 1 so kappa = 1 - W1."""
    measures = {node: measure(graph, node, alpha) for node in graph.nodes}
    curvature: dict[tuple[str, str], float] = {}
    for source, target in graph.edges:
        w1 = wasserstein_distance(graph, measures[source], measures[target])
        curvature[(source, target)] = float("nan") if np.isnan(w1) else 1.0 - w1
    return curvature


def draw_network(graph: nx.Graph, edge_curvature: dict[tuple[str, str], float], gene_curvature: dict[str, float], stem: Path) -> int:
    """Render a curvature-coloured network view, subset to the highest-degree nodes.

    Returns the number of nodes actually drawn (the whole graph if within the cap).
    """
    if graph.number_of_nodes() > DRAW_MAX_NODES:
        chosen = sorted(graph.degree, key=lambda item: item[1], reverse=True)[:DRAW_MAX_NODES]
        drawn = graph.subgraph({node for node, _ in chosen})
    else:
        drawn = graph
    if drawn.number_of_nodes() == 0:
        return 0

    layout = nx.spring_layout(drawn, seed=LAYOUT_SEED, weight="weight")
    edge_values = [edge_curvature.get((u, v), edge_curvature.get((v, u), 0.0)) for u, v in drawn.edges]
    node_values = [gene_curvature.get(node, 0.0) for node in drawn.nodes]
    finite = [value for value in edge_values + node_values if not np.isnan(value)]
    span = max(abs(min(finite)), abs(max(finite))) if finite else 1.0
    span = span if span > 0 else 1.0

    figure, axis = plt.subplots(figsize=(7.4, 6.6))
    nx.draw_networkx_edges(drawn, layout, ax=axis, edge_color=edge_values, edge_cmap=plt.get_cmap(CURVATURE_CMAP), edge_vmin=-span, edge_vmax=span, width=0.7, alpha=0.6)
    nodes = nx.draw_networkx_nodes(drawn, layout, ax=axis, node_color=node_values, cmap=plt.get_cmap(CURVATURE_CMAP), vmin=-span, vmax=span, node_size=26, linewidths=0.3, edgecolors="white")
    axis.set_title("Co-expression Ollivier-Ricci curvature", color=NAVY, fontsize=12, fontweight="bold")
    if drawn.number_of_nodes() < graph.number_of_nodes():
        axis.text(0.5, -0.02, f"Top {drawn.number_of_nodes()} of {graph.number_of_nodes()} nodes by degree", transform=axis.transAxes, ha="center", va="top", color="#697783", fontsize=8)
    axis.set_axis_off()
    colourbar = figure.colorbar(nodes, ax=axis, fraction=0.045, pad=0.02)
    colourbar.set_label("Curvature", color=NAVY)
    figure.tight_layout()
    figure.savefig(f"{stem}.pdf", bbox_inches="tight")
    figure.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(figure)
    return drawn.number_of_nodes()


def draw_histogram(values: list[float], stem: Path) -> None:
    """Render the edge-curvature distribution as a histogram with a zero reference."""
    figure, axis = plt.subplots(figsize=(6.2, 4.4))
    finite = [value for value in values if not np.isnan(value)]
    if finite:
        axis.hist(finite, bins=40, color="#4C78A8", edgecolor="white", linewidth=0.4)
        axis.axvline(0.0, color="#B22222", linestyle="--", linewidth=1.0)
    axis.set_title("Distribution of edge Ollivier-Ricci curvature", color=NAVY, fontsize=12, fontweight="bold")
    axis.set_xlabel("Curvature", color=NAVY)
    axis.set_ylabel("Edges", color=NAVY)
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    figure.savefig(f"{stem}.pdf", bbox_inches="tight")
    figure.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    """Compute co-expression curvature tables and figures for one contrast."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--edges", required=True)
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--contrast-id", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    project = load_project(args.project_config)
    settings = (project.config.get("analysis", {}).get("settings", {}) or {}).get("curvature", {}) or {}
    alpha = float(settings.get("alpha", 0.5))
    top_bridges = int(settings.get("top_bridges", 25))

    outdir = Path(args.outdir)
    tables_dir = outdir / "tables"
    figures_dir = outdir / "figures"

    node_rows = read_tsv(Path(args.nodes))
    edge_rows = read_tsv(Path(args.edges))
    symbol_of = {row["gene_id"]: row.get("gene_symbol", row["gene_id"]) for row in node_rows}
    module_of = {row["gene_id"]: row.get("module", "") for row in node_rows}

    graph = build_graph(edge_rows)
    edge_curvature = compute_edge_curvature(graph, alpha) if graph.number_of_edges() else {}

    # Per-edge table (module-tagged, inter-module flag for bridge selection).
    edge_records: list[dict[str, Any]] = []
    gene_accumulator: dict[str, list[float]] = defaultdict(list)
    module_intra: dict[str, list[float]] = defaultdict(list)
    for (source, target), value in edge_curvature.items():
        source_module = module_of.get(source, "")
        target_module = module_of.get(target, "")
        inter_module = source_module != target_module
        edge_records.append({
            "source": source,
            "target": target,
            "source_symbol": symbol_of.get(source, source),
            "target_symbol": symbol_of.get(target, target),
            "source_module": source_module,
            "target_module": target_module,
            "weight": graph[source][target]["weight"],
            "curvature": value,
            "inter_module": inter_module,
        })
        if not np.isnan(value):
            gene_accumulator[source].append(value)
            gene_accumulator[target].append(value)
            if not inter_module:
                module_intra[source_module].append(value)
    edge_records.sort(key=lambda record: (record["curvature"] if not np.isnan(record["curvature"]) else float("inf")))
    write_tsv(
        tables_dir / "edge_curvature.tsv",
        edge_records,
        ["source", "target", "source_symbol", "target_symbol", "source_module", "target_module", "weight", "curvature", "inter_module"],
    )

    # Gene-level mean curvature.
    gene_records = [
        {
            "gene_id": node,
            "gene_symbol": symbol_of.get(node, node),
            "module": module_of.get(node, ""),
            "degree": graph.degree(node) if node in graph else 0,
            "mean_curvature": float(np.mean(gene_accumulator[node])) if gene_accumulator.get(node) else float("nan"),
        }
        for node in graph.nodes
    ]
    gene_records.sort(key=lambda record: record["gene_symbol"])
    gene_curvature = {record["gene_id"]: record["mean_curvature"] for record in gene_records}
    write_tsv(tables_dir / "gene_curvature.tsv", gene_records, ["gene_id", "gene_symbol", "module", "degree", "mean_curvature"])

    # Per-module robustness scalar (mean intra-module edge curvature).
    module_members: dict[str, int] = defaultdict(int)
    for node in graph.nodes:
        module_members[module_of.get(node, "")] += 1
    module_records = [
        {
            "module": module,
            "n_genes": module_members.get(module, 0),
            "n_intra_edges": len(values),
            "mean_intra_curvature": float(np.mean(values)) if values else float("nan"),
        }
        for module, values in sorted(module_intra.items())
    ]
    module_records.sort(key=lambda record: (record["mean_intra_curvature"] if not np.isnan(record["mean_intra_curvature"]) else float("inf")), reverse=True)
    write_tsv(tables_dir / "module_curvature.tsv", module_records, ["module", "n_genes", "n_intra_edges", "mean_intra_curvature"])

    # Negative-curvature inter-module bridges (bottleneck genes linking modules).
    bridge_records = [
        record for record in edge_records
        if record["inter_module"] and not np.isnan(record["curvature"]) and record["curvature"] < 0
    ][:top_bridges]
    write_tsv(
        tables_dir / "curvature_bridges.tsv",
        bridge_records,
        ["source", "target", "source_symbol", "target_symbol", "source_module", "target_module", "weight", "curvature", "inter_module"],
    )

    figures_dir.mkdir(parents=True, exist_ok=True)
    draw_histogram([record["curvature"] for record in edge_records], figures_dir / "curvature_histogram")
    drawn_nodes = draw_network(graph, edge_curvature, gene_curvature, figures_dir / "curvature_network")

    finite_curvature = [record["curvature"] for record in edge_records if not np.isnan(record["curvature"])]
    most_robust = module_records[0]["module"] if module_records else None
    least_robust = module_records[-1]["module"] if module_records else None
    summary = {
        "contrast_id": args.contrast_id,
        "method": "ollivier_ricci",
        "alpha": alpha,
        "ground_metric": "graph hop distance",
        "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "curvature_failures": sum(1 for record in edge_records if np.isnan(record["curvature"])),
        "mean_curvature": float(np.mean(finite_curvature)) if finite_curvature else None,
        "median_curvature": float(np.median(finite_curvature)) if finite_curvature else None,
        "min_curvature": float(np.min(finite_curvature)) if finite_curvature else None,
        "max_curvature": float(np.max(finite_curvature)) if finite_curvature else None,
        "negative_edges": sum(1 for value in finite_curvature if value < 0),
        "modules": len(module_records),
        "most_robust_module": most_robust,
        "least_robust_module": least_robust,
        "bridges_reported": len(bridge_records),
        "drawn_nodes": drawn_nodes,
    }
    (outdir / "curvature_summary.json").parent.mkdir(parents=True, exist_ok=True)
    (outdir / "curvature_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Computed Ollivier-Ricci curvature for {graph.number_of_edges()} edges over {graph.number_of_nodes()} genes")


if __name__ == "__main__":
    main()
