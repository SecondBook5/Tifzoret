#!/usr/bin/env python3
"""Triangulate regulatory, co-expression, and STRING association layers."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import networkx as nx


def read(path: str) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def symbols(rows: list[dict[str, str]], columns: tuple[str, ...]) -> set[str]:
    """Return normalized symbols from one or more declared table columns."""
    if rows:
        missing = [column for column in columns if column not in rows[0]]
        if missing:
            raise ValueError(f"Missing required column(s): {', '.join(missing)}")
    return {
        value.strip().upper()
        for row in rows
        for column in columns
        if (value := row.get(column, "")) and value.strip()
    }


def write(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grn-edges", required=True)
    parser.add_argument("--wgcna-hubs", required=True)
    parser.add_argument("--string-up", required=True)
    parser.add_argument("--string-down", required=True)
    parser.add_argument("--contrast-id", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()

def main() -> None:
    args = arguments()
    grn = read(args.grn_edges)
    hubs = read(args.wgcna_hubs)
    string = read(args.string_up) + read(args.string_down)
    grn_genes = symbols(grn, ("source", "target"))
    wgcna = symbols(hubs, ("gene_symbol",))
    string_genes = symbols(string, ("source", "target"))

    rows = []
    evidence_layers = (
        ("regulatory", grn_genes),
        ("coexpression", wgcna),
        ("string_association", string_genes),
    )
    for gene in sorted(grn_genes | wgcna | string_genes):
        layers = [name for name, members in evidence_layers if gene in members]
        rows.append(
            {
                "gene_symbol": gene,
                "layer_count": len(layers),
                "layers": ";".join(layers),
                "triangulated": len(layers) >= 2,
            }
        )

    out = Path(args.outdir)
    fields = ["gene_symbol", "layer_count", "layers", "triangulated"]
    write(out / "tables/multilayer_nodes.tsv", rows, fields)
    triangulated = [row for row in rows if row["triangulated"]]
    write(out / "tables/multilayer_triangulated.tsv", triangulated, fields)

    edge_rows = []
    edge_rows += [
        {"source": row["source"].upper(), "target": row["target"].upper(), "layer": "regulatory"}
        for row in grn
    ]
    edge_rows += [
        {"source": row["source"].upper(), "target": row["target"].upper(), "layer": "string_association"}
        for row in string
    ]
    write(out / "tables/multilayer_edges.tsv", edge_rows, ["source", "target", "layer"])

    graph = nx.Graph()
    graph.add_edges_from(
        (row["source"], row["target"], {"layer": row["layer"]})
        for row in edge_rows
    )
    selected = {row["gene_symbol"] for row in triangulated}
    view = graph.subgraph(selected).copy()
    figure, axis = plt.subplots(figsize=(8, 6))
    axis.axis("off")
    axis.set_title("Multilayer triangulation", loc="left", weight="bold", color="#183B56")
    if view.number_of_nodes():
        layer_counts = {row["gene_symbol"]: row["layer_count"] for row in rows}
        positions = nx.spring_layout(view, seed=1)
        nx.draw_networkx_edges(view, positions, alpha=0.25, edge_color="#6C92AE", ax=axis)
        nx.draw_networkx_nodes(
            view,
            positions,
            node_color=[layer_counts[node] for node in view],
            cmap="viridis",
            node_size=100,
            edgecolors="white",
            ax=axis,
        )
        nx.draw_networkx_labels(view, positions, font_size=6, ax=axis)
    else:
        axis.text(
            0.5,
            0.5,
            "No genes occur in at least two layers",
            ha="center",
            va="center",
            transform=axis.transAxes,
        )

    (out / "figures").mkdir(parents=True, exist_ok=True)
    figure.savefig(out / "figures/multilayer_network.pdf", bbox_inches="tight")
    figure.savefig(out / "figures/multilayer_network.png", dpi=300, bbox_inches="tight")
    plt.close(figure)

    summary = {
        "schema_version": 1,
        "contrast_id": args.contrast_id,
        "genes": len(rows),
        "triangulated": len(triangulated),
        "warnings": [
            "Layers encode distinct evidence types: regulatory, co-expression, and STRING association; overlap does not establish causality."
        ],
    }
    (out / "multilayer_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
