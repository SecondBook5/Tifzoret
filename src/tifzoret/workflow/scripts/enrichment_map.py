#!/usr/bin/env python3
"""Enrichment-similarity map: cluster enriched terms by the genes they share.

A per-contrast pathway analysis emits dozens of significant terms, many of which
are redundant -- they are enriched because they share large parts of their gene
membership (e.g. several overlapping inflammation gene sets driven by the same
leading-edge genes). An enrichment map answers "which enriched terms are really
the same biological program?" by drawing a network whose nodes are enriched terms
and whose edges connect terms that overlap in the genes that drove their
enrichment. Community detection then groups the terms into a handful of
non-redundant themes.

Node set: every term called significant by the engine's own pathway analysis --
significant in fgsea (``padj`` below the DE FDR) or in ORA (``adjusted_p_value``
below the DE FDR). No new statistics are computed here. Each term's *signature*
is the union of the genes that drove its enrichment in this contrast: the fgsea
leading-edge genes and/or the ORA overlap genes. Edge weight is the Jaccard
similarity of two terms' signatures; an edge is kept when the similarity is at
least ``min_similarity``.

Communities are found with NetworkX greedy modularity maximisation, which is
deterministic (no random initialisation); nodes and edges are added in sorted
order and communities are relabelled by descending size then first member so the
cluster ids are stable across runs.

Outputs (all under the contrast's ``analyses/enrichment_map/``):

* ``tables/enrichment_map_nodes.tsv`` -- one row per enriched term with its
  signature size, significance, direction, degree, and cluster;
* ``tables/enrichment_map_edges.tsv`` -- term-term Jaccard edges and shared genes;
* ``tables/enrichment_map_clusters.tsv`` -- one row per community (theme);
* ``figures/enrichment_map.{pdf,png}`` -- the term-similarity network;
* ``enrichment_map_summary.json``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from itertools import combinations
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from tifzoret.config import load_project  # noqa: E402

# House palette (matches theme_publication / curvature.py / consensus.py).
NAVY = "#183B56"
MID_GREY = "#697783"
UP_COLOUR = "#C0392B"
DOWN_COLOUR = "#2C6FBB"
# Okabe-Ito qualitative palette for community colours: colourblind-safe by
# construction, so cluster identity survives CVD without secondary encoding.
OKABE_ITO = [
    "#0072B2", "#E69F00", "#009E73", "#CC79A7",
    "#56B4E9", "#D55E00", "#F0E442", "#000000",
]
OTHER_COLOUR = "#B7BEC4"

DEFAULT_MIN_SIMILARITY = 0.25
DEFAULT_TOP_TERMS = 60
DEFAULT_FDR = 0.05
# Cap on how many term labels are drawn (all TABLES are complete; only the drawn
# labels are thinned so the network stays legible).
MAX_DRAWN_LABELS = 24
TINY = 1e-300


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a tab-delimited file into a list of column-keyed dictionaries."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    """Write ``rows`` as a tab-delimited file with a ``fields`` header, creating parents."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def _float(value: str, default: float = float("nan")) -> float:
    """Parse a float, mapping blanks/NA to ``default``."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _genes(value: str) -> set[str]:
    """Split a ``;``-joined gene list into a set, dropping blanks."""
    return {token.strip() for token in (value or "").split(";") if token.strip()}


def collect_terms(
    fgsea_rows: list[dict[str, str]], ora_rows: list[dict[str, str]], fdr: float
) -> dict[str, dict[str, Any]]:
    """Build the node table keyed by term id from significant fgsea and ORA rows.

    A term's signature is the union of its fgsea leading-edge genes and its ORA
    overlap genes; its direction comes from the fgsea NES sign when available,
    otherwise from the more-significant ORA direction.
    """
    terms: dict[str, dict[str, Any]] = {}

    for row in fgsea_rows:
        term = (row.get("pathway") or "").strip()
        if not term:
            continue
        padj = _float(row.get("padj", ""))
        if math.isnan(padj) or padj >= fdr:
            continue
        record = terms.setdefault(term, _blank_term(term, row))
        record["signature"] |= _genes(row.get("leadingEdge", ""))
        record["fgsea_padj"] = padj
        record["nes"] = _float(row.get("NES", ""))
        record["direction"] = (row.get("direction") or "").strip() or record["direction"]

    for row in ora_rows:
        term = (row.get("pathway") or "").strip()
        if not term:
            continue
        adjusted = _float(row.get("adjusted_p_value", ""))
        count = _float(row.get("count", ""), 0.0)
        if math.isnan(adjusted) or adjusted >= fdr or count <= 0:
            continue
        record = terms.setdefault(term, _blank_term(term, row))
        record["signature"] |= _genes(row.get("overlap_genes", ""))
        # Keep the strongest ORA adjusted p across the up/down rows.
        if math.isnan(record["ora_padj"]) or adjusted < record["ora_padj"]:
            record["ora_padj"] = adjusted
            if not record["direction"]:
                record["direction"] = (row.get("direction") or "").strip()

    # A term with no signature genes cannot contribute overlap edges; drop it.
    return {term: record for term, record in terms.items() if record["signature"]}


def _blank_term(term: str, row: dict[str, str]) -> dict[str, Any]:
    """Seed a node record with label/provider carried from the first row that names it."""
    return {
        "term": term,
        "label": (row.get("pathway_label") or term).strip() or term,
        "provider": (row.get("gene_set_source") or "").strip(),
        "direction": "",
        "nes": float("nan"),
        "fgsea_padj": float("nan"),
        "ora_padj": float("nan"),
        "signature": set(),
    }


def best_p(record: dict[str, Any]) -> float:
    """Return the strongest (smallest) adjusted p across fgsea and ORA for ranking."""
    candidates = [p for p in (record["fgsea_padj"], record["ora_padj"]) if not math.isnan(p)]
    return min(candidates) if candidates else 1.0


def build_graph(terms: dict[str, dict[str, Any]], min_similarity: float) -> tuple[nx.Graph, list[dict[str, Any]]]:
    """Build the term-similarity graph and the edge table (Jaccard >= threshold)."""
    graph = nx.Graph()
    for term in sorted(terms):
        graph.add_node(term)
    edge_rows: list[dict[str, Any]] = []
    for a, b in combinations(sorted(terms), 2):
        sig_a, sig_b = terms[a]["signature"], terms[b]["signature"]
        shared = sig_a & sig_b
        if not shared:
            continue
        union = sig_a | sig_b
        jaccard = len(shared) / len(union)
        if jaccard < min_similarity:
            continue
        graph.add_edge(a, b, weight=jaccard)
        edge_rows.append({
            "term_a": a,
            "term_b": b,
            "label_a": terms[a]["label"],
            "label_b": terms[b]["label"],
            "jaccard": jaccard,
            "n_shared": len(shared),
            "shared_genes": ";".join(sorted(shared)),
        })
    edge_rows.sort(key=lambda r: (-r["jaccard"], r["term_a"], r["term_b"]))
    return graph, edge_rows


def detect_communities(graph: nx.Graph, terms: dict[str, dict[str, Any]]) -> dict[str, int]:
    """Assign each node a stable cluster id via deterministic greedy modularity.

    Communities are relabelled by descending size, then by the strongest
    significance among members, then by first term id -- a total order, so the
    ids do not depend on NetworkX's internal community ordering.
    """
    if graph.number_of_nodes() == 0:
        return {}
    communities = nx.community.greedy_modularity_communities(graph, weight="weight")
    ordered = sorted(
        (sorted(group) for group in communities),
        key=lambda members: (-len(members), best_p(terms[members[0]]), members[0]),
    )
    return {term: index for index, members in enumerate(ordered) for term in members}


def draw_map(
    graph: nx.Graph,
    terms: dict[str, dict[str, Any]],
    clusters: dict[str, int],
    n_clusters: int,
    seed: int,
    stem: Path,
) -> int:
    """Render the enrichment map: term nodes coloured by cluster, sized by significance."""
    stem.parent.mkdir(parents=True, exist_ok=True)
    if graph.number_of_nodes() == 0:
        figure, axis = plt.subplots(figsize=(7.0, 5.0))
        axis.text(0.5, 0.5, "No enriched terms above the FDR threshold", ha="center", va="center", color=MID_GREY)
        axis.set_axis_off()
        figure.savefig(f"{stem}.pdf", bbox_inches="tight")
        figure.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
        plt.close(figure)
        return 0

    positions = nx.spring_layout(
        graph, seed=seed, weight="weight",
        k=1.6 / math.sqrt(max(graph.number_of_nodes(), 1)),
    )
    nodes = sorted(graph.nodes())
    cluster_colour = {
        index: (OKABE_ITO[index] if index < len(OKABE_ITO) else OTHER_COLOUR)
        for index in range(n_clusters)
    }
    node_colours = [cluster_colour.get(clusters.get(node, 0), OTHER_COLOUR) for node in nodes]
    ring_colours = [
        UP_COLOUR if terms[node]["direction"] == "up_in_numerator"
        else DOWN_COLOUR if terms[node]["direction"] == "down_in_numerator"
        else MID_GREY
        for node in nodes
    ]
    sizes = [90 + 260 * min(1.0, -math.log10(max(best_p(terms[node]), TINY)) / 6.0) for node in nodes]

    figure, axis = plt.subplots(figsize=(8.4, 6.6))
    if graph.number_of_edges():
        nx.draw_networkx_edges(
            graph, positions, ax=axis, alpha=0.35, edge_color="#9AA6AE",
            width=[0.5 + 3.0 * graph.edges[edge]["weight"] for edge in graph.edges],
        )
    nx.draw_networkx_nodes(
        graph, positions, nodelist=nodes, ax=axis, node_color=node_colours,
        node_size=sizes, edgecolors=ring_colours, linewidths=1.4,
    )
    # Label only the most significant terms to keep the map legible.
    labelled = sorted(nodes, key=lambda node: (best_p(terms[node]), node))[:MAX_DRAWN_LABELS]
    labels = {node: terms[node]["label"] for node in labelled}
    nx.draw_networkx_labels(graph, positions, labels=labels, ax=axis, font_size=6, font_color=NAVY)

    axis.set_title("Enrichment map: shared-gene term similarity", color=NAVY, fontsize=13, fontweight="bold")
    axis.set_axis_off()
    cluster_handles = [
        plt.Line2D([0], [0], marker="o", linestyle="none", markersize=9,
                   markerfacecolor=cluster_colour[index], markeredgecolor="none",
                   label=f"Cluster {index + 1}")
        for index in range(min(n_clusters, len(OKABE_ITO)))
    ]
    direction_handles = [
        plt.Line2D([0], [0], marker="o", linestyle="none", markersize=9,
                   markerfacecolor="white", markeredgecolor=UP_COLOUR, markeredgewidth=1.6, label="Up in numerator"),
        plt.Line2D([0], [0], marker="o", linestyle="none", markersize=9,
                   markerfacecolor="white", markeredgecolor=DOWN_COLOUR, markeredgewidth=1.6, label="Down in numerator"),
    ]
    axis.legend(
        handles=cluster_handles + direction_handles,
        loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7, frameon=False,
    )
    figure.savefig(f"{stem}.pdf", bbox_inches="tight")
    figure.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(figure)
    return len(nodes)


def main() -> None:
    """Build the enrichment-similarity map for one contrast."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--fgsea", required=True, help="pathways fgsea.tsv")
    parser.add_argument("--ora", required=True, help="pathways ora.tsv")
    parser.add_argument("--contrast-id", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    project = load_project(args.project_config)
    analysis = project.config.get("analysis", {}) or {}
    settings = (analysis.get("settings", {}) or {}).get("enrichment_map", {}) or {}
    min_similarity = float(settings.get("min_similarity", DEFAULT_MIN_SIMILARITY))
    top_terms = int(settings.get("top_terms", DEFAULT_TOP_TERMS))
    seed = int(settings.get("seed", analysis.get("random_seed", 1)))
    fdr = float((project.config.get("figures", {}) or {}).get("de", {}).get("fdr", DEFAULT_FDR))

    outdir = Path(args.outdir)
    tables_dir = outdir / "tables"
    figures_dir = outdir / "figures"

    fgsea_rows = read_tsv(Path(args.fgsea))
    ora_rows = read_tsv(Path(args.ora))

    terms = collect_terms(fgsea_rows, ora_rows, fdr)
    # Cap to the most significant terms so the map stays interpretable; all
    # dropped terms are simply absent (there is no separate complete node table).
    kept = sorted(terms, key=lambda term: (best_p(terms[term]), term))[:top_terms]
    terms = {term: terms[term] for term in kept}

    graph, edge_rows = build_graph(terms, min_similarity)
    clusters = detect_communities(graph, terms)
    n_clusters = (max(clusters.values()) + 1) if clusters else 0

    node_fields = [
        "term", "label", "provider", "direction", "nes",
        "fgsea_padj", "ora_padj", "signature_size", "degree", "cluster", "cluster_label",
    ]
    node_rows = []
    for term in sorted(terms):
        record = terms[term]
        cluster = clusters.get(term, 0)
        node_rows.append({
            "term": term,
            "label": record["label"],
            "provider": record["provider"],
            "direction": record["direction"],
            "nes": record["nes"],
            "fgsea_padj": record["fgsea_padj"],
            "ora_padj": record["ora_padj"],
            "signature_size": len(record["signature"]),
            "degree": graph.degree(term) if term in graph else 0,
            "cluster": cluster,
            "cluster_label": f"C{cluster + 1}",
        })
    node_rows.sort(key=lambda r: (r["cluster"], best_p(terms[r["term"]]), r["term"]))
    write_tsv(tables_dir / "enrichment_map_nodes.tsv", node_rows, node_fields)

    write_tsv(
        tables_dir / "enrichment_map_edges.tsv",
        edge_rows,
        ["term_a", "term_b", "label_a", "label_b", "jaccard", "n_shared", "shared_genes"],
    )

    # One row per community (theme).
    cluster_rows: list[dict[str, Any]] = []
    for index in range(n_clusters):
        members = sorted((term for term, c in clusters.items() if c == index), key=lambda t: (best_p(terms[t]), t))
        n_up = sum(1 for term in members if terms[term]["direction"] == "up_in_numerator")
        n_down = sum(1 for term in members if terms[term]["direction"] == "down_in_numerator")
        representative = members[0]
        cluster_rows.append({
            "cluster": index,
            "cluster_label": f"C{index + 1}",
            "size": len(members),
            "representative_term": representative,
            "representative_label": terms[representative]["label"],
            "n_up": n_up,
            "n_down": n_down,
            "member_terms": ";".join(members),
        })
    write_tsv(
        tables_dir / "enrichment_map_clusters.tsv",
        cluster_rows,
        ["cluster", "cluster_label", "size", "representative_term", "representative_label", "n_up", "n_down", "member_terms"],
    )

    drawn = draw_map(graph, terms, clusters, n_clusters, seed, figures_dir / "enrichment_map")

    summary = {
        "schema_version": 1,
        "method": "shared-gene (Jaccard) enrichment map with greedy-modularity clustering",
        "contrast_id": args.contrast_id,
        "fdr": fdr,
        "min_similarity": min_similarity,
        "top_terms": top_terms,
        "terms": len(terms),
        "edges": graph.number_of_edges(),
        "clusters": n_clusters,
        "drawn_terms": drawn,
    }
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "enrichment_map_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Built enrichment map for {args.contrast_id}: {len(terms)} terms, {graph.number_of_edges()} edges, {n_clusters} clusters")


if __name__ == "__main__":
    main()
