#!/usr/bin/env python3
"""Fetch, audit, analyze, and render contrast-specific STRING networks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import networkx as nx

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from bulk_rna_frame.config import load_project  # noqa: E402


API = "https://string-db.org/api/tsv"
STRING_NETWORK_BATCH_SIZE = 900


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def cached_post(
    endpoint: str,
    parameters: dict[str, object],
    cache_dir: Path,
    offline: bool,
    refresh: bool,
) -> str:
    encoded = urllib.parse.urlencode(parameters).encode()
    key = hashlib.sha256(endpoint.encode() + b"\0" + encoded).hexdigest()
    response_path = cache_dir / f"string_{endpoint}_{key}.tsv"
    receipt_path = cache_dir / f"string_{endpoint}_{key}.json"
    if response_path.is_file() and receipt_path.is_file() and not refresh:
        return response_path.read_text(encoding="utf-8")
    if offline:
        raise RuntimeError(f"offline mode: missing STRING cache entry {key}")
    request = urllib.request.Request(
        f"{API}/{endpoint}", data=encoded, headers={"User-Agent": "BulkRNAFrame/0.1"}, method="POST"
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        payload = response.read()
        release = response.headers.get("X-STRING-Version") or response.headers.get("Last-Modified")
    cache_dir.mkdir(parents=True, exist_ok=True)
    response_path.write_bytes(payload)
    receipt = {
        "schema_version": 1,
        "provider": "STRING",
        "endpoint": endpoint,
        "database_release": release or "reported by live API response",
        "retrieval_time_utc": datetime.now(timezone.utc).isoformat(),
        "request_parameters": {key: value for key, value in parameters.items() if key != "identifiers"},
        "requested_identifier_count": len(str(parameters.get("identifiers", "")).split("\r")),
        "license_notice": "STRING data are subject to STRING licensing and attribution terms.",
        "sha256": sha256_bytes(payload),
    }
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return payload.decode("utf-8")


def parse_response(text: str) -> list[dict[str, str]]:
    if not text.strip():
        return []
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def fetch_induced_network(
    identifiers: list[str],
    taxonomy: int,
    required_score: int,
    cache_dir: Path,
    offline: bool,
    refresh: bool,
    batch_size: int = STRING_NETWORK_BATCH_SIZE,
) -> tuple[list[dict[str, str]], int]:
    """Fetch every within- and between-batch edge for an arbitrarily large gene set.

    STRING rejects very large induced-network requests. Pairing deterministic
    chunks covers every possible identifier pair without selecting genes out of
    the submitted set. Repeated within-chunk edges are deduplicated afterward.
    """
    unique = list(dict.fromkeys(identifier for identifier in identifiers if identifier))
    if not unique:
        return [], 0
    chunks = [unique[index : index + batch_size] for index in range(0, len(unique), batch_size)]
    retained: dict[tuple[str, str], dict[str, str]] = {}
    calls = 0
    for left_index, left in enumerate(chunks):
        for right_index in range(left_index, len(chunks)):
            submitted = left if right_index == left_index else left + chunks[right_index]
            response = cached_post(
                "network",
                {
                    "identifiers": "\r".join(submitted),
                    "species": taxonomy,
                    "required_score": required_score,
                    "network_type": "functional",
                },
                cache_dir,
                offline,
                refresh,
            )
            calls += 1
            for row in parse_response(response):
                source = row.get("stringId_A") or row.get("preferredName_A", "")
                target = row.get("stringId_B") or row.get("preferredName_B", "")
                if not source or not target or source == target:
                    continue
                key = tuple(sorted((source, target)))
                previous = retained.get(key)
                if previous is None or float(row.get("score", 0)) > float(previous.get("score", 0)):
                    retained[key] = row
    return list(retained.values()), calls


def empty_plot(path: Path, title: str, message: str) -> None:
    figure, axis = plt.subplots(figsize=(7, 5))
    axis.axis("off")
    axis.set_title(title, loc="left", weight="bold", color="#183B56")
    axis.text(0.5, 0.5, message, ha="center", va="center", color="#697783", transform=axis.transAxes)
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def repelled_label_positions(
    positions: dict[str, tuple[float, float]], labels: list[str], iterations: int = 240
) -> dict[str, tuple[float, float]]:
    """Deterministically separate a small set of network labels in data space."""
    result = {node: [float(positions[node][0]), float(positions[node][1])] for node in labels}
    if len(result) < 2:
        return {node: tuple(value) for node, value in result.items()}
    xs = [float(value[0]) for value in positions.values()]
    ys = [float(value[1]) for value in positions.values()]
    minimum_x = max(max(xs) - min(xs), 1.0) * 0.075
    minimum_y = max(max(ys) - min(ys), 1.0) * 0.050
    for _ in range(iterations):
        movement = {node: [0.0, 0.0] for node in labels}
        for index, left in enumerate(labels):
            for right in labels[index + 1 :]:
                dx = result[left][0] - result[right][0]
                dy = result[left][1] - result[right][1]
                if abs(dx) >= minimum_x or abs(dy) >= minimum_y:
                    continue
                if abs(dx) < 1e-9 and abs(dy) < 1e-9:
                    dx = 1.0 if left < right else -1.0
                    dy = 0.5
                push_x = 0.08 * (minimum_x - abs(dx)) * (1 if dx >= 0 else -1)
                push_y = 0.08 * (minimum_y - abs(dy)) * (1 if dy >= 0 else -1)
                movement[left][0] += push_x
                movement[left][1] += push_y
                movement[right][0] -= push_x
                movement[right][1] -= push_y
        for node in labels:
            movement[node][0] += 0.008 * (float(positions[node][0]) - result[node][0])
            movement[node][1] += 0.008 * (float(positions[node][1]) - result[node][1])
            result[node][0] += movement[node][0]
            result[node][1] += movement[node][1]
    return {node: tuple(value) for node, value in result.items()}


def network_for_direction(
    direction: str,
    genes: list[dict[str, str]],
    taxonomy: int,
    required_score: int,
    max_display: int,
    seed: int,
    cache_dir: Path,
    offline: bool,
    refresh: bool,
    tables: Path,
    figures: Path,
) -> dict[str, Any]:
    prefix = "up" if direction == "up_in_numerator" else "down"
    symbols = list(dict.fromkeys(row["gene_symbol"] for row in genes if row.get("gene_symbol")))
    input_rows = [{"gene_symbol": symbol, "direction": direction, "log2_fold_change": next(row["log2_fold_change"] for row in genes if row["gene_symbol"] == symbol)} for symbol in symbols]
    mapping = parse_response(cached_post(
        "get_string_ids",
        {"identifiers": "\r".join(symbols), "species": taxonomy, "limit": 1, "echo_query": 1},
        cache_dir, offline, refresh,
    )) if symbols else []
    mapping_by_query = {row.get("queryItem", ""): row for row in mapping}
    for row in input_rows:
        mapped = mapping_by_query.get(row["gene_symbol"])
        row["mapped"] = bool(mapped)
        row["string_id"] = mapped.get("stringId", "") if mapped else ""
        row["preferred_name"] = mapped.get("preferredName", "") if mapped else ""

    mapped_identifiers = [str(row["string_id"]) for row in input_rows if row["mapped"]]
    raw_edges, network_api_calls = fetch_induced_network(
        mapped_identifiers,
        taxonomy,
        required_score,
        cache_dir,
        offline,
        refresh,
    )
    graph = nx.Graph()
    effect: dict[str, float] = {}
    for row in input_rows:
        value = float(row["log2_fold_change"])
        effect[row["gene_symbol"]] = value
        if row["preferred_name"]:
            effect[row["preferred_name"]] = value
    edge_rows: list[dict[str, Any]] = []
    for row in raw_edges:
        source, target = row.get("preferredName_A", ""), row.get("preferredName_B", "")
        if not source or not target or source == target:
            continue
        score = float(row.get("score", 0))
        graph.add_edge(source, target, score=score)
        edge_rows.append({
            "source": source, "target": target, "combined_score": score,
            "neighborhood_score": row.get("nscore", ""), "fusion_score": row.get("fscore", ""),
            "cooccurrence_score": row.get("pscore", ""), "coexpression_score": row.get("ascore", ""),
            "experimental_score": row.get("escore", ""), "database_score": row.get("dscore", ""),
            "textmining_score": row.get("tscore", ""),
        })
    communities: dict[str, int] = {}
    if graph.number_of_edges():
        groups = nx.community.louvain_communities(graph, weight="score", seed=seed)
        communities = {node: index + 1 for index, group in enumerate(groups) for node in group}
    degree = nx.degree_centrality(graph) if graph else {}
    betweenness = nx.betweenness_centrality(graph, weight=None) if graph else {}
    node_rows = [{
        "gene_symbol": node,
        "log2_fold_change": effect.get(node, 0.0),
        "degree": graph.degree(node),
        "degree_centrality": degree.get(node, 0.0),
        "betweenness_centrality": betweenness.get(node, 0.0),
        "community": communities.get(node, 0),
    } for node in sorted(graph)]
    connected = set(graph)
    for row in input_rows:
        row["connected"] = row["preferred_name"] in connected or row["gene_symbol"] in connected
    unmapped = [row for row in input_rows if not row["mapped"]]
    unconnected = [row for row in input_rows if row["mapped"] and not row["connected"]]

    selected = sorted(graph, key=lambda node: (-graph.degree(node), node))[:max_display]
    display_graph = graph.subgraph(selected).copy()
    display_nodes = [row for row in node_rows if row["gene_symbol"] in display_graph]
    display_edges = [row for row in edge_rows if row["source"] in display_graph and row["target"] in display_graph]
    write_tsv(tables / f"string_{prefix}_input_genes.tsv", input_rows, ["gene_symbol", "direction", "log2_fold_change", "mapped", "string_id", "preferred_name", "connected"])
    write_tsv(tables / f"string_{prefix}_unmapped_genes.tsv", unmapped, ["gene_symbol", "direction", "log2_fold_change", "mapped", "string_id", "preferred_name", "connected"])
    write_tsv(tables / f"string_{prefix}_unconnected_genes.tsv", unconnected, ["gene_symbol", "direction", "log2_fold_change", "mapped", "string_id", "preferred_name", "connected"])
    write_tsv(tables / f"string_{prefix}_nodes.tsv", node_rows, ["gene_symbol", "log2_fold_change", "degree", "degree_centrality", "betweenness_centrality", "community"])
    write_tsv(tables / f"string_{prefix}_edges.tsv", edge_rows, ["source", "target", "combined_score", "neighborhood_score", "fusion_score", "cooccurrence_score", "coexpression_score", "experimental_score", "database_score", "textmining_score"])
    write_tsv(tables / f"string_{prefix}_nodes_displayed.tsv", display_nodes, ["gene_symbol", "log2_fold_change", "degree", "degree_centrality", "betweenness_centrality", "community"])
    write_tsv(tables / f"string_{prefix}_edges_displayed.tsv", display_edges, ["source", "target", "combined_score", "neighborhood_score", "fusion_score", "cooccurrence_score", "coexpression_score", "experimental_score", "database_score", "textmining_score"])

    stem = figures / f"string_{prefix}_network"
    if not display_graph.number_of_edges():
        empty_plot(stem, f"STRING {prefix}regulated network", "No connected STRING subnetwork")
    else:
        positions = nx.spring_layout(display_graph, seed=seed, weight="score", k=1.4 / math.sqrt(max(len(display_graph), 1)))
        values = [effect.get(node, 0.0) for node in display_graph]
        limit = max(max(abs(value) for value in values), 1e-6)
        figure, axis = plt.subplots(figsize=(8.2, 7.0))
        nx.draw_networkx_edges(display_graph, positions, width=[0.4 + 2.6 * display_graph.edges[edge]["score"] for edge in display_graph.edges], alpha=0.32, edge_color="#6C92AE", ax=axis)
        nodes = nx.draw_networkx_nodes(display_graph, positions, node_color=values, cmap="coolwarm", vmin=-limit, vmax=limit, node_size=[80 + 24 * display_graph.degree(node) for node in display_graph], edgecolors="white", linewidths=0.6, ax=axis)
        labels = sorted(display_graph, key=lambda node: (-display_graph.degree(node), node))[: min(15, len(display_graph))]
        label_positions = repelled_label_positions(positions, labels)
        for node in labels:
            axis.annotate(
                node, xy=positions[node], xytext=label_positions[node], fontsize=6.5,
                color="#183B56", ha="center", va="center",
                arrowprops={"arrowstyle": "-", "color": "#9FB4C3", "lw": 0.35, "alpha": 0.7},
                zorder=5,
            )
        figure.colorbar(nodes, ax=axis, shrink=0.55, label="log2 fold-change")
        figure.suptitle(
            f"STRING {prefix}regulated interaction network", x=0.08, y=0.975,
            ha="left", weight="bold", color="#183B56", fontsize=15,
        )
        figure.text(
            0.08, 0.935,
            f"Full graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges; "
            f"display: {display_graph.number_of_nodes()} nodes; 15 highest-degree labels",
            ha="left", color="#697783", fontsize=8,
        )
        axis.margins(0.14)
        axis.axis("off")
        figure.subplots_adjust(top=0.90)
        figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
        figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
        plt.close(figure)
    return {
        "direction": direction, "input_genes": len(input_rows), "mapped": len(input_rows) - len(unmapped),
        "unmapped": len(unmapped), "unconnected": len(unconnected), "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(), "communities": len(set(communities.values())),
        "display_nodes": len(display_graph), "network_api_calls": network_api_calls,
        "network_batch_size": STRING_NETWORK_BATCH_SIZE,
    }


def enrichment_panel(genes: list[str], taxonomy: int, cache_dir: Path, offline: bool, refresh: bool, tables: Path, figures: Path) -> int:
    rows = parse_response(cached_post(
        "enrichment", {"identifiers": "\r".join(genes), "species": taxonomy}, cache_dir, offline, refresh
    )) if genes else []
    normalized = []
    for row in rows:
        try:
            fdr = float(row.get("fdr", "nan"))
        except ValueError:
            continue
        normalized.append({
            "category": row.get("category", ""), "term": row.get("term", ""),
            "description": row.get("description", ""), "gene_count": row.get("number_of_genes", ""),
            "background_count": row.get("number_of_genes_in_background", ""), "fdr": fdr,
            "negative_log10_fdr": -math.log10(max(fdr, 1e-300)), "input_genes": row.get("inputGenes", ""),
        })
    normalized.sort(key=lambda row: (row["category"], row["fdr"]))
    displayed = []
    by_category: dict[str, list[dict[str, Any]]] = {}
    for row in normalized:
        by_category.setdefault(row["category"], []).append(row)
    # Preserve every returned term in the audit table. The figure uses a
    # deterministic, category-balanced, de-duplicated view so broad STRING
    # ontologies cannot expand a panel into dozens of redundant labels.
    seen_descriptions: set[str] = set()
    category_candidates: dict[str, list[dict[str, Any]]] = {}
    for category, values in by_category.items():
        category_candidates[category] = []
        for row in values:
            description_key = re.sub(r"[^a-z0-9]+", " ", row["description"].strip().casefold()).strip()
            if not description_key or description_key in seen_descriptions:
                continue
            seen_descriptions.add(description_key)
            category_candidates[category].append(row)
    categories = sorted(category_candidates)
    for rank in range(3):
        for category in categories:
            values = category_candidates[category]
            if rank < len(values):
                displayed.append(values[rank])
            if len(displayed) == 24:
                break
        if len(displayed) == 24:
            break
    displayed.sort(key=lambda row: (row["category"], row["fdr"]))
    fields = ["category", "term", "description", "gene_count", "background_count", "fdr", "negative_log10_fdr", "input_genes"]
    write_tsv(tables / "string_enrichment.tsv", normalized, fields)
    write_tsv(tables / "string_enrichment_displayed.tsv", displayed, fields)
    stem = figures / "string_enrichment"
    if not displayed:
        empty_plot(stem, "STRING enrichment", "No STRING enrichment terms returned")
    else:
        labels = [f"[{row['category']}] {row['description'][:56]}" for row in displayed]
        values = [row["negative_log10_fdr"] for row in displayed]
        sizes = [20 + 12 * float(row["gene_count"] or 0) for row in displayed]
        figure, axis = plt.subplots(figsize=(8.5, max(4.8, 0.34 * len(displayed) + 1.8)))
        y = list(range(len(displayed)))
        scatter = axis.scatter(values, y, s=sizes, c=values, cmap="magma", edgecolors="white", linewidths=0.5)
        axis.set_yticks(y, labels)
        axis.invert_yaxis()
        axis.set_xlabel("−log10(FDR)")
        axis.set_title("STRING functional enrichment", loc="left", weight="bold", color="#183B56")
        axis.grid(axis="x", color="#E7ECF0", linewidth=0.6)
        axis.spines[["top", "right", "left"]].set_visible(False)
        figure.colorbar(scatter, ax=axis, shrink=0.55, label="−log10(FDR)")
        figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
        figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
        plt.close(figure)
    return len(normalized)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--de", required=True)
    parser.add_argument("--contrast-id", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--cache-dir", required=True)
    args = parser.parse_args()
    project = load_project(args.project_config)
    outdir = Path(args.outdir).resolve()
    tables, figures = outdir / "tables", outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True); figures.mkdir(parents=True, exist_ok=True)
    de = read_tsv(Path(args.de))
    taxonomy = int(project.config["species"]["taxonomy_id"])
    settings = project.config["analysis"].get("settings", {}).get("networks", {})
    required_score = int(settings.get("required_score", 700))
    max_nodes = int(settings.get("max_nodes", 120))
    seed = int(settings.get("seed", project.config["analysis"].get("random_seed", 1)))
    offline = bool(project.config["resources"].get("offline", False))
    refresh = bool(project.config["resources"].get("refresh", False))
    cache_dir = Path(args.cache_dir).resolve()
    summaries = []
    for direction in ("up_in_numerator", "down_in_numerator"):
        genes = [row for row in de if row.get("direction") == direction]
        summaries.append(network_for_direction(direction, genes, taxonomy, required_score, max_nodes, seed, cache_dir, offline, refresh, tables, figures))
    all_genes = list(dict.fromkeys(row["gene_symbol"] for row in de if row.get("direction") != "not_significant" and row.get("gene_symbol")))
    enrichment_terms = enrichment_panel(all_genes, taxonomy, cache_dir, offline, refresh, tables, figures)
    summary = {
        "schema_version": 1, "contrast_id": args.contrast_id, "taxonomy_id": taxonomy,
        "required_score": required_score, "display_node_limit": max_nodes,
        "selection_policy": "all significant genes are submitted and retained in audit tables; max_nodes limits visualization only",
        "enrichment_display_policy": "up to 24 unique descriptions, selected round-robin across STRING categories with at most three ranks per category; all terms remain in string_enrichment.tsv",
        "directions": summaries, "enrichment_terms": enrichment_terms,
        "random_seed": seed,
        "warnings": ["STRING edges represent functional/physical association evidence and are not regulatory or causal edges."],
    }
    (outdir / "networks_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
