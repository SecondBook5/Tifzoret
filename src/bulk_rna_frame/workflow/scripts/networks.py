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
    """Write ``rows`` as a tab-delimited file with a ``fields`` header, creating parent directories."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a tab-delimited file into a list of column-keyed dictionaries."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_bytes(value: bytes) -> str:
    """Return the hex SHA-256 digest of ``value``."""
    return hashlib.sha256(value).hexdigest()


def cached_post(
    endpoint: str,
    parameters: dict[str, object],
    cache_dir: Path,
    offline: bool,
    refresh: bool,
) -> str:
    """Return the STRING ``endpoint`` response for ``parameters``, caching the TSV
    body and a JSON provenance receipt on disk keyed by a hash of the request;
    serve the cache unless ``refresh`` and raise in ``offline`` mode if it is absent."""
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
    """Parse tab-delimited STRING response text into column-keyed rows (empty text yields none)."""
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
    """Save a titled placeholder figure carrying ``message`` as PDF and PNG at ``path``."""
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


NODE_FIELDS = [
    "gene_symbol", "preferredName", "stringId", "queryItem", "seed_group",
    "log2_fold_change", "padj", "stat", "leading_edge_frequency",
    "degree", "degree_centrality", "betweenness_centrality", "community",
]


def network_for_direction(
    seed_group: str,
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
    """Build the STRING association (not regulatory) network induced on one
    directional seed set: map symbols to STRING IDs, fetch induced edges, compute
    Louvain communities and centralities, write the seed/node/edge tables, and
    render the top-``max_display`` subgraph coloured by log2 fold-change."""
    prefix = seed_group
    de_by_symbol: dict[str, dict[str, str]] = {}
    for row in genes:
        symbol = row.get("gene_symbol", "")
        if symbol:
            de_by_symbol.setdefault(symbol, row)
    symbols = list(dict.fromkeys(row["gene_symbol"] for row in genes if row.get("gene_symbol")))
    input_rows = [
        {
            "gene_symbol": symbol,
            "direction": direction,
            "log2_fold_change": de_by_symbol[symbol].get("log2_fold_change", ""),
            "padj": de_by_symbol[symbol].get("adjusted_p_value", ""),
            "stat": de_by_symbol[symbol].get("statistic", ""),
        }
        for symbol in symbols
    ]
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
    # One node row per MAPPED seed gene (connected or isolated), matching the
    # legacy per-seed node table — isolated mapped seeds get degree 0. Keyed on
    # STRING preferredName so it lines up with the edge endpoints and the graph.
    node_rows = []
    for row in input_rows:
        if not row["mapped"]:
            continue
        name = row["preferred_name"]
        in_graph = name in graph
        node_rows.append({
            "gene_symbol": row["gene_symbol"],
            "preferredName": name,
            "stringId": row["string_id"],
            "queryItem": row["gene_symbol"],
            "seed_group": seed_group,
            "log2_fold_change": row["log2_fold_change"],
            "padj": row["padj"],
            "stat": row["stat"],
            "leading_edge_frequency": "",
            "degree": graph.degree(name) if in_graph else 0,
            "degree_centrality": degree.get(name, 0.0),
            "betweenness_centrality": betweenness.get(name, 0.0),
            "community": communities.get(name, 0),
        })
    node_rows.sort(key=lambda r: (-int(r["degree"]), r["preferredName"]))
    connected = set(graph)
    for row in input_rows:
        row["connected"] = row["preferred_name"] in connected or row["gene_symbol"] in connected
    unmapped = [row for row in input_rows if not row["mapped"]]
    unconnected = [row for row in input_rows if row["mapped"] and not row["connected"]]

    selected = sorted(graph, key=lambda node: (-graph.degree(node), node))[:max_display]
    display_graph = graph.subgraph(selected).copy()
    display_nodes = [row for row in node_rows if row["preferredName"] in display_graph]
    display_edges = [row for row in edge_rows if row["source"] in display_graph and row["target"] in display_graph]
    write_tsv(tables / f"string_{prefix}_input_genes.tsv", input_rows, ["gene_symbol", "direction", "log2_fold_change", "mapped", "string_id", "preferred_name", "connected"])
    write_tsv(tables / f"string_{prefix}_unmapped_genes.tsv", unmapped, ["gene_symbol", "direction", "log2_fold_change", "mapped", "string_id", "preferred_name", "connected"])
    write_tsv(tables / f"string_{prefix}_unconnected_genes.tsv", unconnected, ["gene_symbol", "direction", "log2_fold_change", "mapped", "string_id", "preferred_name", "connected"])
    write_tsv(tables / f"string_{prefix}_nodes.tsv", node_rows, NODE_FIELDS)
    write_tsv(tables / f"string_{prefix}_edges.tsv", edge_rows, ["source", "target", "combined_score", "neighborhood_score", "fusion_score", "cooccurrence_score", "coexpression_score", "experimental_score", "database_score", "textmining_score"])
    write_tsv(tables / f"string_{prefix}_nodes_displayed.tsv", display_nodes, NODE_FIELDS)
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
    """Run STRING functional enrichment on ``genes`` and render the term dot plot.

    Normalise FDR to -log10, write the full and displayed enrichment tables, plot a
    deterministic category-balanced, de-duplicated set of up to 24 terms, and return
    the count of returned terms."""
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


def _num(value: object) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return float("nan")


def de_significance_threshold(config: dict[str, Any]) -> float:
    """Return the DE FDR used to define STRING seed sets.

    Mirrors the legacy `select_string_seed_symbols` cutoff (hard-coded 0.05 in
    ``workflow/stages/ontology/ontology_string.R``). Reads the study's DE FDR
    when present so the seed set tracks the DE significance the figures report.
    """
    for path in (("analysis", "settings", "de", "fdr"), ("figures", "de", "fdr")):
        node: Any = config
        for key in path:
            if isinstance(node, dict) and key in node:
                node = node[key]
            else:
                node = None
                break
        value = _num(node)
        if not math.isnan(value):
            return value
    return 0.05


def select_string_seed_symbols(
    de_rows: list[dict[str, str]],
    fgsea_rows: list[dict[str, str]],
    fdr: float,
    max_nodes: int,
) -> dict[str, Any]:
    """Reproduce the legacy STRING seed policy (top-``max_nodes`` per group).

    Faithful port of ``select_string_seed_symbols`` in
    ``workflow/stages/ontology/ontology_string.R``:

    * up / down  — significant DE genes (``adjusted_p_value < fdr``) of the given
      sign, ordered by adjusted p-value ascending then ``|log2FC|`` descending;
      the |log2FC| threshold is deliberately NOT applied to the seed set.
    * leading edge — union of the leading-edge genes of the six most significant
      database (GMT) fgsea pathways, frequency-ranked (frequency descending, then
      DE adjusted p-value ascending, then ``|statistic|`` descending). Leading-edge
      genes need not be DE-significant. Configured gene programs are excluded so
      only Hallmark/Reactome-style pathways drive the selection, matching legacy.
    """

    def directional(positive: bool) -> list[str]:
        """Return up to ``max_nodes`` significant DE symbols of the requested sign, ordered by adjusted p-value then |log2FC|."""
        ranked: list[tuple[float, float, str]] = []
        for row in de_rows:
            symbol = row.get("gene_symbol", "")
            padj = _num(row.get("adjusted_p_value"))
            lfc = _num(row.get("log2_fold_change"))
            if not symbol or math.isnan(padj) or math.isnan(lfc) or padj >= fdr:
                continue
            if (lfc > 0) != positive:
                continue
            ranked.append((padj, -abs(lfc), symbol))
        ranked.sort(key=lambda item: (item[0], item[1]))
        ordered: list[str] = []
        seen: set[str] = set()
        for _, _, symbol in ranked:
            if symbol not in seen:
                seen.add(symbol)
                ordered.append(symbol)
        return ordered[:max_nodes]

    de_padj: dict[str, float] = {}
    de_stat: dict[str, float] = {}
    for row in de_rows:
        symbol = row.get("gene_symbol", "")
        if symbol:
            de_padj.setdefault(symbol, _num(row.get("adjusted_p_value")))
            de_stat.setdefault(symbol, _num(row.get("statistic")))

    # Restrict the leading-edge pathway pool to the study's assembled MSigDB
    # collections (provider "custom", e.g. Hallmark + Reactome), matching the
    # legacy derivation. Auto-added ontology providers (go, kegg) and configured
    # gene programs are excluded so they cannot displace the curated pathways.
    scored_pathways: list[tuple[float, float, dict[str, str]]] = []
    for row in fgsea_rows:
        if row.get("gene_set_source") != "custom":
            continue
        padj = _num(row.get("padj"))
        if math.isnan(padj) or padj >= fdr:
            continue
        scored_pathways.append((padj, -abs(_num(row.get("NES"))), row))
    scored_pathways.sort(key=lambda item: (item[0], item[1]))

    frequency: dict[str, int] = {}
    for _, _, row in scored_pathways[:6]:
        for symbol in (row.get("leadingEdge", "") or "").split(";"):
            symbol = symbol.strip()
            if symbol:
                frequency[symbol] = frequency.get(symbol, 0) + 1

    def leading_edge_key(symbol: str) -> tuple[int, float, float]:
        """Sort key ranking a leading-edge symbol by frequency, then DE adjusted p-value, then |statistic|."""
        padj = de_padj.get(symbol, float("nan"))
        stat = de_stat.get(symbol, float("nan"))
        return (
            -frequency[symbol],
            padj if not math.isnan(padj) else float("inf"),
            -abs(stat) if not math.isnan(stat) else 0.0,
        )

    leading_edge = sorted(frequency, key=leading_edge_key)[:max_nodes]

    return {
        "up": directional(True),
        "down": directional(False),
        "leading_edge": leading_edge,
        "leading_edge_frequency": {symbol: frequency[symbol] for symbol in leading_edge},
    }


ENRICHMENT_FIELDS = [
    "category", "term", "number_of_genes", "number_of_genes_in_background",
    "ncbiTaxonId", "inputGenes", "preferredNames", "p_value", "fdr",
    "description", "seed_group",
]

SEED_FIELDS = ["gene_symbol", "seed_group", "leading_edge_frequency"]


def map_symbols_to_ids(
    symbols: list[str], taxonomy: int, cache_dir: Path, offline: bool, refresh: bool
) -> list[dict[str, str]]:
    """Map gene symbols to STRING identifiers via the cached ``get_string_ids`` endpoint."""
    if not symbols:
        return []
    return parse_response(cached_post(
        "get_string_ids",
        {"identifiers": "\r".join(symbols), "species": taxonomy, "limit": 1, "echo_query": 1},
        cache_dir, offline, refresh,
    ))


def directional_enrichment(
    seed_symbols: list[str],
    frequency: dict[str, int],
    seed_group: str,
    taxonomy: int,
    cache_dir: Path,
    offline: bool,
    refresh: bool,
    tables: Path,
) -> dict[str, Any]:
    """Run STRING functional enrichment on one directional seed set.

    Emits the enrichment table in the legacy bespoke schema (raw STRING
    ``/enrichment`` columns plus ``seed_group``) so the R renderer can consume
    it unchanged, and an auditable seed-gene table.
    """
    mapping = map_symbols_to_ids(seed_symbols, taxonomy, cache_dir, offline, refresh)
    string_ids = list(dict.fromkeys(row.get("stringId", "") for row in mapping if row.get("stringId")))
    rows = parse_response(cached_post(
        "enrichment", {"identifiers": "\r".join(string_ids), "species": taxonomy}, cache_dir, offline, refresh
    )) if string_ids else []
    for row in rows:
        row["seed_group"] = seed_group
    write_tsv(tables / f"string_{seed_group}_enrichment.tsv", rows, ENRICHMENT_FIELDS)
    seed_rows = [
        {
            "gene_symbol": symbol,
            "seed_group": seed_group,
            "leading_edge_frequency": frequency.get(symbol, "") if frequency else "",
        }
        for symbol in seed_symbols
    ]
    write_tsv(tables / f"string_{seed_group}_seed_genes.tsv", seed_rows, SEED_FIELDS)
    return {
        "seed_group": seed_group,
        "seed_genes": len(seed_symbols),
        "mapped_genes": len(mapping),
        "enrichment_terms": len(rows),
    }


def main() -> None:
    """CLI entry point: build the directional STRING networks, the enrichment panel,
    and per-group directional enrichment for one contrast, then write networks_summary.json."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--de", required=True)
    parser.add_argument("--fgsea", default=None, help="pathways fgsea.tsv for the leading-edge STRING seed set")
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
    # STRING networks and the enrichment facets are built from the SAME capped
    # directional seed sets (legacy parity), so the network figures use the
    # seed-network confidence (STRING default 400) rather than required_score.
    seed_required_score = int(settings.get("seed_required_score", 400))
    seed = int(settings.get("seed", project.config["analysis"].get("random_seed", 1)))
    offline = bool(project.config["resources"].get("offline", False))
    refresh = bool(project.config["resources"].get("refresh", False))
    cache_dir = Path(args.cache_dir).resolve()

    # Compute the directional seed sets first: the up/down STRING networks
    # (Panels B/C) and the three-facet enrichment (Panel A) are both built from
    # these top-N seeds, matching the legacy pipeline that produced the bespoke
    # figures. The seed policy is a faithful port of select_string_seed_symbols.
    seed_fdr = de_significance_threshold(project.config)
    seed_max = int(settings.get("seed_max_nodes", 50))
    fgsea_rows = read_tsv(Path(args.fgsea)) if args.fgsea and Path(args.fgsea).is_file() else []
    seeds = select_string_seed_symbols(de, fgsea_rows, seed_fdr, seed_max)
    de_by_symbol = {}
    for row in de:
        symbol = row.get("gene_symbol", "")
        if symbol:
            de_by_symbol.setdefault(symbol, row)

    summaries = []
    for seed_group, direction in (("up", "up_in_numerator"), ("down", "down_in_numerator")):
        genes = [de_by_symbol[symbol] for symbol in seeds[seed_group] if symbol in de_by_symbol]
        summaries.append(network_for_direction(seed_group, direction, genes, taxonomy, seed_required_score, max_nodes, seed, cache_dir, offline, refresh, tables, figures))
    all_genes = list(dict.fromkeys(row["gene_symbol"] for row in de if row.get("direction") != "not_significant" and row.get("gene_symbol")))
    enrichment_terms = enrichment_panel(all_genes, taxonomy, cache_dir, offline, refresh, tables, figures)

    # Directional STRING enrichment on the same seed sets, emitted per group in
    # the bespoke schema for the R renderer.
    seed_summaries = [
        directional_enrichment(
            seeds[group], seeds["leading_edge_frequency"] if group == "leading_edge" else {},
            group, taxonomy, cache_dir, offline, refresh, tables,
        )
        for group in ("up", "down", "leading_edge")
    ]

    summary = {
        "schema_version": 3, "contrast_id": args.contrast_id, "taxonomy_id": taxonomy,
        "required_score": required_score, "network_required_score": seed_required_score,
        "display_node_limit": max_nodes,
        "selection_policy": (
            f"up/down STRING networks are induced on the top-{seed_max} directional seed sets "
            f"(same seeds as the enrichment facets) at STRING confidence {seed_required_score}; "
            "max_nodes limits visualization only"
        ),
        "enrichment_display_policy": "up to 24 unique descriptions, selected round-robin across STRING categories with at most three ranks per category; all terms remain in string_enrichment.tsv",
        "seed_enrichment_policy": (
            f"directional STRING enrichment on top-{seed_max} seed sets: up/down = "
            f"significant DE genes (adjusted_p_value < {seed_fdr}) by padj then |log2FC|; "
            "leading_edge = union of leading edges of the six most significant GMT fgsea pathways, frequency-ranked"
        ),
        "seed_fdr": seed_fdr, "seed_max_nodes": seed_max,
        "directions": summaries, "enrichment_terms": enrichment_terms,
        "seed_enrichment": seed_summaries,
        "random_seed": seed,
        "warnings": ["STRING edges represent functional/physical association evidence and are not regulatory or causal edges."],
    }
    (outdir / "networks_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
