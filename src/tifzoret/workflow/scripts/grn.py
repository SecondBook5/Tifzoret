#!/usr/bin/env python3
"""Create auditable program-aware regulon views (DoRothEA GRN, Figure 2 Panel E).

DATA layer of the GRN render-seam. This selects the displayed regulator-target
subgraph (top differential regulators + their top DE-annotated targets), assigns
each node to a transcriptional program by anchored network diffusion, runs the
exact target-overlap separation test, and emits the node/edge/sector tables the R
renderer (grn_radial.R) consumes. It also writes an auditable rectangular
matplotlib view; the polished radial panel is rendered downstream in R.

Program assignment mirrors the bespoke publication figure
(the prepare_dorothea_grn_panel routine in the reference figure library): a
personalized PageRank is diffused from each
program's curated anchor set over the likelihood-weighted main component, and a
node is assigned to a program when its epithelial/(epithelial+myogenic) diffusion
probability clears the threshold (else "shared / bridging"). Absent curation, or a
program set that is not the two-program probability model, falls back to the
generic gene-panel majority labelling so the stage still runs for any study.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import networkx as nx

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from tifzoret.config import load_project  # noqa: E402
from tifzoret.figures import normalized_gene_panels  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a tab-delimited file into a list of column-keyed dictionaries."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    """Write ``rows`` as a tab-delimited file with a ``fields`` header, creating parent directories."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def numeric(value: object, default: float) -> float:
    """Coerce ``value`` to a finite float, returning ``default`` on failure or non-finite input."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else default


def program_map(project) -> tuple[dict[str, str], dict[str, str]]:
    """Fallback program labelling from hypothesis gene panels (majority vote)."""
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


def grn_curation(project) -> dict[str, Any] | None:
    """Return the GRN curation block from ``analysis.settings.regulators.grn`` (or None)."""
    settings = ((project.config.get("analysis") or {}).get("settings") or {})
    regulators = settings.get("regulators") or {}
    return regulators.get("grn")


def personalized_pagerank(graph: nx.Graph, seeds: list[str], damping: float) -> dict[str, float] | None:
    """Diffuse a likelihood-weighted personalized PageRank from the ``seeds`` present
    in ``graph`` (uniform seed mass), returning the score map or None when no seed is
    in the graph."""
    present = [seed for seed in seeds if seed in graph]
    if not present:
        return None
    personalization = {node: 0.0 for node in graph.nodes}
    for seed in present:
        personalization[seed] = 1.0 / len(present)
    # High iteration budget / tight tolerance so the power iteration matches the
    # exact solver (igraph PRPACK) the reference figure used closely enough that
    # threshold assignments agree.
    return nx.pagerank(graph, alpha=damping, personalization=personalization, weight="likelihood", max_iter=2000, tol=1.0e-12)


def program_modularity(graph: nx.Graph, labels: dict[str, str]) -> float:
    """Return the weighted modularity of the partition ``labels`` induces on ``graph``
    (0.0 when it has no edges or fewer than two groups)."""
    groups = defaultdict(set)
    for node in graph:
        groups[labels.get(node, "Unassigned")].add(node)
    communities = [nodes for nodes in groups.values() if nodes]
    return nx.community.modularity(graph, communities, weight="weight") if graph.number_of_edges() and len(communities) > 1 else 0.0


def separation_score(jaccard: dict[tuple[str, str], float], group_a: list[str], group_b: list[str]) -> float:
    """Return the target-overlap separation of two groups: mean within-group Jaccard
    minus mean between-group Jaccard."""
    def mean_within(group: list[str]) -> float:
        """Mean pairwise Jaccard within ``group`` (NaN if fewer than two members)."""
        pairs = [jaccard[(a, b)] for a, b in combinations(sorted(group), 2)]
        return sum(pairs) / len(pairs) if pairs else math.nan
    def mean_between(a: list[str], b: list[str]) -> float:
        """Mean Jaccard over all cross pairs between ``a`` and ``b`` (NaN if either is empty)."""
        pairs = [jaccard[(x, y)] for x in a for y in b]
        return sum(pairs) / len(pairs) if pairs else math.nan
    within = [value for value in (mean_within(group_a), mean_within(group_b)) if not math.isnan(value)]
    within_mean = sum(within) / len(within) if within else math.nan
    return within_mean - mean_between(group_a, group_b)


def exact_overlap_test(
    regulators: list[str],
    reg_targets: dict[str, set[str]],
    epithelial: list[str],
    myogenic: list[str],
    cap: int = 3_000_000,
) -> dict[str, Any]:
    """Target-overlap Jaccard separation of the two anchor groups vs all disjoint
    equal-size regulator splits (exact enumeration; Monte Carlo above the cap)."""
    regs = sorted(regulators)
    jaccard: dict[tuple[str, str], float] = {}
    for a, b in combinations(regs, 2):
        sa, sb = reg_targets.get(a, set()), reg_targets.get(b, set())
        union = len(sa | sb)
        value = len(sa & sb) / union if union else 0.0
        jaccard[(a, b)] = value
        jaccard[(b, a)] = value
    observed = separation_score(jaccard, epithelial, myogenic)
    n_a, n_b = len(epithelial), len(myogenic)
    total = math.comb(len(regs), n_a) * math.comb(len(regs) - n_a, n_b) if len(regs) >= n_a + n_b else 0
    eps = math.sqrt(2.2e-16)
    ge = 0
    count = 0
    exact = total <= cap and total > 0
    if exact:
        for candidate_a in combinations(regs, n_a):
            remaining = [reg for reg in regs if reg not in candidate_a]
            for candidate_b in combinations(remaining, n_b):
                if separation_score(jaccard, list(candidate_a), list(candidate_b)) >= observed - eps:
                    ge += 1
                count += 1
    else:
        rng = random.Random(1)
        draws = min(cap, 200_000)
        for _ in range(draws):
            shuffled = regs[:]
            rng.shuffle(shuffled)
            if separation_score(jaccard, shuffled[:n_a], shuffled[n_a:n_a + n_b]) >= observed - eps:
                ge += 1
            count += 1
    p_value = ge / count if count else math.nan
    return {
        "epithelial_anchor_count": n_a,
        "myogenic_anchor_count": n_b,
        "main_component_regulator_count": len(regs),
        "observed_separation": observed,
        "exact_permutations": count,
        "exact_test": exact,
        "exact_p": p_value,
    }


def main() -> None:
    """CLI entry point: select the top differential regulators and their top DE targets,
    assign each node to a program by anchored network diffusion (panel-majority fallback),
    run the program separation test, write the node/edge/sector/separation tables and
    matplotlib views, and emit grn_summary.json."""
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

    # Top regulators by differential activity: adj.P.Val, then |t| (matches the
    # bespoke grn_dorothea.R selection; the |t| tiebreak orders the FDR-tied tail).
    top_n = int(((project.config.get("analysis") or {}).get("settings") or {}).get("regulators", {}).get("top_regulators", 15) or 15)
    ranked = sorted(
        regulator_rows,
        key=lambda row: (numeric(row.get("adj.P.Val"), math.inf), -abs(numeric(row.get("t"), 0.0))),
    )
    top_regulators = [row["regulator"] for row in ranked[:top_n]]
    regulator_activity = {row["regulator"]: numeric(row.get("logFC"), 0.0) for row in regulator_rows}
    regulator_padj = {row["regulator"]: numeric(row.get("adj.P.Val"), math.nan) for row in regulator_rows}

    # Top DE-annotated measured targets per regulator (padj, then |log2FC|, name).
    candidate_edges = [row for row in edges if row["source"] in top_regulators]
    targets_per_regulator = 6
    selected_edges: list[dict[str, Any]] = []
    for regulator in top_regulators:
        current = [row for row in candidate_edges if row["source"] == regulator and row["target"].upper() in de and de[row["target"].upper()].get("adjusted_p_value") not in (None, "", "NA")]
        current.sort(key=lambda row: (
            numeric(de[row["target"].upper()].get("adjusted_p_value"), math.inf),
            -abs(numeric(de[row["target"].upper()].get("log2_fold_change"), 0.0)),
            row["target"],
        ))
        selected_edges.extend(current[:targets_per_regulator])

    # Enriched edge table (mode of regulation, confidence likelihood, target DE).
    edge_rows: list[dict[str, Any]] = []
    reg_targets: dict[str, set[str]] = defaultdict(set)
    for row in selected_edges:
        source, target = row["source"], row["target"]
        mor = numeric(row.get("mor"), 1.0)
        likelihood = numeric(row.get("likelihood"), 1.0)
        de_row = de.get(target.upper(), {})
        edge_rows.append({
            "source": source, "target": target,
            "mode_of_regulation": mor, "likelihood": likelihood,
            "regulation": "Activating" if mor >= 0 else "Repressing",
            "target_log2fc": de_row.get("log2_fold_change", ""),
            "target_padj": de_row.get("adjusted_p_value", ""),
            "weight": max(abs(mor) * abs(likelihood), 0.05),
        })
        reg_targets[source].add(target)

    # Directed graph for the rectangular view + degree; likelihood-weighted
    # undirected graph for diffusion (matches the reference weights = likelihood).
    graph = nx.DiGraph()
    diffusion_graph = nx.Graph()
    for row in edge_rows:
        graph.add_edge(row["source"], row["target"], weight=row["weight"], regulation=row["regulation"])
        if not diffusion_graph.has_edge(row["source"], row["target"]):
            diffusion_graph.add_edge(row["source"], row["target"], likelihood=row["likelihood"])

    curation = grn_curation(project)
    programs = (curation or {}).get("programs") or []
    damping = numeric((curation or {}).get("damping"), 0.85)
    threshold = numeric((curation or {}).get("assignment_threshold"), 0.80)
    shared_label = ((curation or {}).get("shared") or {}).get("label", "Shared / bridging context")
    peripheral_label = ((curation or {}).get("peripheral") or {}).get("label", "Peripheral regulator component")

    labels: dict[str, str] = {}
    probs: dict[str, float] = {}
    overlap: dict[str, Any] = {}
    program_colors: dict[str, str] = {}
    used_diffusion = False

    main_names: set[str] = set()
    if diffusion_graph.number_of_edges():
        components = sorted(nx.connected_components(diffusion_graph), key=len, reverse=True)
        main_names = set(components[0]) if components else set()

    # Anchor-connectivity guarantee. The reference figure asserts (stopifnot)
    # that every curated anchor lands in the diffusion main component -- anchors
    # define the programs, so a disconnected anchor invalidates the assignment.
    # The diffusion graph is the *displayed* top-k subgraph, so whether an anchor
    # connects can hinge on which regulator was selected into the top-N tail: two
    # independent VIPER implementations may swap a near-tied regulator at the
    # boundary, dropping the one whose shared target bridged an anchor in. That
    # strands a genuine anchor (e.g. Tead4) as its own component -- a display
    # artifact, not biology. Reconnect each stranded anchor minimally: add its
    # single strongest shared measured-regulon target to a same-program anchor
    # already in the main component (the bridge node is diffusion-only, never
    # displayed). Non-anchor regulators whose top targets are idiosyncratic
    # (e.g. housekeeping zinc-fingers) are deliberately left peripheral, exactly
    # as in the reference.
    bridged_anchors: list[str] = []
    curation_early = grn_curation(project)
    programs_early = (curation_early or {}).get("programs") or []
    if len(programs_early) == 2 and diffusion_graph.number_of_edges():
        measured_targets: dict[str, dict[str, float]] = defaultdict(dict)
        for row in edges:
            likelihood = numeric(row.get("likelihood"), 0.5)
            if likelihood > measured_targets[row["source"]].get(row["target"], -1.0):
                measured_targets[row["source"]][row["target"]] = likelihood
        for program in programs_early:
            anchors = [a for a in program.get("anchors", []) if a in diffusion_graph]
            for anchor in [a for a in anchors if a not in main_names]:
                partners = [c for c in anchors if c != anchor and c in main_names] or [c for c in anchors if c != anchor]
                best: tuple[float, str, str] | None = None
                for partner in partners:
                    shared = set(measured_targets.get(anchor, {})) & set(measured_targets.get(partner, {}))
                    for target in sorted(shared):
                        score = measured_targets[anchor][target] * measured_targets[partner][target]
                        if best is None or score > best[0]:
                            best = (score, target, partner)
                if best is not None:
                    _, target, partner = best
                    diffusion_graph.add_edge(anchor, target, likelihood=measured_targets[anchor][target])
                    diffusion_graph.add_edge(target, partner, likelihood=measured_targets[partner][target])
                    bridged_anchors.append(anchor)
        if bridged_anchors:
            components = sorted(nx.connected_components(diffusion_graph), key=len, reverse=True)
            main_names = set(components[0]) if components else set()

    main_graph = diffusion_graph.subgraph(main_names).copy() if main_names else nx.Graph()

    if len(programs) == 2 and main_graph.number_of_nodes():
        prog_a, prog_b = programs[0], programs[1]
        score_a = personalized_pagerank(main_graph, prog_a.get("anchors", []), damping)
        score_b = personalized_pagerank(main_graph, prog_b.get("anchors", []), damping)
        if score_a is not None and score_b is not None:
            used_diffusion = True
            program_colors[prog_a["label"]] = prog_a.get("colour", "#34824B")
            program_colors[prog_b["label"]] = prog_b.get("colour", "#D51B70")
            program_colors[shared_label] = ((curation or {}).get("shared") or {}).get("colour", "#9A7B43")
            program_colors[peripheral_label] = ((curation or {}).get("peripheral") or {}).get("colour", "#7F8C96")
            for node in main_graph.nodes:
                a = score_a.get(node, 0.0); b = score_b.get(node, 0.0)
                probability = a / (a + b) if (a + b) > 0 else 0.5
                probs[node] = probability
                if probability >= threshold:
                    labels[node] = prog_a["label"]
                elif probability <= 1.0 - threshold:
                    labels[node] = prog_b["label"]
                else:
                    labels[node] = shared_label
            for node in graph:
                if node not in main_names:
                    labels[node] = peripheral_label
            # Exact target-overlap separation test over the anchor groups actually
            # present as main-component regulators.
            main_regulators = sorted(node for node in main_names if node in top_regulators)
            epi_present = [a for a in prog_a.get("anchors", []) if a in main_regulators]
            myo_present = [a for a in prog_b.get("anchors", []) if a in main_regulators]
            if len(main_regulators) >= len(epi_present) + len(myo_present) and epi_present and myo_present:
                overlap = exact_overlap_test(main_regulators, reg_targets, epi_present, myo_present)

    if not used_diffusion:
        # Generic fallback: gene-panel majority vote.
        mapping, fallback_colors = program_map(project)
        program_colors = fallback_colors
        target_program = {node: mapping.get(node.upper(), "Unassigned") for node in graph if node not in top_regulators}
        labels = dict(target_program)
        for regulator in top_regulators:
            member = [target_program.get(target, "Unassigned") for target in graph.successors(regulator)] if regulator in graph else []
            labels[regulator] = Counter(member).most_common(1)[0][0] if member else "Unassigned"

    fallback_palette = program_map(project)[1]
    def color_of(program: str) -> str:
        """Return the display colour for ``program``, falling back to the panel palette then grey."""
        return program_colors.get(program, fallback_palette.get(program, "#A0A8AE"))

    undirected = graph.to_undirected()
    community_members = nx.community.louvain_communities(undirected, weight="weight", seed=project.config["analysis"].get("random_seed", 1)) if undirected.number_of_edges() else []
    community_map = {node: index + 1 for index, group in enumerate(community_members) for node in group}

    node_rows = []
    for node in graph:
        is_regulator = node in top_regulators
        de_row = de.get(node.upper(), {})
        value = regulator_activity.get(node, "") if is_regulator else de_row.get("log2_fold_change", "")
        padj = regulator_padj.get(node, "") if is_regulator else de_row.get("adjusted_p_value", "")
        node_rows.append({
            "node": node, "node_type": "regulator" if is_regulator else "target",
            "program": labels.get(node, "Unassigned"), "community": community_map.get(node, 0),
            "value": value if value != "" and value is not None else "",
            "padj": padj if padj != "" and padj is not None else "",
            "log2_fold_change": de_row.get("log2_fold_change", ""),
            "adjusted_p_value": de_row.get("adjusted_p_value", ""),
            "activity_logfc": regulator_activity.get(node, ""),
            "diffusion_probability": probs.get(node, ""),
            "degree": undirected.degree(node),
        })
    write_tsv(tables / "grn_nodes_displayed.tsv", node_rows, [
        "node", "node_type", "program", "community", "value", "padj",
        "log2_fold_change", "adjusted_p_value", "activity_logfc", "diffusion_probability", "degree",
    ])
    write_tsv(tables / "grn_edges_displayed.tsv", edge_rows, [
        "source", "target", "mode_of_regulation", "likelihood", "regulation", "target_log2fc", "target_padj", "weight",
    ])

    # Program separation test. Prefer the exact target-overlap statistic used by
    # the reference radial subtitle; otherwise the label-shuffle modularity test.
    if overlap:
        separation = [overlap]
        p_value = overlap["exact_p"]
    else:
        observed = program_modularity(undirected, labels)
        rng = random.Random(project.config["analysis"].get("random_seed", 1))
        nodes = list(undirected); values = [labels.get(node, "Unassigned") for node in nodes]
        permutations = []
        for _ in range(1000):
            shuffled = values[:]; rng.shuffle(shuffled)
            permutations.append(program_modularity(undirected, dict(zip(nodes, shuffled))))
        p_value = (1 + sum(value >= observed for value in permutations)) / (len(permutations) + 1)
        separation = [{"observed_program_modularity": observed, "permutation_p_value": p_value, "permutations": len(permutations), "communities": len(community_members)}]
    write_tsv(tables / "grn_program_separation_test.tsv", separation, list(separation[0]))

    ordered_programs = []
    if used_diffusion:
        for name in (peripheral_label, programs[0]["label"], shared_label, programs[1]["label"]):
            if name not in ordered_programs:
                ordered_programs.append(name)
    for row in node_rows:
        if row["program"] not in ordered_programs:
            ordered_programs.append(row["program"])
    sector_rows = [{
        "program": program, "color": color_of(program),
        "target_count": sum(row["program"] == program and row["node_type"] == "target" for row in node_rows),
        "nodes": sum(row["program"] == program for row in node_rows),
        "regulators": sum(row["program"] == program and row["node_type"] == "regulator" for row in node_rows),
        "targets": sum(row["program"] == program and row["node_type"] == "target" for row in node_rows),
    } for program in ordered_programs]
    sector_rows = [row for row in sector_rows if row["nodes"] > 0]
    write_tsv(tables / "grn_sector_summary.tsv", sector_rows, ["program", "color", "target_count", "nodes", "regulators", "targets"])

    # Rectangular breadth-preserving matplotlib view (auditable secondary).
    regulators = [node for node in graph if node in top_regulators]
    targets = [node for node in graph if node not in top_regulators]
    positions = {node: (0, index) for index, node in enumerate(regulators)}
    positions.update({node: (1, index * max(len(regulators), 1) / max(len(targets), 1)) for index, node in enumerate(targets)})
    figure, axis = plt.subplots(figsize=(10.5, max(6.2, 0.22 * max(len(targets), len(regulators)) + 2)))
    nx.draw_networkx_edges(graph, positions, edge_color=["#D66B5D" if graph.edges[edge]["regulation"] == "Activating" else "#4E88A8" for edge in graph.edges], width=[0.6 + graph.edges[edge]["weight"] for edge in graph.edges], alpha=0.38, arrows=True, arrowsize=8, ax=axis)
    nx.draw_networkx_nodes(graph, positions, nodelist=regulators, node_shape="s", node_color=[color_of(labels.get(node, "Unassigned")) for node in regulators], node_size=140, edgecolors="white", ax=axis)
    nx.draw_networkx_nodes(graph, positions, nodelist=targets, node_shape="o", node_color=[color_of(labels.get(node, "Unassigned")) for node in targets], node_size=75, edgecolors="white", ax=axis)
    nx.draw_networkx_labels(graph, positions, font_size=6.2, font_color="#183B56", ax=axis)
    axis.set_title("Regulator–target network", loc="left", weight="bold", color="#183B56"); axis.axis("off")
    figure.savefig(figures / "grn_rectangular.pdf", bbox_inches="tight"); figure.savefig(figures / "grn_rectangular.png", dpi=300, bbox_inches="tight"); plt.close(figure)

    # Legacy compact radial (matplotlib); the publication radial is grn_radial.R.
    ordered_for_radial = [row["program"] for row in sector_rows] or ["Unassigned"]
    sector_angles = {program: 2 * math.pi * index / max(len(ordered_for_radial), 1) for index, program in enumerate(ordered_for_radial)}
    radial_positions = {}
    for program in ordered_for_radial:
        program_regulators = [node for node in regulators if labels.get(node, "Unassigned") == program]
        program_targets = [node for node in targets if labels.get(node, "Unassigned") == program]
        center = sector_angles[program]
        for index, node in enumerate(program_regulators):
            angle = center + (index - (len(program_regulators) - 1) / 2) * 0.10
            radial_positions[node] = (0.48 * math.cos(angle), 0.48 * math.sin(angle))
        for index, node in enumerate(program_targets):
            angle = center + (index - (len(program_targets) - 1) / 2) * min(0.08, 0.65 / max(len(program_targets), 1))
            radial_positions[node] = (math.cos(angle), math.sin(angle))
    for node in graph:
        radial_positions.setdefault(node, (0.0, 0.0))
    figure, axis = plt.subplots(figsize=(8.4, 8.4))
    nx.draw_networkx_edges(graph, radial_positions, edge_color=["#D66B5D" if graph.edges[edge]["regulation"] == "Activating" else "#4E88A8" for edge in graph.edges], width=0.65, alpha=0.25, arrows=True, arrowsize=7, connectionstyle="arc3,rad=0.08", ax=axis)
    nx.draw_networkx_nodes(graph, radial_positions, nodelist=regulators, node_shape="s", node_color=[color_of(labels.get(node, "Unassigned")) for node in regulators], node_size=130, edgecolors="white", ax=axis)
    nx.draw_networkx_nodes(graph, radial_positions, nodelist=targets, node_shape="o", node_color=[color_of(labels.get(node, "Unassigned")) for node in targets], node_size=65, edgecolors="white", ax=axis)
    label_nodes = regulators + sorted(targets, key=lambda node: -undirected.degree(node))[:25]
    nx.draw_networkx_labels(graph, radial_positions, labels={node: node for node in label_nodes}, font_size=6, font_color="#183B56", ax=axis)
    for program, angle in sector_angles.items():
        axis.text(1.22 * math.cos(angle), 1.22 * math.sin(angle), program.replace("_", " "), ha="center", va="center", color=color_of(program), weight="bold", fontsize=8)
    axis.set_title("Program-aware radial regulator network", loc="left", weight="bold", color="#183B56"); axis.axis("off"); axis.set_aspect("equal")
    figure.savefig(figures / "grn_radial_legacy.pdf", bbox_inches="tight"); figure.savefig(figures / "grn_radial_legacy.png", dpi=300, bbox_inches="tight"); plt.close(figure)

    (outdir / "grn_summary.json").write_text(json.dumps({
        "schema_version": 2, "contrast_id": args.contrast_id, "displayed_regulators": len(regulators),
        "displayed_targets": len(targets), "displayed_edges": len(edge_rows),
        "selection": f"top {top_n} differential regulators (adj.P.Val, |t|) and top {targets_per_regulator} DE-ranked measured targets per regulator",
        "program_assignment": "anchored network diffusion" if used_diffusion else "gene-panel majority",
        "bridged_anchors": bridged_anchors,
        "full_edges": str(Path(args.edges).resolve()),
        "program_separation_p": p_value,
        "programs": ordered_programs,
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
