#!/usr/bin/env python3
"""Cross-contrast consensus: which genes move the same way across contrasts.

A study with two or more pairwise contrasts (e.g. drug-vs-vehicle in several
tissues, or a knockout compared under several conditions) invites the question
the per-contrast tables cannot answer alone: which genes are reproducibly
regulated across contrasts, and in which contrasts does a given program appear?
This study-level reducer reads every pairwise contrast's ``de_results.tsv`` and
reports the consensus.

Direction is taken verbatim from the engine's own per-contrast significance
decision (the ``direction`` column: ``up_in_numerator`` / ``down_in_numerator``
/ ``not_significant``, which already applies the configured FDR and |log2FC|
cutoffs). No thresholding is re-derived here and no direction is inferred from
group names. Each gene is scored by the number of contrasts in which it is
significant and whether those calls agree in sign; a gene significant in
``min_contrasts`` or more contrasts is a consensus gene.

Outputs (all under ``comparison/``):

* ``tables/consensus_membership.tsv`` -- gene x contrast signed-direction matrix
  with per-gene counts, sign-consistency, and a consensus score;
* ``tables/consensus_genes.tsv`` -- the consensus genes only, ranked;
* ``tables/contrast_overlap.tsv`` -- pairwise Jaccard of each contrast's
  significant gene set and the sign-agreement among the shared genes;
* ``tables/consensus_intersections_displayed.tsv`` -- UpSet intersection sizes
  (genes significant in exactly one combination of contrasts);
* ``figures/consensus_upset.{pdf,png}`` -- UpSet-style intersection bars;
* ``figures/consensus_sign_heatmap.{pdf,png}`` -- top consensus genes coloured by
  signed direction across contrasts;
* ``consensus_summary.json``.

Fully deterministic: contrasts are processed in the sorted order supplied on the
command line and every ranking has an explicit total-order tie-break.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from itertools import combinations
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from tifzoret.config import load_project  # noqa: E402

# House palette (matches the R theme_publication constants and curvature.py).
NAVY = "#183B56"
MID_GREY = "#697783"
# Diverging signed-direction encoding: warm = up in numerator, cool = down,
# neutral grey = not significant (a two-hue diverging pair + neutral midpoint).
UP_COLOUR = "#C0392B"
DOWN_COLOUR = "#2C6FBB"
NS_COLOUR = "#E3E7EA"
BAR_COLOUR = "#40354A"
# Cap on how many UpSet intersections and heatmap genes are drawn (all TABLES are
# complete; only the drawn views are truncated to stay legible).
MAX_INTERSECTIONS = 20
DEFAULT_TOP_GENES = 40
DEFAULT_MIN_CONTRASTS = 2

DIRECTION_SIGN = {"up_in_numerator": 1, "down_in_numerator": -1, "not_significant": 0}


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


def _float(value: str) -> float:
    """Parse a float, mapping blanks/NA to NaN so downstream stats can skip them."""
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return result


def load_contrast(path: Path) -> tuple[str, dict[str, dict[str, Any]]]:
    """Return one contrast's id and a per-gene record keyed by gene_id.

    The contrast id is read from the ``contrast_id`` column inside the file (the
    engine writes it into every DE row), falling back to the parent directory
    name so a hand-run file without the column is still usable.
    """
    rows = read_tsv(path)
    contrast_id = ""
    for row in rows:
        contrast_id = (row.get("contrast_id") or "").strip()
        if contrast_id:
            break
    if not contrast_id:
        # .../contrasts/<id>/analyses/de/tables/de_results.tsv
        parts = path.resolve().parts
        contrast_id = parts[parts.index("contrasts") + 1] if "contrasts" in parts else path.stem
    genes: dict[str, dict[str, Any]] = {}
    for row in rows:
        gene_id = (row.get("gene_id") or "").strip()
        if not gene_id:
            continue
        direction = (row.get("direction") or "not_significant").strip()
        genes[gene_id] = {
            "gene_symbol": (row.get("gene_symbol") or gene_id).strip() or gene_id,
            "sign": DIRECTION_SIGN.get(direction, 0),
            "log2_fold_change": _float(row.get("log2_fold_change", "")),
        }
    return contrast_id, genes


def build_upset(
    membership: dict[str, dict[str, int]], contrast_ids: list[str]
) -> list[dict[str, Any]]:
    """Group significant genes by the exact set of contrasts they are significant in.

    Returns one record per observed non-empty combination: the sorted member
    contrasts, a pipe-joined key, and the number of genes significant in exactly
    that set (and no other), sorted by descending size then combination key.
    """
    buckets: dict[tuple[str, ...], int] = {}
    for signs in membership.values():
        present = tuple(cid for cid in contrast_ids if signs.get(cid, 0) != 0)
        if not present:
            continue
        buckets[present] = buckets.get(present, 0) + 1
    records = [
        {
            "contrasts": "|".join(combo),
            "degree": len(combo),
            "gene_count": count,
        }
        for combo, count in buckets.items()
    ]
    records.sort(key=lambda record: (-record["gene_count"], record["contrasts"]))
    return records


def draw_upset(
    intersections: list[dict[str, Any]], contrast_ids: list[str], stem: Path
) -> int:
    """Render an UpSet-style figure: intersection-size bars over a membership dot matrix."""
    shown = intersections[:MAX_INTERSECTIONS]
    figures_dir = stem.parent
    figures_dir.mkdir(parents=True, exist_ok=True)
    if not shown:
        figure, axis = plt.subplots(figsize=(7.2, 4.6))
        axis.text(0.5, 0.5, "No genes significant in any contrast", ha="center", va="center", color=MID_GREY)
        axis.set_axis_off()
        figure.savefig(f"{stem}.pdf", bbox_inches="tight")
        figure.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
        plt.close(figure)
        return 0

    n = len(shown)
    figure, (bar_axis, dot_axis) = plt.subplots(
        2, 1, figsize=(max(6.0, 0.5 * n + 2.2), 5.6), sharex=True,
        gridspec_kw={"height_ratios": [3, max(1.2, 0.42 * len(contrast_ids))], "hspace": 0.05},
    )
    positions = np.arange(n)
    counts = [record["gene_count"] for record in shown]
    bar_axis.bar(positions, counts, color=BAR_COLOUR, width=0.6)
    for x, count in zip(positions, counts):
        bar_axis.text(x, count, str(count), ha="center", va="bottom", fontsize=7, color=NAVY)
    bar_axis.set_ylabel("Genes in intersection", color=NAVY, fontsize=9)
    bar_axis.set_title("Cross-contrast significant-gene intersections", color=NAVY, fontsize=12, fontweight="bold")
    bar_axis.spines[["top", "right"]].set_visible(False)

    member_sets = [set(record["contrasts"].split("|")) for record in shown]
    for row_index, cid in enumerate(contrast_ids):
        y = len(contrast_ids) - 1 - row_index
        for col_index, members in enumerate(member_sets):
            on = cid in members
            dot_axis.plot(
                col_index, y, "o", markersize=9,
                color=NAVY if on else "#D5DADE",
                markeredgecolor="none",
            )
        # Connect the members of each intersection with a vertical line.
    for col_index, members in enumerate(member_sets):
        rows_on = [len(contrast_ids) - 1 - contrast_ids.index(cid) for cid in members]
        if len(rows_on) > 1:
            dot_axis.plot([col_index, col_index], [min(rows_on), max(rows_on)], color=NAVY, linewidth=1.4, zorder=0)
    dot_axis.set_yticks(range(len(contrast_ids)))
    dot_axis.set_yticklabels(list(reversed(contrast_ids)), fontsize=7)
    dot_axis.set_xticks(positions)
    dot_axis.set_xticklabels([])
    dot_axis.set_ylim(-0.5, len(contrast_ids) - 0.5)
    dot_axis.set_xlim(-0.5, n - 0.5)
    dot_axis.spines[["top", "right", "bottom"]].set_visible(False)
    dot_axis.tick_params(axis="x", length=0)
    figure.savefig(f"{stem}.pdf", bbox_inches="tight")
    figure.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(figure)
    return len(shown)


def draw_sign_heatmap(
    genes: list[dict[str, Any]], contrast_ids: list[str], stem: Path
) -> int:
    """Render the top consensus genes as a signed-direction heatmap (gene x contrast)."""
    figures_dir = stem.parent
    figures_dir.mkdir(parents=True, exist_ok=True)
    shown = genes  # already capped to top_genes by the caller
    if not shown:
        figure, axis = plt.subplots(figsize=(6.0, 4.0))
        axis.text(0.5, 0.5, "No consensus genes", ha="center", va="center", color=MID_GREY)
        axis.set_axis_off()
        figure.savefig(f"{stem}.pdf", bbox_inches="tight")
        figure.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
        plt.close(figure)
        return 0
    matrix = np.array([[gene["_signs"][cid] for cid in contrast_ids] for gene in shown], dtype=float)
    cmap = matplotlib.colors.ListedColormap([DOWN_COLOUR, NS_COLOUR, UP_COLOUR])
    norm = matplotlib.colors.BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
    height = max(3.2, 0.24 * len(shown) + 1.2)
    width = max(4.2, 0.7 * len(contrast_ids) + 2.4)
    figure, axis = plt.subplots(figsize=(width, height))
    axis.imshow(matrix, aspect="auto", cmap=cmap, norm=norm)
    axis.set_xticks(range(len(contrast_ids)))
    axis.set_xticklabels(contrast_ids, rotation=45, ha="right", fontsize=7)
    axis.set_yticks(range(len(shown)))
    axis.set_yticklabels([gene["gene_symbol"] for gene in shown], fontsize=6)
    axis.set_title("Consensus genes across contrasts", color=NAVY, fontsize=12, fontweight="bold")
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="none", markersize=9, color=UP_COLOUR, label="Up in numerator"),
        plt.Line2D([0], [0], marker="s", linestyle="none", markersize=9, color=DOWN_COLOUR, label="Down in numerator"),
        plt.Line2D([0], [0], marker="s", linestyle="none", markersize=9, color=NS_COLOUR, label="Not significant"),
    ]
    axis.legend(handles=handles, bbox_to_anchor=(1.02, 1.0), loc="upper left", fontsize=7, frameon=False)
    figure.savefig(f"{stem}.pdf", bbox_inches="tight")
    figure.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(figure)
    return len(shown)


def main() -> None:
    """Compute cross-contrast consensus tables and figures for one study."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--de", required=True, nargs="+", help="de_results.tsv for each pairwise contrast")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    project = load_project(args.project_config)
    settings = (project.config.get("analysis", {}).get("settings", {}) or {}).get("consensus", {}) or {}
    top_genes = int(settings.get("top_genes", DEFAULT_TOP_GENES))
    min_contrasts = int(settings.get("min_contrasts", DEFAULT_MIN_CONTRASTS))

    outdir = Path(args.outdir)
    tables_dir = outdir / "tables"
    figures_dir = outdir / "figures"

    # Load every contrast; order contrasts deterministically by id for stable
    # columns, matrices, and intersection keys regardless of input path order.
    loaded = [load_contrast(Path(path)) for path in args.de]
    loaded.sort(key=lambda item: item[0])
    contrast_ids = [contrast_id for contrast_id, _ in loaded]
    per_contrast = {contrast_id: genes for contrast_id, genes in loaded}

    # Union of genes across contrasts, with a stable symbol per gene.
    symbol_of: dict[str, str] = {}
    for _, genes in loaded:
        for gene_id, record in genes.items():
            symbol_of.setdefault(gene_id, record["gene_symbol"])
    all_gene_ids = sorted(symbol_of)

    # Signed-direction membership matrix and per-gene aggregates.
    membership_signs: dict[str, dict[str, int]] = {}
    membership_rows: list[dict[str, Any]] = []
    for gene_id in all_gene_ids:
        signs = {cid: int(per_contrast[cid].get(gene_id, {}).get("sign", 0)) for cid in contrast_ids}
        membership_signs[gene_id] = signs
        significant = [s for s in signs.values() if s != 0]
        n_up = sum(1 for s in significant if s > 0)
        n_down = sum(1 for s in significant if s < 0)
        n_sig = len(significant)
        lfcs = [
            per_contrast[cid][gene_id]["log2_fold_change"]
            for cid in contrast_ids
            if gene_id in per_contrast[cid]
        ]
        finite_lfcs = [v for v in lfcs if not np.isnan(v)]
        mean_abs_lfc = float(np.mean([abs(v) for v in finite_lfcs])) if finite_lfcs else float("nan")
        sign_consistent = n_sig > 0 and (n_up == 0 or n_down == 0)
        # Consensus score: significant-contrast count, rewarded for sign agreement
        # (a gene up in some contrasts and down in others is genuinely discordant,
        # so its consistent-direction count -- max(n_up, n_down) -- is the score).
        score = max(n_up, n_down)
        row: dict[str, Any] = {
            "gene_id": gene_id,
            "gene_symbol": symbol_of[gene_id],
            "n_significant": n_sig,
            "n_up": n_up,
            "n_down": n_down,
            "sign_consistent": sign_consistent,
            "consensus_score": score,
            "mean_abs_log2fc": mean_abs_lfc,
        }
        for cid in contrast_ids:
            row[cid] = {1: "up", -1: "down", 0: "ns"}[signs[cid]]
        membership_rows.append(row)

    membership_fields = [
        "gene_id", "gene_symbol", "n_significant", "n_up", "n_down",
        "sign_consistent", "consensus_score", "mean_abs_log2fc", *contrast_ids,
    ]
    membership_rows.sort(
        key=lambda r: (-r["consensus_score"], -r["n_significant"], -_safe(r["mean_abs_log2fc"]), r["gene_symbol"])
    )
    write_tsv(tables_dir / "consensus_membership.tsv", membership_rows, membership_fields)

    # Consensus genes: significant (in a consistent direction) in >= min_contrasts.
    consensus_rows = [
        r for r in membership_rows if r["consensus_score"] >= min_contrasts
    ]
    write_tsv(tables_dir / "consensus_genes.tsv", consensus_rows, membership_fields)

    # Pairwise contrast overlap (Jaccard of significant sets + sign agreement).
    overlap_rows: list[dict[str, Any]] = []
    sig_sets = {
        cid: {gid for gid in all_gene_ids if membership_signs[gid][cid] != 0}
        for cid in contrast_ids
    }
    for a, b in combinations(contrast_ids, 2):
        set_a, set_b = sig_sets[a], sig_sets[b]
        shared = set_a & set_b
        union = set_a | set_b
        agree = sum(1 for gid in shared if membership_signs[gid][a] == membership_signs[gid][b])
        overlap_rows.append({
            "contrast_a": a,
            "contrast_b": b,
            "n_a": len(set_a),
            "n_b": len(set_b),
            "intersection": len(shared),
            "union": len(union),
            "jaccard": (len(shared) / len(union)) if union else 0.0,
            "sign_agreement": (agree / len(shared)) if shared else float("nan"),
        })
    overlap_rows.sort(key=lambda r: (-r["jaccard"], r["contrast_a"], r["contrast_b"]))
    write_tsv(
        tables_dir / "contrast_overlap.tsv",
        overlap_rows,
        ["contrast_a", "contrast_b", "n_a", "n_b", "intersection", "union", "jaccard", "sign_agreement"],
    )

    # UpSet intersection sizes.
    intersections = build_upset(membership_signs, contrast_ids)
    write_tsv(
        tables_dir / "consensus_intersections_displayed.tsv",
        intersections,
        ["contrasts", "degree", "gene_count"],
    )

    # Figures.
    drawn_intersections = draw_upset(intersections, contrast_ids, figures_dir / "consensus_upset")
    heatmap_genes = [
        {**r, "_signs": membership_signs[r["gene_id"]]}
        for r in consensus_rows[:top_genes]
    ]
    drawn_genes = draw_sign_heatmap(heatmap_genes, contrast_ids, figures_dir / "consensus_sign_heatmap")

    discordant = sum(1 for r in membership_rows if r["n_up"] > 0 and r["n_down"] > 0)
    summary = {
        "schema_version": 1,
        "method": "cross-contrast direction consensus (engine significance calls)",
        "contrasts": contrast_ids,
        "n_contrasts": len(contrast_ids),
        "genes_considered": len(all_gene_ids),
        "min_contrasts": min_contrasts,
        "consensus_genes": len(consensus_rows),
        "discordant_genes": discordant,
        "intersections": len(intersections),
        "drawn_intersections": drawn_intersections,
        "drawn_consensus_genes": drawn_genes,
    }
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "consensus_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Computed cross-contrast consensus over {len(contrast_ids)} contrasts, {len(consensus_rows)} consensus genes")


def _safe(value: float) -> float:
    """Return ``value`` for sorting, mapping NaN to 0 so it sorts last under negation."""
    return 0.0 if value != value else value


if __name__ == "__main__":
    main()
