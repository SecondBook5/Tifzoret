#!/usr/bin/env python3
"""Render the Tifzoret pipeline diagram for the README.

This is a *curated* view of the workflow — the five phases the README narrates,
with representative edges — drawn in Graphviz so it reads cleanly and matches the
profile-tier palette (a node's color states which profile turns it on). A few
tightly-coupled rules are collapsed into one node (e.g. ``contrast_grn`` and its
radial layout); each node below declares exactly which real rules it stands for.

To keep the picture from silently drifting from the code, the drift guard reads
every ``rule`` declaration in ``src/tifzoret/workflow/`` and asserts that each one
is represented by exactly one node here. A new, renamed, or removed rule therefore
aborts the render until the diagram is updated deliberately — the curation stays
honest without being auto-generated.

Outputs ``docs/dag.svg`` (used by the README) and ``docs/dag.png``.

Usage:  python docs/render_dag.py     (needs graphviz `dot` on PATH)
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RULES_DIR = REPO / "src" / "tifzoret" / "workflow"
SVG_OUT = REPO / "docs" / "dag.svg"
PNG_OUT = REPO / "docs" / "dag.png"

# Profile tiers: fill / stroke / ink. The README legend documents the same palette.
IO = ("#E8EAED", "#555555", "#111111")            # input · resource · provenance
STD = ("#A6CEE3", "#2C6E9B", "#08263B")           # standard
PUB = ("#B7E4C7", "#2D8659", "#0B3320")           # publication
FULL = ("#FAD7A0", "#B9770E", "#5B3A08")          # full

# Each node: (rule-set it represents, tier, shape, HTML-ish label).
# `represents` is the ground truth the drift guard checks against the Snakefile.
# Input *kinds* (BAM/NFC/…) are not rules, so they represent nothing.
def lbl(name, sub=None, bold=True):
    head = f"<b>{name}</b>" if bold else name
    if sub is None:
        return f"<{head}>"
    return f'<{head}<br/><font point-size="9">{sub}</font>>'


# Display names drop the phase-implied `contrast_` prefix (the ③ box supplies it)
# to keep the wide fan compact and consistent with the ④ box's short names — the
# `represents` lists below keep the true rule names, which is what the drift guard
# and the README caption pin to.
NODES = {
    # ① input boundary — the four accepted input kinds (decorative) + the two io rules
    "BAM": ([], IO, "parallelogram", "<aligned BAMs>"),
    "NFC": ([], IO, "parallelogram", "<nf-core/rnaseq>"),
    "ARC": ([], IO, "parallelogram", "<ZIP · TAR archive>"),
    "CNT": ([], IO, "parallelogram", "<integer count matrix>"),
    "MAT": (["materialize_inputs"], IO, "box",
            lbl("materialize_inputs", "featureCounts · GTF-first")),
    "RES": (["resolve_resources"], IO, "box",
            lbl("resolve_resources", "MSigDB · GO · KEGG · STRING")),
    # ② study scope
    "QC": (["study_qc"], STD, "box", lbl("study_qc", "VST · PCA · distances")),
    "BATCH": (["study_batch"], STD, "box", lbl("study_batch", "removeBatchEffect", bold=False)),
    # ③ per contrast (names shown without the shared `contrast_` prefix)
    "DE": (["contrast_de"], STD, "box", lbl("de", "DESeq2 · 5-class")),
    "CONF": (["contrast_de_confirm"], STD, "box", lbl("de_confirm", "edgeR", bold=False)),
    "OMNI": (["contrast_omnibus"], STD, "box", lbl("omnibus", "LRT", bold=False)),
    "PATH": (["contrast_pathways"], STD, "box", lbl("pathways", "fgsea · GSVA · GSEA")),
    "ONT": (["contrast_ontology"], STD, "box", lbl("ontology", "ORA", bold=False)),
    "SPIA": (["contrast_spia"], STD, "box", lbl("spia", bold=False)),
    "EMAP": (["contrast_enrichment_map"], STD, "box", lbl("enrichment_map", bold=False)),
    "COMP": (["contrast_composition"], PUB, "box", lbl("composition", "cell-state · NNLS", bold=False)),
    "REG": (["contrast_regulators"], PUB, "box", lbl("regulators", "DoRothEA / VIPER", bold=False)),
    "NET": (["contrast_networks", "contrast_string_figures"], PUB, "box",
            lbl("networks", "STRING", bold=False)),
    "GRN": (["contrast_grn", "contrast_grn_radial"], PUB, "box",
            lbl("grn → grn_radial", bold=False)),
    "HYP": (["contrast_hypotheses"], PUB, "box", lbl("hypotheses", "evidence scoring")),
    # ④ advanced (full)
    "SVA": (["contrast_sva"], FULL, "box", lbl("sva", bold=False)),
    "WGCNA": (["contrast_wgcna"], FULL, "box", lbl("wgcna", bold=False)),
    "CURV": (["contrast_curvature"], FULL, "box", lbl("curvature", bold=False)),
    "MED": (["contrast_mediation", "contrast_mediation_power"], FULL, "box", lbl("mediation → power", bold=False)),
    "MULTI": (["contrast_multilayer"], FULL, "box", lbl("multilayer", bold=False)),
    "CONS": (["study_consensus"], FULL, "box", lbl("consensus", bold=False)),
    "VP": (["study_variance_partition"], FULL, "box", lbl("variancePartition", bold=False)),
    "DECON": (["study_deconvolution"], FULL, "box", lbl("deconvolution", bold=False)),
    # ⑤ assemble + provenance
    "PUB": (["contrast_publication", "assemble_figure", "assemble_publication"], PUB, "box",
            lbl("publication → assemble_figure", "vector PDF + raster PNG · ×N sets")),
    "FRONT": (["front_door_artifacts"], IO, "box",
              lbl("front_door_artifacts", "promote figures/tables", bold=False)),
    "REPORT": (["report_html"], IO, "box", lbl("report_html", "REPORT.html")),
    "MANIFEST": (["release_manifest"], IO, "box",
                 lbl("release_manifest", "checksums · seeds · git rev")),
}

# Clusters, in flow order: (id, label, border, fill, members, rows).
# `rows` (optional) = (rank_rows, align): pin members into fixed rank rows so a
# wide block wraps into a compact grid instead of one long strip. align=True also
# column-aligns the rows with invisible edges (good for an edge-light grid like
# ④); align=False only pins ranks and lets the real edges route (good for ③,
# which has internal edges that would fight forced column alignment).
CLUSTERS = [
    ("input", "①  Input boundary  ·  choose one kind", "#9aa0a6", "#f4f5f6",
     ["BAM", "NFC", "ARC", "CNT", "MAT", "RES"],
     ([["BAM", "NFC"], ["ARC", "CNT"]], True)),  # stack the four kinds 2×2, not in one wide row
    ("study", "②  Study scope  ·  once per study", "#2C6E9B", "#eef5fb",
     ["QC", "BATCH"], None),
    ("contrast", "③  Per contrast  ·  effect = numerator − denominator", "#6b7280", "#f8f8f7",
     ["DE", "CONF", "OMNI", "PATH", "ONT", "SPIA", "EMAP", "COMP", "REG", "NET", "GRN", "HYP"],
     None),
    ("advanced", "④  Advanced  ·  full profile · flagged exploratory", "#B9770E", "#fdf6ec",
     ["SVA", "WGCNA", "CURV", "MED", "MULTI", "CONS", "VP", "DECON"],
     ([["SVA", "WGCNA", "CURV", "MED"], ["MULTI", "CONS", "VP", "DECON"]], True)),
    ("assemble", "⑤  Assemble + provenance", "#6b7280", "#f6f6f4",
     ["PUB", "FRONT", "REPORT", "MANIFEST"], None),
]

# Invisible edges that only constrain layout (no dependency meaning): seat the
# advanced block *below* the per-contrast block rather than beside it, so the
# DE→advanced→manifest path runs straight down instead of arcing across the canvas.
INVIS_EDGES = [("HYP", "SVA"), ("HYP", "MED")]

# Representative edges. (src, dst, style) — style "" solid, "dashed" conditional/exploratory.
# Compound edges into/out of the advanced cluster clip at its box (lhead/ltail).
EDGES = [
    ("BAM", "MAT", ""), ("NFC", "MAT", ""), ("ARC", "MAT", ""),
    ("CNT", "MAT", "dashed"),                       # count matrix skips featureCounts
    ("MAT", "QC", ""), ("QC", "BATCH", "dashed"),
    ("QC", "DE", ""), ("DE", "CONF", "dashed"), ("DE", "OMNI", "dashed"),
    ("DE", "PATH", ""), ("DE", "ONT", ""), ("DE", "COMP", ""),
    ("DE", "REG", ""), ("DE", "NET", ""), ("DE", "GRN", ""),
    ("RES", "PATH", ""), ("RES", "ONT", ""), ("RES", "REG", ""),
    ("RES", "NET", ""), ("RES", "GRN", ""),
    ("PATH", "SPIA", "dashed"), ("PATH", "EMAP", ""),
    ("DE", "HYP", ""), ("PATH", "HYP", ""), ("REG", "HYP", ""), ("NET", "HYP", ""),
    ("DE", "SVA", "dashed:lhead=cluster_advanced"),        # DE → advanced block
    ("DE", "PUB", ""), ("PATH", "PUB", ""), ("COMP", "PUB", ""), ("REG", "PUB", ""),
    ("NET", "PUB", ""), ("GRN", "PUB", ""), ("HYP", "PUB", ""),
    ("QC", "FRONT", ""), ("PUB", "FRONT", ""),
    ("FRONT", "REPORT", ""), ("FRONT", "MANIFEST", ""),
    ("VP", "MANIFEST", ":ltail=cluster_advanced"),         # advanced block → manifest
]


def declared_rules() -> set[str]:
    """Every `rule <name>:` in the workflow (the ground truth), minus the `all` target."""
    rules: set[str] = set()
    for path in RULES_DIR.rglob("*.smk"):
        rules |= set(re.findall(r"^\s*rule\s+([a-z0-9_]+)\s*:", path.read_text(), re.M))
    snakefile = RULES_DIR / "Snakefile"
    if snakefile.exists():
        rules |= set(re.findall(r"^\s*rule\s+([a-z0-9_]+)\s*:", snakefile.read_text(), re.M))
    rules.discard("all")
    return rules


def check_drift() -> None:
    """Abort unless every workflow rule is represented by exactly one node."""
    declared = declared_rules()
    if not declared:
        sys.exit(f"error: found no rule declarations under {RULES_DIR} — cannot verify the diagram")
    represented: dict[str, str] = {}
    dupes = []
    for node_id, (rules, *_rest) in NODES.items():
        for rule in rules:
            if rule in represented:
                dupes.append(f"{rule} (in {represented[rule]} and {node_id})")
            represented[rule] = node_id
    missing = declared - represented.keys()
    phantom = represented.keys() - declared
    problems = []
    if missing:
        problems.append("rules in the workflow with no node in the diagram (add them): "
                        + ", ".join(sorted(missing)))
    if phantom:
        problems.append("nodes reference rules that no longer exist (remove them): "
                        + ", ".join(sorted(phantom)))
    if dupes:
        problems.append("rules represented by more than one node: " + ", ".join(sorted(dupes)))
    if problems:
        sys.exit("pipeline diagram is out of sync with the workflow:\n  - " + "\n  - ".join(problems))


def build_dot() -> str:
    out = [
        "digraph tifzoret {",
        "  compound=true; rankdir=TB; bgcolor=transparent;",
        "  nodesep=0.28; ranksep=0.50; fontname=\"Helvetica\";",
        '  node[shape=box, style="rounded,filled", fontname="Helvetica", fontsize=10, '
        'penwidth=1.4, margin="0.14,0.08"];',
        '  edge[color="#8f8b80", arrowsize=0.7, penwidth=1.1];',
        "",
    ]
    for cid, clabel, border, fill, members, rows in CLUSTERS:
        out.append(f"  subgraph cluster_{cid} {{")
        out.append(f'    label=<<b>{clabel}</b>>; labeljust="l"; fontsize=12; fontcolor="{border}";')
        out.append(f'    style="rounded,filled,dashed"; color="{border}"; fillcolor="{fill}"; '
                   'penwidth=1.4; margin=12;')
        for nid in members:
            _rules, (bg, stroke, ink), shape, label = NODES[nid]
            out.append(f'    {nid}[label={label}, shape={shape}, fillcolor="{bg}", '
                       f'color="{stroke}", fontcolor="{ink}"];')
        if rows:
            rank_rows, align = rows
            for row in rank_rows:
                out.append("    {rank=same; " + " ".join(f"{n};" for n in row) + "}")
            if align:  # column-align consecutive rows with invisible edges
                for top, bot in zip(rank_rows, rank_rows[1:]):
                    for a, b in zip(top, bot):
                        out.append(f"    {a} -> {b} [style=invis];")
        out.append("  }")
        out.append("")
    for a, b in INVIS_EDGES:
        out.append(f"  {a} -> {b} [style=invis];")
    for src, dst, style in EDGES:
        attrs = []
        if style:
            parts = style.split(":")
            if parts[0] == "dashed":
                attrs.append("style=dashed")
            attrs += [p for p in parts[1:] if p]
        suffix = f" [{', '.join(attrs)}]" if attrs else ""
        out.append(f"  {src} -> {dst}{suffix};")
    out.append("}")
    return "\n".join(out) + "\n"


def main() -> None:
    if shutil.which("dot") is None:
        sys.exit("error: graphviz 'dot' not found on PATH")
    check_drift()
    dot = build_dot()
    subprocess.run(["dot", "-Tsvg", "-o", str(SVG_OUT)], input=dot, text=True, check=True)
    subprocess.run(["dot", "-Tpng", "-Gdpi=150", "-o", str(PNG_OUT)], input=dot, text=True, check=True)
    print(f"wrote {SVG_OUT.relative_to(REPO)} and {PNG_OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
