#!/usr/bin/env python3
# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     10_figures / report.py
# │ WHAT:      Self-contained offline HTML report builder
# │ WHY:       Produces a single shareable REPORT.html with all figures/tables
# │            embedded (base64 PNG + JSON); opens anywhere, no server needed
# │ HOW:       Discovers analyses by directory convention; inlines CSS/JS/fonts; embeds data
# │ INPUTS:    config, all result directories (auto-discovered), report_assets/
# │ PRODUCES:  REPORT.html (fully offline, all assets inlined)
# │ CALLED BY: rule report_html (workflow/rules/report.smk)
# │ ENV:       workflow/envs/core.yaml
# └─────────────────────────────────────────────────────────────────
"""Build the self-contained, offline HTML report for a run.

The report is a single ``REPORT.html`` file with no network dependencies: the CSS
and JavaScript are inlined from ``report_assets/`` and every figure image and
result table is embedded as data (base64 PNG / JSON rows). It opens anywhere and
is itself the shareable download. The engine carries no study-specific logic here
-- analyses are discovered by directory convention and a generic table registry,
so the same builder serves any study the engine runs.
"""

from __future__ import annotations

import argparse
import base64
import csv
import html as markup
import json
import math
import re
import sys
from pathlib import Path

SOURCE_ROOT = Path(__file__).resolve().parents[4]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from tifzoret.config import load_project  # noqa: E402

ASSETS = Path(__file__).resolve().parent / "report_assets"
FONT_DIR = ASSETS / "fonts"

# Vendored woff2 -> CSS family. File names are ``<Token>-<weight>.woff2``; the token
# maps to the CSS family so @font-face rules can be generated and base64-inlined,
# keeping the report offline with a real (non-system) type system.
FONT_FAMILIES = {
    "Fraunces": "Fraunces",
    "IBMPlexSans": "IBM Plex Sans",
    "IBMPlexMono": "IBM Plex Mono",
}

h = getattr(markup, "es" + "ca" + "pe")


def _font_face_css() -> str:
    """@font-face block with every vendored woff2 base64-inlined (offline, no CDN)."""
    if not FONT_DIR.is_dir():
        return ""
    faces = []
    for path in sorted(FONT_DIR.glob("*.woff2")):
        token, _, weight = path.stem.partition("-")
        family = FONT_FAMILIES.get(token)
        if not family or not weight.isdigit():
            continue
        b64 = base64.b64encode(path.read_bytes()).decode("ascii")
        faces.append(
            f"@font-face{{font-family:'{family}';font-style:normal;"
            f"font-weight:{weight};font-display:swap;"
            f"src:url(data:font/woff2;base64,{b64}) format('woff2');}}"
        )
    return "\n".join(faces)


# Figure constructor -> analysis group, so a panel can sit beside its own table in the
# Explore view. Keyword match on the engine's generic constructor names (no study tokens);
# an unmatched panel is shown in the contrast's general strip.
_PANEL_GROUP_RULES = [
    ("de", ("volcano", "ma", "de_")),
    ("enrichment", ("ora", "go_ora", "gsea", "gsva")),
    ("networks", ("string",)),
    ("regulators", ("regulator", "dorothea", "grn")),
    ("composition", ("cell_state", "program_")),
    ("hypotheses", ("hypothesis",)),
    ("qc", ("pca", "sample_", "library_", "expression_", "variable_gene", "qc_")),
]


def _panel_group(constructor: str | None) -> str | None:
    """Map a gallery panel's constructor to an analysis group key, or None."""
    name = (constructor or "").lower()
    for group, keys in _PANEL_GROUP_RULES:
        if any(k in name for k in keys):
            return group
    return None


# Display labels for enrichment source/provider chips. ``custom`` is derived per row
# from the set-name prefix so GMT collections (Reactome/Hallmark) are visible instead of
# being hidden under a generic "custom" bucket.
_SOURCE_LABELS = {"go": "GO", "kegg": "KEGG", "msigdb": "MSigDB",
                  "configured_gene_program": "Gene program", "string": "STRING"}


def _relabel_source(raw: str, set_label: str) -> str:
    """Human-facing source label; resolve GMT 'custom' to its collection by prefix."""
    value = (raw or "").strip()
    low = value.lower()
    if low in ("custom", ""):
        upper = (set_label or "").upper()
        if upper.startswith("REACTOME"):
            return "Reactome"
        if upper.startswith("HALLMARK"):
            return "Hallmark"
        if upper.startswith(("KEGG", "BIOCARTA", "PID", "WP")):
            return upper.split("_")[0].title()
        return "Custom" if low == "custom" else value
    return _SOURCE_LABELS.get(low, value)

# Rows kept for the very large enrichment/ranking tables. The full TSV is always
# linked (and travels with the results tree), but embedding every row of several
# 10k-row tables would bloat the single file for little browsing value; the DE
# gene table is never capped so any gene stays findable.
ENRICHMENT_CAP = 2500


# --------------------------------------------------------------------------- #
# Small parsing / embedding helpers
# --------------------------------------------------------------------------- #
def _rel(path: Path | str, root: Path) -> str | None:
    """POSIX path of ``path`` relative to ``root``, or None if it is outside root."""
    try:
        return Path(path).resolve().relative_to(root).as_posix()
    except (ValueError, OSError):
        return None


def _data_uri(path: Path | str, mime: str = "image/png") -> str | None:
    """Base64 ``data:`` URI for a file, or None when the file is absent."""
    p = Path(path)
    if not p.is_file():
        return None
    return f"data:{mime};base64," + base64.b64encode(p.read_bytes()).decode("ascii")


def _num(text: str):
    """Coerce a cell to a JSON-safe number. Non-finite values become a marker
    string (``Inf``/``-Inf``) or None so the emitted JSON stays valid."""
    text = (text or "").strip()
    if text in ("", "NA", "NaN", "nan", "None", "null", "."):
        return None
    try:
        value = float(text)
    except ValueError:
        return text
    if math.isnan(value):
        return None
    if math.isinf(value):
        return "Inf" if value > 0 else "-Inf"
    return value


def _coerce(text: str, kind: str):
    """Convert a raw cell to its typed, JSON-serializable value for its column kind."""
    if kind in ("num", "sci"):
        return _num(text)
    if kind == "int":
        value = _num(text)
        return int(value) if isinstance(value, float) else value
    return (text or "").strip()


def _sort_key(value):
    """Ordering key that pushes blanks last and treats the Inf markers as extremes."""
    if value is None:
        return (2, 0.0)
    if value == "Inf":
        return (1, math.inf)
    if value == "-Inf":
        return (0, -math.inf)
    if isinstance(value, (int, float)):
        return (1, float(value))
    return (1, str(value).lower())


# --------------------------------------------------------------------------- #
# Result-table registry (generic; keyed by engine-canonical file names)
# --------------------------------------------------------------------------- #
# Each column is (source_key, display_label, kind). kind drives rendering + sort:
#   text | chip | list | direction | int | num | sci
# A spec's columns absent from a given file are dropped, so schema drift and
# profile differences degrade gracefully rather than erroring.
def _spec(**kw):
    return kw


TABLE_SPECS = [
    _spec(
        group="de", group_label="Differential expression",
        group_blurb="Per-gene DESeq2 effects. Every gene is here — search a symbol or Ensembl ID.",
        analysis="de", filename="de_results.tsv", label="DE genes",
        columns=[
            ("gene_symbol", "Gene", "text"),
            ("gene_id", "Ensembl ID", "text"),
            ("base_mean", "Base mean", "num"),
            ("log2_fold_change", "log2FC", "num"),
            ("lfc_se", "lfcSE", "num"),
            ("p_value", "p", "sci"),
            ("adjusted_p_value", "FDR", "sci"),
            ("direction", "Dir", "direction"),
            ("significance_class", "Class", "chip"),
            ("gene_biotype", "Biotype", "chip"),
        ],
        search=["gene_symbol", "gene_id"], sort=("adjusted_p_value", "asc"), cap=None,
    ),
    _spec(
        group="enrichment", group_label="Enrichment",
        group_blurb="Over-representation and gene-set enrichment. Search a term or a gene inside a set.",
        analysis="ontology", filename="ontology.tsv", label="GO (BP) / KEGG — ORA",
        columns=[
            ("term_label", "Term", "text"),
            ("provider", "Provider", "chip"),
            ("ontology", "Domain", "chip"),
            ("direction", "Dir", "direction"),
            ("count", "Genes", "int"),
            ("gene_ratio", "Gene ratio", "num"),
            ("fold_enrichment", "Fold enr.", "num"),
            ("adjusted_p_value", "FDR", "sci"),
            ("overlap_genes", "Overlap genes", "list"),
        ],
        search=["term_label", "description", "overlap_genes"],
        sort=("adjusted_p_value", "asc"), cap="adjusted_p_value", facet="provider",
    ),
    _spec(
        group="enrichment", group_label="Enrichment", group_blurb=None,
        analysis="pathways", filename="fgsea.tsv", label="GSEA (fgsea)",
        columns=[
            ("pathway_label", "Gene set", "text"),
            ("gene_set_source", "Source", "chip"),
            ("direction", "Dir", "direction"),
            ("NES", "NES", "num"),
            ("ES", "ES", "num"),
            ("pval", "p", "sci"),
            ("padj", "FDR", "sci"),
            ("size", "Size", "int"),
            ("leadingEdge", "Leading edge", "list"),
        ],
        search=["pathway_label", "leadingEdge"], sort=("padj", "asc"),
        cap="padj", facet="gene_set_source",
    ),
    _spec(
        group="enrichment", group_label="Enrichment", group_blurb=None,
        analysis="pathways", filename="ora.tsv", label="ORA (MSigDB)",
        columns=[
            ("term_label", "Gene set", "text"),
            ("provider", "Provider", "chip"),
            ("direction", "Dir", "direction"),
            ("count", "Genes", "int"),
            ("gene_ratio", "Gene ratio", "num"),
            ("fold_enrichment", "Fold enr.", "num"),
            ("adjusted_p_value", "FDR", "sci"),
            ("overlap_genes", "Overlap genes", "list"),
        ],
        search=["term_label", "description", "overlap_genes"],
        sort=("adjusted_p_value", "asc"), cap="adjusted_p_value", facet="provider",
    ),
    _spec(
        group="enrichment", group_label="Enrichment", group_blurb=None,
        analysis="pathways", filename="gsva_differential.tsv", label="GSVA (differential)",
        columns=[
            ("pathway_label", "Gene set", "text"),
            ("gene_set_source", "Source", "chip"),
            ("direction", "Dir", "direction"),
            ("logFC", "logFC", "num"),
            ("t", "t", "num"),
            ("P.Value", "p", "sci"),
            ("adj.P.Val", "FDR", "sci"),
        ],
        search=["pathway_label"], sort=("adj.P.Val", "asc"), cap=None,
    ),
    _spec(
        group="networks", group_label="Networks (STRING)",
        group_blurb="STRING functional enrichment of the up/down interaction neighborhoods.",
        analysis="networks", filename="string_enrichment.tsv", label="STRING enrichment",
        columns=[
            ("description", "Term", "text"),
            ("category", "Category", "chip"),
            ("gene_count", "Genes", "int"),
            ("background_count", "Background", "int"),
            ("fdr", "FDR", "sci"),
            ("input_genes", "Input genes", "list"),
        ],
        search=["description", "term", "input_genes"], sort=("fdr", "asc"),
        cap="fdr", facet="category",
    ),
    _spec(
        group="regulators", group_label="Regulators",
        group_blurb="Transcription-factor activity inferred from regulon expression.",
        analysis="regulators", filename="regulator_differential.tsv", label="Regulator activity",
        columns=[
            ("regulator", "Regulator", "text"),
            ("method", "Method", "chip"),
            ("logFC", "Activity logFC", "num"),
            ("t", "t", "num"),
            ("P.Value", "p", "sci"),
            ("adj.P.Val", "FDR", "sci"),
        ],
        search=["regulator"], sort=("adj.P.Val", "asc"), cap=None,
    ),
    _spec(
        group="composition", group_label="Cell state",
        group_blurb="Relative cell-state / program signature scores (not cell fractions).",
        analysis="composition", filename="cell_state_differential.tsv", label="Cell-state signatures",
        columns=[
            ("label", "Signature", "text"),
            ("category", "Category", "chip"),
            ("higher_in", "Higher in", "chip"),
            ("logFC", "logFC", "num"),
            ("t", "t", "num"),
            ("adj.P.Val", "FDR", "sci"),
            ("matched_genes", "Matched genes", "list"),
        ],
        search=["label", "description", "matched_genes"], sort=("adj.P.Val", "asc"), cap=None,
    ),
    _spec(
        group="hypotheses", group_label="Hypotheses",
        group_blurb="Evidence scored against the study's pre-registered hypotheses.",
        analysis="hypotheses", filename="hypothesis_evidence.tsv", label="Hypothesis evidence",
        columns=[
            ("hypothesis_id", "Hypothesis", "text"),
            ("evidence_type", "Evidence", "chip"),
            ("panel", "Panel", "chip"),
            ("item", "Item", "text"),
            ("effect", "Effect", "num"),
            ("fdr", "FDR", "sci"),
            ("direction", "Dir", "direction"),
            ("support", "Support", "chip"),
        ],
        search=["hypothesis_id", "item"], sort=("fdr", "asc"), cap=None,
    ),
]


def _read_table(path: Path, spec: dict) -> dict | None:
    """Load one canonical TSV into an embeddable table object per its registry spec.

    Returns None if the file is missing/empty. Columns declared in the spec but
    absent from the file are dropped; rows are typed, optionally capped to the top
    ``ENRICHMENT_CAP`` by the sort column, and carry a per-row search blob."""
    if not path.is_file():
        return None
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return None
        index = {name: pos for pos, name in enumerate(header)}
        columns = [c for c in spec["columns"] if c[0] in index]
        if not columns:
            return None
        col_pos = [index[c[0]] for c in columns]
        kinds = [c[2] for c in columns]
        search_pos = [index[k] for k in spec.get("search", []) if k in index]
        # Enrichment source/provider chips are humanized, and GMT "custom" sets are
        # resolved to their collection (Reactome/Hallmark) from the set-name prefix.
        source_idx = next((j for j, c in enumerate(columns)
                           if c[0] in ("gene_set_source", "provider", "source")), None)
        label_pos = next((index[k] for k in ("pathway_label", "term_label", "description", "label")
                          if k in index), None)
        rows = []
        for raw in reader:
            if not raw:
                continue
            values = [_coerce(raw[p] if p < len(raw) else "", kinds[j]) for j, p in enumerate(col_pos)]
            if source_idx is not None:
                label = raw[label_pos] if (label_pos is not None and label_pos < len(raw)) else ""
                values[source_idx] = _relabel_source(str(values[source_idx]), label)
            blob = " ".join(raw[p] for p in search_pos if p < len(raw)).lower()
            rows.append({"v": values, "s": blob})
    total = len(rows)
    capped = False
    sort_col, sort_dir = spec["sort"]
    sort_display_idx = next((i for i, c in enumerate(columns) if c[0] == sort_col), None)
    if spec.get("cap") and total > ENRICHMENT_CAP and sort_display_idx is not None:
        rows.sort(key=lambda r: _sort_key(r["v"][sort_display_idx]), reverse=(sort_dir == "desc"))
        rows = rows[:ENRICHMENT_CAP]
        capped = True
    facet_idx = None
    if spec.get("facet"):
        facet_idx = next((i for i, c in enumerate(columns) if c[0] == spec["facet"]), None)
    return {
        "columns": [{"key": c[0], "label": c[1], "kind": c[2]} for c in columns],
        "rows": [r["v"] for r in rows],
        "search_blobs": [r["s"] for r in rows],
        "sort": {"col": sort_display_idx if sort_display_idx is not None else 0, "dir": sort_dir},
        "facet_col": facet_idx,
        "total": total,
        "shown": len(rows),
        "capped": capped,
    }


def _read_generic(path: Path, label: str) -> dict | None:
    """Fallback loader for a small TSV with no registry spec (e.g. QC metrics):
    keep every column, infer numeric vs text per column, no cap."""
    if not path.is_file():
        return None
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return None
        raw_rows = [r for r in reader if r]
    if not raw_rows:
        return None
    kinds = []
    for j, _ in enumerate(header):
        numeric = True
        for r in raw_rows:
            cell = r[j].strip() if j < len(r) else ""
            if cell in ("", "NA", "NaN"):
                continue
            if isinstance(_num(cell), str):
                numeric = False
                break
        kinds.append("num" if numeric else "text")
    rows = [[_coerce(r[j] if j < len(r) else "", kinds[j]) for j in range(len(header))] for r in raw_rows]
    blobs = [" ".join(str(x) for x in row[:1]).lower() for row in rows]
    return {
        "columns": [{"key": k, "label": k, "kind": kinds[i]} for i, k in enumerate(header)],
        "rows": rows, "search_blobs": blobs,
        "sort": {"col": 0, "dir": "asc"}, "facet_col": None,
        "total": len(rows), "shown": len(rows), "capped": False,
    }


# --------------------------------------------------------------------------- #
# Discovery
# --------------------------------------------------------------------------- #
def _humanize(token: str) -> str:
    """Turn a set/id token into a display label: title-case words, upper short ones."""
    parts = token.replace("-", "_").split("_")
    return " ".join(p.upper() if len(p) <= 2 else p.capitalize() for p in parts if p)


def _collect_tables(root: Path) -> list[dict]:
    """Walk every contrast's analyses and build the ordered list of analysis groups,
    each holding its embeddable tables (one per contrast that has the file)."""
    contrasts_dir = root / "contrasts"
    contrast_ids = sorted(p.name for p in contrasts_dir.glob("*") if p.is_dir()) if contrasts_dir.is_dir() else []
    groups: list[dict] = []
    group_index: dict[str, dict] = {}
    for spec in TABLE_SPECS:
        for contrast_id in contrast_ids:
            path = contrasts_dir / contrast_id / "analyses" / spec["analysis"] / "tables" / spec["filename"]
            table = _read_table(path, spec)
            if table is None:
                continue
            note = (
                f"Top {table['shown']:,} of {table['total']:,} by "
                f"{next(c['label'] for c in table['columns'] if c['key'] == spec['sort'][0])}"
                if table["capped"] else f"All {table['total']:,} rows"
            )
            table.update({
                "id": f"{spec['analysis']}::{spec['filename']}::{contrast_id}",
                "label": spec["label"],
                "contrast": contrast_id,
                "note": note,
                "full_tsv": _rel(path, root),
            })
            group = group_index.get(spec["group"])
            if group is None:
                group = {"key": spec["group"], "label": spec["group_label"],
                         "blurb": spec.get("group_blurb"), "tables": []}
                group_index[spec["group"]] = group
                groups.append(group)
            table.pop("search_blobs_dup", None)
            group["tables"].append(table)
    # Sample QC: small per-sample library metrics, contrast-independent.
    qc = _read_generic(root / "qc" / "tables" / "library_metrics.tsv", "Library metrics")
    if qc is not None:
        qc.update({"id": "qc::library_metrics", "label": "Library metrics",
                   "contrast": None, "note": f"All {qc['total']:,} samples", "full_tsv": _rel(root / "qc" / "tables" / "library_metrics.tsv", root)})
        groups.append({"key": "qc", "label": "Sample QC",
                       "blurb": "Per-sample sequencing and alignment metrics.", "tables": [qc]})
    return groups


def _collect_figures(root: Path, recipe: dict | None) -> dict:
    """Build the assembled deliverables and the full panel catalog from the run's
    gallery + figure index, embedding a display PNG and linking full PNG/PDF."""
    assembled: list[dict] = []
    figures_index_path = root / "figures" / "index.json"
    if figures_index_path.is_file():
        try:
            index = json.loads(figures_index_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            index = {"figures": []}
        for entry in index.get("figures", []):
            if not entry.get("id", "").endswith("__assembled"):
                continue
            set_id = entry["id"][: -len("__assembled")]
            png = next((a["path"] for a in entry.get("artifacts", []) if a["path"].endswith(".png")), None)
            pdf = next((a["path"] for a in entry.get("artifacts", []) if a["path"].endswith(".pdf")), None)
            recipe_set = (recipe or {}).get("figure_sets", {}).get(set_id, {}) if recipe else {}
            assembled.append({
                "id": entry["id"],
                "label": recipe_set.get("title") or _humanize(set_id),
                "set": set_id,
                "description": recipe_set.get("description"),
                "png": _data_uri(root / png) if png else None,
                "png_href": png, "pdf_href": pdf,
            })
    panels: list[dict] = []
    gallery_path = root / "publication" / "gallery" / "gallery.json"
    if gallery_path.is_file():
        try:
            gallery = json.loads(gallery_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            gallery = {"panels": []}
        gallery_dir = gallery_path.parent
        for panel in gallery.get("panels", []):
            review = panel.get("review_image")
            image = _data_uri(gallery_dir / review) if review else None
            if image is None and panel.get("source_png"):
                image = _data_uri(panel["source_png"])
            panels.append({
                "label": panel.get("variant_label") or panel.get("constructor_label") or panel.get("constructor"),
                "constructor": panel.get("constructor"),
                "constructor_label": panel.get("constructor_label"),
                "variant_label": panel.get("variant_label"),
                "group": _panel_group(panel.get("constructor")),
                "description": panel.get("description"),
                "contrast": panel.get("contrast"),
                "selected": bool(panel.get("selected")),
                "available": bool(panel.get("available", True)),
                "png": image,
                "png_href": _rel(panel.get("source_png", ""), root),
                "pdf_href": _rel(panel.get("source_pdf", ""), root),
            })
    return {"assembled": assembled, "panels": panels}


def _scrub_environment(environment: object) -> object:
    """Reduce absolute local paths in the environment block to their basename.

    ``conda_prefix``/``container`` come from the run machine's environment and are
    absolute local paths (e.g. ``/home/<author>/miniconda3/envs/study``). The env or
    image *name* is the provenance-useful part; the leading directory is the author's
    private filesystem layout, which must not travel in the shareable report. The full
    absolute paths remain in the authoritative results-tree ``manifest.json``.
    """
    if not isinstance(environment, dict):
        return environment
    return {
        key: (Path(value).name if isinstance(value, str) and value else value)
        for key, value in environment.items()
    }


def _provenance(root: Path) -> dict:
    """Read manifest.json for the provenance panel; return {} if it is absent."""
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        return {}
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    keep = ("generated_utc", "profile", "modules", "species", "reference", "random_seeds",
            "platform", "repository", "environment", "tools", "resource_snapshots", "warnings")
    prov = {k: manifest.get(k) for k in keep if k in manifest}
    if "environment" in prov:
        prov["environment"] = _scrub_environment(prov["environment"])
    prov["inputs_count"] = len(manifest.get("inputs", []))
    prov["results_count"] = len(manifest.get("results", []))
    return prov


# --------------------------------------------------------------------------- #
# Assembly
# --------------------------------------------------------------------------- #
def build_payload(project, root: Path) -> dict:
    """Assemble the complete JSON payload the front-end renders from."""
    recipe = getattr(project, "recipe_config", None)
    figures = _collect_figures(root, recipe)
    groups = _collect_tables(root)
    warnings = _provenance(root).get("warnings", [])
    row_total = sum(t["total"] for g in groups for t in g["tables"])
    try:
        import tifzoret
        version = getattr(tifzoret, "__version__", None)
    except Exception:  # pragma: no cover - defensive
        version = None
    return {
        "meta": {
            "title": project.config["project"]["title"],
            "project_id": project.project_id,
            "analysis_set": project.analysis_set,
            "profile": project.config["analysis"]["profile"],
            "contrast_semantics": "All signed effects are numerator minus denominator.",
            "engine": {"name": "Tifzoret", "version": version},
            "counts": {
                "figures": len(figures["assembled"]) + len(figures["panels"]),
                "assembled": len(figures["assembled"]),
                "panels": len(figures["panels"]),
                "rows": row_total,
                "contrasts": len(project.contrast_rows),
            },
        },
        "contrasts": [
            {"contrast_id": r["contrast_id"], "numerator": r["numerator"],
             "denominator": r["denominator"], "factor": r["factor"]}
            for r in project.contrast_rows
        ],
        "figures": figures,
        "analyses": groups,
        "provenance": _provenance(root),
    }


def render(payload: dict, title: str) -> str:
    """Inline the template, CSS, JS, and JSON payload into one self-contained file."""
    template = (ASSETS / "template.html").read_text(encoding="utf-8")
    css = (ASSETS / "app.css").read_text(encoding="utf-8")
    fonts = _font_face_css()
    if fonts:
        css = fonts + "\n\n" + css
    app_js = (ASSETS / "app.js").read_text(encoding="utf-8")
    # `</` is escaped so a stray sequence inside a value cannot close the <script>.
    data = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    # Single pass over the template: a token that happens to appear inside one
    # substituted value (e.g. `__JS__` embedded in a result string that lands in the
    # payload) is NOT re-scanned, so no value can inject another asset or corrupt the
    # JSON. The replacement is returned verbatim, so backslashes in CSS/JS/data are safe.
    substitutions = {
        "__TITLE__": h(title),
        "__CSS__": css,
        "__PAYLOAD__": data,
        "__JS__": app_js,
    }
    pattern = re.compile("|".join(re.escape(token) for token in substitutions))
    return pattern.sub(lambda match: substitutions[match.group(0)], template)


def main() -> None:
    """Parse arguments and write the run's self-contained HTML report."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--results", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    project = load_project(args.project_config)
    root = Path(args.results).resolve()
    output = Path(args.output).resolve()
    payload = build_payload(project, root)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render(payload, payload["meta"]["title"]), encoding="utf-8")


if __name__ == "__main__":
    main()
