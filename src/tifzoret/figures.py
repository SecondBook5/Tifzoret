"""Registered publication panels and hypothesis-driven review tooling."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable

import yaml
from PIL import Image, ImageDraw, ImageFont

if TYPE_CHECKING:
    from .config import ResolvedProject


@dataclass(frozen=True)
class PanelVariant:
    """One renderable form of a panel: its result-relative source figure stem,
    the displayed-data tables it exposes, the analysis module that must be
    enabled to produce it, and a human-readable label."""

    source: str
    displayed_data: tuple[str, ...]
    required_module: str
    label: str


@dataclass(frozen=True)
class PanelConstructor:
    """A registered publication-panel type in the catalog: its identifier,
    label, and description, whether it is contrast-specific, its default
    variant, and the named variants that render it."""

    id: str
    label: str
    description: str
    contrast_specific: bool
    default_variant: str
    variants: dict[str, PanelVariant]


@dataclass(frozen=True)
class ResolvedPanel:
    """A recipe panel resolved to concrete paths: the chosen constructor,
    variant, and contrast, its label, the source figure path, the displayed-data
    file paths, and the module the workflow must have run to build it."""

    constructor: str | None
    variant: str | None
    contrast: str | None
    label: str
    source: Path
    displayed_data: tuple[Path, ...]
    required_module: str | None


def _variant(source: str, module: str, label: str, *displayed: str) -> PanelVariant:
    return PanelVariant(source, tuple(displayed), module, label)


PANEL_REGISTRY: dict[str, PanelConstructor] = {
    "pca": PanelConstructor(
        "pca", "PCA", "Sample PCA with group-colored covariance ellipses.", False, "default",
        {"default": _variant("qc/figures/pca", "qc", "PCA", "qc/tables/pca_coordinates.tsv", "qc/tables/pca_variance.tsv")},
    ),
    "sample_correlation": PanelConstructor(
        "sample_correlation", "Sample correlation", "Clustered Pearson sample-correlation heatmap.", False, "default",
        {"default": _variant("qc/figures/sample_correlation", "qc", "Sample correlation", "qc/tables/sample_correlation.tsv")},
    ),
    "pca_correlation": PanelConstructor(
        "pca_correlation", "PCA and sample correlation", "Shared-legend PCA and clustered correlation layout.", False, "default",
        {"default": _variant("qc/figures/pca_correlation", "qc", "PCA + correlation", "qc/tables/pca_coordinates.tsv", "qc/tables/pca_variance.tsv", "qc/tables/sample_correlation.tsv", "qc/tables/pca_correlation_layout.json")},
    ),
    "library_metrics": PanelConstructor(
        "library_metrics", "Library metrics", "Per-sample library size and detected-gene sequencing-depth metrics.", False, "default",
        {"default": _variant("qc/figures/library_metrics", "qc", "Library metrics", "qc/tables/library_metrics.tsv")},
    ),
    "detected_genes": PanelConstructor(
        "detected_genes", "Detected genes", "Per-sample detected-gene counts as dots on a data-focused axis with a median reference line.", False, "default",
        {"default": _variant("qc/figures/detected_genes", "qc", "Detected genes", "qc/tables/detected_genes_displayed.tsv")},
    ),
    "expression_density": PanelConstructor(
        "expression_density", "Expression density", "Per-sample log-expression density distributions.", False, "default",
        {"default": _variant("qc/figures/expression_density", "qc", "Expression density", "qc/tables/expression_density_displayed.tsv")},
    ),
    "sample_distance": PanelConstructor(
        "sample_distance", "Sample distance", "Euclidean sample-to-sample distance heatmap on variance-stabilized expression.", False, "default",
        {"default": _variant("qc/figures/sample_distance", "qc", "Sample distance", "qc/tables/sample_distance.tsv")},
    ),
    "variable_gene_heatmap": PanelConstructor(
        "variable_gene_heatmap", "Most-variable-gene heatmap", "Row-scaled expression of the most variable genes across samples.", False, "default",
        {"default": _variant("qc/figures/variable_gene_heatmap", "qc", "Variable-gene heatmap", "qc/tables/variable_gene_heatmap_displayed.tsv")},
    ),
    "sample_outliers": PanelConstructor(
        "sample_outliers", "Sample outlier screen", "PCA-distance outlier screen: per-sample Euclidean distance in the leading PCs against a robust median + MAD threshold (flagged, never removed).", False, "default",
        {"default": _variant("qc/figures/sample_outliers", "qc", "Outlier screen", "qc/tables/outlier_distances.tsv")},
    ),
    "qc_overview": PanelConstructor(
        "qc_overview", "QC overview", "Composite quality-control overview across library, density, PCA, and correlation views.", False, "default",
        {"default": _variant("qc/figures/qc_overview", "qc", "QC overview")},
    ),
    "volcano": PanelConstructor(
        "volcano", "Volcano plot", "Shrunken effect sizes with configured thresholds and labels.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/de/figures/volcano", "de", "Volcano", "contrasts/{contrast}/analyses/de/tables/volcano_displayed.tsv")},
    ),
    "ma": PanelConstructor(
        "ma", "MA plot", "Mean abundance versus shrunken effect size.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/de/figures/ma", "de", "MA", "contrasts/{contrast}/analyses/de/tables/ma_displayed.tsv")},
    ),
    "de_overview": PanelConstructor(
        "de_overview", "Differential-expression overview", "Multi-view differential-expression summary.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/de/figures/de_overview", "de", "DE overview", "contrasts/{contrast}/analyses/de/tables/volcano_displayed.tsv", "contrasts/{contrast}/analyses/de/tables/ma_displayed.tsv")},
    ),
    "de_heatmap": PanelConstructor(
        "de_heatmap", "Differential-expression heatmap", "Top DE heatmap with selectable biological organization.", True, "global_clustered",
        {
            "default": _variant("contrasts/{contrast}/analyses/de/figures/de_heatmap", "de", "DE heatmap", "contrasts/{contrast}/analyses/de/tables/de_heatmap_displayed.tsv"),
            "global_clustered": _variant("contrasts/{contrast}/analyses/publication/figures/de_heatmap_global", "publication", "Global clustering", "contrasts/{contrast}/analyses/publication/tables/de_heatmap_global_displayed.tsv", "contrasts/{contrast}/analyses/publication/tables/de_gene_program_assignments.tsv", "contrasts/{contrast}/analyses/publication/tables/program_definitions.tsv"),
            "program_grouped": _variant("contrasts/{contrast}/analyses/publication/figures/de_heatmap_program_grouped", "publication", "Program-grouped clustering", "contrasts/{contrast}/analyses/publication/tables/de_heatmap_program_grouped_displayed.tsv", "contrasts/{contrast}/analyses/publication/tables/de_gene_program_assignments.tsv", "contrasts/{contrast}/analyses/publication/tables/program_definitions.tsv"),
            "direct_program_labels": _variant("contrasts/{contrast}/analyses/publication/figures/de_heatmap_compact", "publication", "Direct program labels", "contrasts/{contrast}/analyses/publication/tables/de_heatmap_compact_displayed.tsv", "contrasts/{contrast}/analyses/publication/tables/de_gene_program_assignments.tsv", "contrasts/{contrast}/analyses/publication/tables/program_definitions.tsv"),
        },
    ),
    "cell_state_effects": PanelConstructor(
        "cell_state_effects", "Cell-state signature shifts", "Grouped relative signature effects with matched-gene and FDR annotations.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/composition/figures/cell_state_signatures", "composition", "Cell-state signatures", "contrasts/{contrast}/analyses/composition/tables/cell_state_displayed.tsv")},
    ),
    "ora_bidirectional": PanelConstructor(
        "ora_bidirectional", "Bidirectional ORA", "Up- and downregulated terms in one gradient bubble plot.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/pathways/figures/ora_bidirectional", "pathways", "Bidirectional ORA", "contrasts/{contrast}/analyses/pathways/tables/ora_displayed.tsv")},
    ),
    "go_ora": PanelConstructor(
        "go_ora", "GO biological-process ORA", "Combined directional GO Biological Process bubble plot.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/ontology/figures/ontology_bidirectional", "ontology", "GO BP ORA", "contrasts/{contrast}/analyses/ontology/tables/ontology_displayed.tsv")},
    ),
    "gsva_heatmap": PanelConstructor(
        "gsva_heatmap", "GSVA heatmap", "Configured pathway scores across samples.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/pathways/figures/gsva_heatmap", "pathways", "GSVA heatmap", "contrasts/{contrast}/analyses/pathways/tables/gsva_heatmap_displayed.tsv", "contrasts/{contrast}/analyses/pathways/tables/gsva_differential.tsv")},
    ),
    "gsea_multitrack": PanelConstructor(
        "gsea_multitrack", "Advanced GSEA curves", "Enrichment score, hits, ranked metric, NES, FDR, and leading edge.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/pathways/figures/gsea_curves", "pathways", "Advanced GSEA", "contrasts/{contrast}/analyses/pathways/tables/gsea_curves_displayed.tsv", "contrasts/{contrast}/analyses/pathways/tables/fgsea.tsv")},
    ),
    "program_heatmap_effects": PanelConstructor(
        "program_heatmap_effects", "Program heatmap and effects", "Consolidated expression heatmap and gene-level effect estimates.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/publication/figures/program_integrated", "publication", "Program heatmap + effects", "contrasts/{contrast}/analyses/publication/tables/program_integrated_displayed.tsv", "contrasts/{contrast}/analyses/publication/tables/program_definitions.tsv")},
    ),
    "program_violins": PanelConstructor(
        "program_violins", "Program gene distributions", "Consolidated program-shaded violins with adjusted significance brackets.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/publication/figures/program_violins", "publication", "Program violins", "contrasts/{contrast}/analyses/publication/tables/program_violins_displayed.tsv", "contrasts/{contrast}/analyses/publication/tables/program_violins_tests.tsv", "contrasts/{contrast}/analyses/publication/tables/program_definitions.tsv")},
    ),
    "string_enrichment": PanelConstructor(
        "string_enrichment", "STRING functional enrichment", "Three-facet directional STRING enrichment bubble (down / GSEA leading edge / up).", True, "faceted",
        {
            "faceted": _variant("contrasts/{contrast}/analyses/networks/figures/string_enrichment_faceted", "networks", "STRING enrichment", "contrasts/{contrast}/analyses/networks/tables/string_down_enrichment.tsv", "contrasts/{contrast}/analyses/networks/tables/string_leading_edge_enrichment.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_enrichment.tsv"),
            "combined": _variant("contrasts/{contrast}/analyses/networks/figures/string_enrichment", "networks", "STRING enrichment (combined)", "contrasts/{contrast}/analyses/networks/tables/string_enrichment_displayed.tsv"),
        },
    ),
    "string_network": PanelConstructor(
        "string_network", "STRING network", "Direction-specific STRING association network with topology-derived community hulls; auditable full and displayed edges.", True, "upregulated",
        {
            "upregulated": _variant("contrasts/{contrast}/analyses/networks/figures/string_up_network_community", "networks", "Upregulated STRING community network", "contrasts/{contrast}/analyses/networks/tables/string_up_nodes.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_edges.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_input_genes.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_unmapped_genes.tsv"),
            "downregulated": _variant("contrasts/{contrast}/analyses/networks/figures/string_down_network_community", "networks", "Downregulated STRING community network", "contrasts/{contrast}/analyses/networks/tables/string_down_nodes.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_edges.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_input_genes.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_unmapped_genes.tsv"),
            "upregulated_legacy": _variant("contrasts/{contrast}/analyses/networks/figures/string_up_network", "networks", "Upregulated STRING network (matplotlib)", "contrasts/{contrast}/analyses/networks/tables/string_up_nodes_displayed.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_edges_displayed.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_input_genes.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_unmapped_genes.tsv"),
            "downregulated_legacy": _variant("contrasts/{contrast}/analyses/networks/figures/string_down_network", "networks", "Downregulated STRING network (matplotlib)", "contrasts/{contrast}/analyses/networks/tables/string_down_nodes_displayed.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_edges_displayed.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_input_genes.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_unmapped_genes.tsv"),
        },
    ),
    "regulator_activity": PanelConstructor(
        "regulator_activity", "Regulator activity", "Signed regulator-activity heatmap.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/regulators/figures/regulator_activity", "regulators", "Regulator activity", "contrasts/{contrast}/analyses/regulators/tables/regulator_activity_displayed.tsv", "contrasts/{contrast}/analyses/regulators/tables/regulon_edges.tsv")},
    ),
    "dorothea_grn": PanelConstructor(
        "dorothea_grn", "DoRothEA gene-regulatory network", "Program-aware GRN view backed by complete regulon-edge audit data.", True, "radial",
        {
            "rectangular": _variant("contrasts/{contrast}/analyses/regulators/figures/grn_rectangular", "regulators", "Rectangular GRN", "contrasts/{contrast}/analyses/regulators/tables/grn_nodes_displayed.tsv", "contrasts/{contrast}/analyses/regulators/tables/grn_edges_displayed.tsv", "contrasts/{contrast}/analyses/regulators/tables/regulon_edges.tsv"),
            "radial": _variant("contrasts/{contrast}/analyses/regulators/figures/grn_radial", "regulators", "Radial GRN", "contrasts/{contrast}/analyses/regulators/tables/grn_nodes_displayed.tsv", "contrasts/{contrast}/analyses/regulators/tables/grn_edges_displayed.tsv", "contrasts/{contrast}/analyses/regulators/tables/grn_program_separation_test.tsv", "contrasts/{contrast}/analyses/regulators/tables/grn_sector_summary.tsv"),
            "radial_legacy": _variant("contrasts/{contrast}/analyses/regulators/figures/grn_radial_legacy", "regulators", "Radial GRN (matplotlib)", "contrasts/{contrast}/analyses/regulators/tables/grn_nodes_displayed.tsv", "contrasts/{contrast}/analyses/regulators/tables/grn_edges_displayed.tsv", "contrasts/{contrast}/analyses/regulators/tables/grn_sector_summary.tsv"),
        },
    ),
    "wgcna_module_trait": PanelConstructor(
        "wgcna_module_trait", "WGCNA module associations", "Module-trait association heatmap with small-sample warnings.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/advanced/wgcna/figures/wgcna_module_trait", "wgcna", "WGCNA module associations", "contrasts/{contrast}/analyses/advanced/wgcna/tables/wgcna_module_trait.tsv")},
    ),
    "multilayer_network": PanelConstructor(
        "multilayer_network", "Multilayer network", "Typed GRN, co-expression, and STRING evidence integration.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/advanced/multilayer/figures/multilayer_network", "multilayer", "Multilayer network", "contrasts/{contrast}/analyses/advanced/multilayer/tables/multilayer_nodes.tsv", "contrasts/{contrast}/analyses/advanced/multilayer/tables/multilayer_edges.tsv", "contrasts/{contrast}/analyses/advanced/multilayer/tables/multilayer_triangulated.tsv")},
    ),
}


def constructor_catalog() -> list[dict[str, Any]]:
    """Return a JSON-serializable constructor catalog."""
    return [
        {
            "id": constructor.id,
            "label": constructor.label,
            "description": constructor.description,
            "contrast_specific": constructor.contrast_specific,
            "default_variant": constructor.default_variant,
            "variants": {
                key: {"label": value.label, "required_module": value.required_module}
                for key, value in constructor.variants.items()
            },
        }
        for constructor in PANEL_REGISTRY.values()
    ]


def normalized_gene_panels(panel_config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Expose first-class programs through the legacy grouped-panel shape."""
    panels = {key: dict(value) for key, value in panel_config.get("gene_panels", {}).items()}
    for program_id, program in panel_config.get("programs", {}).items():
        label = program.get("label", program_id.replace("_", " ").title())
        panels[program_id] = {
            "description": program.get("description", label),
            "color": program.get("color"),
            "groups": {label: list(program["genes"])},
        }
    return panels


def validate_recipe_contract(
    recipe_config: dict[str, Any], modules: Iterable[str], contrast_ids: Iterable[str]
) -> list[str]:
    """Validate registered constructors beyond the structural JSON schema."""
    enabled = set(modules)
    contrasts = set(contrast_ids)
    errors: list[str] = []
    for figure_set, recipe in recipe_config.get("figure_sets", {}).items():
        for panel in recipe.get("panels", []):
            constructor_id = panel.get("constructor")
            if not constructor_id:
                continue
            constructor = PANEL_REGISTRY.get(constructor_id)
            prefix = f"figure set {figure_set!r} panel {panel.get('id')!r}"
            if constructor is None:
                errors.append(f"{prefix}: unknown constructor {constructor_id!r}")
                continue
            variant_id = panel.get("variant", constructor.default_variant)
            variant = constructor.variants.get(variant_id)
            if variant is None:
                errors.append(
                    f"{prefix}: constructor {constructor_id!r} has no variant {variant_id!r}; "
                    f"choose from {sorted(constructor.variants)}"
                )
                continue
            if variant.required_module not in enabled:
                errors.append(
                    f"{prefix}: constructor {constructor_id!r}/{variant_id!r} requires "
                    f"analysis module {variant.required_module!r}"
                )
            contrast = panel.get("contrast")
            if constructor.contrast_specific:
                if contrast is None and len(contrasts) != 1:
                    errors.append(f"{prefix}: contrast is required when the project has multiple contrasts")
                elif contrast is not None and contrast not in contrasts:
                    errors.append(f"{prefix}: unknown contrast {contrast!r}")
            elif contrast is not None:
                errors.append(f"{prefix}: constructor {constructor_id!r} is study-level and does not accept contrast")
    return errors


def resolve_panel(project: "ResolvedProject", panel: dict[str, Any]) -> ResolvedPanel:
    """Resolve a recipe panel to figure and displayed-data artifact paths."""
    constructor_id = panel.get("constructor")
    if not constructor_id:
        source = Path(panel["source"]).expanduser()
        source = source if source.is_absolute() else project.result_root / source
        source = source.with_suffix("") if source.suffix.lower() in {".pdf", ".png"} else source
        displayed = []
        for value in panel.get("displayed_data", []):
            path = Path(value).expanduser()
            displayed.append(path if path.is_absolute() else project.result_root / path)
        return ResolvedPanel(None, panel.get("variant"), panel.get("contrast"), panel.get("title", "Legacy source panel"), source, tuple(displayed), None)

    constructor = PANEL_REGISTRY[constructor_id]
    variant_id = panel.get("variant", constructor.default_variant)
    variant = constructor.variants[variant_id]
    contrast: str | None = panel.get("contrast")
    if constructor.contrast_specific and contrast is None:
        contrast = project.contrast_rows[0]["contrast_id"]
    values = {"contrast": contrast or ""}
    return ResolvedPanel(
        constructor_id,
        variant_id,
        contrast,
        panel.get("title", variant.label),
        project.result_root / variant.source.format_map(values),
        tuple(project.result_root / value.format_map(values) for value in variant.displayed_data),
        variant.required_module,
    )


def _write_yaml(path: Path, value: dict[str, Any], force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"refusing to overwrite existing publication file: {path}")
    path.write_text(yaml.safe_dump(value, sort_keys=False), encoding="utf-8")


def _read_contrasts(project_path: Path, config: dict[str, Any]) -> list[dict[str, str]]:
    value = config.get("analysis", {}).get("contrasts")
    if not value:
        raise ValueError("project configuration has no analysis.contrasts path")
    path = Path(str(value)).expanduser()
    path = path if path.is_absolute() else project_path.parent / path
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"contrast file contains no rows: {path}")
    return rows


def initialize_figure_workflow(project_path: str | Path, *, force: bool = False) -> tuple[Path, ...]:
    """Scaffold story files and enable the generic publication constructor module."""
    path = Path(project_path).expanduser().resolve()
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or raw.get("version") != 2:
        raise ValueError("figures init requires a version 2 project.yaml")
    contrasts = _read_contrasts(path, raw)
    contrast_ids = [row["contrast_id"] for row in contrasts]
    claims_path = path.parent / "hypotheses.yaml"
    panels_path = path.parent / "hypothesis_panels.yaml"
    recipe_path = path.parent / "figure_recipe.yaml"
    if not force:
        existing = [candidate for candidate in (claims_path, panels_path, recipe_path) if candidate.exists()]
        if existing:
            raise FileExistsError(
                "refusing to overwrite existing publication files: "
                + ", ".join(str(candidate) for candidate in existing)
            )

    claims = {
        "study": raw.get("project", {}).get("title", raw.get("project", {}).get("id", "study")),
        "hypotheses": [
            {
                "id": f"{contrast_id}_working_hypothesis",
                "statement": "Replace with the biological hypothesis this contrast tests.",
                "contrast": contrast_id,
                "expected_direction": "context_dependent",
                "gene_panels": ["program_1"],
                "pathway_panels": ["pathways_1"],
            }
            for contrast_id in contrast_ids
        ],
    }
    panels = {
        "programs": {
            "program_1": {
                "label": "Biological program 1",
                "description": "Replace with a curated, hypothesis-relevant program.",
                "color": "#D97706",
                "genes": ["REPLACE_WITH_GENE_SYMBOLS"],
                "expected_direction": "not_specified",
            }
        },
        "pathway_panels": {
            "pathways_1": {
                "description": "Replace with pathways central to the hypothesis.",
                "pathways": [{"collection": "custom", "pathway": "REPLACE_WITH_PATHWAY_ID"}],
            }
        },
        "gsea_programs": ["program_1"],
        "program_order": ["Biological program 1"],
    }
    first_contrast = contrast_ids[0]
    recipe = {
        "figure_sets": {
            "primary": {
                "title": "Primary publication figure",
                "description": "Edit constructors, variants, and placement; run tifzoret figures gallery to review alternatives.",
                "width": 12,
                "height": 10,
                "units": "in",
                "columns": 2,
                "shared_legends": True,
                "panels": [
                    {"id": "A", "constructor": "pca_correlation", "row": 1, "column": 1, "column_span": 2},
                    {"id": "B", "constructor": "volcano", "contrast": first_contrast, "row": 2, "column": 1},
                    {"id": "C", "constructor": "de_heatmap", "variant": "global_clustered", "contrast": first_contrast, "row": 2, "column": 2},
                ],
            }
        }
    }
    _write_yaml(claims_path, claims, force)
    _write_yaml(panels_path, panels, force)
    _write_yaml(recipe_path, recipe, force)

    raw["hypotheses"] = {"claims": claims_path.name, "panels": panels_path.name}
    raw["publication"] = {"recipe": recipe_path.name}
    raw.setdefault("analysis", {}).setdefault("modules", {})["publication"] = True
    path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    return claims_path, panels_path, recipe_path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _selected_keys(project: "ResolvedProject") -> set[tuple[str, str, str]]:
    selected: set[tuple[str, str, str]] = set()
    for recipe in (project.recipe_config or {}).get("figure_sets", {}).values():
        for panel in recipe.get("panels", []):
            if panel.get("constructor"):
                resolved = resolve_panel(project, panel)
                selected.add((resolved.constructor or "", resolved.variant or "", resolved.contrast or ""))
    return selected


def iter_review_panels(project: "ResolvedProject") -> Iterable[dict[str, Any]]:
    """Yield every registered variant supported by the enabled project modules."""
    selected = _selected_keys(project)
    contrast_ids = [row["contrast_id"] for row in project.contrast_rows]
    for constructor in PANEL_REGISTRY.values():
        targets: list[str | None] = contrast_ids if constructor.contrast_specific else [None]
        for contrast in targets:
            for variant_id, variant in constructor.variants.items():
                if variant.required_module not in project.modules:
                    continue
                values = {"contrast": contrast or ""}
                source = project.result_root / variant.source.format_map(values)
                displayed = [project.result_root / value.format_map(values) for value in variant.displayed_data]
                yield {
                    "constructor": constructor.id,
                    "constructor_label": constructor.label,
                    "description": constructor.description,
                    "variant": variant_id,
                    "variant_label": variant.label,
                    "contrast": contrast,
                    "selected": (constructor.id, variant_id, contrast or "") in selected,
                    "source_pdf": str(source.with_suffix(".pdf")),
                    "source_png": str(source.with_suffix(".png")),
                    "displayed_data": [str(path) for path in displayed],
                    "available": source.with_suffix(".pdf").is_file() and source.with_suffix(".png").is_file(),
                }


def build_gallery(project: "ResolvedProject", output: str | Path | None = None) -> Path:
    """Create an HTML/PNG review gallery of all available registered variants."""
    outdir = Path(output).expanduser().resolve() if output else project.result_root / "publication" / "gallery"
    assets = outdir / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    records = list(iter_review_panels(project))
    available = [record for record in records if record["available"]]
    for index, record in enumerate(available, start=1):
        name = f"{index:03d}_{record['constructor']}_{record['variant']}_{record['contrast'] or 'study'}.png"
        destination = assets / name
        shutil.copy2(record["source_png"], destination)
        record["review_image"] = f"assets/{name}"
        record["sha256"] = _sha256(Path(record["source_png"]))

    width, tile_width, tile_height, columns = 1600, 380, 300, 4
    rows = max(1, (len(available) + columns - 1) // columns)
    sheet = Image.new("RGB", (width, rows * tile_height), "white")
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.load_default(size=16)
    for index, record in enumerate(available):
        row, column = divmod(index, columns)
        x, y = column * tile_width, row * tile_height
        image = Image.open(outdir / record["review_image"]).convert("RGB")
        image.thumbnail((tile_width - 20, tile_height - 55), Image.Resampling.LANCZOS)
        sheet.paste(image, (x + (tile_width - image.width) // 2, y + 35))
        marker = "SELECTED · " if record["selected"] else ""
        label = f"{marker}{record['constructor']} / {record['variant']}"
        draw.text((x + 8, y + 8), label, fill="#14324A", font=font)
    sheet_path = outdir / "contact_sheet.png"
    sheet.save(sheet_path, dpi=(150, 150))

    cards = []
    for record in records:
        status = "selected" if record["selected"] else "available" if record["available"] else "not built"
        image = f'<img src="{html.escape(record.get("review_image", ""))}" alt="review image">' if record["available"] else '<div class="missing">Not built</div>'
        displayed = "".join(f"<li>{html.escape(value)}</li>" for value in record["displayed_data"])
        cards.append(
            f'<article class="card {status.replace(" ", "-")}">{image}'
            f'<h2>{html.escape(record["constructor_label"])} · {html.escape(record["variant_label"])}</h2>'
            f'<p><strong>{html.escape(status.upper())}</strong> · contrast: {html.escape(record["contrast"] or "study-level")}</p>'
            f'<p>{html.escape(record["description"])}</p><details><summary>Displayed data</summary><ul>{displayed}</ul></details></article>'
        )
    page = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Tifzoret figure gallery</title>
<style>body{{font:15px system-ui;margin:2rem;color:#14324A;background:#F7F9FA}}header{{max-width:70rem}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:1rem}}.card{{background:white;border:1px solid #DCE4E8;border-radius:10px;padding:1rem}}.card.selected{{border:3px solid #0F9D78}}img{{width:100%;height:240px;object-fit:contain;background:white}}.missing{{height:240px;display:grid;place-items:center;background:#EEF2F4;color:#6B7C87}}h2{{font-size:1.05rem}}details{{font-size:.8rem;overflow-wrap:anywhere}}</style></head>
<body><header><h1>Tifzoret publication figure gallery</h1><p>{html.escape(project.project_id)} · generated from registered constructors. Statistical results are unchanged; this page reviews presentation variants.</p><p><a href="contact_sheet.png">Open contact sheet</a></p></header><main class="grid">{''.join(cards)}</main></body></html>"""
    index = outdir / "index.html"
    index.write_text(page, encoding="utf-8")
    (outdir / "gallery.json").write_text(json.dumps({"schema_version": 1, "project": project.project_id, "panels": records}, indent=2) + "\n", encoding="utf-8")
    return index
