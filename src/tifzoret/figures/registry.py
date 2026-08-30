"""Panel registry and constructor catalog."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ..config import ResolvedProject


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
        {"default": _variant("contrasts/{contrast}/analyses/pathways/figures/gsva_heatmap", "pathways", "GSVA heatmap", "contrasts/{contrast}/analyses/pathways/tables/gsva_heatmap_displayed.tsv", "contrasts/{contrast}/analyses/pathways/tables/gsva_scores.tsv", "contrasts/{contrast}/analyses/pathways/tables/gsva_differential.tsv")},
    ),
    "gsea_multitrack": PanelConstructor(
        "gsea_multitrack", "Advanced GSEA curves", "Enrichment score, hits, ranked metric, NES, FDR, and leading edge.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/pathways/figures/gsea_curves", "pathways", "Advanced GSEA", "contrasts/{contrast}/analyses/pathways/tables/gsea_curves_displayed.tsv", "contrasts/{contrast}/analyses/pathways/tables/fgsea.tsv")},
    ),
    "program_heatmap_effects": PanelConstructor(
        "program_heatmap_effects", "Program heatmap and effects", "Consolidated expression heatmap and gene-level effect estimates.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/publication/figures/program_integrated", "publication", "Program heatmap + effects", "contrasts/{contrast}/analyses/publication/tables/program_integrated_displayed.tsv", "contrasts/{contrast}/analyses/publication/tables/program_effects_forest_displayed.tsv", "contrasts/{contrast}/analyses/publication/tables/program_definitions.tsv")},
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
            "upregulated": _variant("contrasts/{contrast}/analyses/networks/figures/string_up_network_community", "networks", "Upregulated STRING community network", "contrasts/{contrast}/analyses/networks/tables/string_up_network_displayed.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_nodes.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_edges.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_input_genes.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_unmapped_genes.tsv"),
            "downregulated": _variant("contrasts/{contrast}/analyses/networks/figures/string_down_network_community", "networks", "Downregulated STRING community network", "contrasts/{contrast}/analyses/networks/tables/string_down_network_displayed.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_nodes.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_edges.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_input_genes.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_unmapped_genes.tsv"),
            "upregulated_legacy": _variant("contrasts/{contrast}/analyses/networks/figures/string_up_network", "networks", "Upregulated STRING network (matplotlib)", "contrasts/{contrast}/analyses/networks/tables/string_up_nodes_displayed.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_edges_displayed.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_input_genes.tsv", "contrasts/{contrast}/analyses/networks/tables/string_up_unmapped_genes.tsv"),
            "downregulated_legacy": _variant("contrasts/{contrast}/analyses/networks/figures/string_down_network", "networks", "Downregulated STRING network (matplotlib)", "contrasts/{contrast}/analyses/networks/tables/string_down_nodes_displayed.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_edges_displayed.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_input_genes.tsv", "contrasts/{contrast}/analyses/networks/tables/string_down_unmapped_genes.tsv"),
        },
    ),
    "regulator_activity": PanelConstructor(
        "regulator_activity", "Regulator activity",
        "Regulator-activity heatmap. The `default` variant plots the primary "
        "(first-declared) regulator view — signed DoRothEA under the legacy dual-view "
        "path; the `binding` variant plots the unsuffixed second view produced by a "
        "configured binding prior (unsigned GTRD occupancy).", True, "default",
        {
            "default": _variant("contrasts/{contrast}/analyses/regulators/figures/regulator_activity", "regulators", "Regulator activity", "contrasts/{contrast}/analyses/regulators/tables/regulator_activity_displayed.tsv", "contrasts/{contrast}/analyses/regulators/tables/regulon_edges.tsv"),
            "binding": _variant("contrasts/{contrast}/analyses/regulators/figures/regulator_activity_binding", "regulators", "Regulator activity (binding prior)", "contrasts/{contrast}/analyses/regulators/tables/regulator_activity_displayed_binding.tsv", "contrasts/{contrast}/analyses/regulators/tables/regulon_edges_binding.tsv"),
        },
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
    # ---- Differential-expression diagnostics (de module) -------------------
    # Calibration and structure checks the de stage already renders per signed
    # contrast but that no manuscript recipe has surfaced. A well-behaved test
    # has a flat-with-a-left-spike p-value histogram; anything else (uniform,
    # U-shaped, mid-bump) signals a modeling problem before any biology.
    "de_pvalue_histogram": PanelConstructor(
        "de_pvalue_histogram", "DE p-value histogram", "Raw p-value distribution for the contrast — a null-calibration diagnostic (expect flat with a spike near zero).", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/de/figures/pvalue_distribution", "de", "p-value histogram", "contrasts/{contrast}/analyses/de/tables/pvalue_distribution_displayed.tsv")},
    ),
    "de_lfc_distribution": PanelConstructor(
        "de_lfc_distribution", "DE log-fold-change distribution", "Shrunken effect-size distribution for the contrast, showing the balance and spread of up- versus downregulation.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/de/figures/lfc_distribution", "de", "LFC distribution", "contrasts/{contrast}/analyses/de/tables/lfc_distribution_displayed.tsv")},
    ),
    "de_pca": PanelConstructor(
        "de_pca", "DE-set sample PCA", "PCA of the samples over the contrast's variance-stabilized expression with group covariance ellipses.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/de/figures/de_pca", "de", "DE-set PCA", "contrasts/{contrast}/analyses/de/tables/de_pca_coordinates.tsv", "contrasts/{contrast}/analyses/de/tables/de_pca_ellipses.tsv")},
    ),
    # ---- Omnibus LRT across >2 groups (de module, omnibus contrasts) -------
    # Rendered only for type: omnibus contrasts; the on-disk gate confines these
    # to those contrasts even though the module tag is the shared "de".
    "omnibus_heatmap": PanelConstructor(
        "omnibus_heatmap", "Omnibus multi-group heatmap", "Row-scaled expression of the genes most significant in the omnibus likelihood-ratio test across all groups.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/omnibus/figures/omnibus_heatmap", "de", "Omnibus heatmap", "contrasts/{contrast}/analyses/omnibus/tables/omnibus_heatmap_displayed.tsv", "contrasts/{contrast}/analyses/omnibus/tables/omnibus_results.tsv")},
    ),
    "omnibus_pvalue_histogram": PanelConstructor(
        "omnibus_pvalue_histogram", "Omnibus p-value histogram", "Raw p-value distribution of the omnibus likelihood-ratio test — a null-calibration diagnostic for the multi-group model.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/omnibus/figures/pvalue_distribution", "de", "Omnibus p-value histogram", "contrasts/{contrast}/analyses/omnibus/tables/pvalue_distribution_displayed.tsv")},
    ),
    # ---- Cross-engine confirmatory DE (de_confirm module) ------------------
    "de_concordance": PanelConstructor(
        "de_concordance", "DESeq2 vs edgeR concordance", "Effect-size concordance between the primary DESeq2 fit and a confirmatory edgeR quasi-likelihood fit for the same contrast.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/de_confirm/figures/de_concordance", "de_confirm", "DE concordance", "contrasts/{contrast}/analyses/de_confirm/tables/de_concordance_displayed.tsv", "contrasts/{contrast}/analyses/de_confirm/tables/edger_results.tsv")},
    ),
    # ---- GO domain breadth (ontology module) -------------------------------
    "ontology_domains": PanelConstructor(
        "ontology_domains", "GO domain breadth", "Enriched Gene Ontology terms faceted across the biological-process, cellular-component, and molecular-function domains.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/ontology/figures/ontology_domains", "ontology", "GO domains", "contrasts/{contrast}/analyses/ontology/tables/ontology_domain_displayed.tsv")},
    ),
    # ---- Enrichment-similarity map (enrichment_map module) -----------------
    "enrichment_map": PanelConstructor(
        "enrichment_map", "Enrichment-similarity map", "Jaccard-overlap network of enriched terms, community-clustered to expose redundant term families rather than a flat list.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/enrichment_map/figures/enrichment_map", "enrichment_map", "Enrichment map", "contrasts/{contrast}/analyses/enrichment_map/tables/enrichment_map_nodes.tsv", "contrasts/{contrast}/analyses/enrichment_map/tables/enrichment_map_clusters.tsv")},
    ),
    # ---- Signed pathway-topology impact (spia module) ----------------------
    "spia_impact": PanelConstructor(
        "spia_impact", "SPIA pathway impact", "Two-evidence SPIA plot combining over-representation and signed pathway-topology perturbation for KEGG pathways.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/spia/figures/spia_two_evidence", "spia", "SPIA impact", "contrasts/{contrast}/analyses/spia/tables/spia_displayed.tsv", "contrasts/{contrast}/analyses/spia/tables/spia_pathways.tsv")},
    ),
    # ---- Ollivier-Ricci curvature of the co-expression graph (curvature) ---
    "curvature_network": PanelConstructor(
        "curvature_network", "Co-expression curvature network", "WGCNA co-expression graph colored by Ollivier-Ricci edge curvature, highlighting bridge edges between modules.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/advanced/curvature/figures/curvature_network", "curvature", "Curvature network", "contrasts/{contrast}/analyses/advanced/curvature/tables/edge_curvature.tsv", "contrasts/{contrast}/analyses/advanced/curvature/tables/curvature_bridges.tsv")},
    ),
    "curvature_histogram": PanelConstructor(
        "curvature_histogram", "Curvature distribution", "Distribution of Ollivier-Ricci curvature over the co-expression edges, summarizing graph bottleneck structure.", True, "default",
        {"default": _variant("contrasts/{contrast}/analyses/advanced/curvature/figures/curvature_histogram", "curvature", "Curvature histogram", "contrasts/{contrast}/analyses/advanced/curvature/tables/gene_curvature.tsv", "contrasts/{contrast}/analyses/advanced/curvature/tables/module_curvature.tsv")},
    ),
    # ---- Variance decomposition (variance_partition module, study-level) ---
    "variance_partition": PanelConstructor(
        "variance_partition", "Variance partition", "Fraction of per-gene expression variance attributable to each design covariate versus residual — the factorial's variance budget.", False, "default",
        {"default": _variant("variance_partition/figures/variance_partition", "variance_partition", "Variance partition", "variance_partition/tables/variance_summary.tsv", "variance_partition/tables/variance_fractions.tsv")},
    ),
    # ---- Cross-contrast consensus (consensus module, study-level) ----------
    "consensus_upset": PanelConstructor(
        "consensus_upset", "Cross-contrast consensus (UpSet)", "UpSet plot of shared and contrast-specific differentially expressed genes across every pairwise contrast.", False, "default",
        {"default": _variant("comparison/figures/consensus_upset", "consensus", "Consensus UpSet", "comparison/tables/consensus_intersections_displayed.tsv", "comparison/tables/contrast_overlap.tsv")},
    ),
    "consensus_sign_heatmap": PanelConstructor(
        "consensus_sign_heatmap", "Cross-contrast sign concordance", "Signed direction of the top consensus genes across contrasts, exposing where effects agree, flip, or drop out.", False, "default",
        {"default": _variant("comparison/figures/consensus_sign_heatmap", "consensus", "Consensus sign heatmap", "comparison/tables/consensus_genes.tsv", "comparison/tables/consensus_membership.tsv")},
    ),
    # ---- Cell-fraction deconvolution (deconvolution module, study-level) ---
    "cell_fractions": PanelConstructor(
        "cell_fractions", "Deconvolved cell fractions", "Per-sample cell-type fractions from signature-matrix deconvolution of the bulk expression.", False, "default",
        {"default": _variant("deconvolution/figures/cell_fractions", "deconvolution", "Cell fractions", "deconvolution/tables/cell_fractions.tsv")},
    ),
    # ---- Factorial interaction views (factorial module, study-level) -------
    # Three project-agnostic views that make a crossed 2x2 (or larger) design's
    # interaction legible from the two configured signed-contrast arms, rather
    # than from a single pairwise comparison.
    "effect_vs_effect": PanelConstructor(
        "effect_vs_effect", "Effect-vs-effect quadrant", "Per-gene log2 fold-change in one factorial arm against the other, with the y=x no-interaction diagonal — genes off the diagonal are the interaction.", False, "default",
        {"default": _variant("factorial/figures/effect_vs_effect", "factorial", "Effect vs effect", "factorial/tables/effect_vs_effect_displayed.tsv")},
    ),
    "interaction_profile": PanelConstructor(
        "interaction_profile", "Interaction profile", "Group-mean expression of the top-interaction genes across one factor, one line per level of the other — non-parallel lines are the interaction.", False, "default",
        {"default": _variant("factorial/figures/interaction_profile", "factorial", "Interaction profile", "factorial/tables/interaction_profile_displayed.tsv")},
    ),
    "group_expression": PanelConstructor(
        "group_expression", "Four-group expression", "Per-sample variance-stabilized expression of the top-interaction genes across every crossed factor group, with group means.", False, "default",
        {"default": _variant("factorial/figures/group_expression", "factorial", "Group expression", "factorial/tables/group_expression_displayed.tsv")},
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
