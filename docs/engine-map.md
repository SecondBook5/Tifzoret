# Engine map

This document is the ordered index — the literal "entry → step 1 → step 2 → … → what curvature means" artifact that maps every analysis stage to its purpose, scientific motivation, and place in the pipeline.

## Entry point

`tifzoret run project.yaml` invokes [`cli.py:main`](../src/tifzoret/cli.py) (line 387), which validates the configuration, then calls `_snakemake` (line 71) to build and execute the Snakemake invocation. The workflow is defined in [`workflow/Snakefile`](../src/tifzoret/workflow/Snakefile), which `include:`s rule groups from `workflow/rules/*.smk`. Each rule shells out to a stage script via `{input.script}`.

The dependency graph for a configured project is visualized in [`dag.svg`](dag.svg).

## The ten phases

Tifzoret organizes its 35 stage scripts into ten sequential phases. Every stage is listed below with:
- Its target location in the phase-organized tree (shown as `NN_phase/script`; links point to the current flat location in `workflow/scripts/`)
- What it does (one-line summary)
- Why it exists (the scientific motivation)
- The Snakemake rule that invokes it (and which `.smk` file defines that rule)
- The conda environment it runs in

### Phase 01: Inputs

The input boundary: materializing canonical tables from heterogeneous upstream sources and fetching cached knowledge bases.

| Stage | What | Why | Rule | Env |
|-------|------|-----|------|-----|
| [01_inputs/materialize_inputs.py](../src/tifzoret/workflow/scripts/01_inputs/materialize_inputs.py) | Materialize external sources into canonical input tables | Transforms heterogeneous upstream data (featureCounts, Salmon, nf-core, GEO, archive) into four uniform TSVs (counts, samples, contrasts, annotation); downstream is blind to source | `materialize_inputs` (core.smk) | core.yaml |
| [01_inputs/resources.R](../src/tifzoret/workflow/scripts/01_inputs/resources.R) | Fetch and cache gene sets from external knowledge bases | Provides pathway/ontology/regulon annotations for enrichment tests; pinned versions prevent drift | `resolve_resources` (providers.smk) | r.yaml |
| [01_inputs/export_symbol_map.R](../src/tifzoret/workflow/scripts/01_inputs/export_symbol_map.R) | Generate a frozen gene ID → symbol TSV from annotation contract | Provides a consistent symbol lookup for materialize_inputs.py when annotation.tsv is unavailable at ingestion time | manual helper (not a Snakemake rule) | r.yaml |

### Phase 02: QC

Study-wide quality control, exploratory analysis, and batch diagnostics. These run once per study, not per contrast.

| Stage | What | Why | Rule | Env |
|-------|------|-----|------|-----|
| [02_qc/qc.R](../src/tifzoret/workflow/scripts/02_qc/qc.R) | Study-wide quality control and exploratory analysis | Reveals batch effects, outliers, and baseline variance structure; validates that declared design factors separate in unsupervised space | `study_qc` (core.smk) | r.yaml |
| [02_qc/variancepartition.R](../src/tifzoret/workflow/scripts/02_qc/variancepartition.R) | Variance decomposition across design covariates | Quantifies how much of each gene's expression variance is explained by each covariate; surfaces whether batch dominates signal | `study_variance_partition` (core.smk) | r.yaml |
| [02_qc/batch.R](../src/tifzoret/workflow/scripts/02_qc/batch.R) | Batch-corrected ordination and sample-distance views | Reveals whether nominal batch drives clustering; diagnostic for when to use batch correction in DE models | `study_batch` (core.smk) | r.yaml |
| [02_qc/sva.R](../src/tifzoret/workflow/scripts/02_qc/sva.R) | Surrogate variable analysis (latent technical confounders) | Estimates hidden batch effects; compares DE with/without SVA adjustment to quantify hidden confounder impact | `contrast_sva` (advanced.smk) | r.yaml |

### Phase 03: Differential expression

Per-contrast differential expression engines. These run once per declared contrast.

| Stage | What | Why | Rule | Env |
|-------|------|-----|------|-----|
| [03_differential/family_fit.R](../src/tifzoret/workflow/scripts/03_differential/family_fit.R) | The single Wald fit for one estimand family | Every estimand in a family must share one filter, one set of size factors, one set of dispersions, one sample universe, and one coefficient covariance — structurally, not by convention | `family_fit` (core.smk) | r.yaml |
| [03_differential/estimand.R](../src/tifzoret/workflow/scripts/03_differential/estimand.R) | One estimand extracted from its family's shared Wald fit | Every reported SE is the exact sqrt(c' Sigma c); reparameterized shrinkage allows apeglm on composite contrasts | `family_estimand` (core.smk) | r.yaml |
| [03_differential/de.R](../src/tifzoret/workflow/scripts/03_differential/de.R) | Primary differential expression per contrast (DESeq2 ~condition) | Identifies genes whose expression differs between conditions; produces log2FC + adjusted p-value for every gene | `contrast_de` (core.smk) | r.yaml |
| [03_differential/omnibus.R](../src/tifzoret/workflow/scripts/03_differential/omnibus.R) | Omnibus differential expression (multi-level factor LRT) | Asks whether a gene differs across ANY level of a factor (>2 levels); complements pairwise contrasts for multi-group designs | `contrast_omnibus` (core.smk) | r.yaml |
| [03_differential/factorial.R](../src/tifzoret/workflow/scripts/03_differential/factorial.R) | Factorial interaction exploratory figures (2×2 designs) | Makes interaction (difference-of-differences) legible: genes whose response to treatment A depends on treatment B | `study_factorial` (core.smk) | r.yaml |
| [03_differential/de_confirm.R](../src/tifzoret/workflow/scripts/03_differential/de_confirm.R) | Confirmatory differential expression (edgeR quasi-likelihood) | Two independent NB engines (DESeq2 + edgeR) that agree on direction reduce false discoveries; gives reviewer confidence | `contrast_de_confirm` (modules.smk) | r.yaml |

### Phase 04: Enrichment

Pathway and ontology enrichment analysis: which biological processes/pathways are overrepresented?

| Stage | What | Why | Rule | Env |
|-------|------|-----|------|-----|
| [04_enrichment/pathways.R](../src/tifzoret/workflow/scripts/04_enrichment/pathways.R) | Multi-method pathway enrichment (ORA, GSEA, GSVA, ssGSEA) | Identifies which biological pathways/gene-sets are enriched in the differential gene list; consensus across methods is robust | `contrast_pathways` (modules.smk) | r.yaml |
| [04_enrichment/ontology.R](../src/tifzoret/workflow/scripts/04_enrichment/ontology.R) | Ontology-specific views (GO, KEGG, Reactome ORA results) | Reveals which biological processes, pathways, and cellular components are enriched; structured ontology navigation | `contrast_ontology` (modules.smk) | r.yaml |
| [04_enrichment/spia.R](../src/tifzoret/workflow/scripts/04_enrichment/spia.R) | Topology-aware pathway perturbation (SPIA) | Beyond enrichment: propagates DE fold-changes through KEGG reaction graphs; distinguishes activation vs inhibition | `contrast_spia` (modules.smk) | r.yaml |
| [04_enrichment/enrichment_map.py](../src/tifzoret/workflow/scripts/04_enrichment/enrichment_map.py) | Enrichment-similarity network (redundancy reduction) | Clusters enriched terms by shared leading-edge genes; reveals which pathways are distinct vs. redundant | `contrast_enrichment_map` (modules.smk) | network.yaml |

### Phase 05: Composition

Cell-state signature enrichment and cell-fraction deconvolution: what cell types/states are present?

| Stage | What | Why | Rule | Env |
|-------|------|-----|------|-----|
| [05_composition/composition.R](../src/tifzoret/workflow/scripts/05_composition/composition.R) | Cell-state signature enrichment (ssGSEA) per sample | Quantifies cell-state/marker enrichment as sample-level scores; tests whether composition differs between conditions | `contrast_composition` (modules.smk) | r.yaml |
| [05_composition/deconvolution.R](../src/tifzoret/workflow/scripts/05_composition/deconvolution.R) | Signature-matrix cell-fraction deconvolution (NNLS) | Estimates each bulk sample's cell-type composition from a reference signature matrix; quantifies compositional shifts | `study_deconvolution` (core.smk) | r.yaml |

### Phase 06: Regulators

Transcription-factor activity inference and gene regulatory networks.

| Stage | What | Why | Rule | Env |
|-------|------|-----|------|-----|
| [06_regulators/regulators.R](../src/tifzoret/workflow/scripts/06_regulators/regulators.R) | Transcription-factor activity inference (VIPER / decoupleR) | Infers TF activity from target-gene expression and signed prior networks; identifies master regulators | `contrast_regulators` (modules.smk) | r.yaml |
| [06_regulators/grn.py](../src/tifzoret/workflow/scripts/06_regulators/grn.py) | Gene regulatory network data layer (program-aware regulon subgraph) | Selects top differential TFs + their top DE targets; assigns nodes to expression programs for visualization | `contrast_grn` (modules.smk) | network.yaml |
| [06_regulators/grn_radial.R](../src/tifzoret/workflow/scripts/06_regulators/grn_radial.R) | Gene regulatory network figure layer (publication radial panel) | Shows which regulators drive the contrast's expression programs and their target genes; publication-grade radial layout | `contrast_grn_radial` (modules.smk) | r.yaml |

### Phase 07: Networks

Co-expression networks, protein-protein interaction networks, and multi-layer integration.

| Stage | What | Why | Rule | Env |
|-------|------|-----|------|-----|
| [07_networks/wgcna.R](../src/tifzoret/workflow/scripts/07_networks/wgcna.R) | Weighted gene co-expression network analysis | Discovers co-expression modules (genes with correlated expression); tests which modules associate with traits | `contrast_wgcna` (advanced.smk) | r.yaml |
| [07_networks/networks.py](../src/tifzoret/workflow/scripts/07_networks/networks.py) | STRING protein-protein interaction network fetch + audit | Provides physical/functional interaction context for DE genes; fetches from STRING API with coverage audit | `contrast_networks` (modules.smk) | network.yaml |
| [07_networks/string_figures.R](../src/tifzoret/workflow/scripts/07_networks/string_figures.R) | STRING functional-enrichment bubble visualization | Shows which biological processes are enriched in the STRING subnetworks; reveals functional themes | `contrast_string_figures` (modules.smk) | r.yaml |
| [07_networks/string_network.R](../src/tifzoret/workflow/scripts/07_networks/string_network.R) | STRING community network visualization (publication panel) | Shows how top DE genes cluster into functional communities; igraph community detection + publication-grade layout | `contrast_string_figures` (modules.smk) | r.yaml |
| [07_networks/curvature.py](../src/tifzoret/workflow/scripts/07_networks/curvature.py) | Ollivier-Ricci curvature of the WGCNA co-expression graph | +curvature marks redundant, robust neighbourhoods; −curvature identifies bridging genes between modules | `contrast_curvature` (advanced.smk) | network.yaml |
| [07_networks/multilayer.py](../src/tifzoret/workflow/scripts/07_networks/multilayer.py) | Multi-layer network triangulation (GRN + WGCNA + STRING) | Identifies genes supported by multiple independent evidence types; triangulation is more robust than any single layer | `contrast_multilayer` (advanced.smk) | network.yaml |

### Phase 08: Causal

Causal mediation analysis and power estimation.

| Stage | What | Why | Rule | Env |
|-------|------|-----|------|-----|
| [08_causal/mediation.R](../src/tifzoret/workflow/scripts/08_causal/mediation.R) | Causal mediation analysis (pathway scores as mediators) | Tests whether a mediator pathway accounts for treatment effect on outcome; decomposes total effect into direct + indirect | `contrast_mediation` (advanced.smk) | r.yaml |
| [08_causal/power.py](../src/tifzoret/workflow/scripts/08_causal/power.py) | Mediation power analysis (sample-size requirements) | Estimates per-group N needed for 80% power on mediator and outcome effects; informs future study design | `contrast_mediation_power` (advanced.smk) | core.yaml |

### Phase 09: Synthesis

Cross-contrast synthesis: reproducibly regulated genes and hypothesis evaluation.

| Stage | What | Why | Rule | Env |
|-------|------|-----|------|-----|
| [09_synthesis/consensus.py](../src/tifzoret/workflow/scripts/09_synthesis/consensus.py) | Cross-contrast consensus (reproducibly regulated genes) | Identifies genes regulated the same way across multiple contrasts; reproducibility is strong evidence | `study_consensus` (core.smk) | network.yaml |
| [09_synthesis/hypotheses.py](../src/tifzoret/workflow/scripts/09_synthesis/hypotheses.py) | Hypothesis evaluation engine (claim → evidence → verdict) | Tests configured biological claims (gene/pathway direction) against results; produces structured pass/fail + confidence | `contrast_hypotheses` (modules.smk) | core.yaml |

### Phase 10: Figures

Publication figure assembly, HTML report, and release provenance.

| Stage | What | Why | Rule | Env |
|-------|------|-----|------|-----|
| [10_figures/publication.R](../src/tifzoret/workflow/scripts/10_figures/publication.R) | Publication-grade figure panels (heatmaps, volcano, boxplots) | Renders manuscript-ready figures with study-specific gene panels, hypothesis-driven coloring, and consistent theme | `contrast_publication` (publication.smk) | r.yaml |
| [10_figures/assemble.py](../src/tifzoret/workflow/scripts/10_figures/assemble.py) | Multi-panel figure assembly (PDF + PNG raster preview) | Combines recipe-selected panel PDFs into publication figures with grid layout; renders raster preview for review | `assemble_figure` (publication.smk) | core.yaml |
| [10_figures/report.py](../src/tifzoret/workflow/scripts/10_figures/report.py) | Self-contained offline HTML report builder | Produces a single shareable REPORT.html with all figures/tables embedded; no server or external dependencies | `report_html` (report.smk) | core.yaml |
| [10_figures/manifest.py](../src/tifzoret/workflow/scripts/10_figures/manifest.py) | Release manifest with SHA-256 checksums | Provides auditable provenance for all workflow outputs; checksums detect tampering or corruption | `release_manifest` (report.smk) | core.yaml |
| [10_figures/front_door.py](../src/tifzoret/workflow/scripts/10_figures/front_door.py) | Review-facing artifact promotion (front_door/ populated) | Copies reviewer-critical outputs to a flat, documented front_door/ directory; simplifies external access | `front_door_artifacts` (publication.smk) | core.yaml |

## Shared helpers

These scripts are sourced by pipeline stages rather than invoked as stages themselves.

| Helper | What | Why |
|--------|------|-----|
| [utils.R](../src/tifzoret/workflow/scripts/utils.R) | Common utilities (plotting theme, color scales, I/O) | Ensures visual consistency and reduces duplication across R stages |
| [estimands.R](../src/tifzoret/workflow/scripts/estimands.R) | Cell weights → contrast vector mapping | The single authority on what an estimand means numerically; sourced by family_fit.R, estimand.R, term_test.R, de_confirm.R so DESeq2 and edgeR can never disagree on a contrast's sign |
| [de_render.R](../src/tifzoret/workflow/scripts/03_differential/de_render.R) | Differential expression rendering (volcano, MA, PCA, heatmap) | Shared rendering logic for both de.R and estimand.R; display-sample selection driven by cell membership rather than contrast type |

## How the phases relate

The ten phases form a dependency chain:

- **01_inputs** materializes the canonical boundary and fetches gene sets. Everything downstream depends on this.
- **02_qc** runs once per study to reveal batch effects and variance structure.
- **03_differential** produces per-contrast gene-level statistics. All enrichment and network stages consume these.
- **04_enrichment** and **05_composition** interpret gene lists through pathways and cell-state signatures.
- **06_regulators** infers TF activity and builds GRNs from signed prior networks.
- **07_networks** discovers co-expression modules (WGCNA), fetches PPI from STRING, and triangulates evidence across layers. The **curvature** stage consumes WGCNA output to quantify module topology.
- **08_causal** tests mediation hypotheses; **power** pairs with it to estimate sample-size requirements.
- **09_synthesis** identifies cross-contrast consensus genes and evaluates configured biological hypotheses.
- **10_figures** assembles publication panels, the HTML report, and release provenance.

Every stage is reusable: the same scripts serve all contrasts, and the same engine serves all studies.
