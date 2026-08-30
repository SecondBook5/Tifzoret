# Walkthrough

This document is a single narrated trace of one real run, following the DAG from entry to completion. It shows how the CLI, configuration, and stage scripts fit together to produce the final outputs.

## Starting point: tifzoret run

You have a configured `project.yaml` that declares:
- A project ID, species, and reference metadata
- An input boundary (counts, BAM, nf-core archive, or ZIP/TAR)
- A design formula and contrast table
- An analysis profile (`standard`, `publication`, or `full`)
- Optional: hypothesis claims, figure panels, and a publication recipe

You invoke:

```bash
tifzoret run project.yaml --cores 4
```

The CLI ([`cli.py:main`](../src/tifzoret/cli.py)) validates the configuration, resolves all relative paths, and hands a validated config object to `_snakemake`, which builds the Snakemake command:

```bash
snakemake \
  --snakefile <packaged-workflow>/Snakefile \
  --configfile project.yaml \
  --cores 4 \
  --use-conda \
  --rerun-incomplete \
  --printshellcmds \
  results/<project_id>/<analysis_set>/manifest.json
```

Snakemake reads [`workflow/Snakefile`](../src/tifzoret/workflow/Snakefile), which `include:`s rule definitions from `workflow/rules/*.smk`, builds a dependency graph targeting `manifest.json`, and executes the DAG.

## Phase 1: Materializing the canonical boundary

**Rule: `materialize_inputs`** ([`core.smk`](../src/tifzoret/workflow/rules/core.smk))  
**Script: [`materialize_inputs.py`](../src/tifzoret/workflow/scripts/01_inputs/materialize_inputs.py)**

If your input is `counts`, the script validates and copies your provided tables. If it's BAM (aligned, nf-core, or archive), the script:
1. Runs `featureCounts` on each BAM to produce integer counts
2. Calculates per-gene exon lengths
3. Derives TPM and FPKM abundance from counts and lengths
4. Parses your GTF to build the annotation table

Regardless of source, this stage produces the four canonical files under `results/<project>/<analysis_set>/inputs/`:
- `counts.tsv` — integer counts matrix (genes × samples)
- `samples.tsv` — sample metadata with design covariates
- `contrasts.tsv` — declared contrasts (factor, numerator, denominator)
- `annotation.tsv` — gene metadata (ID, symbol, biotype, chromosome, coordinates)

And an input manifest: `input_manifest.json`.

**Rule: `resolve_resources`** ([`providers.smk`](../src/tifzoret/workflow/rules/providers.smk))  
**Script: [`resources.R`](../src/tifzoret/workflow/scripts/01_inputs/resources.R)**

In parallel, this stage fetches and caches gene sets from external knowledge bases (MSigDB Hallmarks, GO BP/CC/MF, KEGG, Reactome, DoRothEA regulons). Gene sets are filtered to the species and written to `results/<project>/resources/<source>.rds`. These are versioned and pinned so enrichment tests are reproducible.

**Checkpoint: `tifzoret prepare`** stops here. You can inspect the canonical inputs before committing to the full analysis.

## Phase 2: Study-wide QC

**Rule: `study_qc`** ([`core.smk`](../src/tifzoret/workflow/rules/core.smk))  
**Script: [`qc.R`](../src/tifzoret/workflow/scripts/02_qc/qc.R)**

Runs once per study (not per contrast). Produces:
- Normalized count matrices (VST, rlog)
- PCA, MDS, and hierarchical clustering
- Library-size and dispersion diagnostics
- Correlation heatmaps

Output: `results/<project>/<analysis_set>/qc/*.pdf` and `qc_summary.json`.

This reveals whether samples cluster by their biological condition (good) or by batch (bad).

**Optional rule: `study_batch`** ([`core.smk`](../src/tifzoret/workflow/rules/core.smk))  
**Script: [`batch.R`](../src/tifzoret/workflow/scripts/02_qc/batch.R)**

If your samples have a `batch` column, this stage produces batch-corrected PCA and distance views using `limma::removeBatchEffect`. These are diagnostic only — the DE models handle batch as a covariate.

Output: `results/<project>/<analysis_set>/batch/*.pdf`.

## Phase 3: Differential expression (per contrast)

From here on, stages fan out: each declared contrast runs independently through the enabled modules.

**Rule: `contrast_de`** ([`core.smk`](../src/tifzoret/workflow/rules/core.smk))  
**Script: [`de.R`](../src/tifzoret/workflow/scripts/03_differential/de.R)**

For each contrast `<contrast_id>`, this stage:
1. Builds a DESeq2 dataset from counts and samples
2. Estimates size factors and dispersions
3. Fits the design formula (e.g., `~batch + condition`)
4. Extracts the coefficient for `numerator - denominator`
5. Applies shrinkage (apeglm, ashr, normal, or none per config)
6. Writes results with log2FC, adjusted p-value, and base mean

Output: `results/<project>/<analysis_set>/<contrast>/de/de_results.tsv` and summary plots (`volcano.pdf`, `ma_plot.pdf`, `pvalue_hist.pdf`).

**Optional rule: `contrast_de_confirm`** ([`modules.smk`](../src/tifzoret/workflow/rules/modules.smk))  
**Script: [`de_confirm.R`](../src/tifzoret/workflow/scripts/03_differential/de_confirm.R)**

If enabled, re-runs the contrast with edgeR's quasi-likelihood framework. Genes that agree on direction between DESeq2 and edgeR are flagged as high-confidence.

Output: `results/<project>/<analysis_set>/<contrast>/de_confirm/*.tsv` and agreement plots.

## Phase 4: Enrichment (per contrast)

**Rule: `contrast_pathways`** ([`modules.smk`](../src/tifzoret/workflow/rules/modules.smk))  
**Script: [`pathways.R`](../src/tifzoret/workflow/scripts/04_enrichment/pathways.R)**

Runs four enrichment methods (ORA, GSEA via fgsea, GSVA, ssGSEA) against the cached gene sets from Phase 1. Produces:
- Ranked gene-set results with NES, p-value, adjusted p-value
- GSEA curves for top pathways
- Consensus tables (pathways that pass in multiple methods)

Output: `results/<project>/<analysis_set>/<contrast>/pathways/*.tsv` and `gsea_curves.pdf`.

**Rule: `contrast_ontology`** ([`modules.smk`](../src/tifzoret/workflow/rules/modules.smk))  
**Script: [`ontology.R`](../src/tifzoret/workflow/scripts/04_enrichment/ontology.R)**

Slices the ORA results by ontology source (GO BP, KEGG, Reactome) and renders domain-specific views.

Output: `results/<project>/<analysis_set>/<contrast>/ontology/*.pdf`.

## Phase 5: Composition (per contrast)

**Rule: `contrast_composition`** ([`modules.smk`](../src/tifzoret/workflow/rules/modules.smk))  
**Script: [`composition.R`](../src/tifzoret/workflow/scripts/05_composition/composition.R)**

Computes per-sample cell-state signature scores using ssGSEA. If you provided a `composition.signatures` mapping in `project.yaml`, this stage tests whether signature enrichment differs between conditions.

Output: `results/<project>/<analysis_set>/<contrast>/composition/signature_scores.tsv` and boxplots.

## Phase 6: Regulators (per contrast)

**Rule: `contrast_regulators`** ([`modules.smk`](../src/tifzoret/workflow/rules/modules.smk))  
**Script: [`regulators.R`](../src/tifzoret/workflow/scripts/06_regulators/regulators.R)**

Infers transcription-factor activity from gene expression and a signed TF→target prior network (DoRothEA by default). Uses VIPER via decoupleR.

Output: `results/<project>/<analysis_set>/<contrast>/regulators/tf_activity.tsv` and ranked TF plots.

**Rule: `contrast_grn`** ([`modules.smk`](../src/tifzoret/workflow/rules/modules.smk))  
**Script: [`grn.py`](../src/tifzoret/workflow/scripts/06_regulators/grn.py)**

Builds a gene regulatory network (GRN) by selecting top differential TFs and their top DE targets. Assigns nodes to expression programs (up-regulated, down-regulated).

Output: `results/<project>/<analysis_set>/<contrast>/grn/grn_edges.tsv` and `grn_nodes.tsv`.

**Rule: `contrast_grn_radial`** ([`modules.smk`](../src/tifzoret/workflow/rules/modules.smk))  
**Script: [`grn_radial.R`](../src/tifzoret/workflow/scripts/06_regulators/grn_radial.R)**

Renders a publication-grade radial GRN figure: TFs in the center, targets on the periphery, colored by program.

Output: `results/<project>/<analysis_set>/<contrast>/grn/grn_radial.pdf`.

## Phase 7: Networks (per contrast)

**Rule: `contrast_networks`** ([`modules.smk`](../src/tifzoret/workflow/rules/modules.smk))  
**Script: [`networks.py`](../src/tifzoret/workflow/scripts/07_networks/networks.py)**

Fetches protein-protein interaction edges from STRING for the top DE genes. Produces a coverage audit (how many DE genes are present in STRING).

Output: `results/<project>/<analysis_set>/<contrast>/networks/string_edges.tsv` and audit JSON.

**Rule: `contrast_string_figures`** ([`modules.smk`](../src/tifzoret/workflow/rules/modules.smk))  
**Scripts: [`string_figures.R`](../src/tifzoret/workflow/scripts/07_networks/string_figures.R) and [`string_network.R`](../src/tifzoret/workflow/scripts/07_networks/string_network.R)**

From STRING edges:
1. `string_figures.R` runs functional enrichment on STRING communities and renders bubble plots
2. `string_network.R` detects communities via igraph and renders a publication-grade network layout

Output: `results/<project>/<analysis_set>/<contrast>/networks/*.pdf`.

**Optional rule: `contrast_wgcna`** ([`advanced.smk`](../src/tifzoret/workflow/rules/advanced.smk))  
**Script: [`wgcna.R`](../src/tifzoret/workflow/scripts/07_networks/wgcna.R)**

If enabled in `full` profile, discovers weighted co-expression modules. Identifies clusters of genes with correlated expression and tests which modules associate with traits.

Output: `results/<project>/<analysis_set>/<contrast>/wgcna/module_assignment.tsv` and dendrograms.

**Optional rule: `contrast_curvature`** ([`advanced.smk`](../src/tifzoret/workflow/rules/advanced.smk))  
**Script: [`curvature.py`](../src/tifzoret/workflow/scripts/07_networks/curvature.py)**

If WGCNA ran, computes Ollivier-Ricci curvature on the co-expression graph. Positive curvature indicates robust, redundant neighborhoods; negative curvature marks bridging genes between modules.

Output: `results/<project>/<analysis_set>/<contrast>/curvature/curvature.tsv`.

## Phase 8: Synthesis

**Rule: `contrast_hypotheses`** ([`modules.smk`](../src/tifzoret/workflow/rules/modules.smk))  
**Script: [`hypotheses.py`](../src/tifzoret/workflow/scripts/09_synthesis/hypotheses.py)**

If you declared hypotheses in `hypotheses.yaml`, this stage evaluates each claim against the DE results and enrichment scores. Produces structured pass/fail verdicts with confidence levels.

Output: `results/<project>/<analysis_set>/<contrast>/hypotheses/verdicts.json` and summary tables.

**Rule: `study_consensus`** ([`core.smk`](../src/tifzoret/workflow/rules/core.smk))  
**Script: [`consensus.py`](../src/tifzoret/workflow/scripts/09_synthesis/consensus.py)**

Runs once per study. Identifies genes regulated the same direction across multiple contrasts (reproducibly regulated genes).

Output: `results/<project>/<analysis_set>/consensus/consensus_genes.tsv`.

## Phase 9: Publication figures

**Rule: `contrast_publication`** ([`publication.smk`](../src/tifzoret/workflow/rules/publication.smk))  
**Script: [`publication.R`](../src/tifzoret/workflow/scripts/10_figures/publication.R)**

If you provided a `publication.recipe`, this stage renders the requested figure panels (heatmaps, volcano plots, boxplots) with hypothesis-driven gene highlighting and consistent publication theme.

Output: `results/<project>/<analysis_set>/publication/<figure_set>/panels/*.pdf`.

**Rule: `assemble_figure`** ([`publication.smk`](../src/tifzoret/workflow/rules/publication.smk))  
**Script: [`assemble.py`](../src/tifzoret/workflow/scripts/10_figures/assemble.py)**

Combines recipe-selected panel PDFs into multi-panel publication figures. Renders both PDF and PNG raster previews.

Output: `results/<project>/<analysis_set>/publication/<figure_set>/<figure_name>.pdf` and `.png`.

## Phase 10: Report and provenance

**Rule: `report_html`** ([`report.smk`](../src/tifzoret/workflow/rules/report.smk))  
**Script: [`report.py`](../src/tifzoret/workflow/scripts/10_figures/report.py)**

Collects all figures, tables, and summaries into a single self-contained HTML report. No server or external dependencies required — open it in any browser.

Output: `results/<project>/<analysis_set>/REPORT.html`.

**Rule: `release_manifest`** ([`report.smk`](../src/tifzoret/workflow/rules/report.smk))  
**Script: [`manifest.py`](../src/tifzoret/workflow/scripts/10_figures/manifest.py)**

The final stage. Walks the entire result tree, computes SHA-256 checksums for every output, and writes a JSON manifest. This is the workflow's target: when `manifest.json` exists and passes validation, the run is complete.

Output: `results/<project>/<analysis_set>/manifest.json`.

**Optional rule: `front_door_artifacts`** ([`publication.smk`](../src/tifzoret/workflow/rules/publication.smk))  
**Script: [`front_door.py`](../src/tifzoret/workflow/scripts/10_figures/front_door.py)**

If configured, copies reviewer-critical outputs (assembled figures, REPORT.html, manifest) to a flat `front_door/` directory for simplified external access.

Output: `results/<project>/<analysis_set>/front_door/`.

## What you get

After the run completes, you have:

- **`results/<project>/<analysis_set>/`** — the entire result tree
  - `inputs/` — canonical count, sample, contrast, and annotation tables
  - `qc/` — study-wide quality-control figures
  - `<contrast>/` — per-contrast results
    - `de/` — differential expression results and volcano/MA plots
    - `pathways/`, `ontology/`, `composition/`, `regulators/`, `grn/`, `networks/` — enrichment and network outputs
    - `hypotheses/` — hypothesis evaluation verdicts
  - `publication/` — assembled multi-panel figures (if configured)
  - `REPORT.html` — navigable HTML report with all results
  - `manifest.json` — checksummed provenance for every output

Snakemake caches all intermediate outputs. Re-running with the same config resumes from cached results; changing only one contrast re-runs only that contrast's downstream stages.

## Inspecting the plan

Before committing to a full run, you can:

- **`tifzoret dry-run project.yaml`** — prints the DAG and planned commands without executing
- **`tifzoret prepare project.yaml`** — runs only to canonical inputs, then stops
- **`tifzoret validate project.yaml`** — checks config validity without invoking Snakemake

These let you verify the plan before spinning up the full analysis.
