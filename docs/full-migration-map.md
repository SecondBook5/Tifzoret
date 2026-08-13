# BulkRNAFrame full migration map

Status date: 2026-08-11

Inventory baseline:

- biological source HEAD: `f48aec9180c01eee2b13ac68914c0797e1af2830`
  with an intentionally preserved dirty worktree containing current biological
  analyses, v2 consumer configs, results, and manuscript assets;
- BulkRNAFrame implementation HEAD before this planning update:
  `881b2c29f7e4a7665809626747b572bc8cc5c3f5`.

This document is the authoritative map for separating the reusable bulk RNA-seq
workflow from `lymphatic-flow-homeostasis` into the public BulkRNAFrame tool.
It records where every major source capability belongs, what has already been
implemented, what has actually passed parity checks, and what must remain in
the biological consumer repository.

The migration is not complete merely because an implementation exists. A
capability is safe to retire from the source repository only after its declared
acceptance evidence has passed against an external golden reference.

## Status language

| Status | Meaning |
|---|---|
| **Implemented** | Reusable code and an output contract exist in BulkRNAFrame. |
| **Parity verified** | The relevant candidate artifacts passed a declared comparison with the biological golden reference. |
| **Partial** | A usable implementation exists, but either legacy breadth or verification coverage is incomplete. |
| **Planned** | The capability has a destination and acceptance contract but has not been completed. |
| **Consumer-owned** | The content is biological study configuration or evidence and must not move into the tool implementation. |
| **Retire after parity** | The source implementation is duplicated and may be removed only after its migration gates pass. |
| **Reference only** | The artifact is evidence for comparison, not reusable source code. |

## Repository ownership boundary

### BulkRNAFrame owns

- schemas, configuration resolution, CLI behavior, profiles, dependency
  resolution, and workflow orchestration;
- BAM, nf-core/rnaseq, archive, and count-matrix input adapters;
- canonical inputs, QC, differential expression, enrichment, ontology,
  composition scoring, regulator inference, networks, hypotheses, advanced
  systems analyses, and collection analyses;
- species/resource provider interfaces and deterministic resource caching;
- generic publication-grade figure constructors and recipe-driven assembly;
- reports, manifests, verification logic, tests, package metadata, locked
  environments, containers, and public documentation.

### The biological repository owns

- sample sheets, cohort inclusion/exclusion, study-specific covariates, and
  contrast declarations;
- GTF/reference selections for each run and locations of upstream BAM or
  nf-core/rnaseq products;
- biological claims, named programs, selected genes/pathways, expected
  directions, program colors, panel choices, and figure layouts;
- study results and immutable golden-reference snapshots;
- manuscript-facing captions, narrative, and journal-specific assembly or
  submission files;
- the invocation of nf-core/rnaseq from FASTQ, because BulkRNAFrame begins at
  aligned BAMs or validated counts.

### Generated data rule

Generated study results never become package data. Small synthetic fixtures and
license-compatible provider fixtures may live in BulkRNAFrame tests; real CAPE,
ligation, aging, and RelA results remain external consumer evidence.

## Current evidence and its limit

The five current migration reports in
`lymphatic-flow-homeostasis/results/bulk_rna_frame_v2/migration/` are green for:

- CAPE primary;
- CAPE full;
- ligation versus sham;
- old versus young with the sex covariate;
- RelA versus WT through the nf-core input boundary.

Those reports currently prove only:

1. exact canonical count equality;
2. configured DE numerical tolerances and identical significance/direction
   decisions; and
3. validity of required PDF/PNG figure pairs and front-door indexes.

They do **not** yet prove enrichment-universe or leading-edge parity, GSVA score
parity, regulator/network audit-table parity, displayed-data equality, report or
manifest completeness, or bounded perceptual similarity. The current reports
must therefore be described as **core migration passes**, not full workflow
parity passes.

Current tool evidence includes a passing local test suite (34 passed and one
scheduled live-provider integration test skipped) and a successful full-profile
smoke run. Publication to PyPI/OCI, fresh-machine execution, human parity, and
live-provider release testing remain release gates.

Any consumer documentation that currently calls the legacy `--scope all`
reports “complete” must be corrected when the expanded verifier lands. The
files are valid evidence for their three implemented gates, but their scope
name is broader than their present behavior.

## Target project and output contracts

Every consumer project is expressed as:

```text
my-study/
├── project.yaml
├── samples.tsv
├── contrasts.tsv
├── hypotheses.yaml          # optional
├── hypothesis_panels.yaml   # optional
└── figure_recipe.yaml       # optional
```

All downstream modules consume the canonical files:

```text
counts.tsv
samples.tsv
annotation.tsv
contrasts.tsv
input_manifest.json
```

`gene_id` and `sample_id` are unique keys. Every signed effect is explicitly
`numerator - denominator`; sample order and group names never determine
direction.

The stable result layout is:

```text
results/<project>/<analysis_set>/
├── REPORT.html
├── manifest.json
├── figures/
├── tables/
├── inputs/
├── qc/
├── contrasts/<contrast_id>/analyses/
│   ├── de/
│   ├── pathways/
│   ├── ontology/
│   ├── composition/
│   ├── regulators/
│   ├── networks/
│   ├── hypotheses/
│   ├── publication/
│   └── advanced/
├── publication/<figure_set>/
│   ├── panels/
│   └── assembled/
└── .cache/
    ├── logs/
    └── resources/
```

## Configuration, orchestration, and shared-library mapping

| Source in `lymphatic-flow-homeostasis` | BulkRNAFrame destination | Disposition | Current status and remaining gate |
|---|---|---|---|
| `config/config.schema.yaml`, `config/defaults.yaml` | `src/bulk_rna_frame/schemas/project.schema.yaml`, `config.py`, input templates | Replace study-oriented schema with layered v2 contracts | **Implemented**; schemas and migration tests pass. Real consumer configurations still need full parity gates. |
| `config/deconvolution_signatures.yaml` | `signatures.schema.yaml` plus consumer-owned signature file referenced by `project.yaml` | Separate neutral schema from curated content | **Implemented** interface; the actual lymphatic signatures remain **consumer-owned**. |
| `workflow/Snakefile` | packaged `src/bulk_rna_frame/workflow/Snakefile` | Split orchestration into installable rule groups | **Implemented**. |
| `workflow/rules/core.smk`, `workflow/rules/publication.smk` | packaged `rules/core.smk`, `providers.smk`, `modules.smk`, `advanced.smk`, `publication.smk`, `report.smk` | Replace monolithic study orchestration with profile-aware modules | **Implemented**; DAG/unit and smoke evidence exists. |
| `workflow/lib/analysis_set.py` | `config.py`, canonical result-root resolution, CLI | Preserve explicit analysis-set selection | **Implemented**. |
| `workflow/lib/configmerge.py` | `config.py` schema defaults, normalization, module resolution | Replace ad hoc merge behavior | **Implemented**. |
| `workflow/lib/annotate.R` | `materialize_inputs.py`, `utils.R`, provider annotation helpers | Preserve GTF-first annotation and unique gene keys | **Implemented**; counts/core parity exists, annotation-table parity still needs an explicit gate. |
| `workflow/lib/contrast.R` | `config.py`, `utils.R`, every module CLI | Preserve numerator-minus-denominator semantics | **Implemented** and core parity verified across the current consumer studies. |
| `workflow/lib/evidence.R` | `hypotheses.py`, displayed-data contracts, report evidence sections | Generalize evidence records | **Implemented**, but hypothesis-table parity is not yet verified. |
| `workflow/lib/io.R`, `packages.R` | `utils.R`, module-local checked imports, Python utilities | Consolidate stable I/O and dependency behavior | **Implemented**; retire source helpers after all dependent modules pass. |
| `workflow/lib/viz.R`, `workflow/stages/common/plot_style_utils.R` | `utils.R` and generic constructors in `qc.R`, `de.R`, `pathways.R`, `publication.R` | Replace study styling with configuration-driven themes | **Implemented/partial**; displayed-data and visual-regression gates remain. |
| `workflow/envs/*.yaml` | packaged `workflow/envs/{core,r,network}.yaml`, lock files, OCI image | Consolidate rule-specific environments | **Implemented** environments; lock/fresh-machine/container release gates remain. |
| `workflow/report/analysis.qmd` | `report.py`, `front_door.py`, recipe metadata | Replace study-bound Quarto front door with portable HTML | **Implemented/partial**; link, content, warning, and manifest/report parity gates remain. |

## Input, preflight, counting, QC, and DE mapping

| Source capability | BulkRNAFrame destination | Required preserved behavior | Status |
|---|---|---|---|
| `workflow/stages/common/adapter.py` | `materialize_inputs.py` | Normalize all input boundaries into one canonical contract | **Implemented** for counts, BAM, nf-core/rnaseq, and archive. |
| `validate_study.py`, `validate_hypotheses.py` | layered schemas and `config.py` semantic validation | Fail before expensive work; validate IDs, contrasts, dependencies, and optional story files | **Implemented**. |
| `preflight_bams.py` | `materialize_inputs.py` | BAM existence/readability, sample mapping, headers, indexes, integrity, reference compatibility | **Implemented**; fresh real-BAM integration gate remains. |
| `strand_test.py` | BAM materialization/counting path | Declared or inferred strandedness with auditable evidence | **Implemented/partial**; exact legacy inference parity needs a real fixture/reference gate. |
| `count_genes.py` | `materialize_inputs.py`, `rules/core.smk` | featureCounts options, GTF checksum, assignment metrics, deterministic column mapping | **Implemented**; exact counts passed for current consumers. |
| archive handling added during generalization | `materialize_inputs.py` | Reject traversal/unsafe members; extract supported BAM archives into a controlled directory | **Implemented** and unit tested. |
| historical comparison and QC gate scripts | `verification.py`, `qc.R`, report/manifest warnings | Separate analytical QC from migration comparison | **Implemented/partial**; QC threshold/warning parity still needs explicit table checks. |
| `workflow/stages/qc/stage.R` | `qc.R` | library sizes, detected genes, zero and mitochondrial fractions, expression density, PCA, sample correlations/distances, variable-gene heatmap, QC summaries | **Implemented**; individual QC table/displayed-data parity not yet verified. |
| `workflow/stages/core/deseq2.R` and all `workflow/stages/de/*.R` | `de.R`, `rules/core.smk` | arbitrary design formulas/covariates, explicit contrasts, DESeq2, apeglm shrinkage, full result tables, MA/volcano/distribution/heatmap/PCA figures | **Implemented**; count and DE decisions/numerics are core-parity verified. Figure data/layout parity remains. |
| `limma_sensitivity.R` | sensitivity output in DE/advanced configuration | Preserve optional method-sensitivity analysis | **Partial**; declare as an optional sensitivity module or document retirement after confirming it is not required by any consumer. |
| `release_manifest.py` | `manifest.py` | checksum inputs/results, resolved config, contrast semantics, environments, providers, seeds, warnings, revision | **Implemented/partial**; manifest completeness has not yet been compared with a required-field acceptance fixture. |

## Discovery and biological-context module mapping

### Pathways and ontology

| Source | Destination | Preserved contract | Status |
|---|---|---|---|
| `workflow/stages/pathways/pathway_gene_sets.R` | `resources.R`, provider receipts, custom GMT interface | Resolve species-backed or user-provided gene sets with version/hash provenance | **Implemented/partial**; cached provider tests exist, live release integration and exact snapshot parity remain. |
| `pathway_rankings.R`, `pathway_fgsea.R` | `pathways.R` | deterministic ranking, seeded fgsea, NES direction, FDR, leading edge | **Implemented**; ranking universe/NES/leading-edge parity not yet verified. |
| `pathway_ora.R` | `pathways.R`, `ontology.R` | tested universe, directional up/down membership, adjusted p-values and gene counts | **Implemented**; universe/membership parity pending. |
| `pathway_gsva.R`, `pathway_gsva_tables.R`, `pathway_gsva_plots.R` | `pathways.R`, `publication.R` | per-sample scores, contrast summaries, selected heatmaps, displayed tables | **Implemented**; score and displayed-selection parity pending. |
| remaining pathway I/O/theme/table/plot fragments | `pathways.R`, `utils.R` | complete machine-readable tables plus PDF/PNG/displayed data | **Implemented/partial**; artifact completeness and visual gates pending. |
| `ontology_go.R` | `ontology.R`, GO provider | GO BP enrichment with explicit ID mapping and universe | **Implemented**; mapping/universe parity pending. |
| `ontology_kegg.R` | `ontology.R`, KEGG provider | species-aware KEGG enrichment and displayed data | **Implemented**; mouse snapshot parity and scheduled live tests pending. |
| `ontology_string.R` | `networks.py`, STRING provider | Keep functional enrichment and interaction queries under the STRING provider rather than conflating ontology and networks | **Implemented/partial**; exact mapped/unmapped and term-table parity pending. |
| ontology plots/tables/theme | `ontology.R`, `publication.R` | complete tables and configurable displayed selections | **Implemented/partial**. |

### Composition and cell-state scoring

| Source | Destination | Preserved contract | Status |
|---|---|---|---|
| `workflow/stages/composition/deconvolution_signatures.R` | signature schema, provider/custom signature loader | Validate and map curated signatures | **Implemented** interface; signatures remain consumer-owned. |
| composition I/O/tables/plots/stage | `composition.R` | per-sample relative signature scores, contrast effects, matched-gene counts, audit tables, publication plot | **Implemented/partial**; output must be labeled relative state/composition scoring, never cell fractions. Exact scores and selected states have not yet passed parity. |

### Regulators and GRNs

| Source | Destination | Preserved contract | Status |
|---|---|---|---|
| `regulators_prior.R` and GTRD-backed target sets | `resources.R`, custom regulon edges, provider receipt | Maintain unsigned target-program scoring separately from signed regulation | **Partial**; custom GTRD snapshots are accepted, but a general live GTRD provider is not implemented. |
| `viper.R`, `viper_dorothea.R` | `regulators.R` | signed DoRothEA/VIPER inference, confidence levels, per-sample activity, differential activity | **Implemented/partial**; full activity and regulon-edge parity pending. |
| `grn_dorothea.R` | `grn.py` | export complete retained regulon edges before any display selection; preserve sign, regulator, target, evidence | **Implemented/partial**; node/edge completeness and deterministic-layout parity pending. |
| regulator tables/plots/theme | `regulators.R`, `grn.py`, `publication.R` | activity heatmaps and rectangular/radial GRN alternatives with displayed-data metadata | **Implemented/partial**; figure and audit-table gates pending. |

### STRING and integrated networks

| Source | Destination | Preserved contract | Status |
|---|---|---|---|
| `ontology_string.R` and network input preparation | `resources.R`, `networks.py` | submit complete configured seed breadth; retain mapped and unmapped inputs | **Implemented/partial**; mapping tables need parity checks. |
| `network_integrate.R` | `networks.py`, `multilayer.py` | connected subnetworks, community detection, typed evidence layers, deterministic layout | **Implemented/partial**; breadth, community, and layout parity pending. |
| network tables | `networks.py` | complete nodes/edges, centrality/community fields, direction, evidence type | **Implemented/partial**; centrality breadth from legacy scripts still needs completion. |
| network plots/theme | `networks.py`, `publication.R` | enrichment summary, upregulated and downregulated STRING networks, clean program/community annotations | **Implemented/partial**; displayed node/edge and perceptual gates pending. |

### Hypothesis engine

| Source | Destination | Preserved contract | Status |
|---|---|---|---|
| `workflow/stages/hypothesis/*.R` | `hypotheses.py`, `hypotheses.schema.yaml`, `hypothesis_panels.schema.yaml` | arbitrary claims, expected direction, gene/pathway/program evidence, auditable supporting/contradicting lines | **Implemented/partial**; study claims are consumer-owned and evidence-table parity pending. |

## Publication figure migration map

All generic constructors move into BulkRNAFrame. CAPE names, group labels,
genes, pathways, biological program assignments, program colors, panel letters,
dimensions, and layout decisions stay in consumer YAML/TSV files.

The generic constructor registry and `bulk-rna figures init|catalog|build|gallery`
workflow are now implemented. Constructor/variant/module/contrast validation,
self-contained staged panels, displayed-data copies, hashes, and review galleries
are available; CAPE displayed-data and visual parity remain acceptance gates.

Every final constructor must emit:

- vector PDF;
- review PNG;
- exact displayed-data TSV;
- selection and layout metadata JSON.

SVG and TIFF produced by the exploratory CAPE notebooks are optional export
formats, not v1 guarantees. Notebook HTML previews and LaTeX auxiliary files
are reference-only and must not become package outputs.

| CAPE source function/artifact family | Generic BulkRNAFrame destination | Status/gate |
|---|---|---|
| `make_pca_plot`, PCA requested edit, PCA coordinates/variance tables | `qc.R`: group-colored points, group-colored ellipses, sample labels, displayed coordinates and variance | **Implemented/partial**; exact displayed coordinate and figure-layout parity pending. |
| `make_correlation_heatmap` | `qc.R`: clustered sample correlation with configured condition colors and optional shared legend | **Implemented/partial**. |
| `make_pca_correlation` | `publication.R`/assembly recipe: shared-legend composite | **Implemented/partial**; assembly geometry gate pending. |
| `make_cell_state`, `make_cell_state_hybrid` | `composition.R` and `publication.R`: grouped state effects with matched-gene encoding and FDR annotation | **Implemented/partial**; exact state ordering/data pending. |
| original volcano and `make_volcano` | `de.R` plus recipe-selected original/updated variant | **Implemented**; DE data parity exists, displayed label/geometry parity pending. |
| `make_de_heatmap` and option 1 global clustering | `publication.R`: global row/column clustering and one-column program legend | **Implemented/partial**. |
| option 2 program-grouped clustering | `publication.R`: program blocks with within-program clustering | **Implemented/partial**. |
| option 3 compact direct program labels | `publication.R`: direct program labels without wide legend | **Implemented/partial**. |
| gene-program assignment/audit tables | consumer `hypothesis_panels.yaml`/TSV plus displayed-data output | Interface **implemented**; CAPE assignments remain **consumer-owned**. |
| combined bidirectional GO ORA bubble plot | `pathways.R`/`publication.R`: all configured up/down terms, staggered axes, magnitude gradient, gene-count size | **Implemented/partial**; exact terms/statistics pending. |
| Hallmark/Reactome GSVA heatmap | `pathways.R`/`publication.R`: per-sample or contrast heatmap with configured selections | **Implemented/partial**; scores/selections pending. |
| advanced multi-track GSEA curves | `pathways.R`/`publication.R`: ES curve, hits, ranked metric, NES/FDR and leading edge | **Implemented/partial**; curve coordinates/leading edges pending. |
| integrated program heatmap and effects | `publication.R`: program-colored heatmap plus effect summaries | **Implemented/partial**. |
| consolidated program violins | `publication.R`: one program-grouped figure, program shading, sample points/boxes, tests and significance brackets | **Implemented/partial**; exact tests and annotation-placement parity pending. |
| STRING enrichment summary | `networks.py`: configured down/leading-edge/up enrichment display | **Implemented/partial**. |
| STRING upregulated network | `networks.py`: unpruned audit network plus clean displayed subnetwork/community labels | **Implemented/partial**. |
| STRING downregulated network | `networks.py`: same contract with direction preserved as downregulated in numerator, not relabeled as a different contrast | **Implemented/partial**. |
| regulator activity heatmap | `regulators.R`/`publication.R`: signed activity with condition annotation | **Implemented/partial**. |
| rectangular DoRothEA GRN | `grn.py`: program/community-separated rectangular option | **Implemented/partial**. |
| radial DoRothEA GRN | `grn.py`: compact radial option retaining complete edge audit outside the display | **Implemented/partial**. |
| TeX assemblies for Figures 5 and 6 | `figure_recipe.yaml` plus `assemble.py` | **Implemented/partial**; exact panel dimensions, letters, shared legends, placements, and assembled PDF/PNG parity pending. |
| CAPE Rmd/notebook scripts | consumer migration record and golden outputs | **Reference only** after constructors and displayed data pass. |

## Advanced systems-biology and collection mapping

| Source | Destination | Current disposition |
|---|---|---|
| `sva_check.R` | `sva.R` | **Implemented/partial**: surrogate variables, covariate correlations, and DE sensitivity exist; compare exact summary and decision changes. |
| `wgcna_analysis.R` | `wgcna.R` | **Implemented/partial**: module construction/eigengenes exist; preserve power diagnostics, module membership, composition correlations, and module GO enrichment where required. |
| `wgcna_hub_genes.R` | `wgcna.R` | **Partial**: hub table exists; verify intramodular connectivity definitions and ranked hub parity. |
| `wgcna_hypothesis_overlap.py` | `hypotheses.py` or `multilayer.py` | **Planned/partial**: retain an explicit module-program overlap table rather than burying this evidence in narrative. |
| `mediation_analysis.R` | `mediation.R` | **Implemented/partial**: configurable mediation runs and warnings exist; model/coefficient/interval parity pending. |
| `mediation_analysis_rela_regulator.R` | generalized mediator definitions in `mediation.R` | **Partial**: regulator-mediated analysis must be configuration rather than a RelA-specific script. |
| `mediation_power_analysis.py` | `power.py` | **Implemented/partial**; assumptions and numerical parity pending. |
| `string_centrality.py` | `networks.py` or `multilayer.py` | **Partial**: restore complete within-study and collection-level centrality summaries if they remain required. |
| `multilayer_network_integration.py` | `multilayer.py` | **Implemented/partial**: typed GRN/WGCNA/STRING integration exists; triangulated-gene parity and evidence-line completeness pending. |
| `meta_leave_one_out.py` | `collection.py` | **Implemented/partial**: signed weighted Stouffer, BH, and leave-one-study-out direction stability exist; no real four-study consumer collection run is yet recorded. |
| `meta_grn_candidate_list.py` | `collection.py` plus hypothesis/multilayer output | **Planned/partial**: generalize candidate intersection/ranking, never embed CAPE/RelA gene choices. |
| `export_full_regulon_targets.R` | `regulators.R`/`grn.py` complete regulon-edge outputs | **Implemented/partial**; verify that display pruning never alters the full audit table. |
| generated advanced TSV/JSON outputs | immutable external golden references | **Reference only**; never ship real biological outputs with the tool. |

Small-sample advanced analyses may run when enabled, but every affected table,
report, and manifest must contain a conspicuous suitability warning. A warning
is not a parity failure; a missing warning is.

## Consumer-study migration map

| Consumer | Input boundary and design | Config ownership | Evidence already passed | Evidence still required |
|---|---|---|---|---|
| CAPE primary | aligned BAM/count-derived canonical input; CAPE versus Control | biological repo `studies/cape/*` and publication recipe files | exact counts, DE values/decisions, PDF/PNG format contract | QC tables, enrichment/GSVA/GSEA, composition, regulators, STRING/GRN, hypotheses, all displayed data, assembly geometry, report/manifest |
| CAPE full | same cohort with full profile/modules | biological repo | exact counts, DE values/decisions, PDF/PNG format contract | every publication/full module plus advanced warnings and audit tables |
| Ligation versus sham | aligned BAM boundary and explicit contrast | biological repo `studies/ligation/*` | exact counts, DE values/decisions, PDF/PNG format contract | enrichment/context/advanced/report/manifest parity as configured |
| Old versus young | aligned BAM/count boundary with design `~ sex + age` and explicit age contrast | biological repo `studies/aging/*` | exact counts, covariate-aware DE values/decisions, PDF/PNG format contract | covariate diagnostics and configured downstream parity |
| RelA versus WT | nf-core/rnaseq BAM-output adapter and explicit genotype contrast | biological repo `studies/rela/*` | exact counts, DE values/decisions, PDF/PNG format contract | adapter provenance plus configured downstream parity |

The current source workflow remains runnable until every enabled module in each
consumer has passed its declared gate. A study may be migrated at one profile
without claiming that its full profile has passed.

## Path-level disposition index

This index covers the remaining repository families that are easy to miss in a
capability-only map. Generated files are grouped by contract rather than listed
one image at a time.

| Source path or family | Destination/owner | Disposition |
|---|---|---|
| `analysis/bulk_rna_frame/compare_counts.py` | `src/bulk_rna_frame/verification.py` | Reusable comparison logic moves to the CLI verifier; consumer-specific reference paths stay in the biological repository. **Implemented/partial** pending expanded gates. |
| `analysis/bulk_rna_frame/compare_de.py` | `verification.py` | Same treatment; field tolerances and decision equality are already implemented and core-parity verified. |
| `analysis/bulk_rna_frame/runner.py` | `src/bulk_rna_frame/cli.py` and installed `bulk-rna` command | Replace sibling-repository orchestration with the installed public CLI. **Implemented**; remove helper after clean-checkout migration. |
| `analysis/bulk_rna_frame/export_mouse_gene_sets.R` | `resources.R`, cached provider snapshots, or an explicitly versioned consumer GMT | Generalize resource resolution; retain the script only as a migration-record utility until provider snapshot parity passes. |
| `studies/*/config.yaml` | biological repository | Legacy golden-run configuration. **Consumer-owned/reference only**; freeze until replacement passes. |
| `studies/*/bulk_rna_frame.yaml` and `bulk_rna_frame_full.yaml` | biological repository | v2 consumer configuration. **Consumer-owned** and retained permanently. |
| `studies/*/samples.tsv`, `contrasts.tsv` | biological repository | Cohort/design truth. **Consumer-owned** and retained permanently. |
| `studies/*/hypotheses.yaml`, `hypothesis_panels.yaml`, `figure_recipe.yaml` | biological repository | Story and presentation choices. **Consumer-owned** and retained permanently. |
| `studies/shared/cell_state_signatures.yaml` | biological repository referenced through the generic signatures schema | Curated biology. **Consumer-owned**. |
| `studies/rela_ko_vs_wt/fastq_inventory.tsv` | biological repository/upstream nf-core record | Upstream provenance; BulkRNAFrame consumes nf-core BAM output and does not own FASTQ processing. |
| `studies/migration_references.yaml` | biological repository | External golden-reference registry. **Consumer-owned**; may be read by verification but never packaged as a default. |
| `tests/test_analysis_set.py`, `tests/test_defaults_merge.py` | `tests/test_config.py`, `tests/test_cli.py`, schema fixtures | Preserve behavior through neutral configuration tests. **Implemented**, with consumer integration tests still needed. |
| `tests/r/test_annotate.R`, `test_contrast.R`, `test_evidence.R`, `test_io.R`, `test_viz.R` | Python unit tests plus R script integration/fixture tests in BulkRNAFrame | Port assertions for public contracts, not source-file identity. **Partial**: config/adapter/manifest tests exist; deeper R statistical/figure fixtures should be added. |
| source `tests/test_repository.py` | destination `tests/test_repository.py` | Preserve repository hygiene and package-data assertions. **Implemented**. |
| `analysis/cape_xizhao_edits/*.R`, `*.Rmd` | registered constructors plus biological figure recipe | Constructor registry, story scaffold, build, gallery, staging, and assembly are **implemented**; notebooks remain **reference only** until panel parity passes. |
| `analysis/cape_xizhao_edits/gene_programs.tsv`, `input_provenance.tsv`, `deliverables.tsv` | biological story files, result provenance, and migration records | **Consumer-owned/reference only**. |
| `analysis/cape_xizhao_edits/figures/**/*.tsv` | immutable displayed-data golden references | **Reference only** and required inputs to expanded verification. |
| `analysis/cape_xizhao_edits/figures/**/*.{pdf,png}` | immutable figure golden references | **Reference only** for dimensions/labels/perceptual gates. |
| `analysis/cape_xizhao_edits/figures/**/*.{svg,tiff}` | biological archive | Optional source exports; not required BulkRNAFrame v1 formats. |
| `analysis/cape_xizhao_edits/**/*.nb.html`, previews | biological archive | Review aids only; do not migrate into tool outputs. |
| `analysis/cape_xizhao_edits/final_figure_*/assemble*.tex` | `figure_recipe.yaml` and `assemble.py` | Replace TeX-specific placement with declared recipe geometry after assembly parity. |
| `analysis/cape_xizhao_edits/finalized_panels_pdf/` | biological golden reference/deliverable archive | **Reference only**; selected variants are named by the consumer recipe. |
| `manuscript_figures/systems_biology_analysis/*.R`, `*.py` | advanced scripts and collection modules listed above | **Retire after parity** on a script-by-script basis. |
| `manuscript_figures/systems_biology_analysis/**/*.tsv`, `*.json` | immutable advanced-analysis golden references | **Reference only**; use to verify numerical breadth and warnings. |
| source `results/`, workflow caches, logs, R objects, rendered reports | biological repository/external archival storage | Never migrate as package data. Register golden subsets and checksums rather than copying bulk results into the public tool. |
| Python/R `__pycache__`, notebook auxiliary, TeX `.aux/.log`, temporary previews | neither repository | Generated housekeeping artifacts; ignore and remove only through safe, scoped cleanup. |

When source fragments within a family implement materially different behavior,
their acceptance evidence must be separate even if BulkRNAFrame consolidates
them into one script. Consolidation is an implementation choice, not permission
to drop outputs.

## Verification framework expansion

`bulk-rna verify` currently supports generic TSV comparison and a legacy mode
whose `all` scope checks counts, DE, and PDF/PNG validity. Before it can be used
as the final migration authority, add named gates with machine-readable output:

| Gate | Required comparison |
|---|---|
| canonical inputs | exact counts/sample/contrast keys; annotation field equality or declared mapping; input-manifest options and hashes |
| QC | metric values, flags/warnings, selected samples, PCA coordinates up to sign, correlations/distances, displayed tables |
| DE | field-specific tolerances, identical significance and direction decisions, tested gene universe, shrinkage method |
| fgsea/GSEA | ranking universe/order, pathway universe, NES direction/value tolerance, FDR decisions, hit positions and leading-edge membership |
| ORA/ontology | tested background, input membership, mapped/unmapped genes, term IDs, gene counts, adjusted p-values, displayed selection |
| GSVA/ssGSEA | per-sample score tolerance, contrast summaries, pathway ordering, displayed selection |
| composition | signature definitions/checksums, mapped genes, per-sample scores, effect statistics and selected states |
| regulators | regulon snapshot/checksum, full signed edges, per-sample activity, differential activity, displayed regulators |
| networks/GRN | complete input/mapped/unmapped lists, node and edge sets, sign/evidence types, communities, centrality, deterministic layout seed, displayed subset |
| hypotheses | claim definitions/checksums, evidence lines, support/contradiction classification, displayed panels |
| advanced | SVA/WGCNA/mediation/power/multilayer tables, method settings, suitability warnings |
| publication data | exact displayed-data TSVs and selection metadata before image comparison |
| figures | PDF validity, PNG validity, dimensions, expected labels/legends/panel count, bounded perceptual difference after displayed-data equality |
| assembly | panel identity/variant, placement coordinates, shared legends, letters, page size, assembled PDF/PNG validity |
| report | complete internal links, promoted artifact inventory, contrast semantics, method/resource/warning sections |
| manifest | all declared input/result checksums, resolved configuration, environment/container versions, provider receipts, seeds, warnings, repository revision |

Verification results must distinguish `not_run`, `passed`, `failed`, and
`not_applicable`. A missing enabled-module reference may not silently pass.

## Resource and species mapping

The provider API is species-neutral. Analysis code receives species metadata;
it does not embed mouse constants.

Mouse v1 must have tested mappings for Ensembl/GTF annotation, taxonomy 10090,
mouse MSigDB, `org.Mm.eg.db`, `dorothea_mm`, STRING, GO, and KEGG. Custom GMT
and regulon-edge files remain valid offline/species-independent boundaries.

Each fetched resource must record provider, organism, release, retrieval time,
request parameters, license notice, and checksum. Cached offline reuse must be
deterministic, and missing required offline snapshots must fail clearly.

Human support is the next provider milestone: taxonomy 9606,
`org.Hs.eg.db`, `dorothea_hs`, human MSigDB/STRING/GO/KEGG fixtures, and parity
tests must be added without changing module interfaces. Current schema/provider
acceptance of human metadata is not equivalent to human analytical parity.

## Prioritized gap register

### P0: required before legacy retirement

1. Expand `bulk-rna verify` to all gates listed above.
2. Capture immutable golden references for every enabled module and displayed
   panel in all four studies, with CAPE primary and CAPE full treated separately.
3. Run the expanded verifier and resolve every unexplained difference.
4. Reproduce the selected CAPE Figures 5 and 6 solely from consumer
   configuration, with no CAPE names, genes, pathways, or colors in tool code.
5. Validate report and manifest completeness, including warnings and resource
   receipts.
6. Keep the source workflow intact until these gates pass.

### P1: required for complete current-capability parity

1. Complete legacy WGCNA diagnostics, intramodular connectivity, module GO,
   composition correlation, and hypothesis-overlap outputs where used.
2. Generalize regulator-mediated mediation and meta-GRN candidate ranking.
3. Restore full within-study and cross-study STRING centrality summaries.
4. Record a real collection config and four-study meta-analysis/leave-one-out
   run in the biological repository.
5. Decide explicitly whether limma sensitivity is a supported module or a
   retired exploratory analysis.
6. Exercise all required live providers on schedule and pin/cache accepted
   resource snapshots.

### P2: required for public release milestones

1. Finish human provider fixtures and parity tests.
2. Run fresh-machine tests through locked Conda environments and both Docker
   and Apptainer-compatible OCI execution.
3. Complete public examples, methods/configuration references, license notices,
   citation metadata, semantic versioning, migration guide, and changelog.
4. Publish the Python distribution and versioned OCI image.

## Execution sequence and exit criteria

### Phase 0: freeze and map the baseline

- Preserve the current source repository and results without rewriting its
  dirty biological worktree.
- Record source revisions, configuration, output roots, and golden artifact
  checksums.
- Maintain this document as the capability ledger.

Exit: every reusable source path has a destination or an explicit retirement
decision, and every biological artifact is marked consumer-owned/reference-only.

### Phase 1: complete verification before more migration claims

- Implement named module gates in `verification.py` and CLI scope selection.
- Add fixtures/unit tests for pass, fail, missing-reference, tolerance, sign,
  set equality, warning, and perceptual-comparison behavior.
- Generate structured migration reports for every consumer/profile.

Exit: `all` means all enabled modules and cannot pass on counts/DE/format alone.

### Phase 2: core and enrichment parity

- Close canonical annotation, preflight, strandedness, QC, DE-figure,
  fgsea/GSEA, ORA, GSVA, GO, and KEGG differences.
- Validate cached provider receipts and offline replay.

Exit: all standard-profile gates pass for CAPE, ligation, aging, and RelA.

### Phase 3: biological-context and advanced parity

- Close composition, GTRD/custom targets, DoRothEA/VIPER, STRING, GRN,
  hypothesis, SVA, WGCNA, mediation, power, centrality, and multilayer gaps.
- Run the real collection analysis and leave-one-out checks.

Exit: every configured publication/full module passes, including full audit
tables and visible suitability warnings.

### Phase 4: publication parity

- Express all CAPE biological choices in consumer story and figure recipes.
- Generate the three DE heatmap options and selected option, combined ORA,
  GSVA, advanced GSEA, program effect/heatmap, consolidated violins, STRING
  panels, regulator heatmap, GRN alternatives, and assembled figures.
- Compare exact displayed data before bounded image comparison.

Exit: selected Figures 5 and 6 regenerate from BulkRNAFrame with no study code
inside the package and pass artifact/data/layout checks.

### Phase 5: migrate consumers and retire duplication

- Point each study front door to its BulkRNAFrame v2 configuration and result
  root.
- Preserve migration reports and immutable golden references.
- Tag/archive the last legacy workflow revision before removing duplicated
  reusable workflow code from the biological repository.
- Retain study configs, hypotheses, panel recipes, manuscript text, golden
  references, and migration records.

Exit: a clean consumer checkout can reproduce accepted results through an
installed BulkRNAFrame release, and no manuscript analysis depends on deleted
source code.

### Phase 6: public mouse v1 and human milestone

- Run clean Conda/container installations, scheduled provider checks, and
  public examples.
- Publish mouse v1 artifacts and documentation.
- Add human resource/provider parity without changing analytical contracts.

Exit: mouse v1 satisfies the full release checklist; human support receives its
own declared parity evidence rather than inheriting the mouse claim.

## Safe legacy retirement procedure

1. Freeze/tag the final legacy revision and checksum external golden results.
2. Run both workflows into separate, explicit output roots.
3. Require expanded verifier passes for each study/profile.
4. Switch consumer documentation and automation to an installed BulkRNAFrame
   version, not a sibling-source import.
5. Re-run from a clean checkout/environment.
6. Remove only duplicated reusable workflow code; do not delete biological
   configuration, claims, recipes, results, or migration evidence.
7. Record the BulkRNAFrame version and source tag that replaced each legacy
   module.

## Definition of done

BulkRNAFrame mouse v1 is complete only when all of the following are true:

- all supported input boundaries produce validated canonical inputs;
- exact counts and declared DE equivalence pass for all consumer studies;
- enrichment universes, directions, leading edges, ORA membership, and GSVA
  scores pass;
- regulator and network node/edge completeness, sign semantics, mapped and
  unmapped inputs, communities/centralities, and deterministic layouts pass;
- every selected publication panel has PDF, PNG, exact displayed data, and
  selection/layout metadata;
- the selected CAPE assembled figures regenerate from consumer configuration;
- reports and manifests contain complete links, checksums, environments,
  providers, seeds, contrast semantics, and warnings;
- standard, publication, full, and collection behavior have synthetic tests;
- all four biological consumers pass their enabled-module gates;
- fresh-machine Conda and container execution pass;
- no study-specific biology exists in BulkRNAFrame implementation defaults;
- the biological repository can retire duplicated code without losing
  reproducibility.
