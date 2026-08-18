# Publication-Grade Consolidation of bulk-rna-frame

- **Date:** 2026-08-18
- **Status:** Design — awaiting user review before writing-plans
- **Author:** SecondBook5 (with Claude)
- **Repos in scope:** `bulk-rna-frame` (the package/engine, primary), `lymphatic-flow-homeostasis` (the driver repo hosting the Cape and PTEN studies)

## 1. Goal & context

`bulk-rna-frame` is being prepared as a publication-grade bulk RNA-seq engine
(BAM → featureCounts → DESeq2 → publication figures, any design). The CAPE
thoracic-duct study is going to *Nature Communications*; the package itself
should read as if it were the subject of a methods software paper — tight,
consolidated, reproducible, and defensible.

The engine is **already the consolidated core**. It has a single shared
primitives file (`utils.R`: theme, dual PDF+PNG save, direction resolver,
palette), every stage reads `contrasts.tsv`, and the paper's figures already
exist as engine **recipes**: `figure_recipe.yaml` builds Figure 1 (panels A–I)
and Figure 2 (panels A–E) from 20 registered panel constructors, and
`figure_1.pdf` / `figure_2.pdf` already assemble for Cape.

The remaining work is therefore **not** to build figures or re-consolidate
engine internals. It is to:

1. **Prove and close panel-level parity** between the engine's assembled
   figures and the exact panels that are in the paper
   (`analysis/cape_xizhao_edits/`), fixing any gap *in the engine constructor*.
2. **Retire the sprawl** in the driver repo that duplicates the engine
   (4 Cape configs, 26 ad-hoc `manuscript_figures/*` scripts, the
   `cape_xizhao_edits` working scripts, the superseded old pipeline).
3. **Generalize the same figure system to PTEN and other datasets.**
4. **Harden reproducibility, professionalism, and scientific defensibility** to
   methods-paper bar.

### The scope decision (from the user)

- Cape is the **primary** use, but the package must handle **other datasets and
  more complex requests** — the PTEN 2×2 factorial study is the forcing example.
- **The figures in `analysis/cape_xizhao_edits/` are the *only* figures in the
  Nature Communications paper**, plus **QC figures**. That set — and nothing
  else — is the figure deliverable.
- All 14 engine modules (core + the 6 advanced: regulators, networks, sva,
  wgcna, mediation, multilayer) are **kept**. Tightness means killing config and
  script sprawl and reproducibility gaps, not analytical capability.

### Honesty constraint (non-negotiable)

Any paper-facing material must be grounded in actual results. The PTEN study's
pre-registered "enhanced antitumor immunity" hypothesis has **zero** supporting
evidence lines; PTEN figures anchor on the proliferation-vs-ECM rewiring story
(`studies/lec_pten_tumor/findings_for_xizhao.md`), never the aspirational
hypothesis.

## 2. The deliverable — exact figure inventory

The paper's figures are enumerated by `analysis/cape_xizhao_edits/`. Each maps
1:1 to an engine constructor already registered in the catalog. The parity gate
(§5.1) is defined against exactly this table.

### Figure 1 — nine panels (A–I)

| Panel | Reference (`finalized_panels_pdf/`) | Engine constructor | Variant |
|---|---|---|---|
| A | PCA and Sample Correlation | `pca_correlation` | default |
| B | Cell-State Signatures | `cell_state_effects` | default |
| C | **Original Volcano Plot** | `volcano` | **5-class (see §5.1)** |
| D | Global Clustering / Program Grouped / Direct Program Labels | `de_heatmap` | `global_clustered`, `program_grouped`, `direct_program_labels` |
| E | GO Biological Process ORA | `go_ora` | default |
| F | Hallmark GSVA Heatmap | `gsva_heatmap` | default |
| G | Curated Program GSEA | `gsea_multitrack` | default |
| H | Integrated Heatmap and Effects | `program_heatmap_effects` | default |
| I | Consolidated Violins | `program_violins` | default |

### Figure 2 — five panels (A–E)

| Panel | Reference | Engine constructor | Variant |
|---|---|---|---|
| A | STRING functional enrichment | `string_enrichment` | default |
| B | STRING network (up) | `string_network` | `upregulated` |
| C | STRING network (down) | `string_network` | `downregulated` |
| D | Regulator activity | `regulator_activity` | default |
| E | DoRothEA GRN | `dorothea_grn` | `radial` |

### QC figures (paper-grade, likely supplementary)

The engine already emits eight QC figures at publication bar and must keep them
first-class: `pca`, `sample_correlation`, `pca_correlation`, `library_metrics`,
`expression_density`, `sample_distance`, `variable_gene_heatmap`, `qc_overview`.
(`pca_correlation` doubles as Figure 1 Panel A.)

Output format is **PDF + PNG per panel** (`save_plot_pair` already does exactly
this). Individual panels are the deliverable; the LaTeX assembly is convenience,
not a requirement.

## 3. Non-goals

- No removal of any of the 14 modules.
- No rewrite of engine internals that already work (DESeq2 core, apeglm
  shrinkage, the `resolve_contrast` direction model, the recipe/assembly system).
- No touching the unrelated `nfkb/` (rela / NF-κB) study sprawl — different study.
- No new figure types beyond the inventory in §2.
- No Nextflow migration (engine stays Snakemake).

## 4. Architecture — the engine as the single source

The engine already expresses the right architecture; the consolidation makes it
the *only* path and documents it as the stable public contract:

- **A study = three files**: `bulk_rna_frame.yaml` (config), `contrasts.tsv`
  (the contrasts, pairwise and/or coefficient), and `figure_recipe.yaml` (which
  named constructors compose which panels). Nothing else is needed to go from
  BAMs to paper figures.
- **The constructor catalog is the stable API.** `bulk-rna figures catalog`
  lists the 20 constructors and their variants. New datasets author a recipe
  against this catalog; they do not write figure code. This is what "handles
  other datasets / more complex requests" means concretely.
- **Complex designs** (factorial) are already first-class: a coefficient
  contrast (e.g. `genotypePTEN.tumoryes` from `~ genotype * tumor`) runs DE +
  pathways; the `PAIRWISE_CONTRAST_IDS` gate scopes the two-group modules to
  pairwise contrasts. This is the mechanism that lets one package serve both a
  simple 2-group study (Cape) and a 2×2 factorial (PTEN).

## 5. Workstreams

### 5.1 Panel-parity gate (the centerpiece)

Produce a **parity report** — one row per panel in §2 — comparing the engine's
assembled panel against the `cape_xizhao_edits` reference (content, palette,
typography, layout). Every gap is fixed **in the engine constructor**, never in
a side script. The reference styling library is
`cape_xizhao_edits/publication_figure_functions.R`; the old pipeline
(`workflow/stages/`) holds additional ground-truth styling that the paper
actually used and must be mined before it is retired.

**Confirmed gap #1 — the volcano (Panel C).** Panel C is documented as "unchanged
volcano plot from the validated primary analysis output" and uses the old
pipeline's **5-class** `SIGNIFICANCE_PALETTE`
(`workflow/stages/de/de_packages.R:129`):

| Class | Meaning | Color |
|---|---|---|
| `significant_up` | FDR<0.05 & log2FC ≥ +1 | `#B22222` |
| `significant_down` | FDR<0.05 & log2FC ≤ −1 | `#2166AC` |
| `padj_only` | FDR<0.05 but \|log2FC\| < 1 | `#2A9D8F` |
| `lfc_only` | \|log2FC\| ≥ 1 but not FDR-significant | `#F4A261` (the amber "yellow") |
| `ns` | neither | `#C7CDD4` |

The engine's `volcano` constructor is currently **3-class** (up / not-significant
/ down, blue-red) and collapses `padj_only` and `lfc_only` into "not
significant."

**Confirmed gap #2 — the MA plot.** The old pipeline's MA plot
(`workflow/stages/de/de_plot_main.R:902`) uses the *same* 5-class
`SIGNIFICANCE_PALETTE`. The MA plot is not itself a Figure 1/2 paper panel, but
the engine emits `ma.pdf` and the combined `de_overview.pdf`, so it must share
the volcano's classification for internal consistency.

**Fix (both gaps, one change):** introduce a single shared 5-class significance
classifier + palette in the engine (thresholds from `figures.de.fdr` /
`figures.de.abs_log2fc`), and have both the `volcano` and `ma` constructors — and
therefore `de_overview` — consume it. One palette definition, DRY, matching the
paper's DE styling exactly, and it must continue to work for PTEN and any other
dataset.

**Acceptance:** every panel row in §2 is marked `MATCH` in the parity report, or
`MATCH-after-fix` with the engine change referenced. No panel remains sourced
from an ad-hoc script.

### 5.2 Retire the sprawl

- **Cape configs 4 → 1.** Collapse `config.yaml` (stale old-pipeline),
  `bulk_rna_frame.yaml`, `bulk_rna_frame_full.yaml`, and
  `bulk_rna_frame_publication.yaml` into a single canonical
  `bulk_rna_frame.yaml`. The `analysis_set` (primary / full) and profile
  variants become a documented switch (config field + CLI override), not three
  near-duplicate files. `config.yaml` is deleted (belongs to the retired old
  pipeline).
- **Retire ad-hoc figure scripts** (`manuscript_figures/cape/*`,
  `cape_xizhao_edits/*_requested_edit.R`, `figure2_*_panel_*.R`) **only after**
  the corresponding panels pass the §5.1 parity gate. Keep the final reference
  PDFs (`finalized_panels_pdf/`, `final_figure_1/`, `final_figure_2/`) plus a
  one-page provenance note as the archival record.
- **Old pipeline** (`workflow/stages/`, driver `workflow/Snakefile`): mine for
  paper-used styling (§5.1), then retire per `docs/bulkrnaframe_migration.md`.
  **Deletion is the last step and is gated on parity** — see Decision D1.
- Leave `nfkb/` untouched.

### 5.3 Generalize to PTEN and other datasets

- PTEN's `figure_recipe.yaml` already exists and mirrors Cape's structure
  (Figure 1 A–I, Figure 2 A–E), anchored on `pten_vs_wt_tumor`. Verify it
  assembles end-to-end and passes the same panel-quality bar (including the
  5-class volcano from §5.1).
- Enforce the honesty guard: PTEN figure titles/anchors reflect the
  proliferation-vs-ECM rewiring story; the interaction contrast
  (`pten_x_tumor_interaction`, de + pathways only) supplies the
  difference-in-differences evidence.
- Document the "author a recipe for a new dataset" path (config + contrasts +
  recipe against the constructor catalog) so the generalization is a documented
  capability, not tribal knowledge.

### 5.4 Reproducibility & correctness

- **Direction convention — consolidate to one code path.** Verified consistent
  today: `de.R` (via `resolve_contrast`), `composition.R:48`, and
  `regulators.R:101` all relevel to `denominator` and read the `numerator`
  coefficient (positive = up in numerator); `networks.py` consumes DE's own
  `up_in_numerator` / `down_in_numerator` labels. This is **not** an active bug
  (unlike the old pipeline). The task is to route `composition.R` and
  `regulators.R` through the shared `resolve_contrast()` that `de.R`/`pathways.R`
  use, eliminating three inline re-derivations — identical results, one code
  path. Defensibility/DRY closure.
- **Verify-and-fix in the engine** (these were old-pipeline audit findings; their
  status in `bulk-rna-frame` is unconfirmed and must be checked, then fixed if
  present): deterministic `set.seed` before every Monte-Carlo enrichment (fgsea,
  gseGO/gseKEGG, any permutation ORA); QC metrics (library size, detected genes,
  zero/mt fraction) computed on the **raw pre-filter** matrix; ORA/GO background
  universe = all tested genes, not only non-NA-padj genes.
- **Pinned environment + version manifest.** Pin the conda/R env used for the
  manuscript run (lockfile or pinned `envs/*.yaml`) and record tool versions in a
  run manifest alongside the figures.
- **Methods ↔ code provenance.** The methods text and the engine parameters
  (thresholds, shrinkage, GSEA/ORA settings) must derive from one source so they
  cannot drift.

### 5.5 Software professionalism

Docstrings/type hints on the Python surface (`figures.py`, `assemble.py`,
`config.py`, `materialize_inputs.py`, providers); a pytest suite with meaningful
coverage of config validation, contrast resolution (pairwise + coefficient), and
panel resolution; lint + CI; a clean README and docs so the workflow reads
end-to-end; release hygiene (version, `CITATION.cff`, MIT `LICENSE` — already
added). Clarify the two-repo story in docs (engine = package; studies = driver).

### 5.6 Scientific defensibility

Audit that FDR control, apeglm shrinkage (done), DE thresholds, and GSEA/ORA
settings are consistent across stages and documented in methods. Confirm the
background-universe and seeding fixes from §5.4 land. Confirm the coefficient
(interaction) contrast is correctly specified and reported.

## 6. Decisions

- **D1 — Old pipeline: mine-then-retire, not delete-blind.** The old pipeline
  holds ground-truth styling that is *in the paper* (the 5-class volcano). It is
  mined for paper-used styling, ported into engine constructors, parity is
  proven, and only then is it deleted. Deletion is the last step, gated on the
  §5.1 parity report. *(Revised from the initial "delete now" recommendation
  after the volcano finding.)*
- **D2 — Spec and engine hardening live in `bulk-rna-frame`; driver-repo cleanup
  is tracked here and executed in `lymphatic-flow-homeostasis`.** The package
  published as a methods paper is the engine; the Nature Comms analysis is the
  Cape study in the driver repo. Two repos, one spec.

## 7. Definition of done

1. The engine reproduces every panel in §2 for Cape, each marked MATCH in the
   parity report; a single shared 5-class significance scheme is implemented and
   consumed by the `volcano` and `ma` constructors (and `de_overview`).
2. The same recipe machinery produces the PTEN figures end-to-end, on the honest
   story, at the same bar.
3. Cape has exactly one config file; the retired ad-hoc scripts and old pipeline
   are gone (post-parity), with reference PDFs archived.
4. Direction resolution is a single code path; seeds/QC-matrix/background-universe
   are verified/fixed; the manuscript run is reproducible from a pinned env with
   a version manifest.
5. Engine has tests + CI + docs such that a reader can go BAM → paper figures by
   authoring only config + contrasts + recipe.

## 8. Risks & open questions

- **Panel parity beyond the volcano.** Other panels may have subtle reference
  deviations (typography, legend placement, program colors). The §5.1 report
  will surface them; each is an engine-constructor fix. Risk: a reference panel
  used a bespoke transform not expressible by the current constructor — handled
  by extending the constructor, not by re-introducing a side script.
- **Config consolidation and existing results paths.** Collapsing analysis_set
  variants must not silently change the canonical `results/` layout that current
  artifacts live under; the switch must preserve `primary` / `full` separation.
- **Old-pipeline deletion breadth.** Confirm nothing still-referenced (docs,
  Makefile targets) depends on `workflow/stages/` before deletion.
- **Reproducibility findings may be no-ops.** If seeds/QC-matrix/background are
  already correct in the engine, those tasks close as verified — reported
  honestly, not padded.

## 9. Repo/scope split

| Work | Repo |
|---|---|
| Volcano 5-class + any constructor parity fixes | `bulk-rna-frame` |
| Direction consolidation, seeds, QC matrix, background universe | `bulk-rna-frame` |
| Pinned env, version manifest, tests, CI, docs | `bulk-rna-frame` |
| Cape config 4→1, retire ad-hoc scripts, retire old pipeline | `lymphatic-flow-homeostasis` |
| PTEN recipe verification + honesty guard | `lymphatic-flow-homeostasis` |
| Parity report (references live in driver repo) | `lymphatic-flow-homeostasis`, informing engine fixes |
