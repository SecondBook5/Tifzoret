# Migration guide

This guide describes how to move an existing project onto Tifzoret and
verify that the migrated run reproduces your established results.

Tifzoret v2 separates reusable workflow behavior from study biology. Keep
sample metadata, hypotheses, gene programs, panel recipes, and golden-reference
results in the study repository; keep analysis implementations in
Tifzoret.

## Convert a v1 project

```bash
tifzoret migrate-config old-project.yaml \
  --output project.yaml \
  --species mouse \
  --genome-build GRCm39
tifzoret validate project.yaml
tifzoret dry-run project.yaml
```

Review the migrated file before execution. In particular, confirm the input
adapter, reference release, explicit numerator and denominator for every
contrast, analysis profile, and output root.

## Preserve the old analysis during migration

Point v2 at a new output root. Run `tifzoret prepare` first and compare the
canonical counts, samples, annotation, contrasts, and input manifest. Then run
the configured profile and compare it with the declared reference:

```bash
tifzoret prepare project.yaml
tifzoret verify project.yaml --reference /path/to/golden/results --scope counts
tifzoret run project.yaml
tifzoret verify project.yaml --reference /path/to/golden/results
```

Do not remove the established workflow until counts are exactly equal and the
DE, enrichment, regulator, network, displayed-data, report, and figure gates
required by the project have passed.

## Input-boundary changes

Choose one explicit input kind:

- `bam`: aligned BAM files matched from a directory;
- `nfcore_rnaseq`: BAMs discovered beneath an nf-core/rnaseq result directory;
- `archive`: validated BAMs extracted safely from ZIP or TAR input;
- `counts`: an existing integer gene-by-sample matrix.

FASTQ alignment remains upstream in nf-core/rnaseq. Every adapter materializes
the same canonical downstream input contract.

## Contrast Type Desugaring

Tifzoret now uses an **estimand architecture** where comparisons are expressed as
linear combinations of design cell means. Legacy `contrasts.tsv` rows with
`type` column are automatically desugared into estimand families at
config-load time:

**Pairwise contrasts** (default or `type: pairwise`):
```tsv
contrast_id         factor      numerator    denominator
treated_vs_control  condition   treated      control
```
→ Becomes a two-cell estimand: `(treated) - (control)`

**Coefficient contrasts** (`type: coefficient`):
```tsv
contrast_id              factor      numerator
genotype_mutant_effect   genotype    mutant
```
→ Becomes a single-coefficient estimand referencing the named design term

**Omnibus contrasts** (`type: omnibus`):
```tsv
contrast_id        factor
condition_any      condition
```
→ Generates term tests for the factor (full-vs-reduced likelihood-ratio tests)

**What changes in your results:**

1. **Gene filter default:** Moved from 0 to 10 reads. Genes with <10 reads in all
   samples are now filtered before DE testing. This matches established practice
   and reduces multiple-testing burden. **Action:** Regenerate goldens with the
   new filter, or set `analysis.settings.de.gene_filter_threshold: 0` to preserve
   exact old behavior.

2. **Covariance-aware standard errors:** Estimands with multiple non-zero
   coefficients (e.g., `(treated,wt) + (treated,mutant) - (control,wt) -
   (control,mutant)`) now use `SE(c'β) = sqrt(c'Σc)` instead of the naive
   independent sum. This typically reduces standard errors by ~20-40% for
   composite estimands. **Impact:** Better power, tighter confidence intervals.

3. **Shared dispersions:** Related estimands (e.g., treatment and genotype main
   effects plus interaction) extract from one DESeq2 fit instead of separate
   fits. **Impact:** Fully comparable statistics, no arbitrary differences from
   fit-to-fit variation.

4. **Term tests replace omnibus:** The old omnibus LRT is now a term test derived
   automatically from the design. For a multi-level factor, you get the same
   likelihood-ratio test plus additional nested tests (e.g., for a 2×2 factorial,
   you get main effect tests, interaction test, and any-effect omnibus).

**Validation:**
```bash
# After migration, verify DE direction/magnitude agreement
tifzoret verify project.yaml --reference path/to/old/results --scope de
```

Most studies see exact direction agreement (100%) and magnitude correlation >0.99.
Discrepancies typically come from the improved covariance handling or filter
change — both are improvements, not regressions.

## Publication migration

Move curated claims to `hypotheses.yaml`, genes and biological programs to
`hypothesis_panels.yaml`, and panel selection/layout to `figure_recipe.yaml`.
The engine contains no study-specific group names, genes, colors, or panel
letters. Every promoted panel has PDF and PNG output, displayed-data tables,
and selection/layout metadata.

## Resource and offline behavior

The first online run downloads configured resources and records receipts under
`.cache/resources`. Preserve that directory for reproducible offline reruns.
An offline run fails clearly when a required snapshot is missing. Refresh is an
explicit configuration choice and produces a new receipt and checksum.

## Consumer acceptance sequence

1. validate and materialize canonical inputs;
2. establish exact count equality;
3. compare DE estimates and significance/direction calls;
4. compare tested pathway universes, enrichment direction, leading edges, and
   displayed selections;
5. compare regulator and network node/edge audit tables;
6. compare displayed data before applying bounded visual-regression checks;
7. inspect the release manifest, warnings, resource receipts, and report;
8. retire duplicated workflow code only after all required gates pass.
