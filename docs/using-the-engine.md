# Using the engine

This document answers the practical operator questions: how to configure inputs, skip or choose analysis modules, resume from partial runs, and understand where outputs land.

## What is the input?

Tifzoret defines **one canonical boundary**: four validated TSV files (counts, samples, contrasts, annotation). Every upstream source is normalized to this boundary by the input adapter, so all downstream stages are blind to the original format.

**Four input adapters are provided:**

1. **`counts`** — You already have validated integer counts and annotation. You provide:
   - `inputs.counts` — path to counts TSV (genes × samples)
   - `inputs.annotation` — path to annotation TSV (gene metadata)
   - `inputs.samples` — path to samples TSV (sample metadata with design covariates)
   - `analysis.contrasts` — path to contrasts TSV (declared comparisons)

2. **`bam`** — You have aligned BAM files. You provide:
   - `inputs.bam_root` — directory containing BAMs
   - `inputs.gtf` — path to GTF annotation
   - `inputs.samples` — samples TSV with a `bam` column (BAM filename per sample)
   - `counting` block (threads, feature type, strandedness, etc.)

   The adapter runs `featureCounts` on each BAM, calculates exon lengths, derives TPM/FPKM, and parses the GTF to build the annotation table.

3. **`nfcore_rnaseq`** — You ran nf-core/rnaseq and have its output tree. You provide:
   - `inputs.root` — nf-core output directory
   - `inputs.gtf` — GTF annotation
   - `inputs.bam_pattern` — template for BAM paths (e.g., `star_salmon/{sample_id}/{sample_id}.markdup.sorted.bam`)
   - `inputs.samples` — samples TSV (metadata columns interpolated into `bam_pattern`)
   - `counting` block

   The adapter resolves each sample's BAM from the pattern and proceeds as with `bam`.

4. **`archive`** — You have a ZIP or TAR.GZ archive containing BAMs. You provide:
   - `inputs.archive` — path to archive
   - `inputs.gtf` — GTF annotation
   - `inputs.samples` — samples TSV with a `bam` column (archive-relative BAM path)
   - Optional `inputs.member_root` — prefix to strip from archive member paths
   - `counting` block

   The adapter extracts BAMs on demand and runs `featureCounts`.

**Result:** Regardless of which adapter you use, the workflow produces the same four canonical files under `results/<project>/<analysis_set>/inputs/`. Everything downstream sees the same interface.

## How do I change the input?

**Option 1: Edit the `input` block in `project.yaml`.**

Only the `input.boundary` key and its associated paths change. The rest of the configuration (design, contrasts, profiles, modules) remains identical. For example, switching from `counts` to `bam`:

```yaml
input:
  boundary: bam
  bam_root: data/bams
  gtf: resources/genome.gtf
  samples: metadata/samples.tsv

counting:
  threads: 4
  feature_type: exon
  gtf_attribute: gene_id
  paired_end: true
  strandedness: infer
  # ... rest of counting config
```

**Option 2: Re-scaffold with `tifzoret init`.**

To generate a fresh template for a different input boundary:

```bash
tifzoret init my-study --input counts
tifzoret init my-study --input bam
tifzoret init my-study --input nfcore-rnaseq
tifzoret init my-study --input archive
```

Each scaffold contains the correct `input` block structure and placeholder paths. Copy your existing `analysis`, `resources`, and `figures` blocks into the new scaffold.

## How do I skip or choose analysis modules?

Tifzoret uses **profiles** to control which stages run. Each profile is a superset of the one before:

- **`standard`** — Core analysis only: QC, differential expression, pathway enrichment, ontology views.
- **`publication`** — Adds composition, regulators, networks, GRN, hypotheses, and figure assembly.
- **`full`** — Adds exploratory advanced modules: SVA, WGCNA, curvature, multilayer, mediation/power.

Set the profile in `project.yaml`:

```yaml
analysis:
  profile: publication
```

**Fine-grained control:** Override individual modules with `analysis.modules`. Each module is a boolean flag:

```yaml
analysis:
  profile: publication
  modules:
    # Disable GRN radial figure in the publication profile
    grn_radial: false
    # Enable WGCNA even though we're not on full profile
    wgcna: true
```

Available module flags (see [`configuration.md`](configuration.md) for the full list):
- `de_confirm`, `sva`, `batch`, `variance_partition` — QC and confirmatory DE
- `pathways`, `ontology`, `spia`, `enrichment_map` — enrichment layers
- `composition`, `deconvolution` — cell-state scoring
- `regulators`, `grn`, `grn_radial` — TF activity and GRN
- `networks`, `wgcna`, `curvature`, `multilayer` — network layers
- `mediation`, `power` — causal analysis
- `consensus`, `hypotheses` — synthesis

**Conditional stages:** Some stages activate automatically when their trigger data is present:
- `study_batch` — runs only if `samples.tsv` has a `batch` column
- `contrast_de_confirm` — runs only if `analysis.modules.de_confirm: true`
- `contrast_omnibus` — runs only for contrasts with `type: omnibus` in `contrasts.tsv`
- `contrast_spia` — runs only if `analysis.modules.spia: true` and KEGG is available
- `study_deconvolution` — runs only if `composition.signature_matrix` is provided

If the trigger isn't present, Snakemake skips the stage without error.

## How do I run part of the workflow or resume from a checkpoint?

**Run only to canonical inputs:**

```bash
tifzoret prepare project.yaml --cores 4
```

This stops after materializing the input boundary. You can inspect `results/<project>/<analysis_set>/inputs/` before committing to the full analysis.

**Dry-run to see the plan:**

```bash
tifzoret dry-run project.yaml
```

Prints the dependency graph and planned shell commands without executing. Use this to verify which stages will run.

**Resume from cached intermediates:**

Snakemake tracks all intermediate outputs. If you run the workflow, it fails partway, and you fix the problem:

```bash
tifzoret run project.yaml --cores 4
```

Snakemake resumes from the last successful checkpoint. Already-completed stages are skipped; only failed or downstream stages re-run.

**Run only to the HTML report:**

```bash
tifzoret report project.yaml --cores 4
```

Targets `REPORT.html`. Runs any missing upstream stages, but stops before publication figure assembly (if configured).

**Run only publication figures:**

```bash
tifzoret assemble project.yaml --cores 4
# or for a specific figure set:
tifzoret figures build project.yaml --figure-set Fig1 --cores 4
```

Targets the configured `publication.recipe`. Runs any missing upstream stages required by the recipe's panel constructors.

## Where do outputs land?

All outputs land under the resolved result root: `results/<project_id>/<analysis_set>/`. The directory structure mirrors the ten-phase organization:

```
results/<project_id>/<analysis_set>/
├── inputs/
│   ├── counts.tsv
│   ├── samples.tsv
│   ├── contrasts.tsv
│   ├── annotation.tsv
│   └── input_manifest.json
├── qc/
│   ├── pca.pdf
│   ├── sample_distances.pdf
│   └── qc_summary.json
├── batch/  (if batch column present)
├── <contrast_id>/
│   ├── de/
│   │   ├── de_results.tsv
│   │   ├── volcano.pdf
│   │   └── ma_plot.pdf
│   ├── pathways/
│   ├── ontology/
│   ├── composition/
│   ├── regulators/
│   ├── grn/
│   ├── networks/
│   └── hypotheses/
├── consensus/  (study-wide)
├── publication/
│   └── <figure_set>/
│       ├── panels/
│       │   └── <panel_id>_<variant>.pdf
│       └── <figure_name>.pdf
├── REPORT.html
└── manifest.json
```

**Key paths:**

- **Canonical inputs:** `results/<project>/<analysis_set>/inputs/*.tsv`
- **Study-wide QC:** `results/<project>/<analysis_set>/qc/*.pdf`
- **Per-contrast DE results:** `results/<project>/<analysis_set>/<contrast>/de/de_results.tsv`
- **Pathway enrichment:** `results/<project>/<analysis_set>/<contrast>/pathways/*.tsv`
- **GRN radial figure:** `results/<project>/<analysis_set>/<contrast>/grn/grn_radial.pdf`
- **Assembled publication figures:** `results/<project>/<analysis_set>/publication/<figure_set>/*.pdf`
- **HTML report:** `results/<project>/<analysis_set>/REPORT.html`
- **Release manifest:** `results/<project>/<analysis_set>/manifest.json` (checksummed provenance)

**Changing the result root:**

By default, `output.root: results` in `project.yaml`. You can override this to write results elsewhere:

```yaml
output:
  root: /scratch/analysis_results
```

Then results land in `/scratch/analysis_results/<project_id>/<analysis_set>/`.

**Optional front-door layout:**

If you configure `publication.front_door: true`, the workflow copies reviewer-critical outputs to a flat `front_door/` directory:

```
results/<project>/<analysis_set>/front_door/
├── REPORT.html
├── manifest.json
├── Fig1.pdf
├── Fig2.pdf
└── ...
```

This simplifies external access — reviewers don't need to navigate the full result tree.

## How do I verify results match a reference?

Use `tifzoret verify` to compare a candidate run against a reference run within numerical tolerances:

```bash
tifzoret verify project.yaml \
  --reference results_reference/<project>/<analysis_set> \
  --candidate results_candidate/<project>/<analysis_set> \
  --scope core \
  --output verification.json
```

**Scopes:**
- `counts` — verifies only canonical inputs
- `core` — verifies inputs + QC + DE
- `all` — verifies the entire result tree

**Tolerances:**
- `--atol` — absolute tolerance (default 1e-8)
- `--rtol` — relative tolerance (default 1e-6)

Writes a structured verification report and exits 0 if all checks pass, 1 otherwise. See [`docs/migration.md`](migration.md) for using this to validate that a migrated project reproduces legacy results.

## Common CLI patterns

**Validate configuration without running:**

```bash
tifzoret validate project.yaml
```

Checks config validity, resolves all paths, and prints a JSON summary. Exits 2 on validation error.

**List available figure constructors:**

```bash
tifzoret figures catalog
# or machine-readable:
tifzoret figures catalog --json
```

**Initialize hypothesis-driven publication files:**

```bash
tifzoret figures init project.yaml
```

Scaffolds `hypotheses.yaml`, `hypothesis_panels.yaml`, and `figure_recipe.yaml` beside the project.

**Build a panel variant review gallery:**

```bash
tifzoret figures gallery project.yaml --output review.html
```

Produces an HTML contact sheet with all built panel variants for side-by-side comparison.

**Run a multi-study collection (meta-analysis):**

```bash
tifzoret collection validate collection.yaml
tifzoret collection run collection.yaml
```

See [`docs/configuration.md`](configuration.md) for collection schema details.
