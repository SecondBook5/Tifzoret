# Authoring a study (start here)

This guide is for the scientist who ran the experiment. You do **not** need to
read the code, know Python or R, or understand Snakemake to use this engine. You
describe your experiment and your hypotheses in one plain-text file, run one
command, and get back a publication-grade report with figures and tables.

If you are using the authoring UI, it writes this same file for you — reading
this guide tells you what every choice in the UI means and how to sanity-check
what it produced.

---

## 1. The one rule that governs every result

Every comparison in this engine is written as **numerator vs denominator**, and
a positive number *always* means **higher in the numerator**.

> A log2 fold-change of `+1.5` for gene *X* in the contrast `treated vs control`
> means gene *X* is higher in **treated**. A `-1.5` means it is higher in
> **control**.

The engine never guesses direction from group names, alphabetical order, or the
order of your samples. **You** state it, once, in the contrast table (below), and
every figure, table, and sentence in the report obeys it. If a volcano plot looks
"flipped" from what you expected, the fix is always the contrast definition — not
the plot.

Keep this rule in mind for the rest of the guide; it is the single most common
source of confusion when reading results.

---

## 2. What you provide

A study is a small folder. In the fullest case it holds five files, but three of
them can be folded into the main file so that **the entire study is one
self-contained file** (this is what the UI produces):

```text
my-study/
├── project.yaml        ← the main file: what/where the data is, and what to run
├── samples.tsv         ← one row per sample
├── contrasts.tsv       ← the comparisons you want (defines direction)
├── hypotheses.yaml         ← optional: your biological claims  ┐ can be written
├── hypothesis_panels.yaml  ← optional: the gene/pathway panels ├ directly inside
└── figure_recipe.yaml      ← optional: which figure to build   ┘ project.yaml
```

- **`samples.tsv` and `contrasts.tsv`** are always separate spreadsheets — they
  are tables, and tables belong in tables.
- **`project.yaml`** is the main configuration file, described below.
- The **three publication files** are only needed if you want the hypothesis
  evidence and the assembled multi-panel figure. Each can either sit beside
  `project.yaml` as its own file **or** be written straight into `project.yaml`.
  One self-contained file is the recommended form.

You never edit anything else. You never touch the analysis code.

---

## 3. `samples.tsv` — your samples

A plain tab-separated spreadsheet, one row per sample. The first column must be
`sample_id` (a unique name). Every other column is a property of that sample that
you might want to compare on or correct for — condition, genotype, sex, batch,
timepoint, and so on. Add as many as you like.

```tsv
sample_id	condition	sex	batch
ctrl_1	control	F	1
ctrl_2	control	M	1
ctrl_3	control	F	2
trt_1	treated	F	1
trt_2	treated	M	1
trt_3	treated	M	2
```

If your inputs are aligned BAM files, add a `bam` column pointing at each file.
You can keep several experiments in one sheet and pick one with an optional
`analysis_set` column — but for a single study you can ignore that.

---

## 4. `contrasts.tsv` — the comparisons (and their direction)

This is where you state, explicitly, what is compared with what. One row per
comparison. Four required columns:

```tsv
contrast_id	factor	numerator	denominator
treated_vs_control	condition	treated	control
```

- **`contrast_id`** — a short name you choose; it labels the outputs.
- **`factor`** — which column of `samples.tsv` to compare on (`condition` here).
- **`numerator`** and **`denominator`** — the two levels of that factor. Positive
  results are higher in the **numerator** (see §1).

You can list several contrasts. The factor must appear in your design formula
(next section), and both levels must actually be present in the samples you're
analyzing.

---

## 5. `project.yaml` — the main file

This describes your data and what analyses to run. You do not write it from
scratch — run `bulk-rna init my-study --input counts` (or `--input bam`) to get a
filled-in scaffold, then edit the values. The important sections, in plain terms:

```yaml
project:
  id: my_study                 # a short, file-safe name
  title: My treated-vs-control study

species: {provider: mouse, scientific_name: Mus musculus, taxonomy_id: 10090}
reference: {genome_build: GRCm39, annotation_release: 107}

inputs:
  kind: counts                 # 'counts' if you have a count matrix; 'bam' for aligned reads
  counts: counts.tsv           # your gene-by-sample count matrix
  annotation: annotation.tsv   # gene IDs → symbols  (for 'bam' inputs, give bam_root + gtf instead)
  samples: samples.tsv

analysis:
  design: "~ condition"        # the statistical model; add covariates like "~ batch + condition"
  contrasts: contrasts.tsv
  profile: publication         # how much to run — see below
  random_seed: 1               # makes every run reproducible

figures:
  group: condition             # which sample property colours the plots
  palette: {control: "#A6CEE3", treated: "#F4A6A6"}   # one colour per level
  de: {fdr: 0.05, abs_log2fc: 1.0}                    # significance thresholds

output: {root: ../results}
```

Two things worth understanding:

- **`design`** is a statistical model written in R's formula style. `~ condition`
  compares conditions directly. `~ batch + condition` compares conditions *after*
  correcting for a batch effect. List the things you want to account for; put the
  factor you're contrasting last.
- **`profile`** decides how much runs:
  - `standard` — QC, differential expression, pathways, ontology, and a report.
  - `publication` — everything in `standard` plus cell-state scoring, regulator
    activity, networks, your hypothesis evidence, and the assembled figure.
  - `full` — everything in `publication` plus advanced sensitivity analyses.

  Start with `standard` to see your data, move to `publication` for the figures.

---

## 6. Stating a hypothesis (publication profile)

This is the heart of the engine: you write your biological claim in words, list
the genes and pathways that would be evidence for it, and the engine tests it and
reports how the evidence fell. This lives in `hypotheses.claims` (a file, or
inline in `project.yaml`):

```yaml
hypotheses:
  claims:
    study: my_study
    hypotheses:
      - id: emt_activation
        statement: >
          Treatment induces an epithelial-to-mesenchymal transition — loss of
          epithelial identity genes and gain of mesenchymal genes.
        contrast: treated_vs_control       # which comparison tests this
        expected_direction: increased_emt  # what you predict
        gene_panels: [epithelial, mesenchymal]   # the evidence gene sets
        pathway_panels: [emt_pathways]
        regulators: [emt_drivers]
```

The **gene panels** referenced here are defined once in `hypotheses.panels` (again
a file or inline) — a named list of genes grouped biologically:

```yaml
hypotheses:
  panels:
    gene_panels:
      epithelial:
        description: Epithelial identity markers.
        groups:
          identity: [CDH1, EPCAM, KRT8, KRT18]
      mesenchymal:
        description: Mesenchymal markers gained during EMT.
        groups:
          identity: [CDH2, VIM, FN1, ACTA2]
```

You write the biology; the engine handles the statistics, the multiple-testing
correction, the figures, and the wording — all obeying the direction rule from §1.

---

## 7. Running it

Each command is one line you type in a terminal. Run them in order the first time;
afterwards you only re-run the steps you need.

```bash
bulk-rna validate my-study/project.yaml      # check every file is consistent before doing work
bulk-rna dry-run  my-study/project.yaml      # show the plan without running it
bulk-rna run      my-study/project.yaml --cores 4   # run the full analysis
bulk-rna report   my-study/project.yaml      # build the browsable HTML report
```

For the publication profile, also:

```bash
bulk-rna figures build   my-study/project.yaml --cores 4   # render the panels
bulk-rna figures gallery my-study/project.yaml             # side-by-side variant comparison
```

`validate` is your friend — it catches typos, a colour missing from the palette, a
contrast level that isn't in the samples, or a gene panel referenced but never
defined, and it tells you exactly which file and line to fix. Nothing heavy runs
until `validate` is clean.

---

## 8. Reading the output

Everything lands under `results/<project.id>/<analysis_set>/`:

- **`REPORT.html`** — open this in a browser. It's the guided tour: QC, the
  differential-expression results, pathway and ontology enrichment, your
  hypothesis verdicts, and every figure, with captions.
- **`figures/`** — each panel as both a vector **PDF** (for the manuscript) and a
  **PNG** (for quick viewing).
- **`tables/`** — the displayed numbers behind every figure, as spreadsheets.
- **`manifest.json`** — a complete provenance record: your exact settings, data
  checksums, tool versions, and random seeds. This is what makes a result
  reproducible and defensible for a reviewer.

When you read a fold-change, a volcano plot, or a "higher in …" sentence, apply
§1: positive means higher in the numerator you named in `contrasts.tsv`.

---

## 9. Common mistakes (and what `validate` says)

| Symptom | Cause | Fix |
|---|---|---|
| Plot looks "flipped" | numerator/denominator swapped | edit the row in `contrasts.tsv` |
| "level not found in samples" | a typo in `numerator`/`denominator` | match the exact spelling in `samples.tsv` |
| "palette missing level X" | a sample group has no colour | add every level of `figures.group` to `palette` |
| "panel referenced but not defined" | a hypothesis lists a gene panel that isn't in `hypotheses.panels` | define it, or fix the name |
| "factor not in design" | contrast factor absent from `analysis.design` | add it to the formula |

---

## 10. Where to go next

- [configuration.md](configuration.md) — every field of `project.yaml`, in full.
- [figures.md](figures.md) — the figure recipe: which panels exist and how to
  arrange them into a multi-panel figure.
- [methods.md](methods.md) — the statistical methods, in the words you'd put in a
  paper's Methods section.
- [troubleshooting.md](troubleshooting.md) — when a run fails.
