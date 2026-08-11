# BulkRNAFrame

BulkRNAFrame is a configuration-driven downstream bulk RNA-seq workflow. It
starts from aligned BAM files, standardizes counting and metadata, performs
reproducible differential-expression and pathway analyses, and exports both
publication-ready figures and the exact data displayed in them.

nf-core/rnaseq remains an upstream FASTQ-to-BAM workflow. BulkRNAFrame can point
directly at an nf-core/rnaseq result bundle, locate its BAMs with a declared
pattern, and then enter the same downstream counting contract as any other BAM
study. A validated integer count matrix is also supported as an explicit bypass
when quantification has already been completed.

## Supported input boundaries

- `bam`: sample metadata names BAMs beneath a configured root directory;
- `nfcore_rnaseq`: sample metadata and a path pattern locate BAMs in an
  nf-core/rnaseq result directory;
- `counts`: an integer gene-count matrix and gene annotation bypass counting.

For BAM-based inputs, BulkRNAFrame checks BAM integrity and sort order, can infer
strandedness, runs featureCounts using fully declared settings, and records the
source and counting provenance. It does not duplicate read QC, alignment, or
other upstream responsibilities already handled by nf-core/rnaseq.

## Current scope

- arbitrary sample and condition names;
- arbitrary design covariates supported by DESeq2 formulas;
- multiple named contrasts from one project;
- contrast direction standardized as numerator minus denominator;
- PCA with group-colored ellipses and a clustered correlation heatmap;
- DESeq2 tables, volcano plots, and globally clustered DE heatmaps;
- combined up/down ORA bubble plots;
- preranked fgsea tables and multi-track enrichment curves;
- per-sample ssGSEA/GSVA scores and pathway heatmaps;
- declared PDF and PNG outputs plus the exact displayed-data tables;
- checksummed release manifests.

Species-backed external-resource providers such as MSigDB, KEGG, STRING, and
DoRothEA remain planned adapters rather than hidden assumptions in the core.

## Quick start

```bash
python -m pip install -e .
bulk-rna init my_project
bulk-rna validate my_project/project.yaml
bulk-rna prepare my_project/project.yaml --cores 4
bulk-rna dry-run my_project/project.yaml
bulk-rna run my_project/project.yaml --cores 4
```

An executable synthetic count-matrix project is included as the built-in
template. Running `bulk-rna init <directory>` copies it to a new working
directory.

```bash
bulk-rna init synthetic_project
bulk-rna validate synthetic_project/project.yaml
bulk-rna dry-run synthetic_project/project.yaml
bulk-rna run synthetic_project/project.yaml --cores 2
```

For real data, create a source-specific scaffold and then replace its declared
paths and example metadata:

```bash
bulk-rna init my_bam_study --input bam
bulk-rna init my_nfcore_study --input nfcore-rnaseq
```

### BAM input

```yaml
inputs:
  kind: bam
  bam_root: /data/aligned_bams
  samples: samples.tsv
  gtf: /references/genes.gtf
counting:
  threads: 8
  feature_type: exon
  attribute: gene_id
  paired_end: true
  count_read_pairs: true
  require_both_ends_aligned: true
  exclude_chimeric_fragments: true
  strandedness: infer
  strand_test_modes: [1, 2]
  strand_min_dominance: 0.80
```

The sample sheet includes a path relative to `bam_root`:

```text
sample_id  bam                         condition
control_1  control_1.sorted.bam        control
treated_1  batch_2/treated_1.bam       treated
```

### nf-core/rnaseq result input

```yaml
inputs:
  kind: nfcore_rnaseq
  root: /data/nfcore-rnaseq/star_salmon
  samples: samples.tsv
  gtf: /references/genes.gtf
  bam_pattern: "{sample_id}.markdup.sorted.bam"
```

`bam_pattern` may reference any column in `samples.tsv`. The same `counting`
block shown above is required, making nf-core and delivered-BAM studies
comparable at the downstream boundary.

## Output contract

```text
results/<project_id>/
├── inputs/
│   ├── counts.tsv
│   ├── samples.tsv
│   ├── annotation.tsv
│   └── input_manifest.json
├── qc/
│   ├── figures/
│   ├── tables/
│   └── objects/
├── contrasts/<contrast_id>/
│   ├── figures/
│   ├── tables/
│   └── objects/
└── manifest.json
```

Every figure is a declared Snakemake output. Each panel has a matching table
containing the exact values displayed in the figure. The input manifest records
which adapter ran, selected samples, BAM metadata, GTF checksum, resolved strand
mode, counting settings, and tool versions; the release manifest checksums the
complete result bundle.

`bulk-rna prepare` stops after creating the canonical input contract. This is
useful when migrating an established analysis because its counts can be checked
against the previous workflow before any downstream results are accepted.

See [architecture](docs/architecture.md) and [migration plan](docs/migration.md).
