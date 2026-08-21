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
