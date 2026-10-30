# Troubleshooting

Failures are grouped by the stage that raises them. Most contract violations are
caught by `tifzoret validate` (exit code `2`, message on stderr) before any rule
runs — run it first.

## Configuration and validation

- **Unset variable.** Set every `${NAME}` referenced by `project.yaml`;
  validation names each unresolved variable. Export it in the environment or
  define it under the config's variable block.
- **Schema violation.** `validate` reports the offending key and the expected
  type/enum. The machine-readable schemas live in `src/tifzoret/schemas/`; the
  key reference is [configuration.md](configuration.md).
- **Profile vs. module conflict.** A `publication`/`full` profile expects the
  companion files (or their inline equivalents); a missing `publication.recipe`
  is reported when you `assemble` or `figures build`.

## Inputs and quantification

- **Archive error.** BAM member paths are relative to `member_root`; absolute
  paths, symlinks, and `..` traversal are rejected. Repackage with clean
  relative paths.
- **"BAM and GTF share no chromosome names".** `validate` compares the contig names in the
  first BAM header with the sequence names in the GTF, because a BAM aligned to UCSC-style names
  (`chr1`) counted against an Ensembl-style GTF (`1`) assigns no reads. The message says whether
  they differ only in a `chr` prefix. Use an annotation built for the same reference as the BAMs,
  or add or remove the prefix in the GTF (for example `awk 'BEGIN{FS=OFS="\t"} !/^#/{$1="chr"$1} 1'`
  for Ensembl-to-UCSC, keeping `MT` as `chrM`).
- **A rerun counted again, or did not.** The raw featureCounts result is cached under
  `<results>/.cache/counting/`, keyed on the BAMs (path, size, modification time), the GTF, the
  `counting` settings and the featureCounts version. Editing anything else in `project.yaml`
  (design, contrasts, modules, figures, hypotheses) reuses it and finishes in seconds. Changing
  any keyed input, or touching a BAM, recounts. Set `TIFZORET_RECOUNT=1` to force a recount.
- **featureCounts failed.** The stage surfaces featureCounts' own stdout/stderr
  and the strand mode it was running. Confirm the BAMs are coordinate-sorted,
  the annotation GTF matches the alignment genome, and paired-end/strandedness
  settings are declared correctly.
- **Strand inference looks wrong.** Fix strandedness explicitly in the counting
  config rather than relying on inference when forward/reverse rates are close.
- **No abundance emitted.** Validated count-matrix inputs carry no exon lengths,
  so TPM/FPKM are intentionally absent — start from BAMs if you need abundance.

## Resources and providers

- **Offline resource error.** Run once with network access to populate the
  cache, or point `resources.cache` at an already-populated cache. Provider
  snapshots are then reused deterministically.
- **No gene sets remain.** Verify symbol species/case, the collection selection,
  and the `min_size`/`max_size` bounds — an over-tight size window can empty the
  set list.
- **STRING returns nothing.** The network stage calls the STRING REST API; check
  connectivity and that input symbols map to the configured species. Unmapped
  inputs are reported rather than dropped silently.
- **SPIA/pathway-topology skipped.** SPIA needs `SPIA`, `graphite`, the
  species-matched `org.*.eg.db`, and a reachable KEGG topology database; when any
  is missing the module writes a documented skip instead of failing the run.

## Differential expression and contrasts

- **Coefficient cannot be resolved.** Confirm the contrast factor is in the
  design formula and the denominator is a real, selected level. For factorial
  `coefficient` contrasts, the named interaction term must exist after
  releveling.
- **Shrinkage prior rejected.** `shrinkage` must be one of `apeglm`, `ashr`,
  `normal`, `none`; other values are rejected at validation.
- **Downstream module skipped for a contrast.** `coefficient`/`omnibus`
  contrasts and weighted estimands (interactions, main effects, composites) run
  DE, pathways and ontology. The other two-group modules (composition,
  regulators, networks, hypotheses, publication, SPIA, enrichment maps, SVA,
  WGCNA) expand over pairwise contrasts only; a multi-factor family yields no
  pairwise rows, so they do not run for it.

## Small-sample and exploratory modules

- **Small-sample warning.** The module ran as requested, but its result should
  remain exploratory. The warning is machine-readable and retained in the
  report and the manifest — it is not an error.
- **variancePartition produced empty outputs.** Expected when no covariate is
  usable; the stage emits well-formed empty results rather than failing.

## Environments and execution

- **Missing executable (`snakemake` not found).** Exit code `127`. Install the
  workflow extra (`pip install -e '.[workflow]'`) or pass `--snakemake` with an
  explicit path.
- **Conda solve fails / dependency missing under `--no-conda`.** The single-env
  `environment.yaml` is a complete superset of the per-rule envs; recreate it, or
  drop `--no-conda` to use the pinned per-rule locks.
- **Relative paths don't resolve.** Workflow commands run from the config file's
  directory. Keep `samples.tsv`, `contrasts.tsv`, and data pointers relative to
  `project.yaml`, exactly as `validate` resolves them.

## Migration and verification

- **Migration/verify mismatch.** Run `tifzoret prepare` first and confirm exact
  counts, then compare DE and displayed-data tables with `tifzoret verify`.
  Verification uses field-specific tolerances (`--atol`/`--rtol`) and requires
  identical threshold/direction decisions; a mismatch returns exit code `1` and
  is itemized in the report file.
- **Wrong scope compared.** Narrow with `--scope counts` to isolate the counting
  boundary before comparing `core` or `all`.
