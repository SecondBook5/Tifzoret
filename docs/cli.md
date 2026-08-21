# CLI reference

Tifzoret's public interface is the `tifzoret` command plus the `project.yaml`
configuration contract. There is no Python import API; everything a study needs
is driven through the CLI below and the keys documented in
[configuration.md](configuration.md). The registered figure constructors — the
figure "API" — are listed by `tifzoret figures catalog` and described in
[figures.md](figures.md).

```text
tifzoret <command> [options]
```

Run `tifzoret <command> --help` for the argparse-generated summary of any
command.

## Exit codes

| Code | Meaning |
|------|---------|
| `0` | Success. |
| `1` | `verify` ran and the comparison did not pass within tolerances. |
| `2` | Configuration/contract violation (`validation error:` on stderr). |
| `127` | A required executable (e.g. `snakemake`) was not found. |
| other | For workflow commands, Snakemake's own exit code is passed through unchanged. |

## Common workflow flags

Every command that drives the Snakemake workflow (`prepare`, `dry-run`, `run`,
`report`, `assemble`, `figures build`) accepts:

| Flag | Default | Meaning |
|------|---------|---------|
| `--cores N` | `1` | Cores handed to Snakemake. |
| `--snakemake PATH` | `snakemake` | Snakemake executable to invoke. |
| `--no-conda` | off | Skip per-rule conda provisioning and run in the current environment (the single-env `environment.yaml` path). Omit it to use `--use-conda` with the pinned per-rule locks. |

These commands always validate the project first, then run Snakemake from the
configuration file's directory so relative paths resolve exactly as they do
under `validate`.

## Lifecycle commands

### `init`

```text
tifzoret init <directory> [--input counts|bam|nfcore-rnaseq|archive]
```

Scaffold a new project directory from an input-boundary template (`--input`
default `counts`). The destination must be empty.

### `validate`

```text
tifzoret validate <project>
```

Validate the project and every tabular input, printing a resolved JSON summary.
Exits `2` on any contract violation.

### `prepare`

```text
tifzoret prepare <project> [--cores N] [--snakemake PATH] [--no-conda]
```

Run only up to the canonical inputs (counting/annotation) and stop, targeting
`inputs/input_manifest.json`, so the data boundary can be inspected before the
full analysis.

### `dry-run`

```text
tifzoret dry-run <project> [--cores N] [--snakemake PATH] [--no-conda]
```

Resolve and print the Snakemake DAG without executing any rule.

### `run`

```text
tifzoret run <project> [--cores N] [--snakemake PATH] [--no-conda]
```

Execute the full analysis for the resolved profile, targeting the release
`manifest.json`.

### `report`

```text
tifzoret report <project> [--cores N] [--snakemake PATH] [--no-conda]
```

Build the navigable `REPORT.html`, running any missing upstream steps.

### `assemble`

```text
tifzoret assemble <project> [--cores N] [--snakemake PATH] [--no-conda]
```

Assemble the configured multi-panel publication figure(s). Requires
`publication.recipe` (a file or inline mapping).

## Figure commands

### `figures init`

```text
tifzoret figures init <project> [--force]
```

Scaffold `hypotheses.yaml`, `hypothesis_panels.yaml`, and `figure_recipe.yaml`
beside the project and validate the recipe. Existing files are preserved unless
`--force` is given.

### `figures catalog`

```text
tifzoret figures catalog [--json]
```

List the registered figure constructors, their contrast-vs-study scope, and
available variants. `--json` emits the machine-readable catalog.

### `figures build`

```text
tifzoret figures build <project> [--figure-set NAME] [--cores N] [--snakemake PATH] [--no-conda]
```

Render the recipe's panels and assemble its figure set(s). With `--figure-set`
only that set is built; without it, every configured set is assembled. Requires
`publication.recipe`.

### `figures gallery`

```text
tifzoret figures gallery <project> [--output PATH]
```

Build an HTML + contact-sheet review of every built panel variant so a reviewer
can compare variants before choosing the final panel. `--output` overrides the
destination.

## Verification

### `verify`

```text
tifzoret verify <project> --reference REF [--candidate DIR] \
  [--scope counts|core|all] [--output verification.json] \
  [--atol 1e-8] [--rtol 1e-6]
```

Compare a project's result against a reference run within declared tolerances.

| Flag | Default | Meaning |
|------|---------|---------|
| `--reference` | required | Reference (golden) result directory. |
| `--candidate` | resolved result | Override the candidate directory checked. |
| `--scope` | `all` | Restrict the comparison: `counts`, `core`, or `all`. |
| `--output` | `verification.json` | Report path. |
| `--atol` | `1e-8` | Absolute tolerance. |
| `--rtol` | `1e-6` | Relative tolerance. |

Returns `0` when the comparison passes, `1` otherwise.

## Migration

### `migrate-config`

```text
tifzoret migrate-config <project> --output OUT \
  [--species mouse|human|custom] [--scientific-name NAME] [--taxonomy-id INT] \
  [--genome-build STR] [--annotation-release STR] [--force]
```

Convert a development version-1 configuration into the public version-2
contract, applying the requested species/reference metadata (which v1 predates)
and re-validating the written file in place. Refuses to overwrite an existing
`--output` unless `--force`.

## Collections

### `collection validate`

```text
tifzoret collection validate <collection>
```

Validate a multi-study collection and print its resolved study summary.

### `collection run`

```text
tifzoret collection run <collection>
```

Run the cross-study (meta-analysis) collection over already materialized project
results and report the output path.

## A typical session

```bash
tifzoret init my-study --input bam           # scaffold
# edit project.yaml, samples.tsv, contrasts.tsv
tifzoret validate my-study/project.yaml      # check the contract
tifzoret prepare  my-study/project.yaml --cores 8   # materialize counts, inspect
tifzoret run      my-study/project.yaml --cores 8   # full analysis
tifzoret report   my-study/project.yaml      # REPORT.html
tifzoret figures init    my-study/project.yaml      # publication companions
tifzoret figures build   my-study/project.yaml --cores 8
tifzoret figures gallery my-study/project.yaml      # review variants
```
