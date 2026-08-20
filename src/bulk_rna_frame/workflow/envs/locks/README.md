# Frozen environment locks

The `*.yaml` files one level up are the human-readable, version-pinned Conda
specifications that Snakemake solves under `--use-conda`. The `*.lock.txt`
files here are the **byte-exact, fully transitive** solutions of those specs —
one `conda create --file`-installable line per package, with exact build strings
and checksums.

| Lock | Solved from | Packages | Role |
|------|-------------|----------|------|
| `core.lock.txt`    | `../core.yaml`    | 90  | BAM validation, featureCounts, assembly, front-door, report |
| `network.lock.txt` | `../network.yaml` | 152 | Python network/GRN figure rendering (networkx + scipy) |
| `r.lock.txt`       | `../r.yaml`       | 414 | DESeq2, fgsea, GSVA, ComplexHeatmap, decoupleR/dorothea/viper, WGCNA |

## Platform

All three locks are `linux-64`. They are the exact environments that produced
the reference figure set and are the authoritative record of what the published
numbers were computed with. They are **not** portable to other platforms —
regenerate per platform (see below).

## How they were produced

Captured from the Snakemake-solved environments after a successful end-to-end
run (Snakemake writes each solved env under the run's
`.snakemake/conda/<hash>_/`):

```sh
conda list -p <prefix> --explicit > <name>.lock.txt
```

Nothing here is hand-authored — every pin is the real solved version.

## Recreating an environment exactly

```sh
conda create --name bulk-rna-frame-r --file r.lock.txt
```

## Relationship to the YAML specs

The YAML specs are the source of truth Snakemake consumes; edit those to change
the environment, then re-solve and refresh the matching lock here. The lock is
the frozen witness of one solve — if a YAML pin and its lock disagree, the YAML
was changed without refreshing the lock.
