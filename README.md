# BulkRNAFrame

BulkRNAFrame is a configuration-driven, end-to-end downstream bulk RNA-seq workflow. It accepts aligned BAMs, nf-core/rnaseq BAM outputs, ZIP/TAR archives containing BAMs, or validated integer count matrices. FASTQ processing remains upstream in nf-core/rnaseq.

The stable v1 target is mouse. The analysis contracts are species-neutral, and the human provider implementation is already exercised at the schema boundary while parity fixtures are completed.

## What it runs

- canonical counting and GTF-first annotation;
- full sample QC, DESeq2 with apeglm shrinkage, and explicit numerator-minus-denominator contrasts;
- cached MSigDB, GO BP, KEGG, STRING, and DoRothEA provider stages;
- fgsea, bidirectional ORA, GSVA/ssGSEA, and multi-track GSEA curves;
- relative cell-state signature scoring, regulator activity, STRING networks, and GRNs;
- configurable hypothesis evidence and publication panels;
- vector PDF/raster PNG figure assembly and a navigable HTML report;
- optional SVA, WGCNA, mediation/power, multilayer integration, and collection meta-analysis;
- complete checksummed provenance and numerical run verification.

Cell-state results are relative signature/composition scores, not estimated cell fractions. Regulatory, co-expression, and STRING association edges remain explicitly distinguished.

## Install and start

```bash
python -m pip install -e '.[workflow]'
bulk-rna init my-study --input counts
bulk-rna validate my-study/project.yaml
bulk-rna prepare my-study/project.yaml --cores 4
bulk-rna dry-run my-study/project.yaml
bulk-rna run my-study/project.yaml --cores 4
bulk-rna report my-study/project.yaml
bulk-rna figures init my-study/project.yaml
bulk-rna figures build my-study/project.yaml --cores 4
bulk-rna figures gallery my-study/project.yaml
```

Input scaffolds are available for `counts`, `bam`, `nfcore-rnaseq`, and `archive`.

## Project contract

```text
my-study/
├── project.yaml
├── samples.tsv
├── contrasts.tsv
├── hypotheses.yaml          # publication/full profiles
├── hypothesis_panels.yaml   # publication/full profiles
└── figure_recipe.yaml       # publication/full profiles
```

Positive effects always mean `numerator - denominator`. No rule infers direction from group names or sample order.

Profiles:

- `standard`: inputs, QC, DE, pathways, ontology, report, and manifest;
- `publication`: standard plus composition, regulators, networks, hypotheses, publication panels, and assembly;
- `full`: publication plus SVA, WGCNA, mediation/power, and multilayer integration.

Individual modules can be overridden under `analysis.modules`. Small-sample advanced analyses continue when enabled, but warnings are written into summaries, reports, and the release manifest.

## Commands

```text
bulk-rna init
bulk-rna validate
bulk-rna dry-run
bulk-rna prepare
bulk-rna run
bulk-rna report
bulk-rna assemble
bulk-rna figures init|catalog|build|gallery
bulk-rna verify
bulk-rna migrate-config
bulk-rna collection validate|run
```

`prepare` stops at the canonical inputs. `verify` compares candidate TSV artifacts with an external reference using declared tolerances. `assemble` reads `figure_recipe.yaml` and creates a vector multi-panel PDF, review PNG, and placement metadata.

```bash
bulk-rna verify my-study/project.yaml --reference /path/to/golden/results --scope counts
bulk-rna verify my-study/project.yaml --reference /path/to/golden/results
```

Legacy migration mode recognizes the established featureCounts and DESeq2
cache layout. It requires exact counts, uses field-specific DE tolerances,
requires identical threshold/direction decisions, and validates PDF/PNG pairs.

## Output contract

```text
results/<project>/<analysis_set>/
├── REPORT.html
├── manifest.json
├── figures/                         # promoted review-facing figures + metadata
├── tables/                          # promoted review-facing displayed data
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
├── publication/<figure_set>/assembled/
└── .cache/{logs,resources}/
```

Every rendered constructor declares PDF, PNG, and displayed-data artifacts. The release manifest records normalized configuration, inputs/results checksums, contrast semantics, provider receipts, random seeds, warnings, environment/container information, and repository revision.

See [configuration](docs/configuration.md), [architecture](docs/architecture.md), [methods](docs/methods.md), [hypothesis-driven figures](docs/figures.md), [migration](docs/migration.md), the [full source-to-tool migration map](docs/full-migration-map.md), and [troubleshooting](docs/troubleshooting.md).
