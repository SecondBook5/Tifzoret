# Native gene-level RSEM boundary

The `rsem` input kind imports one `.genes.results` file per biological sample
using the actual tximport package. It is independent of study biology.

```yaml
inputs:
  kind: rsem
  samples: samples.tsv
  annotation: annotation.tsv
  counts_from_abundance: "no"
  min_total_count: 0
  count_threshold_operator: ge
  exclude_gene_regex: null
  zero_length_policy: error
```

`samples.tsv` requires `sample_id` and `rsem`, plus declared biological design
columns. RSEM paths resolve relative to project.yaml, with tilde/environment
expansion. Files and genes must be unique; every file must share ordered gene
identifiers and finite nonnegative quantification values. Annotation must cover
all source gene IDs. The adapter deliberately supports gene-level RSEM only.
`scaledTPM` and `lengthScaledTPM` require transcript-level estimates and are
rejected at this boundary.

## Contract and normalization

Preparation emits canonical integer counts and a typed `tximport.rds` object
containing the original estimated counts, imported TPM abundance, effective
length matrices, and countsFromAbundance semantics. Additional TSVs expose the
estimated counts, lengths, TPM and feature-retention decisions. The import
receipt records actual R/tximport versions and options; the input manifest hashes
the source files and additional quantification artifacts.

QC and family fits use a shared checked constructor calling
DESeqDataSetFromTximport. Family subsets retain their aligned length matrices.
Missing bundles, changed canonical counts and identity/order mismatches fail
explicitly. Integer count/BAM/archive paths retain their original constructor.
Existing fitted-object consumers (estimands, LRTs and SVA) inherit the fitted
object's normalization information.

Filtering is explicit at import: rounded row sums are compared using `gt` or
`ge` against `min_total_count`; an optional regex excludes declared features.
The default `zero_length_policy: error` rejects retained zero lengths. An explicitly
justified `one` policy replaces exactly zero with one for model construction;
original lengths and a cell-level adjustment audit are exported separately. Defaults impose no count
threshold or exclusions. Set a family filter explicitly when the imported gene
universe must be preserved. Input filtering is distinct from DESeq2 independent
filtering and from a family's tested-gene filter.

Count-only edgeR confirmation and deconvolution are currently rejected for RSEM
studies. They require separately tested adapters; they must not silently drop
length information. This restriction does not affect existing count/BAM studies.

## Commands

```
tifzoret init my-rsem-study --input rsem
tifzoret validate my-rsem-study/project.yaml
tifzoret prepare my-rsem-study/project.yaml --cores 1
```

The scaffold contains synthetic quantifications for contract testing; replace
them with the study's actual data and metadata. `prepare` fits no DE models.
A dedicated rsem.yaml runtime supplies Python, tximport and jsonlite. Its YAML
pins are provided; a byte-exact conda lock is not yet supplied. Existing lockfiles
are preserved. Import receipts record the actual runtime used, which may differ
from the environment specification in preprovisioned `--no-conda` execution.

Boundary-only studies may leave the custom GMT empty when pathways, ontology
and hypotheses are disabled. Explicit estimand families may use a header-only
legacy contrasts.tsv; the resolved synthesized contrasts drive the DAG.

## Validation boundaries

Tests compare normalized expression from the shared constructor to an
independent direct tximport/DESeq2 import with deliberately varying effective
lengths. They also demonstrate that plain rounded-count normalization differs,
exercise family subsetting, reject corrupted inputs and check DAG dependencies.
External golden validation belongs to the consumer study: input agreement does
not establish complete DE reconstruction or equivalence of package versions.

Method reference:
https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
