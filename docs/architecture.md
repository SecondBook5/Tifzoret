# Architecture

## Design rule

The engine may discover and quantify evidence, but it must never silently
invent a biological story. Automated results and curated story choices are
separate inputs with separate provenance.

## Layers

The packaged Snakefile resolves configuration and output contracts, then loads
separate `core`, `providers`, `modules`, `advanced`, `publication`, and `report`
rule groups. Collection analyses operate across already materialized project
results through the collection CLI rather than nesting project DAGs.

### Input adapters

Adapters convert supported upstream products into the canonical contract:

- integer counts keyed by `gene_id`;
- sample metadata keyed by `sample_id`;
- `gene_id` to `gene_symbol` annotation;
- named contrasts with explicit numerator and denominator.

The primary adapter accepts aligned BAMs, validates them, resolves or infers
strandedness, and runs featureCounts with declared settings. The
`nfcore_rnaseq` adapter treats nf-core/rnaseq as an upstream FASTQ-to-BAM
producer and locates its BAMs without reimplementing upstream processing. The
`archive` adapter safely extracts BAMs from ZIP or TAR archives before applying
the same validation and counting contract. The `counts` adapter is an explicit
bypass for previously quantified data.

All adapters materialize the same canonical files beneath the result directory.
Downstream rules never branch on the upstream source type.

### Statistical core

The core owns sample validation, normalization, QC, design matrices, DESeq2,
and directionality. Every signed effect is numerator minus denominator.

### Discovery modules

Discovery modules consume canonical DE and expression contracts. Custom GMT
gene sets are species-independent. Species-backed providers resolve MSigDB,
GO, KEGG, STRING, and regulator resources through explicit organism metadata,
cache each response, and emit provenance receipts. Offline execution requires
the corresponding cache entries rather than silently substituting resources.

### Figure library

Figure constructors consume tables, not hidden global state. Every constructor
emits PDF, PNG, and a displayed-data TSV. Figure selection thresholds and
layout values are configuration, not source-code constants.

### Story recipes

Optional recipes define curated gene programs, program colors, highlighted
pathways, network anchors, panel variants, and figure assembly. They are not
required for exploratory analysis.

### Release layer

The release manifest records project configuration, input and result hashes,
tool versions, and contrast semantics. BAM provenance uses file metadata,
whole-file checksums, samtools header checksums, integrity status, the GTF
checksum, counting options, and assignment rates. Resource receipts record the
provider, organism, release, request, retrieval time, license notice, and hash.

## Initial repository relationship

`lymphatic-flow-homeostasis` remains a biological project and golden-reference
consumer. Code is promoted into this repository only after its inputs and
outputs have been expressed through neutral contracts.
