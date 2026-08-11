# Architecture

## Design rule

The engine may discover and quantify evidence, but it must never silently
invent a biological story. Automated results and curated story choices are
separate inputs with separate provenance.

## Layers

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
`counts` adapter is an explicit bypass for previously quantified data.

All adapters materialize the same canonical files beneath the result directory.
Downstream rules never branch on the upstream source type.

### Statistical core

The core owns sample validation, normalization, QC, design matrices, DESeq2,
and directionality. Every signed effect is numerator minus denominator.

### Discovery modules

Discovery modules consume canonical DE and expression contracts. Custom GMT
gene sets are species-independent and therefore form the MVP pathway backend.
Species-backed providers are adapters layered on top later.

### Figure library

Figure constructors consume tables, not hidden global state. Every constructor
emits PDF, PNG, and a displayed-data TSV. Figure selection thresholds and
layout values are configuration, not source-code constants.

### Story recipes

Optional recipes will define curated gene programs, program colors, highlighted
pathways, network anchors, panel variants, and figure assembly. They are not
required for exploratory analysis.

### Release layer

The release manifest records project configuration, input and result hashes,
tool versions, and contrast semantics. BAM provenance uses file metadata,
samtools header checksums, integrity status, the GTF checksum, counting options,
and assignment rates without rehashing multi-gigabyte BAM bodies on every run.

## Initial repository relationship

`lymphatic-flow-homeostasis` remains a biological project and golden-reference
consumer. Code is promoted into this repository only after its inputs and
outputs have been expressed through neutral contracts.
