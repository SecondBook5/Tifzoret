# nf-core/rnaseq result scaffold

This scaffold begins downstream of nf-core/rnaseq. Replace the example rows in
`samples.tsv`, set `NFCORE_RNASEQ_ROOT` and `BULK_RNA_GTF` (or replace those
values in `project.yaml`), and adjust `bam_pattern` to the selected nf-core
aligner output. BulkRNAFrame will uniformly recount the located BAMs.
