# Normalized-expression program observations

`src/tifzoret/program_observations.py` provides `normalized_program_observations(expression, annotation, programs, *, scale)` for explicitly labeled CPM or TPM matrices. This observation-only boundary computes equal-gene mean log1p abundance scores and gene-coverage records. It does not normalize libraries, import integer counts, run differential expression, estimate cell fractions or modify the primary pipeline DAG.

Inputs must contain finite nonnegative abundances and unique gene/library identifiers. Each prespecified target must map uniquely to a measured gene; missing and ambiguous targets are rejected rather than silently dropped. Programs must not overlap. Preserve source processing, annotation provenance and sample identities in the calling study.

The first consumer is a downstream temporal-program analysis that lives in its own study repository. Effects and uncertainty belong to the calling study. Scores are RNA observations, not direct measurements of biochemical activity or kinetic states.
