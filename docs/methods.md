# Methods

This reference describes what each analysis stage computes and the tool behind
it, in the order the workflow runs. It is written to seed a manuscript methods
section: every stated default is the engine's actual default, and every named
tool is one the workflow actually invokes. Cite the specific method/database
packages for the modules a study enables (see [`CITATION.cff`](../CITATION.cff)).

Two conventions hold throughout:

- **Direction.** Every effect is `numerator − denominator`, taken from the
  explicit contrast definition. No stage infers direction from group names or
  sample order.
- **Multiplicity.** Significance uses Benjamini–Hochberg adjusted *p*-values
  (DESeq2's default `padj`), thresholded by the study's configured
  `figures.de.fdr` and `figures.de.abs_log2fc`.

Randomization seeds, provider snapshots, and environment/container identity are
recorded in the release manifest (see [Reproducibility](#reproducibility-and-provenance)).

## 1. Input quantification and annotation

Aligned, coordinate-sorted BAMs are checked with samtools and quantified with
`featureCounts` (subread) using fully declared paired-end and strandedness
settings; when strandedness is not fixed, it is inferred by comparing forward
and reverse assignment rates across the configured test modes. GTF annotations
supply stable gene identifiers, symbols, coordinates, strand, and biotype
(GTF-first, so no identifier is guessed). Each counting boundary also emits
per-gene exon lengths and length-normalized TPM and FPKM matrices alongside the
integer counts. Validated count-matrix inputs carry no exon lengths and
therefore no abundance.

## 2. Quality control

QC operates on DESeq2 variance-stabilized (VST) expression and reports library
size, detected genes, zero-count and mitochondrial fractions, expression
densities, PCA, Pearson correlations, Euclidean sample distances, and
variable-gene heatmaps. PCA group envelopes use the within-group covariance when
at least three samples are available; for nearly collinear three-point groups
the envelope uses a documented minimum 0.20 minor-to-major axis ratio and
expands to contain every observed point. When a batch covariate is declared,
`limma::removeBatchEffect` produces batch-corrected PCA and sample-distance views
over the VST expression, with raw and corrected states shown side by side and
the per-batch variance explained reported.

## 3. Differential expression

DESeq2 fits arbitrary design formulas and resolves explicit contrasts.
Effect-size shrinkage is selectable — `apeglm`, `ashr`, `normal`, or `none` —
and FDR-significant genes are split into a five-class scheme
(`significant_up`, `significant_down`, `padj_only`, `lfc_only`,
`not_significant`) by the configured FDR and |log2FC| thresholds.

Three contrast types are first-class:

- **pairwise** — a named two-level `factor_numerator_vs_denominator` coefficient
  (the historical behavior, reproduced exactly);
- **coefficient** — an interaction/difference-in-differences coefficient from a
  factorial design (relevel + refit + shrink), which displays all design groups
  colored by `figures.group` rather than the two-level factor;
- **omnibus** — a DESeq2 likelihood-ratio test of a multi-level factor.

Any pairwise contrast can be independently confirmed with an edgeR
quasi-likelihood fit whose fold changes are reported next to the DESeq2 result.
When several pairwise contrasts are declared, a cross-contrast **consensus** step
tabulates each gene's per-contrast significance and signed direction, reports the
genes that agree in direction across a configurable number of contrasts, and
renders an UpSet-style intersection plot with a signed-direction heatmap. An
optional **variancePartition** step decomposes each gene's VST expression over
the design covariates (categorical terms as random effects, numeric as fixed)
across the most-variable genes and summarizes the median variance fraction per
covariate, emitting well-formed empty outputs when no covariate is usable.

Coefficient and omnibus contrasts are deliberately scoped to DE and pathway
stages only; the downstream two-group modules (composition, regulators,
networks, hypotheses, publication) expand over pairwise contrasts, which they
assume.

## 4. Functional enrichment and ontology

Enrichment combines hypergeometric over-representation against the tested gene
universe with `fgsea` (`fgseaMultilevel`) preranked enrichment whose leading-edge
genes are exported per pathway, ssGSEA/GSVA scored with `limma` contrasts, and
multi-track enrichment curves. GO defaults to the biological-process domain and
optionally adds cellular-component and molecular-function breadth; Reactome
pathways (the MSigDB CP:REACTOME subcollection) can be enabled alongside KEGG
(resolved via `KEGGREST`). The bespoke GO-BP figures are unchanged when the extra
domains are off — additional domains appear only in an additive faceted view.

When KEGG is available, an optional signaling-pathway-impact analysis (`SPIA`
over `graphite` KEGG reaction graphs, with Ensembl→Entrez mapping through the
species `org.*.eg.db`) combines over-representation with a topology-based
perturbation signal to rank pathways by directional impact; it degrades to a
documented skip when the topology database or a required package is unavailable.
An optional enrichment map clusters significant terms into an unweighted-Jaccard
gene-overlap network with deterministic greedy-modularity communities. Focused
curves may be drawn from configured gene panels and are labeled as configured
programs, not database pathways.

## 5. Cell-state composition

Curated signatures are scored and reported as **relative** state/composition
evidence, not cell fractions. When an external signature matrix is supplied — or
one of the curated marker-panel presets shipped with the package is named —
non-negative least-squares deconvolution (`nnls` / `scipy.optimize`)
additionally estimates per-sample cell fractions and reports each sample's
reconstruction correlation. These absolute fractions are kept distinct from the
relative signature scores, and the binary presets are documented as a relative
screen, not a calibrated quantitative reference.

## 6. Regulatory inference

Regulator outputs separate **unsigned** target programs (including imported GTRD
binding snapshots) from **signed** DoRothEA/VIPER inference. The primary signed
engine is `viper::viper(method = "scale")` over DoRothEA regulons weighted by
confidence tier (A/B/C/D/E → 1.0/0.75/0.5/0.25/0.1; Garcia-Alonso et al. 2019);
`decoupleR::run_viper` is the fallback when the `viper` package is unavailable.
GRN outputs retain the complete regulon edge set while separately recording the
displayed selection.

## 7. Networks

Protein–protein association edges are retrieved from the STRING REST API
(`string-db.org/api`) for the mapped input genes; the stage exports mapped and
unmapped inputs, the complete association edge set, detected communities, and
centrality. Association edges (STRING), regulatory edges (DoRothEA/VIPER), and
co-expression edges (WGCNA) are kept explicitly distinct throughout.

## 8. Advanced and exploratory modules

Optional SVA, WGCNA, Ollivier-Ricci co-expression curvature, mediation, power,
and multilayer integration are exploratory and flagged as such when sample sizes
are small. Curvature runs on the WGCNA co-expression graph under the classic
graph hop-distance ground metric and reports per-edge and per-gene curvature, a
per-module robustness scalar, and the negative-curvature inter-module bridge
genes that bottleneck the network. Every small-sample caveat is written as a
machine-readable warning and surfaced in the report and the release manifest.

## Reproducibility and provenance

Every run records random seeds, provider snapshots, and the resolved
configuration. The release `manifest.json` records normalized configuration,
inputs/results checksums, contrast semantics, provider receipts, warnings,
environment/container information, and repository revision. Environments are
pinned as byte-exact per-rule conda lockfiles (shipped in the wheel), and
`tifzoret verify` compares a run against a golden reference within declared
absolute/relative tolerances. See [migration.md](migration.md) for the
verification workflow and [architecture.md](architecture.md) for how the stages
compose.
