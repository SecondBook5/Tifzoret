# Changelog

## 0.1.0 - unreleased mouse v1

- Configuration v2, profiles, four input adapters, canonical contracts, and migration command.
- QC, DE, enrichment, ontology, composition, regulators, STRING/GRN, hypotheses, publication, and systems modules.
- Report, vector assembly, collection meta-analysis, verification, provider caching, and provenance manifests.
- DE shrinkage is selectable (`apeglm`/`ashr`/`normal`/`none`); likelihood-ratio `omnibus` contrasts test multi-level factors; an optional edgeR quasi-likelihood fit confirms pairwise contrasts.
- Length-normalized TPM/FPKM abundance and per-gene exon lengths are emitted on every counting boundary.
- Opt-in modules (no profile, `analysis.modules` only): `batch` (limma batch-corrected PCA/distance), `deconvolution` (signature-matrix NNLS cell fractions), `curvature` (Ollivier-Ricci curvature of the WGCNA co-expression graph), `consensus` (cross-contrast significance/sign concordance with an UpSet intersection plot and signed heatmap), `spia` (SPIA/graphite pathway-topology impact over KEGG, with a graceful skip when the topology database is unavailable), `variance_partition` (per-covariate expression-variance decomposition), and `enrichment_map` (Jaccard gene-overlap network of enriched terms with greedy-modularity communities).
- Enrichment breadth: the `go` provider resolves configurable GO domains (`BP` by default; opt-in `CC`/`MF`) and a `reactome` provider adds MSigDB CP:REACTOME pathways; the bespoke GO-BP outputs are unchanged when the extra domains stay off. fgsea leading-edge genes are exported per pathway.
- `deconvolution` accepts `resources.deconvolution_preset`, a curated marker-panel reference matrix shipped with the package (documented as a relative screen), as an alternative to a supplied signature.
- The single-env `environment.yaml` (`--no-conda` path) is now a complete superset of the per-rule envs — it gained `scipy` and every R package the optional modules need (`viper`, `ComplexHeatmap`, `circlize`, `dendextend`, `igraph`, `ggforce`, `colorspace`, `ashr`, `edgeR`, `nnls`, `variancePartition`, `matrixStats`, `SPIA`, `graphite`) so a full-feature run under `--no-conda` no longer fails on a missing dependency; a test guards the superset invariant.
- Release hardening: the pinned conda lockfiles (`workflow/envs/locks/*.lock.txt`) now ship in the built wheel as well as the sdist, so a `pip install` of the package can reconstruct the exact solved environments. The project-agnostic guard was widened to catch every planned reference-study fingerprint across the whole packaged root (not just the workflow directory), and orphaned relocated bytecode was removed. Execution coverage was added for the pure-Python opt-in modules (`consensus`, `enrichment_map`, `curvature` — the last gated on `scipy`) and for the `coefficient`/`omnibus` contrast validators, DE `shrinkage` prior rejection, and bounded opt-in GO domains.

