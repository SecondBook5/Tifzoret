# Changelog

## 0.1.1

- Code style: `ruff check` is clean on `src/`. One-line `a; b` statements were split onto separate lines (verified by comparing the parsed syntax tree before and after), unused imports and variables were removed, and two lambda assignments became functions; behaviour is unchanged.
- `analysis.settings.de.shrinkage` and the new `analysis.settings.de.filter` now reach the families desugared from `contrasts.tsv`; since the `de.R` retirement nothing read `shrinkage`, and `normal` was accepted by the schema but rejected by the compiler.
- An estimand that is already a single coefficient of the canonical fit is shrunk in place. Re-deriving it through a releveled refit failed the identity assertion on real data (genes with zero counts in a group have no unique MLE), which silently replaced apeglm with ashr for the whole table.
- DE summaries and `de_results.tsv` record FDR-significant genes whose displayed (shrunken) effect lies outside their unshrunk 95% interval (`shrinkage_outside_ci`), the signature of an apeglm optimiser failure on zero-count genes, and surface the count as a caveat.
- The study-wide QC variance-stabilising transform sees only genes with at least 10 reads in total (31,845 of 56,748 genes in a GRCm39 matrix were constant or near-empty), through `blind_vst`, which uses DESeq2's `vst()` when it has 1,000 informative genes and `varianceStabilizingTransformation()` otherwise. PCA variance, sample distance, composition scores and pathway scores change; library metrics do not.
- A curated GSVA panel (`figures.pathways.gsva_sets`) scores and tests every measured set of the collection it was curated from, so the adjusted P values are over the collection and the listed sets only choose what is drawn. A listed set that is not measured is an error.
- The GO over-representation plot derives its colour and size scales from the displayed terms; fixed limits silently removed terms with -log10(FDR) below 4 or more than 80 genes.
- A static source panel addressed with `..` out of the results root now resolves lexically, so a fresh clone no longer fails at DAG build because the results directory does not exist yet.
- Documentation no longer describes `gene_filter_threshold`, which never existed, and states the real change of filter default.

## 0.1.0 - unreleased mouse v1

- Pathways and ontology now run for every estimand, not only pairwise rows: the GSVA differential reuses the family's fitted design and the DE job's contrast vector (`estimand_score_contrast`), so interactions, main effects and weighted composites get pathway and ontology results.
- `condition_palette` falls back to a colour-blind-safe palette (with a warning) for levels the configured palette does not name, instead of returning `NA` colours.
- `validate` now fails early, in about a second, when BAMs and the GTF use different chromosome names (a `chr` prefix mismatch), instead of after hours of counting.
- The raw featureCounts result is cached per run (keyed on BAMs, GTF, counting settings, tool version), so editing non-counting config no longer recounts the BAMs; `TIFZORET_RECOUNT=1` forces a recount.
- Estimand DE figures no longer show a dangling "versus" or "−", and the volcano legend no longer overflows the panel.
- GSVA heatmaps show GO term names instead of GO ids; ontology figures draw every term and no longer overlap their labels; `ontology_domain_displayed.tsv` no longer splits rows on embedded line breaks.
- The report shows the per-analysis figures when no publication recipe is configured.
- Claims and figure recipes may reference estimands declared in `analysis.families` (families are compiled before claim validation).

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

