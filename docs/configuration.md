# Configuration v2

`project.yaml` is validated before Snakemake builds a DAG. Relative paths are resolved from the configuration directory and unset environment variables fail with their variable names.

The required sections are `project`, `species`, `reference`, `inputs`, `analysis`, `resources`, `figures`, and `output`. BAM, nf-core, and archive boundaries also require `counting`.

`samples.tsv` is keyed by unique `sample_id` and can contain arbitrary design covariates. An optional comma-delimited `analysis_set` column allows one sheet to define multiple selected cohorts. `contrasts.tsv` requires `contrast_id`, `factor`, `numerator`, and `denominator`; the factor must occur in the design formula and both levels must occur in the selected samples.

Publication profiles require claims, panels, signatures, and a figure recipe. Their schemas are installed with the package and cross-file references are checked during `validate`.

Legacy development configurations can be converted without enabling new modules:

```bash
tifzoret migrate-config old.yaml --output project.yaml \
  --species mouse --genome-build GRCm39 --annotation-release 107
```

## Project and species

- `project.id`: filesystem-safe stable identifier; `project.title` and
  `project.description` are display metadata.
- `species.provider`: `mouse`, `human`, or `custom`; pair mouse/human with the
  correct scientific name and taxonomy identifier.
- `reference.genome_build`, `annotation_release`, and optional
  `expected_contigs` describe the reference used upstream.

## Inputs

Every input declares `samples` and may declare `analysis_set`.

- `counts`: add `counts` and `annotation`.
- `bam`: add `bam_root`, `gtf`, and a `bam` column in `samples.tsv`.
- `nfcore_rnaseq`: add `root`, `gtf`, and `bam_pattern`; placeholders are
  resolved from metadata columns.
- `archive`: add `archive`, `gtf`, optional `member_root`, and a `bam` column
  containing safe archive-relative member paths.

BAM boundaries require `counting`: threads, feature type, GTF attribute,
paired-end flags, strandedness (`infer`, `unstranded`, `forward`, or `reverse`),
strand test modes, and minimum inference dominance. Every counting boundary also
writes per-gene exon lengths and length-normalized TPM/FPKM matrices next to the
integer counts; the `counts` input carries no exon lengths and so emits no
abundance.

## Analysis

`analysis.design` is an R formula and can include arbitrary covariates.
`analysis.contrasts` points to the explicit contrast table; each row's optional
`type` column selects `pairwise` (default), `coefficient` (a named design
coefficient), or `omnibus` (a DESeq2 likelihood-ratio test across all levels of
a multi-level factor, which carries no numerator/denominator). `profile` selects
`standard`, `publication`, or `full`; `analysis.modules` can override each
resolved module. Batch correction, edgeR confirmation, signature-matrix
deconvolution, Ollivier-Ricci curvature, cross-contrast consensus, SPIA
pathway-topology impact, variancePartition variance decomposition, and
enrichment-map term clustering belong to no profile and run only when switched
on under `analysis.modules`; batch correction additionally reads
`analysis.batch`, a `samples.tsv` column, curvature requires `wgcna`, consensus
requires at least two pairwise contrasts, and SPIA requires
`resources.providers.kegg`. `random_seed` is inherited by deterministic
selections and layouts.

Optional settings live under `analysis.settings.<module>` (for example
`analysis.settings.de.shrinkage` or `analysis.settings.deconvolution.min_genes`),
namespaced by module:

- composition: minimum matched genes;
- regulators: confidence classes, minimum targets, and display count;
- networks: STRING score, display node cap, and layout seed;
- de: `shrinkage` (`apeglm`, `ashr`, `normal`, or `none`) and, for the optional
  confirmatory fit, `confirm_method` (`edger`);
- deconvolution: `method` (`nnls`) and `min_genes`;
- curvature: lazy-random-walk `alpha` and the reported `top_bridges` count;
- consensus: `min_contrasts` (the count in which a gene must agree in direction)
  and displayed `top_genes`;
- spia: DE-subset `fdr` (defaults to `figures.de.fdr`) and displayed
  `top_pathways`;
- variance_partition: `covariates` (defaults to every fixed design term) and
  `top_variable_genes`;
- enrichment_map: edge `min_similarity` (Jaccard) and `top_terms` per direction;
- SVA/WGCNA/mediation: method parameters and recommended-sample warning
  thresholds. WGCNA `network_neighbors` sets the per-gene edge count of the
  co-expression graph exported for curvature.

Dependency errors are reported at validation time. Enabling a module never
silently enables a missing resource or changes contrast direction.

## Resources

`resources.gene_sets` always declares a custom GMT and size bounds; optional
collection identifiers select live MSigDB content. Provider switches control
MSigDB, GO, KEGG, Reactome, STRING, DoRothEA, and snapshot-backed GTRD.
`resources.go_domains` selects the GO domains the `go` provider resolves; it
defaults to `["BP"]` so existing GO-BP studies are byte-identical, and `CC`
and/or `MF` add cellular-component and molecular-function breadth. Reactome
pathways are drawn from the MSigDB `CP:REACTOME` subcollection, so the provider
needs no dependency beyond MSigDB. Publication composition also requires
`cell_state_signatures`; custom regulator analysis may use `regulon_edges`.
Signature-matrix deconvolution requires either `resources.deconvolution_signature`,
a gene-by-cell-type reference table (a gene column and at least two cell-type
columns), or `resources.deconvolution_preset` naming a reference matrix shipped
with the package; the two are mutually exclusive. Presets are curated binary
marker panels (see the deconvolution data directory's provenance note), so their
fractions are a relative screen rather than calibrated abundances.

Provider resources use explicit species metadata. Mouse is taxonomy 10090;
human is 9606. Cached resource receipts include provider, organism, release,
retrieval time, parameters, upstream license notice, and checksum.
`resources.offline: true` reuses cache entries and fails clearly when one is
absent; `refresh: true` forces retrieval.

GTRD is intentionally snapshot-backed: set `resources.providers.gtrd: true`
and point `resources.regulon_edges` to an exported `source`/`target` table
obtained under the applicable upstream terms. Its checksum is part of the cache
key and resource receipt. Because binding does not establish activation or
repression, edges without `mor` are analyzed as unsigned target-program
evidence.

## Figures and publication files

`figures.group` selects the displayed metadata field and `palette` must cover
every selected level. PCA ellipse confidence, QC selection, DE thresholds and
z-score limit, and pathway display counts are explicit configuration.

Publication/full projects declare:

- `hypotheses.claims`: statements, expected directions, and evidence panels;
- `hypotheses.panels`: first-class biological `programs` (preferred),
  legacy-compatible grouped `gene_panels`, selected pathways, and expected
  effects. Optional `program_annotations` assigns a mutually exclusive
  display program to each top-DE heatmap gene, while `program_colors` and
  `program_order` control its study-owned legend and grouped layout.
  Optional `gsea_programs` is an ordered list of gene-panel identifiers; each
  referenced panel is tested as a user-configured preranked gene program and
  receives a focused multi-track GSEA curve. Database pathways remain the
  automatic fallback when this list is absent;
- `publication.recipe`: registered constructors or custom panel sources,
  variants, contrast IDs, letters, grid positions, dimensions, and legend
  policy.

Each of `hypotheses.claims`, `hypotheses.panels`, and `publication.recipe` may
be given either as a path to a sibling file or as an inline mapping written
directly in `project.yaml`. Inlining all three yields a single self-contained
study file — the form an authoring UI produces — and the engine validates and
runs it identically to the multi-file form.

Run `tifzoret figures init project.yaml` to scaffold all three publication
files. Constructor-based recipes are validated against the enabled modules and
declared contrasts. Run `tifzoret figures catalog` to inspect the available
constructors and `tifzoret figures gallery project.yaml` to compare built
variants. See [hypothesis-driven publication figures](figures.md) for the full
contract.

The recipe selects already declared workflow artifacts. It cannot execute
arbitrary code or make an undeclared analysis direction.

## Output

`output.root` resolves relative to `project.yaml`. A run writes beneath
`<root>/<project.id>/<analysis_set>/`, allowing primary, full, and sensitivity
analyses to coexist without overwriting one another.
