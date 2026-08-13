# Configuration v2

`project.yaml` is validated before Snakemake builds a DAG. Relative paths are resolved from the configuration directory and unset environment variables fail with their variable names.

The required sections are `project`, `species`, `reference`, `inputs`, `analysis`, `resources`, `figures`, and `output`. BAM, nf-core, and archive boundaries also require `counting`.

`samples.tsv` is keyed by unique `sample_id` and can contain arbitrary design covariates. An optional comma-delimited `analysis_set` column allows one sheet to define multiple selected cohorts. `contrasts.tsv` requires `contrast_id`, `factor`, `numerator`, and `denominator`; the factor must occur in the design formula and both levels must occur in the selected samples.

Publication profiles require claims, panels, signatures, and a figure recipe. Their schemas are installed with the package and cross-file references are checked during `validate`.

Legacy development configurations can be converted without enabling new modules:

```bash
bulk-rna migrate-config old.yaml --output project.yaml \
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
strand test modes, and minimum inference dominance.

## Analysis

`analysis.design` is an R formula and can include arbitrary covariates.
`analysis.contrasts` points to the explicit contrast table. `profile` selects
`standard`, `publication`, or `full`; `analysis.modules` can override each
resolved module. `random_seed` is inherited by deterministic selections and
layouts.

Optional settings are namespaced by module:

- composition: minimum matched genes;
- regulators: confidence classes, minimum targets, and display count;
- networks: STRING score, display node cap, and layout seed;
- SVA/WGCNA/mediation: method parameters and recommended-sample warning
  thresholds.

Dependency errors are reported at validation time. Enabling a module never
silently enables a missing resource or changes contrast direction.

## Resources

`resources.gene_sets` always declares a custom GMT and size bounds; optional
collection identifiers select live MSigDB content. Provider switches control
MSigDB, GO, KEGG, STRING, DoRothEA, and snapshot-backed GTRD. Publication
composition also requires `cell_state_signatures`; custom regulator analysis
may use `regulon_edges`.

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

Run `bulk-rna figures init project.yaml` to scaffold all three publication
files. Constructor-based recipes are validated against the enabled modules and
declared contrasts. Run `bulk-rna figures catalog` to inspect the available
constructors and `bulk-rna figures gallery project.yaml` to compare built
variants. See [hypothesis-driven publication figures](figures.md) for the full
contract.

The recipe selects already declared workflow artifacts. It cannot execute
arbitrary code or make an undeclared analysis direction.

## Output

`output.root` resolves relative to `project.yaml`. A run writes beneath
`<root>/<project.id>/<analysis_set>/`, allowing primary, full, and sensitivity
analyses to coexist without overwriting one another.
