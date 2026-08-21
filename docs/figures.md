# Hypothesis-driven publication figures

Tifzoret separates statistical discovery from presentation choices. Analysis
modules produce complete tables and standard figures. A publication project then
declares biological programs, expected effects, registered panel constructors,
variants, and assembly geometry without modifying analysis code.

```text
analysis results
      ↓
hypotheses.yaml + hypothesis_panels.yaml
      ↓
registered constructors and variants
      ↓
figure_recipe.yaml
      ↓
staged panels, review gallery, assembled PDF/PNG
```

Hypotheses guide which evidence is displayed; they never alter statistical
results, significance thresholds, or contrast direction.

## Start a publication configuration

From an existing version 2 project:

```bash
tifzoret figures init project.yaml
```

This creates `hypotheses.yaml`, `hypothesis_panels.yaml`, and
`figure_recipe.yaml`, adds their references to `project.yaml`, and enables the
publication constructor module. It does not enable composition, regulators, or
networks automatically; constructors backed by those modules become available
when the corresponding analysis modules and resources are enabled.

The generated biological program and pathway identifiers are deliberate
placeholders. Replace them before building. Use `--force` only when
intentionally replacing all three publication files.

## Declare programs and expectations

Prefer the first-class `programs` contract for new projects:

```yaml
programs:
  contractile:
    label: Contractile / cytoskeletal
    description: Smooth-muscle contractile machinery
    color: "#D97706"
    genes: [ACTA2, CNN1, MYL9, TPM1, TPM2, TNNT2]
    expected_direction: increased
    contrast: treated_vs_control

  epithelial:
    label: Epithelial identity / junctions
    color: "#E9A300"
    genes: [EPCAM, CLDN7, KRT5, KRT18, GATA3]

pathway_panels:
  remodeling:
    description: Pathways central to the remodeling hypothesis
    pathways:
      - collection: hallmark
        pathway: HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
        expected_direction: increased

expected_effects:
  contractile_in_treatment:
    contrast: treated_vs_control
    target_type: program
    target: contractile
    direction: increased
    rationale: Prespecified mechanistic expectation.

gsea_programs: [contractile, epithelial]
program_order: [Contractile / cytoskeletal, Epithelial identity / junctions]
```

Legacy `gene_panels` remain valid. Internally, both forms resolve to the same
auditable grouped-gene contract.

`tables/program_definitions.tsv` records every requested gene, its resolved
measured symbol, mapping status, program color, expected direction, and expected
contrast. Unmapped genes are retained in this audit table rather than silently
disappearing.

## Select registered constructors

List the installed catalog:

```bash
tifzoret figures catalog
tifzoret figures catalog --json
```

Recipes should use a constructor instead of a raw result path:

```yaml
figure_sets:
  primary:
    title: Primary publication figure
    width: 12
    height: 10
    units: in
    columns: 2
    shared_legends: true
    panels:
      - id: A
        constructor: pca_correlation
        row: 1
        column: 1
        column_span: 2

      - id: B
        constructor: de_heatmap
        contrast: treated_vs_control
        variant: program_grouped
        row: 2
        column: 1

      - id: C
        constructor: gsea_multitrack
        contrast: treated_vs_control
        row: 2
        column: 2
        options:
          scale: 0.90
```

Raw `source` paths remain supported for custom external panels. Such panels
should declare `displayed_data` explicitly; otherwise their staged metadata
contains a visible audit warning.

## Constructor catalog

| Constructor | Variants |
|---|---|
| `pca`, `sample_correlation`, `pca_correlation` | `default` |
| `library_metrics`, `expression_density`, `sample_distance` | `default` |
| `variable_gene_heatmap`, `sample_outliers`, `qc_overview` | `default` |
| `volcano`, `ma`, `de_overview` | `default` |
| `de_heatmap` | `default`, `global_clustered`, `program_grouped`, `direct_program_labels` |
| `cell_state_effects` | `default` |
| `ora_bidirectional`, `go_ora` | `default` |
| `gsva_heatmap`, `gsea_multitrack` | `default` |
| `program_heatmap_effects`, `program_violins` | `default` |
| `string_enrichment` | `faceted`, `combined` |
| `string_network` | `upregulated`, `downregulated` |
| `regulator_activity` | `default` |
| `dorothea_grn` | `rectangular`, `radial` |
| `wgcna_module_trait`, `multilayer_network` | `default` |

This table lists the primary variant of each constructor. A few constructors
also carry a `*_legacy` matplotlib-rendered fallback (for example
`string_network`'s `upregulated_legacy`); the complete, always-current
enumeration is `tifzoret figures catalog`, and every registered variant is
rendered in the review gallery. Omitting `variant` selects the constructor's
default (shown first above where a constructor has several).

Configuration validation rejects unknown constructors/variants, missing
contrast IDs, and constructors whose analysis module is disabled.

## Build and review

```bash
tifzoret figures build project.yaml --cores 8
tifzoret figures gallery project.yaml
```

Build one configured figure set with:

```bash
tifzoret figures build project.yaml --figure-set primary --cores 8
```

The gallery contains `index.html`, `contact_sheet.png`, `gallery.json`, and
review copies of every available registered PNG beneath
`publication/gallery/`. It shows all built constructor variants supported by
the enabled modules and marks those selected in `figure_recipe.yaml`.

## Auditable panel output

Each selected panel is staged beneath:

```text
publication/<figure_set>/panels/<panel_id>/
├── panel.pdf
├── panel.png
├── panel.json
└── displayed_data/
    └── <copied displayed tables and selection metadata>
```

`panel.json` records the constructor, variant, contrast, caption, merged
constructor/panel options, source paths, SHA-256 hashes, displayed-data hashes,
and warnings. `panels/index.json` inventories the complete figure set.

The assembled output remains:

```text
publication/<figure_set>/assembled/
├── <figure_set>.pdf
├── <figure_set>.png
└── assembly.json
```

Constructor defaults may be shared in `hypothesis_panels.yaml` and overridden
per panel:

```yaml
constructor_defaults:
  de_heatmap:
    scale: 0.94
    show_panel_label: true
```

Currently implemented assembly options are `scale` (0.1–1.0),
`show_panel_label`, and `panel_label`. Biological selections and heatmap/GSEA
variants belong in the program/pathway definitions and `variant` field rather
than being hidden in layout options.
