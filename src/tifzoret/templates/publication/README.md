# Publication profile scaffold

This template demonstrates the **publication** analysis profile, which includes advanced composition, regulator, and publication figure modules.

## Quick start

```bash
tifzoret validate project.yaml
tifzoret run project.yaml --cores 4
tifzoret figures build project.yaml --cores 4
```

## What's included

**Core inputs:**
- `counts.tsv`, `samples.tsv`, `annotation.tsv` — synthetic gene expression data
- `contrasts.tsv` — single treated vs. control comparison
- `gene_sets.gmt` — custom pathway definitions

**Publication requirements:**
- `cell_state_signatures.yaml` — cell-state composition signatures
- `regulon_edges.tsv` — transcription-factor regulatory network
- `hypotheses.yaml` — scientific claims
- `hypothesis_panels.yaml` — biological programs and evidence panels
- `figure_recipe.yaml` — publication figure layout

## Modules enabled

This scaffold runs **offline** with the following analysis modules:

- ✓ **qc** — quality control and PCA
- ✓ **de** — differential expression (DESeq2)
- ✓ **pathways** — pathway enrichment (ORA, GSEA, GSVA)
- ✓ **composition** — cell-state composition (ssGSEA)
- ✓ **regulators** — transcription factor activity (VIPER)
- ✓ **publication** — hypothesis testing and figure assembly
- ✗ **networks** — disabled (requires network connectivity)

The **networks** module is disabled because it requires:
1. A live STRING database fetch (network connectivity)
2. `species.taxonomy_id` set to a real organism (e.g., 10090 for mouse)
3. `resources.providers.string: true` in the configuration

To enable networks, set a real taxonomy ID in `project.yaml`:

```yaml
species:
  taxonomy_id: 10090  # mouse
resources:
  providers:
    string: true
analysis:
  modules:
    networks: true  # or remove the override
```

## File formats

### cell_state_signatures.yaml

Defines cell-type or cell-state marker gene sets for composition analysis:

```yaml
signatures:
  - id: signature_name
    label: Display label
    category: cell_state
    description: Optional description
    genes:
      - GeneA
      - GeneB
```

Each signature requires at least 2 genes. Gene symbols must match those in `annotation.tsv`.

### regulon_edges.tsv

Defines transcription factor regulatory relationships:

```
source	target	mor
TF1	GeneA	1
TF1	GeneB	-1
TF2	GeneC	1
```

- `source` — transcription factor (gene symbol)
- `target` — regulated gene (gene symbol)
- `mor` — mode of regulation: `1` for activation, `-1` for repression

All symbols must exist in `annotation.tsv`.

## Next steps

1. Replace synthetic data with your real counts/samples/contrasts
2. Update `gene_sets.gmt` with relevant pathways
3. Define cell-state signatures for your tissue/system
4. Provide a regulatory network (custom edges or enable DoRothEA)
5. Write hypotheses and define biological programs
6. Design your figure layout in `figure_recipe.yaml`
7. Run the full pipeline

See `tifzoret --help` and the documentation for details.
