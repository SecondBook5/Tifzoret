# Deconvolution reference presets — provenance and honest caveats

These files are **binary marker panels**, not quantitative expression
signatures. Each column is a cell type; each row is a gene; a `1` means the gene
is a canonical, textbook marker of that cell type and a `0` means it is not.
They ship so a study can enable the `deconvolution` module without first
sourcing a signature matrix, by setting `resources.deconvolution_preset: <name>`
instead of `resources.deconvolution_signature: <path>`.

## What these presets ARE

- Curated lists of **well-established, widely-cited lineage markers** for common
  mouse cell types, encoded as a presence/absence (0/1) matrix keyed by mouse
  gene symbol.
- A convenient, dependency-free starting point for an **exploratory, relative,
  qualitative** read on how bulk composition shifts between conditions.

## What these presets ARE NOT — read before interpreting results

- **Not quantitative signatures.** A real deconvolution reference (e.g.
  CIBERSORTx LM22, or a signature derived from sorted-cell or single-cell RNA-seq
  on comparable tissue and platform) carries per-gene *expression magnitudes*
  that encode how strongly each gene distinguishes cell types. A binary panel
  discards all of that. NNLS against a 0/1 matrix therefore yields **rough,
  relative** fractions, not calibrated absolute abundances.
- **Not tissue- or platform-matched.** The markers are generic lineage markers,
  not tuned to any particular tissue, disease state, or protocol.
- **Not exhaustive.** Only a handful of markers per cell type are included, and
  cell types absent from the panel are silently unmodelled — their signal is
  redistributed across the included types.
- **Sensitive to marker overlap and co-expression.** Genes co-expressed across
  lineages, or absent from a given dataset, degrade identifiability.

**Recommendation:** treat preset-based fractions as a screen. For any claim that
reaches a figure or the text, supply a proper `deconvolution_signature` matched
to the tissue and platform, and report the per-sample reconstruction correlation
(the module already writes it) so a weak fit is visible.

## Files

- `mouse_immune.tsv` — 6 immune cell types
  (T cell, B cell, NK cell, macrophage, neutrophil, dendritic cell).
- `mouse_immune_stromal.tsv` — the 6 immune types above plus 3 structural/stromal
  types (fibroblast, endothelial, epithelial).

Marker sources are the standard immunology/cell-biology lineage markers found in
canonical references (e.g. CD3/CD19/NCR1 lineage genes, F4/80·Csf1r for
macrophages, Ly6g·S100a8/9 for neutrophils, Pecam1·Cdh5 for endothelium,
Col1a1·Pdgfra for fibroblasts, Epcam·Krt8/18 for epithelium). They are common
knowledge in the field rather than proprietary to any single dataset.

## Adding a preset

Drop a `<name>.tsv` here with a `gene_symbol` (or gene-id) first column followed
by one column per cell type. It is discovered automatically —
`resources.deconvolution_preset: <name>` will resolve to it, and the name will
appear in the validation error listing available presets. Keep this file honest:
document what the numbers mean and how the panel was built.
