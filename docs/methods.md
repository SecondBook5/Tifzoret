# Methods implemented

Aligned, coordinate-sorted BAM files are checked with samtools and quantified with featureCounts using fully declared paired-end and strandedness settings. Strand inference compares forward and reverse assignment rates. GTF annotations supply stable gene identifiers, symbols, coordinates, strand, and biotype.

QC uses DESeq2 variance stabilization and reports library size, detected genes, zero-count and mitochondrial fractions, expression densities, PCA, Pearson correlations, Euclidean distances, and variable-gene heatmaps. DESeq2 models arbitrary formulas and explicit contrasts; apeglm provides shrunken effect sizes.

Enrichment includes hypergeometric ORA with the tested gene universe, fgseaMultilevel preranked enrichment, ssGSEA/GSVA with limma contrasts, GO BP and KEGG views, and multi-track enrichment curves. Random seeds and provider snapshots are recorded.

Cell-state analysis scores curated signatures and is reported as relative state/composition evidence rather than cell fractions. Regulator outputs separate unsigned target programs, including imported GTRD binding snapshots, from signed DoRothEA/VIPER inference. STRING exports mapped/unmapped inputs, complete association edges, communities, and centrality. GRN outputs retain complete regulon edges while separately recording displayed selections.

Optional SVA, WGCNA, mediation, power, and multilayer results are exploratory when sample sizes are small. Their warnings are machine-readable and surfaced in the report and manifest.
