#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(patchwork)
})

args <- parse_cli(c("project-config", "counts", "samples", "annotation", "outdir"))
cfg <- read_project(args[["project-config"]])
cfg$.counts <- normalizePath(args$counts, mustWork = TRUE)
cfg$.samples <- normalizePath(args$samples, mustWork = TRUE)
cfg$.annotation <- normalizePath(args$annotation, mustWork = TRUE)
dirs <- ensure_output_dirs(args$outdir)

counts <- read_counts_contract(cfg$.counts)
metadata <- readr::read_tsv(cfg$.samples, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
annotation <- read_annotation_contract(cfg$.annotation)
metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id

design_formula <- stats::as.formula(cfg$design$formula)
for (field in all.vars(design_formula)) metadata[[field]] <- factor(metadata[[field]])
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = design_formula)
dds <- DESeq2::estimateSizeFactors(dds)
vst <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
saveRDS(vst, file.path(dirs$objects, "vst.rds"))

vst_matrix <- SummarizedExperiment::assay(vst)
symbol_matrix <- matrix_to_symbols(vst_matrix, annotation)
readr::write_tsv(
  as.data.frame(symbol_matrix, check.names = FALSE) %>% tibble::rownames_to_column("gene_symbol"),
  file.path(dirs$tables, "vst_expression.tsv")
)

pca <- stats::prcomp(t(vst_matrix), center = TRUE, scale. = FALSE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_table <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(metadata %>% tibble::rownames_to_column("metadata_row") %>% select(-metadata_row), by = "sample_id")
variance_table <- data.frame(component = paste0("PC", seq_along(variance)), variance_percent = variance)
readr::write_tsv(pca_table, file.path(dirs$tables, "pca_coordinates.tsv"))
readr::write_tsv(variance_table, file.path(dirs$tables, "pca_variance.tsv"))

group_col <- cfg$figures$group
groups <- unique(as.character(pca_table[[group_col]]))
palette <- condition_palette(cfg, groups)
ellipse_table <- ellipse_coordinates(pca_table, group_col, cfg$figures$pca$ellipse_level)

pca_plot <- ggplot(pca_table, aes(PC1, PC2, colour = .data[[group_col]])) +
  {if (nrow(ellipse_table)) geom_path(
    data = ellipse_table,
    aes(PC1, PC2, colour = ellipse_group, group = ellipse_group),
    inherit.aes = FALSE,
    linewidth = 0.85,
    alpha = 0.95
  )} +
  geom_point(size = 3.2, alpha = 0.95) +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 2.6, show.legend = FALSE, max.overlaps = Inf) +
  scale_colour_manual(values = palette, drop = FALSE) +
  labs(
    title = "Principal-component analysis",
    subtitle = "Ellipses summarize within-group covariance when at least three samples are available",
    x = sprintf("PC1 (%.1f%%)", variance[[1]]),
    y = sprintf("PC2 (%.1f%%)", variance[[2]]),
    colour = NULL
  ) +
  theme_publication(9.2) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
save_plot_pair(pca_plot, file.path(dirs$figures, "pca"), 6.2, 5.2)

correlation <- stats::cor(vst_matrix, method = "pearson")
correlation_order <- colnames(correlation)[stats::hclust(stats::as.dist(1 - correlation), method = "average")$order]
correlation_table <- as.data.frame(correlation, check.names = FALSE) %>% tibble::rownames_to_column("sample_id")
readr::write_tsv(correlation_table, file.path(dirs$tables, "sample_correlation.tsv"))

correlation_long <- correlation_table %>%
  pivot_longer(-sample_id, names_to = "sample_id_y", values_to = "pearson_r") %>%
  mutate(
    sample_id = factor(sample_id, levels = correlation_order),
    sample_id_y = factor(sample_id_y, levels = rev(correlation_order))
  )
annotation_plot <- metadata %>%
  tibble::rownames_to_column("sample_key") %>%
  mutate(sample_key = factor(sample_key, levels = correlation_order)) %>%
  ggplot(aes(sample_key, 1, fill = .data[[group_col]])) +
  geom_tile() +
  scale_fill_manual(values = palette, drop = FALSE) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 6, 0, 38))
correlation_plot <- ggplot(correlation_long, aes(sample_id, sample_id_y, fill = pearson_r)) +
  geom_tile(colour = "white", linewidth = 0.22) +
  scale_fill_gradientn(colours = c("#355C7D", "#A6CEE3", "#F7F4EE", "#F4A6A6", "#B55252"), limits = c(min(correlation), 1), name = "Pearson r") +
  labs(title = "Sample correlation", subtitle = "Average-linkage ordering of 1 − Pearson correlation", x = NULL, y = NULL) +
  theme_publication(8.2) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks = element_blank())
combined_correlation <- annotation_plot / correlation_plot + patchwork::plot_layout(heights = c(0.055, 1))
save_plot_pair(combined_correlation, file.path(dirs$figures, "sample_correlation"), 6.2, 5.7)

write_json_file(
  list(
    project_id = cfg$project$id,
    design = cfg$design$formula,
    samples = ncol(counts),
    genes = nrow(counts),
    pca_variance_pc1 = variance[[1]],
    pca_variance_pc2 = variance[[2]],
    ellipse_groups_drawn = unique(ellipse_table$ellipse_group)
  ),
  file.path(args$outdir, "qc_summary.json")
)
