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
display_group <- cfg$figures$group
configured_group_order <- names(unlist(cfg$figures$palette, use.names = TRUE))
observed_group_order <- configured_group_order[configured_group_order %in% as.character(metadata[[display_group]])]
metadata[[display_group]] <- factor(metadata[[display_group]], levels = observed_group_order)

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

gene_symbols <- annotation$gene_symbol[match(rownames(counts), annotation$gene_id)]
gene_symbols[is.na(gene_symbols)] <- rownames(counts)[is.na(gene_symbols)]
mitochondrial <- grepl("^mt-", gene_symbols, ignore.case = TRUE)
library_sizes <- colSums(counts)
library_metrics <- data.frame(
  sample_id = colnames(counts),
  library_size = as.numeric(library_sizes),
  detected_genes = as.numeric(colSums(counts > 0)),
  zero_fraction = as.numeric(colMeans(counts == 0)),
  mitochondrial_fraction = if (any(mitochondrial)) {
    as.numeric(colSums(counts[mitochondrial, , drop = FALSE]) / pmax(library_sizes, 1))
  } else {
    rep(0, ncol(counts))
  },
  stringsAsFactors = FALSE
) %>%
  left_join(metadata %>% tibble::rownames_to_column("metadata_row") %>% select(-metadata_row), by = "sample_id")
readr::write_tsv(library_metrics, file.path(dirs$tables, "library_metrics.tsv"))

log_cpm <- log2(t(t(counts) / pmax(library_sizes, 1) * 1e6) + 1)
density_table <- bind_rows(lapply(seq_len(ncol(log_cpm)), function(index) {
  estimate <- stats::density(log_cpm[, index], from = 0, n = 256)
  data.frame(sample_id = colnames(log_cpm)[[index]], log2_cpm = estimate$x, density = estimate$y)
})) %>%
  left_join(library_metrics %>% select(sample_id, all_of(cfg$figures$group)), by = "sample_id")
readr::write_tsv(density_table, file.path(dirs$tables, "expression_density_displayed.tsv"))

pca <- stats::prcomp(t(vst_matrix), center = TRUE, scale. = FALSE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_table <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(metadata %>% tibble::rownames_to_column("metadata_row") %>% select(-metadata_row), by = "sample_id")
variance_table <- data.frame(component = paste0("PC", seq_along(variance)), variance_percent = variance)
readr::write_tsv(pca_table, file.path(dirs$tables, "pca_coordinates.tsv"))
readr::write_tsv(variance_table, file.path(dirs$tables, "pca_variance.tsv"))

group_col <- cfg$figures$group
groups <- levels(droplevels(pca_table[[group_col]]))
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

sample_distance <- as.matrix(stats::dist(t(vst_matrix), method = "euclidean"))
distance_order <- colnames(sample_distance)[stats::hclust(stats::as.dist(sample_distance), method = "average")$order]
distance_table <- as.data.frame(sample_distance, check.names = FALSE) %>% tibble::rownames_to_column("sample_id")
readr::write_tsv(distance_table, file.path(dirs$tables, "sample_distance.tsv"))

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
combined_pca_correlation <- pca_plot | combined_correlation
save_plot_pair(combined_pca_correlation, file.path(dirs$figures, "pca_correlation"), 12.4, 5.8)
write_json_file(list(
  panels = list(
    list(
      id = "PCA",
      displayed_data = "tables/pca_coordinates.tsv",
      ellipses = unique(ellipse_table$ellipse_group),
      ellipse_method = "covariance ellipse with a 0.20 minimum minor-to-major axis ratio"
    ),
    list(id = "correlation", displayed_data = "tables/sample_correlation.tsv", clustering = "average linkage on 1 - Pearson correlation")
  ),
  shared_group_legend = TRUE,
  correlation_condition_legend = FALSE
), file.path(dirs$tables, "pca_correlation_layout.json"))

metrics_long <- library_metrics %>%
  select(sample_id, all_of(group_col), library_size, detected_genes, zero_fraction, mitochondrial_fraction) %>%
  pivot_longer(c(library_size, detected_genes, zero_fraction, mitochondrial_fraction), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("library_size", "detected_genes", "zero_fraction", "mitochondrial_fraction"), labels = c("Library size", "Detected genes", "Zero-count fraction", "Mitochondrial fraction")))
metrics_plot <- ggplot(metrics_long, aes(sample_id, value, fill = .data[[group_col]])) +
  geom_col(width = 0.72) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = palette, drop = FALSE) +
  scale_y_continuous(labels = scales::label_number()) +
  labs(title = "Library-level quality metrics", x = NULL, y = NULL, fill = NULL) +
  theme_publication(8.3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_plot_pair(metrics_plot, file.path(dirs$figures, "library_metrics"), 7.2, 5.6)

density_plot <- ggplot(density_table, aes(log2_cpm, density, colour = .data[[group_col]], group = sample_id)) +
  geom_line(linewidth = 0.65, alpha = 0.85) +
  scale_colour_manual(values = palette, drop = FALSE) +
  labs(title = "Expression distributions", subtitle = "Per-sample kernel densities of log2(CPM + 1)", x = "log2(CPM + 1)", y = "Density", colour = NULL) +
  theme_publication(8.8) + theme(legend.position = "top")
save_plot_pair(density_plot, file.path(dirs$figures, "expression_density"), 6.6, 4.8)

distance_long <- distance_table %>%
  pivot_longer(-sample_id, names_to = "sample_id_y", values_to = "euclidean_distance") %>%
  mutate(sample_id = factor(sample_id, levels = distance_order), sample_id_y = factor(sample_id_y, levels = rev(distance_order)))
distance_plot <- ggplot(distance_long, aes(sample_id, sample_id_y, fill = euclidean_distance)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "magma", name = "Euclidean\ndistance") +
  labs(title = "Sample distance", subtitle = "Euclidean distance on variance-stabilized expression", x = NULL, y = NULL) +
  theme_publication(8.1) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks = element_blank())
save_plot_pair(distance_plot, file.path(dirs$figures, "sample_distance"), 6.2, 5.5)

top_variable_count <- if (is.null(cfg$figures$qc$top_variable_genes)) 50L else as.integer(cfg$figures$qc$top_variable_genes)
variances <- apply(vst_matrix, 1, stats::var)
variable_ids <- names(sort(variances, decreasing = TRUE))[seq_len(min(top_variable_count, length(variances)))]
variable_matrix <- vst_matrix[variable_ids, , drop = FALSE]
rownames(variable_matrix) <- gene_symbols[match(variable_ids, rownames(counts))]
variable_matrix <- variable_matrix[!duplicated(rownames(variable_matrix)), , drop = FALSE]
variable_z <- row_zscore(variable_matrix, 2)
variable_rows <- rownames(variable_z)[stats::hclust(stats::dist(variable_z), method = "complete")$order]
variable_columns <- colnames(variable_z)[stats::hclust(stats::as.dist(1 - stats::cor(variable_z)), method = "average")$order]
variable_heatmap <- tile_heatmap(variable_z, variable_rows, variable_columns, legend_title = "Row z-score", base_size = 7.1)
readr::write_tsv(variable_heatmap$table, file.path(dirs$tables, "variable_gene_heatmap_displayed.tsv"))
variable_heatmap$plot <- variable_heatmap$plot + labs(title = "Most variable genes", subtitle = paste0("Top ", nrow(variable_z), " genes by variance-stabilized variance"))
save_plot_pair(variable_heatmap$plot, file.path(dirs$figures, "variable_gene_heatmap"), 7.2, max(6.2, 0.14 * nrow(variable_z) + 2))

qc_overview <- (pca_plot | metrics_plot) / (density_plot | distance_plot) +
  patchwork::plot_annotation(title = "Bulk RNA-seq quality-control overview")
save_plot_pair(qc_overview, file.path(dirs$figures, "qc_overview"), 13.2, 10.2)

write_json_file(
  list(
    project_id = cfg$project$id,
    design = cfg$design$formula,
    samples = ncol(counts),
    genes = nrow(counts),
    mitochondrial_genes = sum(mitochondrial),
    median_library_size = stats::median(library_sizes),
    median_detected_genes = stats::median(library_metrics$detected_genes),
    pca_variance_pc1 = variance[[1]],
    pca_variance_pc2 = variance[[2]],
    ellipse_groups_drawn = unique(ellipse_table$ellipse_group),
    ellipse_minimum_minor_major_ratio = 0.20
  ),
  file.path(args$outdir, "qc_summary.json")
)
