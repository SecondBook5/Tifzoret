#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(patchwork)
})

args <- parse_cli(c("project-config", "counts", "samples", "annotation", "contrasts", "contrast-id", "outdir"))
cfg <- read_project(args[["project-config"]])
cfg$.counts <- normalizePath(args$counts, mustWork = TRUE)
cfg$.samples <- normalizePath(args$samples, mustWork = TRUE)
cfg$.annotation <- normalizePath(args$annotation, mustWork = TRUE)
cfg$.contrasts <- normalizePath(args$contrasts, mustWork = TRUE)
dirs <- ensure_output_dirs(args$outdir)

counts <- read_counts_contract(cfg$.counts)
metadata <- readr::read_tsv(cfg$.samples, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
annotation <- read_annotation_contract(cfg$.annotation)
contrasts <- readr::read_tsv(cfg$.contrasts, show_col_types = FALSE, progress = FALSE)
contrast <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
if (nrow(contrast) != 1L) stop("Could not resolve exactly one contrast: ", args[["contrast-id"]], call. = FALSE)

resolved <- resolve_contrast(contrast, cfg$design$formula)
factor_name <- resolved$factor_name
numerator <- resolved$numerator
denominator <- resolved$denominator
design_formula <- resolved$design_formula
metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id
for (field in all.vars(design_formula)) metadata[[field]] <- factor(metadata[[field]])
for (relevel_factor in names(resolved$reference_levels)) {
  metadata[[relevel_factor]] <- stats::relevel(factor(metadata[[relevel_factor]]), ref = resolved$reference_levels[[relevel_factor]])
}

keep <- rowSums(counts) >= 10L
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE], colData = metadata, design = design_formula)
dispersion_fit <- "parametric"
dds <- tryCatch(
  DESeq2::DESeq(dds, fitType = "parametric", quiet = TRUE),
  error = function(error) {
    if (!grepl("all gene-wise dispersion estimates are within", conditionMessage(error), fixed = TRUE)) {
      stop(error)
    }
    message("Parametric dispersion trend unavailable; using DESeq2 gene-wise dispersion estimates.")
    dispersion_fit <<- "gene-wise"
    fallback <- DESeq2::estimateSizeFactors(dds)
    fallback <- DESeq2::estimateDispersionsGeneEst(fallback, quiet = TRUE)
    DESeq2::dispersions(fallback) <- S4Vectors::mcols(fallback)$dispGeneEst
    DESeq2::nbinomWaldTest(fallback, quiet = TRUE)
  }
)
coefficient <- resolved$coefficient_name
if (!coefficient %in% DESeq2::resultsNames(dds)) {
  stop(
    "Could not resolve coefficient ", coefficient,
    "; available: ", paste(DESeq2::resultsNames(dds), collapse = ", "),
    call. = FALSE
  )
}
raw <- DESeq2::results(dds, name = coefficient, alpha = cfg$figures$de$fdr)
shrunk <- DESeq2::lfcShrink(dds, coef = coefficient, type = "apeglm")
saveRDS(dds, file.path(dirs$objects, "deseq2.rds"))

result_table <- data.frame(
  gene_id = rownames(raw),
  base_mean = raw$baseMean,
  log2_fold_change_raw = raw$log2FoldChange,
  log2_fold_change = shrunk$log2FoldChange,
  lfc_se = shrunk$lfcSE,
  statistic = raw$stat,
  p_value = raw$pvalue,
  adjusted_p_value = raw$padj,
  stringsAsFactors = FALSE
) %>%
  left_join(annotation, by = "gene_id") %>%
  mutate(
    gene_symbol = ifelse(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol),
    safe_p_value = ifelse(
      is.na(p_value),
      NA_real_,
      pmax(p_value, .Machine$double.xmin)
    ),
    negative_log10_p = -log10(safe_p_value),
    direction = case_when(
      !is.na(adjusted_p_value) & adjusted_p_value < cfg$figures$de$fdr & log2_fold_change >= cfg$figures$de$abs_log2fc ~ "up_in_numerator",
      !is.na(adjusted_p_value) & adjusted_p_value < cfg$figures$de$fdr & log2_fold_change <= -cfg$figures$de$abs_log2fc ~ "down_in_numerator",
      TRUE ~ "not_significant"
    ),
    # Paper's 5-class scheme for volcano/MA colouring. Kept SEPARATE from the
    # `direction` column above (consumed byte-identically by networks.py);
    # thresholds are the same config-driven cfg$figures$de$fdr / abs_log2fc.
    significance_class = classify_significance(
      adjusted_p_value, log2_fold_change, cfg$figures$de$fdr, cfg$figures$de$abs_log2fc
    ),
    contrast_id = args[["contrast-id"]],
    numerator = numerator,
    denominator = denominator
  ) %>%
  arrange(adjusted_p_value, desc(abs(log2_fold_change)))
readr::write_tsv(result_table, file.path(dirs$tables, "de_results.tsv"), na = "NA")

# Volcano gene labels. Default (`pooled`) = top-N significant genes by
# (padj, |log2FC|). With `figures.de.label_balance: directional` the top_labels
# budget is split into up- and down-regulated halves, each ranked by -log10(p)
# (raw p, tie-broken by |log2FC|). This mirrors the primary analysis, which
# labels the top 10 up + top 10 down significant genes by -log10(p): set
# `top_labels: 20` to reproduce it. Backfill keeps the total at top_labels when
# one direction has fewer significant genes than its half of the budget.
label_balance <- cfg$figures$de$label_balance
if (is.null(label_balance)) label_balance <- "pooled"
label_pool <- result_table %>%
  filter(direction != "not_significant")
if (identical(label_balance, "directional")) {
  n_total <- cfg$figures$de$top_labels
  n_up <- ceiling(n_total / 2)
  n_down <- n_total - n_up
  up_lab <- label_pool %>% filter(log2_fold_change > 0) %>% arrange(desc(negative_log10_p), desc(abs(log2_fold_change))) %>% slice_head(n = n_up)
  down_lab <- label_pool %>% filter(log2_fold_change < 0) %>% arrange(desc(negative_log10_p), desc(abs(log2_fold_change))) %>% slice_head(n = n_down)
  label_table <- bind_rows(up_lab, down_lab)
  if (nrow(label_table) < n_total) {
    backfill <- label_pool %>%
      filter(!gene_id %in% label_table$gene_id) %>%
      arrange(desc(negative_log10_p), desc(abs(log2_fold_change))) %>%
      slice_head(n = n_total - nrow(label_table))
    label_table <- bind_rows(label_table, backfill)
  }
} else {
  label_table <- label_pool %>%
    arrange(adjusted_p_value, desc(abs(log2_fold_change))) %>%
    slice_head(n = cfg$figures$de$top_labels)
}
volcano_table <- result_table %>%
  mutate(label = ifelse(gene_id %in% label_table$gene_id, gene_symbol, NA_character_))
readr::write_tsv(volcano_table, file.path(dirs$tables, "volcano_displayed.tsv"), na = "NA")

volcano_plot <- ggplot(volcano_table, aes(log2_fold_change, negative_log10_p, colour = significance_class)) +
  annotate("rect", xmin = cfg$figures$de$abs_log2fc, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F4A6A6", alpha = 0.10) +
  annotate("rect", xmin = -Inf, xmax = -cfg$figures$de$abs_log2fc, ymin = -Inf, ymax = Inf, fill = "#A6CEE3", alpha = 0.10) +
  geom_vline(xintercept = c(-cfg$figures$de$abs_log2fc, cfg$figures$de$abs_log2fc), colour = "#9AA4AC", linetype = 3, linewidth = 0.35) +
  geom_point(size = 1.15, alpha = 0.78) +
  ggrepel::geom_label_repel(
    data = dplyr::filter(volcano_table, !is.na(label)),
    aes(label = label), size = 2.35, label.size = 0.15, fill = scales::alpha("white", 0.90),
    box.padding = 0.26, point.padding = 0.18, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_colour_manual(
    values = SIGNIFICANCE_PALETTE,
    breaks = SIGNIFICANCE_CLASSES,
    drop = FALSE,
    labels = c(
      significant_up = paste0("Upregulated in ", numerator),
      significant_down = paste0("Downregulated in ", numerator),
      padj_only = sprintf("FDR < %.2g only", cfg$figures$de$fdr),
      lfc_only = sprintf("|log2FC| ≥ %.2g only", cfg$figures$de$abs_log2fc),
      ns = "Not significant"
    )
  ) +
  labs(
    title = paste0(numerator, " versus ", denominator),
    subtitle = sprintf("DESeq2 with apeglm shrinkage; FDR < %.2g and |log2FC| ≥ %.2g", cfg$figures$de$fdr, cfg$figures$de$abs_log2fc),
    x = paste0("Shrunken log2 fold-change (", numerator, " − ", denominator, ")"),
    y = expression(-log[10](italic(p))),
    colour = NULL
  ) +
  theme_publication(9.2) +
  theme(legend.position = "top")
save_plot_pair(volcano_plot, file.path(dirs$figures, "volcano"), 6.4, 5.4)

ma_table <- result_table %>%
  select(gene_id, gene_symbol, base_mean, log2_fold_change, adjusted_p_value, direction, significance_class, contrast_id)
readr::write_tsv(ma_table, file.path(dirs$tables, "ma_displayed.tsv"), na = "NA")
ma_plot <- ggplot(ma_table, aes(base_mean, log2_fold_change, colour = significance_class)) +
  geom_hline(yintercept = 0, colour = "#8B979F", linewidth = 0.35) +
  geom_point(size = 1.05, alpha = 0.72) +
  scale_x_log10(labels = scales::label_number()) +
  scale_colour_manual(values = SIGNIFICANCE_PALETTE, drop = FALSE, guide = "none") +
  labs(title = "MA plot", subtitle = "Shrunken effects versus mean normalized abundance", x = "Mean normalized count", y = "Shrunken log2 fold-change") +
  theme_publication(8.8)
save_plot_pair(ma_plot, file.path(dirs$figures, "ma"), 6.2, 4.9)

histogram_table <- function(values, bins = 40L) {
  values <- values[is.finite(values)]
  estimate <- graphics::hist(values, breaks = bins, plot = FALSE)
  data.frame(bin_left = head(estimate$breaks, -1), bin_right = tail(estimate$breaks, -1), bin_midpoint = estimate$mids, count = estimate$counts)
}
pvalue_histogram <- histogram_table(result_table$p_value, 40L)
lfc_histogram <- histogram_table(result_table$log2_fold_change, 50L)
readr::write_tsv(pvalue_histogram, file.path(dirs$tables, "pvalue_distribution_displayed.tsv"))
readr::write_tsv(lfc_histogram, file.path(dirs$tables, "lfc_distribution_displayed.tsv"))
pvalue_plot <- ggplot(pvalue_histogram, aes(bin_midpoint, count)) +
  geom_col(width = stats::median(pvalue_histogram$bin_right - pvalue_histogram$bin_left), fill = "#6C92AE", colour = "white", linewidth = 0.15) +
  labs(title = "P-value distribution", x = "Raw p-value", y = "Genes") + theme_publication(8.8)
lfc_plot <- ggplot(lfc_histogram, aes(bin_midpoint, count, fill = bin_midpoint >= 0)) +
  geom_col(width = stats::median(lfc_histogram$bin_right - lfc_histogram$bin_left), colour = "white", linewidth = 0.12) +
  scale_fill_manual(values = c(`FALSE` = "#7EAFCB", `TRUE` = "#D98282"), guide = "none") +
  labs(title = "Effect-size distribution", x = "Shrunken log2 fold-change", y = "Genes") + theme_publication(8.8)
save_plot_pair(pvalue_plot, file.path(dirs$figures, "pvalue_distribution"), 5.8, 4.5)
save_plot_pair(lfc_plot, file.path(dirs$figures, "lfc_distribution"), 5.8, 4.5)

selected_genes <- result_table %>%
  filter(!is.na(adjusted_p_value), is.finite(log2_fold_change)) %>%
  arrange(adjusted_p_value, desc(abs(log2_fold_change))) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  slice_head(n = cfg$figures$de$top_heatmap_genes)
expression_transform <- if (dispersion_fit == "parametric") "variance_stabilizing" else "log2_normalized"
vst_contrast <- if (dispersion_fit == "parametric") {
  DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
} else {
  DESeq2::normTransform(dds)
}
expression <- SummarizedExperiment::assay(vst_contrast)
# Pairwise contrasts display only their two factor levels (unchanged). A
# coefficient (interaction) contrast spans every group in the design, so it
# displays all samples and colours by the palette's canonical grouping
# (figures.group) rather than the label factor, which has no palette entries.
if (identical(resolved$type, "pairwise")) {
  display_samples <- metadata[[factor_name]] %in% c(denominator, numerator)
  display_group_col <- factor_name
  display_palette <- condition_palette(cfg, c(denominator, numerator))
  display_subtitle <- paste(numerator, "and", denominator, "samples")
} else {
  display_samples <- rep(TRUE, nrow(metadata))
  display_group_col <- cfg$figures$group
  display_palette <- condition_palette(cfg, levels(factor(metadata[[display_group_col]])))
  display_subtitle <- "All design groups (interaction contrast)"
}
de_pca <- stats::prcomp(t(expression[, display_samples, drop = FALSE]), center = TRUE, scale. = FALSE)
de_pca_variance <- 100 * de_pca$sdev^2 / sum(de_pca$sdev^2)
de_pca_table <- as.data.frame(de_pca$x[, 1:2, drop = FALSE]) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(metadata %>% tibble::rownames_to_column("metadata_row") %>% select(-metadata_row), by = "sample_id")
de_pca_ellipses <- ellipse_coordinates(de_pca_table, display_group_col, cfg$figures$pca$ellipse_level)
readr::write_tsv(de_pca_table, file.path(dirs$tables, "de_pca_coordinates.tsv"))
readr::write_tsv(de_pca_ellipses, file.path(dirs$tables, "de_pca_ellipses.tsv"))
de_pca_plot <- ggplot(de_pca_table, aes(PC1, PC2, colour = .data[[display_group_col]])) +
  {if (nrow(de_pca_ellipses)) geom_path(data = de_pca_ellipses, aes(PC1, PC2, colour = ellipse_group, group = ellipse_group), inherit.aes = FALSE, linewidth = 0.8)} +
  geom_point(size = 3.1) +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 2.5, show.legend = FALSE, max.overlaps = Inf) +
  scale_colour_manual(values = display_palette, drop = FALSE) +
  labs(title = "Contrast PCA", subtitle = display_subtitle, x = sprintf("PC1 (%.1f%%)", de_pca_variance[[1]]), y = sprintf("PC2 (%.1f%%)", de_pca_variance[[2]]), colour = NULL) +
  theme_publication(8.8) + theme(legend.position = "top")
save_plot_pair(de_pca_plot, file.path(dirs$figures, "de_pca"), 6.1, 5.0)
expression <- expression[selected_genes$gene_id, display_samples, drop = FALSE]
rownames(expression) <- selected_genes$gene_symbol[match(rownames(expression), selected_genes$gene_id)]
z <- row_zscore(expression, cfg$figures$de$z_limit)
row_order <- rownames(z)[stats::hclust(stats::dist(z), method = "complete")$order]
column_distance <- stats::as.dist(1 - stats::cor(z, method = "pearson"))
column_order <- colnames(z)[stats::hclust(column_distance, method = "average")$order]
heatmap <- tile_heatmap(z, row_order, column_order, legend_title = "Row z-score", base_size = 7.7)
heatmap_table <- heatmap$table %>%
  mutate(
    condition = as.character(metadata[as.character(sample_id), factor_name]),
    contrast_id = args[["contrast-id"]]
  )
readr::write_tsv(heatmap_table, file.path(dirs$tables, "de_heatmap_displayed.tsv"))

annotation_plot <- data.frame(
  sample_id = factor(column_order, levels = column_order),
  condition = metadata[column_order, display_group_col]
) %>%
  ggplot(aes(sample_id, 1, fill = condition)) +
  geom_tile() +
  scale_fill_manual(values = display_palette, drop = FALSE) +
  theme_void() +
  theme(legend.position = "top", plot.margin = margin(0, 55, 0, 35))
heatmap_subtitle <- if (identical(resolved$type, "pairwise")) {
  paste0("Displayed samples are limited to ", denominator, " and ", numerator, "; row z-scores clipped at ±", cfg$figures$de$z_limit)
} else {
  paste0("All design groups shown; row z-scores clipped at ±", cfg$figures$de$z_limit)
}
heatmap$plot <- heatmap$plot +
  labs(
    title = "Top DE genes with global hierarchical clustering",
    subtitle = heatmap_subtitle
  )
combined_heatmap <- annotation_plot / heatmap$plot + patchwork::plot_layout(heights = c(0.07, 1))
save_plot_pair(combined_heatmap, file.path(dirs$figures, "de_heatmap"), 7.3, max(6.0, 0.18 * nrow(z) + 2.2))

de_overview <- (volcano_plot | ma_plot) / (pvalue_plot | lfc_plot) +
  patchwork::plot_annotation(title = paste0("Differential-expression overview: ", numerator, " versus ", denominator))
save_plot_pair(de_overview, file.path(dirs$figures, "de_overview"), 12.6, 9.6)

write_json_file(
  list(
    project_id = cfg$project$id,
    contrast_id = args[["contrast-id"]],
    factor = factor_name,
    numerator = numerator,
    denominator = denominator,
    design = resolved$design_text,
    contrast_type = resolved$type,
    coefficient = coefficient,
    dispersion_fit = dispersion_fit,
    expression_transform = expression_transform,
    samples = ncol(dds),
    genes_tested = nrow(result_table),
    significant_up = sum(result_table$direction == "up_in_numerator"),
    significant_down = sum(result_table$direction == "down_in_numerator")
  ),
  file.path(args$outdir, "de_summary.json")
)
