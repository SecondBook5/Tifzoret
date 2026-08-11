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

factor_name <- contrast$factor[[1]]
numerator <- contrast$numerator[[1]]
denominator <- contrast$denominator[[1]]
metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id
design_formula <- stats::as.formula(cfg$design$formula)
for (field in all.vars(design_formula)) metadata[[field]] <- factor(metadata[[field]])
metadata[[factor_name]] <- stats::relevel(factor(metadata[[factor_name]]), ref = denominator)

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
coefficient_pattern <- paste0("^", factor_name, "_", numerator, "_vs_", denominator, "$")
coefficient <- grep(coefficient_pattern, DESeq2::resultsNames(dds), value = TRUE)
if (length(coefficient) != 1L) {
  stop(
    "Could not resolve coefficient ", coefficient_pattern,
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
    contrast_id = args[["contrast-id"]],
    numerator = numerator,
    denominator = denominator
  ) %>%
  arrange(adjusted_p_value, desc(abs(log2_fold_change)))
readr::write_tsv(result_table, file.path(dirs$tables, "de_results.tsv"), na = "NA")

label_table <- result_table %>%
  filter(direction != "not_significant") %>%
  arrange(adjusted_p_value, desc(abs(log2_fold_change))) %>%
  slice_head(n = cfg$figures$de$top_labels)
volcano_table <- result_table %>%
  mutate(label = ifelse(gene_id %in% label_table$gene_id, gene_symbol, NA_character_))
readr::write_tsv(volcano_table, file.path(dirs$tables, "volcano_displayed.tsv"), na = "NA")

direction_colors <- c(
  down_in_numerator = "#39799C",
  not_significant = "#B8C1C8",
  up_in_numerator = "#B55252"
)
volcano_plot <- ggplot(volcano_table, aes(log2_fold_change, negative_log10_p, colour = direction)) +
  annotate("rect", xmin = cfg$figures$de$abs_log2fc, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F4A6A6", alpha = 0.10) +
  annotate("rect", xmin = -Inf, xmax = -cfg$figures$de$abs_log2fc, ymin = -Inf, ymax = Inf, fill = "#A6CEE3", alpha = 0.10) +
  geom_vline(xintercept = c(-cfg$figures$de$abs_log2fc, cfg$figures$de$abs_log2fc), colour = "#9AA4AC", linetype = 3, linewidth = 0.35) +
  geom_hline(yintercept = -log10(cfg$figures$de$fdr), colour = "#9AA4AC", linetype = 3, linewidth = 0.35) +
  geom_point(size = 1.15, alpha = 0.78) +
  ggrepel::geom_label_repel(
    data = dplyr::filter(volcano_table, !is.na(label)),
    aes(label = label), size = 2.35, label.size = 0.15, fill = scales::alpha("white", 0.90),
    box.padding = 0.26, point.padding = 0.18, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_colour_manual(values = direction_colors, breaks = names(direction_colors), labels = c(
    down_in_numerator = paste0("Downregulated in ", numerator),
    not_significant = "Not significant",
    up_in_numerator = paste0("Upregulated in ", numerator)
  )) +
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
display_samples <- metadata[[factor_name]] %in% c(denominator, numerator)
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

palette <- condition_palette(cfg, c(denominator, numerator))
annotation_plot <- data.frame(
  sample_id = factor(column_order, levels = column_order),
  condition = metadata[column_order, factor_name]
) %>%
  ggplot(aes(sample_id, 1, fill = condition)) +
  geom_tile() +
  scale_fill_manual(values = palette, drop = FALSE) +
  theme_void() +
  theme(legend.position = "top", plot.margin = margin(0, 55, 0, 35))
heatmap$plot <- heatmap$plot +
  labs(
    title = "Top DE genes with global hierarchical clustering",
    subtitle = paste0("Displayed samples are limited to ", denominator, " and ", numerator, "; row z-scores clipped at ±", cfg$figures$de$z_limit)
  )
combined_heatmap <- annotation_plot / heatmap$plot + patchwork::plot_layout(heights = c(0.07, 1))
save_plot_pair(combined_heatmap, file.path(dirs$figures, "de_heatmap"), 7.3, max(6.0, 0.18 * nrow(z) + 2.2))

write_json_file(
  list(
    project_id = cfg$project$id,
    contrast_id = args[["contrast-id"]],
    factor = factor_name,
    numerator = numerator,
    denominator = denominator,
    design = cfg$design$formula,
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
