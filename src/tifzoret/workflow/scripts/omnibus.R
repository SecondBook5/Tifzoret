#!/usr/bin/env Rscript

# Omnibus (analysis-of-deviance) differential expression for a factor with >2
# levels. Unlike the signed pairwise/coefficient DE in de.R, an omnibus contrast
# asks a single question -- "does this gene differ across ANY level of the
# factor?" -- via a DESeq2 likelihood-ratio test of the full design against a
# reduced design that drops the tested factor. There is no numerator/denominator
# and no shrunken log2 fold-change, so this stage emits an LRT statistic and
# per-gene p-values rather than a volcano/MA. The reduced formula is supplied by
# the contrast's `reduced` column and is validated (must drop the factor) in
# config.load_project and parsed by resolve_contrast in utils.R.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
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
if (!identical(resolved$type, "omnibus")) {
  stop("contrast ", args[["contrast-id"]], " is not an omnibus contrast (type = ", resolved$type, ")", call. = FALSE)
}
factor_name <- resolved$factor_name
design_formula <- resolved$design_formula
reduced_formula <- resolved$reduced_formula

# Every design/reduced variable must be a factor with fixed reference levels so
# the LRT is well defined and reproducible (mirrors de.R's coercion + relevel).
model_vars <- union(all.vars(design_formula), all.vars(reduced_formula))
metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id
for (field in model_vars) metadata[[field]] <- factor(metadata[[field]])
for (relevel_factor in names(resolved$reference_levels)) {
  metadata[[relevel_factor]] <- stats::relevel(factor(metadata[[relevel_factor]]), ref = resolved$reference_levels[[relevel_factor]])
}
factor_levels <- levels(metadata[[factor_name]])

keep <- rowSums(counts) >= 10L
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE], colData = metadata, design = design_formula)
# LRT of the full design against the reduced design (which drops the tested
# factor). Fall back to gene-wise dispersions when the parametric trend fails,
# exactly as de.R does, so small-n multi-level designs stay robust.
dispersion_fit <- "parametric"
dds <- tryCatch(
  DESeq2::DESeq(dds, test = "LRT", reduced = reduced_formula, fitType = "parametric", quiet = TRUE),
  error = function(error) {
    if (!grepl("all gene-wise dispersion estimates are within", conditionMessage(error), fixed = TRUE)) {
      stop(error)
    }
    message("Parametric dispersion trend unavailable; using DESeq2 gene-wise dispersion estimates.")
    dispersion_fit <<- "gene-wise"
    fallback <- DESeq2::estimateSizeFactors(dds)
    fallback <- DESeq2::estimateDispersionsGeneEst(fallback, quiet = TRUE)
    DESeq2::dispersions(fallback) <- S4Vectors::mcols(fallback)$dispGeneEst
    DESeq2::nbinomLRT(fallback, reduced = reduced_formula, quiet = TRUE)
  }
)
saveRDS(dds, file.path(dirs$objects, "deseq2_lrt.rds"))

fdr <- cfg$figures$de$fdr
lrt <- DESeq2::results(dds, alpha = fdr)
result_table <- data.frame(
  gene_id = rownames(lrt),
  base_mean = lrt$baseMean,
  lrt_statistic = lrt$stat,
  p_value = lrt$pvalue,
  adjusted_p_value = lrt$padj,
  stringsAsFactors = FALSE
) %>%
  left_join(annotation, by = "gene_id") %>%
  mutate(
    gene_symbol = ifelse(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol),
    safe_p_value = ifelse(is.na(p_value), NA_real_, pmax(p_value, .Machine$double.xmin)),
    negative_log10_p = -log10(safe_p_value),
    significant = !is.na(adjusted_p_value) & adjusted_p_value < fdr,
    contrast_id = args[["contrast-id"]],
    factor = factor_name
  ) %>%
  arrange(adjusted_p_value, desc(lrt_statistic))
readr::write_tsv(result_table, file.path(dirs$tables, "omnibus_results.tsv"), na = "NA")

# P-value distribution -- the primary omnibus diagnostic (a spike near zero over
# a uniform background is the expected shape when the factor drives real signal).
histogram_table <- function(values, bins = 40L) {
  values <- values[is.finite(values)]
  estimate <- graphics::hist(values, breaks = bins, plot = FALSE)
  data.frame(bin_left = head(estimate$breaks, -1), bin_right = tail(estimate$breaks, -1), bin_midpoint = estimate$mids, count = estimate$counts)
}
pvalue_histogram <- histogram_table(result_table$p_value, 40L)
readr::write_tsv(pvalue_histogram, file.path(dirs$tables, "pvalue_distribution_displayed.tsv"))
pvalue_plot <- ggplot(pvalue_histogram, aes(bin_midpoint, count)) +
  geom_col(width = stats::median(pvalue_histogram$bin_right - pvalue_histogram$bin_left), fill = "#6C92AE", colour = "white", linewidth = 0.15) +
  labs(title = "Omnibus p-value distribution", subtitle = sprintf("LRT across %d levels of %s", length(factor_levels), factor_name), x = "Raw p-value", y = "Genes") +
  theme_publication(8.8)
save_plot_pair(pvalue_plot, file.path(dirs$figures, "pvalue_distribution"), 5.8, 4.5)

# Top-gene heatmap: VST expression of the most-significant genes across every
# sample, columns grouped by the tested factor's levels so the multi-level
# pattern that the LRT detects is visible. VST is built from this stage's own
# dds (independent of the qc module), matching de.R's fallback-aware transform.
expression_transform <- if (dispersion_fit == "parametric") "variance_stabilizing" else "log2_normalized"
vst_object <- if (dispersion_fit == "parametric") {
  DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
} else {
  DESeq2::normTransform(dds)
}
expression_all <- SummarizedExperiment::assay(vst_object)
top_genes <- result_table %>%
  filter(!is.na(adjusted_p_value)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  slice_head(n = cfg$figures$de$top_heatmap_genes)
if (nrow(top_genes) >= 2L) {
  column_order <- rownames(metadata)[order(metadata[[factor_name]], rownames(metadata))]
  expression <- expression_all[top_genes$gene_id, column_order, drop = FALSE]
  rownames(expression) <- top_genes$gene_symbol[match(rownames(expression), top_genes$gene_id)]
  z <- row_zscore(expression, cfg$figures$de$z_limit)
  row_order <- rownames(z)[stats::hclust(stats::dist(z), method = "complete")$order]
  heatmap <- tile_heatmap(z, row_order, column_order, legend_title = "Row z-score", base_size = 7.7)
  heatmap_table <- heatmap$table %>%
    mutate(
      condition = as.character(metadata[as.character(sample_id), factor_name]),
      contrast_id = args[["contrast-id"]]
    )
  readr::write_tsv(heatmap_table, file.path(dirs$tables, "omnibus_heatmap_displayed.tsv"))
  display_palette <- condition_palette(cfg, factor_levels)
  annotation_plot <- data.frame(
    sample_id = factor(column_order, levels = column_order),
    condition = factor(metadata[column_order, factor_name], levels = factor_levels)
  ) %>%
    ggplot(aes(sample_id, 1, fill = condition)) +
    geom_tile() +
    scale_fill_manual(values = display_palette, drop = FALSE) +
    theme_void() +
    theme(legend.position = "top", plot.margin = margin(0, 55, 0, 35))
  heatmap$plot <- heatmap$plot +
    labs(
      title = "Top omnibus genes across factor levels",
      subtitle = paste0("Samples grouped by ", factor_name, "; row z-scores clipped at ±", cfg$figures$de$z_limit)
    )
  combined_heatmap <- annotation_plot / heatmap$plot + patchwork::plot_layout(heights = c(0.07, 1))
} else {
  # Fewer than two significant genes -> no informative heatmap; keep the output
  # contract (both files present) with an explicit placeholder + empty table.
  readr::write_tsv(
    data.frame(feature = character(0), sample_id = character(0), value = numeric(0), condition = character(0), contrast_id = character(0)),
    file.path(dirs$tables, "omnibus_heatmap_displayed.tsv")
  )
  combined_heatmap <- empty_plot("Top omnibus genes across factor levels", "No genes passed the display threshold")
}
save_plot_pair(combined_heatmap, file.path(dirs$figures, "omnibus_heatmap"), 7.3, max(6.0, 0.18 * nrow(top_genes) + 2.2))

write_json_file(
  list(
    project_id = cfg$project$id,
    contrast_id = args[["contrast-id"]],
    factor = factor_name,
    levels = factor_levels,
    contrast_type = resolved$type,
    design = resolved$design_text,
    reduced = resolved$reduced_text,
    test = "LRT",
    dispersion_fit = dispersion_fit,
    expression_transform = expression_transform,
    samples = ncol(dds),
    genes_tested = nrow(result_table),
    significant = sum(result_table$significant)
  ),
  file.path(args$outdir, "omnibus_summary.json")
)
