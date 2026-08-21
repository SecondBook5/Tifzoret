#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Publication "clipped" volcano (opt-in via figures.de.volcano_style: clipped).
# ---------------------------------------------------------------------------
# Reproduces the primary-analysis volcano exactly: a display-clipped y-axis with
# off-scale genes drawn as open up-triangles at the top edge (so a handful of
# extreme hits do not stretch the axis and flatten the main cloud), a four-class
# "Gene class" legend (the non-significant cloud is drawn but omitted from the
# legend), and a tiered label selector that ranks by a salience score
# (significance x effect size) rather than adjusted p-value alone. An optional
# priority-gene list is labelled first so a study's mechanism genes are shown
# even in weak-signal contrasts. Display-only: no statistic, threshold, ranking,
# or table value is altered -- only which points are clamped and which are
# labelled. Ported from the validated old pipeline (workflow/stages/de/*.R).

# Symmetric display cap: max(min_cap, quantile(|values|, prob)). Display only.
volcano_display_cap <- function(values, min_cap, prob = 0.995) {
  finite <- values[is.finite(values)]
  if (!length(finite)) return(as.numeric(min_cap))
  cap <- as.numeric(stats::quantile(abs(finite), probs = prob, na.rm = TRUE))
  if (!is.finite(cap)) return(as.numeric(min_cap))
  max(as.numeric(min_cap), cap)
}

# Tiered volcano-label selector operating on engine result-table columns.
# Tier order: priority (study mechanism genes, only when FDR-significant) ->
# strict up -> strict down -> other FDR-significant -> strongest remaining. The
# up/down budget left after the priority tier is split half-and-half; unused
# slots spill forward. Ranks by a salience score and deduplicates on the visible
# symbol so repeated symbols never waste slots. Returns selected gene_ids.
select_clipped_volcano_labels <- function(tbl, top_n, fdr, lfc, priority_symbols = character(0)) {
  top_n <- as.integer(top_n)
  if (is.na(top_n) || top_n <= 0L) return(character(0))
  candidates <- tbl %>%
    mutate(
      plot_label = as.character(gene_symbol),
      has_symbol = !is.na(gene_symbol) & nzchar(trimws(as.character(gene_symbol))),
      abs_lfc = abs(log2_fold_change),
      neglogp = ifelse(!is.na(negative_log10_p) & is.finite(negative_log10_p), negative_log10_p, 0),
      basemean_score = ifelse(!is.na(base_mean) & is.finite(base_mean) & base_mean >= 0, log10(base_mean + 1), 0),
      label_score = (neglogp * pmax(abs_lfc, 0.25)) + (0.35 * abs_lfc) + (0.05 * basemean_score),
      is_strict_up = !is.na(adjusted_p_value) & is.finite(adjusted_p_value) & adjusted_p_value < fdr & is.finite(log2_fold_change) & log2_fold_change >= lfc,
      is_strict_down = !is.na(adjusted_p_value) & is.finite(adjusted_p_value) & adjusted_p_value < fdr & is.finite(log2_fold_change) & log2_fold_change <= -lfc,
      is_other_fdr = !is.na(adjusted_p_value) & is.finite(adjusted_p_value) & adjusted_p_value < fdr & !(is_strict_up | is_strict_down)
    ) %>%
    filter(
      !is.na(gene_id), nzchar(trimws(as.character(gene_id))),
      !is.na(plot_label), nzchar(trimws(plot_label)),
      is.finite(log2_fold_change), is.finite(label_score)
    )
  if (!nrow(candidates)) return(character(0))

  take_rows <- function(pool, n, excl_ids, excl_labels, tier) {
    if (is.na(n) || n <= 0L || !nrow(pool)) return(pool[0, , drop = FALSE])
    pool %>%
      filter(!(gene_id %in% excl_ids), !(plot_label %in% excl_labels)) %>%
      arrange(desc(has_symbol), desc(label_score), adjusted_p_value, desc(abs_lfc), desc(base_mean)) %>%
      distinct(plot_label, .keep_all = TRUE) %>%
      arrange(desc(label_score), adjusted_p_value, desc(abs_lfc), desc(base_mean)) %>%
      slice_head(n = as.integer(n)) %>%
      mutate(tier = tier)
  }

  priority_symbols <- toupper(trimws(as.character(priority_symbols)))
  priority_symbols <- priority_symbols[nzchar(priority_symbols)]
  sel_priority <- candidates[0, , drop = FALSE]
  if (length(priority_symbols)) {
    priority_pool <- candidates %>%
      filter(toupper(plot_label) %in% priority_symbols, is_strict_up | is_strict_down | is_other_fdr)
    if (nrow(priority_pool)) {
      sel_priority <- take_rows(priority_pool, min(nrow(priority_pool), top_n), character(0), character(0), "priority")
    }
  }
  budget <- max(0L, top_n - nrow(sel_priority))
  up_target <- ceiling(budget / 2)
  down_target <- budget - up_target
  sel_up <- take_rows(filter(candidates, is_strict_up), up_target, sel_priority$gene_id, sel_priority$plot_label, "strict_up")
  sel_down <- take_rows(
    filter(candidates, is_strict_down), down_target,
    c(sel_priority$gene_id, sel_up$gene_id), c(sel_priority$plot_label, sel_up$plot_label), "strict_down"
  )
  selected <- bind_rows(sel_priority, sel_up, sel_down)
  remaining <- top_n - nrow(selected)
  if (remaining > 0L) {
    sel_other <- take_rows(filter(candidates, is_other_fdr), remaining, selected$gene_id, selected$plot_label, "other_fdr")
    selected <- bind_rows(selected, sel_other)
    remaining <- top_n - nrow(selected)
  }
  if (remaining > 0L) {
    sel_overall <- take_rows(candidates, remaining, selected$gene_id, selected$plot_label, "overall")
    selected <- bind_rows(selected, sel_overall)
  }
  selected <- selected %>%
    mutate(tier = factor(tier, levels = c("priority", "strict_up", "strict_down", "other_fdr", "overall"))) %>%
    arrange(tier, desc(label_score), adjusted_p_value, desc(abs_lfc), desc(base_mean)) %>%
    slice_head(n = top_n)
  as.character(selected$gene_id)
}

# Render the clipped publication volcano. Returns the ggplot plus the displayed
# (x-filtered, y-clamped, label-annotated) table written as volcano_displayed.tsv.
build_clipped_volcano <- function(result_table, cfg, numerator, denominator, priority_symbols = character(0)) {
  fdr <- cfg$figures$de$fdr
  lfc <- cfg$figures$de$abs_log2fc
  point_size <- 1.4
  vt <- result_table %>% filter(is.finite(log2_fold_change), is.finite(negative_log10_p))
  x_cap <- volcano_display_cap(vt$log2_fold_change, min_cap = max(2.5, lfc * 2.5))
  y_cap <- volcano_display_cap(vt$negative_log10_p, min_cap = 6)
  n_fdr_sig <- sum(!is.na(vt$adjusted_p_value) & vt$adjusted_p_value < fdr)
  n_strict <- sum(vt$significance_class %in% c("significant_up", "significant_down"))
  # Clamp (never drop) the strongest genes to the top edge so they stay visible
  # and labelable; the x-cap filters the few genes beyond the fold-change window.
  visible <- vt %>%
    filter(abs(log2_fold_change) <= x_cap) %>%
    mutate(y_offscale = negative_log10_p > y_cap, neg_log10_p_display = pmin(negative_log10_p, y_cap))
  label_ids <- select_clipped_volcano_labels(visible, cfg$figures$de$top_labels, fdr, lfc, priority_symbols)
  label_rows <- visible %>% filter(gene_id %in% label_ids) %>% distinct(gene_symbol, .keep_all = TRUE)
  labeled_points <- visible %>% filter(gene_id %in% label_ids)
  displayed <- visible %>% mutate(label = ifelse(gene_id %in% label_ids, gene_symbol, NA_character_))
  subtitle <- paste0(
    numerator, " vs ", denominator,
    " · ", n_fdr_sig, " genes at FDR < ", formatC(fdr, format = "f", digits = 2),
    " · ", n_strict, " also |log2FC| ≥ ", formatC(lfc, format = "f", digits = 1)
  )
  plot <- ggplot(visible, aes(log2_fold_change, neg_log10_p_display)) +
    annotate("rect", xmin = -Inf, xmax = -lfc, ymin = 0, ymax = Inf, fill = "#E8EEF6", alpha = 0.28) +
    annotate("rect", xmin = lfc, xmax = Inf, ymin = 0, ymax = Inf, fill = "#F7E8E8", alpha = 0.28) +
    geom_point(
      data = dplyr::filter(visible, significance_class == "ns"),
      colour = SIGNIFICANCE_PALETTE[["ns"]], size = point_size, alpha = 0.32, na.rm = TRUE
    ) +
    geom_point(
      data = dplyr::filter(visible, significance_class != "ns", !y_offscale),
      aes(colour = significance_class), size = point_size + 0.1, alpha = 0.82, na.rm = TRUE
    ) +
    geom_point(
      data = dplyr::filter(visible, significance_class != "ns", y_offscale),
      aes(colour = significance_class), shape = 24, fill = "white",
      size = point_size + 0.6, stroke = 0.5, alpha = 0.95, na.rm = TRUE, show.legend = FALSE
    ) +
    geom_point(
      data = labeled_points, aes(colour = significance_class),
      size = point_size + 1.0, alpha = 0.96, stroke = 0.2, na.rm = TRUE, show.legend = FALSE
    ) +
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dotted", colour = "#6B7280", linewidth = 0.45) +
    {if (nrow(label_rows)) ggrepel::geom_label_repel(
      data = label_rows, aes(log2_fold_change, neg_log10_p_display, label = gene_symbol),
      size = 3.15, fontface = "bold", colour = "#111827", fill = scales::alpha("white", 0.96),
      box.padding = 0.5, point.padding = 0.25,
      label.padding = grid::unit(0.14, "lines"), label.r = grid::unit(0.10, "lines"), label.size = 0.18,
      segment.alpha = 0.78, segment.size = 0.32, segment.colour = "#4B5563",
      max.overlaps = Inf, force = 3.2, force_pull = 0.7, min.segment.length = 0,
      direction = "both", seed = 42, show.legend = FALSE, inherit.aes = FALSE
    )} +
    scale_colour_manual(
      values = SIGNIFICANCE_PALETTE,
      breaks = c("significant_up", "significant_down", "padj_only", "lfc_only"),
      labels = c(
        paste0("FDR<", fdr, ", up"),
        paste0("FDR<", fdr, ", down"),
        paste0("FDR<", fdr, " only"),
        paste0("|LFC|≥", lfc, " only")
      ),
      name = "Gene class"
    ) +
    labs(
      title = paste0(numerator, ": volcano plot"),
      subtitle = subtitle,
      x = paste0("Shrunken log2 fold change (", numerator, " vs ", denominator, ")"),
      y = expression(-log[10](raw ~ pvalue))
    ) +
    coord_cartesian(xlim = c(-x_cap, x_cap), ylim = c(0, y_cap)) +
    theme_publication(9.2) +
    theme(legend.position = "right")
  list(plot = plot, displayed = displayed)
}

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
# lfcShrink estimator for the reported log2 fold-change. apeglm is the validated
# default; ashr/normal are alternative shrinkage priors; none reports the raw MLE
# (shrunk <- raw) for callers that want unshrunken effects. v1 configs and studies
# that omit the knob keep apeglm, so this is fully backward-compatible.
shrinkage <- cfg$analysis$settings$de$shrinkage
if (is.null(shrinkage)) shrinkage <- "apeglm"
shrunk <- if (identical(shrinkage, "none")) raw else DESeq2::lfcShrink(dds, coef = coefficient, type = shrinkage)
shrinkage_label <- if (identical(shrinkage, "none")) "no shrinkage (raw MLE)" else paste0(shrinkage, " shrinkage")
effect_axis <- if (identical(shrinkage, "none")) "log2 fold-change" else "Shrunken log2 fold-change"
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

# Volcano gene labels + plot. Two styles, both opt-in and backward-compatible:
#
#   figures.de.volcano_style: clipped  -> the validated primary-analysis volcano
#     (display-clipped y-axis with off-scale up-triangles, four-class legend,
#     tiered salience label selector). Optional figures.de.priority_genes labels
#     a study's mechanism genes first. This reproduces the published Panel C.
#
#   (default / standard) -> the house-style volcano. Labels default to top-N by
#     (padj, |log2FC|); figures.de.label_balance: directional splits the
#     top_labels budget into up/down halves ranked by -log10(p).
volcano_style <- cfg$figures$de$volcano_style
if (is.null(volcano_style)) volcano_style <- "standard"
priority_symbols <- cfg$figures$de$priority_genes
if (is.null(priority_symbols)) priority_symbols <- character(0)
priority_symbols <- as.character(unlist(priority_symbols, use.names = FALSE))

if (identical(volcano_style, "clipped")) {
  clipped <- build_clipped_volcano(result_table, cfg, numerator, denominator, priority_symbols)
  volcano_plot <- clipped$plot
  readr::write_tsv(clipped$displayed, file.path(dirs$tables, "volcano_displayed.tsv"), na = "NA")
} else {
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
      subtitle = sprintf("DESeq2 with %s; FDR < %.2g and |log2FC| ≥ %.2g", shrinkage_label, cfg$figures$de$fdr, cfg$figures$de$abs_log2fc),
      x = paste0(effect_axis, " (", numerator, " − ", denominator, ")"),
      y = expression(-log[10](italic(p))),
      colour = NULL
    ) +
    theme_publication(9.2) +
    theme(legend.position = "top")
}
save_plot_pair(volcano_plot, file.path(dirs$figures, "volcano"), 6.4, 5.4)

ma_table <- result_table %>%
  select(gene_id, gene_symbol, base_mean, log2_fold_change, adjusted_p_value, direction, significance_class, contrast_id)
readr::write_tsv(ma_table, file.path(dirs$tables, "ma_displayed.tsv"), na = "NA")
ma_plot <- ggplot(ma_table, aes(base_mean, log2_fold_change, colour = significance_class)) +
  geom_hline(yintercept = 0, colour = "#8B979F", linewidth = 0.35) +
  geom_point(size = 1.05, alpha = 0.72) +
  scale_x_log10(labels = scales::label_number()) +
  scale_colour_manual(values = SIGNIFICANCE_PALETTE, drop = FALSE, guide = "none") +
  labs(title = "MA plot", subtitle = sprintf("%s versus mean normalized abundance", if (identical(shrinkage, "none")) "Maximum-likelihood effects" else "Shrunken effects"), x = "Mean normalized count", y = effect_axis) +
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
  labs(title = "Effect-size distribution", x = effect_axis, y = "Genes") + theme_publication(8.8)
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
    shrinkage = shrinkage,
    dispersion_fit = dispersion_fit,
    expression_transform = expression_transform,
    samples = ncol(dds),
    genes_tested = nrow(result_table),
    significant_up = sum(result_table$direction == "up_in_numerator"),
    significant_down = sum(result_table$direction == "down_in_numerator")
  ),
  file.path(args$outdir, "de_summary.json")
)
