#!/usr/bin/env Rscript

# variancePartition: how much of each gene's expression variance is attributable
# to each design covariate. QC (qc.R) shows how samples cluster; this module
# quantifies WHY -- for every gene it decomposes the variance-stabilized
# expression into the fraction explained by each covariate (condition, batch,
# sex, ...) and an unexplained Residual. Study-level and display/diagnostic only:
# it does not feed DE. Categorical covariates enter as random effects ((1|x), the
# variancePartition convention) and continuous covariates as fixed effects; the
# per-gene fractions sum to 1. A covariate whose median fraction is large is a
# major structured source of variance and a candidate for the DE design or the
# batch term.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(variancePartition)
})

args <- parse_cli(c("project-config", "vst", "samples", "outdir"))
cfg <- read_project(args[["project-config"]])
cfg$.samples <- normalizePath(args$samples, mustWork = TRUE)
dirs <- ensure_output_dirs(args$outdir)

settings <- cfg$analysis$settings$variance_partition
if (is.null(settings)) settings <- list()
top_variable_genes <- if (is.null(settings$top_variable_genes)) 2000L else as.integer(settings$top_variable_genes)

vst <- readRDS(normalizePath(args$vst, mustWork = TRUE))
vst_matrix <- SummarizedExperiment::assay(vst)
metadata <- readr::read_tsv(cfg$.samples, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
metadata <- metadata[match(colnames(vst_matrix), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id

# Covariates to decompose: the configured list if given, else the biological
# grouping plus the technical batch (when set) -- the two structured sources a
# study almost always cares about.
requested <- settings$covariates
if (is.null(requested) || length(requested) == 0L) {
  requested <- unique(c(cfg$figures$group, if (!is.null(cfg$analysis$batch) && nzchar(cfg$analysis$batch)) cfg$analysis$batch))
}
requested <- unique(as.character(requested))

warnings <- character(0)
usable <- character(0)
for (cov in requested) {
  if (!cov %in% colnames(metadata)) {
    warnings <- c(warnings, sprintf("covariate '%s' is not a samples.tsv column; skipped", cov))
    next
  }
  column <- metadata[[cov]]
  if (is.numeric(column)) {
    if (length(unique(column[!is.na(column)])) < 2L) {
      warnings <- c(warnings, sprintf("covariate '%s' is constant; skipped", cov))
      next
    }
  } else {
    metadata[[cov]] <- factor(column)
    if (nlevels(metadata[[cov]]) < 2L) {
      warnings <- c(warnings, sprintf("covariate '%s' has a single level; skipped", cov))
      next
    }
  }
  usable <- c(usable, cov)
}

fraction_fields <- c("gene_id", usable, "Residuals")
if (length(usable) == 0L) {
  # Nothing to decompose; emit empty-but-well-formed outputs and a placeholder
  # figure so the pipeline surfaces the reason rather than hard-failing.
  warnings <- c(warnings, "no usable covariates (need >=1 present, non-constant covariate); nothing decomposed")
  readr::write_tsv(
    stats::setNames(data.frame(matrix(character(0), nrow = 0, ncol = length(fraction_fields))), fraction_fields),
    file.path(dirs$tables, "variance_fractions.tsv")
  )
  readr::write_tsv(
    data.frame(variable = character(0), median_fraction = numeric(0), mean_fraction = numeric(0), median_percent = numeric(0)),
    file.path(dirs$tables, "variance_summary.tsv")
  )
  save_plot_pair(empty_plot("variancePartition", "No usable covariates to decompose"), file.path(dirs$figures, "variance_partition"), 6.4, 4.6)
  write_json_file(
    list(project_id = cfg$project$id, method = "variancePartition::fitExtractVarPartModel", covariates = list(), genes = 0L, warnings = warnings),
    file.path(args$outdir, "variance_partition_summary.json")
  )
  quit(save = "no", status = 0)
}

# Restrict to the most variable genes: variancePartition is per-gene and the
# tail of invariant genes is uninformative and slow to fit.
row_variance <- matrixStats::rowVars(vst_matrix)
n_keep <- min(top_variable_genes, nrow(vst_matrix))
kept_genes <- rownames(vst_matrix)[order(row_variance, decreasing = TRUE)][seq_len(n_keep)]
expression <- vst_matrix[kept_genes, , drop = FALSE]

# Categorical covariates as random effects ((1|x)); continuous as fixed effects.
term_for <- function(cov) if (is.numeric(metadata[[cov]])) cov else sprintf("(1|%s)", cov)
formula <- stats::as.formula(paste("~", paste(vapply(usable, term_for, character(1)), collapse = " + ")))

var_part <- suppressWarnings(variancePartition::fitExtractVarPartModel(expression, formula, metadata))
fractions <- as.data.frame(var_part, check.names = FALSE)
# variancePartition names random-effect columns by the bare covariate; keep the
# canonical column order (covariates then Residuals) for a stable table.
ordered_cols <- c(usable[usable %in% colnames(fractions)], "Residuals")
fractions <- fractions[, ordered_cols, drop = FALSE]
fractions_out <- tibble::rownames_to_column(fractions, "gene_id")
readr::write_tsv(fractions_out, file.path(dirs$tables, "variance_fractions.tsv"))

summary_table <- data.frame(
  variable = ordered_cols,
  median_fraction = vapply(ordered_cols, function(col) stats::median(fractions[[col]], na.rm = TRUE), numeric(1)),
  mean_fraction = vapply(ordered_cols, function(col) mean(fractions[[col]], na.rm = TRUE), numeric(1)),
  stringsAsFactors = FALSE
)
summary_table$median_percent <- 100 * summary_table$median_fraction
# Order covariates by median fraction (descending), Residuals always last.
covariate_order <- summary_table$variable[summary_table$variable != "Residuals"]
covariate_order <- covariate_order[order(-summary_table$median_fraction[match(covariate_order, summary_table$variable)])]
display_levels <- c(covariate_order, "Residuals")
summary_table <- summary_table[match(display_levels, summary_table$variable), , drop = FALSE]
readr::write_tsv(summary_table, file.path(dirs$tables, "variance_summary.tsv"))

long <- fractions_out %>%
  tidyr::pivot_longer(-gene_id, names_to = "variable", values_to = "fraction") %>%
  dplyr::mutate(percent = 100 * fraction, variable = factor(variable, levels = display_levels))

# Colours: Okabe-Ito (colourblind-safe) for covariates, neutral grey for Residuals.
okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#000000")
fill_values <- stats::setNames(rep(MID_GREY, length(display_levels)), display_levels)
covariate_levels <- display_levels[display_levels != "Residuals"]
if (length(covariate_levels)) {
  fill_values[covariate_levels] <- okabe_ito[(seq_along(covariate_levels) - 1L) %% length(okabe_ito) + 1L]
}

var_plot <- ggplot(long, aes(variable, percent, fill = variable)) +
  geom_violin(scale = "width", colour = "#3A4750", linewidth = 0.3, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", colour = NAVY, linewidth = 0.3) +
  scale_fill_manual(values = fill_values, guide = "none") +
  labs(
    title = "Variance partition across samples",
    subtitle = sprintf("Per-gene expression variance explained by each covariate (top %d variable genes)", n_keep),
    x = NULL, y = "% variance explained"
  ) +
  theme_publication(8.6) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot_pair(var_plot, file.path(dirs$figures, "variance_partition"), max(5.0, 1.1 * length(display_levels) + 2.0), 4.8)

write_json_file(
  list(
    project_id = cfg$project$id,
    method = "variancePartition::fitExtractVarPartModel",
    formula = paste(deparse(formula), collapse = " "),
    covariates = as.list(usable),
    genes = n_keep,
    median_percent = as.list(stats::setNames(summary_table$median_percent, summary_table$variable)),
    warnings = warnings
  ),
  file.path(args$outdir, "variance_partition_summary.json")
)
