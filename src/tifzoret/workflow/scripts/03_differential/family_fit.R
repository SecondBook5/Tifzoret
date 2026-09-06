#!/usr/bin/env Rscript
# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     03_differential / family_fit.R
# │ WHAT:      The single Wald fit for one estimand family
# │ WHY:       Every estimand in a family must share one filter, one set of size
# │            factors, one set of dispersions, one sample universe, and one
# │            coefficient covariance -- structurally, not by convention
# │ HOW:       DESeq2 NB GLM over the family design; exports Sigma_beta, the
# │            compiled contrast matrix, the tested gene universe, diagnostics
# │ INPUTS:    inputs/{counts,samples,families,estimands,estimand_cell_weights}.tsv
# │ PRODUCES:  families/<id>/{objects/deseq2.rds, tables/*, figures/*, family_summary.json}
# │ CALLED BY: rule family_fit (workflow/rules/core.smk)
# │ ENV:       workflow/envs/r.yaml
# └─────────────────────────────────────────────────────────────────

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "utils.R"), local = FALSE)
source(file.path(dirname(script_path), "..", "estimands.R"), local = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(patchwork)
})

args <- parse_cli(c("project-config", "counts", "samples", "families", "estimands",
                    "cell-weights", "family-id", "outdir"))
cfg <- read_project(args[["project-config"]])
dirs <- ensure_output_dirs(args$outdir)

counts <- read_counts_contract(normalizePath(args$counts, mustWork = TRUE))
metadata <- readr::read_tsv(normalizePath(args$samples, mustWork = TRUE),
                            show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id

family <- read_family_spec(args$families, args[["family-id"]])
design_formula <- stats::as.formula(family$design)
for (field in all.vars(design_formula)) metadata[[field]] <- factor(metadata[[field]])
for (relevel_factor in names(family$reference_levels)) {
  metadata[[relevel_factor]] <- stats::relevel(
    factor(metadata[[relevel_factor]]), ref = family$reference_levels[[relevel_factor]])
}

# ---------------------------------------------------------------------------
# Gene universe. ONE filter for the family, exported as an artifact so DESeq2,
# edgeR, the pathway universe and the ORA universe all consume the same genes
# (spec §12). design_aware is the default because filterByExpr derives its
# thresholds from the design's group sizes rather than from a flat total.
# ---------------------------------------------------------------------------
design_matrix <- stats::model.matrix(design_formula, data = metadata)
keep <- switch(
  family$filter,
  design_aware = {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("filter: design_aware requires edgeR", call. = FALSE)
    }
    edgeR::filterByExpr(edgeR::DGEList(counts = counts), design_matrix)
  },
  total_count = rowSums(counts) >= 10L,
  none = rep(TRUE, nrow(counts)),
  stop("unknown filter: ", family$filter, call. = FALSE)
)
readr::write_tsv(
  data.frame(gene_id = rownames(counts), filter = family$filter, retained = unname(keep)),
  file.path(dirs$tables, "tested_gene_universe.tsv")
)

# ---------------------------------------------------------------------------
# The canonical fit. DESeq() (not the manual estimate* sequence) so Cook's
# outlier filtering/replacement behaves consistently; estimand.R reuses THIS
# object's dispersions and replaced counts.
# ---------------------------------------------------------------------------
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts[keep, , drop = FALSE], colData = metadata, design = design_formula)
fitted <- fit_deseq_with_dispersion_fallback(dds)
dds <- fitted$dds
dispersion_fit <- fitted$dispersion_fit
saveRDS(dds, file.path(dirs$objects, "deseq2.rds"))
coefficient_names <- DESeq2::resultsNames(dds)

# ---------------------------------------------------------------------------
# Compiled contrast matrix C. One row per estimand, carrying BOTH the source
# expression and the numeric vector so provenance is auditable from the file.
# ---------------------------------------------------------------------------
estimand_ids <- read_estimand_ids(args$estimands, family$family_id)
contrast_rows <- lapply(estimand_ids, function(estimand_id) {
  spec <- read_estimand_spec(args$estimands, estimand_id)
  weights <- read_cell_weights(args[["cell-weights"]], estimand_id)
  vector <- compile_contrast_vector(
    c(family, spec), weights, metadata, design_formula, coefficient_names)
  row <- data.frame(
    estimand_id = estimand_id,
    label = spec$label,
    role = spec$role,
    atom_kind = spec$atom_kind,
    expression = spec$expression,
    cell_weights = paste(sprintf("%s=%g", names(weights), weights), collapse = ";"),
    stringsAsFactors = FALSE
  )
  for (index in seq_along(coefficient_names)) row[[coefficient_names[[index]]]] <- vector[[index]]
  row
})
readr::write_tsv(dplyr::bind_rows(contrast_rows), file.path(dirs$tables, "contrast_matrix.tsv"))

# ---------------------------------------------------------------------------
# Coefficient covariance. Sigma_beta = (X' W X)^-1 per gene, from DESeq2's
# stored per-gene weights. Exported long-form for the top genes by base mean so
# a downstream table's SE can be recomputed from the file (spec §11 test 3).
# ---------------------------------------------------------------------------
covariance_limit <- cfg$analysis$settings$family$covariance_genes
if (is.null(covariance_limit)) covariance_limit <- 2000L
base_mean <- S4Vectors::mcols(dds)$baseMean
selected <- utils::head(order(base_mean, decreasing = TRUE), as.integer(covariance_limit))
gene_ids <- rownames(dds)
covariance_of <- function(index) {
  # DESeq2's Wald SEs come from the same (X' W X)^-1 with W the NB IRLS weights.
  mu <- SummarizedExperiment::assays(dds)[["mu"]][index, ]
  dispersion <- DESeq2::dispersions(dds)[[index]]
  if (any(is.na(mu)) || is.na(dispersion)) return(NULL)
  weight <- 1 / (1 / mu + dispersion)
  scaled <- design_matrix * sqrt(weight)
  solved <- tryCatch(solve(t(scaled) %*% scaled), error = function(e) NULL)
  if (is.null(solved)) return(NULL)
  # Convert from the natural-log scale DESeq2 fits on to log2, matching lfcSE.
  solved / (log(2)^2)
}
covariance_frames <- lapply(selected, function(index) {
  sigma <- covariance_of(index)
  if (is.null(sigma)) return(NULL)
  colnames(sigma) <- coefficient_names
  rownames(sigma) <- coefficient_names
  expand <- expand.grid(
    coefficient_row = coefficient_names, coefficient_col = coefficient_names,
    stringsAsFactors = FALSE)
  expand$gene_id <- gene_ids[[index]]
  expand$covariance <- mapply(
    function(r, c) sigma[r, c], expand$coefficient_row, expand$coefficient_col)
  expand[, c("gene_id", "coefficient_row", "coefficient_col", "covariance")]
})
readr::write_tsv(
  dplyr::bind_rows(covariance_frames), file.path(dirs$tables, "coefficient_covariance.tsv"))

# ---------------------------------------------------------------------------
# Design diagnostics. WARN, never fail: the hard fails already ran at config
# time in config/estimands.py, before any compute (spec §6).
# ---------------------------------------------------------------------------
singular <- svd(design_matrix)$d
rank <- qr(design_matrix)$rank
condition_number <- if (min(singular) > 0) max(singular) / min(singular) else Inf
hat <- diag(design_matrix %*% MASS_safe_solve(design_matrix) %*% t(design_matrix))
residual_df <- nrow(design_matrix) - rank
cell_keys <- apply(metadata[, family$cells, drop = FALSE], 1,
                   function(row) paste(row, collapse = CELL_KEY_SEPARATOR))
cell_sizes <- as.data.frame(table(cell_keys), stringsAsFactors = FALSE)
names(cell_sizes) <- c("label", "value")
leverage_threshold <- 3 * rank / nrow(design_matrix)
warnings <- character(0)
if (nrow(cell_sizes) && min(cell_sizes$value) > 0 &&
    max(cell_sizes$value) / min(cell_sizes$value) >= 3) {
  warnings <- c(warnings, sprintf(
    "cell sizes are imbalanced (max/min = %.1f)",
    max(cell_sizes$value) / min(cell_sizes$value)))
}
if (is.finite(condition_number) && condition_number > 30) {
  warnings <- c(warnings, sprintf("design condition number is %.1f", condition_number))
}
high_leverage <- names(hat)[hat > leverage_threshold]
if (length(high_leverage)) {
  warnings <- c(warnings, sprintf("high-leverage samples (h_ii > 3p/n): %s",
                                  paste(high_leverage, collapse = ", ")))
}
if (residual_df < 4) {
  warnings <- c(warnings, sprintf("residual degrees of freedom is only %d", residual_df))
}
diagnostics <- dplyr::bind_rows(
  data.frame(metric = "rank", label = "", value = rank, stringsAsFactors = FALSE),
  data.frame(metric = "residual_df", label = "", value = residual_df, stringsAsFactors = FALSE),
  data.frame(metric = "condition_number", label = "", value = condition_number, stringsAsFactors = FALSE),
  data.frame(metric = "leverage_threshold", label = "", value = leverage_threshold, stringsAsFactors = FALSE),
  data.frame(metric = "singular_value", label = as.character(seq_along(singular)),
             value = singular, stringsAsFactors = FALSE),
  data.frame(metric = "leverage", label = rownames(design_matrix), value = unname(hat),
             stringsAsFactors = FALSE),
  data.frame(metric = "cell_size", label = cell_sizes$label, value = as.numeric(cell_sizes$value),
             stringsAsFactors = FALSE)
)
readr::write_tsv(diagnostics, file.path(dirs$tables, "design_diagnostics.tsv"))

cell_plot <- ggplot2::ggplot(cell_sizes, ggplot2::aes(label, value)) +
  ggplot2::geom_col(fill = "#6C92AE", colour = "white", linewidth = 0.15) +
  ggplot2::labs(title = "Design cell sizes", x = NULL, y = "Samples") +
  theme_publication(8.6) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
leverage_plot <- ggplot2::ggplot(
  data.frame(sample_id = rownames(design_matrix), leverage = unname(hat)),
  ggplot2::aes(stats::reorder(sample_id, leverage), leverage)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_hline(yintercept = leverage_threshold, linetype = "dashed", colour = "#B22222") +
  ggplot2::labs(title = "Per-sample leverage", subtitle = "Dashed line marks 3p/n",
       x = NULL, y = expression(h[ii])) +
  theme_publication(8.6) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
save_plot_pair(cell_plot / leverage_plot, file.path(dirs$figures, "design_diagnostics"), 7.0, 6.4)

write_json_file(
  list(
    project_id = cfg$project$id,
    family_id = family$family_id,
    design = family$design,
    cells = as.list(family$cells),
    declared = family$declared,
    shrinkage = family$shrinkage,
    filter = family$filter,
    coefficients = as.list(coefficient_names),
    estimands = as.list(estimand_ids),
    samples = ncol(dds),
    genes_tested = nrow(dds),
    genes_filtered_out = sum(!keep),
    rank = rank,
    residual_df = residual_df,
    condition_number = condition_number,
    dispersion_fit = dispersion_fit,
    covariance_genes = length(covariance_frames),
    warnings = as.list(warnings)
  ),
  file.path(args$outdir, "family_summary.json")
)
