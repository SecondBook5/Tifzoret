#!/usr/bin/env Rscript
# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     03_differential / estimand.R
# │ WHAT:      One estimand extracted from its family's shared Wald fit
# │ WHY:       Every reported SE is the exact sqrt(c' Sigma c); the naive
# │            independent-sum formula it replaces inflated SE ~1.8x
# │ HOW:       results(dds, contrast = c) for inference; apeglm on a
# │            reparameterized refit for the display effect only
# │ INPUTS:    families/<id>/objects/deseq2.rds, inputs/*, family contrast matrix
# │ PRODUCES:  contrasts/<estimand_id>/analyses/de/{tables,figures}/*, de_summary.json
# │ CALLED BY: rule family_estimand (workflow/rules/core.smk)
# │ ENV:       workflow/envs/r.yaml
# └─────────────────────────────────────────────────────────────────

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "utils.R"), local = FALSE)
source(file.path(dirname(script_path), "..", "estimands.R"), local = FALSE)
source(file.path(dirname(script_path), "de_render.R"), local = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(patchwork)
})

args <- parse_cli(c("project-config", "counts", "samples", "annotation", "families",
                    "estimands", "cell-weights", "family-dir", "estimand-id", "outdir"))
cfg <- read_project(args[["project-config"]])
dirs <- ensure_output_dirs(args$outdir)
annotation <- read_annotation_contract(normalizePath(args$annotation, mustWork = TRUE))

spec <- read_estimand_spec(args$estimands, args[["estimand-id"]])
family <- read_family_spec(args$families, spec$family_id)
weights <- read_cell_weights(args[["cell-weights"]], spec$estimand_id)
design_formula <- stats::as.formula(family$design)

dds <- readRDS(file.path(args[["family-dir"]], "objects", "deseq2.rds"))
metadata <- as.data.frame(SummarizedExperiment::colData(dds))
coefficient_names <- DESeq2::resultsNames(dds)
contrast_vector <- compile_contrast_vector(
  c(family, spec), weights, metadata, design_formula, coefficient_names)

fdr <- cfg$figures$de$fdr
raw <- DESeq2::results(dds, contrast = as.numeric(contrast_vector), alpha = fdr)

# ---------------------------------------------------------------------------
# Display shrinkage. apeglm refuses a contrast ("only for use with 'coef'"), but
# an estimand that is a contrast under one parameterization is usually a plain
# coefficient under another, and releveling does not change the model's column
# space -- dispersions, size factors, fitted means and deviance are identical,
# only the coefficient basis moves. So: find a releveling that turns c into a
# single coefficient, re-run ONLY nbinomWaldTest on the same object (dispersions
# and replaced counts inherited), shrink that coefficient, and assert the
# reparameterized MLE reproduces the canonical contrast before trusting it.
#
# Inference NEVER comes from the reparameterized fit. If the identity assertion
# fails we fall back to ashr (which accepts contrasts) rather than silently
# reporting a display effect from a model we could not verify.
# ---------------------------------------------------------------------------
single_coefficient <- function(vector) {
  nonzero <- which(abs(vector) > 1e-9)
  if (length(nonzero) == 1L && abs(abs(vector[[nonzero]]) - 1) < 1e-9) {
    return(list(name = names(vector)[[nonzero]], sign = sign(vector[[nonzero]])))
  }
  NULL
}

reparameterize <- function(cells, cell_keys) {
  # Try each cell as the reference for every cell factor; the first combination
  # that renders the estimand a single +/-1 coefficient wins. Deterministic.
  candidates <- list()
  parts <- strsplit(cell_keys, CELL_KEY_SEPARATOR, fixed = TRUE)
  for (piece in parts) {
    references <- stats::setNames(as.list(piece), cells)
    candidates[[length(candidates) + 1L]] <- references
  }
  for (references in candidates) {
    local_metadata <- metadata
    for (column in names(references)) {
      local_metadata[[column]] <- stats::relevel(
        factor(local_metadata[[column]]), ref = references[[column]])
    }
    candidate_dds <- dds
    SummarizedExperiment::colData(candidate_dds) <- S4Vectors::DataFrame(local_metadata)
    DESeq2::design(candidate_dds) <- design_formula
    if (!identical(DESeq2::dispersions(candidate_dds), DESeq2::dispersions(dds))) next
    candidate_dds <- DESeq2::nbinomWaldTest(candidate_dds, quiet = TRUE)
    names_here <- DESeq2::resultsNames(candidate_dds)
    vector <- compile_contrast_vector(
      c(family, spec), weights, local_metadata, design_formula, names_here)
    single <- single_coefficient(vector)
    if (!is.null(single)) {
      return(list(dds = candidate_dds, coefficient = single$name, sign = single$sign))
    }
  }
  NULL
}

requested <- family$shrinkage
shrinkage_applied <- "none"
reparameterized <- FALSE
identity_gap_abs <- 0
identity_gap_relative <- 0
shrunk_lfc <- raw$log2FoldChange

if (identical(requested, "apeglm")) {
  attempt <- if (identical(spec$atom_kind, "coef")) {
    single <- single_coefficient(contrast_vector)
    if (is.null(single)) NULL else list(dds = dds, coefficient = single$name, sign = single$sign)
  } else {
    reparameterize(family$cells, names(weights))
  }
  if (!is.null(attempt)) {
    check <- DESeq2::results(attempt$dds, name = attempt$coefficient)
    aligned <- attempt$sign * check$log2FoldChange
    # Identity assertion: the reparameterized MLE must reproduce the canonical
    # contrast to within numerical precision. We check both absolute difference
    # and difference in units of SE; the latter is scale-free and the binding
    # constraint. A thousandth of a standard error cannot change inference, and
    # a genuine mis-mapping is ~1 SE, so this catches errors with thousandfold
    # margin while tolerating IRLS convergence noise.
    identity_gap_abs <- max(abs(aligned - raw$log2FoldChange), na.rm = TRUE)
    valid_se <- is.finite(raw$lfcSE) & raw$lfcSE > 0
    if (any(valid_se)) {
      relative_diff <- abs(aligned - raw$log2FoldChange)[valid_se] / raw$lfcSE[valid_se]
      identity_gap_relative <- max(relative_diff, na.rm = TRUE)
    } else {
      identity_gap_relative <- Inf
    }
    if (is.finite(identity_gap_relative) && identity_gap_relative < 1e-3) {
      shrunk <- DESeq2::lfcShrink(
        attempt$dds, coef = attempt$coefficient, type = "apeglm", quiet = TRUE)
      shrunk_lfc <- attempt$sign * shrunk$log2FoldChange
      shrinkage_applied <- "apeglm"
      reparameterized <- TRUE
    } else {
      requested <- "ashr"
    }
  } else {
    requested <- "ashr"
  }
}
if (identical(requested, "ashr") && !identical(shrinkage_applied, "apeglm")) {
  if (requireNamespace("ashr", quietly = TRUE)) {
    shrunk <- DESeq2::lfcShrink(
      dds, contrast = as.numeric(contrast_vector), type = "ashr", quiet = TRUE)
    shrunk_lfc <- shrunk$log2FoldChange
    shrinkage_applied <- "ashr"
  } else {
    tz_degrade(cfg, "ashr unavailable; reporting the unshrunken MLE effect")
  }
}
shrinkage_label <- if (identical(shrinkage_applied, "none")) {
  "no shrinkage (raw MLE)"
} else {
  paste0(shrinkage_applied, " shrinkage")
}
effect_axis <- if (identical(shrinkage_applied, "none")) {
  "log2 fold-change"
} else {
  "Shrunken log2 fold-change"
}

result_table <- data.frame(
  gene_id = rownames(raw),
  base_mean = raw$baseMean,
  log2_fold_change_raw = raw$log2FoldChange,
  log2_fold_change = shrunk_lfc,
  lfc_se = raw$lfcSE,
  statistic = raw$stat,
  p_value = raw$pvalue,
  adjusted_p_value = raw$padj,
  ci_low = raw$log2FoldChange - stats::qnorm(0.975) * raw$lfcSE,
  ci_high = raw$log2FoldChange + stats::qnorm(0.975) * raw$lfcSE,
  shrinkage = shrinkage_applied,
  family_id = family$family_id,
  estimand_id = spec$estimand_id,
  estimand_role = spec$role,
  stringsAsFactors = FALSE
) %>%
  left_join(annotation, by = "gene_id") %>%
  mutate(
    gene_symbol = ifelse(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol),
    safe_p_value = ifelse(is.na(p_value), NA_real_, pmax(p_value, .Machine$double.xmin)),
    negative_log10_p = -log10(safe_p_value),
    direction = case_when(
      !is.na(adjusted_p_value) & adjusted_p_value < fdr &
        log2_fold_change >= cfg$figures$de$abs_log2fc ~ "up_in_numerator",
      !is.na(adjusted_p_value) & adjusted_p_value < fdr &
        log2_fold_change <= -cfg$figures$de$abs_log2fc ~ "down_in_numerator",
      TRUE ~ "not_significant"
    ),
    significance_class = classify_significance(
      adjusted_p_value, log2_fold_change, fdr, cfg$figures$de$abs_log2fc),
    contrast_id = spec$estimand_id,
    numerator = spec$label,
    denominator = ""
  ) %>%
  arrange(adjusted_p_value, desc(abs(log2_fold_change)))
readr::write_tsv(result_table, file.path(dirs$tables, "de_results.tsv"), na = "NA")

display_samples <- display_samples_for_cells(metadata, family$cells, names(weights))
display_group_col <- if (length(family$cells)) family$cells[[1]] else cfg$figures$group
rendered <- render_de_outputs(
  result_table, cfg, dds, metadata, dirs,
  numerator = spec$label, denominator = "",
  display_samples = display_samples,
  display_group_col = display_group_col,
  display_subtitle = sprintf("Samples in the %d referenced design cell(s)", length(weights)),
  shrinkage_label = shrinkage_label, effect_axis = effect_axis,
  dispersion_fit = "parametric", contrast_id = spec$estimand_id
)

write_json_file(
  list(
    project_id = cfg$project$id,
    contrast_id = spec$estimand_id,
    estimand_id = spec$estimand_id,
    family_id = family$family_id,
    label = spec$label,
    role = spec$role,
    expression = spec$expression,
    contrast_vector = as.list(stats::setNames(as.numeric(contrast_vector), coefficient_names)),
    contrast_type = "estimand",
    shrinkage_requested = family$shrinkage,
    shrinkage_applied = shrinkage_applied,
    reparameterized = reparameterized,
    identity_max_abs_difference = identity_gap_abs,
    identity_max_relative_difference = identity_gap_relative,
    samples = ncol(dds),
    genes_tested = nrow(result_table),
    significant_up = sum(result_table$direction == "up_in_numerator"),
    significant_down = sum(result_table$direction == "down_in_numerator")
  ),
  file.path(args$outdir, "de_summary.json")
)
