#!/usr/bin/env Rscript
# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     03_differential / de.R
# │ WHAT:      Primary differential expression per contrast (DESeq2 ~condition)
# │ WHY:       Identifies genes whose expression differs between conditions;
# │            the signed log2 fold change and FDR are the study's core result
# │ HOW:       DESeq2 negative-binomial GLM + Wald test; optional LFC shrinkage
# │            (ashr/normal/none); direction resolved via shared resolve_contrast()
# │ INPUTS:    inputs/{counts.tsv, samples.tsv, contrasts.tsv, annotation.tsv}
# │ PRODUCES:  de/tables/de_results.tsv, de/objects/deseq2.rds, volcano/MA figures
# │ CALLED BY: rule contrast_de (workflow/rules/core.smk)
# │ ENV:       workflow/envs/r.yaml
# └─────────────────────────────────────────────────────────────────

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "utils.R"), local = FALSE)
source(file.path(dirname(script_path), "de_render.R"), local = FALSE)

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

# Pairwise contrasts display only their two factor levels. A coefficient
# (interaction) contrast spans every group in the design, so it displays
# all samples and colours by the palette's canonical grouping (figures.group)
# rather than the label factor, which has no palette entries.
if (identical(resolved$type, "pairwise")) {
  display_samples <- metadata[[factor_name]] %in% c(denominator, numerator)
  display_group_col <- factor_name
  display_subtitle <- paste(numerator, "and", denominator, "samples")
} else {
  display_samples <- rep(TRUE, nrow(metadata))
  display_group_col <- cfg$figures$group
  display_subtitle <- "All design groups (interaction contrast)"
}

rendered <- render_de_outputs(
  result_table, cfg, dds, metadata, dirs,
  numerator = numerator, denominator = denominator,
  display_samples = display_samples,
  display_group_col = display_group_col,
  display_subtitle = display_subtitle,
  shrinkage_label = shrinkage_label, effect_axis = effect_axis,
  dispersion_fit = dispersion_fit, contrast_id = args[["contrast-id"]]
)
expression_transform <- rendered$expression_transform

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
