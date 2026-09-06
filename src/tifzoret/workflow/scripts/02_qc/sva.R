#!/usr/bin/env Rscript
# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     02_qc / sva.R
# │ WHAT:      Surrogate variable analysis (latent technical confounders)
# │ WHY:       Estimates hidden batch effects; compares DE with/without SVA adjustment
# │            to quantify sensitivity (exploratory, per-contrast, opt-in)
# │ HOW:       sva::sva on VST expression; refit DESeq2 design augmented with SVs
# │ INPUTS:    families/<id>/objects/deseq2.rds, qc/objects/vst.rds, contrasts.tsv, config
# │ PRODUCES:  advanced/sva/{surrogate_variables.tsv, sva_de_sensitivity.tsv}
# │ CALLED BY: rule contrast_sva (workflow/rules/advanced.smk)
# │ ENV:       workflow/envs/r.yaml
# └─────────────────────────────────────────────────────────────────

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "utils.R"), local = FALSE)
suppressPackageStartupMessages({library(DESeq2); library(sva)})

args <- parse_cli(c("project-config", "dds", "vst", "contrasts", "contrast-id", "outdir"))
cfg <- read_project(args[["project-config"]]); dirs <- ensure_output_dirs(args$outdir)
dds <- readRDS(args$dds); vst <- readRDS(args$vst); metadata <- as.data.frame(SummarizedExperiment::colData(dds))
contrasts <- readr::read_tsv(args$contrasts, show_col_types = FALSE); contrast <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
factor_name <- contrast$factor[[1]]; numerator <- contrast$numerator[[1]]; denominator <- contrast$denominator[[1]]
metadata[[factor_name]] <- stats::relevel(factor(metadata[[factor_name]]), ref = denominator)
full_formula <- stats::as.formula(cfg$design$formula)
null_terms <- setdiff(all.vars(full_formula), factor_name)
null_formula <- stats::as.formula(if (length(null_terms)) paste("~", paste(null_terms, collapse = " + ")) else "~ 1")
mod <- stats::model.matrix(full_formula, metadata); mod0 <- stats::model.matrix(null_formula, metadata)
expression <- SummarizedExperiment::assay(vst)
set.seed(if (is.null(cfg$analysis$random_seed)) 1L else as.integer(cfg$analysis$random_seed))
n_sv <- tryCatch(sva::num.sv(expression, mod, method = "be"), error = function(error) 0L)
minimum_n <- if (is.null(cfg$analysis$settings$sva$minimum_recommended_samples)) 10L else as.integer(cfg$analysis$settings$sva$minimum_recommended_samples)
warnings <- list(); if (nrow(metadata) < minimum_n) warnings <- c(warnings, sprintf("SVA n=%d is below the configured recommended minimum n=%d; sensitivity results are exploratory.", nrow(metadata), minimum_n))
sv_table <- data.frame(sample_id = rownames(metadata))
comparison <- data.frame()
if (n_sv > 0L) {
  estimate <- sva::sva(expression, mod, mod0, n.sv = n_sv)
  sv_table <- cbind(sv_table, as.data.frame(estimate$sv)); names(sv_table)[-1] <- paste0("SV", seq_len(n_sv))
  augmented <- cbind(metadata, estimate$sv); names(augmented)[(ncol(metadata) + 1):ncol(augmented)] <- paste0("SV", seq_len(n_sv))
  design_terms <- paste(c(sub("^~", "", cfg$design$formula), paste0("SV", seq_len(n_sv))), collapse = " + ")
  refit <- dds; SummarizedExperiment::colData(refit) <- S4Vectors::DataFrame(augmented); DESeq2::design(refit) <- stats::as.formula(paste("~", design_terms))
  # Share the engine's dispersion fallback rather than calling DESeq() bare. The
  # SV-augmented refit hits the same parametric-trend failure the canonical fit
  # guards against ("all gene-wise dispersion estimates are within 2 orders of
  # magnitude..."), and because this is an exploratory sensitivity module, an
  # unguarded refit aborted the ENTIRE workflow on data the main DE path handled
  # fine.
  refit <- fit_deseq_with_dispersion_fallback(refit)$dds
  # Resolve the comparison by contrast rather than by coefficient name. The
  # family fit keeps ONE reference for every estimand it serves, so a coefficient
  # named "<factor>_<numerator>_vs_<denominator>" need not exist -- only the
  # retired de.R, which releveled per contrast, could rely on that. This form is
  # exact under any reference level and carries the requested direction.
  contrast_spec <- c(factor_name, numerator, denominator)
  original <- DESeq2::results(dds, contrast = contrast_spec); adjusted <- DESeq2::results(refit, contrast = contrast_spec)
  comparison <- data.frame(gene_id = rownames(original), log2_fold_change_original = original$log2FoldChange, adjusted_p_value_original = original$padj, log2_fold_change_sva = adjusted$log2FoldChange, adjusted_p_value_sva = adjusted$padj)
}
readr::write_tsv(sv_table, file.path(dirs$tables, "surrogate_variables.tsv")); readr::write_tsv(comparison, file.path(dirs$tables, "sva_de_sensitivity.tsv"), na = "NA")
write_json_file(list(contrast_id = args[["contrast-id"]], samples = nrow(metadata), surrogate_variables = n_sv, warnings = warnings), file.path(args$outdir, "sva_summary.json"))

