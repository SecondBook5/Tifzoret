#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)
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
  refit <- dds; SummarizedExperiment::colData(refit) <- S4Vectors::DataFrame(augmented); DESeq2::design(refit) <- stats::as.formula(paste("~", design_terms)); refit <- DESeq2::DESeq(refit, quiet = TRUE)
  pattern <- paste0("^", factor_name, "_", numerator, "_vs_", denominator, "$")
  original_name <- grep(pattern, DESeq2::resultsNames(dds), value = TRUE)[1]; refit_name <- grep(pattern, DESeq2::resultsNames(refit), value = TRUE)[1]
  original <- DESeq2::results(dds, name = original_name); adjusted <- DESeq2::results(refit, name = refit_name)
  comparison <- data.frame(gene_id = rownames(original), log2_fold_change_original = original$log2FoldChange, adjusted_p_value_original = original$padj, log2_fold_change_sva = adjusted$log2FoldChange, adjusted_p_value_sva = adjusted$padj)
}
readr::write_tsv(sv_table, file.path(dirs$tables, "surrogate_variables.tsv")); readr::write_tsv(comparison, file.path(dirs$tables, "sva_de_sensitivity.tsv"), na = "NA")
write_json_file(list(contrast_id = args[["contrast-id"]], samples = nrow(metadata), surrogate_variables = n_sv, warnings = warnings), file.path(args$outdir, "sva_summary.json"))

