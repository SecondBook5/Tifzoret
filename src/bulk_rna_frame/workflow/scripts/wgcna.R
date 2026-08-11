#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

args <- parse_cli(c("project-config", "vst", "samples", "annotation", "contrasts", "contrast-id", "outdir"))
cfg <- read_project(args[["project-config"]]); dirs <- ensure_output_dirs(args$outdir)
expression <- SummarizedExperiment::assay(readRDS(args$vst)); annotation <- read_annotation_contract(args$annotation)
metadata <- readr::read_tsv(args$samples, show_col_types = FALSE) %>% as.data.frame(); rownames(metadata) <- metadata$sample_id
contrasts <- readr::read_tsv(args$contrasts, show_col_types = FALSE); contrast <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
factor_name <- contrast$factor[[1]]; numerator <- contrast$numerator[[1]]; denominator <- contrast$denominator[[1]]
top_n <- if (is.null(cfg$analysis$settings$wgcna$top_variable_genes)) 5000L else as.integer(cfg$analysis$settings$wgcna$top_variable_genes)
min_module <- if (is.null(cfg$analysis$settings$wgcna$minimum_module_size)) 15L else as.integer(cfg$analysis$settings$wgcna$minimum_module_size)
minimum_n <- if (is.null(cfg$analysis$settings$wgcna$minimum_recommended_samples)) 15L else as.integer(cfg$analysis$settings$wgcna$minimum_recommended_samples)
warnings <- list(); if (ncol(expression) < minimum_n) warnings <- c(warnings, sprintf("WGCNA n=%d is below the configured recommended minimum n=%d; modules are exploratory.", ncol(expression), minimum_n))
keep <- head(order(apply(expression, 1, stats::var), decreasing = TRUE), min(top_n, nrow(expression))); datExpr <- t(expression[keep, , drop = FALSE])
good <- WGCNA::goodSamplesGenes(datExpr, verbose = 0); datExpr <- datExpr[good$goodSamples, good$goodGenes, drop = FALSE]
powers <- 1:20; fit <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)$fitIndices
candidate <- fit$Power[fit$SFT.R.sq >= 0.8]; power <- if (length(candidate)) min(candidate) else fit$Power[which.max(fit$SFT.R.sq)]
net <- WGCNA::blockwiseModules(datExpr, power = power, networkType = "signed", TOMType = "signed", minModuleSize = min_module, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, maxBlockSize = ncol(datExpr) + 1, verbose = 0)
colors <- WGCNA::labels2colors(net$colors); eigengenes <- WGCNA::orderMEs(WGCNA::moduleEigengenes(datExpr, colors)$eigengenes)
trait <- as.numeric(factor(metadata[rownames(datExpr), factor_name], levels = c(denominator, numerator))) - 1
correlation <- stats::cor(eigengenes, trait, use = "pairwise.complete.obs"); pvalues <- WGCNA::corPvalueStudent(correlation, nrow(datExpr))
trait_table <- data.frame(module = rownames(correlation), correlation = correlation[, 1], p_value = pvalues[, 1]); readr::write_tsv(trait_table, file.path(dirs$tables, "wgcna_module_trait.tsv"))
membership <- data.frame(gene_id = colnames(datExpr), gene_symbol = annotation$gene_symbol[match(colnames(datExpr), annotation$gene_id)], module = colors)
kme <- WGCNA::signedKME(datExpr, eigengenes); membership$kME <- apply(abs(kme), 1, max, na.rm = TRUE); readr::write_tsv(membership, file.path(dirs$tables, "wgcna_membership.tsv"))
hubs <- membership %>% filter(module != "grey") %>% group_by(module) %>% slice_max(kME, n = 20, with_ties = FALSE) %>% ungroup(); readr::write_tsv(hubs, file.path(dirs$tables, "wgcna_hubs.tsv"))
plot_table <- trait_table %>% mutate(module = factor(module, levels = rev(module)))
plot <- ggplot(plot_table, aes("Condition", module, fill = correlation)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f\nP %.2g", correlation, p_value)), size = 3) + scale_fill_gradient2(low = "#39799C", mid = "white", high = "#B55252", midpoint = 0, limits = c(-1, 1)) + labs(title = "WGCNA module–condition associations", x = NULL, y = NULL, fill = "Correlation") + theme_publication(8)
save_plot_pair(plot, file.path(dirs$figures, "wgcna_module_trait"), 6.2, max(4.5, 0.35 * nrow(trait_table) + 2))
write_json_file(list(contrast_id = args[["contrast-id"]], samples = nrow(datExpr), genes = ncol(datExpr), soft_power = power, modules = length(setdiff(unique(colors), "grey")), warnings = warnings), file.path(args$outdir, "wgcna_summary.json"))

