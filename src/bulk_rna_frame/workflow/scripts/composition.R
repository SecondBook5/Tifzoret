#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(GSVA)
  library(limma)
})

args <- parse_cli(c("project-config", "samples", "annotation", "contrasts", "contrast-id", "vst", "signatures", "outdir"))
cfg <- read_project(args[["project-config"]])
dirs <- ensure_output_dirs(args$outdir)
metadata <- readr::read_tsv(args$samples, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
annotation <- read_annotation_contract(args$annotation)
contrasts <- readr::read_tsv(args$contrasts, show_col_types = FALSE, progress = FALSE)
contrast <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
if (nrow(contrast) != 1L) stop("Could not resolve contrast", call. = FALSE)
# Route direction resolution through the shared resolver used by de.R/pathways.R
# so the relevel reference and numerator identity come from one code path. This
# stage only ever runs on pairwise contrasts (PAIRWISE_CONTRAST_IDS gate), so the
# guard is a documented assertion and cannot fire in practice.
resolved <- resolve_contrast(contrast, cfg$design$formula)
if (!identical(resolved$type, "pairwise")) stop("composition stage supports pairwise contrasts only", call. = FALSE)
factor_name <- resolved$factor_name; numerator <- resolved$numerator; denominator <- resolved$denominator

definition <- yaml::read_yaml(args$signatures)$signatures
signature_table <- bind_rows(lapply(definition, function(item) data.frame(
  signature_id = item$id, label = item$label, category = item$category,
  description = if (is.null(item$description)) "" else item$description,
  gene_symbol = unlist(item$genes), stringsAsFactors = FALSE
)))
vst <- readRDS(args$vst)
expression <- matrix_to_symbols(SummarizedExperiment::assay(vst), annotation)
metadata <- metadata[match(colnames(expression), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id
symbol_lookup <- setNames(rownames(expression), toupper(rownames(expression)))
signature_table <- signature_table %>%
  mutate(configured_gene_symbol = gene_symbol, gene_symbol = unname(symbol_lookup[toupper(gene_symbol)])) %>%
  filter(!is.na(gene_symbol))
min_genes <- if (is.null(cfg$analysis$settings$composition$min_genes)) 3L else as.integer(cfg$analysis$settings$composition$min_genes)
gene_sets <- split(signature_table$gene_symbol, signature_table$signature_id)
gene_sets <- lapply(gene_sets, function(genes) intersect(unique(genes), rownames(expression)))
matched <- data.frame(signature_id = names(gene_sets), matched_genes = lengths(gene_sets))
gene_sets <- gene_sets[lengths(gene_sets) >= min_genes]
if (!length(gene_sets)) stop("No cell-state signatures retain the configured minimum number of measured genes", call. = FALSE)

parameter <- GSVA::ssgseaParam(exprData = expression, geneSets = gene_sets, minSize = min_genes, maxSize = Inf, normalize = TRUE)
scores <- GSVA::gsva(parameter, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
score_table <- as.data.frame(scores, check.names = FALSE) %>% tibble::rownames_to_column("signature_id")
readr::write_tsv(score_table, file.path(dirs$tables, "cell_state_scores.tsv"))

metadata$contrast_group <- stats::relevel(factor(metadata[[factor_name]]), ref = resolved$reference_levels[[factor_name]])
formula_text <- gsub(paste0("\\b", factor_name, "\\b"), "contrast_group", cfg$design$formula)
design <- stats::model.matrix(stats::as.formula(formula_text), metadata)
coefficient <- grep(paste0("^contrast_group", make.names(numerator), "$"), colnames(design), value = TRUE)
if (length(coefficient) != 1L) stop("Could not resolve cell-state model coefficient", call. = FALSE)
fit <- limma::eBayes(limma::lmFit(scores, design))
differential <- limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("signature_id") %>%
  left_join(signature_table %>% select(signature_id, label, category, description) %>% distinct(), by = "signature_id") %>%
  left_join(matched, by = "signature_id") %>%
  mutate(
    contrast_id = args[["contrast-id"]], numerator = numerator, denominator = denominator,
    higher_in = ifelse(logFC >= 0, numerator, denominator),
    fdr_label = paste0("FDR ", formatC(adj.P.Val, format = "g", digits = 2)),
    text_hjust = ifelse(logFC >= 0, -0.12, 1.12)
  ) %>% arrange(category, logFC)
readr::write_tsv(differential, file.path(dirs$tables, "cell_state_differential.tsv"), na = "NA")
readr::write_tsv(differential, file.path(dirs$tables, "cell_state_displayed.tsv"), na = "NA")

palette <- condition_palette(cfg, c(denominator, numerator))
differential <- differential %>% mutate(label = factor(label, levels = unique(label)))
plot <- ggplot(differential, aes(logFC, label, colour = higher_in, size = matched_genes)) +
  geom_vline(xintercept = 0, colour = "#87939C", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = logFC, yend = label), linewidth = 0.75, show.legend = FALSE) +
  geom_point(alpha = 0.96) +
  geom_text(aes(label = fdr_label, hjust = text_hjust), colour = MID_GREY, size = 2.5, show.legend = FALSE) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = palette, breaks = c(denominator, numerator), drop = FALSE, name = "Higher in") +
  scale_size_continuous(range = c(2.8, 7), name = "Matched genes") +
  scale_x_continuous(expand = expansion(mult = c(0.22, 0.22))) +
  labs(title = "Cell-state signature shifts", subtitle = "Relative signature scores; positive effects are higher in the contrast numerator", x = "Signature score log2 fold-change", y = NULL) +
  theme_publication(8.3) +
  theme(legend.position = "top", strip.placement = "outside", strip.background = element_rect(fill = "#F2F5F7", colour = NA), strip.text.y.left = element_text(face = "bold"), panel.grid.major.y = element_blank())
save_plot_pair(plot, file.path(dirs$figures, "cell_state_signatures"), 8.6, max(5.4, 0.42 * nrow(differential) + 2.5))

write_json_file(list(
  contrast_id = args[["contrast-id"]], signatures_tested = nrow(differential),
  warnings = list("Cell-state outputs are relative signature/composition scores and must not be interpreted as cell fractions.")
), file.path(args$outdir, "composition_summary.json"))
