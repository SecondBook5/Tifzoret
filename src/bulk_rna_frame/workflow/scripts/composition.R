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

# ---------------------------------------------------------------------------
# Panel B presentation parity. Replicates the finalized manuscript cell-state
# panel (reference make_cell_state_hybrid) onto the engine's already-computed
# `differential` columns. No statistics are refit and the cell_state_displayed.tsv
# table written above is untouched; only the ggplot construction is styled to
# match the reference: the two-tone lollipop (light fill + dark ink ring/stem),
# matched-gene point sizing, per-point FDR annotations, left-switched category
# facets, palette, legend and typography. The panel letter ("B") is drawn at
# assembly time by assemble.py's _label_overlay, so no plot.tag is baked in here
# (mirrors ontology.R). The reference's hand-picked ink hexes (CONTROL_INK
# #39799C / CAPE_INK #B55252) are study-editorial and are not carried by the
# engine's single configured fill palette, so a hue-preserving darkened
# companion is derived per condition to reproduce the two-tone appearance.
condition_fill <- condition_palette(cfg, c(denominator, numerator))
condition_ink <- vapply(condition_fill, function(hex) {
  hsv_values <- grDevices::rgb2hsv(grDevices::col2rgb(hex))
  grDevices::hsv(h = hsv_values[1L, ], s = pmin(1, hsv_values[2L, ] * 2), v = hsv_values[3L, ] * 0.68)
}, character(1))
differential <- differential %>%
  mutate(label = factor(label, levels = unique(label)), category_display = clean_term(category))
plot <- ggplot(differential, aes(logFC, label)) +
  geom_vline(xintercept = 0, linewidth = 0.48, colour = "#87939D") +
  geom_segment(aes(x = 0, xend = logFC, yend = label, colour = higher_in), linewidth = 0.8, lineend = "round") +
  geom_point(aes(fill = higher_in, colour = higher_in, size = matched_genes), shape = 21, stroke = 0.6) +
  geom_text(aes(label = paste0("FDR ", formatC(adj.P.Val, format = "e", digits = 1)), hjust = text_hjust),
            size = 2.2, colour = MID_GREY, show.legend = FALSE) +
  facet_grid(category_display ~ ., scales = "free_y", space = "free_y", switch = "y", drop = TRUE) +
  scale_colour_manual(values = condition_ink, breaks = c(denominator, numerator), drop = FALSE) +
  scale_fill_manual(values = condition_fill, breaks = c(denominator, numerator), drop = FALSE) +
  scale_size_continuous(range = c(2.5, 5.2), breaks = c(9, 10, 11)) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  labs(
    title = "Differential cell-state signatures",
    subtitle = sprintf("Positive scores are higher in %s; negative scores are higher in %s", numerator, denominator),
    x = "Signature score log2 fold-change", y = NULL,
    fill = "Higher in", colour = "Higher in", size = "Matched genes"
  ) +
  theme_publication(8.5) +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.spacing.x = grid::unit(2.4, "mm"),
    legend.margin = margin(1, 0, 6, 0),
    legend.box.margin = margin(0, 0, 2, 0),
    legend.title = element_text(size = rel(0.82)),
    legend.text = element_text(size = rel(0.78)),
    axis.title = element_text(size = rel(0.92)),
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.text = element_text(size = rel(0.82)),
    axis.text.y = element_text(margin = margin(r = 4)),
    panel.grid.major.y = element_blank(),
    panel.spacing.y = grid::unit(2.6, "mm"),
    strip.placement = "outside",
    strip.background = element_rect(fill = "#EEF2F5", colour = NA),
    strip.text.y.left = element_text(angle = 0, face = "bold", colour = NAVY, size = 7.1, margin = margin(3, 5, 3, 5)),
    plot.title = element_text(size = rel(1.13), margin = margin(b = 3)),
    plot.subtitle = element_text(margin = margin(b = 8)),
    plot.margin = margin(10, 14, 10, 12)
  ) +
  coord_cartesian(clip = "off") +
  guides(
    colour = "none",
    fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 3.4)),
    size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#D3DBE1", colour = NAVY))
  )
save_plot_pair(plot, file.path(dirs$figures, "cell_state_signatures"), 8.6, max(5.4, 0.42 * nrow(differential) + 2.5))

write_json_file(list(
  contrast_id = args[["contrast-id"]], signatures_tested = nrow(differential),
  warnings = list("Cell-state outputs are relative signature/composition scores and must not be interpreted as cell fractions.")
), file.path(args$outdir, "composition_summary.json"))
