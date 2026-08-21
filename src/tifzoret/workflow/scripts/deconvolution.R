#!/usr/bin/env Rscript

# Signature-matrix cell-fraction deconvolution. Given a user-supplied signature
# matrix (genes x reference cell types, resources.deconvolution_signature),
# estimate each bulk sample's cell-type composition by non-negative least squares
# (NNLS): solve min_x ||S x - b||^2 subject to x >= 0 for every sample's linear
# expression vector b against the signature columns S, then normalise the
# non-negative coefficients to fractions summing to one. This is the classic
# Abbas/DeconRNASeq baseline -- deterministic, dependency-light, and exploratory.
# It assumes the signature and the mixture are on comparable LINEAR scales; we
# feed it CPM (library-size-normalised counts, linear), and report per-sample
# reconstruction correlation so a weak fit is visible rather than hidden.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(nnls)
})

# Okabe-Ito qualitative palette (colourblind-safe); cell types are assigned hues
# in signature-column order and cycled only if a signature declares more than
# eight cell types (rare for a curated signature).
CELL_TYPE_PALETTE <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#7A7A7A"
)

args <- parse_cli(c("project-config", "counts", "samples", "annotation", "signature", "outdir"))
cfg <- read_project(args[["project-config"]])
cfg$.counts <- normalizePath(args$counts, mustWork = TRUE)
cfg$.samples <- normalizePath(args$samples, mustWork = TRUE)
cfg$.annotation <- normalizePath(args$annotation, mustWork = TRUE)
dirs <- ensure_output_dirs(args$outdir)

settings <- cfg$analysis$settings$deconvolution
min_genes <- if (is.null(settings$min_genes)) 2L else as.integer(settings$min_genes)

# Emit empty-but-well-formed outputs plus a placeholder figure and a summary that
# records why, so an unusable signature (too little gene overlap, a malformed
# matrix, or degenerate library sizes) degrades gracefully instead of aborting
# the whole run -- matching the graceful-skip contract of the other opt-in
# modules (see spia.R). deconvolution is an opt-in, exploratory screen; it must
# never take the rest of the pipeline down with it.
write_skip <- function(reason) {
  warning(reason, call. = FALSE)
  empty_fractions <- data.frame(
    sample_id = character(0), cell_type = character(0),
    fraction = numeric(0), reconstruction_r = numeric(0),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(empty_fractions, file.path(dirs$tables, "cell_fractions.tsv"), na = "NA")
  save_plot_pair(
    empty_plot("Signature-matrix cell-fraction deconvolution", reason),
    file.path(dirs$figures, "cell_fractions"), 7.2, 4.8
  )
  write_json_file(
    list(project_id = cfg$project$id, method = "nnls", samples = 0L,
         shared_genes = 0L, skipped = TRUE, reason = reason),
    file.path(args$outdir, "deconvolution_summary.json")
  )
  quit(save = "no", status = 0)
}

counts <- read_counts_contract(cfg$.counts)
annotation <- read_annotation_contract(cfg$.annotation)
metadata <- readr::read_tsv(cfg$.samples, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]

# Library-size-normalised linear expression (CPM). NNLS deconvolution works in
# linear space; log/VST would distort the additive mixing model.
library_size <- colSums(counts)
if (any(library_size == 0)) write_skip("sample(s) with zero total counts; cannot deconvolve")
cpm <- sweep(counts, 2, library_size, "/") * 1e6

signature <- readr::read_tsv(normalizePath(args$signature, mustWork = TRUE), show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
if (ncol(signature) < 3L) {
  write_skip("signature needs a gene column plus >= 2 cell-type columns")
}
gene_key <- names(signature)[1]
cell_types <- names(signature)[-1]
signature_genes <- as.character(signature[[gene_key]])

# The signature may be keyed by gene symbol or by gene id; collapse CPM to symbol
# space (highest-mean row per symbol) and also keep the id-keyed matrix, then use
# whichever keying overlaps the signature more.
symbol_cpm <- matrix_to_symbols(cpm, annotation)
overlap_symbol <- intersect(rownames(symbol_cpm), signature_genes)
overlap_id <- intersect(rownames(cpm), signature_genes)
if (length(overlap_id) > length(overlap_symbol)) {
  expression <- cpm
  shared <- overlap_id
  key_space <- "gene_id"
} else {
  expression <- symbol_cpm
  shared <- overlap_symbol
  key_space <- "gene_symbol"
}
if (length(shared) < min_genes) {
  write_skip(sprintf(
    "only %d signature genes overlap the expression matrix (%s space); need >= %d",
    length(shared), key_space, min_genes
  ))
}

signature_matrix <- as.matrix(signature[match(shared, signature_genes), cell_types, drop = FALSE])
rownames(signature_matrix) <- shared
mode(signature_matrix) <- "numeric"
expression_shared <- expression[shared, , drop = FALSE]

# Per-sample NNLS: coefficients (>= 0) normalised to fractions; reconstruction
# Pearson r between the fitted mixture S x and the observed sample b flags poor fits.
sample_ids <- colnames(expression_shared)
fractions <- matrix(NA_real_, nrow = length(sample_ids), ncol = length(cell_types), dimnames = list(sample_ids, cell_types))
reconstruction_r <- stats::setNames(rep(NA_real_, length(sample_ids)), sample_ids)
for (sample in sample_ids) {
  observed <- expression_shared[, sample]
  solution <- nnls::nnls(signature_matrix, observed)
  coefficients <- solution$x
  total <- sum(coefficients)
  if (total > 0) fractions[sample, ] <- coefficients / total
  fitted <- as.numeric(signature_matrix %*% coefficients)
  if (stats::sd(fitted) > 0 && stats::sd(observed) > 0) {
    reconstruction_r[sample] <- suppressWarnings(stats::cor(fitted, observed))
  }
}

fractions_long <- as.data.frame(fractions) %>%
  tibble::rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "cell_type", values_to = "fraction") %>%
  left_join(metadata %>% select(sample_id, dplyr::any_of(cfg$figures$group)), by = "sample_id") %>%
  mutate(reconstruction_r = reconstruction_r[sample_id])
readr::write_tsv(fractions_long, file.path(dirs$tables, "cell_fractions.tsv"), na = "NA")

# Order samples by biological group (palette order), cell types in signature
# column order, so the stacked bars group conditions and stack consistently.
group_col <- cfg$figures$group
configured_group_order <- names(unlist(cfg$figures$palette, use.names = TRUE))
if (!is.null(group_col) && group_col %in% names(metadata)) {
  observed_group_order <- configured_group_order[configured_group_order %in% as.character(metadata[[group_col]])]
  sample_order <- metadata$sample_id[order(match(as.character(metadata[[group_col]]), observed_group_order), metadata$sample_id)]
} else {
  sample_order <- sort(metadata$sample_id)
}
plot_table <- fractions_long %>%
  filter(!is.na(fraction)) %>%
  mutate(
    sample_id = factor(sample_id, levels = sample_order),
    cell_type = factor(cell_type, levels = cell_types)
  )
cell_type_colours <- stats::setNames(
  CELL_TYPE_PALETTE[(seq_along(cell_types) - 1L) %% length(CELL_TYPE_PALETTE) + 1L],
  cell_types
)

if (nrow(plot_table) > 0) {
  fractions_plot <- ggplot(plot_table, aes(sample_id, fraction, fill = cell_type)) +
    geom_col(width = 0.82, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = cell_type_colours, name = "Cell type") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = "Signature-matrix cell-fraction deconvolution",
      subtitle = sprintf("NNLS on %d shared signature genes (%s)", length(shared), key_space),
      x = NULL, y = "Estimated fraction"
    ) +
    theme_publication(8.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
} else {
  fractions_plot <- empty_plot("Signature-matrix cell-fraction deconvolution")
}
save_plot_pair(fractions_plot, file.path(dirs$figures, "cell_fractions"), 7.2, 4.8)

write_json_file(
  list(
    project_id = cfg$project$id,
    method = "nnls",
    key_space = key_space,
    shared_genes = length(shared),
    cell_types = cell_types,
    samples = length(sample_ids),
    min_reconstruction_r = if (all(is.na(reconstruction_r))) NA_real_ else min(reconstruction_r, na.rm = TRUE),
    median_reconstruction_r = if (all(is.na(reconstruction_r))) NA_real_ else stats::median(reconstruction_r, na.rm = TRUE)
  ),
  file.path(args$outdir, "deconvolution_summary.json")
)
