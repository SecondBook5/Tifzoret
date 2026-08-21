#!/usr/bin/env Rscript

# Batch-corrected ordination views. QC (qc.R) draws the PCA and sample-distance
# heatmap on the raw variance-stabilized expression; when a study carries a known
# technical batch (analysis.batch, a samples.tsv column), this module removes that
# batch effect with limma::removeBatchEffect -- preserving the biological grouping
# (figures.group) via the retained design -- and redraws the PCA (before vs after)
# and the corrected sample-distance heatmap. Display/diagnostic only: the
# corrected matrix is NOT fed back into DE (DESeq2 models the batch term in its
# design instead); this shows whether the nominal batch actually drives
# clustering and how much collapses once it is removed.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(limma)
})

args <- parse_cli(c("project-config", "vst", "samples", "outdir"))
cfg <- read_project(args[["project-config"]])
cfg$.samples <- normalizePath(args$samples, mustWork = TRUE)
dirs <- ensure_output_dirs(args$outdir)

vst <- readRDS(normalizePath(args$vst, mustWork = TRUE))
vst_matrix <- SummarizedExperiment::assay(vst)
metadata <- readr::read_tsv(cfg$.samples, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
metadata <- metadata[match(colnames(vst_matrix), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id

batch_variable <- cfg$analysis$batch
if (is.null(batch_variable) || !nzchar(batch_variable)) {
  stop("batch module requires analysis.batch (a samples.tsv column)", call. = FALSE)
}
group_col <- cfg$figures$group
configured_group_order <- names(unlist(cfg$figures$palette, use.names = TRUE))
observed_group_order <- configured_group_order[configured_group_order %in% as.character(metadata[[group_col]])]
metadata[[group_col]] <- factor(metadata[[group_col]], levels = observed_group_order)
metadata[[batch_variable]] <- factor(metadata[[batch_variable]])
batch_levels <- levels(metadata[[batch_variable]])

# Preserve the biological grouping while removing the batch effect. When batch is
# single-level or confounded with the group (rank-deficient retained design),
# correction is undefined; keep the raw matrix and record why so the manifest's
# warning collector surfaces it rather than emitting a silently-wrong figure.
warnings <- character(0)
corrected_matrix <- vst_matrix
correction_applied <- FALSE
if (length(batch_levels) < 2L) {
  warnings <- c(warnings, sprintf("analysis.batch %s has a single level; no batch correction applied", batch_variable))
} else {
  preserve_design <- stats::model.matrix(~ metadata[[group_col]])
  corrected_matrix <- tryCatch(
    {
      out <- limma::removeBatchEffect(vst_matrix, batch = metadata[[batch_variable]], design = preserve_design)
      correction_applied <- TRUE
      out
    },
    error = function(error) {
      warnings <<- c(warnings, paste0("removeBatchEffect failed (batch likely confounded with ", group_col, "): ", conditionMessage(error)))
      vst_matrix
    }
  )
}
readr::write_tsv(
  as.data.frame(corrected_matrix, check.names = FALSE) %>% tibble::rownames_to_column("gene_id"),
  file.path(dirs$tables, "batch_corrected_expression.tsv")
)

palette <- condition_palette(cfg, observed_group_order)
# Batch identity is encoded by point shape (open shapes so the group fill reads
# through); cap at the number of distinct shapes ggplot draws cleanly.
shape_pool <- c(21, 22, 24, 23, 25)
batch_shapes <- stats::setNames(shape_pool[(seq_along(batch_levels) - 1L) %% length(shape_pool) + 1L], batch_levels)

pca_for <- function(mat, state) {
  pca <- stats::prcomp(t(mat), center = TRUE, scale. = FALSE)
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  table <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  table$sample_id <- rownames(table)
  idx <- match(table$sample_id, metadata$sample_id)
  table[[group_col]] <- metadata[[group_col]][idx]
  table[[batch_variable]] <- metadata[[batch_variable]][idx]
  table$state <- state
  # Proportion of PC1/PC2 variance explained by batch (R^2 of PC ~ batch); a drop
  # from "before" to "after" quantifies how much clustering the batch drove.
  batch_r2 <- function(component) {
    if (length(batch_levels) < 2L) return(NA_real_)
    summary(stats::lm(table[[component]] ~ table[[batch_variable]]))$r.squared
  }
  list(table = table, variance = variance, pc1_batch_r2 = batch_r2("PC1"), pc2_batch_r2 = batch_r2("PC2"))
}
before <- pca_for(vst_matrix, "Before correction")
after <- pca_for(corrected_matrix, "After correction")
pca_table <- bind_rows(before$table, after$table) %>%
  mutate(state = factor(state, levels = c("Before correction", "After correction")))
readr::write_tsv(pca_table, file.path(dirs$tables, "batch_pca_coordinates.tsv"))

pca_plot <- ggplot(pca_table, aes(PC1, PC2)) +
  geom_point(aes(fill = .data[[group_col]], shape = .data[[batch_variable]]), size = 3.2, stroke = 0.7, colour = "#3A4750") +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 2.3, colour = NAVY, max.overlaps = Inf, show.legend = FALSE) +
  facet_wrap(~ state, nrow = 1, scales = "free") +
  scale_fill_manual(values = palette, drop = FALSE, labels = cond_display) +
  scale_shape_manual(values = batch_shapes) +
  labs(
    title = "Batch-corrected principal-component analysis",
    subtitle = sprintf("limma::removeBatchEffect on %s, preserving %s", batch_variable, group_col),
    x = "PC1", y = "PC2", fill = NULL, shape = batch_variable
  ) +
  theme_publication(8.5) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)))
save_plot_pair(pca_plot, file.path(dirs$figures, "batch_pca"), 9.6, 5.2)

# Corrected sample-distance heatmap (Euclidean on the corrected matrix), matching
# the qc.R sample-distance panel's magma ramp and average-linkage ordering.
sample_distance <- as.matrix(stats::dist(t(corrected_matrix), method = "euclidean"))
distance_order <- colnames(sample_distance)[stats::hclust(stats::as.dist(sample_distance), method = "average")$order]
distance_table <- as.data.frame(sample_distance, check.names = FALSE) %>% tibble::rownames_to_column("sample_id")
readr::write_tsv(distance_table, file.path(dirs$tables, "batch_sample_distance.tsv"))
distance_long <- distance_table %>%
  pivot_longer(-sample_id, names_to = "sample_id_y", values_to = "euclidean_distance") %>%
  mutate(sample_id = factor(sample_id, levels = distance_order), sample_id_y = factor(sample_id_y, levels = rev(distance_order)))
distance_plot <- ggplot(distance_long, aes(sample_id, sample_id_y, fill = euclidean_distance)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "magma", name = "Euclidean\ndistance") +
  labs(title = "Batch-corrected sample distance", subtitle = "Euclidean distance on batch-corrected expression", x = NULL, y = NULL) +
  theme_publication(8.1) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks = element_blank())
save_plot_pair(distance_plot, file.path(dirs$figures, "batch_sample_distance"), 6.2, 5.5)

write_json_file(
  list(
    project_id = cfg$project$id,
    batch_variable = batch_variable,
    batches = batch_levels,
    group = group_col,
    method = "limma::removeBatchEffect",
    correction_applied = correction_applied,
    samples = ncol(vst_matrix),
    pc1_batch_r2_before = before$pc1_batch_r2,
    pc2_batch_r2_before = before$pc2_batch_r2,
    pc1_batch_r2_after = after$pc1_batch_r2,
    pc2_batch_r2_after = after$pc2_batch_r2,
    warnings = warnings
  ),
  file.path(args$outdir, "batch_summary.json")
)
