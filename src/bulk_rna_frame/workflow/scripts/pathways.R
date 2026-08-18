#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(fgsea)
  library(GSVA)
  library(limma)
  library(patchwork)
})

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  fields <- strsplit(lines, "\t", fixed = TRUE)
  sets <- lapply(fields, function(row) unique(row[-c(1, 2)]))
  names(sets) <- vapply(fields, `[[`, character(1), 1L)
  sets
}

run_ora <- function(gene_sets, selected, universe, direction) {
  selected <- unique(intersect(selected, universe))
  universe <- unique(universe)
  rows <- lapply(names(gene_sets), function(name) {
    members <- unique(intersect(gene_sets[[name]], universe))
    overlap <- intersect(selected, members)
    n_selected <- length(selected)
    n_members <- length(members)
    count <- length(overlap)
    p_value <- if (n_selected == 0L || n_members == 0L) 1 else stats::phyper(
      count - 1L, n_members, length(universe) - n_members, n_selected,
      lower.tail = FALSE
    )
    data.frame(
      pathway = name,
      direction = direction,
      count = count,
      selected_genes = n_selected,
      pathway_genes = n_members,
      universe_genes = length(universe),
      gene_ratio = if (n_selected) count / n_selected else 0,
      fold_enrichment = if (n_selected && n_members) (count / n_selected) / (n_members / length(universe)) else 0,
      p_value = p_value,
      overlap_genes = paste(overlap, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

running_enrichment <- function(ranks, genes) {
  ranks <- sort(ranks, decreasing = TRUE)
  hit <- names(ranks) %in% genes
  n_hits <- sum(hit)
  n_misses <- length(ranks) - n_hits
  hit_weight <- if (n_hits) abs(ranks) / sum(abs(ranks[hit])) else rep(0, length(ranks))
  increments <- ifelse(hit, hit_weight, if (n_misses) -1 / n_misses else 0)
  data.frame(
    rank = seq_along(ranks),
    gene_symbol = names(ranks),
    metric = unname(ranks),
    hit = hit,
    running_es = cumsum(increments),
    stringsAsFactors = FALSE
  )
}

build_curve <- function(curve_table, pathway_row) {
  colour <- if (pathway_row$NES[[1]] >= 0) "#B55252" else "#39799C"
  title <- stringr::str_wrap(pathway_row$pathway_label[[1]], width = 46)
  leading_count <- sum(curve_table$leading_edge)
  subtitle <- sprintf("NES %.2f   FDR %s   leading edge %d genes", pathway_row$NES[[1]], formatC(pathway_row$padj[[1]], format = "g", digits = 2), leading_count)
  peak_index <- if (pathway_row$NES[[1]] >= 0) which.max(curve_table$running_es) else which.min(curve_table$running_es)
  peak_rank <- curve_table$rank[[peak_index]]
  zero_candidates <- which(diff(sign(curve_table$metric)) != 0)
  zero_rank <- if (length(zero_candidates)) curve_table$rank[[zero_candidates[[1]]]] else NA_integer_
  es_plot <- ggplot(curve_table, aes(rank, running_es)) +
    geom_hline(yintercept = 0, colour = "#AEB8BF", linewidth = 0.3) +
    geom_vline(xintercept = peak_rank, colour = colour, linetype = 3, linewidth = 0.3) +
    geom_line(colour = colour, linewidth = 0.9) +
    annotate("point", x = peak_rank, y = curve_table$running_es[[peak_index]], colour = colour, size = 1.8) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "ES") +
    theme_publication(7.5) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
  hits_plot <- ggplot(curve_table %>% filter(hit), aes(rank, 0, xend = rank, yend = 1, colour = leading_edge)) +
    geom_segment(linewidth = 0.30) +
    scale_colour_manual(values = c(`FALSE` = NAVY, `TRUE` = colour), guide = "none") +
    scale_y_continuous(NULL, breaks = NULL) +
    labs(x = NULL) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  strip_plot <- ggplot(curve_table, aes(rank, 1, fill = metric)) +
    geom_tile() +
    scale_fill_gradient2(low = "#39799C", mid = "#F7F4EE", high = "#B55252", midpoint = 0, guide = "none") +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
  metric_plot <- ggplot(curve_table, aes(rank, metric)) +
    geom_hline(yintercept = 0, colour = "#AEB8BF", linewidth = 0.25) +
    geom_area(data = curve_table %>% filter(metric >= 0), fill = "#F4A6A6", alpha = 0.8) +
    geom_area(data = curve_table %>% filter(metric < 0), fill = "#A6CEE3", alpha = 0.8) +
    geom_line(colour = MID_GREY, linewidth = 0.25) +
    {if (is.finite(zero_rank)) geom_vline(xintercept = zero_rank, linetype = 3, colour = "#697783", linewidth = 0.3)} +
    {if (is.finite(zero_rank)) annotate("text", x = zero_rank, y = 0, label = paste0(" zero cross: ", zero_rank), hjust = 0, vjust = -0.5, size = 2.1, colour = MID_GREY)} +
    labs(x = "Rank in ordered dataset", y = "Metric") +
    theme_publication(7.2) +
    theme(panel.grid.major.x = element_blank())
  es_plot / hits_plot / strip_plot / metric_plot + patchwork::plot_layout(heights = c(2.0, 0.35, 0.20, 1.0))
}

args <- parse_cli(c("project-config", "samples", "annotation", "contrasts", "gmt", "resource-table", "contrast-id", "vst", "de", "outdir"))
cfg <- read_project(args[["project-config"]])
cfg$.samples <- normalizePath(args$samples, mustWork = TRUE)
cfg$.annotation <- normalizePath(args$annotation, mustWork = TRUE)
cfg$.contrasts <- normalizePath(args$contrasts, mustWork = TRUE)
cfg$.gmt <- normalizePath(args$gmt, mustWork = TRUE)
dirs <- ensure_output_dirs(args$outdir)
set.seed(cfg$figures$pathways$seed)

de <- readr::read_tsv(args$de, show_col_types = FALSE, progress = FALSE)
metadata <- readr::read_tsv(cfg$.samples, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
contrasts <- readr::read_tsv(cfg$.contrasts, show_col_types = FALSE, progress = FALSE)
contrast <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
if (nrow(contrast) != 1L) stop("Could not resolve contrast: ", args[["contrast-id"]], call. = FALSE)
resolved <- resolve_contrast(contrast, cfg$design$formula)
factor_name <- resolved$factor_name
numerator <- resolved$numerator
denominator <- resolved$denominator

vst <- readRDS(args$vst)
expression <- SummarizedExperiment::assay(vst)
annotation <- read_annotation_contract(cfg$.annotation)
expression <- matrix_to_symbols(expression, annotation)
metadata <- metadata[match(colnames(expression), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id

resource_catalog <- readr::read_tsv(args[["resource-table"]], show_col_types = FALSE, progress = FALSE) %>%
  filter(!is.na(term), term != "") %>%
  mutate(
    description = ifelse(
      is.na(description) | description == "" | grepl("^MSigDB mouse-native via msigdbr", description),
      clean_term(term),
      description
    ),
    .provider_priority = ifelse(provider == "custom", 2L, 1L)
  ) %>%
  arrange(term, .provider_priority) %>%
  distinct(term, .keep_all = TRUE)
gene_set_labels <- setNames(resource_catalog$description, resource_catalog$term)
gene_set_sources <- setNames(resource_catalog$provider, resource_catalog$term)

gene_sets <- read_gmt(cfg$.gmt)
gene_sets <- lapply(gene_sets, intersect, y = rownames(expression))
gene_sets <- gene_sets[lengths(gene_sets) >= cfg$gene_sets$min_size & lengths(gene_sets) <= cfg$gene_sets$max_size]

configured_curve_ids <- character()
if (!is.null(args$panels) && nzchar(args$panels)) {
  panel_cfg <- yaml::read_yaml(args$panels)
  requested_programs <- panel_cfg$gsea_programs
  symbol_lookup <- setNames(rownames(expression), toupper(rownames(expression)))
  configured_panels <- configured_gene_panels(panel_cfg)
  if (length(requested_programs)) {
    for (panel_id in requested_programs) {
      panel <- configured_panels[[panel_id]]
      configured_id <- paste0("CONFIGURED_PROGRAM_", toupper(gsub("[^A-Za-z0-9]+", "_", panel_id)))
      configured_genes <- unique(unname(symbol_lookup[toupper(unlist(panel$groups, use.names = FALSE))]))
      configured_genes <- configured_genes[!is.na(configured_genes)]
      if (length(configured_genes) < 3L) {
        warning("Configured GSEA program '", panel_id, "' has fewer than three measured genes and will be omitted.")
        next
      }
      gene_sets[[configured_id]] <- configured_genes
      gene_set_labels[[configured_id]] <- clean_term(panel_id)
      gene_set_sources[[configured_id]] <- "configured_gene_program"
      configured_curve_ids <- c(configured_curve_ids, configured_id)
    }
  }
}
if (!length(gene_sets)) stop("No gene sets remain after measurement and size filtering.", call. = FALSE)

label_gene_set <- function(ids) {
  labels <- unname(gene_set_labels[ids])
  missing <- is.na(labels) | labels == ""
  labels[missing] <- clean_term(ids[missing])
  labels
}
source_gene_set <- function(ids) {
  sources <- unname(gene_set_sources[ids])
  sources[is.na(sources) | sources == ""] <- "gmt"
  sources
}
effective_min_size <- if (length(configured_curve_ids)) 3L else cfg$gene_sets$min_size

rank_table <- de %>%
  filter(!is.na(statistic), is.finite(statistic), !is.na(gene_symbol), gene_symbol != "") %>%
  arrange(desc(abs(statistic))) %>%
  distinct(gene_symbol, .keep_all = TRUE)
ranks <- rank_table$statistic
names(ranks) <- rank_table$gene_symbol
ranks <- sort(ranks, decreasing = TRUE)

fgsea_table <- suppressWarnings(fgsea::fgseaMultilevel(
  pathways = gene_sets,
  stats = ranks,
  minSize = effective_min_size,
  maxSize = cfg$gene_sets$max_size,
  eps = 0,
  BPPARAM = BiocParallel::SerialParam()
)) %>%
  as.data.frame() %>%
  mutate(
    leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"),
    pathway_label = label_gene_set(pathway),
    gene_set_source = source_gene_set(pathway),
    direction = ifelse(NES >= 0, "up_in_numerator", "down_in_numerator"),
    contrast_id = args[["contrast-id"]]
  ) %>%
  arrange(padj, desc(abs(NES)))
readr::write_tsv(fgsea_table, file.path(dirs$tables, "fgsea.tsv"), na = "NA")

universe <- unique(de$gene_symbol[!is.na(de$gene_symbol) & de$gene_symbol != ""])
up_genes <- de$gene_symbol[de$direction == "up_in_numerator"]
down_genes <- de$gene_symbol[de$direction == "down_in_numerator"]
ora_table <- bind_rows(
  run_ora(gene_sets, up_genes, universe, "up_in_numerator"),
  run_ora(gene_sets, down_genes, universe, "down_in_numerator")
) %>%
  mutate(
    adjusted_p_value = p.adjust(p_value, method = "BH"),
    negative_log10_adjusted_p = -log10(pmax(adjusted_p_value, .Machine$double.xmin)),
    pathway_label = label_gene_set(pathway),
    gene_set_source = source_gene_set(pathway),
    contrast_id = args[["contrast-id"]]
  ) %>%
  arrange(direction, adjusted_p_value, desc(fold_enrichment))
readr::write_tsv(ora_table, file.path(dirs$tables, "ora.tsv"))
ora_displayed <- ora_table %>%
  filter(count > 0) %>%
  group_by(direction) %>%
  slice_min(order_by = adjusted_p_value, n = cfg$figures$pathways$top_ora_terms, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    direction_label = ifelse(
      direction == "up_in_numerator",
      paste0("Up in ", numerator),
      paste0("Down in ", numerator)
    ),
    pathway_label = stringr::str_wrap(pathway_label, width = 44)
  )
readr::write_tsv(ora_displayed, file.path(dirs$tables, "ora_displayed.tsv"))

if (nrow(ora_displayed)) {
  ora_displayed <- ora_displayed %>%
    arrange(direction_label, gene_ratio, adjusted_p_value) %>%
    mutate(pathway_label = factor(pathway_label, levels = unique(pathway_label)))
  ora_plot <- ggplot(ora_displayed, aes(gene_ratio, pathway_label, size = count, colour = negative_log10_adjusted_p)) +
    geom_point(alpha = 0.94) +
    facet_grid(direction_label ~ ., scales = "free_y", space = "free_y") +
    scale_colour_viridis_c(option = "magma", direction = -1, name = expression(-log[10](FDR))) +
    scale_size_continuous(range = c(2.2, 7.5), name = "Gene count") +
    labs(
      title = "Gene-set over-representation",
      subtitle = "Directions are defined relative to the contrast numerator",
      x = "Gene ratio",
      y = NULL
    ) +
    theme_publication(8.3) +
    theme(strip.text.y = element_text(face = "bold", colour = NAVY), panel.grid.major.y = element_blank())
} else {
  ora_plot <- empty_plot("Gene-set over-representation")
}
save_plot_pair(ora_plot, file.path(dirs$figures, "ora_bidirectional"), 8.1, max(4.8, 0.34 * nrow(ora_displayed) + 2.2))

ssgsea_parameter <- GSVA::ssgseaParam(
  exprData = expression,
  geneSets = gene_sets,
  minSize = effective_min_size,
  maxSize = cfg$gene_sets$max_size,
  normalize = TRUE
)
gsva_matrix <- GSVA::gsva(ssgsea_parameter, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
gsva_wide <- as.data.frame(gsva_matrix, check.names = FALSE) %>% tibble::rownames_to_column("pathway")
readr::write_tsv(gsva_wide, file.path(dirs$tables, "gsva_scores.tsv"))

if (identical(resolved$type, "pairwise")) {
  metadata$contrast_group <- stats::relevel(factor(metadata[[factor_name]]), ref = denominator)
  gsva_formula_text <- gsub(paste0("\\b", factor_name, "\\b"), "contrast_group", cfg$design$formula)
  gsva_design <- stats::model.matrix(stats::as.formula(gsva_formula_text), data = metadata)
  coefficient_pattern <- paste0("^contrast_group", make.names(numerator), "$")
  gsva_coefficient <- grep(coefficient_pattern, colnames(gsva_design), value = TRUE)
  if (length(gsva_coefficient) != 1L) {
    stop("Could not resolve GSVA coefficient; available: ", paste(colnames(gsva_design), collapse = ", "), call. = FALSE)
  }
} else {
  # Interaction / named-coefficient contrast: build the per-row design directly
  # and extract the named coefficient (model.matrix renames ":" to "." so match
  # on the make.names-normalized column set, then fit by column index).
  for (field in all.vars(resolved$design_formula)) metadata[[field]] <- factor(metadata[[field]])
  for (relevel_factor in names(resolved$reference_levels)) {
    metadata[[relevel_factor]] <- stats::relevel(factor(metadata[[relevel_factor]]), ref = resolved$reference_levels[[relevel_factor]])
  }
  gsva_design <- stats::model.matrix(resolved$design_formula, data = metadata)
  matched <- which(make.names(colnames(gsva_design)) == make.names(resolved$coefficient_name))
  if (length(matched) != 1L) {
    stop("Could not resolve GSVA coefficient ", resolved$coefficient_name, "; available: ", paste(colnames(gsva_design), collapse = ", "), call. = FALSE)
  }
  gsva_coefficient <- matched
}
gsva_fit <- limma::eBayes(limma::lmFit(gsva_matrix, gsva_design))
gsva_diff <- limma::topTable(gsva_fit, coef = gsva_coefficient, number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("pathway") %>%
  mutate(
    pathway_label = label_gene_set(pathway),
    gene_set_source = source_gene_set(pathway),
    direction = ifelse(logFC >= 0, "up_in_numerator", "down_in_numerator"),
    contrast_id = args[["contrast-id"]],
    numerator = numerator,
    denominator = denominator
  )
readr::write_tsv(gsva_diff, file.path(dirs$tables, "gsva_differential.tsv"), na = "NA")

selected_pathways <- gsva_diff %>%
  arrange(adj.P.Val, desc(abs(logFC))) %>%
  slice_head(n = cfg$figures$pathways$top_gsva_terms) %>%
  pull(pathway)
if (identical(resolved$type, "pairwise")) {
  display_samples <- metadata[[factor_name]] %in% c(denominator, numerator)
  display_group_col <- factor_name
  display_palette <- condition_palette(cfg, c(denominator, numerator))
  display_breaks <- c(denominator, numerator)
} else {
  display_samples <- rep(TRUE, nrow(metadata))
  display_group_col <- cfg$figures$group
  display_breaks <- levels(factor(metadata[[display_group_col]]))
  display_palette <- condition_palette(cfg, display_breaks)
}
gsva_display <- gsva_matrix[selected_pathways, display_samples, drop = FALSE]
gsva_z <- row_zscore(gsva_display, 1.5)
row_order <- rownames(gsva_z)[stats::hclust(stats::dist(gsva_z), method = "complete")$order]
column_order <- colnames(gsva_z)[stats::hclust(stats::as.dist(1 - stats::cor(gsva_z)), method = "average")$order]
display_labels <- setNames(make.unique(label_gene_set(rownames(gsva_z))), rownames(gsva_z))
rownames(gsva_z) <- unname(display_labels[rownames(gsva_z)])
row_order_clean <- unname(display_labels[row_order])
gsva_heatmap <- tile_heatmap(gsva_z, row_order_clean, column_order, legend_title = "Row-scaled\nssGSEA score", base_size = 7.7)
gsva_heatmap$plot <- gsva_heatmap$plot + scale_y_discrete(labels = function(value) stringr::str_wrap(value, width = 43))
gsva_displayed <- gsva_heatmap$table %>%
  mutate(
    pathway_id = names(display_labels)[match(as.character(feature), display_labels)],
    condition = metadata[as.character(sample_id), display_group_col],
    contrast_id = args[["contrast-id"]]
  )
readr::write_tsv(gsva_displayed, file.path(dirs$tables, "gsva_heatmap_displayed.tsv"))
annotation_plot <- data.frame(
  sample_id = factor(column_order, levels = column_order),
  condition = metadata[column_order, display_group_col]
) %>%
  ggplot(aes(sample_id, 1, fill = condition)) +
  geom_tile() +
  scale_fill_manual(values = display_palette, breaks = display_breaks, drop = FALSE) +
  theme_void() +
  theme(legend.position = "top", plot.margin = margin(0, 55, 0, 35))
gsva_heatmap$plot <- gsva_heatmap$plot +
  labs(title = "Pathway activity heatmap", subtitle = "Top differential ssGSEA programs; row-scaled within displayed samples")
combined_gsva <- annotation_plot / gsva_heatmap$plot + patchwork::plot_layout(heights = c(0.07, 1))
save_plot_pair(combined_gsva, file.path(dirs$figures, "gsva_heatmap"), 7.7, max(5.2, 0.30 * nrow(gsva_z) + 2.2))

curve_n <- cfg$figures$pathways$gsea_curves_per_direction
if (length(configured_curve_ids)) {
  selected_gsea <- fgsea_table %>%
    filter(pathway %in% configured_curve_ids) %>%
    mutate(.configured_order = match(pathway, configured_curve_ids)) %>%
    arrange(.configured_order) %>%
    select(-.configured_order)
} else {
  selected_gsea <- bind_rows(
    fgsea_table %>% filter(NES >= 0) %>% slice_head(n = curve_n),
    fgsea_table %>% filter(NES < 0) %>% slice_head(n = curve_n)
  ) %>% distinct(pathway, .keep_all = TRUE)
}
curve_tables <- list()
curve_plots <- list()
for (index in seq_len(nrow(selected_gsea))) {
  pathway_row <- selected_gsea[index, , drop = FALSE]
  curve <- running_enrichment(ranks, gene_sets[[pathway_row$pathway[[1]]]]) %>%
    mutate(
      leading_edge = gene_symbol %in% strsplit(pathway_row$leadingEdge[[1]], ";", fixed = TRUE)[[1]],
      pathway = pathway_row$pathway[[1]],
      pathway_label = pathway_row$pathway_label[[1]],
      gene_set_source = pathway_row$gene_set_source[[1]],
      NES = pathway_row$NES[[1]],
      adjusted_p_value = pathway_row$padj[[1]],
      contrast_id = args[["contrast-id"]]
    )
  curve_tables[[index]] <- curve
  curve_plots[[index]] <- build_curve(curve, pathway_row)
}
gsea_displayed <- bind_rows(curve_tables)
readr::write_tsv(gsea_displayed, file.path(dirs$tables, "gsea_curves_displayed.tsv"), na = "NA")
if (length(curve_plots)) {
  gsea_plot <- patchwork::wrap_plots(curve_plots, ncol = min(2L, length(curve_plots))) +
    patchwork::plot_annotation(
      title = "Preranked gene-set enrichment",
      subtitle = paste0("Positive NES is enriched toward ", numerator, "; negative NES toward ", denominator),
      theme = theme(plot.title = element_text(face = "bold", colour = NAVY, size = 12), plot.subtitle = element_text(colour = MID_GREY, size = 8.5))
    )
} else {
  gsea_plot <- empty_plot("Preranked gene-set enrichment")
}
save_plot_pair(gsea_plot, file.path(dirs$figures, "gsea_curves"), 10.2, max(5.4, 4.2 * ceiling(length(curve_plots) / 2)))

write_json_file(
  list(
    project_id = cfg$project$id,
    contrast_id = args[["contrast-id"]],
    numerator = numerator,
    denominator = denominator,
    gene_sets_tested = length(gene_sets),
    fgsea_significant = sum(fgsea_table$padj < cfg$figures$de$fdr, na.rm = TRUE),
    ora_terms = nrow(ora_table),
    gsva_programs = nrow(gsva_diff),
    gsea_curves = nrow(selected_gsea),
    configured_gsea_programs_requested = length(configured_curve_ids),
    configured_gsea_programs_displayed = sum(selected_gsea$gene_set_source == "configured_gene_program")
  ),
  file.path(args$outdir, "pathways_summary.json")
)
