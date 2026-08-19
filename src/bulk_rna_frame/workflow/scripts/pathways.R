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

build_curve <- function(curve_table, pathway_row, numerator = "numerator", denominator = "denominator") {
  # Presentation parity with the finalized publication panel (Figure 1G): a
  # stacked GSEA multitrack -- running-enrichment ribbon, leading-edge hit rug,
  # ranked-metric colour strip, and signed ranking-statistic profile -- composed
  # with patchwork. Colours, sizes, annotations and layout replicate the
  # reference styling library; every plotted value is read from the already
  # computed curve table, and no statistic is recomputed here.
  cape_ink <- "#B55252"      # numerator-enriched ink (reference CAPE_INK)
  control_ink <- "#39799C"   # denominator-enriched ink (reference CONTROL_INK)

  nes <- pathway_row$NES[[1]]
  n_genes <- nrow(curve_table)
  line_col <- if (nes >= 0) cape_ink else control_ink

  peak_index <- if (nes >= 0) which.max(curve_table$running_es) else which.min(curve_table$running_es)
  peak_rank <- curve_table$rank[[peak_index]]
  es <- curve_table$running_es[[peak_index]]
  zero_cross <- which.min(abs(curve_table$metric))

  mapped_genes <- sum(curve_table$hit)
  leading_count <- sum(curve_table$leading_edge)
  collection_label <- if (identical(pathway_row$gene_set_source[[1]], "configured_gene_program")) {
    "Custom manuscript program"
  } else {
    "MSigDB gene set"
  }
  direction_label <- paste0("Enriched toward ", cond_display(if (nes >= 0) numerator else denominator))
  title <- stringr::str_wrap(pathway_row$pathway_label[[1]], width = 46)
  subtitle <- sprintf(
    "%s · %s  |  NES %.2f  |  FDR %s  |  %d mapped genes  |  %d leading-edge",
    collection_label, direction_label, nes,
    formatC(pathway_row$padj[[1]], format = "e", digits = 1),
    mapped_genes, leading_count
  )

  hits <- curve_table %>% filter(hit)
  gene_labs <- hits %>% filter(leading_edge) %>% slice_max(abs(metric), n = 4)

  es_plot <- ggplot(curve_table, aes(rank, running_es)) +
    annotate("rect",
             xmin = if (nes >= 0) 1 else peak_rank,
             xmax = if (nes >= 0) peak_rank else n_genes,
             ymin = -Inf, ymax = Inf, fill = alpha(line_col, 0.08)) +
    geom_hline(yintercept = 0, colour = "#9DA7AF", linewidth = 0.35) +
    geom_vline(xintercept = peak_rank, colour = alpha(line_col, 0.45), linewidth = 0.35, linetype = 3) +
    geom_ribbon(aes(ymin = 0, ymax = running_es), fill = alpha(line_col, 0.17)) +
    geom_line(colour = line_col, linewidth = 0.8) +
    geom_point(data = curve_table[peak_index, , drop = FALSE], colour = line_col,
               fill = "white", shape = 21, size = 2.2, stroke = 0.65) +
    annotate("text", x = peak_rank, y = es, label = sprintf(" ES %.2f", es),
             hjust = if (es > 0) 0 else 1, vjust = if (es > 0) -0.8 else 1.4,
             size = 2.35, fontface = "bold", colour = line_col) +
    labs(title = title, subtitle = subtitle, y = "Running enrichment score", x = NULL) +
    scale_x_continuous(limits = c(1, n_genes), expand = c(0, 0)) +
    theme_publication(8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(5, 5, 0, 5), panel.grid.major.x = element_blank())
  hits_plot <- ggplot(hits, aes(rank, 0, colour = leading_edge)) +
    geom_segment(aes(xend = rank, y = -0.42, yend = 0.42), linewidth = 0.38) +
    ggrepel::geom_text_repel(
      data = gene_labs, aes(label = gene_symbol), y = 0.48, size = 1.9,
      direction = "x", angle = 45, hjust = 0, vjust = 0, min.segment.length = 0,
      segment.colour = "#AAB3BA", max.overlaps = Inf
    ) +
    scale_colour_manual(values = c(`FALSE` = "#111111", `TRUE` = line_col)) +
    scale_x_continuous(limits = c(1, n_genes), expand = c(0, 0)) +
    coord_cartesian(ylim = c(-0.48, 1.15), clip = "off") +
    labs(y = "Gene hits", x = NULL) +
    theme_publication(7.4) +
    theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), plot.margin = margin(0, 5, 0, 5))
  strip_plot <- ggplot(curve_table, aes(rank, 1, fill = metric)) +
    geom_raster() +
    scale_fill_gradient2(low = control_ink, mid = "#F7F4EE", high = cape_ink,
                         midpoint = 0, limits = range(curve_table$metric), guide = "none") +
    geom_vline(xintercept = zero_cross, colour = "#5E6871", linewidth = 0.28, linetype = 3) +
    annotate("text", x = n_genes * 0.02, y = 1, label = paste0(cond_display(numerator), "-correlated"),
             hjust = 0, size = 2.05, colour = "white", fontface = "bold") +
    annotate("text", x = n_genes * 0.98, y = 1, label = paste0(cond_display(denominator), "-correlated"),
             hjust = 1, size = 2.05, colour = "white", fontface = "bold") +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(1, n_genes), expand = FALSE) +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 5, 0, 5))
  metric_plot <- ggplot(curve_table, aes(rank, metric)) +
    geom_hline(yintercept = 0, colour = "#9DA7AF", linewidth = 0.3) +
    geom_vline(xintercept = zero_cross, colour = "#6F7B85", linewidth = 0.35, linetype = 3) +
    geom_ribbon(aes(ymin = 0, ymax = pmax(metric, 0)), fill = alpha(cape_ink, 0.5)) +
    geom_ribbon(aes(ymin = pmin(metric, 0), ymax = 0), fill = alpha(control_ink, 0.5)) +
    geom_line(colour = "#65737E", linewidth = 0.25) +
    annotate("text", x = zero_cross, y = 0,
             label = paste0("Zero cross: rank ", comma(zero_cross)),
             vjust = -0.6, hjust = 0.5, size = 2, colour = MID_GREY) +
    scale_x_continuous(limits = c(1, n_genes), expand = c(0, 0), labels = comma) +
    labs(x = "Rank in ordered gene list", y = "Wald statistic") +
    theme_publication(7.4) +
    theme(panel.grid.major.x = element_blank(), plot.margin = margin(0, 5, 4, 5))
  es_plot / hits_plot / strip_plot / metric_plot + patchwork::plot_layout(heights = c(3.2, 0.9, 0.27, 1.25))
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
  program_labels <- panel_cfg$gsea_program_labels
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
      # Prefer the study's editorial display label (gsea_program_labels, shared
      # with Panels H/I) so the curve title reads "WNT/PCP valve" rather than the
      # de-underscored fallback "Wnt pcp valve"; keep clean_term for the absent case.
      program_label <- if (!is.null(program_labels[[panel_id]])) program_labels[[panel_id]] else clean_term(panel_id)
      gene_sets[[configured_id]] <- configured_genes
      gene_set_labels[[configured_id]] <- program_label
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

# GSVA panel data choice (Figure 1F). Method and the displayed set list are
# curation-driven so a study can reproduce an exact published panel: CAPE uses
# Hänzelmann GSVA over an explicit, effect-ordered Hallmark list; the generic
# default stays ssGSEA over every measured set. Restricting the set list here
# does not change other sets' scores (each set is scored independently), and the
# audit table below records exactly what was scored.
gsva_method <- if (is.null(cfg$figures$pathways$gsva_method)) "ssgsea" else cfg$figures$pathways$gsva_method
gsva_set_ids <- cfg$figures$pathways$gsva_sets
if (!is.null(gsva_set_ids)) {
  gsva_input_sets <- gene_sets[intersect(as.character(gsva_set_ids), names(gene_sets))]
} else {
  gsva_input_sets <- gene_sets[setdiff(names(gene_sets), configured_curve_ids)]
}
if (!length(gsva_input_sets)) stop("No gene sets remain for the GSVA panel after curation filtering.", call. = FALSE)
gsva_param <- if (identical(gsva_method, "gsva")) {
  GSVA::gsvaParam(exprData = expression, geneSets = gsva_input_sets, minSize = effective_min_size, maxSize = cfg$gene_sets$max_size, kcdf = "Gaussian")
} else {
  GSVA::ssgseaParam(exprData = expression, geneSets = gsva_input_sets, minSize = effective_min_size, maxSize = cfg$gene_sets$max_size, normalize = TRUE)
}
gsva_matrix <- GSVA::gsva(gsva_param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
gsva_score_label <- if (identical(gsva_method, "gsva")) "GSVA" else "ssGSEA"
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

if (!is.null(gsva_set_ids)) {
  # Curated panel: display exactly the configured sets, in the exact order they
  # are listed (that list is the authoritative published row order -- see the
  # row-ordering block below).
  selected_pathways <- intersect(as.character(gsva_set_ids), rownames(gsva_matrix))
} else {
  selected_pathways <- gsva_diff %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    slice_head(n = cfg$figures$pathways$top_gsva_terms) %>%
    pull(pathway)
}
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
gsva_cluster <- if (is.null(cfg$figures$pathways$gsva_cluster)) TRUE else isTRUE(cfg$figures$pathways$gsva_cluster)
gsva_order <- if (is.null(cfg$figures$pathways$gsva_order)) "differential" else cfg$figures$pathways$gsva_order
if (!is.null(gsva_set_ids)) {
  # Curated panel: the configured gsva_sets list IS the authoritative row order.
  # It is pre-ordered by the published numerator-minus-denominator effect (e.g.
  # the CAPE Panel F order), so honour it verbatim rather than re-deriving from
  # the engine's own GSVA scores -- those differ slightly from the bespoke
  # pipeline and flip near-tied adjacent rows, which would not reproduce the
  # paper exactly. Columns stay in canonical denominator-then-numerator order.
  row_order <- selected_pathways
  order_condition <- metadata[colnames(gsva_z), display_group_col]
  column_order <- colnames(gsva_z)[order(factor(order_condition, levels = display_breaks), colnames(gsva_z))]
} else if (gsva_cluster) {
  row_order <- rownames(gsva_z)[stats::hclust(stats::dist(gsva_z), method = "complete")$order]
  column_order <- colnames(gsva_z)[stats::hclust(stats::as.dist(1 - stats::cor(gsva_z)), method = "average")$order]
} else {
  # Unclustered auto-selected panel: rows by numerator-minus-denominator effect
  # (ascending); columns in canonical denominator-then-numerator sample order.
  if (identical(resolved$type, "pairwise") && identical(gsva_order, "effect")) {
    display_condition <- metadata[colnames(gsva_display), factor_name]
    num_cols <- colnames(gsva_display)[display_condition == numerator]
    den_cols <- colnames(gsva_display)[display_condition == denominator]
    effect <- rowMeans(gsva_display[, num_cols, drop = FALSE]) - rowMeans(gsva_display[, den_cols, drop = FALSE])
    row_order <- names(sort(effect))
  } else {
    row_order <- rownames(gsva_z)
  }
  order_condition <- metadata[colnames(gsva_z), display_group_col]
  column_order <- colnames(gsva_z)[order(factor(order_condition, levels = display_breaks), colnames(gsva_z))]
}
# --- Panel F: Hallmark GSVA activity heatmap (ComplexHeatmap) -----------------
# Mirror the bespoke make_gsva_heatmap(): an unclustered heatmap with rows in the
# configured (numerator-minus-denominator effect) order, a condition strip as the
# top annotation (its own legend suppressed -- the shared condition key rides on
# the neighbouring panels), left-side row names capped at 43 mm, display column
# labels, the #356D9A/#F8F7F3/#C94F4F ramp anchored at +/-1.5 (row_zscore already
# clamps there), white hairline cell borders, and a Row z-score scale on the
# right. Pathway ids stay the matrix rownames; display labels ride on row_labels.
display_labels <- setNames(make.unique(prettify_gene_set_label(rownames(gsva_z))), rownames(gsva_z))
gsva_ordered <- gsva_z[row_order, column_order, drop = FALSE]
row_labels_display <- unname(display_labels[rownames(gsva_ordered)])
column_labels_display <- sample_display_labels(colnames(gsva_ordered), metadata[colnames(gsva_ordered), display_group_col])

gsva_displayed <- as.data.frame(gsva_ordered, check.names = FALSE) %>%
  tibble::rownames_to_column("pathway_id") %>%
  tidyr::pivot_longer(-pathway_id, names_to = "sample_id", values_to = "value") %>%
  mutate(
    feature = unname(display_labels[pathway_id]),
    condition = metadata[as.character(sample_id), display_group_col],
    contrast_id = args[["contrast-id"]]
  )
readr::write_tsv(gsva_displayed, file.path(dirs$tables, "gsva_heatmap_displayed.tsv"))

# A fully-Hallmark curated panel is the published Figure 1F: title it "Hallmark
# GSVA activity", matching the paper. Any other panel keeps the generic title.
gsva_is_hallmark <- !is.null(gsva_set_ids) && length(selected_pathways) > 0 &&
  all(startsWith(selected_pathways, "HALLMARK_"))
gsva_title <- if (gsva_is_hallmark) "Hallmark GSVA activity" else "Pathway activity heatmap"

gsva_condition_strip <- cond_display(as.character(metadata[colnames(gsva_ordered), display_group_col]))
names(gsva_condition_strip) <- colnames(gsva_ordered)
gsva_condition_colors <- setNames(unname(display_palette), cond_display(names(display_palette)))
gsva_top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Condition = gsva_condition_strip,
  col = list(Condition = gsva_condition_colors),
  show_legend = FALSE,
  simple_anno_size = grid::unit(3.8, "mm"),
  annotation_name_side = "left",
  annotation_name_gp = grid::gpar(fontface = "bold", fontsize = 7.5, col = NAVY)
)
gsva_ht <- ComplexHeatmap::Heatmap(
  gsva_ordered,
  name = "Row\nz-score",
  col = circlize::colorRamp2(c(-1.5, 0, 1.5), c("#356D9A", "#F8F7F3", "#C94F4F")),
  top_annotation = gsva_top_annotation,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_labels = row_labels_display,
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 6.7, col = NAVY),
  row_names_max_width = grid::unit(43, "mm"),
  column_labels = column_labels_display,
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 7.1, col = NAVY),
  column_title = gsva_title,
  column_title_gp = grid::gpar(fontface = "bold", fontsize = 10, col = NAVY),
  border = FALSE,
  rect_gp = grid::gpar(col = "white", lwd = 0.35),
  use_raster = FALSE,
  heatmap_legend_param = list(
    at = c(-1.5, 0, 1.5),
    labels = c("−1.5", "0", "1.5"),
    title_gp = grid::gpar(fontface = "bold", fontsize = 7.5, col = NAVY),
    labels_gp = grid::gpar(fontsize = 7, col = NAVY),
    legend_height = grid::unit(24, "mm")
  )
)
save_complexheatmap_pair(gsva_ht, file.path(dirs$figures, "gsva_heatmap"), 7.2, max(5.2, 0.30 * nrow(gsva_ordered) + 2.2))

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
  curve_plots[[index]] <- build_curve(curve, pathway_row, numerator, denominator)
}
gsea_displayed <- bind_rows(curve_tables)
readr::write_tsv(gsea_displayed, file.path(dirs$tables, "gsea_curves_displayed.tsv"), na = "NA")
if (length(curve_plots)) {
  # No overall title/subtitle: each curve already carries its program title and
  # a self-describing subtitle, and the assembled figure supplies the panel
  # letter. Matches the finalized Panel G, which has no umbrella heading.
  gsea_plot <- patchwork::wrap_plots(curve_plots, ncol = min(2L, length(curve_plots)))
} else {
  gsea_plot <- empty_plot("Preranked gene-set enrichment")
}
save_plot_pair(gsea_plot, file.path(dirs$figures, "gsea_curves"), 13.5, max(5.4, 4.2 * ceiling(length(curve_plots) / 2)))

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
