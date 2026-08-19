#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

suppressPackageStartupMessages(library(patchwork))

args <- parse_cli(c("project-config", "samples", "annotation", "contrasts", "contrast-id", "vst", "de", "panels", "outdir"))
cfg <- read_project(args[["project-config"]])
dirs <- ensure_output_dirs(args$outdir)
metadata <- readr::read_tsv(args$samples, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
annotation <- read_annotation_contract(args$annotation)
contrasts <- readr::read_tsv(args$contrasts, show_col_types = FALSE, progress = FALSE)
contrast <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
factor_name <- contrast$factor[[1]]; numerator <- contrast$numerator[[1]]; denominator <- contrast$denominator[[1]]
de <- readr::read_tsv(args$de, show_col_types = FALSE, progress = FALSE)
panel_cfg <- yaml::read_yaml(args$panels)
expression <- matrix_to_symbols(SummarizedExperiment::assay(readRDS(args$vst)), annotation)
metadata <- metadata[match(colnames(expression), metadata$sample_id), , drop = FALSE]; rownames(metadata) <- metadata$sample_id
display_samples <- metadata[[factor_name]] %in% c(denominator, numerator)
expression <- expression[, display_samples, drop = FALSE]
symbol_lookup <- setNames(rownames(expression), toupper(rownames(expression)))

OTHER_PROGRAM <- "Other / poorly characterized"
OTHER_PROGRAM_COLOR <- "#969696"
panel_rows <- list(); panel_colors <- c(`Other / poorly characterized` = OTHER_PROGRAM_COLOR)
default_colors <- c("#D97706", "#E9A300", "#0F9D78", "#C76C9E", "#167BB5", "#53A7D8", "#7A5195", "#8C8C8C")
panel_index <- 0L
configured_panels <- configured_gene_panels(panel_cfg)
for (panel_id in names(configured_panels)) {
  panel_index <- panel_index + 1L
  panel <- configured_panels[[panel_id]]
  color <- if (is.null(panel$color)) default_colors[[1 + (panel_index - 1L) %% length(default_colors)]] else panel$color
  for (group in names(panel$groups)) {
    panel_colors[[group]] <- color
    panel_rows[[length(panel_rows) + 1L]] <- data.frame(
      panel = panel_id,
      program = group,
      gene_symbol = unlist(panel$groups[[group]]),
      program_color = color,
      expected_direction = if (is.null(panel$expected_direction)) "not_specified" else panel$expected_direction,
      expected_contrast = if (is.null(panel$contrast)) NA_character_ else panel$contrast
    )
  }
}
panel_definitions <- bind_rows(panel_rows) %>% distinct(gene_symbol, .keep_all = TRUE)
if (!is.null(panel_cfg$program_colors)) {
  configured_program_colors <- unlist(panel_cfg$program_colors, use.names = TRUE)
  panel_colors[names(configured_program_colors)] <- configured_program_colors
}
panel_definitions <- panel_definitions %>%
  mutate(program_color = ifelse(program %in% names(panel_colors), unname(panel_colors[program]), program_color))
configured_program_definitions <- panel_definitions %>%
  mutate(
    configured_gene_symbol = gene_symbol,
    measured_gene_symbol = unname(symbol_lookup[toupper(gene_symbol)]),
    measured = !is.na(measured_gene_symbol)
  )
readr::write_tsv(
  configured_program_definitions,
  file.path(dirs$tables, "program_definitions.tsv"),
  na = "NA"
)
heatmap_definitions <- panel_definitions
if (!is.null(panel_cfg$program_annotations)) {
  annotations <- data.frame(gene_symbol = names(panel_cfg$program_annotations), program = unlist(panel_cfg$program_annotations), stringsAsFactors = FALSE)
  missing_colors <- setdiff(unique(annotations$program), names(panel_colors))
  if (length(missing_colors)) panel_colors[missing_colors] <- scales::hue_pal()(length(missing_colors))
  heatmap_definitions <- bind_rows(
    annotations %>% mutate(panel = "program_annotations", program_color = unname(panel_colors[program])),
    panel_definitions %>% filter(!gene_symbol %in% annotations$gene_symbol)
  )
}
heatmap_definitions <- heatmap_definitions %>%
  mutate(configured_gene_symbol = gene_symbol, gene_symbol = unname(symbol_lookup[toupper(gene_symbol)])) %>%
  filter(!is.na(gene_symbol)) %>% distinct(gene_symbol, .keep_all = TRUE)
panel_definitions <- panel_definitions %>%
  mutate(configured_gene_symbol = gene_symbol, gene_symbol = unname(symbol_lookup[toupper(gene_symbol)])) %>%
  filter(!is.na(gene_symbol)) %>% distinct(gene_symbol, .keep_all = TRUE)
# Optional single-direction selection for the top-DE heatmap (Panel D shows the
# top CAPE-up genes). Defaults to both directions when the knob is absent.
heatmap_direction <- cfg$figures$de$heatmap_direction
if (is.null(heatmap_direction)) heatmap_direction <- "both"
top <- de %>% filter(!is.na(adjusted_p_value), gene_symbol %in% rownames(expression))
if (identical(heatmap_direction, "up")) top <- top %>% filter(log2_fold_change > 0)
if (identical(heatmap_direction, "down")) top <- top %>% filter(log2_fold_change < 0)
top <- top %>% arrange(adjusted_p_value, desc(abs(log2_fold_change))) %>% distinct(gene_symbol, .keep_all = TRUE) %>% slice_head(n = cfg$figures$de$top_heatmap_genes)
assignments <- top %>% select(gene_symbol, gene_id, log2_fold_change, adjusted_p_value) %>%
  left_join(heatmap_definitions, by = "gene_symbol") %>%
  mutate(
    program = ifelse(is.na(program), OTHER_PROGRAM, program),
    program_color = dplyr::coalesce(as.character(program_color), OTHER_PROGRAM_COLOR)
  )
readr::write_tsv(assignments, file.path(dirs$tables, "de_gene_program_assignments.tsv"), na = "NA")
z <- row_zscore(expression[assignments$gene_symbol, , drop = FALSE], cfg$figures$de$z_limit)
column_order <- colnames(z)[stats::hclust(stats::as.dist(1 - stats::cor(z)), method = "average")$order]
global_order <- rownames(z)[stats::hclust(stats::dist(z), method = "complete")$order]
program_levels <- if (is.null(panel_cfg$program_order)) unique(assignments$program) else {
  c(intersect(panel_cfg$program_order, unique(assignments$program)), setdiff(unique(assignments$program), panel_cfg$program_order))
}
program_order <- unlist(lapply(program_levels, function(program) {
  genes <- assignments$gene_symbol[assignments$program == program]
  if (length(genes) < 2L) genes else genes[stats::hclust(stats::dist(z[genes, , drop = FALSE]), method = "complete")$order]
}))

condition_plot <- data.frame(
  sample_id = factor(column_order, levels = column_order),
  condition = factor(metadata[column_order, factor_name], levels = c(denominator, numerator))
) %>%
  ggplot(aes(sample_id, 1, fill = condition)) + geom_tile() +
  geom_text(aes(label = sample_id), angle = 45, hjust = 0, nudge_y = -0.06, size = 2.5, colour = NAVY) +
  scale_fill_manual(values = condition_palette(cfg, c(denominator, numerator)), breaks = c(denominator, numerator), drop = FALSE) +
  coord_cartesian(clip = "off") + theme_void() + theme(legend.position = "top", plot.margin = margin(2, 25, 25, 0))

heatmap_variant <- function(row_order, stem, title, direct_labels = FALSE) {
  table <- as.data.frame(z[row_order, column_order, drop = FALSE], check.names = FALSE) %>%
    tibble::rownames_to_column("gene_symbol") %>% pivot_longer(-gene_symbol, names_to = "sample_id", values_to = "row_z_score") %>%
    left_join(assignments %>% select(gene_symbol, program, program_color, log2_fold_change, adjusted_p_value), by = "gene_symbol") %>%
    mutate(gene_symbol = factor(gene_symbol, levels = rev(row_order)), sample_id = factor(sample_id, levels = column_order), variant = stem, condition = metadata[as.character(sample_id), factor_name])
  readr::write_tsv(table, file.path(dirs$tables, paste0(stem, "_displayed.tsv")), na = "NA")
  heat <- ggplot(table, aes(sample_id, gene_symbol, fill = row_z_score)) + geom_tile(colour = "white", linewidth = 0.18) +
    scale_fill_gradient2(low = "#4C78A8", mid = "#F7F4EE", high = "#D95F5F", midpoint = 0, limits = c(-cfg$figures$de$z_limit, cfg$figures$de$z_limit), name = "Row z-score") +
    labs(title = title, x = NULL, y = NULL) + theme_publication(7.5) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "right")
  strip_data <- assignments %>% filter(gene_symbol %in% row_order) %>% mutate(gene_symbol = factor(gene_symbol, levels = rev(row_order)))
  if (direct_labels) {
    label_data <- strip_data %>%
      mutate(row_number = match(as.character(gene_symbol), levels(gene_symbol))) %>%
      group_by(program) %>%
      summarize(mid = mean(row_number), color = as.character(program_color[[1]]), .groups = "drop")
    left <- ggplot(label_data, aes(1, mid, label = stringr::str_replace_all(program, "_", " "), colour = program)) + geom_text(hjust = 1, fontface = "bold", size = 3) +
      scale_colour_manual(values = setNames(label_data$color, label_data$program), guide = "none") + scale_y_continuous(limits = c(0.5, length(row_order) + 0.5)) + coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(0, 8, 0, 0))
  } else {
    left <- ggplot(strip_data, aes(1, gene_symbol, fill = program)) + geom_tile() +
      scale_fill_manual(values = panel_colors, guide = guide_legend(ncol = 1, title = "Program")) +
      labs(x = NULL, y = NULL) + theme_void() + theme(legend.position = "right", legend.text = element_text(size = 7), plot.margin = margin(0, 4, 0, 0))
  }
  body <- left | heat
  body <- body + plot_layout(widths = if (direct_labels) c(0.24, 1) else c(0.06, 1))
  top <- plot_spacer() | condition_plot
  top <- top + plot_layout(widths = if (direct_labels) c(0.24, 1) else c(0.06, 1))
  combined <- top / body + plot_layout(heights = c(0.10, 1))
  save_plot_pair(combined, file.path(dirs$figures, stem), if (direct_labels) 9.2 else 10.2, max(7.0, 0.19 * length(row_order) + 2.5))
}

heatmap_variant(global_order, "de_heatmap_global", "Top DE genes with global hierarchical clustering")
heatmap_variant(program_order, "de_heatmap_program_grouped", "Top DE genes grouped by biological program")
heatmap_variant(program_order, "de_heatmap_compact", "Top DE genes grouped by biological program", TRUE)

configured_genes <- unique(panel_definitions$gene_symbol)
configured_genes <- configured_genes[configured_genes %in% rownames(expression)]
program_data <- panel_definitions %>% filter(gene_symbol %in% configured_genes) %>% distinct(gene_symbol, .keep_all = TRUE)
program_order_genes <- unlist(lapply(unique(program_data$program), function(program) program_data$gene_symbol[program_data$program == program]))
if (!length(program_order_genes)) {
  stop(
    "No configured biological-program genes matched measured gene symbols. Review tables/program_definitions.tsv and hypothesis_panels.yaml.",
    call. = FALSE
  )
}
# ---------------------------------------------------------------------------
# Panel H: integrated program heatmap + differential-effect forest.
# ---------------------------------------------------------------------------
# The curated gsea_programs are drawn as facet blocks (a gene may appear in more
# than one program -- e.g. Mylk in both contractile and calcium), rows keep
# their configured order, samples keep their natural denominator -> numerator
# order (no clustering), and each gene carries its shrunken-log2FC effect + CI.
gsea_program_ids <- panel_cfg$gsea_programs
if (is.null(gsea_program_ids)) gsea_program_ids <- names(configured_panels)
program_label_overrides <- panel_cfg$gsea_program_labels
program_color_overrides <- panel_cfg$gsea_program_colors
program_default_palette <- c("#D56A24", "#B66A9C", "#2E9C76", "#397EAF", "#B5892E", "#6F7D87")

cond_display <- function(value) {
  ifelse(!is.na(value) & nchar(value) > 1 & value == toupper(value), value, stringr::str_to_title(value))
}

program_rows <- list()
program_h_colors <- character(0)
program_h_levels <- character(0)
program_idx <- 0L
for (program_id in gsea_program_ids) {
  panel <- configured_panels[[program_id]]
  if (is.null(panel)) next
  genes <- unique(unlist(panel$groups, use.names = FALSE))
  measured <- unname(symbol_lookup[toupper(genes)])
  measured <- measured[!is.na(measured) & measured %in% rownames(expression)]
  if (!length(measured)) next
  program_idx <- program_idx + 1L
  label <- if (!is.null(program_label_overrides[[program_id]])) program_label_overrides[[program_id]] else clean_term(program_id)
  color <- if (!is.null(program_color_overrides[[program_id]])) program_color_overrides[[program_id]] else program_default_palette[[1 + (program_idx - 1L) %% length(program_default_palette)]]
  program_h_colors[[label]] <- color
  program_h_levels <- c(program_h_levels, label)
  program_rows[[length(program_rows) + 1L]] <- data.frame(
    program = label, gene_index = seq_along(measured), gene_symbol = measured, stringsAsFactors = FALSE
  )
}
if (!length(program_rows)) {
  stop("No configured gsea_programs matched measured gene symbols for Panel H.", call. = FALSE)
}
program_h_def <- dplyr::bind_rows(program_rows) %>%
  mutate(program = factor(program, levels = program_h_levels), row_id = paste(program, gene_symbol, sep = "__")) %>%
  arrange(program, gene_index)
program_h_levels_row <- program_h_def$row_id                                     # top -> bottom
program_h_row_labels <- setNames(as.character(program_h_def$gene_symbol), program_h_def$row_id)

# Natural sample order: denominator replicates first, then numerator replicates.
h_samples <- colnames(expression)
h_sample_order <- h_samples[order(factor(metadata[h_samples, factor_name], levels = c(denominator, numerator)), h_samples)]
h_sample_labels <- setNames(sample_display_labels(h_sample_order, metadata[h_sample_order, factor_name]), h_sample_order)
h_fill <- condition_palette(cfg, c(denominator, numerator))

# Per-gene row z-scores over the displayed samples (unclamped; the fill scale
# squishes to the configured z limit for display).
h_genes <- unique(program_h_def$gene_symbol)
h_z <- t(scale(t(expression[h_genes, h_sample_order, drop = FALSE])))
h_z[!is.finite(h_z)] <- 0
program_h_long <- as.data.frame(h_z, check.names = FALSE) %>% tibble::rownames_to_column("gene_symbol") %>%
  pivot_longer(-gene_symbol, names_to = "sample_id", values_to = "row_z_score") %>%
  dplyr::inner_join(program_h_def, by = "gene_symbol", relationship = "many-to-many") %>%
  mutate(row_id = factor(row_id, levels = rev(program_h_levels_row)), sample_id = factor(sample_id, levels = h_sample_order), condition = metadata[as.character(sample_id), factor_name])
readr::write_tsv(program_h_long, file.path(dirs$tables, "program_integrated_displayed.tsv"), na = "NA")

# Differential-effect forest data (one row per program-gene).
program_h_forest <- program_h_def %>%
  left_join(de %>% distinct(gene_symbol, .keep_all = TRUE) %>% select(gene_symbol, log2_fold_change, lfc_se, adjusted_p_value), by = "gene_symbol") %>%
  mutate(
    row_id = factor(row_id, levels = rev(program_h_levels_row)),
    ci_low = log2_fold_change - 1.96 * lfc_se, ci_high = log2_fold_change + 1.96 * lfc_se,
    point_class = dplyr::case_when(
      !is.na(adjusted_p_value) & adjusted_p_value < cfg$figures$de$fdr & log2_fold_change >= 0 ~ "numerator",
      !is.na(adjusted_p_value) & adjusted_p_value < cfg$figures$de$fdr & log2_fold_change < 0 ~ "denominator",
      TRUE ~ "ns"
    )
  )
forest_fill <- c(numerator = unname(h_fill[numerator]), denominator = unname(h_fill[denominator]), ns = "white")

# Condition bar over the heatmap columns.
h_condition_bar <- data.frame(
  sample_id = factor(h_sample_order, levels = h_sample_order),
  condition = factor(metadata[h_sample_order, factor_name], levels = c(denominator, numerator)),
  sample_label = unname(h_sample_labels[h_sample_order])
) %>%
  ggplot(aes(sample_id, 1, fill = condition)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sample_label), colour = NAVY, size = 2.45, fontface = "bold") +
  scale_fill_manual(values = h_fill, breaks = c(denominator, numerator), drop = FALSE) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(title = "Sample expression", x = NULL, y = NULL) +
  theme_void(base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 9, colour = NAVY, hjust = 0.5), legend.position = "none", plot.margin = margin(0, 6, 1, 0))

# Left programme colour strip with wrapped block labels.
program_strip_labels <- program_h_def %>% group_by(program) %>% slice(ceiling(dplyr::n() / 2)) %>% ungroup() %>%
  mutate(row_id = factor(row_id, levels = rev(program_h_levels_row)), program_label = stringr::str_replace(as.character(program), " ", "\n"))
program_strip_plot <- program_h_def %>% mutate(row_id = factor(row_id, levels = rev(program_h_levels_row))) %>%
  ggplot(aes(1, row_id, fill = program)) +
  geom_tile(width = 1, height = 1, colour = NA) +
  geom_text(data = program_strip_labels, aes(label = program_label, colour = program), lineheight = 0.92, fontface = "bold", size = 2.65) +
  facet_grid(rows = vars(program), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = setNames(scales::alpha(program_h_colors, 0.11), names(program_h_colors))) +
  scale_colour_manual(values = program_h_colors) +
  scale_y_discrete(expand = expansion(add = 0)) +
  coord_cartesian(xlim = c(0.5, 1.5), expand = FALSE) +
  theme_void(base_family = "sans") +
  theme(legend.position = "none", panel.spacing.y = grid::unit(1.8, "mm"), strip.text.y = element_blank(), strip.background.y = element_blank(), plot.margin = margin(0, 0, 28, 0))

program_heat <- ggplot(program_h_long, aes(sample_id, row_id, fill = row_z_score)) +
  geom_tile(colour = "white", linewidth = 0.28) +
  facet_grid(rows = vars(program), scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = h_sample_labels, expand = c(0, 0)) +
  scale_y_discrete(labels = program_h_row_labels, expand = expansion(add = 0)) +
  scale_fill_gradient2(low = "#356D9A", mid = "#F8F7F3", high = "#C94F4F", midpoint = 0, limits = c(-cfg$figures$de$z_limit, cfg$figures$de$z_limit), oob = scales::squish, name = "Row\nz-score") +
  labs(x = NULL, y = NULL) +
  theme_publication(8.2) +
  theme(panel.grid = element_blank(), panel.spacing.y = grid::unit(1.8, "mm"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(face = "italic", size = 7.0), strip.background.y = element_blank(), strip.text.y = element_blank(), legend.position = "bottom", legend.direction = "horizontal") +
  guides(fill = guide_colourbar(title.position = "left", title.hjust = 0.5, barheight = grid::unit(3, "mm"), barwidth = grid::unit(32, "mm"), ticks.colour = NAVY, frame.colour = NA))

forest_title <- paste(cond_display(numerator), "vs", cond_display(denominator))
program_forest_header <- ggplot() +
  annotate("text", x = 0, y = 1, label = forest_title, hjust = 0.5, vjust = 0.5, fontface = "bold", size = 3.2, colour = NAVY) +
  xlim(-1, 1) + ylim(0, 2) + theme_void() + theme(plot.margin = margin(0, 5, 0, 0))

program_forest <- ggplot(program_h_forest, aes(y = row_id)) +
  geom_vline(xintercept = 0, colour = "#84919B", linewidth = 0.45) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = row_id), linewidth = 0.7, colour = "#52636F") +
  geom_point(aes(x = log2_fold_change, fill = point_class), shape = 21, size = 2.5, stroke = 0.55, colour = "#33434E") +
  facet_grid(rows = vars(program), scales = "free_y", space = "free_y") +
  scale_y_discrete(expand = expansion(add = 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.04))) +
  scale_fill_manual(values = forest_fill, guide = "none") +
  labs(x = "Shrunken log2 fold-change (95% CI)", y = NULL) +
  theme_publication(8.2) +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = LIGHT_GREY, linewidth = 0.3), panel.spacing.y = grid::unit(1.8, "mm"), axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text.y = element_blank(), strip.background.y = element_blank(), legend.position = "none", plot.margin = margin(0, 5, 28, 0))

program_h_header <- (plot_spacer() | h_condition_bar | program_forest_header) + plot_layout(widths = c(0.35, 1.35, 1))
program_h_body <- (program_strip_plot | program_heat | program_forest) + plot_layout(widths = c(0.35, 1.35, 1))
integrated <- (program_h_header / program_h_body) +
  plot_layout(heights = c(0.05, 1)) +
  plot_annotation(
    title = "Curated gene programs: expression patterns and differential effects",
    subtitle = "Genes are grouped by program; filled effect-size points denote FDR < 0.05",
    theme = theme(plot.title = element_text(face = "bold", size = 13, colour = NAVY), plot.subtitle = element_text(size = 9, colour = MID_GREY))
  )
save_plot_pair(integrated, file.path(dirs$figures, "program_integrated"), 10.4, 12.4)

# ---------------------------------------------------------------------------
# Panel I: consolidated per-gene violins, bounded to the curated gsea_programs.
# ---------------------------------------------------------------------------
# A biological program is the atomic unit of the panel, never the individual
# gene: the gene set and its order are identical to Panel H (program_h_def), so
# the panel can never grow past (curated programs x their measured genes) no
# matter how many genes pass FDR. Gene repeats across programs are kept (e.g.
# Mylk in both contractile and calcium). Significance is the DESeq2 adjusted
# p-value -- the same statistic as the volcano and heatmaps -- not a separately
# re-run Wilcoxon test.
violin_def <- program_h_def %>%
  left_join(
    de %>% distinct(gene_symbol, .keep_all = TRUE) %>% select(gene_symbol, adjusted_p_value),
    by = "gene_symbol"
  ) %>%
  mutate(
    program_color = unname(program_h_colors[as.character(program)]),
    significance = case_when(
      is.na(adjusted_p_value) ~ NA_character_,
      adjusted_p_value < 0.001 ~ "***",
      adjusted_p_value < 0.01 ~ "**",
      adjusted_p_value < 0.05 ~ "*",
      TRUE ~ NA_character_
    )
  )

violin <- as.data.frame(t(expression[unique(violin_def$gene_symbol), h_sample_order, drop = FALSE]), check.names = FALSE) %>%
  tibble::rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "gene_symbol", values_to = "expression") %>%
  inner_join(violin_def, by = "gene_symbol", relationship = "many-to-many") %>%
  mutate(condition = factor(metadata[sample_id, factor_name], levels = c(denominator, numerator)))

y_positions <- violin %>% group_by(row_id) %>%
  summarize(y = max(expression, na.rm = TRUE) + 0.10 * diff(range(expression, na.rm = TRUE)), .groups = "drop")
tests <- violin_def %>% left_join(y_positions, by = "row_id") %>% mutate(x1 = 1, x2 = 2)

# Apply the shared program->gene display order (top-left -> bottom-right).
violin$row_id <- factor(violin$row_id, levels = program_h_levels_row)
tests$row_id <- factor(tests$row_id, levels = program_h_levels_row)
facet_bg <- violin_def %>% distinct(row_id, program_color) %>%
  mutate(row_id = factor(row_id, levels = program_h_levels_row))
facet_labels <- setNames(as.character(violin_def$gene_symbol), as.character(violin_def$row_id))

readr::write_tsv(violin, file.path(dirs$tables, "program_violins_displayed.tsv"), na = "NA")
readr::write_tsv(tests, file.path(dirs$tables, "program_violins_tests.tsv"), na = "NA")

# Brackets are drawn only where the DESeq2 FDR is significant; ns slots carry no
# bracket rather than an empty line with no annotation above it.
sig_tests <- dplyr::filter(tests, !is.na(significance))
violin_plot <- ggplot(violin, aes(condition, expression, fill = condition)) +
  geom_rect(data = facet_bg, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = scales::alpha(facet_bg$program_color, 0.09)) +
  geom_violin(trim = FALSE, alpha = 0.72, colour = "white", linewidth = 0.3) + geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.45, linewidth = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 1.6, fill = "white", colour = NAVY, stroke = 0.5) +
  geom_segment(data = sig_tests, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, colour = NAVY, linewidth = 0.35) +
  geom_segment(data = sig_tests, aes(x = x1, xend = x1, y = y, yend = y - 0.03), inherit.aes = FALSE, colour = NAVY, linewidth = 0.35) +
  geom_segment(data = sig_tests, aes(x = x2, xend = x2, y = y, yend = y - 0.03), inherit.aes = FALSE, colour = NAVY, linewidth = 0.35) +
  geom_text(data = sig_tests, aes(x = 1.5, y = y, label = significance), inherit.aes = FALSE, vjust = -0.25, colour = NAVY, size = 2.8) +
  facet_wrap(~ row_id, scales = "free_y", ncol = 5, labeller = as_labeller(facet_labels)) +
  scale_fill_manual(values = condition_palette(cfg, c(denominator, numerator)), breaks = c(denominator, numerator), drop = FALSE) +
  labs(title = "Curated program genes", subtitle = paste0(forest_title, "; DESeq2 FDR:  * < 0.05   ** < 0.01   *** < 0.001"), x = NULL, y = "Variance-stabilized expression", fill = NULL) +
  theme_publication(7.4) + theme(legend.position = "top", strip.background = element_rect(fill = "#F3F6F8", colour = NA), axis.text.x = element_text(angle = 35, hjust = 1))
n_violin_slots <- nrow(violin_def)
save_plot_pair(violin_plot, file.path(dirs$figures, "program_violins"), 13.0, max(7.2, 2.25 * ceiling(n_violin_slots / 5)))

write_json_file(list(
  contrast_id = args[["contrast-id"]],
  hypothesis_panel_source = normalizePath(args$panels, mustWork = TRUE),
  de_heatmap_variants = c("global", "program_grouped", "compact_direct_labels"),
  registered_constructor_ids = c("de_heatmap", "program_heatmap_effects", "program_violins"),
  configured_program_genes = length(program_order_genes),
  configured_program_genes_requested = nrow(configured_program_definitions),
  configured_program_genes_measured = sum(configured_program_definitions$measured),
  expected_effects = panel_cfg$expected_effects,
  constructor_defaults = panel_cfg$constructor_defaults,
  figures = list(
    de_heatmap_global = list(displayed_data = "tables/de_heatmap_global_displayed.tsv"),
    de_heatmap_program_grouped = list(displayed_data = "tables/de_heatmap_program_grouped_displayed.tsv"),
    de_heatmap_compact = list(displayed_data = "tables/de_heatmap_compact_displayed.tsv"),
    program_integrated = list(displayed_data = "tables/program_integrated_displayed.tsv"),
    program_violins = list(displayed_data = "tables/program_violins_displayed.tsv", tests = "tables/program_violins_tests.tsv")
  )
), file.path(args$outdir, "publication_summary.json"))
