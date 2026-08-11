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
for (panel_id in names(panel_cfg$gene_panels)) {
  panel_index <- panel_index + 1L
  panel <- panel_cfg$gene_panels[[panel_id]]
  color <- if (is.null(panel$color)) default_colors[[1 + (panel_index - 1L) %% length(default_colors)]] else panel$color
  for (group in names(panel$groups)) {
    panel_colors[[group]] <- color
    panel_rows[[length(panel_rows) + 1L]] <- data.frame(panel = panel_id, program = group, gene_symbol = unlist(panel$groups[[group]]), program_color = color)
  }
}
panel_definitions <- bind_rows(panel_rows) %>% distinct(gene_symbol, .keep_all = TRUE)
if (!is.null(panel_cfg$program_colors)) {
  configured_program_colors <- unlist(panel_cfg$program_colors, use.names = TRUE)
  panel_colors[names(configured_program_colors)] <- configured_program_colors
}
panel_definitions <- panel_definitions %>%
  mutate(program_color = ifelse(program %in% names(panel_colors), unname(panel_colors[program]), program_color))
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
top <- de %>% filter(!is.na(adjusted_p_value), gene_symbol %in% rownames(expression)) %>% arrange(adjusted_p_value, desc(abs(log2_fold_change))) %>% distinct(gene_symbol, .keep_all = TRUE) %>% slice_head(n = cfg$figures$de$top_heatmap_genes)
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
program_z <- row_zscore(expression[program_order_genes, column_order, drop = FALSE], cfg$figures$de$z_limit)
program_long <- as.data.frame(program_z, check.names = FALSE) %>% tibble::rownames_to_column("gene_symbol") %>% pivot_longer(-gene_symbol, names_to = "sample_id", values_to = "row_z_score") %>%
  left_join(program_data, by = "gene_symbol") %>% left_join(de %>% select(gene_symbol, log2_fold_change, lfc_se, adjusted_p_value), by = "gene_symbol") %>%
  mutate(gene_symbol = factor(gene_symbol, levels = rev(program_order_genes)), sample_id = factor(sample_id, levels = column_order), condition = metadata[as.character(sample_id), factor_name])
readr::write_tsv(program_long, file.path(dirs$tables, "program_integrated_displayed.tsv"), na = "NA")
program_heat <- ggplot(program_long, aes(sample_id, gene_symbol, fill = row_z_score)) + geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#4C78A8", mid = "#F7F4EE", high = "#D95F5F", midpoint = 0, limits = c(-cfg$figures$de$z_limit, cfg$figures$de$z_limit), name = "Row z-score") +
  labs(title = "Configured biological programs", x = NULL, y = NULL) + theme_publication(7.2) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks = element_blank())
effects <- program_long %>% select(gene_symbol, program, program_color, log2_fold_change, lfc_se, adjusted_p_value) %>% distinct()
effect_plot <- ggplot(effects, aes(log2_fold_change, gene_symbol, colour = program)) + geom_vline(xintercept = 0, colour = "#8D989F", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = log2_fold_change - 1.96 * lfc_se, xmax = log2_fold_change + 1.96 * lfc_se), height = 0, alpha = 0.6) + geom_point(size = 2) +
  scale_colour_manual(values = panel_colors, guide = guide_legend(ncol = 1)) + labs(title = "DE effect", x = "Shrunken log2FC", y = NULL, colour = "Program") + theme_publication(7.2) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
integrated <- program_heat | effect_plot
integrated <- integrated + plot_layout(widths = c(1, 0.55))
save_plot_pair(integrated, file.path(dirs$figures, "program_integrated"), 12.2, max(7.2, 0.18 * length(program_order_genes) + 2.5))

violin <- as.data.frame(t(expression[program_order_genes, , drop = FALSE]), check.names = FALSE) %>% tibble::rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "gene_symbol", values_to = "expression") %>% left_join(program_data, by = "gene_symbol") %>%
  mutate(condition = metadata[sample_id, factor_name], condition = factor(condition, levels = c(denominator, numerator)), gene_symbol = factor(gene_symbol, levels = program_order_genes))
tests <- violin %>% group_by(gene_symbol, program, program_color) %>% summarize(
  p_value = tryCatch(stats::wilcox.test(expression[condition == numerator], expression[condition == denominator], exact = FALSE)$p.value, error = function(error) NA_real_),
  y = max(expression, na.rm = TRUE) + 0.10 * diff(range(expression, na.rm = TRUE)), .groups = "drop"
) %>% mutate(adjusted_p_value = p.adjust(p_value, method = "BH"), significance = case_when(adjusted_p_value < 0.001 ~ "***", adjusted_p_value < 0.01 ~ "**", adjusted_p_value < 0.05 ~ "*", TRUE ~ "ns"), x1 = 1, x2 = 2)
readr::write_tsv(violin, file.path(dirs$tables, "program_violins_displayed.tsv"))
readr::write_tsv(tests, file.path(dirs$tables, "program_violins_tests.tsv"), na = "NA")
background <- tests %>% select(gene_symbol, program, program_color) %>% distinct()
violin_plot <- ggplot(violin, aes(condition, expression, fill = condition)) +
  geom_rect(data = background, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = scales::alpha(background$program_color, 0.07)) +
  geom_violin(trim = FALSE, alpha = 0.72, colour = "white", linewidth = 0.3) + geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.45, linewidth = 0.3) +
  geom_segment(data = tests, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, colour = NAVY, linewidth = 0.35) +
  geom_segment(data = tests, aes(x = x1, xend = x1, y = y, yend = y - 0.03), inherit.aes = FALSE, colour = NAVY, linewidth = 0.35) +
  geom_segment(data = tests, aes(x = x2, xend = x2, y = y, yend = y - 0.03), inherit.aes = FALSE, colour = NAVY, linewidth = 0.35) +
  geom_text(data = tests, aes(x = 1.5, y = y, label = significance), inherit.aes = FALSE, vjust = -0.25, colour = NAVY, size = 2.8) +
  facet_wrap(~ gene_symbol, scales = "free_y", ncol = 5) + scale_fill_manual(values = condition_palette(cfg, c(denominator, numerator)), breaks = c(denominator, numerator), drop = FALSE) +
  labs(title = "Configured program genes", subtitle = "Wilcoxon tests with BH correction; program colors are shown as light facet shading", x = NULL, y = "Variance-stabilized expression", fill = NULL) +
  theme_publication(7.4) + theme(legend.position = "top", strip.background = element_rect(fill = "#F3F6F8", colour = NA), axis.text.x = element_text(angle = 35, hjust = 1))
save_plot_pair(violin_plot, file.path(dirs$figures, "program_violins"), 13.0, max(7.2, 2.25 * ceiling(length(program_order_genes) / 5)))

write_json_file(list(
  contrast_id = args[["contrast-id"]],
  de_heatmap_variants = c("global", "program_grouped", "compact_direct_labels"),
  configured_program_genes = length(program_order_genes),
  figures = list(
    de_heatmap_global = list(displayed_data = "tables/de_heatmap_global_displayed.tsv"),
    de_heatmap_program_grouped = list(displayed_data = "tables/de_heatmap_program_grouped_displayed.tsv"),
    de_heatmap_compact = list(displayed_data = "tables/de_heatmap_compact_displayed.tsv"),
    program_integrated = list(displayed_data = "tables/program_integrated_displayed.tsv"),
    program_violins = list(displayed_data = "tables/program_violins_displayed.tsv", tests = "tables/program_violins_tests.tsv")
  )
), file.path(args$outdir, "publication_summary.json"))
