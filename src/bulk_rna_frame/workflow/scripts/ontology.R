#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

args <- parse_cli(c("project-config", "ora", "resource-table", "contrast-id", "outdir"))
cfg <- read_project(args[["project-config"]])
dirs <- ensure_output_dirs(args$outdir)
contrasts <- readr::read_tsv(cfg$.contrasts, show_col_types = FALSE, progress = FALSE)
contrast <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
if (nrow(contrast) != 1L) stop("Could not resolve contrast", call. = FALSE)
numerator <- contrast$numerator[[1]]
ora <- readr::read_tsv(args$ora, show_col_types = FALSE, progress = FALSE)
resources <- readr::read_tsv(args[["resource-table"]], show_col_types = FALSE, progress = FALSE) %>%
  select(term, description, provider) %>% distinct()
ontology <- ora %>%
  inner_join(resources, by = c("pathway" = "term")) %>%
  filter(provider %in% c("go", "kegg")) %>%
  mutate(
    ontology = ifelse(provider == "go", "GO Biological Process", "KEGG"),
    direction_label = ifelse(
      direction == "up_in_numerator",
      paste0("Upregulated in ", numerator),
      paste0("Downregulated in ", numerator)
    ),
    term_label = ifelse(is.na(description) | description == "", clean_term(pathway), description)
  ) %>%
  arrange(ontology, direction, adjusted_p_value, desc(gene_ratio))
readr::write_tsv(ontology, file.path(dirs$tables, "ontology.tsv"))

eligible <- ontology %>% filter(ontology == "GO Biological Process", count > 0)
significant <- eligible %>% filter(adjusted_p_value < cfg$figures$de$fdr)
if (!nrow(significant)) significant <- eligible
displayed <- significant %>%
  group_by(direction) %>%
  slice_min(adjusted_p_value, n = cfg$figures$pathways$top_ora_terms, with_ties = FALSE) %>%
  ungroup()
displayed <- displayed %>% arrange(ontology, direction_label, gene_ratio, adjusted_p_value) %>% mutate(term_label = factor(term_label, levels = unique(term_label)))
readr::write_tsv(displayed, file.path(dirs$tables, "ontology_displayed.tsv"))
if (nrow(displayed)) {
  # Presentation parity with the finalized manuscript panel E (combined
  # directional GO BP ORA bubble plot). This replicates the reference
  # constructors prepare_combined_ora()/make_combined_ora() onto the engine's
  # already-computed columns; no statistics are refit and the displayed table
  # written above is untouched. The panel letter ("E") is drawn at assembly
  # time by assemble.py's _label_overlay, so no plot.tag is baked in here.
  #
  # Column mapping: engine gene_ratio = count / selected_genes (a fraction),
  # so gene_ratio_pct = 100 * gene_ratio matches the reference's
  # 100 * Count / <GeneRatio denominator>; engine negative_log10_adjusted_p is
  # the reference's `score` = -log10(FDR).
  plot_data <- displayed %>%
    mutate(
      score = negative_log10_adjusted_p,
      gene_ratio_pct = 100 * gene_ratio,
      direction_facet = factor(
        ifelse(
          direction == "up_in_numerator",
          paste0("Upregulated\nin ", numerator),
          paste0("Downregulated\nin ", numerator)
        ),
        levels = c(
          paste0("Upregulated\nin ", numerator),
          paste0("Downregulated\nin ", numerator)
        )
      ),
      # Reference wraps the pathway term at width 31 for the y-axis labels.
      term_display = str_wrap(as.character(term_label), width = 31),
      # Disambiguate identical terms appearing in both directions.
      term_key = paste(direction, pathway, sep = "::")
    )
  # Order terms within each direction by ascending score (reference term_levels).
  term_levels <- plot_data %>% arrange(direction_facet, score) %>% pull(term_key)
  plot_data <- plot_data %>% mutate(term_key = factor(term_key, levels = term_levels))
  score_limit <- ceiling(max(plot_data$score))
  ratio_limit <- ceiling(max(plot_data$gene_ratio_pct) / 2) * 2

  plot <- ggplot(plot_data, aes(gene_ratio_pct, term_key)) +
    geom_point(
      aes(size = count, fill = score),
      shape = 21, colour = "#40354A", stroke = 0.72
    ) +
    facet_wrap(
      vars(direction_facet),
      ncol = 1,
      scales = "free_y",
      strip.position = "top"
    ) +
    scale_y_discrete(labels = setNames(plot_data$term_display, plot_data$term_key)) +
    scale_x_continuous(
      limits = c(0, ratio_limit),
      breaks = seq(0, ratio_limit, by = 2),
      labels = scales::label_number(suffix = "%", accuracy = 1),
      expand = expansion(mult = c(0.01, 0.04))
    ) +
    scale_fill_viridis_c(
      option = "magma",
      begin = 0.14,
      end = 0.94,
      limits = c(4, score_limit),
      breaks = c(4, 8, 12, 16),
      name = expression(-log[10](FDR)),
      guide = guide_colourbar(
        title.position = "top",
        barheight = grid::unit(24, "mm"),
        barwidth = grid::unit(3.6, "mm")
      )
    ) +
    scale_size_area(
      max_size = 8.5,
      limits = c(0, 80),
      breaks = c(10, 30, 50, 70),
      name = "Genes"
    ) +
    labs(
      title = "GO biological process enrichment",
      subtitle = sprintf(
        "Comparative over-representation analysis of %s-responsive gene sets",
        numerator
      ),
      x = "Gene ratio",
      y = NULL
    ) +
    # Reuse the engine theme helper (theme_publication) at the reference base
    # size (9), then reconstruct theme_pub()'s specific values plus the
    # make_combined_ora() theme overrides so the typography matches exactly.
    theme_publication(9) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.13), margin = margin(b = 2.5)),
      plot.subtitle = element_text(colour = MID_GREY, size = rel(0.86), margin = margin(b = 5)),
      axis.title = element_text(face = "bold", size = rel(0.92)),
      legend.title = element_text(face = "bold", size = rel(0.82)),
      legend.text = element_text(size = rel(0.78)),
      legend.position = "right",
      legend.box = "vertical",
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = LIGHT_GREY, linewidth = 0.30),
      strip.background = element_rect(fill = "#F3F4F5", colour = LIGHT_GREY, linewidth = 0.35),
      strip.text = element_text(
        face = "bold", colour = NAVY, size = 8.5, lineheight = 0.95,
        margin = margin(3, 3, 3, 3)
      ),
      axis.text.y = element_text(size = 7.1, lineheight = 0.93),
      axis.text.x = element_text(size = 7.4),
      axis.title.x = element_text(margin = margin(t = 6)),
      panel.spacing.y = grid::unit(3.5, "mm"),
      legend.key.height = grid::unit(4.6, "mm"),
      plot.margin = margin(5, 6, 5, 6)
    )
} else {
  plot <- empty_plot("GO biological process enrichment", "No GO BP terms passed the configured criteria")
}
# Fixed panel geometry matches the reference's stacked (ncol = 1) layout,
# saved by save_pub() at 7.25 x 5.7 in.
save_plot_pair(plot, file.path(dirs$figures, "ontology_bidirectional"), 7.25, 5.7)
write_json_file(list(
  contrast_id = args[["contrast-id"]],
  terms_tested = nrow(ontology),
  terms_displayed = nrow(displayed),
  significant_terms = sum(ontology$adjusted_p_value < cfg$figures$de$fdr, na.rm = TRUE)
), file.path(args$outdir, "ontology_summary.json"))
