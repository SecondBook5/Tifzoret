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
# Map each term to its ontology domain from the provider and the term-id prefix.
# GO ids carry the domain (GO_BP_/GO_CC_/GO_MF_); the GO_BP_ label is identical
# to the previous BP-only behaviour, so the bespoke panel below is unchanged.
domain_label <- function(provider, term) {
  dplyr::case_when(
    provider == "kegg" ~ "KEGG",
    provider == "reactome" ~ "Reactome",
    provider == "go" & grepl("^GO_BP_", term) ~ "GO Biological Process",
    provider == "go" & grepl("^GO_CC_", term) ~ "GO Cellular Component",
    provider == "go" & grepl("^GO_MF_", term) ~ "GO Molecular Function",
    provider == "go" ~ "GO",
    TRUE ~ provider
  )
}
ontology <- ora %>%
  inner_join(resources, by = c("pathway" = "term")) %>%
  filter(provider %in% c("go", "kegg", "reactome")) %>%
  mutate(
    ontology = domain_label(provider, pathway),
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

# Additive breadth view: when extra GO domains (CC/MF) or Reactome are enabled,
# summarise the top enriched terms per domain and direction in one faceted bar
# chart. Byte-identical for BP-only studies except that the (single-facet) figure
# and table are also emitted; the bespoke panel above is never touched.
per_domain <- min(10L, cfg$figures$pathways$top_ora_terms)
domain_top <- ontology %>%
  filter(count > 0, !is.na(adjusted_p_value)) %>%
  group_by(ontology, direction) %>%
  slice_min(adjusted_p_value, n = per_domain, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    score = negative_log10_adjusted_p,
    direction_display = ifelse(direction == "up_in_numerator", paste0("Up in ", numerator), paste0("Down in ", numerator)),
    term_key = paste(ontology, direction, pathway, sep = "::"),
    term_display = str_wrap(as.character(term_label), width = 34)
  )
domain_top <- domain_top %>% arrange(ontology, direction, score) %>% mutate(term_key = factor(term_key, levels = unique(term_key)))
readr::write_tsv(domain_top, file.path(dirs$tables, "ontology_domain_displayed.tsv"))

direction_colours <- c("up" = "#C0392B", "down" = "#2C6FBB")
domain_fill <- stats::setNames(
  ifelse(grepl("^Up", sort(unique(domain_top$direction_display))), direction_colours[["up"]], direction_colours[["down"]]),
  sort(unique(domain_top$direction_display))
)
if (nrow(domain_top)) {
  domains_plot <- ggplot(domain_top, aes(score, term_key, fill = direction_display)) +
    geom_col(width = 0.72, colour = "#40354A", linewidth = 0.2) +
    facet_grid(rows = vars(ontology), scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = setNames(domain_top$term_display, domain_top$term_key)) +
    scale_fill_manual(values = domain_fill, name = NULL) +
    labs(
      title = "Enrichment across ontology domains",
      subtitle = sprintf("Top %d over-represented terms per domain and direction", per_domain),
      x = expression(-log[10](FDR)), y = NULL
    ) +
    theme_publication(8.4) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "#F3F4F5", colour = LIGHT_GREY, linewidth = 0.35),
      strip.text.y = element_text(face = "bold", colour = NAVY, size = 7.4, angle = 0),
      axis.text.y = element_text(size = 6.6, lineheight = 0.92),
      legend.position = "bottom"
    )
} else {
  domains_plot <- empty_plot("Enrichment across ontology domains", "No over-represented terms to display")
}
domain_count <- length(unique(domain_top$ontology))
save_plot_pair(domains_plot, file.path(dirs$figures, "ontology_domains"), 7.4, max(4.4, 1.6 * max(domain_count, 1) + 1.6))

write_json_file(list(
  contrast_id = args[["contrast-id"]],
  terms_tested = nrow(ontology),
  terms_displayed = nrow(displayed),
  significant_terms = sum(ontology$adjusted_p_value < cfg$figures$de$fdr, na.rm = TRUE),
  domains = as.list(sort(unique(ontology$ontology))),
  domain_terms_displayed = nrow(domain_top)
), file.path(args$outdir, "ontology_summary.json"))
