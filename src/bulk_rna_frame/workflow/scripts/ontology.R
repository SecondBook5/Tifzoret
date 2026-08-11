#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

args <- parse_cli(c("project-config", "ora", "resource-table", "contrast-id", "outdir"))
cfg <- read_project(args[["project-config"]])
dirs <- ensure_output_dirs(args$outdir)
ora <- readr::read_tsv(args$ora, show_col_types = FALSE, progress = FALSE)
resources <- readr::read_tsv(args[["resource-table"]], show_col_types = FALSE, progress = FALSE) %>%
  select(term, description, provider) %>% distinct()
ontology <- ora %>%
  inner_join(resources, by = c("pathway" = "term")) %>%
  filter(provider %in% c("go", "kegg")) %>%
  mutate(
    ontology = ifelse(provider == "go", "GO Biological Process", "KEGG"),
    direction_label = ifelse(direction == "up_in_numerator", "Upregulated in numerator", "Downregulated in numerator"),
    term_label = ifelse(is.na(description) | description == "", clean_term(pathway), description)
  ) %>%
  arrange(ontology, direction, adjusted_p_value, desc(gene_ratio))
readr::write_tsv(ontology, file.path(dirs$tables, "ontology.tsv"))

displayed <- ontology %>% filter(count > 0, adjusted_p_value < cfg$figures$de$fdr)
if (!nrow(displayed) && nrow(ontology)) {
  displayed <- ontology %>% filter(count > 0) %>% group_by(ontology, direction) %>% slice_min(adjusted_p_value, n = cfg$figures$pathways$top_ora_terms, with_ties = FALSE) %>% ungroup()
}
displayed <- displayed %>% arrange(ontology, direction_label, gene_ratio, adjusted_p_value) %>% mutate(term_label = factor(term_label, levels = unique(term_label)))
readr::write_tsv(displayed, file.path(dirs$tables, "ontology_displayed.tsv"))
if (nrow(displayed)) {
  plot <- ggplot(displayed, aes(gene_ratio, term_label, size = count, colour = negative_log10_adjusted_p)) +
    geom_point(alpha = 0.94) +
    facet_grid(ontology ~ direction_label, scales = "free_y", space = "free_y") +
    scale_colour_viridis_c(option = "magma", direction = -1, name = expression(-log[10](FDR))) +
    scale_size_continuous(range = c(2.2, 8), name = "Gene count") +
    labs(title = "Ontology and pathway over-representation", subtitle = "Every significant GO BP and KEGG term is retained; directions refer to the numerator", x = "Gene ratio", y = NULL) +
    theme_publication(8.1) + theme(panel.grid.major.y = element_blank(), strip.text = element_text(face = "bold"))
} else {
  plot <- empty_plot("Ontology and pathway over-representation", "No GO BP or KEGG terms passed the configured criteria")
}
save_plot_pair(plot, file.path(dirs$figures, "ontology_bidirectional"), 10.2, max(5.4, 0.27 * nrow(displayed) + 2.7))
write_json_file(list(
  contrast_id = args[["contrast-id"]],
  terms_tested = nrow(ontology),
  terms_displayed = nrow(displayed),
  significant_terms = sum(ontology$adjusted_p_value < cfg$figures$de$fdr, na.rm = TRUE)
), file.path(args$outdir, "ontology_summary.json"))

