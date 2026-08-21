#!/usr/bin/env Rscript

# Publication-grade STRING figures for the functional/regulatory network figure.
# Consumes the directional STRING enrichment tables emitted by networks.py and
# renders the three-facet functional-enrichment bubble (Figure 2, Panel A).
#
# Faithful port of make_string_enrichment_bubble()/prepare_string_enrichment_bubble()
# from the reference bespoke figure library. Editorial text is derived from the
# contrast numerator so the renderer is study-agnostic (e.g. the up-regulated
# facet is titled "Upregulated in <numerator>"), while still reproducing the
# original published panel byte-for-byte on the same inputs.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

args <- parse_cli(c("up", "down", "leading-edge", "contrasts", "contrast-id", "outdir"))
dirs <- ensure_output_dirs(args$outdir)

# Number of nonredundant terms shown per facet (bespoke default).
N_PER_GROUP <- 8L
# STRING GO categories retained by the bespoke bubble (BP/MF/CC).
BUBBLE_CATEGORIES <- c("Process", "Function", "Component")

contrasts <- readr::read_tsv(args$contrasts, show_col_types = FALSE, progress = FALSE)
contrast <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
numerator <- if (nrow(contrast) && !is.na(contrast$numerator[[1]])) contrast$numerator[[1]] else "numerator"

down_label <- sprintf("Downregulated\nin %s", numerator)
le_label <- "GSEA leading\nedge"
up_label <- sprintf("Upregulated\nin %s", numerator)

read_enrichment <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  tbl <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  if (!nrow(tbl)) NULL else tbl
}

string_up_enr <- read_enrichment(args$up)
string_down_enr <- read_enrichment(args$down)
string_le_enr <- read_enrichment(args[["leading-edge"]])

prepare_group <- function(x, direction_label, direction_index) {
  if (is.null(x) || !nrow(x)) {
    return(NULL)
  }
  x %>%
    dplyr::filter(category %in% BUBBLE_CATEGORIES) %>%
    # Deterministic ranking: FDR ascending, then favour the better-supported term
    # (more mapped genes) and finally alphabetical, so tied-FDR terms at the
    # slice_head() boundary resolve reproducibly rather than by STRING's raw
    # response order (which differs run-to-run). Matches the bespoke reference panel.
    dplyr::arrange(fdr, dplyr::desc(number_of_genes), description) %>%
    dplyr::distinct(description, .keep_all = TRUE) %>%
    dplyr::slice_head(n = N_PER_GROUP) %>%
    dplyr::transmute(
      direction_label = direction_label,
      direction_index = direction_index,
      category = factor(category, levels = BUBBLE_CATEGORIES),
      term = term,
      description = description,
      mapped_genes = number_of_genes,
      fdr = fdr,
      score = -log10(pmax(fdr, .Machine$double.xmin)),
      preferredNames = preferredNames
    )
}

bubble <- dplyr::bind_rows(
  prepare_group(string_down_enr, down_label, 1L),
  prepare_group(string_le_enr, le_label, 2L),
  prepare_group(string_up_enr, up_label, 3L)
)

stem <- file.path(dirs$figures, "string_enrichment_faceted")

if (is.null(bubble) || !nrow(bubble)) {
  save_plot_pair(
    empty_plot("STRING functional enrichment", "No STRING enrichment terms in the displayed categories"),
    stem, width = 9.2, height = 4.4
  )
  message("string_figures: no enrichment terms to display; wrote placeholder panel")
  quit(save = "no", status = 0)
}

bubble <- bubble %>%
  dplyr::mutate(
    direction_label = factor(direction_label, levels = c(down_label, le_label, up_label)),
    term_label = str_wrap(str_to_sentence(description), width = 27),
    term_key = paste(direction_index, dplyr::row_number(), sep = "__")
  ) %>%
  dplyr::group_by(direction_label) %>%
  # Order rows up the y-axis by significance; break tied scores by mapped-gene
  # count (better-supported term sits higher) then alphabetically, so the
  # vertical order is deterministic instead of inheriting STRING's row order.
  dplyr::arrange(score, mapped_genes, description, .by_group = TRUE) %>%
  dplyr::ungroup()
bubble$term_key <- factor(bubble$term_key, levels = unique(bubble$term_key))

term_labels <- setNames(as.character(bubble$term_label), as.character(bubble$term_key))
score_limit <- ceiling(max(bubble$score) / 10) * 10

panel <- ggplot(bubble, aes(score, term_key)) +
  geom_point(
    aes(size = mapped_genes, fill = score, shape = category),
    colour = "#3E4650", stroke = 0.65
  ) +
  facet_wrap(vars(direction_label), nrow = 1, scales = "free") +
  scale_y_discrete(labels = term_labels, expand = expansion(add = 0.55)) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.10)), breaks = pretty_breaks(n = 4)) +
  scale_fill_viridis_c(
    option = "magma", begin = 0.12, end = 0.94,
    limits = c(0, score_limit), breaks = c(4, 10, 30, 50, 80), oob = squish,
    name = expression(-log[10](FDR)),
    guide = guide_colourbar(order = 1, title.position = "top", barheight = unit(24, "mm"), barwidth = unit(3.6, "mm"))
  ) +
  scale_size_area(max_size = 8.2, limits = c(0, 50), breaks = c(10, 30, 50), name = "Mapped\ngenes",
                  guide = guide_legend(order = 2)) +
  scale_shape_manual(
    values = c(Process = 21, Function = 22, Component = 24), drop = FALSE,
    name = "STRING\ncategory",
    guide = guide_legend(order = 3, override.aes = list(size = 3.2, fill = "white", colour = "#3E4650", stroke = 0.7))
  ) +
  labs(
    title = "STRING functional enrichment",
    subtitle = sprintf("Top nonredundant ontology terms from each %s-responsive seed set", numerator),
    x = expression(-log[10](FDR)), y = NULL
  ) +
  theme_publication(8.6) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = LIGHT_GREY, linewidth = 0.30),
    strip.background = element_rect(fill = "#F1F3F5", colour = LIGHT_GREY, linewidth = 0.35),
    strip.text = element_text(face = "bold", colour = NAVY, size = 8.2, lineheight = 0.92, margin = margin(3, 4, 3, 4)),
    axis.text.y = element_text(size = 6.7, lineheight = 0.90),
    axis.text.x = element_text(size = 7.0),
    axis.title.x = element_text(margin = margin(t = 5)),
    panel.spacing.x = unit(4.0, "mm"),
    legend.key.height = unit(4.4, "mm"),
    plot.margin = margin(5, 6, 5, 7)
  )

save_plot_pair(panel, stem, width = 9.4, height = 4.4)
message(sprintf("string_figures: wrote %s.{pdf,png} (%d terms across %d facets)",
                stem, nrow(bubble), dplyr::n_distinct(bubble$direction_label)))
