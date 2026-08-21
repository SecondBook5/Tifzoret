#!/usr/bin/env Rscript

# Publication-grade DoRothEA radial regulon map (Figure 2, Panel E).
#
# FIGURE layer of the GRN render-seam. Consumes the program-labelled node/edge
# tables emitted by grn.py (contrast_grn) -- where each node's transcriptional
# program was assigned by anchored network diffusion and the target-overlap
# separation P was computed -- and reproduces the bespoke radial panel: measured
# targets on an outer ring, regulators on an inner ring, sectors sized by program
# target count, ggforce bezier edges coloured by mode of regulation and weighted
# by prior likelihood, nodes filled by a diverging DE effect scale, rotated target
# labels and repelled regulator labels.
#
# Faithful, study-agnostic port of prepare_dorothea_grn_radial() /
# make_dorothea_grn_radial_panel() from the reference bespoke figure library.
# All figure semantics that were hardcoded in the original (the program
# levels/labels/colours, the regulator display order, the diverging fill limits,
# the edge palette, sector geometry) are study-specific curation and are read
# from the project config
# (analysis.settings.regulators.grn). With no curation the renderer writes a
# placeholder and defers to the auditable matplotlib views (grn_rectangular,
# grn_radial_legacy) so the stage still completes for any study.

suppressPackageStartupMessages({
  library(ggforce)
  library(ggrepel)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

args <- parse_cli(c("nodes", "edges", "separation", "project-config", "contrasts", "contrast-id", "outdir"))
dirs <- ensure_output_dirs(args$outdir)
stem <- file.path(dirs$figures, "grn_radial")

cfg <- read_project(args[["project-config"]])
grn <- cfg$analysis$settings$regulators$grn

# Editorial text follows the contrast numerator/denominator so the panel is
# study-agnostic.
contrasts <- readr::read_tsv(args$contrasts, show_col_types = FALSE, progress = FALSE)
contrast_row <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
numerator <- if (nrow(contrast_row) && !is.na(contrast_row$numerator[[1]])) contrast_row$numerator[[1]] else "numerator"
denominator <- if (nrow(contrast_row) && !is.na(contrast_row$denominator[[1]])) contrast_row$denominator[[1]] else "denominator"

# Curation is required for the bespoke radial (it supplies the program identity
# the diffusion diffuses toward). Absent it, defer to the matplotlib legacy view.
if (is.null(grn) || length(grn$programs) != 2L) {
  save_plot_pair(
    empty_plot("DoRothEA radial regulon map",
               "No regulators.grn curation configured; see grn_radial_legacy / grn_rectangular"),
    stem, width = 9.4, height = 8.6
  )
  message("grn_radial: no regulators.grn curation; wrote placeholder")
  quit(save = "no", status = 0)
}

# --- curation -> program levels, labels, colours, regulator order -------------
prog_a <- grn$programs[[1]]          # epithelial-analogue (numerator of the ratio)
prog_b <- grn$programs[[2]]          # myogenic-analogue
shared_label <- grn$shared$label %||% "Shared / bridging context"
peripheral_label <- grn$peripheral$label %||% "Peripheral regulator component"
program_levels <- c(peripheral_label, prog_a$label, shared_label, prog_b$label)
regulator_order <- unlist(grn$regulator_order %||% list(), use.names = FALSE)

program_colours <- stats::setNames(
  c(grn$peripheral$colour %||% "#7F8C96", prog_a$colour %||% "#34824B",
    grn$shared$colour %||% "#9A7B43", prog_b$colour %||% "#D51B70"),
  program_levels
)
program_labels <- stats::setNames(
  c(grn$peripheral$short %||% "Peripheral components", prog_a$short %||% prog_a$label,
    grn$shared$short %||% "Shared / bridging", prog_b$short %||% prog_b$label),
  program_levels
)

start_angle <- as.numeric(grn$start_angle %||% 150)
sector_gap <- as.numeric(grn$sector_gap %||% 7)
fill_limits <- unlist(grn$fill_limits %||% list(-10.5, 10.5), use.names = FALSE)
fill_breaks <- unlist(grn$fill_breaks %||% list(-10, -5, 0, 5, 10), use.names = FALSE)
fill_low <- grn$fill_colours$low %||% "#326D9B"
fill_mid <- grn$fill_colours$mid %||% "#F7F4EE"
fill_high <- grn$fill_colours$high %||% "#D9543F"
edge_activating <- grn$edge_colours$Activating %||% "#5C9BB3"
edge_repressing <- grn$edge_colours$Repressing %||% "#9565A1"
lw_range <- unlist(grn$likelihood$range %||% list(0.28, 0.72), use.names = FALSE)
lw_limits <- unlist(grn$likelihood$limits %||% list(0.5, 0.75), use.names = FALSE)
lw_breaks <- unlist(grn$likelihood$breaks %||% list(0.5, 0.75), use.names = FALSE)
repel_seed <- as.integer(grn$seed %||% 814L)

# --- data from grn.py ---------------------------------------------------------
node_program <- readr::read_tsv(args$nodes, show_col_types = FALSE, progress = FALSE) %>%
  dplyr::transmute(
    name = node, node_type,
    program = factor(program, levels = program_levels),
    value = suppressWarnings(as.numeric(value)),
    padj = suppressWarnings(as.numeric(padj))
  )
grn_edges <- readr::read_tsv(args$edges, show_col_types = FALSE, progress = FALSE) %>%
  dplyr::transmute(
    from = source, to = target, regulation,
    likelihood = suppressWarnings(as.numeric(likelihood)),
    mode_of_regulation = suppressWarnings(as.numeric(mode_of_regulation)),
    target_log2fc = suppressWarnings(as.numeric(target_log2fc)),
    target_padj = suppressWarnings(as.numeric(target_padj))
  )
sep_tbl <- readr::read_tsv(args$separation, show_col_types = FALSE, progress = FALSE)
exact_p <- if ("exact_p" %in% names(sep_tbl) && nrow(sep_tbl)) suppressWarnings(as.numeric(sep_tbl$exact_p[[1]])) else NA_real_

target_nodes <- node_program %>% dplyr::filter(node_type == "target")
if (!nrow(target_nodes)) {
  save_plot_pair(empty_plot("DoRothEA radial regulon map", "No displayed regulator targets"), stem, width = 9.4, height = 8.6)
  message("grn_radial: no displayed targets; wrote placeholder")
  quit(save = "no", status = 0)
}
regulator_program <- node_program %>%
  dplyr::filter(node_type == "regulator") %>%
  dplyr::select(regulator = name, regulator_program = program)

# --- target ordering within each program sector (primary regulator, ranked) ---
target_order_key <- grn_edges %>%
  dplyr::select(regulator = from, target = to) %>%
  dplyr::inner_join(target_nodes %>% dplyr::select(target = name, target_program = program), by = "target") %>%
  dplyr::left_join(regulator_program, by = "regulator") %>%
  dplyr::mutate(
    same_program = as.character(regulator_program) == as.character(target_program),
    regulator_rank = match(regulator, regulator_order)
  ) %>%
  dplyr::group_by(target, target_program) %>%
  dplyr::arrange(dplyr::desc(same_program), regulator_rank, regulator, .by_group = TRUE) %>%
  dplyr::slice_head(n = 1L) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(name = target, program = target_program, primary_regulator = regulator, primary_regulator_rank = regulator_rank)

target_key <- target_nodes %>%
  dplyr::select(name, program) %>%
  dplyr::left_join(target_order_key, by = c("name", "program")) %>%
  dplyr::mutate(primary_regulator_rank = dplyr::coalesce(primary_regulator_rank, Inf)) %>%
  dplyr::arrange(program, primary_regulator_rank, primary_regulator, name)

# --- proportional sector geometry (angles sweep clockwise from start_angle) ---
program_counts <- target_key %>%
  dplyr::count(program, name = "target_count", .drop = FALSE) %>%
  dplyr::filter(target_count > 0) %>%
  dplyr::mutate(
    program = factor(as.character(program), levels = program_levels),
    available_angle = 360 - sector_gap * dplyr::n(),
    sector_span = available_angle * target_count / sum(target_count),
    sector_start = NA_real_, sector_end = NA_real_
  ) %>%
  dplyr::arrange(program)
current_angle <- start_angle
for (i in seq_len(nrow(program_counts))) {
  program_counts$sector_start[[i]] <- current_angle
  program_counts$sector_end[[i]] <- current_angle - program_counts$sector_span[[i]]
  current_angle <- program_counts$sector_end[[i]] - sector_gap
}
program_counts <- program_counts %>% dplyr::mutate(sector_mid = (sector_start + sector_end) / 2)

target_key <- target_key %>%
  dplyr::left_join(program_counts, by = "program") %>%
  dplyr::group_by(program) %>%
  dplyr::mutate(
    target_index = dplyr::row_number(),
    angle = sector_start - sector_span * (target_index - 0.5) / target_count
  ) %>%
  dplyr::ungroup()

regulator_angles <- grn_edges %>%
  dplyr::select(regulator = from, target = to) %>%
  dplyr::inner_join(target_key %>% dplyr::select(target = name, target_angle = angle), by = "target") %>%
  dplyr::group_by(regulator) %>%
  dplyr::summarise(angle = stats::median(target_angle), .groups = "drop")
regulator_nodes <- node_program %>%
  dplyr::filter(node_type == "regulator") %>%
  dplyr::left_join(regulator_angles, by = c("name" = "regulator")) %>%
  dplyr::left_join(program_counts %>% dplyr::select(program, sector_mid), by = "program") %>%
  dplyr::mutate(angle = dplyr::coalesce(angle, sector_mid))

radial_nodes <- dplyr::bind_rows(
  target_key %>%
    dplyr::select(name, program, angle) %>%
    dplyr::left_join(node_program %>% dplyr::select(name, node_type, value, padj), by = "name") %>%
    dplyr::mutate(radius = 1.00),
  regulator_nodes %>%
    dplyr::select(name, program, angle, node_type, value, padj) %>%
    dplyr::mutate(radius = 0.58)
) %>%
  dplyr::mutate(
    angle_radians = angle * pi / 180,
    x = radius * cos(angle_radians),
    y = radius * sin(angle_radians),
    angle_normalized = angle %% 360,
    right_side = cos(angle_radians) >= 0,
    label_angle = dplyr::if_else(right_side, angle_normalized, angle_normalized + 180),
    label_hjust = dplyr::if_else(right_side, 0, 1),
    label_x = dplyr::if_else(node_type == "target", 1.075, 0.50) * cos(angle_radians),
    label_y = dplyr::if_else(node_type == "target", 1.075, 0.50) * sin(angle_radians),
    node_size = dplyr::if_else(node_type == "regulator", 5.2, 3.0)
  )

# --- bezier edge control paths (regulator inner ring -> target outer ring) ----
radial_edges <- grn_edges %>%
  dplyr::left_join(radial_nodes %>% dplyr::select(from = name, from_angle = angle, from_radius = radius), by = "from") %>%
  dplyr::left_join(radial_nodes %>% dplyr::select(to = name, to_angle = angle, to_radius = radius), by = "to") %>%
  dplyr::filter(!is.na(from_angle), !is.na(to_angle)) %>%
  dplyr::mutate(
    edge_id = dplyr::row_number(),
    angular_delta = ((to_angle - from_angle + 180) %% 360) - 180,
    target_control_radius = dplyr::if_else(to_radius < 0.8, 0.34, 0.72),
    target_endpoint_radius = dplyr::if_else(to_radius < 0.8, 0.61, 0.975)
  )
bezier_point <- function(edges, k, angle, radius) {
  tibble::tibble(
    edge_id = edges$edge_id, point_order = k,
    x = radius * cos(angle * pi / 180), y = radius * sin(angle * pi / 180),
    regulation = edges$regulation,
    # In-data squish to the width-scale limits: bespoke uses limits c(0.5, 0.75),
    # and A-tier edges (likelihood 1.0) would otherwise be censored to NA (dropped)
    # by the scale. Clamping here maps them to the top of the width band and is
    # portable across ggplot2 versions (scale_linewidth_continuous predates `oob`).
    likelihood = pmin(pmax(edges$likelihood, lw_limits[1]), lw_limits[2])
  )
}
edge_paths <- dplyr::bind_rows(
  bezier_point(radial_edges, 1L, radial_edges$from_angle, radial_edges$from_radius + 0.025),
  bezier_point(radial_edges, 2L, radial_edges$from_angle, 0.34),
  bezier_point(radial_edges, 3L, radial_edges$from_angle + radial_edges$angular_delta, radial_edges$target_control_radius),
  bezier_point(radial_edges, 4L, radial_edges$to_angle, radial_edges$target_endpoint_radius)
) %>%
  dplyr::arrange(edge_id, point_order)

# --- sector arcs + labels -----------------------------------------------------
sector_paths <- do.call(rbind, lapply(seq_len(nrow(program_counts)), function(i) {
  program_name <- as.character(program_counts$program[[i]])
  angle <- seq(program_counts$sector_start[[i]], program_counts$sector_end[[i]], length.out = 100)
  data.frame(
    program = program_name, angle = angle,
    x = 1.035 * cos(angle * pi / 180), y = 1.035 * sin(angle * pi / 180),
    colour = unname(program_colours[program_name]), stringsAsFactors = FALSE
  )
}))
sector_labels <- program_counts %>%
  dplyr::mutate(
    program = as.character(program),
    angle_radians = sector_mid * pi / 180,
    x = 0.80 * cos(angle_radians), y = 0.80 * sin(angle_radians),
    label = unname(program_labels[program]), colour = unname(program_colours[program])
  )

# --- render (make_dorothea_grn_radial_panel) ----------------------------------
target_ring <- radial_nodes %>% dplyr::filter(node_type == "target")
regulator_ring <- radial_nodes %>% dplyr::filter(node_type == "regulator")
fill_title <- paste0(numerator, " − ", denominator, "\n(target log2FC / regulator activity)")
subtitle <- paste0(
  "All ", nrow(radial_nodes), " nodes and ", nrow(radial_edges),
  " signed prior edges; regulators occupy the inner ring and all measured targets the outer ring",
  if (!is.na(exact_p)) paste0(
    "\nSector ordering follows network-diffusion programs (target-overlap permutation P = ",
    formatC(exact_p, format = "f", digits = 4), ")"
  ) else "\nSector ordering follows anchored network-diffusion programs"
)

p <- ggplot()
for (program_name in unique(sector_paths$program)) {
  program_path <- sector_paths[sector_paths$program == program_name, , drop = FALSE]
  p <- p + geom_path(
    data = program_path, aes(x, y, group = program),
    colour = unique(program_path$colour), linewidth = 4.2, alpha = 0.22, lineend = "round"
  )
}
panel <- p +
  ggforce::geom_bezier(
    data = edge_paths, aes(x, y, group = edge_id, colour = regulation, linewidth = likelihood),
    alpha = 0.48, lineend = "round", n = 45, arrow = arrow(length = unit(0.95, "mm"), type = "closed")
  ) +
  geom_point(
    data = target_ring, aes(x, y, fill = value, shape = node_type, size = node_size),
    colour = "#526978", stroke = 0.42
  ) +
  geom_point(
    data = regulator_ring, aes(x, y, fill = value, shape = node_type, size = node_size),
    colour = "#304C60", stroke = 0.75
  ) +
  geom_text(
    data = target_ring, aes(label_x, label_y, label = name, angle = label_angle, hjust = label_hjust),
    colour = NAVY, size = 2.15
  ) +
  geom_label_repel(
    data = regulator_ring, aes(x, y, label = name), seed = repel_seed,
    colour = NAVY, fill = scales::alpha("white", 0.90), label.size = 0, fontface = "bold",
    size = 2.45, box.padding = 0.30, point.padding = 0.18, force = 1.6, force_pull = 0.4,
    min.segment.length = 0, segment.colour = "#9FAAB3", segment.size = 0.22, max.overlaps = Inf, max.time = 5
  ) +
  geom_label(
    data = sector_labels, aes(x, y, label = label), colour = sector_labels$colour,
    fill = scales::alpha("white", 0.90), fontface = "bold", size = 2.6, linewidth = 0,
    label.padding = unit(1.0, "mm")
  ) +
  scale_fill_gradient2(
    low = fill_low, mid = fill_mid, high = fill_high, midpoint = 0,
    limits = fill_limits, breaks = fill_breaks, oob = scales::squish, name = fill_title
  ) +
  scale_shape_manual(
    values = c(regulator = 22, target = 21),
    labels = c(regulator = "Regulator", target = "Measured target"), name = "Node type"
  ) +
  scale_size_identity(guide = "none") +
  scale_colour_manual(values = c(Activating = edge_activating, Repressing = edge_repressing), name = "DoRothEA edge") +
  # Likelihood is squished into lw_limits in-data (see bezier_point) so A-tier
  # (likelihood 1.0) edges render at the top of the width band, not dropped.
  scale_linewidth_continuous(range = lw_range, limits = lw_limits, breaks = lw_breaks, name = "Prior\nlikelihood") +
  guides(
    fill = guide_colourbar(order = 1, title.position = "top", direction = "horizontal", barwidth = unit(36, "mm"), barheight = unit(3.2, "mm")),
    shape = guide_legend(order = 2, override.aes = list(size = c(4.0, 3.0), fill = "white", colour = "#405A6C")),
    colour = guide_legend(order = 3),
    linewidth = guide_legend(order = 4)
  ) +
  labs(title = "DoRothEA radial regulon map", subtitle = subtitle) +
  coord_equal(xlim = c(-1.43, 1.43), ylim = c(-1.43, 1.43), clip = "off", expand = FALSE) +
  theme_void(base_size = 8.2) +
  theme(
    text = element_text(colour = NAVY),
    plot.title = element_text(face = "bold", size = 12.3, margin = margin(b = 2)),
    plot.subtitle = element_text(colour = MID_GREY, size = 8.0, lineheight = 1.08, margin = margin(b = 5)),
    legend.position = "bottom", legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 7.2), legend.text = element_text(size = 6.7),
    legend.key.width = unit(5.0, "mm"), plot.margin = margin(8, 12, 6, 10)
  )

save_plot_pair(panel, stem, width = 9.4, height = 8.6)
message(sprintf("grn_radial: wrote %s.{pdf,png} (%d nodes, %d edges, %d sectors)",
                stem, nrow(radial_nodes), nrow(radial_edges), nrow(program_counts)))
