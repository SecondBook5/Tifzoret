#!/usr/bin/env Rscript

# Publication-grade STRING community network (Figure 2, Panels B/C).
#
# Consumes the directional STRING node/edge tables emitted by networks.py
# (built from the top-N directional seed sets) and reproduces the bespoke
# community-network panel: igraph Louvain communities (seeded) over a
# Fruchterman-Reingold layout, ggforce community hulls, degree-sized points
# filled by shrunken log2FC, and repelled gene labels.
#
# Faithful, study-agnostic port of prepare_string_{up,down}_network() /
# make_string_{up,down}_network_panel() from the reference bespoke figure
# library. Community NAMES, colours, the log2FC gradient, and the point-size
# scale are study-specific
# curation and are read from the project config
# (analysis.settings.networks.community_curation.{up,down}); with no curation
# the renderer falls back to topology-derived community labels, a default
# categorical palette, and a data-ranged diverging gradient so it produces a
# sensible network for any study.

suppressPackageStartupMessages({
  library(igraph)
  library(ggforce)
  library(ggrepel)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

args <- parse_cli(c("nodes", "edges", "project-config", "contrasts", "contrast-id", "direction", "outdir"))
direction <- args$direction
stopifnot(direction %in% c("up", "down"))
dirs <- ensure_output_dirs(args$outdir)

cfg <- read_project(args[["project-config"]])
curation <- cfg$analysis$settings$networks$community_curation[[direction]]
if (is.null(curation)) curation <- list()

# Editorial text follows the contrast numerator so the panel is study-agnostic.
contrasts <- readr::read_tsv(args$contrasts, show_col_types = FALSE, progress = FALSE)
contrast_row <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
numerator <- if (nrow(contrast_row) && !is.na(contrast_row$numerator[[1]])) contrast_row$numerator[[1]] else "numerator"
direction_word <- if (direction == "up") "upregulated" else "downregulated"

seed <- as.integer(curation$seed %||% if (direction == "up") 241L else 242L)
# Display policy: "main_component" shows only the largest connected component
# (bespoke up); "connected" shows every node with at least one retained edge
# (bespoke down). Default matches the bespoke up panel.
display_policy <- curation$display %||% if (direction == "up") "main_component" else "connected"

stem <- file.path(dirs$figures, sprintf("string_%s_network_community", direction))

nodes_raw <- readr::read_tsv(args$nodes, show_col_types = FALSE, progress = FALSE)
edges_raw <- readr::read_tsv(args$edges, show_col_types = FALSE, progress = FALSE)

no_network <- is.null(nodes_raw) || !nrow(nodes_raw) || is.null(edges_raw) || !nrow(edges_raw)
if (no_network) {
  save_plot_pair(
    empty_plot(sprintf("STRING network of genes %s in %s", direction_word, numerator),
               "No connected STRING subnetwork for the directional seed set"),
    stem, width = 8.2, height = 6.2
  )
  message(sprintf("string_network[%s]: no connected subnetwork; wrote placeholder", direction))
  quit(save = "no", status = 0)
}

nodes <- nodes_raw %>%
  dplyr::distinct(preferredName, .keep_all = TRUE) %>%
  dplyr::transmute(
    name = preferredName,
    string_id = stringId,
    log2fc = suppressWarnings(as.numeric(log2_fold_change)),
    padj = suppressWarnings(as.numeric(padj)),
    de_statistic = suppressWarnings(as.numeric(stat))
  )
edges <- edges_raw %>%
  dplyr::transmute(from = source, to = target, string_score = suppressWarnings(as.numeric(combined_score))) %>%
  dplyr::distinct(from, to, .keep_all = TRUE)

graph_all <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes %>% dplyr::select(name))

# Restrict to the displayed vertex set per the display policy before detecting
# communities, matching the bespoke derivation (communities are computed on the
# graph that is actually drawn).
if (display_policy == "main_component") {
  comp <- igraph::components(graph_all)
  keep <- igraph::V(graph_all)$name[comp$membership == which.max(comp$csize)]
} else {
  keep <- igraph::V(graph_all)$name[igraph::degree(graph_all) > 0]
}
graph_disp <- igraph::induced_subgraph(graph_all, vids = keep)

# Community detection. igraph::cluster_louvain is sensitive to edge/vertex order
# and to the igraph RNG version, so a single seed is not reproducible across
# environments. The default is a DETERMINISTIC max-modularity multi-start: run
# many restarts, keep the highest-modularity partition (canonical-signature
# tie-break), which is order- and version-independent and recovers the reported
# optimum. A study may instead PIN an exact partition to reproduce a specific
# published figure; pinned membership is authoritative and every displayed node
# must be covered exactly once.
partition_signature <- function(memb) {
  blocks <- split(names(memb), memb)
  paste(sort(vapply(blocks, function(b) paste(sort(b), collapse = ","), character(1))), collapse = " | ")
}

disp_names <- igraph::V(graph_disp)$name
pinned <- curation$pinned_partition
if (!is.null(pinned) && length(pinned)) {
  groups <- lapply(pinned, function(g) intersect(unlist(g, use.names = FALSE), disp_names))
  assigned <- unlist(groups, use.names = FALSE)
  missing_nodes <- setdiff(disp_names, assigned)
  dup_nodes <- unique(assigned[duplicated(assigned)])
  if (length(missing_nodes) || length(dup_nodes)) {
    stop(sprintf("string_network[%s]: pinned_partition must cover every displayed node exactly once; unassigned={%s} duplicated={%s}",
                 direction, paste(missing_nodes, collapse = ","), paste(dup_nodes, collapse = ",")))
  }
  membership <- stats::setNames(integer(length(disp_names)), disp_names)
  for (i in seq_along(groups)) membership[groups[[i]]] <- i
  message(sprintf("string_network[%s]: using pinned partition (%d communities)", direction, length(groups)))
} else {
  resolution <- as.numeric(curation$resolution %||% 1.0)
  restarts <- as.integer(curation$restarts %||% 200L)
  best_mod <- -Inf; membership <- NULL; best_sig <- NULL
  for (s in seq_len(restarts)) {
    set.seed(s)
    fit <- igraph::cluster_louvain(graph_disp, weights = igraph::E(graph_disp)$string_score, resolution = resolution)
    m <- igraph::membership(fit)
    mod <- igraph::modularity(graph_disp, m, weights = igraph::E(graph_disp)$string_score, resolution = resolution)
    sig <- partition_signature(m)
    if (mod > best_mod + 1e-9 || (abs(mod - best_mod) <= 1e-9 && (is.null(best_sig) || sig < best_sig))) {
      best_mod <- mod; membership <- m; best_sig <- sig
    }
  }
  message(sprintf("string_network[%s]: max-modularity Louvain over %d restarts (resolution=%.2f) -> %d communities, modularity=%.6f",
                  direction, restarts, resolution, length(unique(membership)), best_mod))
}

# Resolve each topology community to a curated label if any member gene appears
# in a curation rule's gene set (first match wins, in config order); otherwise a
# stable "Topology community N" label. Mirrors the bespoke case_when().
rules <- curation$communities %||% list()
rule_label <- function(member_names) {
  for (rule in rules) {
    genes <- unlist(rule$genes, use.names = FALSE)
    if (length(genes) && any(member_names %in% genes)) return(rule$label)
  }
  NA_character_
}
community_ids <- sort(unique(as.integer(membership)))
community_key <- do.call(rbind, lapply(community_ids, function(cid) {
  member_names <- names(membership)[as.integer(membership) == cid]
  label <- rule_label(member_names)
  if (is.na(label)) label <- paste("Topology community", cid)
  data.frame(community_id = cid, community = label, stringsAsFactors = FALSE)
}))

set.seed(seed)
xy <- igraph::layout_with_fr(graph_disp, weights = igraph::E(graph_disp)$string_score, niter = 5000, grid = "nogrid")
xy[, 1] <- scales::rescale(xy[, 1], to = c(-1, 1))
xy[, 2] <- scales::rescale(xy[, 2], to = c(-0.92, 0.92))

disp_nodes <- tibble::tibble(
  name = igraph::V(graph_disp)$name,
  x = xy[, 1],
  y = xy[, 2],
  community_id = as.integer(membership[igraph::V(graph_disp)$name]),
  graph_degree = as.numeric(igraph::degree(graph_disp))
) %>%
  dplyr::left_join(community_key, by = "community_id") %>%
  dplyr::left_join(nodes, by = "name")

edge_plot <- igraph::as_data_frame(graph_disp, what = "edges") %>%
  dplyr::transmute(from, to, string_score) %>%
  dplyr::left_join(disp_nodes %>% dplyr::select(from = name, x, y), by = "from") %>%
  dplyr::left_join(disp_nodes %>% dplyr::select(to = name, xend = x, yend = y), by = "to")

# --- scales: curated when provided, else sensible data-ranged defaults --------
present_communities <- unique(disp_nodes$community)
curated_colours <- stats::setNames(
  vapply(rules, function(r) r$colour %||% NA_character_, character(1)),
  vapply(rules, function(r) r$label, character(1))
)
DEFAULT_COMMUNITY_PALETTE <- c(
  "#C94F5E", "#536FA8", "#278B8A", "#D49A27", "#8567AF",
  "#A26493", "#738F3B", "#D18932", "#2B8C7E", "#6C8EBF"
)
community_colours <- stats::setNames(rep(NA_character_, length(present_communities)), present_communities)
for (nm in present_communities) {
  if (!is.na(curated_colours[nm]) && nm %in% names(curated_colours)) community_colours[nm] <- curated_colours[[nm]]
}
open_slots <- is.na(community_colours)
if (any(open_slots)) {
  fillers <- rep(DEFAULT_COMMUNITY_PALETTE, length.out = sum(open_slots))
  community_colours[open_slots] <- fillers
}

grad <- curation$fill_gradient
if (!is.null(grad)) {
  fill_scale <- scale_fill_gradientn(
    colours = unlist(grad$colours, use.names = FALSE),
    values = scales::rescale(unlist(grad$values, use.names = FALSE)),
    limits = unlist(grad$limits, use.names = FALSE),
    breaks = unlist(grad$breaks, use.names = FALSE),
    oob = scales::squish, name = "Shrunken\nlog2FC"
  )
} else if (direction == "up") {
  hi <- max(c(disp_nodes$log2fc, 1), na.rm = TRUE)
  fill_scale <- scale_fill_gradientn(
    colours = c("#FFF8F6", "#F4A6A6", "#D96363", "#8E2530"),
    limits = c(0, hi), oob = scales::squish, name = "Shrunken\nlog2FC"
  )
} else {
  lo <- min(c(disp_nodes$log2fc, -1), na.rm = TRUE)
  fill_scale <- scale_fill_gradientn(
    colours = c("#214B75", "#4F83AC", "#A6CEE3", "#F7FBFD"),
    limits = c(lo, 0), oob = scales::squish, name = "Shrunken\nlog2FC"
  )
}

size_cur <- curation$size
size_range <- if (!is.null(size_cur$range)) unlist(size_cur$range, use.names = FALSE) else c(2.8, 7.2)
size_limits <- if (!is.null(size_cur$limits)) unlist(size_cur$limits, use.names = FALSE) else c(1, max(c(disp_nodes$graph_degree, 1), na.rm = TRUE))
size_breaks <- if (!is.null(size_cur$breaks)) unlist(size_cur$breaks, use.names = FALSE) else waiver()

panel <- ggplot() +
  ggforce::geom_mark_hull(
    data = disp_nodes, aes(x, y, group = community, colour = community),
    fill = "#F7F9FA", alpha = 0.55, expand = unit(2.7, "mm"), concavity = 4,
    radius = unit(1.5, "mm"), linewidth = 0.55, linetype = "22", show.legend = FALSE
  ) +
  geom_segment(
    data = edge_plot, aes(x, y, xend = xend, yend = yend, linewidth = string_score),
    colour = "#B8C9D4", alpha = 0.68, lineend = "round"
  ) +
  geom_point(
    data = disp_nodes, aes(x, y, fill = log2fc, size = graph_degree, colour = community),
    shape = 21, stroke = 0.85
  ) +
  geom_text_repel(
    data = disp_nodes, aes(x, y, label = name), seed = seed, colour = NAVY, size = 2.6,
    box.padding = 0.30, point.padding = 0.22, min.segment.length = 0,
    segment.colour = "#AEB8C0", segment.size = 0.25, max.overlaps = Inf, max.time = 3
  ) +
  fill_scale +
  scale_colour_manual(values = community_colours, name = "Topology-derived\ncommunity") +
  scale_size_continuous(range = size_range, limits = size_limits, breaks = size_breaks, name = "Network\ndegree") +
  scale_linewidth_continuous(range = c(0.25, 1.65), breaks = c(0.5, 0.7, 0.9), name = "STRING\nscore") +
  guides(
    colour = guide_legend(order = 1, override.aes = list(shape = 21, fill = "white", size = 3.6, stroke = 1)),
    fill = guide_colourbar(order = 2, barheight = unit(22, "mm"), barwidth = unit(3.2, "mm")),
    size = guide_legend(order = 3),
    linewidth = guide_legend(order = 4)
  ) +
  coord_equal(clip = "off", expand = TRUE) +
  labs(
    title = sprintf("STRING network of genes %s in %s", direction_word, numerator),
    subtitle = sprintf("%d connected genes and %d interactions; communities are detected from weighted network topology",
                       nrow(disp_nodes), nrow(edge_plot))
  ) +
  theme_void(base_size = 8.2) +
  theme(
    text = element_text(colour = NAVY),
    plot.title = element_text(face = "bold", colour = NAVY, size = 12.2, margin = margin(b = 2)),
    plot.subtitle = element_text(colour = MID_GREY, size = 8.2, margin = margin(b = 2)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 7.4),
    legend.text = element_text(size = 6.8),
    legend.key.height = unit(4.2, "mm"),
    plot.margin = margin(6, 6, 4, 8)
  )

save_plot_pair(panel, stem, width = 8.2, height = 6.2)
message(sprintf("string_network[%s]: wrote %s.{pdf,png} (%d nodes, %d edges, %d communities)",
                direction, stem, nrow(disp_nodes), nrow(edge_plot), length(present_communities)))
