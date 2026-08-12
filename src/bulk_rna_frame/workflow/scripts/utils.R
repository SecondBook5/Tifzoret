suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(yaml)
  library(jsonlite)
  library(scales)
  library(stringr)
})

NAVY <- "#183B56"
MID_GREY <- "#697783"
LIGHT_GREY <- "#E7ECF0"

parse_cli <- function(required) {
  args <- commandArgs(trailingOnly = TRUE)
  result <- list()
  index <- 1L
  while (index <= length(args)) {
    key <- sub("^--", "", args[[index]])
    if (index == length(args)) stop("Missing value for --", key, call. = FALSE)
    result[[key]] <- args[[index + 1L]]
    index <- index + 2L
  }
  missing <- setdiff(required, names(result))
  if (length(missing)) stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
  result
}

resolve_path <- function(base, value) {
  if (grepl("^/", value)) normalizePath(value, mustWork = FALSE) else normalizePath(file.path(base, value), mustWork = FALSE)
}

read_project <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  cfg <- yaml::read_yaml(path)
  base <- dirname(path)
  if (identical(as.integer(cfg$version), 2L)) {
    # Downstream scripts consume one stable internal shape while public v2
    # keeps analysis and resource policy grouped explicitly.
    cfg$design <- list(formula = cfg$analysis$design)
    cfg$contrasts <- cfg$analysis$contrasts
    cfg$gene_sets <- cfg$resources$gene_sets
  } else {
    cfg$analysis <- list(profile = "standard", random_seed = cfg$figures$pathways$seed)
    cfg$resources <- list(
      cache = "~/.cache/bulk-rna-frame/resources",
      offline = FALSE,
      refresh = FALSE,
      gene_sets = cfg$gene_sets,
      providers = list()
    )
    cfg$species <- list(provider = "custom", scientific_name = "unspecified", taxonomy_id = NULL)
    cfg$reference <- list(genome_build = "unspecified", annotation_release = NULL)
  }
  cfg$.config_path <- path
  cfg$.base <- base
  cfg$.samples <- resolve_path(base, cfg$inputs$samples)
  if (identical(cfg$inputs$kind, "counts")) {
    cfg$.counts <- resolve_path(base, cfg$inputs$counts)
    cfg$.annotation <- resolve_path(base, cfg$inputs$annotation)
  }
  cfg$.contrasts <- resolve_path(base, cfg$contrasts)
  cfg$.gmt <- resolve_path(base, cfg$gene_sets$gmt)
  cfg
}

ensure_output_dirs <- function(outdir) {
  dirs <- file.path(outdir, c("figures", "tables", "objects", "logs"))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  list(
    root = outdir,
    figures = file.path(outdir, "figures"),
    tables = file.path(outdir, "tables"),
    objects = file.path(outdir, "objects")
  )
}

theme_publication <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(colour = NAVY),
      plot.title = element_text(face = "bold", size = rel(1.18), margin = margin(b = 2)),
      plot.subtitle = element_text(colour = MID_GREY, size = rel(0.86), margin = margin(b = 6)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = NAVY),
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(colour = LIGHT_GREY, linewidth = 0.28),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(6, 7, 6, 7)
    )
}

save_plot_pair <- function(plot, stem, width, height, dpi = 300) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  invisible(paste0(stem, c(".pdf", ".png")))
}

read_counts_contract <- function(path) {
  tab <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  ids <- tab$gene_id
  mat <- as.matrix(tab[, setdiff(names(tab), "gene_id"), drop = FALSE])
  if (all(is.na(mat))) stop("count matrix contains no non-NA values", call. = FALSE)
  if (max(mat, na.rm = TRUE) > 2147483647) stop("count matrix exceeds 32-bit integer range; refusing silent NA coercion", call. = FALSE)
  storage.mode(mat) <- "integer"
  rownames(mat) <- ids
  mat
}

read_annotation_contract <- function(path) {
  tab <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  tab <- tab[!duplicated(tab$gene_id), , drop = FALSE]
  tab$gene_symbol[is.na(tab$gene_symbol) | tab$gene_symbol == ""] <- tab$gene_id[is.na(tab$gene_symbol) | tab$gene_symbol == ""]
  tab
}

matrix_to_symbols <- function(mat, annotation) {
  mapping <- annotation$gene_symbol[match(rownames(mat), annotation$gene_id)]
  mapping[is.na(mapping) | mapping == ""] <- rownames(mat)[is.na(mapping) | mapping == ""]
  means <- rowMeans(mat, na.rm = TRUE)
  keep <- !duplicated(mapping[order(-means)])
  chosen <- order(-means)[keep]
  out <- mat[chosen, , drop = FALSE]
  rownames(out) <- mapping[chosen]
  out
}

condition_palette <- function(cfg, observed) {
  values <- unlist(cfg$figures$palette, use.names = TRUE)
  values[observed]
}

ellipse_coordinates <- function(data, group_col, level = 0.80, points = 120L) {
  groups <- split(data, data[[group_col]])
  rows <- lapply(names(groups), function(group_name) {
    group <- groups[[group_name]]
    if (nrow(group) < 3L) return(NULL)
    covariance <- stats::cov(group[, c("PC1", "PC2")])
    if (any(!is.finite(covariance))) return(NULL)
    eig <- eigen(covariance, symmetric = TRUE)
    if (max(eig$values) <= 1e-12) return(NULL)
    # Three-point groups can yield an almost rank-one covariance estimate. A
    # literal ellipse then renders as a line and no longer communicates a
    # group envelope. Preserve its orientation and major axis while enforcing
    # a small, documented visual minor axis.
    minimum_eigenvalue <- max(eig$values) * 0.04
    stabilized_values <- pmax(eig$values, minimum_eigenvalue)
    theta <- seq(0, 2 * pi, length.out = points)
    centered <- sweep(as.matrix(group[, c("PC1", "PC2")]), 2, c(mean(group$PC1), mean(group$PC2)))
    rotated <- centered %*% eig$vectors
    observed_radius <- max(sqrt(rowSums(sweep(rotated^2, 2, stabilized_values, "/"))))
    radius <- max(sqrt(stats::qchisq(level, df = 2)), observed_radius * 1.08)
    circle <- rbind(cos(theta), sin(theta))
    coords <- t(matrix(c(mean(group$PC1), mean(group$PC2)), nrow = 2, ncol = points) +
      radius * eig$vectors %*% diag(sqrt(stabilized_values), 2) %*% circle)
    data.frame(PC1 = coords[, 1], PC2 = coords[, 2], ellipse_group = group_name)
  })
  dplyr::bind_rows(rows)
}

clean_term <- function(value) {
  value %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_sentence()
}

row_zscore <- function(mat, limit = 1.5) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z[z < -limit] <- -limit
  z[z > limit] <- limit
  z
}

tile_heatmap <- function(matrix, row_order = rownames(matrix), column_order = colnames(matrix),
                         low = "#4C78A8", mid = "#F7F4EE", high = "#D95F5F",
                         legend_title = "Row z-score", base_size = 8) {
  table <- as.data.frame(matrix, check.names = FALSE) %>%
    tibble::rownames_to_column("feature") %>%
    tidyr::pivot_longer(-feature, names_to = "sample_id", values_to = "value") %>%
    mutate(
      feature = factor(feature, levels = rev(row_order)),
      sample_id = factor(sample_id, levels = column_order)
    )
  plot <- ggplot(table, aes(sample_id, feature, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    scale_fill_gradient2(low = low, mid = mid, high = high, midpoint = 0, name = legend_title) +
    labs(x = NULL, y = NULL) +
    theme_publication(base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks = element_blank(),
      legend.position = "right"
    )
  list(plot = plot, table = table)
}

write_json_file <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(value, path, pretty = TRUE, auto_unbox = TRUE, na = "null")
}

empty_plot <- function(title, subtitle = "No displayable results for the configured thresholds") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = subtitle, colour = MID_GREY, size = 3.5) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = title) +
    theme_void(base_size = 9) +
    theme(
      text = element_text(colour = NAVY),
      plot.title = element_text(face = "bold", size = 11),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}
