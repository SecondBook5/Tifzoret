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

# ---------------------------------------------------------------------------
# Differential-expression significance scheme (shared by de.R volcano + MA).
# ---------------------------------------------------------------------------
# Five-class palette copied byte-for-byte from the old pipeline's
# SIGNIFICANCE_PALETTE (workflow/stages/de/de_packages.R) so the engine's DE
# plots match the published figures exactly. Single source of truth: both the
# volcano and MA constructors in de.R (and therefore de_overview) consume this
# palette together with classify_significance() below.
SIGNIFICANCE_PALETTE <- c(
  significant_up = "#B22222",
  significant_down = "#2166AC",
  padj_only = "#2A9D8F",
  lfc_only = "#F4A261",
  ns = "#C7CDD4"
)

# Canonical class order, reused for the factor levels and the legend breaks.
SIGNIFICANCE_CLASSES <- c("significant_up", "significant_down", "padj_only", "lfc_only", "ns")

# Classify each gene into one of the five significance classes from its adjusted
# p-value and shrunken log2 fold-change, given the FDR and |log2FC| thresholds
# (cfg$figures$de$fdr / cfg$figures$de$abs_log2fc). Mirrors the paper's de_class
# rule exactly: FDR-significant genes split by fold-change sign into up/down,
# FDR-significant but sub-threshold fold-change is padj_only, genes clearing the
# fold-change cutoff without FDR significance (or with NA padj) are lfc_only, and
# everything else is ns. Vectorised; returns a factor in the canonical order.
classify_significance <- function(adjusted_p_value, log2_fold_change, fdr, abs_log2fc) {
  classes <- dplyr::case_when(
    !is.na(adjusted_p_value) & adjusted_p_value < fdr & !is.na(log2_fold_change) & log2_fold_change >= abs_log2fc ~ "significant_up",
    !is.na(adjusted_p_value) & adjusted_p_value < fdr & !is.na(log2_fold_change) & log2_fold_change <= -abs_log2fc ~ "significant_down",
    !is.na(adjusted_p_value) & adjusted_p_value < fdr ~ "padj_only",
    !is.na(log2_fold_change) & abs(log2_fold_change) >= abs_log2fc ~ "lfc_only",
    TRUE ~ "ns"
  )
  factor(classes, levels = SIGNIFICANCE_CLASSES)
}

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

# Safely read one cell from a single-row contrast tibble. Returns "" when the
# column is absent (optional columns do not exist in pairwise-only studies) or
# when the value is NA/blank, so callers can treat "absent" and "empty" alike.
contrast_field <- function(contrast_row, column) {
  if (!column %in% names(contrast_row)) return("")
  value <- contrast_row[[column]][[1]]
  if (is.null(value) || (length(value) == 1L && is.na(value))) return("")
  trimws(as.character(value))
}

# Parse "factorA=level1;factorB=level2" into a named list keyed by factor. An
# empty string yields an empty list. Each entry must be a single factor=level.
parse_reference_levels <- function(text) {
  text <- trimws(text)
  if (!nzchar(text)) return(list())
  pieces <- trimws(strsplit(text, ";", fixed = TRUE)[[1]])
  pieces <- pieces[nzchar(pieces)]
  result <- list()
  for (piece in pieces) {
    kv <- strsplit(piece, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2L || !nzchar(trimws(kv[[1]])) || !nzchar(trimws(kv[[2]]))) {
      stop("reference_levels entry must be 'factor=level': ", piece, call. = FALSE)
    }
    result[[trimws(kv[[1]])]] <- trimws(kv[[2]])
  }
  result
}

# Single source of truth for interpreting one contrast row, shared by de.R and
# pathways.R. Pairwise rows (the default, and every row in a study without the
# optional columns) reproduce historical behavior exactly: the global design,
# relevel `factor` to `denominator`, extract `factor_numerator_vs_denominator`.
# Coefficient rows carry a per-row design, explicit reference levels, and a
# named resultsNames() coefficient (a difference-in-differences interaction).
resolve_contrast <- function(contrast_row, global_design) {
  type <- contrast_field(contrast_row, "type")
  if (!nzchar(type)) type <- "pairwise"
  if (!type %in% c("pairwise", "coefficient")) {
    stop("contrast 'type' must be 'pairwise' or 'coefficient', got: ", type, call. = FALSE)
  }
  factor_name <- contrast_field(contrast_row, "factor")
  numerator <- contrast_field(contrast_row, "numerator")
  denominator <- contrast_field(contrast_row, "denominator")

  design_text <- contrast_field(contrast_row, "design")
  if (!nzchar(design_text)) design_text <- global_design
  design_formula <- stats::as.formula(design_text)

  reference_levels <- parse_reference_levels(contrast_field(contrast_row, "reference_levels"))
  coefficient_name <- contrast_field(contrast_row, "coefficient")

  if (identical(type, "pairwise")) {
    # Encode today's implicit relevel (factor -> denominator) so both branches
    # share one relevel path, and reconstruct the DESeq2 coefficient name.
    if (is.null(reference_levels[[factor_name]])) {
      reference_levels[[factor_name]] <- denominator
    }
    coefficient_name <- paste0(factor_name, "_", numerator, "_vs_", denominator)
  } else if (!nzchar(coefficient_name)) {
    stop("coefficient contrast requires a non-empty 'coefficient' column", call. = FALSE)
  }

  list(
    type = type,
    factor_name = factor_name,
    numerator = numerator,
    denominator = denominator,
    design_formula = design_formula,
    design_text = design_text,
    reference_levels = reference_levels,
    coefficient_name = coefficient_name
  )
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
                         low = "#356D9A", mid = "#F8F7F3", high = "#C94F4F",
                         legend_title = "Row\nz-score", base_size = 8, zlimit = NULL) {
  table <- as.data.frame(matrix, check.names = FALSE) %>%
    tibble::rownames_to_column("feature") %>%
    tidyr::pivot_longer(-feature, names_to = "sample_id", values_to = "value") %>%
    mutate(
      feature = factor(feature, levels = rev(row_order)),
      sample_id = factor(sample_id, levels = column_order)
    )
  plot <- ggplot(table, aes(sample_id, feature, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = low, mid = mid, high = high, midpoint = 0,
      limits = if (is.null(zlimit)) NULL else c(-zlimit, zlimit),
      breaks = if (is.null(zlimit)) ggplot2::waiver() else c(-zlimit, 0, zlimit),
      oob = scales::squish, name = legend_title
    ) +
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

# MSigDB-style gene-set id -> display label, reproducing the published Hallmark
# casing (HALLMARK_MYC_TARGETS_V1 -> "MYC Targets V1", HALLMARK_DNA_REPAIR ->
# "Dna Repair", HALLMARK_WNT_BETA_CATENIN_SIGNALING -> "WNT/beta-catenin
# Signaling"). Non-Hallmark ids are simply de-underscored + title-cased.
prettify_gene_set_label <- function(ids) {
  as.character(ids) %>%
    stringr::str_remove("^HALLMARK_") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_title() %>%
    stringr::str_replace_all(c(
      "Myc" = "MYC",
      "Mtorc1" = "mTORC1",
      "E2f" = "E2F",
      "G2m" = "G2M",
      "Tgf Beta" = "TGF-beta",
      "Wnt Beta Catenin" = "WNT/beta-catenin"
    ))
}

# Sample id -> publication display label: the condition (kept upper-case when it
# is already all-caps, e.g. CAPE; otherwise title-cased) plus the trailing
# replicate number -- control1 -> "Control 1", Cape2 -> "CAPE 2". Falls back to
# the raw id when no condition vector is supplied.
sample_display_labels <- function(ids, conditions = NULL) {
  ids <- as.character(ids)
  if (is.null(conditions)) return(ids)
  conditions <- as.character(conditions)
  condition_display <- ifelse(
    !is.na(conditions) & nchar(conditions) > 1 & conditions == toupper(conditions),
    conditions,
    stringr::str_to_title(conditions)
  )
  replicate <- stringr::str_extract(ids, "[0-9]+$")
  ifelse(is.na(replicate), condition_display, paste0(condition_display, " ", replicate))
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

configured_gene_panels <- function(panel_cfg) {
  panels <- panel_cfg$gene_panels
  if (is.null(panels)) panels <- list()
  programs <- panel_cfg$programs
  if (!is.null(programs)) {
    for (program_id in names(programs)) {
      program <- programs[[program_id]]
      label <- if (is.null(program$label)) clean_term(program_id) else program$label
      panels[[program_id]] <- list(
        description = if (is.null(program$description)) label else program$description,
        color = program$color,
        groups = stats::setNames(list(unlist(program$genes, use.names = FALSE)), label),
        expected_direction = program$expected_direction,
        contrast = program$contrast
      )
    }
  }
  panels
}
