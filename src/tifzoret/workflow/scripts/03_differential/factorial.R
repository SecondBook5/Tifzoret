#!/usr/bin/env Rscript
# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     03_differential / factorial.R
# │ WHAT:      Factorial interaction exploratory figures (2×2 designs)
# │ WHY:       Makes interaction (difference-of-differences) legible: genes whose
# │            response to one factor depends on the level of the other
# │ HOW:       Three views: LFC_x vs LFC_y scatter, group-mean profiles, per-sample expression
# │ INPUTS:    qc/vst_expression.tsv, samples.tsv, de/{effect_x, effect_y}_results.tsv, config
# │ PRODUCES:  factorial/figures/{effect_vs_effect, interaction_profile, group_expression}.png
# │ CALLED BY: rule study_factorial (workflow/rules/core.smk)
# │ ENV:       workflow/envs/r.yaml
# └─────────────────────────────────────────────────────────────────

# factorial: three project-agnostic views that make a crossed factorial design's
# INTERACTION legible, rather than reading it off a single pairwise contrast.
# A 2x2 study (two two-level factors, e.g. genotype x tumor) asks "does the
# effect of one factor depend on the level of the other?" -- a difference of
# differences. This module renders that question three ways from the two
# configured signed-contrast "arms" (effect_x, effect_y) and the QC
# variance-stabilized expression:
#
#   1. effect_vs_effect  -- each gene's log2 fold-change in arm X against arm Y,
#      with the y=x no-interaction diagonal. Genes ON the diagonal respond the
#      same in both arms; genes OFF it are the interaction.
#   2. interaction_profile -- group-mean expression of the top-interaction genes
#      across factor B (x-axis), one line per level of factor A. Parallel lines =
#      no interaction; non-parallel lines = interaction.
#   3. group_expression -- per-sample expression of those genes across every
#      crossed group, with group means, as the raw distribution behind (2).
#
# Study-level and display-only: it adds no statistics. It reuses the two arms'
# DESeq2 fold-changes exactly as computed; the "interaction" it highlights is the
# descriptive difference of those two shrunken effects, NOT a refit interaction
# coefficient (that is the coefficient-contrast DE table, rendered elsewhere).

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "utils.R"), local = FALSE)

args <- parse_cli(c("project-config", "vst-expression", "samples", "de-x", "de-y", "outdir"))
cfg <- read_project(args[["project-config"]])
dirs <- ensure_output_dirs(args$outdir)

settings <- cfg$analysis$settings$factorial
if (is.null(settings)) {
  stop("factorial module requires analysis.settings.factorial (factors, effect_x, effect_y)", call. = FALSE)
}
factors <- as.character(unlist(settings$factors, use.names = FALSE))
if (length(factors) != 2L) {
  stop("factorial requires exactly two factors, got: ", paste(factors, collapse = ", "), call. = FALSE)
}
effect_x_id <- as.character(settings$effect_x)
effect_y_id <- as.character(settings$effect_y)
explicit_genes <- if (is.null(settings$genes)) character(0) else as.character(unlist(settings$genes, use.names = FALSE))
top_genes <- if (is.null(settings$top_genes)) 16L else as.integer(settings$top_genes)

fdr <- if (is.null(cfg$figures$de$fdr)) 0.05 else as.numeric(cfg$figures$de$fdr)
abs_log2fc <- if (is.null(cfg$figures$de$abs_log2fc)) 1.0 else as.numeric(cfg$figures$de$abs_log2fc)

warnings <- character(0)

# ---------------------------------------------------------------------------
# Inputs. vst_expression.tsv is symbol-keyed (first column gene_symbol, one
# column per sample); the two DE tables carry gene_symbol, log2_fold_change and
# adjusted_p_value. Nothing here needs the annotation map.
# ---------------------------------------------------------------------------
vst_tbl <- readr::read_tsv(normalizePath(args[["vst-expression"]], mustWork = TRUE),
                           show_col_types = FALSE, progress = FALSE)
symbol_col <- names(vst_tbl)[1]
sample_cols <- setdiff(names(vst_tbl), symbol_col)

samples <- readr::read_tsv(normalizePath(args$samples, mustWork = TRUE),
                           show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
for (fac in factors) {
  if (!fac %in% colnames(samples)) {
    stop("factorial factor '", fac, "' is not a samples.tsv column", call. = FALSE)
  }
}
# Only samples that appear as VST columns can contribute expression.
samples <- samples[samples$sample_id %in% sample_cols, , drop = FALSE]
if (!nrow(samples)) stop("no samples overlap the VST expression columns", call. = FALSE)

# Grouping column for the four-group view: the study's configured figures.group
# when present (so the palette applies), else the crossing of the two factors.
group_col <- cfg$figures$group
if (is.null(group_col) || !group_col %in% colnames(samples)) {
  group_col <- ".factorial_group"
  samples[[group_col]] <- paste(samples[[factors[1]]], samples[[factors[2]]], sep = "_")
  warnings <- c(warnings, "figures.group absent; grouping by the crossed factors")
}

read_de <- function(path) {
  de <- readr::read_tsv(normalizePath(path, mustWork = TRUE), show_col_types = FALSE, progress = FALSE)
  if (!"gene_symbol" %in% names(de)) stop("DE table lacks gene_symbol: ", path, call. = FALSE)
  num_col <- function(name) if (name %in% names(de)) suppressWarnings(as.numeric(de[[name]])) else NA_real_
  tibble::tibble(
    gene_symbol = as.character(de$gene_symbol),
    lfc = num_col("log2_fold_change"),
    se = num_col("lfc_se"),
    base_mean = num_col("base_mean"),
    padj = num_col("adjusted_p_value")
  ) %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)
}
de_x <- read_de(args[["de-x"]])
de_y <- read_de(args[["de-y"]])

# ---------------------------------------------------------------------------
# View 1 -- effect vs effect. Every gene with a finite fold-change in both arms.
# delta = LFC_y - LFC_x is the descriptive interaction (how much more the gene
# moves in arm Y than in arm X); interaction_z = delta / sqrt(se_x^2 + se_y^2)
# is the Wald statistic for that difference of two coefficients -- large only
# when the difference is both sizeable AND well estimated, so a huge fold-change
# resting on near-zero counts (large SE) does not masquerade as an interaction.
# ---------------------------------------------------------------------------
merged <- dplyr::inner_join(
  dplyr::rename(de_x, lfc_x = lfc, se_x = se, base_mean_x = base_mean, padj_x = padj),
  dplyr::rename(de_y, lfc_y = lfc, se_y = se, base_mean_y = base_mean, padj_y = padj),
  by = "gene_symbol"
) %>%
  dplyr::filter(is.finite(lfc_x), is.finite(lfc_y)) %>%
  dplyr::mutate(
    delta = lfc_y - lfc_x,
    base_mean = dplyr::coalesce(base_mean_x, base_mean_y),
    se_combined = sqrt(dplyr::coalesce(se_x, 0)^2 + dplyr::coalesce(se_y, 0)^2),
    interaction_z = ifelse(se_combined > 0, delta / se_combined, NA_real_),
    sig_either = (!is.na(padj_x) & padj_x < fdr) | (!is.na(padj_y) & padj_y < fdr),
    category = dplyr::case_when(
      sig_either & abs(delta) >= abs_log2fc ~ "Interaction",
      sig_either ~ "Concordant",
      TRUE ~ "Not significant"
    )
  )
category_levels <- c("Interaction", "Concordant", "Not significant")
merged$category <- factor(merged$category, levels = category_levels)

# ---------------------------------------------------------------------------
# Gene selection for the profile / expression views. An explicit list wins
# (the analyst's biological curation -- gene-name filtering is deliberately NOT
# in the engine, which is species/study-agnostic). Otherwise, among genes
# significant in either arm and present in the VST matrix, rank by the
# interaction Wald |z| when standard errors are available, and fall back to
# |delta| above the study-median base_mean (so low-count fold-change artifacts
# do not dominate) when they are not.
# ---------------------------------------------------------------------------
vst_symbols <- vst_tbl[[symbol_col]]
if (length(explicit_genes)) {
  selected <- explicit_genes[explicit_genes %in% vst_symbols]
  dropped <- setdiff(explicit_genes, vst_symbols)
  if (length(dropped)) warnings <- c(warnings, sprintf("configured genes absent from expression: %s", paste(dropped, collapse = ", ")))
  selection_method <- "explicit configured gene list"
} else {
  pool <- merged %>% dplyr::filter(sig_either, gene_symbol %in% vst_symbols)
  if (any(is.finite(pool$interaction_z))) {
    pool <- pool %>%
      dplyr::filter(is.finite(interaction_z)) %>%
      dplyr::arrange(dplyr::desc(abs(interaction_z)))
    selection_method <- "top interaction Wald |z| = |LFC_y - LFC_x| / sqrt(se_x^2 + se_y^2) among genes significant in either arm"
  } else {
    expressed <- merged$base_mean[is.finite(merged$base_mean) & merged$base_mean > 0]
    floor_bm <- if (length(expressed)) stats::median(expressed) else 0
    pool <- pool %>%
      dplyr::filter(!is.finite(base_mean) | base_mean >= floor_bm) %>%
      dplyr::arrange(dplyr::desc(abs(delta)))
    selection_method <- sprintf("top |LFC_y - LFC_x| among genes significant in either arm and expressed >= study-median base_mean (%.1f)", floor_bm)
    warnings <- c(warnings, "lfc_se absent from the DE tables; ranked interaction genes by |delta| over an expression floor rather than the Wald z")
  }
  selected <- utils::head(pool$gene_symbol, top_genes)
}
selected <- unique(selected)
if (!length(selected)) {
  warnings <- c(warnings, "no genes selected for the profile/expression views (none significant in either arm)")
}

# Long, per-sample expression for the selected genes, joined to sample metadata.
if (length(selected)) {
  expr_long <- vst_tbl %>%
    dplyr::filter(.data[[symbol_col]] %in% selected) %>%
    tidyr::pivot_longer(dplyr::all_of(sample_cols), names_to = "sample_id", values_to = "expression") %>%
    dplyr::rename(gene_symbol = !!symbol_col) %>%
    dplyr::inner_join(samples, by = "sample_id") %>%
    dplyr::mutate(gene_symbol = factor(gene_symbol, levels = selected))
} else {
  expr_long <- NULL
}

# ---------------------------------------------------------------------------
# Palettes. Factor-A levels get a colourblind-safe qualitative palette; the four
# groups reuse the study palette (figures.palette) where it covers them.
# ---------------------------------------------------------------------------
okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#000000")
level_order <- function(column) unique(as.character(samples[[column]]))
levels_a <- level_order(factors[1])
levels_b <- level_order(factors[2])

# Humanize a contrast id into an axis-ready effect label using the study's
# contrasts table: "numerator vs denominator" with underscores turned to spaces
# and the original casing kept (so WT stays WT, not "Wt"). Falls back to the
# de-underscored contrast id when the row or the table is unavailable, so the
# renderer stays project-agnostic and never fails on a missing lookup.
humanize_effect <- function(effect_id) {
  cid <- as.character(effect_id)
  label <- gsub("_", " ", cid)
  contrasts_path <- cfg$.contrasts
  if (!is.null(contrasts_path) && file.exists(contrasts_path)) {
    ct <- tryCatch(readr::read_tsv(contrasts_path, show_col_types = FALSE, progress = FALSE),
                   error = function(e) NULL)
    if (!is.null(ct) && "contrast_id" %in% names(ct)) {
      row <- ct[ct$contrast_id == cid, , drop = FALSE]
      if (nrow(row) && all(c("numerator", "denominator") %in% names(row))) {
        num <- trimws(as.character(row$numerator[[1]]))
        den <- trimws(as.character(row$denominator[[1]]))
        if (!is.na(num) && nzchar(num) && !is.na(den) && nzchar(den)) {
          label <- sprintf("%s vs %s", gsub("_", " ", num), gsub("_", " ", den))
        }
      }
    }
  }
  label
}

# Factor-A level -> colour. Prefer a study-palette group whose name contains the
# level (case-insensitive) so each factor level inherits the study's hue for the
# group it names (the last match is usually the "strong"/challenged group, e.g.
# the treated arm); fall back to a colourblind-safe qualitative palette when
# nothing matches.
factor_level_palette <- function(levels) {
  pal_values <- unlist(cfg$figures$palette, use.names = TRUE)
  pal_groups <- names(pal_values)
  out <- stats::setNames(rep(NA_character_, length(levels)), levels)
  for (i in seq_along(levels)) {
    # Match the level as a literal substring (case-insensitive), never a regex:
    # factor levels are arbitrary sample values, so an unbalanced bracket or paren
    # (e.g. a truncated "dose[hi" or "Cre(neo" annotation) would abort the render
    # with a TRE pattern-compilation error, and a metacharacter like "." would
    # silently match the WRONG group. fixed = TRUE makes the match literal on both
    # counts (a balanced "[KO]"/"(+)" compiles but still matches by regex meaning).
    hits <- pal_groups[nzchar(pal_groups) & grepl(tolower(levels[[i]]), tolower(pal_groups), fixed = TRUE)]
    if (length(hits)) out[[i]] <- unname(pal_values[[hits[[length(hits)]]]])
  }
  missing <- which(is.na(out))
  if (length(missing)) out[missing] <- okabe_ito[(seq_along(missing) - 1L) %% length(okabe_ito) + 1L]
  out
}
palette_a <- factor_level_palette(levels_a)

observed_groups <- unique(as.character(samples[[group_col]]))
# Order the crossed groups by the study palette's declared order when it covers
# them (its intended reading order); genes not in the palette fall to the end in
# first-appearance order.
palette_order <- names(unlist(cfg$figures$palette, use.names = TRUE))
if (length(palette_order)) {
  observed_groups <- observed_groups[order(match(observed_groups, palette_order), seq_along(observed_groups))]
}
group_palette <- condition_palette(cfg, observed_groups)
# Fill any group the study palette does not cover -- including a palette-free
# config, where condition_palette returns an empty vector -- with the
# colourblind-safe fallback, mirroring factor_level_palette. setdiff() is robust
# to a zero-length palette where is.na()-indexing was not (it aborted the group
# view with "Insufficient values in manual scale").
missing_group_colours <- setdiff(observed_groups, names(group_palette)[!is.na(group_palette)])
if (length(missing_group_colours)) {
  group_palette[missing_group_colours] <- okabe_ito[(seq_along(missing_group_colours) - 1L) %% length(okabe_ito) + 1L]
}
group_palette <- group_palette[observed_groups]

category_palette <- c(
  Interaction = unname(SIGNIFICANCE_PALETTE["significant_up"]),
  Concordant = unname(SIGNIFICANCE_PALETTE["padj_only"]),
  "Not significant" = unname(SIGNIFICANCE_PALETTE["ns"])
)

# ---------------------------------------------------------------------------
# Write the effect-vs-effect displayed table (all genes) and render the scatter.
# ---------------------------------------------------------------------------
effect_table <- merged %>%
  dplyr::transmute(gene_symbol, base_mean, lfc_x, se_x, padj_x, lfc_y, se_y, padj_y,
                   delta, interaction_z, category = as.character(category)) %>%
  dplyr::arrange(dplyr::desc(abs(dplyr::coalesce(interaction_z, delta))))
readr::write_tsv(effect_table, file.path(dirs$tables, "effect_vs_effect_displayed.tsv"))

if (nrow(merged)) {
  # Axis window is driven by the DATA DENSITY, not by the labelled genes. The
  # top-interaction genes routinely carry a very large fold-change in one arm
  # and ~0 in the other; anchoring the axes to them (the old behaviour) stretched
  # the range until the informative core -- where the great majority of genes
  # live -- collapsed into a dot. Instead: size the window to the 99th percentile
  # of |log2FC|, and CLAMP the handful of off-window genes to the border, drawn
  # as outward triangles so they stay visible and honest without rescaling.
  finite_lfc <- c(merged$lfc_x, merged$lfc_y)
  finite_lfc <- finite_lfc[is.finite(finite_lfc)]
  core <- if (length(finite_lfc)) stats::quantile(abs(finite_lfc), 0.99, names = FALSE) else 1
  axis_limit <- max(3, ceiling(core * 2.2 * 2) / 2)   # >= 3, rounded up to 0.5

  plot_data <- merged %>%
    dplyr::mutate(
      off_scale = pmax(abs(lfc_x), abs(lfc_y)) > axis_limit,
      x_plot = pmax(pmin(lfc_x, axis_limit), -axis_limit),
      y_plot = pmax(pmin(lfc_y, axis_limit), -axis_limit)
    ) %>%
    # Draw order: NS first (a faint background haze), then Concordant, then the
    # Interaction genes on top so the story sits above the crowd.
    dplyr::arrange(category == "Concordant", category == "Interaction")
  label_data <- plot_data %>% dplyr::filter(gene_symbol %in% selected)
  n_off <- sum(plot_data$off_scale, na.rm = TRUE)

  subtitle <- str_wrap(sprintf(
    "Each gene's log2 fold-change in the two contrasts. Genes on the dashed y = x line respond identically in both arms; distance from it is the interaction (top %d by Wald z labelled)%s.",
    length(selected),
    if (n_off) sprintf("; %d gene%s beyond ±%.1f clamped to the border (▲)",
                       n_off, ifelse(n_off == 1L, "", "s"), axis_limit) else ""
  ), 94)

  effect_plot <- ggplot(plot_data, aes(x_plot, y_plot)) +
    geom_hline(yintercept = 0, colour = LIGHT_GREY, linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = LIGHT_GREY, linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, colour = MID_GREY, linetype = "dashed", linewidth = 0.4) +
    geom_point(aes(colour = category, size = category, alpha = category, shape = off_scale)) +
    scale_colour_manual(values = category_palette, drop = FALSE, name = NULL) +
    scale_size_manual(values = c(Interaction = 2.0, Concordant = 1.5, "Not significant" = 0.5),
                      drop = FALSE, guide = "none") +
    scale_alpha_manual(values = c(Interaction = 0.9, Concordant = 0.75, "Not significant" = 0.14),
                       drop = FALSE, guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    coord_equal(xlim = c(-axis_limit, axis_limit), ylim = c(-axis_limit, axis_limit), clip = "on") +
    labs(
      title = "Effect vs effect: the interaction view",
      subtitle = subtitle,
      x = sprintf("log2 fold-change\n%s", humanize_effect(effect_x_id)),
      y = sprintf("log2 fold-change\n%s", humanize_effect(effect_y_id))
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2.3, alpha = 1, shape = 16))) +
    theme_publication(8.6) +
    theme(legend.position = "bottom")
  if (nrow(label_data)) {
    effect_plot <- effect_plot +
      ggrepel::geom_text_repel(
        data = label_data, aes(x_plot, y_plot, label = gene_symbol),
        size = 2.5, colour = NAVY, fontface = "plain",
        min.segment.length = 0, segment.colour = "#9FAAB3", segment.size = 0.22,
        box.padding = 0.4, point.padding = 0.2, max.overlaps = Inf, seed = 814
      )
  }
} else {
  effect_plot <- empty_plot("Effect vs effect: the interaction view", "No genes shared by the two effect contrasts")
}
save_plot_pair(effect_plot, file.path(dirs$figures, "effect_vs_effect"), 6.0, 6.4)

# ---------------------------------------------------------------------------
# View 2 -- interaction profile (group means across factor B, one line per
# level of factor A, faceted by gene).
# ---------------------------------------------------------------------------
profile_fields <- c("gene_symbol", factors[1], factors[2], "mean_expression", "sem", "n")
if (length(selected)) {
  profile <- expr_long %>%
    dplyr::mutate(.fa = factor(.data[[factors[1]]], levels = levels_a),
                  .fb = factor(.data[[factors[2]]], levels = levels_b)) %>%
    dplyr::group_by(gene_symbol, .fa, .fb) %>%
    dplyr::summarise(
      mean_expression = mean(expression),
      sem = if (dplyr::n() > 1L) stats::sd(expression) / sqrt(dplyr::n()) else 0,
      n = dplyr::n(),
      .groups = "drop"
    )
  profile_table <- profile %>%
    dplyr::transmute(gene_symbol = as.character(gene_symbol),
                     !!factors[1] := as.character(.fa),
                     !!factors[2] := as.character(.fb),
                     mean_expression, sem, n)
  readr::write_tsv(profile_table, file.path(dirs$tables, "interaction_profile_displayed.tsv"))

  profile_plot <- ggplot(profile, aes(.fb, mean_expression, colour = .fa, group = .fa)) +
    geom_line(linewidth = 0.7) +
    geom_errorbar(aes(ymin = mean_expression - sem, ymax = mean_expression + sem), width = 0.14, linewidth = 0.4) +
    geom_point(size = 2.0) +
    scale_colour_manual(values = palette_a, name = cond_display(factors[1]), labels = function(x) cond_display(x)) +
    scale_x_discrete(labels = function(x) cond_display(x)) +
    facet_wrap(~ gene_symbol, scales = "free_y") +
    labs(
      title = sprintf("Interaction profile: %s response by %s", cond_display(factors[2]), cond_display(factors[1])),
      subtitle = str_wrap(sprintf(
        "Group-mean VST expression (±SEM) of the top interaction genes, one line per %s across %s. Parallel lines mean no interaction; non-parallel lines are the interaction.",
        cond_display(factors[1]), cond_display(factors[2])), 92),
      x = cond_display(factors[2]), y = "VST expression"
    ) +
    theme_publication(8.4) +
    theme(legend.position = "bottom", panel.grid.major.x = element_blank())
  n_facets <- length(selected)
  facet_cols <- max(1, ceiling(sqrt(n_facets)))
  facet_rows <- ceiling(n_facets / facet_cols)
} else {
  readr::write_tsv(
    stats::setNames(data.frame(matrix(character(0), nrow = 0, ncol = length(profile_fields))), profile_fields),
    file.path(dirs$tables, "interaction_profile_displayed.tsv")
  )
  profile_plot <- empty_plot("Interaction profile", "No genes selected")
  facet_cols <- 1; facet_rows <- 1
}
save_plot_pair(profile_plot, file.path(dirs$figures, "interaction_profile"),
               max(4.5, 2.1 * facet_cols + 0.8), max(3.8, 2.0 * facet_rows + 1.0))

# ---------------------------------------------------------------------------
# View 3 -- four-group expression (per-sample points + group mean, faceted by
# gene), the raw distribution behind the profile.
# ---------------------------------------------------------------------------
expression_fields <- c("gene_symbol", "sample_id", "group", "expression")
if (length(selected)) {
  expression_table <- expr_long %>%
    dplyr::transmute(gene_symbol = as.character(gene_symbol), sample_id,
                     group = as.character(.data[[group_col]]), expression)
  readr::write_tsv(expression_table, file.path(dirs$tables, "group_expression_displayed.tsv"))

  expr_plot_data <- expr_long %>%
    dplyr::mutate(group = factor(as.character(.data[[group_col]]), levels = observed_groups))
  expr_plot <- ggplot(expr_plot_data, aes(group, expression, colour = group)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.58, linewidth = 0.55, colour = NAVY) +
    geom_jitter(position = position_jitter(width = 0.14, height = 0, seed = 814L), size = 1.6, alpha = 0.9) +
    scale_colour_manual(values = group_palette, guide = "none") +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
    facet_wrap(~ gene_symbol, scales = "free_y") +
    labs(
      title = sprintf("Group expression across the %d-group design", length(observed_groups)),
      subtitle = "Per-sample VST expression in each crossed group; the bar marks the group mean.",
      x = NULL, y = "VST expression"
    ) +
    theme_publication(8.4) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.major.x = element_blank())
  facet_cols3 <- max(1, ceiling(sqrt(length(selected))))
  facet_rows3 <- ceiling(length(selected) / facet_cols3)
} else {
  readr::write_tsv(
    stats::setNames(data.frame(matrix(character(0), nrow = 0, ncol = length(expression_fields))), expression_fields),
    file.path(dirs$tables, "group_expression_displayed.tsv")
  )
  expr_plot <- empty_plot("Group expression", "No genes selected")
  facet_cols3 <- 1; facet_rows3 <- 1
}
save_plot_pair(expr_plot, file.path(dirs$figures, "group_expression"),
               max(4.8, 2.2 * facet_cols3 + 0.8), max(3.8, 2.0 * facet_rows3 + 1.0))

# ---------------------------------------------------------------------------
# Summary receipt.
# ---------------------------------------------------------------------------
write_json_file(
  list(
    project_id = cfg$project$id,
    method = "descriptive difference of two shrunken DESeq2 effects (no refit)",
    factors = as.list(factors),
    effect_x = effect_x_id,
    effect_y = effect_y_id,
    group_column = group_col,
    fdr = fdr,
    abs_log2fc = abs_log2fc,
    genes_compared = nrow(merged),
    selection_method = selection_method,
    selected_genes = as.list(selected),
    warnings = warnings
  ),
  file.path(args$outdir, "factorial_summary.json")
)
