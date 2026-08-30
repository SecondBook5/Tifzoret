#!/usr/bin/env Rscript
# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     06_regulators / regulators.R
# │ WHAT:      Transcription-factor activity inference (VIPER / decoupleR)
# │ WHY:       Infers TF activity from target-gene expression and signed prior networks;
# │            identifies which regulators drive the observed transcriptional changes
# │ HOW:       viper::viper (aREA scale method) on VST expression + regulon priors; limma contrast test on activities
# │ INPUTS:    qc/vst.rds, samples.tsv, annotation.tsv, contrasts.tsv, regulon priors (resources)
# │ PRODUCES:  modules/regulators/tables/{scores, differential}.tsv + regulator figures
# │ CALLED BY: rule contrast_regulators (workflow/rules/modules.smk)
# │ ENV:       workflow/envs/r.yaml
# └─────────────────────────────────────────────────────────────────

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(limma)
  library(patchwork)
})

args <- parse_cli(c("project-config", "samples", "annotation", "contrasts", "contrast-id", "vst", "outdir"))
cfg <- read_project(args[["project-config"]])
dirs <- ensure_output_dirs(args$outdir)
metadata <- readr::read_tsv(args$samples, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
annotation <- read_annotation_contract(args$annotation)
contrasts <- readr::read_tsv(args$contrasts, show_col_types = FALSE, progress = FALSE)
contrast <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
if (nrow(contrast) != 1L) stop("Could not resolve contrast", call. = FALSE)
# Route direction resolution through the shared resolver used by de.R/pathways.R
# so the relevel reference and numerator identity come from one code path. This
# stage only ever runs on pairwise contrasts (PAIRWISE_CONTRAST_IDS gate), so the
# guard is a documented assertion and cannot fire in practice.
resolved <- resolve_contrast(contrast, cfg$design$formula)
if (!identical(resolved$type, "pairwise")) stop("regulators stage supports pairwise contrasts only", call. = FALSE)
factor_name <- resolved$factor_name; numerator <- resolved$numerator; denominator <- resolved$denominator
expression <- matrix_to_symbols(SummarizedExperiment::assay(readRDS(args$vst)), annotation)
metadata <- metadata[match(colnames(expression), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id
symbol_lookup <- setNames(rownames(expression), toupper(rownames(expression)))

signed_fallback <- function(edges, matrix, signed = TRUE) {
  rows <- lapply(split(edges, edges$source), function(group) {
    weights <- if (signed) group$mor * group$likelihood else abs(group$likelihood)
    values <- matrix[group$target, , drop = FALSE]
    score <- colSums(values * weights) / sum(abs(weights))
    data.frame(source = group$source[[1]], condition = colnames(matrix), score = score)
  })
  bind_rows(rows)
}

# Regulator activity. Primary engine = viper::viper(method="scale"), the canonical
# aREA implementation (Alvarez et al. 2016) and the exact activity model the
# publication panels were built on. A signed regulon is passed as a named list of
# {tfmode = mor, likelihood} per source. decoupleR::run_viper and the signed
# weighted target score are ordered fallbacks so the stage still runs where the
# viper package is unavailable.
build_viper_regulons <- function(edges) {
  lapply(split(edges, edges$source), function(group) {
    group <- group[!duplicated(group$target), , drop = FALSE]
    list(
      tfmode = stats::setNames(as.numeric(group$mor), group$target),
      likelihood = stats::setNames(as.numeric(group$likelihood), group$target)
    )
  })
}

# ---------------------------------------------------------------------------
# One regulator-activity view.
# ---------------------------------------------------------------------------
# Normalize a source/target regulon, score per-sample activity (VIPER ->
# decoupleR -> signed weighted target-score fallback), fit the contrast with
# limma, and render the publication panel. Every output path is supplied by
# `out`, so the SAME code produces both the primary (signed) view and, when a
# binding prior is configured, an unsigned second view under distinct filenames.
#   * force_unsigned = TRUE collapses the prior's target modes to +1 -- a binding
#     prior establishes occupancy, not activation/repression, so it is scored as
#     unsigned target-program evidence (the same reasoning the config docs give
#     for GTRD edges without `mor`).
#   * provider_label, when set, tags the method string (e.g. VIPER ->
#     VIPER_unsigned_GTRD) so the two panels' provenance is never confused.
# With force_unsigned = FALSE and provider_label = NULL the behavior and every
# byte of output are identical to the historical single-view engine.
run_regulator_view <- function(regulon, warnings, out, force_unsigned = FALSE, provider_label = NULL) {
  if (!"source" %in% names(regulon) && "tf" %in% names(regulon)) regulon$source <- regulon$tf
  if (!all(c("source", "target") %in% names(regulon))) {
    stop("Regulon edges require source (or tf) and target columns", call. = FALSE)
  }
  if (!"mor" %in% names(regulon)) regulon$mor <- 1
  if (!"confidence" %in% names(regulon)) regulon$confidence <- "custom"
  if (!"likelihood" %in% names(regulon)) regulon$likelihood <- NA_real_
  regulon$configured_target <- regulon$target
  regulon$target <- unname(symbol_lookup[toupper(regulon$target)])
  confidence <- cfg$analysis$settings$regulators$confidence
  if (is.null(confidence)) confidence <- c("A", "B", "C")
  if ("confidence" %in% names(regulon) && any(regulon$confidence %in% c("A", "B", "C", "D", "E"))) regulon <- regulon[regulon$confidence %in% confidence, , drop = FALSE]
  mor_original <- as.character(regulon$mor)
  regulon$mor <- suppressWarnings(as.numeric(regulon$mor))
  # A non-numeric mode-of-regulation (e.g. "up"/"down") would coerce silently to
  # NA and propagate NA activity scores; fail loud so a malformed custom regulon
  # is caught rather than silently mis-scored.
  coercion_lost <- is.na(regulon$mor) & !is.na(mor_original) & nzchar(trimws(mor_original))
  if (any(coercion_lost)) {
    stop("Regulon 'mor' has non-numeric values (e.g. '", trimws(mor_original[which(coercion_lost)[1]]),
         "'); mode-of-regulation must be numeric (for example +1 / -1).", call. = FALSE)
  }
  regulon$likelihood <- as.numeric(regulon$likelihood)
  # A binding prior records occupancy, not direction: collapse every target mode
  # to +1 so the view is unsigned target-program evidence rather than signed aREA.
  if (force_unsigned) regulon$mor <- 1
  # Confidence-tier edge likelihood (Garcia-Alonso et al. 2019 DoRothEA weights):
  # A/B/C/D/E -> 1.0/0.75/0.5/0.25/0.1. The dorothea package ships no likelihood
  # column, so a genuine signed VIPER regulon derives edge likelihood from the
  # confidence tier (higher-confidence edges weigh more in aREA). An explicit
  # likelihood on a custom regulon is preserved; anything still missing -> 1.
  DOROTHEA_CONFIDENCE_LIKELIHOOD <- c(A = 1.0, B = 0.75, C = 0.5, D = 0.25, E = 0.1)
  tier_weight <- unname(DOROTHEA_CONFIDENCE_LIKELIHOOD[as.character(regulon$confidence)])
  regulon$likelihood <- ifelse(is.na(regulon$likelihood), tier_weight, regulon$likelihood)
  regulon$likelihood[is.na(regulon$likelihood)] <- 1
  regulon$measured <- regulon$target %in% rownames(expression)
  readr::write_tsv(regulon, out$edges, na = "NA")

  min_targets <- if (is.null(cfg$analysis$settings$regulators$min_targets)) 5L else as.integer(cfg$analysis$settings$regulators$min_targets)
  max_targets <- if (is.null(cfg$analysis$settings$regulators$max_targets)) Inf else as.numeric(cfg$analysis$settings$regulators$max_targets)
  measured <- regulon %>% filter(measured) %>% distinct(source, target, .keep_all = TRUE)
  target_counts <- table(measured$source)
  eligible <- names(target_counts[target_counts >= min_targets & target_counts <= max_targets])
  measured <- measured %>% filter(source %in% eligible)
  if (!nrow(measured)) stop("No regulators retain the configured measured-target window", call. = FALSE)

  method <- "signed weighted target score"
  signed_long <- NULL
  if (requireNamespace("viper", quietly = TRUE)) {
    signed_long <- tryCatch({
      regulons <- build_viper_regulons(measured)
      scores <- as.matrix(viper::viper(
        expression, regulons, method = "scale", minsize = min_targets,
        eset.filter = FALSE, cores = 1, verbose = FALSE
      ))
      as.data.frame(scores, check.names = FALSE) %>%
        tibble::rownames_to_column("source") %>%
        tidyr::pivot_longer(-source, names_to = "condition", values_to = "score")
    }, error = function(error) {
      # viper IS installed here (requireNamespace was TRUE), so a scoring error is
      # a real environment/data fault, not an absent-package downgrade. Fail loud:
      # silently falling back would publish different, non-VIPER numbers.
      stop("viper::viper is installed but failed to score regulator activity: ",
           conditionMessage(error), call. = FALSE)
    })
    if (!is.null(signed_long)) method <- "VIPER"
  }
  if (is.null(signed_long) && requireNamespace("decoupleR", quietly = TRUE)) {
    signed_long <- tryCatch(
      decoupleR::run_viper(
        mat = expression, network = measured, .source = "source", .target = "target",
        .mor = "mor", .likelihood = "likelihood", minsize = min_targets,
        eset.filter = FALSE, pleiotropy = TRUE, verbose = FALSE
      ) %>% select(source, condition, score),
      error = function(error) {
        # decoupleR is installed but failed to score. Lenient: record it and fall
        # through to the deterministic proxy. Strict: fail loud (tz_degrade stops).
        warnings <<- tz_degrade(cfg, paste0("decoupleR VIPER failed; used signed weighted target score: ", conditionMessage(error)), warnings)
        NULL
      }
    )
    if (!is.null(signed_long)) method <- "VIPER"
  }
  if (is.null(signed_long)) {
    # Reached only when neither viper nor decoupleR is installed. Record it so the
    # summary's audit trail is unambiguous that these are NOT VIPER aREA scores. In
    # strict (publication) mode this is a hard failure: a published regulator panel
    # must be built on canonical VIPER aREA scoring, never the proxy (tz_degrade stops).
    warnings <- tz_degrade(cfg, paste0(
      "No VIPER backend (viper/decoupleR) is installed; regulator activity used the ",
      "deterministic signed weighted target-score fallback -- these are NOT VIPER aREA ",
      "scores. Install viper to reproduce the canonical scoring."
    ), warnings)
    signed_long <- signed_fallback(measured, expression, TRUE)
  }
  unsigned_long <- signed_fallback(measured, expression, FALSE)
  # A binding-prior view tags its method with the provider so the signed and
  # unsigned panels are never confused (VIPER -> VIPER_unsigned_GTRD).
  final_method <- if (force_unsigned && !is.null(provider_label)) paste0(method, "_unsigned_", toupper(provider_label)) else method

  signed_wide <- signed_long %>% pivot_wider(names_from = condition, values_from = score)
  unsigned_wide <- unsigned_long %>% pivot_wider(names_from = condition, values_from = score)
  readr::write_tsv(signed_wide, out$signed_scores)
  readr::write_tsv(unsigned_wide, out$unsigned_scores)

  signed_matrix <- signed_wide %>% tibble::column_to_rownames("source") %>% as.matrix()
  signed_matrix <- signed_matrix[, metadata$sample_id, drop = FALSE]
  metadata$contrast_group <- stats::relevel(factor(metadata[[factor_name]]), ref = resolved$reference_levels[[factor_name]])
  formula_text <- gsub(paste0("\\b", factor_name, "\\b"), "contrast_group", cfg$design$formula)
  design <- stats::model.matrix(stats::as.formula(formula_text), metadata)
  coefficient <- grep(paste0("^contrast_group", make.names(numerator), "$"), colnames(design), value = TRUE)
  if (length(coefficient) != 1L) stop("Could not resolve regulator model coefficient", call. = FALSE)
  fit <- limma::eBayes(limma::lmFit(signed_matrix, design))
  differential <- limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("regulator") %>%
    mutate(contrast_id = args[["contrast-id"]], numerator = numerator, denominator = denominator, method = final_method)
  readr::write_tsv(differential, out$differential, na = "NA")

  top_n <- if (is.null(cfg$analysis$settings$regulators$top_regulators)) 15L else as.integer(cfg$analysis$settings$regulators$top_regulators)
  selected <- head(differential$regulator, min(top_n, nrow(differential)))
  display <- row_zscore(signed_matrix[selected, metadata[[factor_name]] %in% c(denominator, numerator), drop = FALSE], 1.5)
  row_order <- rev(selected)
  column_order <- colnames(display)[order(metadata[colnames(display), factor_name])]
  heatmap <- tile_heatmap(display, row_order, column_order, legend_title = "Row-scaled\nactivity", base_size = 7.6)
  displayed <- heatmap$table %>% mutate(condition = metadata[as.character(sample_id), factor_name], contrast_id = args[["contrast-id"]], method = final_method)
  readr::write_tsv(displayed, out$displayed)
  # ---------------------------------------------------------------------------
  # Publication panel (Figure 2, Panel D): regulator activity.
  # ---------------------------------------------------------------------------
  # Presentation parity with the finalized paper panel. The reference styling
  # library (the prepare_regulator_activity_panel / make_regulator_activity_panel
  # routines in the reference bespoke figure functions) renders this as a two-part
  # patchwork composite rather than a plain
  # tile heatmap: a diverging tile heatmap of row-scaled activity (with a sample
  # header strip and per-row direction chips) beside a diverging effect-size
  # lollipop of differential activity with FDR annotations. We replicate that
  # construction here using the engine's ALREADY-COMPUTED variables -- `display`
  # (row z-scores, clamped +/-1.5 by row_zscore) and `differential` (limma logFC /
  # adj.P.Val). No statistic is recomputed and the audit table written above is
  # untouched. Direction is keyed off logFC sign so the numerator maps to the warm
  # styling and the denominator to the cool styling.
  # Reference diverging palette (not present in utils.R, defined locally):
  DENOMINATOR_FILL <- "#A6CEE3"; NUMERATOR_FILL <- "#F4A6A6"
  DENOMINATOR_INK <- "#39799C"; NUMERATOR_INK <- "#B55252"
  fill_by_dir <- stats::setNames(c(DENOMINATOR_FILL, NUMERATOR_FILL), c(denominator, numerator))
  ink_by_dir <- stats::setNames(c(DENOMINATOR_INK, NUMERATOR_INK), c(denominator, numerator))

  # Per-regulator effect-size / significance table (mirrors regulator_key).
  reg_stats <- differential %>%
    filter(regulator %in% selected) %>%
    transmute(
      regulator, logFC, adj.P.Val,
      higher_in = if_else(logFC >= 0, numerator, denominator)
    )

  # Order each direction block by hierarchical clustering of the row z-scores
  # (reference clustered_ids), numerator block on top then denominator block.
  cluster_block <- function(direction) {
    ids <- selected[selected %in% reg_stats$regulator[reg_stats$higher_in == direction]]
    if (length(ids) <= 2L) return(ids)
    fit <- stats::hclust(stats::dist(display[ids, , drop = FALSE]), method = "complete")
    ids[fit$order]
  }
  display_ids <- c(cluster_block(numerator), cluster_block(denominator))
  n_regulators <- length(display_ids)
  n_denominator <- sum(reg_stats$higher_in == denominator)
  boundary_y <- n_denominator + 0.5

  reg_key <- reg_stats %>%
    mutate(
      display_order = match(regulator, display_ids),
      y = n_regulators + 1L - display_order,
      regulator_display = stringr::str_replace_all(regulator, "_", "/"),
      fdr_label = if_else(
        adj.P.Val < 0.001, "q < 0.001",
        paste0("q = ", formatC(adj.P.Val, format = "f", digits = 3))
      ),
      fdr_fontface = if_else(adj.P.Val < 0.05, "bold", "plain")
    ) %>%
    arrange(display_order)

  # Sample header key: denominator columns first, then numerator.
  sample_key <- tibble::tibble(sample_id = colnames(display)) %>%
    mutate(condition = as.character(metadata[sample_id, factor_name])) %>%
    arrange(factor(condition, levels = c(denominator, numerator)), sample_id) %>%
    group_by(condition) %>%
    mutate(sample_number = row_number()) %>%
    ungroup() %>%
    mutate(sample_label = paste0(condition, "\n", sample_number), sample_x = row_number())
  n_samples <- nrow(sample_key)

  heat_long <- as.data.frame(display, check.names = FALSE) %>%
    tibble::rownames_to_column("regulator") %>%
    tidyr::pivot_longer(-regulator, names_to = "sample_id", values_to = "activity_z") %>%
    left_join(select(reg_key, regulator, y, higher_in), by = "regulator") %>%
    left_join(select(sample_key, sample_id, sample_x, condition), by = "sample_id")

  y_limits <- c(0.45, n_regulators + 1.65)
  activity_label <- if (grepl("VIPER", final_method, fixed = TRUE)) "Row-scaled VIPER activity" else "Row-scaled activity"

  heatmap_plot <- ggplot() +
    geom_tile(
      data = heat_long, aes(sample_x, y, fill = activity_z),
      width = 0.98, height = 0.94, colour = "white", linewidth = 0.42
    ) +
    geom_tile(
      data = filter(sample_key, condition == denominator),
      aes(sample_x, y = n_regulators + 1.05),
      width = 0.98, height = 0.78, fill = DENOMINATOR_FILL, colour = "white", linewidth = 0.42
    ) +
    geom_tile(
      data = filter(sample_key, condition == numerator),
      aes(sample_x, y = n_regulators + 1.05),
      width = 0.98, height = 0.78, fill = NUMERATOR_FILL, colour = "white", linewidth = 0.42
    ) +
    geom_text(
      data = sample_key, aes(sample_x, y = n_regulators + 1.05, label = sample_label),
      colour = NAVY, size = 2.3, lineheight = 0.88, fontface = "bold"
    ) +
    geom_tile(
      data = filter(reg_key, higher_in == numerator),
      aes(x = 0.28, y = y), width = 0.18, height = 0.94, fill = NUMERATOR_FILL
    ) +
    geom_tile(
      data = filter(reg_key, higher_in == denominator),
      aes(x = 0.28, y = y), width = 0.18, height = 0.94, fill = DENOMINATOR_FILL
    ) +
    geom_hline(yintercept = boundary_y, colour = "white", linewidth = 2.0) +
    scale_fill_gradient2(
      low = DENOMINATOR_INK, mid = "#F7F4EE", high = NUMERATOR_INK, midpoint = 0,
      limits = c(-1.5, 1.5), breaks = c(-1.5, 0, 1.5), oob = scales::squish,
      name = activity_label
    ) +
    scale_x_continuous(limits = c(0.12, n_samples + 0.52), expand = c(0, 0)) +
    scale_y_continuous(
      limits = y_limits, breaks = reg_key$y, labels = reg_key$regulator_display,
      expand = c(0, 0)
    ) +
    guides(fill = guide_colourbar(
      title.position = "top", title.hjust = 0.5, direction = "horizontal",
      barwidth = grid::unit(34, "mm"), barheight = grid::unit(3.3, "mm")
    )) +
    labs(x = NULL, y = NULL) +
    theme_publication(8.0) +
    theme(
      axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
      axis.text.y = element_text(face = "bold", size = 7.4, margin = margin(r = 4)),
      panel.grid = element_blank(),
      legend.position = "bottom", legend.margin = margin(t = 1),
      legend.title = element_text(size = 6.56), legend.text = element_text(size = 6.24),
      plot.margin = margin(2, 4, 2, 5)
    )

  # The reference hardcodes study-tuned x-limits c(-5.8, 8.2); the engine derives the
  # equivalent proportions from the data so the lollipop generalizes across
  # studies without clipping (data domain + a right-hand FDR-label column).
  lfc_abs <- max(abs(reg_key$logFC), na.rm = TRUE)
  if (!is.finite(lfc_abs) || lfc_abs <= 0) lfc_abs <- 1
  data_extent <- lfc_abs * 1.08
  fdr_pad <- lfc_abs * 0.62
  fdr_x <- data_extent + fdr_pad * 0.95
  x_break <- signif(lfc_abs * 0.9, 1)
  x_breaks <- unique(c(-x_break, 0, x_break))

  effect_plot <- ggplot(reg_key, aes(y = y)) +
    annotate(
      "rect", xmin = -Inf, xmax = 0, ymin = 0.5, ymax = n_regulators + 0.5,
      fill = DENOMINATOR_FILL, alpha = 0.10
    ) +
    annotate(
      "rect", xmin = 0, xmax = Inf, ymin = 0.5, ymax = n_regulators + 0.5,
      fill = NUMERATOR_FILL, alpha = 0.10
    ) +
    geom_vline(xintercept = 0, colour = "#87939D", linewidth = 0.55) +
    geom_hline(yintercept = boundary_y, colour = "white", linewidth = 2.0) +
    geom_segment(
      aes(x = 0, xend = logFC, yend = y, colour = higher_in),
      linewidth = 1.05, lineend = "round"
    ) +
    geom_point(
      aes(x = logFC, fill = higher_in, colour = higher_in),
      shape = 21, size = 3.8, stroke = 0.75
    ) +
    geom_text(
      aes(x = fdr_x, label = fdr_label, fontface = fdr_fontface),
      hjust = 1, colour = MID_GREY, size = 2.35
    ) +
    annotate(
      "text", x = fdr_x, y = n_regulators + 1.05, label = "FDR", hjust = 1,
      colour = NAVY, fontface = "bold", size = 2.65
    ) +
    scale_colour_manual(values = ink_by_dir, guide = "none") +
    scale_fill_manual(values = fill_by_dir, guide = "none") +
    scale_x_continuous(limits = c(-data_extent, data_extent + fdr_pad), breaks = x_breaks, expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits, breaks = NULL, expand = c(0, 0)) +
    labs(
      title = "Differential activity",
      x = paste0("Activity logFC (", numerator, " − ", denominator, ")"), y = NULL
    ) +
    theme_publication(8.0) +
    theme(
      plot.title = element_text(size = 8.3, margin = margin(b = 4)),
      axis.title.x = element_text(size = 7.36),
      axis.text.x = element_text(size = 6.56),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = LIGHT_GREY, linewidth = 0.30),
      plot.margin = margin(2, 5, 17, 4)
    )

  panel <- heatmap_plot + effect_plot +
    patchwork::plot_layout(widths = c(1.34, 1.0)) +
    patchwork::plot_annotation(
      title = "Regulator activity",
      subtitle = paste0("Top ", n_regulators, " regulators by differential activity; row-scaled ", final_method, " scores"),
      theme = theme(
        plot.title = element_text(face = "bold", colour = NAVY, size = 12.2, margin = margin(b = 2)),
        plot.subtitle = element_text(colour = MID_GREY, size = 8.2, margin = margin(b = 5)),
        plot.margin = margin(6, 6, 4, 8)
      )
    )
  save_plot_pair(panel, out$figure_stem, 8.6, max(4.4, 0.42 * n_regulators + 1.8))
  write_json_file(list(
    contrast_id = args[["contrast-id"]], method = final_method, regulators_tested = nrow(differential),
    regulon_edges = nrow(regulon), measured_edges = nrow(measured), warnings = warnings
  ), out$summary)
  invisible(differential)
}

# ---------------------------------------------------------------------------
# Output-path scheme for one view. The FIRST (primary) view keeps the canonical
# filenames the downstream GRN/hypothesis stages and figure recipes read; every
# later view writes `<file>_<name>` variants. With name = "binding" the non-primary
# scheme reproduces the legacy binding-prior filenames exactly, so the generalized
# views path and the built-in dual-view path stay filename-compatible. This naming
# is mirrored by contrast_regulators' declared outputs in modules.smk -- the two
# MUST agree or Snakemake reports a missing output.
# ---------------------------------------------------------------------------
regulator_out <- function(dirs, outdir, name, primary) {
  if (primary) {
    list(
      edges = file.path(dirs$tables, "regulon_edges.tsv"),
      signed_scores = file.path(dirs$tables, "dorothea_activity_scores.tsv"),
      unsigned_scores = file.path(dirs$tables, "regulator_target_program_scores.tsv"),
      differential = file.path(dirs$tables, "regulator_differential.tsv"),
      displayed = file.path(dirs$tables, "regulator_activity_displayed.tsv"),
      figure_stem = file.path(dirs$figures, "regulator_activity"),
      summary = file.path(outdir, "regulators_summary.json")
    )
  } else {
    list(
      edges = file.path(dirs$tables, paste0("regulon_edges_", name, ".tsv")),
      signed_scores = file.path(dirs$tables, paste0("regulator_activity_scores_", name, ".tsv")),
      unsigned_scores = file.path(dirs$tables, paste0("regulator_target_program_scores_", name, ".tsv")),
      differential = file.path(dirs$tables, paste0("regulator_differential_", name, ".tsv")),
      displayed = file.path(dirs$tables, paste0("regulator_activity_displayed_", name, ".tsv")),
      figure_stem = file.path(dirs$figures, paste0("regulator_activity_", name)),
      summary = file.path(outdir, paste0("regulators_", name, "_summary.json"))
    )
  }
}

# ---------------------------------------------------------------------------
# Generalized named-views driver (analysis.settings.regulators.views). Score an
# arbitrary list of source/target edge tables, each as its own view, in one run.
# `signed = FALSE` collapses the view's target modes to +1 (occupancy, not
# direction) and tags the method with `provider_label` (defaults to the view
# name). The first view is the primary. Replaces the built-in primary(+binding)
# driver below; the config layer enforces the two are mutually exclusive.
# ---------------------------------------------------------------------------
run_regulator_views <- function(views) {
  for (i in seq_along(views)) {
    view <- views[[i]]
    signed <- isTRUE(view$signed)
    label <- if (is.null(view$provider_label)) view$name else view$provider_label
    edges_path <- resolve_path(cfg$.base, view$edges)
    edges <- readr::read_tsv(edges_path, show_col_types = FALSE, progress = FALSE)
    edges$provider <- label
    view_warnings <- list()
    if (!signed) {
      view_warnings <- c(view_warnings, paste0(
        "Regulator view '", view$name, "': edges carry no direction; target modes set to +1 ",
        "and scored as unsigned target-program activity (", label, ")."
      ))
    } else if (!"mor" %in% names(edges)) {
      view_warnings <- c(view_warnings, paste0(
        "Regulator view '", view$name, "': signed edges have no mor column; activity uses +1 ",
        "and should be interpreted as unsigned target-program activity."
      ))
    }
    run_regulator_view(
      edges, view_warnings,
      out = regulator_out(dirs, args$outdir, view$name, primary = (i == 1L)),
      force_unsigned = !signed,
      provider_label = if (!signed) label else NULL
    )
  }
}

views <- cfg$analysis$settings$regulators$views
if (!is.null(views) && length(views) > 0) {
  run_regulator_views(views)
} else {

# ---------------------------------------------------------------------------
# Primary view (signed): a custom regulon or DoRothEA. Byte-identical to the
# historical single-view engine.
# ---------------------------------------------------------------------------
warnings <- list()
if (!is.null(cfg$resources$regulon_edges)) {
  regulon_path <- resolve_path(cfg$.base, cfg$resources$regulon_edges)
  regulon <- readr::read_tsv(regulon_path, show_col_types = FALSE, progress = FALSE)
  if (!"mor" %in% names(regulon)) {
    regulon$mor <- 1
    warnings <- c(warnings, "Custom regulon edges have no mor column; signed activity uses +1 and should be interpreted as unsigned target-program activity.")
  }
  if (!"confidence" %in% names(regulon)) regulon$confidence <- "custom"
  if (!"likelihood" %in% names(regulon)) regulon$likelihood <- 1
  regulon$provider <- if (isTRUE(cfg$resources$providers$gtrd)) "gtrd" else "custom"
} else {
  if (!requireNamespace("dorothea", quietly = TRUE)) stop("DoRothEA provider requires the dorothea package", call. = FALSE)
  object_name <- if (cfg$species$provider == "mouse") "dorothea_mm" else if (cfg$species$provider == "human") "dorothea_hs" else stop("DoRothEA provider supports mouse or human", call. = FALSE)
  data(list = object_name, package = "dorothea", envir = environment())
  regulon <- get(object_name, envir = environment()) %>% as.data.frame()
  regulon$provider <- "dorothea"
}
run_regulator_view(
  regulon, warnings,
  out = list(
    edges = file.path(dirs$tables, "regulon_edges.tsv"),
    signed_scores = file.path(dirs$tables, "dorothea_activity_scores.tsv"),
    unsigned_scores = file.path(dirs$tables, "regulator_target_program_scores.tsv"),
    differential = file.path(dirs$tables, "regulator_differential.tsv"),
    displayed = file.path(dirs$tables, "regulator_activity_displayed.tsv"),
    figure_stem = file.path(dirs$figures, "regulator_activity"),
    summary = file.path(args$outdir, "regulators_summary.json")
  )
)

# ---------------------------------------------------------------------------
# Optional binding-prior view (unsigned): e.g. a GTRD occupancy snapshot scored
# in the SAME run, so a single `tifzoret run` regenerates both regulator panels
# (the signed DoRothEA-style view and the unsigned binding-prior view). The
# binding table needs only source/target; direction is intentionally discarded.
# ---------------------------------------------------------------------------
if (!is.null(cfg$resources$binding_prior_edges)) {
  binding_path <- resolve_path(cfg$.base, cfg$resources$binding_prior_edges)
  binding <- readr::read_tsv(binding_path, show_col_types = FALSE, progress = FALSE)
  binding_label <- if (isTRUE(cfg$resources$providers$gtrd)) "gtrd" else "binding"
  binding$provider <- binding_label
  binding_warnings <- list(paste0(
    "Binding-prior view: edges carry no direction; target modes set to +1 and scored as unsigned target-program activity (",
    binding_label, ")."
  ))
  run_regulator_view(
    binding, binding_warnings,
    out = list(
      edges = file.path(dirs$tables, "regulon_edges_binding.tsv"),
      signed_scores = file.path(dirs$tables, "regulator_activity_scores_binding.tsv"),
      unsigned_scores = file.path(dirs$tables, "regulator_target_program_scores_binding.tsv"),
      differential = file.path(dirs$tables, "regulator_differential_binding.tsv"),
      displayed = file.path(dirs$tables, "regulator_activity_displayed_binding.tsv"),
      figure_stem = file.path(dirs$figures, "regulator_activity_binding"),
      summary = file.path(args$outdir, "regulators_binding_summary.json")
    ),
    force_unsigned = TRUE, provider_label = binding_label
  )
}

}
