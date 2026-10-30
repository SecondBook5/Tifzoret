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
# Differential-expression significance scheme (shared by estimand.R volcano + MA).
# ---------------------------------------------------------------------------
# Five-class palette copied byte-for-byte from the old pipeline's
# SIGNIFICANCE_PALETTE (workflow/stages/de/de_packages.R) so the engine's DE
# plots match the published figures exactly. Single source of truth: both the
# volcano and MA constructors in estimand.R (and therefore de_overview) consume this
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
  # Match config.py::_resolve so validate-time and run-time agree: expand a
  # leading ~ (or ~user) and any $VAR / ${VAR} from the environment BEFORE the
  # relative-vs-absolute decision. An unset variable is a hard error (never
  # silently resolve to a wrong or empty path). Longest tokens substitute first
  # so "$VAR" cannot corrupt a "$VARLONG" occurrence.
  value <- path.expand(value)
  tokens <- regmatches(value, gregexpr("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}|\\$[A-Za-z_][A-Za-z0-9_]*", value))[[1]]
  tokens <- unique(tokens)
  for (token in tokens[order(-nchar(tokens))]) {
    name <- gsub("^\\$\\{?|\\}$", "", token)
    resolved <- Sys.getenv(name, unset = NA_character_)
    if (is.na(resolved)) stop("required environment variable is not set: ", name, call. = FALSE)
    value <- gsub(token, resolved, value, fixed = TRUE)
  }
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
      cache = "~/.cache/tifzoret/resources",
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

# Single source of truth for interpreting one contrast row, shared by de_confirm.R,
# pathways.R, composition.R, and regulators.R. Pairwise rows (the default, and every
# row in a study without the optional columns) reproduce historical behavior exactly:
# the global design, relevel `factor` to `denominator`, extract `factor_numerator_vs_denominator`.
# Coefficient rows carry a per-row design, explicit reference levels, and a
# named resultsNames() coefficient (a difference-in-differences interaction).
# Estimand rows (a weighted cell contrast from a declared family) carry no design
# of their own: the family's fitted model and cell weights define them, so only
# callers that read those (allow_estimand = TRUE) may resolve one.
resolve_contrast <- function(contrast_row, global_design, allow_estimand = FALSE) {
  type <- contrast_field(contrast_row, "type")
  if (!nzchar(type)) type <- "pairwise"
  allowed <- c("pairwise", "coefficient", "omnibus", if (allow_estimand) "estimand")
  if (!type %in% allowed) {
    stop("contrast 'type' must be ", paste0("'", allowed[-length(allowed)], "'", collapse = ", "),
         " or '", allowed[[length(allowed)]], "', got: ", type, call. = FALSE)
  }
  factor_name <- contrast_field(contrast_row, "factor")
  numerator <- contrast_field(contrast_row, "numerator")
  denominator <- contrast_field(contrast_row, "denominator")

  design_text <- contrast_field(contrast_row, "design")
  if (!nzchar(design_text)) design_text <- global_design
  design_formula <- stats::as.formula(design_text)

  reference_levels <- parse_reference_levels(contrast_field(contrast_row, "reference_levels"))
  coefficient_name <- contrast_field(contrast_row, "coefficient")

  # Omnibus rows carry a reduced design formula for the DESeq2 LRT (test="LRT");
  # the reduced model drops the tested factor. Non-omnibus rows leave these empty.
  reduced_text <- ""
  reduced_formula <- NULL
  if (identical(type, "pairwise")) {
    # Encode today's implicit relevel (factor -> denominator) so both branches
    # share one relevel path, and reconstruct the DESeq2 coefficient name.
    if (is.null(reference_levels[[factor_name]])) {
      reference_levels[[factor_name]] <- denominator
    }
    coefficient_name <- paste0(factor_name, "_", numerator, "_vs_", denominator)
  } else if (identical(type, "coefficient")) {
    if (!nzchar(coefficient_name)) {
      stop("coefficient contrast requires a non-empty 'coefficient' column", call. = FALSE)
    }
  } else if (identical(type, "estimand")) {
    # Nothing to resolve here: the caller compiles the contrast from the family fit.
  } else {
    reduced_text <- contrast_field(contrast_row, "reduced")
    if (!nzchar(reduced_text)) {
      stop("omnibus contrast requires a non-empty 'reduced' column", call. = FALSE)
    }
    reduced_formula <- stats::as.formula(reduced_text)
  }

  list(
    type = type,
    factor_name = factor_name,
    numerator = numerator,
    denominator = denominator,
    design_formula = design_formula,
    design_text = design_text,
    reference_levels = reference_levels,
    coefficient_name = coefficient_name,
    reduced_formula = reduced_formula,
    reduced_text = reduced_text
  )
}

# Fit a DESeqDataSet, degrading to gene-wise dispersions when the parametric
# trend cannot be fitted.
#
# DESeq()'s parametric dispersion-trend fit fails HARD -- it errors, it does not
# warn -- when every gene-wise estimate sits within two orders of magnitude of
# the minimum. That happens whenever biological variability is tight and uniform:
# synthetic fixtures, small templates, tightly controlled designs. DESeq2's own
# error text prescribes the remedy (use the gene-wise estimates as final), which
# is exactly what the fallback below does.
#
# Any script that fits or REFITS a model needs this, not just the canonical fit:
# sva.R refits the design augmented with surrogate variables and aborted the
# whole workflow on this error until it shared this path. Only that one DESeq2
# message is caught -- every other error still propagates, so a genuinely broken
# fit is never silently downgraded.
#
# Returns a list of the fitted object and which trend was used, so callers can
# report it (the choice changes downstream transform selection: see estimand.R,
# which prefers log2 over VST under a gene-wise fit).
fit_deseq_with_dispersion_fallback <- function(dds) {
  tryCatch(
    list(dds = DESeq2::DESeq(dds, fitType = "parametric", quiet = TRUE),
         dispersion_fit = "parametric"),
    error = function(error) {
      if (!grepl("all gene-wise dispersion estimates are within",
                 conditionMessage(error), fixed = TRUE)) {
        stop(error)
      }
      message("Parametric dispersion trend unavailable; using gene-wise dispersion estimates.")
      fallback <- DESeq2::estimateSizeFactors(dds)
      fallback <- DESeq2::estimateDispersionsGeneEst(fallback, quiet = TRUE)
      DESeq2::dispersions(fallback) <- S4Vectors::mcols(fallback)$dispGeneEst
      list(dds = DESeq2::nbinomWaldTest(fallback, quiet = TRUE),
           dispersion_fit = "gene-wise")
    }
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

# ComplexHeatmap analogue of save_plot_pair. ComplexHeatmap draws through grid
# rather than ggplot's device abstraction, so we open a cairo PDF and a cairo
# PNG explicitly and call draw() into each (same width/height/dpi contract as
# save_plot_pair). Returns the drawn HeatmapList so callers can recover the
# realized row/column order for their displayed-data table. legend sides mirror
# ComplexHeatmap::draw().
save_complexheatmap_pair <- function(ht, stem, width, height, dpi = 300,
                                     heatmap_legend_side = "right",
                                     annotation_legend_side = "right",
                                     merge_legends = FALSE) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  render <- function() ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = heatmap_legend_side,
    annotation_legend_side = annotation_legend_side,
    merge_legends = merge_legends
  )
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width, height = height, bg = "white")
  drawn <- render(); grDevices::dev.off()
  grDevices::png(paste0(stem, ".png"), width = width, height = height,
                 units = "in", res = dpi, bg = "white", type = "cairo")
  render(); grDevices::dev.off()
  invisible(drawn)
}

# Recover the realized top-to-bottom (or left-to-right) index order from a drawn
# ComplexHeatmap. Split heatmaps return a per-slice list in slice order; flatten
# it so the caller sees one contiguous display order.
complexheatmap_order <- function(order_value) {
  if (is.list(order_value)) unlist(order_value, use.names = FALSE) else order_value
}

read_counts_contract <- function(path) {
  tab <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  ids <- tab$gene_id
  mat <- as.matrix(tab[, setdiff(names(tab), "gene_id"), drop = FALSE])
  if (all(is.na(mat))) stop("count matrix contains no non-NA values", call. = FALSE)
  if (max(mat, na.rm = TRUE) > 2147483647) stop("count matrix exceeds 32-bit integer range; refusing silent NA coercion", call. = FALSE)
  if (any(mat != round(mat), na.rm = TRUE)) {
    warning("non-integer counts detected; rounding to nearest integer (matches DESeqDataSetFromTximport)", call. = FALSE)
    mat <- round(mat)
  }
  storage.mode(mat) <- "integer"
  rownames(mat) <- ids
  mat
}

# Read a typed quantification bundle only when the project declares RSEM inputs.
# Reject absent bundles and count/bundle identity mismatches rather than falling
# back to count-only normalization. Non-RSEM inputs retain their original path.
read_tximport_contract <- function(cfg, counts_path, counts) {
  # Leave existing integer-count and BAM studies unchanged.
  if (!identical(cfg$inputs$kind, "rsem")) return(NULL)
  # Locate the declared canonical bundle alongside counts.tsv.
  path <- file.path(dirname(counts_path), "tximport.rds")
  # Refuse missing length-aware information.
  if (!file.exists(path)) stop("RSEM tximport bundle is missing: ", path, call. = FALSE)
  # Load the serialized object produced by the actual tximport package.
  txi <- readRDS(path)
  # Require supported import semantics and all three matrices.
  if (!is.list(txi) || !identical(txi$countsFromAbundance, "no")) {
    stop("unsupported tximport countsFromAbundance", call. = FALSE)
  }
  # Check every quantification matrix against the canonical count universe.
  for (name in c("counts", "abundance", "length")) {
    # Require complete matrices with finite nonnegative values.
    if (!is.matrix(txi[[name]]) || any(!is.finite(txi[[name]])) || any(txi[[name]] < 0) ||
        !identical(dimnames(txi[[name]]), dimnames(counts))) {
      stop("invalid tximport matrix alignment: ", name, call. = FALSE)
    }
  }
  # Confirm that counts.tsv has not been replaced independently of its bundle.
  if (!isTRUE(all.equal(round(txi$counts), counts, check.attributes = TRUE))) {
    stop("canonical counts differ from the tximport bundle", call. = FALSE)
  }
  # Require usable length offsets in estimated-count mode.
  if (identical(txi$countsFromAbundance, "no") && any(txi$length <= 0)) {
    stop("nonpositive tximport effective lengths", call. = FALSE)
  }
  # Return the validated bundle without transforming estimates.
  txi
}

# Construct one DESeq2 object from the appropriate neutral quantification contract.
# Rows may be filtered for a family, but sample and feature alignment are exact.
make_deseq_dataset <- function(cfg, counts_path, counts, metadata, design_formula) {
  # Preserve the original count-only constructor for all existing studies.
  if (!identical(cfg$inputs$kind, "rsem")) {
    return(DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metadata,
                                        design = design_formula))
  }
  # Read the complete counts to verify bundle identity before family filtering.
  complete_counts <- read_counts_contract(counts_path)
  # Validate the typed bundle against its canonical count artifact.
  txi <- read_tximport_contract(cfg, counts_path, complete_counts)
  # Require exact sample ordering and a valid subset of the imported gene universe.
  if (!identical(colnames(counts), rownames(metadata)) ||
      !identical(colnames(counts), colnames(complete_counts)) ||
      anyDuplicated(rownames(counts)) || any(!rownames(counts) %in% rownames(complete_counts))) {
    stop("invalid RSEM model sample/feature alignment", call. = FALSE)
  }
  # Match the family's gene universe while preserving effective lengths.
  for (name in c("counts", "abundance", "length")) {
    txi[[name]] <- txi[[name]][rownames(counts), , drop = FALSE]
  }
  # Reject independently modified family counts before constructing its model.
  if (!isTRUE(all.equal(round(txi$counts), counts, check.attributes = TRUE))) {
    stop("family counts differ from the tximport bundle", call. = FALSE)
  }
  # Let DESeq2 preserve avgTxLength and apply its supported normalization semantics.
  DESeq2::DESeqDataSetFromTximport(txi, colData = metadata, design = design_formula)
}


read_annotation_contract <- function(path) {
  tab <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  tab <- tab[!duplicated(tab$gene_id), , drop = FALSE]
  tab$gene_symbol[is.na(tab$gene_symbol) | tab$gene_symbol == ""] <- tab$gene_id[is.na(tab$gene_symbol) | tab$gene_symbol == ""]
  tab
}

# Study-wide variance-stabilising transform for unsupervised use (PCA, sample
# distance, sample-level scores). It is blind to the design on purpose: a
# transform that knows the groups would flatter their separation. DESeq2's vst()
# fits the dispersion trend on a 1,000-gene subset and is its recommended
# transform; it refuses data with fewer than 1,000 genes of mean normalised
# count above 5 and points to varianceStabilizingTransformation(), which fits
# every gene and so serves small matrices. The choice is made by that stated
# condition, not by a setting.
blind_vst <- function(dds) {
  informative <- sum(rowMeans(DESeq2::counts(dds, normalized = TRUE)) > 5)
  if (informative >= 1000L) {
    DESeq2::vst(dds, blind = TRUE)
  } else {
    DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  }
}

# Colour and size scales for the over-representation bubble plot, derived from the
# displayed terms. A ggplot continuous scale turns every value outside its limits
# into NA and removes the point, so fixed limits silently erase terms; these always
# span the data. Breaks fall on whole numbers, stepping by the first of 1, 2, 4, 5,
# 10, 20 that leaves at most five intervals; the gene-count scale reaches at least 80.
ora_scale_limits <- function(scores, counts) {
  lower <- floor(min(scores))
  upper <- max(ceiling(max(scores)), lower + 1)
  steps <- c(1, 2, 4, 5, 10, 20)
  step <- steps[which((upper - lower) / steps <= 5)[1]]
  if (is.na(step)) step <- 20
  count_upper <- max(80, 10 * ceiling(max(counts) / 10))
  list(
    score = c(lower, upper),
    score_breaks = seq(lower, upper, by = step),
    count = c(0, count_upper),
    count_breaks = seq(10, count_upper, by = 20)
  )
}

# The multiplicity family of a curated GSVA panel: every measured gene set of the
# collection(s) the listed sets belong to. A collection is the id prefix before the
# first underscore (HALLMARK_, REACTOME_, GO_BP_ -> GO). Scoring and testing only the
# listed sets would run Benjamini-Hochberg over a list already filtered by effect
# size, so the listed ids choose what is drawn, never what is tested. A listed id
# that is not a measured set is an error: dropping it silently would hide a typo or
# a set that fell out of the size limits.
gsva_family_sets <- function(gene_sets, requested_ids) {
  missing <- setdiff(requested_ids, names(gene_sets))
  if (length(missing)) {
    stop("gsva_sets lists gene sets that are not measured in this study: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  collection <- function(ids) sub("_.*$", "", ids)
  names(gene_sets)[collection(names(gene_sets)) %in% collection(requested_ids)]
}

# Marks displayed (shrunken) effects of the genes under test that the unshrunk data
# rule out. Shrinkage moves an estimate toward zero but leaves it compatible with the
# unshrunk Wald 95% interval; a significant gene whose displayed effect lies outside
# that interval is an optimiser failure -- apeglm's posterior is bimodal for a gene
# with zero counts in one group, and its default optimiser can settle on the mode
# near zero. Only `eligible` genes (the FDR-significant ones) are judged: a
# non-significant gene is meant to be shrunk hard, so a large pull is not evidence of
# failure there. Missing or non-finite inputs are never flagged.
shrinkage_outside_ci <- function(shrunk, raw, se, eligible = TRUE, level = 0.95) {
  half_width <- stats::qnorm(1 - (1 - level) / 2) * se
  flagged <- abs(shrunk - raw) > half_width & eligible
  flagged[is.na(flagged) | !is.finite(shrunk) | !is.finite(raw) | !is.finite(se)] <- FALSE
  flagged
}

# Foreground genes of an over-representation test. Significance is the adjusted P of the
# Wald test and the effect is the UNSHRUNK estimate: that is the quantity the test
# inferred, whereas the shrunken effect is a display value that depends on whether an
# optimiser converged. A gene is "up" when its adjusted P is below `fdr` and its unshrunk
# log2 fold change is at least `min_abs_log2fc` (and positive); "down" mirrors that.
# A zero threshold classifies by sign alone. Genes without an adjusted P (filtered out
# before testing, or independent-filtered) are never in the foreground.
ora_foreground <- function(de, fdr, min_abs_log2fc) {
  tested <- !is.na(de$adjusted_p_value) & de$adjusted_p_value < fdr & !is.na(de$log2_fold_change_raw)
  raw <- de$log2_fold_change_raw
  list(
    up = unique(de$gene_symbol[tested & raw > 0 & raw >= min_abs_log2fc]),
    down = unique(de$gene_symbol[tested & raw < 0 & raw <= -min_abs_log2fc])
  )
}

# Hypergeometric over-representation of `selected` genes in every set of `gene_sets`.
# The background is per collection: the tested genes that belong to at least one set of
# the set's own collection (`set_source`, a named vector term -> collection). A gene that
# no set of the collection contains cannot be a member of any of them, so counting it as
# a non-member inflates every enrichment; restricting to annotated genes is
# clusterProfiler's definition and what a reader of "GeneRatio k/N" expects. Members are
# counted only among tested genes, and the foreground is restricted to the same background.
# The background is built from `annotation_sets` (the collection before size filtering).
run_ora <- function(gene_sets, selected, universe, direction, set_source,
                    annotation_sets = gene_sets) {
  universe <- unique(universe[!is.na(universe) & universe != ""])
  members <- lapply(gene_sets, function(genes) unique(intersect(genes, universe)))
  collection <- unname(set_source[names(gene_sets)])
  collection[is.na(collection)] <- "unspecified"
  # `annotation_sets` are the collection's sets BEFORE any size window: the window decides
  # which terms are tested, not which genes the collection can annotate.
  annotation_members <- lapply(annotation_sets, function(genes) unique(intersect(genes, universe)))
  annotation_collection <- unname(set_source[names(annotation_sets)])
  annotation_collection[is.na(annotation_collection)] <- "unspecified"
  background <- lapply(split(annotation_members, annotation_collection),
                       function(sets) unique(unlist(sets, use.names = FALSE)))
  rows <- lapply(seq_along(gene_sets), function(index) {
    name <- names(gene_sets)[[index]]
    pool <- background[[collection[[index]]]]
    set_genes <- members[[index]]
    chosen <- unique(intersect(selected, pool))
    overlap <- intersect(chosen, set_genes)
    n_chosen <- length(chosen)
    n_members <- length(set_genes)
    count <- length(overlap)
    p_value <- if (n_chosen == 0L || n_members == 0L) 1 else stats::phyper(
      count - 1L, n_members, length(pool) - n_members, n_chosen, lower.tail = FALSE
    )
    data.frame(
      pathway = name,
      direction = direction,
      count = count,
      selected_genes = n_chosen,
      pathway_genes = n_members,
      universe_genes = length(pool),
      gene_ratio = if (n_chosen) count / n_chosen else 0,
      fold_enrichment = if (n_chosen && n_members) (count / n_chosen) / (n_members / length(pool)) else 0,
      p_value = p_value,
      overlap_genes = paste(overlap, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

# The two level names of an estimand that is exactly one level minus another within a
# single factor, or NULL for anything else (a composite, a weighted contrast, an
# interaction). Labels and figure axes name "<numerator> vs <denominator>" from these.
pairwise_levels <- function(cells, weights) {
  if (length(cells) != 1L || length(weights) != 2L) return(NULL)
  positive <- names(weights)[abs(weights - 1) < 1e-9]
  negative <- names(weights)[abs(weights + 1) < 1e-9]
  if (length(positive) != 1L || length(negative) != 1L) return(NULL)
  list(numerator = positive, denominator = negative)
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

# Colours for the observed group levels, in order. Configured colours are used as given.
# A level the palette does not name (e.g. a pairwise contrast on a factor whose levels
# differ from the configured plotting group) gets a deterministic colour-blind-safe
# Okabe-Ito fallback, with a warning, so figures still draw instead of failing on NA.
condition_palette <- function(cfg, observed) {
  values <- unlist(cfg$figures$palette, use.names = TRUE)
  result <- values[observed]
  missing <- is.na(result)
  if (any(missing)) {
    pool <- setdiff(grDevices::palette.colors(9, "Okabe-Ito")[-1], values)
    if (!length(pool)) pool <- grDevices::palette.colors(9, "Okabe-Ito")[-1]
    unnamed <- unique(observed[missing])
    fallback <- stats::setNames(rep_len(unname(pool), length(unnamed)), unnamed)
    warning("figures.palette has no colour for: ", paste(unnamed, collapse = ", "),
            "; using a fallback palette", call. = FALSE)
    result[missing] <- fallback[observed[missing]]
  }
  stats::setNames(unname(result), observed)
}

# Derive a darker, slightly more saturated "ink" of a fill colour, preserving
# hue. Used for per-animal point outlines that sit on top of pale condition
# fills (a darker dot of the same hue reads as the same category). No external
# colour package required -- pure grDevices HSV arithmetic.
darken_hex <- function(hex, value_mult = 0.70, sat_add = 0.30) {
  hsv_mat <- grDevices::rgb2hsv(grDevices::col2rgb(hex))
  grDevices::hsv(
    h = hsv_mat["h", ],
    s = pmin(1, hsv_mat["s", ] + sat_add),
    v = pmax(0, pmin(1, hsv_mat["v", ] * value_mult))
  )
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

# Condition -> publication display label: keep it upper-case when it is already
# all-caps (an acronym such as WT, KO), otherwise title-case it
# (control -> "Control", "old vs young" -> "Old Vs Young"). Vectorised; the
# single source of truth for condition casing across every renderer.
cond_display <- function(x) {
  x <- as.character(x)
  ifelse(!is.na(x) & nchar(x) > 1L & x == toupper(x), x, stringr::str_to_title(x))
}

# Sample id -> publication display label: the condition (kept upper-case when it
# is already all-caps, e.g. KO; otherwise title-cased) plus the trailing
# replicate number -- control1 -> "Control 1", Ko2 -> "KO 2". Falls back to
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

# True when the project selects strict (publication) execution mode. Reads
# execution.strict from the config read by read_project(); absent => FALSE, so
# a config that predates the key behaves exactly as before.
tz_strict <- function(cfg) isTRUE(cfg$execution$strict)

# Strict-mode degrade gate. Every stage that would silently fall back to a
# degraded result (a non-canonical scoring proxy, a skipped step, a coarser
# statistic) routes that decision through here instead of appending to `warnings`
# directly. In lenient mode (the default) it returns `c(warnings, message)` --
# byte-identical to the historical `warnings <- c(warnings, message)` pattern, so
# existing outputs are unchanged. In strict mode it stops the run so a publication
# build fails loudly rather than emitting the degraded result. Use as:
#   warnings <- tz_degrade(cfg, "…", warnings)
tz_degrade <- function(cfg, message, warnings = NULL) {
  if (tz_strict(cfg)) stop("strict mode (execution.strict): ", message, call. = FALSE)
  c(warnings, message)
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

# (X'X)^-1 for a design matrix, used for hat-value diagnostics. Falls back to a
# pseudo-inverse via SVD when the crossproduct is singular so a rank-deficient
# design still reports leverage rather than aborting the diagnostics.
MASS_safe_solve <- function(design_matrix) {
  crossproduct <- t(design_matrix) %*% design_matrix
  inverse <- tryCatch(solve(crossproduct), error = function(e) NULL)
  if (!is.null(inverse)) return(inverse)
  decomposition <- svd(crossproduct)
  positive <- decomposition$d > max(decomposition$d) * 1e-12
  decomposition$v[, positive, drop = FALSE] %*%
    diag(1 / decomposition$d[positive], sum(positive)) %*%
    t(decomposition$u[, positive, drop = FALSE])
}
