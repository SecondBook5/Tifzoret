#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

suppressPackageStartupMessages(library(limma))

args <- parse_cli(c("project-config", "samples", "annotation", "contrasts", "contrast-id", "vst", "outdir"))
cfg <- read_project(args[["project-config"]])
dirs <- ensure_output_dirs(args$outdir)
metadata <- readr::read_tsv(args$samples, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
annotation <- read_annotation_contract(args$annotation)
contrasts <- readr::read_tsv(args$contrasts, show_col_types = FALSE, progress = FALSE)
contrast <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
if (nrow(contrast) != 1L) stop("Could not resolve contrast", call. = FALSE)
factor_name <- contrast$factor[[1]]; numerator <- contrast$numerator[[1]]; denominator <- contrast$denominator[[1]]
expression <- matrix_to_symbols(SummarizedExperiment::assay(readRDS(args$vst)), annotation)
metadata <- metadata[match(colnames(expression), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id
symbol_lookup <- setNames(rownames(expression), toupper(rownames(expression)))

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
regulon$configured_target <- regulon$target
regulon$target <- unname(symbol_lookup[toupper(regulon$target)])
confidence <- cfg$analysis$settings$regulators$confidence
if (is.null(confidence)) confidence <- c("A", "B", "C")
if ("confidence" %in% names(regulon) && any(regulon$confidence %in% c("A", "B", "C", "D", "E"))) regulon <- regulon[regulon$confidence %in% confidence, , drop = FALSE]
regulon$measured <- regulon$target %in% rownames(expression)
regulon$mor <- as.numeric(regulon$mor)
regulon$likelihood <- as.numeric(regulon$likelihood)
readr::write_tsv(regulon, file.path(dirs$tables, "regulon_edges.tsv"), na = "NA")

min_targets <- if (is.null(cfg$analysis$settings$regulators$min_targets)) 5L else as.integer(cfg$analysis$settings$regulators$min_targets)
measured <- regulon %>% filter(measured) %>% distinct(source, target, .keep_all = TRUE)
target_counts <- table(measured$source)
eligible <- names(target_counts[target_counts >= min_targets])
measured <- measured %>% filter(source %in% eligible)
if (!nrow(measured)) stop("No regulators retain the configured minimum measured targets", call. = FALSE)

signed_fallback <- function(edges, matrix, signed = TRUE) {
  rows <- lapply(split(edges, edges$source), function(group) {
    weights <- if (signed) group$mor * group$likelihood else abs(group$likelihood)
    values <- matrix[group$target, , drop = FALSE]
    score <- colSums(values * weights) / sum(abs(weights))
    data.frame(source = group$source[[1]], condition = colnames(matrix), score = score)
  })
  bind_rows(rows)
}

method <- "signed weighted target score"
signed_long <- NULL
if (requireNamespace("decoupleR", quietly = TRUE)) {
  signed_long <- tryCatch(
    decoupleR::run_viper(
      mat = expression, network = measured, .source = "source", .target = "target",
      .mor = "mor", .likelihood = "likelihood", minsize = min_targets,
      eset_filter = FALSE, pleiotropy = TRUE, verbose = FALSE
    ) %>% select(source, condition, score),
    error = function(error) {
      warnings <<- c(warnings, paste0("VIPER execution failed; used signed weighted target score: ", conditionMessage(error)))
      NULL
    }
  )
  if (!is.null(signed_long)) method <- "VIPER"
}
if (is.null(signed_long)) signed_long <- signed_fallback(measured, expression, TRUE)
unsigned_long <- signed_fallback(measured, expression, FALSE)

signed_wide <- signed_long %>% pivot_wider(names_from = condition, values_from = score)
unsigned_wide <- unsigned_long %>% pivot_wider(names_from = condition, values_from = score)
readr::write_tsv(signed_wide, file.path(dirs$tables, "dorothea_activity_scores.tsv"))
readr::write_tsv(unsigned_wide, file.path(dirs$tables, "regulator_target_program_scores.tsv"))

signed_matrix <- signed_wide %>% tibble::column_to_rownames("source") %>% as.matrix()
signed_matrix <- signed_matrix[, metadata$sample_id, drop = FALSE]
metadata$contrast_group <- stats::relevel(factor(metadata[[factor_name]]), ref = denominator)
formula_text <- gsub(paste0("\\b", factor_name, "\\b"), "contrast_group", cfg$design$formula)
design <- stats::model.matrix(stats::as.formula(formula_text), metadata)
coefficient <- grep(paste0("^contrast_group", make.names(numerator), "$"), colnames(design), value = TRUE)
if (length(coefficient) != 1L) stop("Could not resolve regulator model coefficient", call. = FALSE)
fit <- limma::eBayes(limma::lmFit(signed_matrix, design))
differential <- limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("regulator") %>%
  mutate(contrast_id = args[["contrast-id"]], numerator = numerator, denominator = denominator, method = method)
readr::write_tsv(differential, file.path(dirs$tables, "regulator_differential.tsv"), na = "NA")

top_n <- if (is.null(cfg$analysis$settings$regulators$top_regulators)) 15L else as.integer(cfg$analysis$settings$regulators$top_regulators)
selected <- head(differential$regulator, min(top_n, nrow(differential)))
display <- row_zscore(signed_matrix[selected, metadata[[factor_name]] %in% c(denominator, numerator), drop = FALSE], 1.5)
row_order <- rev(selected)
column_order <- colnames(display)[order(metadata[colnames(display), factor_name])]
heatmap <- tile_heatmap(display, row_order, column_order, legend_title = "Row-scaled\nactivity", base_size = 7.6)
displayed <- heatmap$table %>% mutate(condition = metadata[as.character(sample_id), factor_name], contrast_id = args[["contrast-id"]], method = method)
readr::write_tsv(displayed, file.path(dirs$tables, "regulator_activity_displayed.tsv"))
heatmap$plot <- heatmap$plot + labs(title = "Regulator activity", subtitle = paste0(method, "; top regulators by differential activity"))
save_plot_pair(heatmap$plot, file.path(dirs$figures, "regulator_activity"), 7.6, max(5.2, 0.34 * length(selected) + 2.0))
write_json_file(list(
  contrast_id = args[["contrast-id"]], method = method, regulators_tested = nrow(differential),
  regulon_edges = nrow(regulon), measured_edges = nrow(measured), warnings = warnings
), file.path(args$outdir, "regulators_summary.json"))
