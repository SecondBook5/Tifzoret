#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

args <- parse_cli(c("project-config", "scores", "samples", "contrasts", "contrast-id", "outdir"))
cfg <- read_project(args[["project-config"]]); dirs <- ensure_output_dirs(args$outdir)
settings <- cfg$analysis$settings$mediation
metadata <- readr::read_tsv(args$samples, show_col_types = FALSE) %>% as.data.frame(); rownames(metadata) <- metadata$sample_id
contrasts <- readr::read_tsv(args$contrasts, show_col_types = FALSE); contrast <- contrasts[contrasts$contrast_id == args[["contrast-id"]], , drop = FALSE]
factor_name <- contrast$factor[[1]]; numerator <- contrast$numerator[[1]]; denominator <- contrast$denominator[[1]]
scores <- readr::read_tsv(args$scores, show_col_types = FALSE)
warnings <- list(); result <- list(contrast_id = args[["contrast-id"]], status = "not_configured")
input <- data.frame()
if (is.null(settings$mediator_pathway) || is.null(settings$outcome_pathways)) {
  warnings <- c(warnings, "Mediation was enabled without mediator_pathway and outcome_pathways; no model was fit.")
} else {
  requested <- c(settings$mediator_pathway, unlist(settings$outcome_pathways))
  available <- intersect(requested, scores$pathway)
  if (length(available) != length(requested)) {
    warnings <- c(warnings, paste0("Mediation pathways missing from GSVA scores: ", paste(setdiff(requested, available), collapse = ", ")))
    result$status <- "missing_pathways"
  } else {
    wide <- scores %>% filter(pathway %in% requested) %>% pivot_longer(-pathway, names_to = "sample_id", values_to = "score") %>% pivot_wider(names_from = pathway, values_from = score) %>% left_join(metadata %>% tibble::rownames_to_column("metadata_row") %>% select(-metadata_row), by = "sample_id")
    wide$treatment <- as.numeric(factor(wide[[factor_name]], levels = c(denominator, numerator))) - 1
    wide$mediator <- wide[[settings$mediator_pathway]]
    wide$outcome <- rowMeans(wide[, unlist(settings$outcome_pathways), drop = FALSE], na.rm = TRUE)
    input <- wide
    minimum_n <- if (is.null(settings$minimum_recommended_samples)) 20L else as.integer(settings$minimum_recommended_samples)
    if (nrow(wide) < minimum_n) warnings <- c(warnings, sprintf("Mediation n=%d is below the configured recommended minimum n=%d; estimates are exploratory.", nrow(wide), minimum_n))
    if (requireNamespace("mediation", quietly = TRUE)) {
      set.seed(if (is.null(cfg$analysis$random_seed)) 1L else as.integer(cfg$analysis$random_seed))
      mediator_model <- stats::lm(mediator ~ treatment, data = wide); outcome_model <- stats::lm(outcome ~ treatment + mediator, data = wide)
      simulations <- if (is.null(settings$simulations)) 2000L else as.integer(settings$simulations)
      fit <- mediation::mediate(mediator_model, outcome_model, treat = "treatment", mediator = "mediator", boot = FALSE, sims = simulations)
      summary <- summary(fit)
      result <- list(contrast_id = args[["contrast-id"]], status = "fit", mediator_pathway = settings$mediator_pathway, outcome_pathways = unlist(settings$outcome_pathways), samples = nrow(wide), simulations = simulations, acme = unname(summary$d.avg), acme_p = unname(summary$d.avg.p), ade = unname(summary$z.avg), ade_p = unname(summary$z.avg.p), total_effect = unname(summary$tau.coef), total_effect_p = unname(summary$tau.p), proportion_mediated = unname(summary$n.avg), proportion_mediated_p = unname(summary$n.avg.p))
    } else {
      result$status <- "package_unavailable"; warnings <- c(warnings, "The mediation package was unavailable; no model was fit.")
    }
  }
}
readr::write_tsv(input, file.path(dirs$tables, "mediation_inputs.tsv")); readr::write_tsv(as.data.frame(result), file.path(dirs$tables, "mediation_results.tsv"), na = "NA")
result$warnings <- warnings; write_json_file(result, file.path(args$outdir, "mediation_summary.json"))

