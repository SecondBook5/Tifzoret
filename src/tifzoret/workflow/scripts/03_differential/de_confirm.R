#!/usr/bin/env Rscript
# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     03_differential / de_confirm.R
# │ WHAT:      Confirmatory differential expression (edgeR quasi-likelihood)
# │ WHY:       Two independent NB engines (DESeq2 + edgeR) that agree on direction
# │            and top hits increase confidence results are not method artefacts
# │ HOW:       edgeR glmQLFit + glmQLFTest on same counts/design; concordance scatter
# │ INPUTS:    inputs/{counts.tsv, samples.tsv, annotation.tsv}, contrasts.tsv, de/de_results.tsv
# │ PRODUCES:  modules/de_confirm/{tables/edger_results+concordance.tsv, figures/concordance.png}
# │ CALLED BY: rule contrast_de_confirm (workflow/rules/modules.smk)
# │ ENV:       workflow/envs/r.yaml
# └─────────────────────────────────────────────────────────────────

# Confirmatory second differential-expression engine: edgeR quasi-likelihood
# (glmQLFit + glmQLFTest) run on the SAME counts, design, and contrast that
# DESeq2 used in the family path, then cross-checked for concordance. This is a
# robustness control, not a replacement -- two independent negative-binomial
# engines that agree on direction and on the top hits give a reviewer confidence
# the calls are not an artefact of one method's shrinkage or dispersion model.
#
# Generalized (spec §10): accepts arbitrary estimands via the shared family path.
# edgeR's QL model accepts a general design matrix and a general contrast vector,
# so there is no statistical reason to restrict this to two-group comparisons.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "utils.R"), local = FALSE)
source(file.path(dirname(script_path), "..", "estimands.R"), local = FALSE)

suppressPackageStartupMessages({
  library(edgeR)
})

# Four-class concordance palette (both engines / one only / neither), reusing the
# engine's DE class colours so the scatter reads like the volcano.
CONCORDANCE_PALETTE <- c(
  both = "#B22222",
  deseq2_only = "#2166AC",
  edger_only = "#F4A261",
  neither = "#C7CDD4"
)
CONCORDANCE_CLASSES <- c("both", "deseq2_only", "edger_only", "neither")

args <- parse_cli(c("project-config", "counts", "samples", "annotation", "families", "estimands", "cell-weights", "family-dir", "estimand-id", "de", "outdir"))
cfg <- read_project(args[["project-config"]])
cfg$.counts <- normalizePath(args$counts, mustWork = TRUE)
cfg$.samples <- normalizePath(args$samples, mustWork = TRUE)
cfg$.annotation <- normalizePath(args$annotation, mustWork = TRUE)
dirs <- ensure_output_dirs(args$outdir)

counts <- read_counts_contract(cfg$.counts)
metadata <- readr::read_tsv(cfg$.samples, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
annotation <- read_annotation_contract(cfg$.annotation)

# Read family_id from the family directory's summary
family_summary_path <- file.path(args[["family-dir"]], "family_summary.json")
if (!file.exists(family_summary_path)) {
  stop("family summary not found: ", family_summary_path, call. = FALSE)
}
family_summary <- jsonlite::fromJSON(family_summary_path)
family_id <- family_summary$family_id

# Read family and estimand specifications from the shared TSV manifests
spec <- read_estimand_spec(args$estimands, args[["estimand-id"]])
family <- read_family_spec(args$families, spec$family_id)
weights <- read_cell_weights(args[["cell-weights"]], spec$estimand_id)
design_formula <- stats::as.formula(family$design)

metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id
for (field in all.vars(design_formula)) metadata[[field]] <- factor(metadata[[field]])

# edgeR's QL model accepts a general design matrix and a general contrast vector,
# so there is no statistical reason to restrict this to a two-group comparison.
# The contrast vector comes from the SAME estimands.R function the DESeq2 path
# used -- two implementations of that mapping is how two engines end up
# disagreeing on sign while both look plausible.
design <- stats::model.matrix(design_formula, data = metadata)
contrast_vector <- compile_contrast_vector(
  c(family, spec), weights, metadata, design_formula, make.names(colnames(design)))

# The family's exported universe, so concordance is computed on identical genes
# rather than on an inner join of two different filters.
universe <- readr::read_tsv(
  file.path(args[["family-dir"]], "tables", "tested_gene_universe.tsv"),
  show_col_types = FALSE, progress = FALSE)
keep <- rownames(counts) %in% universe$gene_id[universe$retained]
dge <- edgeR::DGEList(counts = counts[keep, , drop = FALSE])
dge <- edgeR::calcNormFactors(dge)
dge <- edgeR::estimateDisp(dge, design)
fit <- edgeR::glmQLFit(dge, design)
qlf <- edgeR::glmQLFTest(fit, contrast = as.numeric(contrast_vector))
edger <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table

fdr <- cfg$figures$de$fdr
lfc <- cfg$figures$de$abs_log2fc
edger_table <- data.frame(
  gene_id = rownames(edger),
  log_cpm = edger$logCPM,
  edger_log2_fold_change = edger$logFC,
  edger_f = edger$F,
  edger_p_value = edger$PValue,
  edger_fdr = edger$FDR,
  stringsAsFactors = FALSE
) %>%
  left_join(annotation, by = "gene_id") %>%
  mutate(
    gene_symbol = ifelse(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol),
    edger_significant = !is.na(edger_fdr) & edger_fdr < fdr & abs(edger_log2_fold_change) >= lfc,
    contrast_id = args[["estimand-id"]]
  ) %>%
  arrange(edger_fdr, desc(abs(edger_log2_fold_change)))
readr::write_tsv(edger_table, file.path(dirs$tables, "edger_results.tsv"), na = "NA")

# Concordance against the DESeq2 result table. Inner-join on gene_id so only
# genes tested by both engines contribute (edgeR's filterByExpr and DESeq2's
# rowSums>=10 filters differ slightly); classify each gene by which engine calls
# it significant at the shared thresholds.
deseq2 <- readr::read_tsv(args$de, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    gene_id,
    deseq2_log2_fold_change = log2_fold_change,
    deseq2_adjusted_p_value = adjusted_p_value,
    deseq2_significant = !is.na(adjusted_p_value) & adjusted_p_value < fdr & abs(log2_fold_change) >= lfc
  )
concordance <- edger_table %>%
  inner_join(deseq2, by = "gene_id") %>%
  mutate(
    concordance_class = factor(
      case_when(
        deseq2_significant & edger_significant ~ "both",
        deseq2_significant & !edger_significant ~ "deseq2_only",
        !deseq2_significant & edger_significant ~ "edger_only",
        TRUE ~ "neither"
      ),
      levels = CONCORDANCE_CLASSES
    )
  )
readr::write_tsv(
  concordance %>% select(gene_id, gene_symbol, deseq2_log2_fold_change, edger_log2_fold_change,
                         deseq2_adjusted_p_value, edger_fdr, deseq2_significant, edger_significant,
                         concordance_class, contrast_id),
  file.path(dirs$tables, "de_concordance_displayed.tsv"),
  na = "NA"
)

finite_lfc <- concordance %>% filter(is.finite(deseq2_log2_fold_change), is.finite(edger_log2_fold_change))
safe_cor <- function(x, y, method) {
  if (length(x) < 3L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x, y, method = method))
}
spearman_lfc <- safe_cor(finite_lfc$deseq2_log2_fold_change, finite_lfc$edger_log2_fold_change, "spearman")
pearson_lfc <- safe_cor(finite_lfc$deseq2_log2_fold_change, finite_lfc$edger_log2_fold_change, "pearson")
n_deseq2_sig <- sum(concordance$deseq2_significant)
n_edger_sig <- sum(concordance$edger_significant)
n_both <- sum(concordance$concordance_class == "both")
n_either <- sum(concordance$deseq2_significant | concordance$edger_significant)
jaccard <- if (n_either > 0L) n_both / n_either else NA_real_
# Sign agreement among genes significant in either engine (does the direction of
# effect agree even when significance calls differ?).
sig_either <- concordance %>% filter(deseq2_significant | edger_significant)
sign_concordance <- if (nrow(sig_either) > 0L) {
  mean(sign(sig_either$deseq2_log2_fold_change) == sign(sig_either$edger_log2_fold_change))
} else {
  NA_real_
}

subtitle <- sprintf(
  "Spearman rho = %s · %d/%d significant genes shared (Jaccard %s)",
  ifelse(is.na(spearman_lfc), "NA", formatC(spearman_lfc, format = "f", digits = 3)),
  n_both, n_either,
  ifelse(is.na(jaccard), "NA", formatC(jaccard, format = "f", digits = 3))
)
concordance_plot <- ggplot(finite_lfc, aes(deseq2_log2_fold_change, edger_log2_fold_change, colour = concordance_class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#8B979F", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "#D5DBE0", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#D5DBE0", linewidth = 0.3) +
  geom_point(size = 1.1, alpha = 0.72, na.rm = TRUE) +
  scale_colour_manual(
    values = CONCORDANCE_PALETTE, drop = FALSE,
    breaks = CONCORDANCE_CLASSES,
    labels = c(both = "Both engines", deseq2_only = "DESeq2 only", edger_only = "edgeR only", neither = "Neither"),
    name = "Significant in"
  ) +
  labs(
    title = paste0("DESeq2 vs edgeR concordance: ", spec$label),
    subtitle = subtitle,
    x = "DESeq2 log2 fold-change",
    y = "edgeR log2 fold-change"
  ) +
  theme_publication(9.0) +
  theme(legend.position = "top")
save_plot_pair(concordance_plot, file.path(dirs$figures, "de_concordance"), 6.2, 5.4)

write_json_file(
  list(
    project_id = cfg$project$id,
    contrast_id = args[["estimand-id"]],
    family_id = family$family_id,
    estimand_id = args[["estimand-id"]],
    method = "edger_quasi_likelihood",
    genes_compared = nrow(concordance),
    deseq2_significant = n_deseq2_sig,
    edger_significant = n_edger_sig,
    significant_both = n_both,
    significant_either = n_either,
    jaccard = jaccard,
    spearman_log2fc = spearman_lfc,
    pearson_log2fc = pearson_lfc,
    sign_concordance = sign_concordance
  ),
  file.path(args$outdir, "de_confirm_summary.json")
)
