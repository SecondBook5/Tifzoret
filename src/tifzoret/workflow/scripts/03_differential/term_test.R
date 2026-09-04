#!/usr/bin/env Rscript
# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     03_differential / term_test.R
# │ WHAT:      One full-vs-reduced likelihood-ratio test within a family
# │ WHY:       Model-based corroboration of the Wald estimands: does the gene
# │            need this term at all? Generalizes the old multi-level omnibus
# │ HOW:       DESeq2 nbinomLRT against the reduced design, reusing the family's
# │            dispersions (estimated on the FULL design, as DESeq2 prescribes)
# │ INPUTS:    families/<id>/{objects/deseq2.rds, tables/tested_gene_universe.tsv}
# │ PRODUCES:  families/<id>/term_tests/<tid>/{tables,figures}/*, term_test_summary.json
# │ CALLED BY: rule family_term_test (workflow/rules/core.smk)
# │ ENV:       workflow/envs/r.yaml
# └─────────────────────────────────────────────────────────────────

# Wald is the primary inference; this is corroboration. The engine deliberately
# never intersects the two into a combined significance rule -- that would be an
# arbitrary conjunction, not a test.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "utils.R"), local = FALSE)
source(file.path(dirname(script_path), "..", "estimands.R"), local = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
})

args <- parse_cli(c("project-config", "annotation", "families", "term-tests",
                    "family-dir", "family-id", "term-test-id", "outdir"))
cfg <- read_project(args[["project-config"]])
dirs <- ensure_output_dirs(args$outdir)
annotation <- read_annotation_contract(normalizePath(args$annotation, mustWork = TRUE))

family <- read_family_spec(args$families, args[["family-id"]])
tests <- read_term_tests(args[["term-tests"]], args[["family-id"]])
row <- tests[tests$term_test_id == args[["term-test-id"]], , drop = FALSE]
if (nrow(row) != 1L) {
  stop("could not resolve exactly one term test: ", args[["term-test-id"]], call. = FALSE)
}
reduced_formula <- stats::as.formula(row$reduced[[1]])
declared_df <- as.integer(row$df[[1]])

dds <- readRDS(file.path(args[["family-dir"]], "objects", "deseq2.rds"))
# Dispersions come from the family's full-design fit, which is exactly what
# DESeq2's test="LRT" guidance prescribes, so the LRT and the Wald estimands are
# comparable gene-for-gene over one shared universe.
lrt_dds <- DESeq2::nbinomLRT(dds, reduced = reduced_formula, quiet = TRUE)
fdr <- cfg$figures$de$fdr
lrt <- DESeq2::results(lrt_dds, alpha = fdr)

result_table <- data.frame(
  gene_id = rownames(lrt),
  base_mean = lrt$baseMean,
  lrt_statistic = lrt$stat,
  df = declared_df,
  p_value = lrt$pvalue,
  adjusted_p_value = lrt$padj,
  family_id = family$family_id,
  term_test_id = args[["term-test-id"]],
  stringsAsFactors = FALSE
) %>%
  left_join(annotation, by = "gene_id") %>%
  mutate(
    gene_symbol = ifelse(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol),
    safe_p_value = ifelse(is.na(p_value), NA_real_, pmax(p_value, .Machine$double.xmin)),
    negative_log10_p = -log10(safe_p_value),
    significant = !is.na(adjusted_p_value) & adjusted_p_value < fdr
  ) %>%
  arrange(adjusted_p_value, desc(lrt_statistic))
readr::write_tsv(result_table, file.path(dirs$tables, "lrt_results.tsv"), na = "NA")

histogram_table <- function(values, bins = 40L) {
  values <- values[is.finite(values)]
  estimate <- graphics::hist(values, breaks = bins, plot = FALSE)
  data.frame(bin_left = head(estimate$breaks, -1), bin_right = tail(estimate$breaks, -1),
             bin_midpoint = estimate$mids, count = estimate$counts)
}
pvalue_histogram <- histogram_table(result_table$p_value, 40L)
readr::write_tsv(pvalue_histogram, file.path(dirs$tables, "pvalue_distribution_displayed.tsv"))
pvalue_plot <- ggplot(pvalue_histogram, aes(bin_midpoint, count)) +
  geom_col(width = stats::median(pvalue_histogram$bin_right - pvalue_histogram$bin_left),
           fill = "#6C92AE", colour = "white", linewidth = 0.15) +
  labs(title = paste0("Term test: ", args[["term-test-id"]]),
       subtitle = sprintf("LRT of %s against %s (%d df)", family$design, row$reduced[[1]], declared_df),
       x = "Raw p-value", y = "Genes") +
  theme_publication(8.8)
save_plot_pair(pvalue_plot, file.path(dirs$figures, "pvalue_distribution"), 5.8, 4.5)

write_json_file(
  list(
    project_id = cfg$project$id,
    family_id = family$family_id,
    term_test_id = args[["term-test-id"]],
    design = family$design,
    reduced = row$reduced[[1]],
    df = declared_df,
    test = "LRT",
    genes_tested = nrow(result_table),
    significant = sum(result_table$significant)
  ),
  file.path(args$outdir, "term_test_summary.json")
)
