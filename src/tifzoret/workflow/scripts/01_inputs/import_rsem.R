#!/usr/bin/env Rscript
# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     01_inputs / import_rsem.R
# │ WHAT:      Import gene-level RSEM estimates with tximport into canonical integer counts
# │            plus a typed tximport bundle (counts, TPM, effective lengths)
# │ WHY:       DESeq2 must see RSEM's effective lengths, not rounded counts alone, so the
# │            import keeps length information and records every filtering decision
# │ HOW:       tximport(type = "rsem", gene level), declared count/regex filter, explicit
# │            zero-length policy; writes audit tables and an import receipt; fits no model
# │ INPUTS:    source sample manifest (sample_id, rsem), import-settings JSON
# │ PRODUCES:  counts.tsv, tximport.rds, rsem_{expected_counts,effective_lengths,tpm,
# │            model_lengths,length_adjustments,gene_universe}.tsv, rsem_import_receipt.json
# │ CALLED BY: materialize_inputs.py (subprocess) when inputs.kind is rsem
# │ ENV:       workflow/envs/rsem.yaml
# └─────────────────────────────────────────────────────────────────

# Import gene-level RSEM estimates with the actual tximport package.
# This boundary performs quantification import and filtering only; no DE fit.

# Require the canonical manifest and an explicit import settings document.
args <- commandArgs(trailingOnly = TRUE)
# Fail rather than guessing paths or import semantics.
if (length(args) != 4L) stop("expected samples, settings, counts and bundle paths")
# Load the sole required input-boundary packages.
if (!requireNamespace("tximport", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
  stop("RSEM import requires tximport and jsonlite")
}
# Read ordered sample identities and resolved RSEM file paths.
samples <- read.delim(args[[1]], check.names = FALSE, stringsAsFactors = FALSE)
# Read user-declared import and filtering choices.
settings <- jsonlite::fromJSON(args[[2]])
# Require unique, nonempty biological sample identifiers.
if (anyDuplicated(samples$sample_id) || anyNA(samples$sample_id) || any(!nzchar(samples$sample_id))) {
  stop("RSEM samples must have unique nonempty sample_id values")
}
# Name the quantification files in canonical sample order.
files <- setNames(samples$rsem, samples$sample_id)
# Invoke tximport's gene-level RSEM adapter without an invented count substitute.
txi <- tximport::tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE,
                         countsFromAbundance = settings$counts_from_abundance)
# Require complete, aligned matrices returned by the package.
for (name in c("counts", "abundance", "length")) {
  # Reject non-finite or negative quantification values.
  if (!is.matrix(txi[[name]]) || any(!is.finite(txi[[name]])) || any(txi[[name]] < 0)) {
    stop("invalid tximport matrix: ", name)
  }
  # Preserve exact sample and feature identity across all import matrices.
  if (!identical(colnames(txi[[name]]), samples$sample_id) ||
      !identical(rownames(txi[[name]]), rownames(txi$counts))) {
    stop("tximport matrix identity/order differs: ", name)
  }
}
# Match DESeqDataSetFromTximport's rounding only at the canonical count boundary.
rounded <- round(txi$counts)
# Refuse silent coercion of counts beyond DESeq2's integer range.
if (any(rounded > .Machine$integer.max)) stop("RSEM counts exceed 32-bit integer range")
# Apply the declared threshold to rounded counts, preserving its exact comparator.
keep <- if (settings$count_threshold_operator == "gt") {
  rowSums(rounded) > settings$min_total_count
} else {
  rowSums(rounded) >= settings$min_total_count
}
# Apply only an explicitly declared feature exclusion expression.
if (!is.null(settings$exclude_gene_regex) && nzchar(settings$exclude_gene_regex)) {
  keep <- keep & !grepl(settings$exclude_gene_regex, rownames(rounded), perl = TRUE)
}
# Fail on empty retained data rather than creating unusable downstream artifacts.
if (!any(keep)) stop("RSEM input filtering retained no genes")
# Record the original feature universe and deterministic retention decisions.
universe <- data.frame(gene_id = rownames(rounded), retained = keep)
# Subset every matrix together so length information cannot lose its alignment.
for (name in c("counts", "abundance", "length")) txi[[name]] <- txi[[name]][keep, , drop = FALSE]
# Preserve the imported lengths before any explicitly requested repair.
raw_lengths <- txi$length
# Locate every retained zero-length cell for the adjustment receipt.
zero_indices <- which(raw_lengths == 0, arr.ind = TRUE)
# Require a deliberate zero-length policy instead of silently inventing offsets.
if (nrow(zero_indices) > 0 && settings$zero_length_policy == "error") {
  stop("retained genes have zero effective lengths; declare zero_length_policy: one only with processing justification")
}
# Apply the documented one-base replacement only when explicitly configured.
if (settings$zero_length_policy == "one") txi$length[txi$length == 0] <- 1
# Record the original and used length for every modified gene/sample cell.
adjustments <- data.frame(gene_id = rownames(raw_lengths)[zero_indices[, 1]],
                          sample_id = colnames(raw_lengths)[zero_indices[, 2]],
                          original_length = rep(0, nrow(zero_indices)),
                          model_length = rep(1, nrow(zero_indices)))
# Serialize all tables with sufficient precision for numerical validation.
options(digits = 17)
# Persist the cell-level transformation audit, including an empty table if none.
write.table(adjustments, file.path(dirname(args[[4]]), "rsem_length_adjustments.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
# Keep model lengths distinct from the original imported effective lengths.
write.table(data.frame(gene_id = rownames(txi$length), txi$length, check.names = FALSE),
            file.path(dirname(args[[4]]), "rsem_model_lengths.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# Emit integer canonical counts while preserving fractional estimates in the bundle.
counts <- round(txi$counts)
# Write the canonical count matrix with explicit stable identifiers.
write.table(data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
            args[[3]], sep = "\t", quote = FALSE, row.names = FALSE)
# Preserve the full tximport object for DESeqDataSetFromTximport downstream.
saveRDS(txi, args[[4]])
# Export auditable fractional estimates, effective lengths and imported TPMs.
for (name in c("counts", "length", "abundance")) {
  # Keep descriptive filenames distinct from integer counts.tsv.
  output_name <- switch(name, counts = "rsem_expected_counts.tsv", length = "rsem_effective_lengths.tsv", abundance = "rsem_tpm.tsv")
  # Write matrices without changing feature or sample order.
  write.table(data.frame(gene_id = rownames(txi[[name]]), if (name == "length") raw_lengths else txi[[name]], check.names = FALSE),
              file.path(dirname(args[[4]]), output_name), sep = "\t", quote = FALSE, row.names = FALSE)
}
# Retain feature-filter provenance independently of the fitted family universe.
write.table(universe, file.path(dirname(args[[4]]), "rsem_gene_universe.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
# Record import settings and actual package versions; make no reproduction claim.
jsonlite::write_json(list(import_type = "rsem", txIn = FALSE, txOut = FALSE,
                         settings = settings, tximport_version = as.character(utils::packageVersion("tximport")),
                         R_version = R.version.string, zero_length_cells = nrow(zero_indices), genes_before = nrow(universe), genes_after = sum(keep)),
                     file.path(dirname(args[[4]]), "rsem_import_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
