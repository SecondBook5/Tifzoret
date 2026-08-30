#!/usr/bin/env Rscript
# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     01_inputs / export_symbol_map.R
# │ WHAT:      Generate a frozen gene ID → symbol TSV from annotation contract
# │ WHY:       Provides a consistent symbol lookup for materialize_inputs.py when
# │            annotation.tsv is unavailable at ingestion time
# │ HOW:       Reads annotation.tsv, drops NA/blank symbols, writes gene_id → gene_symbol
# │ INPUTS:    annotation.tsv (from prior run or external source)
# │ PRODUCES:  symbol_map.tsv (gene_id, gene_symbol)
# │ CALLED BY: manual helper — not a Snakemake rule (see materialize_inputs.py::_load_symbol_map)
# │ ENV:       workflow/envs/r.yaml
# └─────────────────────────────────────────────────────────────────

# Offline, author-run helper -- NOT a pipeline stage.
#
# Export a pinned gene-id -> symbol map from a locally installed Bioconductor
# annotation package (an OrgDb such as org.Mm.eg.db / org.Hs.eg.db) into the
# two-column TSV that the symbol resolver consumes as its `secondary_map`
# (see materialize_inputs.py::_load_symbol_map and the
# `reference.symbol_resolution` config block).
#
# Why a frozen export rather than a live query at run time:
#   * Determinism -- the map is a version-pinned snapshot committed alongside the
#     study, so a reviewer's rerun is byte-identical regardless of when they run it.
#   * Boundary -- the engine ships this GENERATOR (no study fingerprint, no bundled
#     reference data); the resulting map is reference data and belongs in the paper
#     repo, never redistributed by the engine.
#   * No runtime R/network in materialize -- the materialize step runs in a Python
#     env; it reads this pre-built TSV offline.
#
# The output header carries a `#` provenance line (package + version + keytype)
# that the resolver skips. Multi-mapping ids resolve to the lexically smallest
# symbol -- the same deterministic tie-break the resolver applies -- and ids with
# no symbol are dropped (the resolver keeps the accession for those).
#
# Usage:
#   Rscript export_symbol_map.R \
#     --orgdb org.Mm.eg.db \
#     --output resources/org_mm_eg_db_symbols.tsv \
#     [--keytype ENSEMBL] [--symbol-column SYMBOL] \
#     [--gene-ids annotation.tsv]   # optional: restrict to ids in column 1
#
# The produced TSV is wired into a study via:
#   reference:
#     symbol_resolution: {enabled: true, secondary_map: resources/org_mm_eg_db_symbols.tsv}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

args <- parse_cli(c("orgdb", "output"))
orgdb <- args[["orgdb"]]
keytype <- if (is.null(args[["keytype"]])) "ENSEMBL" else args[["keytype"]]
symbol_column <- if (is.null(args[["symbol-column"]])) "SYMBOL" else args[["symbol-column"]]

if (!requireNamespace(orgdb, quietly = TRUE)) {
  stop(sprintf(
    "annotation package '%s' is not installed. Install it (e.g. BiocManager::install('%s')) and rerun; the engine intentionally does not bundle organism annotation databases.",
    orgdb, orgdb
  ), call. = FALSE)
}
suppressPackageStartupMessages({library(orgdb, character.only = TRUE); library(AnnotationDbi)})
db <- get(orgdb)

available <- AnnotationDbi::keytypes(db)
if (!keytype %in% available) {
  stop(sprintf("keytype '%s' is not available in %s. Available: %s", keytype, orgdb, paste(available, collapse = ", ")), call. = FALSE)
}
if (!symbol_column %in% AnnotationDbi::columns(db)) {
  stop(sprintf("symbol column '%s' is not a column of %s. Available: %s", symbol_column, orgdb, paste(AnnotationDbi::columns(db), collapse = ", ")), call. = FALSE)
}

# Restrict to the ids we actually need, when an id list is supplied. Version
# suffixes (ENSMUSG...\.\d+) are stripped so the keys match the org db's bare ids.
keys <- AnnotationDbi::keys(db, keytype = keytype)
if (!is.null(args[["gene-ids"]])) {
  id_table <- readr::read_tsv(args[["gene-ids"]], show_col_types = FALSE, progress = FALSE)
  wanted <- sub("\\.\\d+$", "", as.character(id_table[[1]]))
  keys <- intersect(keys, unique(wanted))
}

# suppressMessages: mapIds warns on 1:many; we resolve that deterministically below.
mapped <- suppressMessages(AnnotationDbi::mapIds(
  db, keys = keys, column = symbol_column, keytype = keytype,
  multiVals = function(values) sort(values[!is.na(values)])[1]
))

result <- data.frame(gene_id = names(mapped), symbol = unname(mapped), stringsAsFactors = FALSE)
result <- result[!is.na(result$symbol) & nzchar(result$symbol), , drop = FALSE]
result <- result[order(result$gene_id), , drop = FALSE]

output <- args[["output"]]
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
provenance <- sprintf(
  "# tifzoret secondary symbol map | orgdb=%s version=%s keytype=%s symbol_column=%s rows=%d",
  orgdb, as.character(utils::packageVersion(orgdb)), keytype, symbol_column, nrow(result)
)
writeLines(provenance, output)
readr::write_tsv(result, output, append = TRUE, col_names = TRUE)

message(sprintf("[export-symbol-map] wrote %d id->symbol rows to %s (%s %s)",
                nrow(result), output, orgdb, as.character(utils::packageVersion(orgdb))))
