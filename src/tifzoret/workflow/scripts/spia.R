#!/usr/bin/env Rscript

# SPIA: Signaling Pathway Impact Analysis. Over-representation (fgsea/ORA) asks
# "are the DE genes enriched in this pathway's member list?"; SPIA additionally
# asks "given WHERE the DE genes sit in the pathway's wiring and HOW they change,
# is the pathway perturbed, and activated or inhibited?". It combines two
# independent lines of evidence -- the over-representation p (pNDE) and a
# topology-aware perturbation p (pPERT) computed by propagating the DE log2 fold
# changes through the KEGG reaction graph -- into a single global p (pG) and a
# net-accumulation sign (tA > 0 activated, < 0 inhibited).
#
# Per-contrast, opt-in (off in every profile), and requires the KEGG provider
# (resources.providers.kegg). Topology comes from graphite's KEGG graphs; if the
# graphite/SPIA machinery or the pathway database is unavailable, the module
# warns and emits empty-but-well-formed outputs rather than failing the run.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

args <- parse_cli(c("project-config", "de", "contrast-id", "outdir"))
cfg <- read_project(args[["project-config"]])
dirs <- ensure_output_dirs(args$outdir)

settings <- cfg$analysis$settings$spia
if (is.null(settings)) settings <- list()
fdr <- if (is.null(cfg$figures$de$fdr)) 0.05 else as.numeric(cfg$figures$de$fdr)
if (!is.null(settings$fdr)) fdr <- as.numeric(settings$fdr)
top_pathways <- if (is.null(settings$top_pathways)) 25L else as.integer(settings$top_pathways)

PATHWAY_FIELDS <- c("Name", "pSize", "NDE", "pNDE", "tA", "pPERT", "pG", "pGFdr", "pGFWER", "Status", "direction", "contrast_id")

# Emit empty-but-well-formed outputs plus a placeholder figure and a summary that
# records why, so an unavailable topology database degrades gracefully.
write_skip <- function(reason) {
  warning(reason, call. = FALSE)
  readr::write_tsv(
    stats::setNames(data.frame(matrix(character(0), nrow = 0, ncol = length(PATHWAY_FIELDS))), PATHWAY_FIELDS),
    file.path(dirs$tables, "spia_pathways.tsv")
  )
  readr::write_tsv(
    stats::setNames(data.frame(matrix(character(0), nrow = 0, ncol = length(PATHWAY_FIELDS))), PATHWAY_FIELDS),
    file.path(dirs$tables, "spia_displayed.tsv")
  )
  save_plot_pair(empty_plot("SPIA pathway impact", reason), file.path(dirs$figures, "spia_two_evidence"), 6.6, 5.0)
  write_json_file(
    list(project_id = cfg$project$id, contrast_id = args[["contrast-id"]], method = "SPIA (graphite KEGG topology)",
         fdr = fdr, pathways = 0L, significant = 0L, skipped = TRUE, reason = reason),
    file.path(args$outdir, "spia_summary.json")
  )
  quit(save = "no", status = 0)
}

if (!requireNamespace("SPIA", quietly = TRUE) || !requireNamespace("graphite", quietly = TRUE)) {
  write_skip("SPIA and/or graphite not installed; skipped pathway-topology impact analysis")
}

provider <- cfg$species$provider
orgdb_package <- if (provider == "mouse") "org.Mm.eg.db" else if (provider == "human") "org.Hs.eg.db" else NA_character_
graphite_species <- if (provider == "mouse") "mmusculus" else if (provider == "human") "hsapiens" else NA_character_
if (is.na(orgdb_package) || !requireNamespace(orgdb_package, quietly = TRUE)) {
  write_skip(sprintf("SPIA requires the species-matched org.*.eg.db package (provider '%s')", provider))
}
orgdb <- get(orgdb_package, asNamespace(orgdb_package))

de <- readr::read_tsv(normalizePath(args$de, mustWork = TRUE), show_col_types = FALSE, progress = FALSE) %>% as.data.frame()

# Map the DE table's Ensembl gene ids to Entrez (SPIA/KEGG identifiers). Genes
# that fail to map are dropped from both the DE set and the measured universe.
entrez <- suppressWarnings(AnnotationDbi::mapIds(orgdb, keys = de$gene_id, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"))
de$entrez <- unname(entrez[de$gene_id])
de <- de[!is.na(de$entrez) & de$entrez != "", , drop = FALSE]
de <- de[!duplicated(de$entrez), , drop = FALSE]
if (nrow(de) == 0L) write_skip("no DE genes could be mapped to Entrez ids; skipped SPIA")

all_entrez <- de$entrez
significant <- de[de$direction %in% c("up_in_numerator", "down_in_numerator"), , drop = FALSE]
if (nrow(significant) == 0L) write_skip("no significant DE genes for this contrast; skipped SPIA")
de_lfc <- stats::setNames(as.numeric(significant$log2_fold_change), significant$entrez)

# Build the graphite KEGG topology and the SPIA data files under a private
# working directory (prepareSPIA writes "<db>SPIA.RData" into getwd()). Any
# failure here -- missing data package, no network for the graphite database --
# degrades to a graceful skip.
result <- tryCatch({
  db_name <- paste0("tifzoret_kegg_", args[["contrast-id"]])
  spia_dir <- file.path(dirs$objects, "spia")
  dir.create(spia_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- setwd(spia_dir)
  on.exit(setwd(old_wd), add = TRUE)

  pathway_db <- graphite::pathways(graphite_species, "kegg")
  pathway_db <- graphite::convertIdentifiers(pathway_db, "ENTREZID")
  graphite::prepareSPIA(pathway_db, db_name)

  # graphite tags nodes as "ENTREZID:<id>"; SPIA de/all names must match.
  de_named <- stats::setNames(de_lfc, paste0("ENTREZID:", names(de_lfc)))
  all_named <- paste0("ENTREZID:", all_entrez)
  graphite::runSPIA(de = de_named, all = all_named, db_name)
}, error = function(error) error)

if (inherits(result, "error")) {
  write_skip(paste0("graphite/SPIA topology unavailable: ", conditionMessage(result)))
}
if (is.null(result) || nrow(result) == 0L) {
  write_skip("SPIA returned no pathways for this contrast")
}

spia_table <- result
# SPIA Status is "Activated"/"Inhibited"; express it in the engine's numerator
# convention so it reads the same direction as the DE tables.
spia_table$direction <- ifelse(spia_table$Status == "Activated", "up_in_numerator",
                        ifelse(spia_table$Status == "Inhibited", "down_in_numerator", "not_significant"))
spia_table$contrast_id <- args[["contrast-id"]]
present_fields <- intersect(PATHWAY_FIELDS, colnames(spia_table))
spia_table <- spia_table[order(spia_table$pGFdr, -abs(spia_table$tA)), c(present_fields, setdiff(colnames(spia_table), present_fields)), drop = FALSE]
readr::write_tsv(spia_table, file.path(dirs$tables, "spia_pathways.tsv"), na = "NA")

displayed <- utils::head(spia_table, top_pathways)
readr::write_tsv(displayed, file.path(dirs$tables, "spia_displayed.tsv"), na = "NA")

# Two-evidence plot: the SPIA signature view. x = over-representation evidence
# (-log10 pNDE), y = perturbation evidence (-log10 pPERT); fill marks global
# significance (pGFdr < FDR) and outline colour marks activated vs inhibited.
plot_table <- spia_table
plot_table$neg_log10_pNDE <- -log10(pmax(plot_table$pNDE, .Machine$double.xmin))
plot_table$neg_log10_pPERT <- -log10(pmax(plot_table$pPERT, .Machine$double.xmin))
plot_table$significant <- ifelse(!is.na(plot_table$pGFdr) & plot_table$pGFdr < fdr, "Significant", "Not significant")
status_colours <- c(Activated = "#C0392B", Inhibited = "#2C6FBB", "No change" = MID_GREY)
plot_table$status_display <- ifelse(plot_table$Status %in% names(status_colours), plot_table$Status, "No change")
labelled <- utils::head(plot_table[order(plot_table$pGFdr), , drop = FALSE], min(12L, nrow(plot_table)))

if (nrow(plot_table) == 0L) {
  spia_plot <- empty_plot("SPIA pathway impact", "No pathways to display")
} else {
  spia_plot <- ggplot(plot_table, aes(neg_log10_pNDE, neg_log10_pPERT)) +
    geom_point(aes(colour = status_display, fill = status_display, shape = significant, size = pSize), stroke = 0.7, alpha = 0.9) +
    ggrepel::geom_text_repel(data = labelled, aes(label = Name), size = 2.2, colour = NAVY, max.overlaps = Inf, min.segment.length = 0) +
    scale_colour_manual(values = status_colours, name = "Perturbation") +
    scale_fill_manual(values = status_colours, guide = "none") +
    scale_shape_manual(values = c("Significant" = 21, "Not significant" = 1), name = sprintf("pG FDR < %.2g", fdr)) +
    scale_size_continuous(range = c(1.8, 6.5), name = "Pathway size") +
    labs(
      title = "SPIA pathway-topology impact",
      subtitle = sprintf("Over-representation vs perturbation evidence · %s", args[["contrast-id"]]),
      x = expression(-log[10]~"p"[NDE]~" (over-representation)"),
      y = expression(-log[10]~"p"[PERT]~" (perturbation)")
    ) +
    theme_publication(8.6)
}
save_plot_pair(spia_plot, file.path(dirs$figures, "spia_two_evidence"), 7.4, 5.6)

write_json_file(
  list(
    project_id = cfg$project$id,
    contrast_id = args[["contrast-id"]],
    method = "SPIA (graphite KEGG topology)",
    species = graphite_species,
    fdr = fdr,
    de_genes = nrow(significant),
    measured_genes = length(all_entrez),
    pathways = nrow(spia_table),
    significant = sum(!is.na(spia_table$pGFdr) & spia_table$pGFdr < fdr),
    activated = sum(spia_table$Status == "Activated", na.rm = TRUE),
    inhibited = sum(spia_table$Status == "Inhibited", na.rm = TRUE),
    skipped = FALSE
  ),
  file.path(args$outdir, "spia_summary.json")
)
