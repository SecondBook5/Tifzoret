#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(KEGGREST)
  library(digest)
})

args <- parse_cli(c("project-config", "custom-gmt", "gmt", "table", "receipt"))
cfg <- read_project(args[["project-config"]])
providers <- cfg$resources$providers
if (is.null(providers)) providers <- list()
enabled <- function(name) isTRUE(providers[[name]])
# GO ontology domains to resolve (Biological Process by default; Cellular
# Component and Molecular Function are opt-in breadth). Kept BP-only by default so
# the ORA universe -- and therefore the bespoke GO-BP figure's FDR values -- is
# unchanged for existing studies.
go_domains <- unlist(cfg$resources$go_domains)
if (is.null(go_domains) || length(go_domains) == 0L) go_domains <- c("BP")
go_domains <- toupper(as.character(go_domains))
gtrd_snapshot <- NULL
if (enabled("gtrd")) {
  if (is.null(cfg$resources$regulon_edges)) {
    stop("GTRD provider requires resources.regulon_edges with an exported source/target snapshot", call. = FALSE)
  }
  gtrd_snapshot <- resolve_path(cfg$.base, cfg$resources$regulon_edges)
  if (!file.exists(gtrd_snapshot)) stop("GTRD regulon snapshot does not exist: ", gtrd_snapshot, call. = FALSE)
}

cache_root <- path.expand(if (is.null(cfg$resources$cache)) "~/.cache/tifzoret/resources" else cfg$resources$cache)
cache_key <- digest::digest(list(
  cache_schema_version = 3L,
  species = cfg$species,
  reference = cfg$reference,
  providers = providers,
  go_domains = sort(go_domains),
  collections = cfg$resources$gene_sets$collections,
  gtrd_snapshot = if (is.null(gtrd_snapshot)) NULL else digest::digest(file = gtrd_snapshot, algo = "sha256"),
  custom_gmt = digest::digest(file = args[["custom-gmt"]], algo = "sha256")
), algo = "sha256")
cache_dir <- file.path(cache_root, cache_key)
cache_gmt <- file.path(cache_dir, "gene_sets.gmt")
cache_table <- file.path(cache_dir, "gene_sets.tsv")
cache_receipt <- file.path(cache_dir, "receipt.json")
dir.create(dirname(args$gmt), recursive = TRUE, showWarnings = FALSE)

copy_cache <- function() {
  file.copy(cache_gmt, args$gmt, overwrite = TRUE)
  file.copy(cache_table, args$table, overwrite = TRUE)
  file.copy(cache_receipt, args$receipt, overwrite = TRUE)
}
if (file.exists(cache_gmt) && file.exists(cache_table) && file.exists(cache_receipt) && !isTRUE(cfg$resources$refresh)) {
  copy_cache()
  quit(save = "no", status = 0)
}
if (isTRUE(cfg$resources$offline) && any(vapply(names(providers), enabled, logical(1)))) {
  stop("Offline mode requested, but the required resource cache entry is absent: ", cache_dir, call. = FALSE)
}

read_gmt_long <- function(path, provider = "custom") {
  lines <- readLines(path, warn = FALSE)
  rows <- lapply(lines[nzchar(lines)], function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 4) return(NULL)
    data.frame(term = fields[[1]], description = fields[[2]], gene_symbol = fields[-c(1, 2)], provider = provider)
  })
  bind_rows(rows)
}

sets <- read_gmt_long(args[["custom-gmt"]])
provider_versions <- list(custom = list(source = normalizePath(args[["custom-gmt"]], mustWork = TRUE)))

if (enabled("msigdb")) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("MSigDB provider requires the msigdbr package", call. = FALSE)
  database_species <- if (cfg$species$provider == "mouse") "MM" else "HS"
  msig <- msigdbr::msigdbr(db_species = database_species, species = cfg$species$scientific_name)
  collection_key <- ifelse(is.na(msig$gs_subcollection) | msig$gs_subcollection == "", msig$gs_collection, paste(msig$gs_collection, msig$gs_subcollection, sep = ":"))
  requested <- unlist(cfg$resources$gene_sets$collections)
  if (cfg$species$provider == "mouse") {
    requested <- sub("^H$", "MH", requested)
    requested <- sub("^C([0-9])", "M\\1", requested)
  }
  if (length(requested)) msig <- msig[collection_key %in% requested, , drop = FALSE]
  sets <- bind_rows(sets, data.frame(term = msig$gs_name, description = msig$gs_description, gene_symbol = msig$gene_symbol, provider = "msigdb"))
  provider_versions$msigdb <- list(package = as.character(utils::packageVersion("msigdbr")), database_release = unique(msig$db_version))
}

orgdb <- NULL
orgdb_package <- NULL
if (cfg$species$provider == "mouse" && requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  orgdb <- get("org.Mm.eg.db", asNamespace("org.Mm.eg.db")); orgdb_package <- "org.Mm.eg.db"
}
if (cfg$species$provider == "human" && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  orgdb <- get("org.Hs.eg.db", asNamespace("org.Hs.eg.db")); orgdb_package <- "org.Hs.eg.db"
}

if (enabled("go")) {
  if (is.null(orgdb)) stop("GO provider requires the species-matched org.*.eg.db package", call. = FALSE)
  go <- AnnotationDbi::select(orgdb, keys = AnnotationDbi::keys(orgdb, keytype = "ENTREZID"), columns = c("SYMBOL", "GOALL", "ONTOLOGYALL"), keytype = "ENTREZID")
  go <- go[!is.na(go$SYMBOL) & !is.na(go$GOALL) & go$ONTOLOGYALL %in% go_domains, , drop = FALSE]
  go_terms <- if (requireNamespace("GO.db", quietly = TRUE)) AnnotationDbi::mapIds(get("GO.db", asNamespace("GO.db")), keys = unique(go$GOALL), column = "TERM", keytype = "GOID", multiVals = "first") else setNames(unique(go$GOALL), unique(go$GOALL))
  # Term id carries the domain (GO_BP_/GO_CC_/GO_MF_) so ontology.R can facet by
  # domain; GO_BP_ ids stay byte-identical to the BP-only default.
  sets <- bind_rows(sets, data.frame(term = paste0("GO_", go$ONTOLOGYALL, "_", go$GOALL), description = unname(go_terms[go$GOALL]), gene_symbol = go$SYMBOL, provider = "go"))
  provider_versions$go <- list(orgdb = orgdb_package, package = as.character(utils::packageVersion(orgdb_package)), domains = sort(go_domains))
}

if (enabled("kegg")) {
  if (is.null(orgdb)) stop("KEGG provider requires the species-matched org.*.eg.db package", call. = FALSE)
  organism <- if (cfg$species$provider == "mouse") "mmu" else if (cfg$species$provider == "human") "hsa" else stop("KEGG provider supports mouse or human", call. = FALSE)
  links <- KEGGREST::keggLink("pathway", organism)
  link_table <- data.frame(entrez = sub(paste0("^", organism, ":"), "", names(links)), pathway = sub("^path:", "", unname(links)))
  symbols <- AnnotationDbi::mapIds(orgdb, keys = unique(link_table$entrez), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  names_table <- KEGGREST::keggList("pathway", organism)
  names(names_table) <- sub("^path:", "", names(names_table))
  link_table$gene_symbol <- unname(symbols[link_table$entrez])
  link_table$description <- unname(names_table[link_table$pathway])
  link_table <- link_table[!is.na(link_table$gene_symbol), , drop = FALSE]
  sets <- bind_rows(sets, data.frame(term = paste0("KEGG_", link_table$pathway), description = link_table$description, gene_symbol = link_table$gene_symbol, provider = "kegg"))
  provider_versions$kegg <- list(package = as.character(utils::packageVersion("KEGGREST")), organism = organism, retrieval = "live KEGG REST")
}

# Reactome curated pathways, sourced from MSigDB's C2:CP:REACTOME subcollection
# (M2:CP:REACTOME for mouse) so the identifiers and gene symbols match the rest
# of the MSigDB-backed resources. Opt-in breadth alongside KEGG.
if (enabled("reactome")) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Reactome provider requires the msigdbr package", call. = FALSE)
  database_species <- if (cfg$species$provider == "mouse") "MM" else "HS"
  reactome_key <- if (cfg$species$provider == "mouse") "M2:CP:REACTOME" else "C2:CP:REACTOME"
  react <- msigdbr::msigdbr(db_species = database_species, species = cfg$species$scientific_name)
  react_collection_key <- ifelse(is.na(react$gs_subcollection) | react$gs_subcollection == "", react$gs_collection, paste(react$gs_collection, react$gs_subcollection, sep = ":"))
  react <- react[react_collection_key %in% reactome_key, , drop = FALSE]
  if (nrow(react)) sets <- bind_rows(sets, data.frame(term = react$gs_name, description = react$gs_description, gene_symbol = react$gene_symbol, provider = "reactome"))
  provider_versions$reactome <- list(package = as.character(utils::packageVersion("msigdbr")), subcollection = reactome_key, database_release = unique(react$db_version))
}

if (enabled("string")) provider_versions$string <- list(api = "STRING", taxonomy_id = cfg$species$taxonomy_id, status = "resolved by the network module")
if (enabled("dorothea")) provider_versions$dorothea <- list(dataset = if (cfg$species$provider == "mouse") "dorothea_mm" else "dorothea_hs", status = "resolved by the regulator module")
if (enabled("gtrd")) provider_versions$gtrd <- list(
  source = normalizePath(gtrd_snapshot, mustWork = TRUE),
  sha256 = digest::digest(file = gtrd_snapshot, algo = "sha256"),
  status = "user-supplied exported snapshot; interpreted as unsigned binding evidence"
)

sets <- sets %>% filter(!is.na(term), !is.na(gene_symbol), term != "", gene_symbol != "") %>% distinct(term, gene_symbol, .keep_all = TRUE) %>% arrange(term, gene_symbol)
if (!nrow(sets)) stop("Resource resolution produced no gene sets", call. = FALSE)
readr::write_tsv(sets, args$table)
gmt_lines <- lapply(split(sets, sets$term), function(group) paste(c(group$term[[1]], group$description[[1]], sort(unique(group$gene_symbol))), collapse = "\t"))
writeLines(unname(unlist(gmt_lines)), args$gmt)

receipt <- list(
  schema_version = 1,
  provider = "Tifzoret gene-set resolver",
  organism = cfg$species,
  database_release = provider_versions,
  retrieval_time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  request_parameters = list(collections = cfg$resources$gene_sets$collections, providers = providers),
  license_notice = "Individual provider terms apply; users are responsible for compliance with upstream licenses.",
  gene_sets = length(unique(sets$term)),
  mappings = nrow(sets),
  sha256 = digest::digest(file = args$gmt, algo = "sha256")
)
write_json_file(receipt, args$receipt)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(args$gmt, cache_gmt, overwrite = TRUE)
file.copy(args$table, cache_table, overwrite = TRUE)
file.copy(args$receipt, cache_receipt, overwrite = TRUE)
