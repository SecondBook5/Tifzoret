# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     estimands.R (shared helper, not a stage)
# │ WHAT:      Cell weights -> model-coefficient contrast vector
# │ WHY:       The single authority on what an estimand means numerically, so
# │            DESeq2 and edgeR can never disagree on a contrast's sign
# │ HOW:       Evaluate model.matrix at one synthetic row per referenced cell
# │            with non-cell variables held fixed, then take the weighted sum
# │ INPUTS:    inputs/{families,estimands,estimand_cell_weights}.tsv
# │ PRODUCES:  nothing on disk; sourced by the 03_differential stages
# │ CALLED BY: family_fit.R, estimand.R, term_test.R, de_confirm.R
# │ ENV:       workflow/envs/r.yaml
# └─────────────────────────────────────────────────────────────────

# Expressions are parsed ONCE, in Python (config/estimands.py), and reach this
# script as cell weights. Nothing here parses expression syntax; nothing in
# Python interprets model.matrix coding. That split is deliberate: a second
# implementation of coefficient semantics is how silent sign errors are made.
#
# The mapping is c = sum_k w_k * x_k, where x_k is the model-matrix row at cell
# k with every non-cell variable HELD at a fixed value. Because the weights sum
# to zero, those held columns cancel exactly -- which is why sum-to-zero is a
# hard requirement upstream. Averaging the real model-matrix rows of each cell's
# samples instead would, under an unbalanced covariate, leave nonzero weight on
# that covariate and silently fold a nuisance difference into the estimand.

CELL_KEY_SEPARATOR <- "|"

# Which samples the heatmap/PCA display. de.R branched on contrast type: a
# pairwise contrast showed its two factor levels; a coefficient contrast showed
# everything. An estimand generalizes both -- it shows exactly the samples in
# the cells its expression references, which reproduces the pairwise behaviour
# for a two-cell estimand and the all-groups behaviour for an interaction.
display_samples_for_cells <- function(metadata, cells, cell_keys) {
  if (!length(cells) || !length(cell_keys)) return(rep(TRUE, nrow(metadata)))
  keys <- apply(metadata[, cells, drop = FALSE], 1,
                function(row) paste(row, collapse = CELL_KEY_SEPARATOR))
  referenced <- keys %in% cell_keys
  if (!any(referenced)) return(rep(TRUE, nrow(metadata)))
  referenced
}

read_tsv_plain <- function(path) {
  utils::read.delim(
    normalizePath(path, mustWork = TRUE),
    stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character"
  )
}

read_family_spec <- function(families_path, family_id) {
  table <- read_tsv_plain(families_path)
  row <- table[table$family_id == family_id, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("could not resolve exactly one family: ", family_id, call. = FALSE)
  }
  cells <- strsplit(row$cells[[1]], CELL_KEY_SEPARATOR, fixed = TRUE)[[1]]
  cells <- cells[nzchar(cells)]
  references <- list()
  raw <- trimws(row$reference_levels[[1]])
  if (nzchar(raw)) {
    for (piece in strsplit(raw, ";", fixed = TRUE)[[1]]) {
      kv <- strsplit(piece, "=", fixed = TRUE)[[1]]
      if (length(kv) == 2L) references[[trimws(kv[[1]])]] <- trimws(kv[[2]])
    }
  }
  list(
    family_id = family_id,
    design = row$design[[1]],
    cells = cells,
    reference_levels = references,
    replicate_unit = if (nzchar(row$replicate_unit[[1]])) row$replicate_unit[[1]] else NULL,
    shrinkage = row$shrinkage[[1]],
    filter = row$filter[[1]],
    declared = identical(row$declared[[1]], "true")
  )
}

read_estimand_spec <- function(estimands_path, estimand_id) {
  table <- read_tsv_plain(estimands_path)
  row <- table[table$estimand_id == estimand_id, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("could not resolve exactly one estimand: ", estimand_id, call. = FALSE)
  }
  list(
    family_id = row$family_id[[1]],
    estimand_id = estimand_id,
    label = row$label[[1]],
    role = row$role[[1]],
    expression = row$expression[[1]],
    atom_kind = row$atom_kind[[1]]
  )
}

read_estimand_ids <- function(estimands_path, family_id) {
  table <- read_tsv_plain(estimands_path)
  table$estimand_id[table$family_id == family_id]
}

read_cell_weights <- function(weights_path, estimand_id) {
  table <- read_tsv_plain(weights_path)
  rows <- table[table$estimand_id == estimand_id, , drop = FALSE]
  if (!nrow(rows)) stop("no cell weights for estimand: ", estimand_id, call. = FALSE)
  stats::setNames(as.numeric(rows$weight), rows$cell_key)
}

read_term_tests <- function(term_tests_path, family_id) {
  table <- read_tsv_plain(term_tests_path)
  table[table$family_id == family_id, , drop = FALSE]
}

# One synthetic row per referenced cell. Cell variables take that cell's levels;
# every other design variable is held at a fixed value (first level / zero).
# Factor levels are taken from the FITTING metadata so model.matrix produces the
# same columns it produced for the fit.
cell_newdata <- function(cells, cell_keys, metadata, design_formula) {
  variables <- all.vars(design_formula)
  frame <- data.frame(row.names = seq_along(cell_keys))
  parts <- strsplit(cell_keys, CELL_KEY_SEPARATOR, fixed = TRUE)
  # Arity check: every cell key must split into exactly length(cells) pieces.
  # Python should already reject malformed keys, but the sign authority defends itself.
  arities <- vapply(parts, length, integer(1))
  expected <- length(cells)
  if (any(arities != expected)) {
    bad_indices <- which(arities != expected)
    bad_keys <- cell_keys[bad_indices]
    bad_arities <- arities[bad_indices]
    stop(
      "cell key arity mismatch: expected ", expected, " level(s) but got ",
      paste(bad_arities, collapse = ", "), " for key(s): ",
      paste(bad_keys, collapse = ", "), call. = FALSE
    )
  }
  for (index in seq_along(cells)) {
    column <- cells[[index]]
    values <- vapply(parts, function(piece) piece[[index]], character(1))
    reference <- metadata[[column]]
    if (is.null(reference)) stop("cell column absent from metadata: ", column, call. = FALSE)
    known <- levels(factor(reference))
    unknown <- setdiff(values, known)
    if (length(unknown)) {
      stop("unknown level(s) for ", column, ": ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    frame[[column]] <- factor(values, levels = levels(factor(reference)))
  }
  for (column in setdiff(variables, cells)) {
    reference <- metadata[[column]]
    if (is.null(reference)) next
    frame[[column]] <- if (is.numeric(reference)) {
      rep(0, length(cell_keys))
    } else {
      factor(rep(levels(factor(reference))[[1]], length(cell_keys)),
             levels = levels(factor(reference)))
    }
  }
  frame
}

# Align a synthetic model matrix to DESeq2's resultsNames() POSITIONALLY.
#
# DESeq2 renames main-effect columns (factor_aa2 -> factor_a_a2_vs_a1) but never
# reorders them, so model.matrix(design, newdata) yields columns in the same order
# as resultsNames(dds). A name-based alignment would abort on every family with a
# main effect because make.names("factor_aa2") != "factor_a_a2_vs_a1".
#
# For interactions, make.names already agrees (factor_aa2:factor_bb2 ->
# factor_aa2.factor_bb2 = DESeq2's naming), so we sanity-check those. A mismatch
# there means the positional assumption has broken and we must stop.
#
# A future reader "fixing" this to match by name would break every main-effect family.
align_to_coefficients <- function(model_matrix, coefficient_names) {
  if (ncol(model_matrix) != length(coefficient_names)) {
    stop(
      "model matrix column count does not match coefficient count: ",
      "matrix has ", ncol(model_matrix), " columns [",
      paste(colnames(model_matrix), collapse = ", "),
      "]; coefficients has ", length(coefficient_names), " names [",
      paste(coefficient_names, collapse = ", "), "]",
      call. = FALSE
    )
  }
  munged <- make.names(colnames(model_matrix))
  # Sanity check: where make.names already agrees (interaction case), confirm
  for (i in seq_along(munged)) {
    if (munged[[i]] == coefficient_names[[i]]) {
      # Already matches — positional assumption holding
    } else if (munged[[i]] %in% coefficient_names) {
      # Name exists but in wrong position — positional assumption broken
      stop(
        "model matrix column order mismatch: position ", i,
        " has ", colnames(model_matrix)[[i]], " (make.names: ", munged[[i]],
        ") but coefficients has ", coefficient_names[[i]],
        "; positional alignment assumption violated",
        call. = FALSE
      )
    }
  }
  # Align by position, name with coefficient_names
  result <- model_matrix[, seq_along(coefficient_names), drop = FALSE]
  colnames(result) <- coefficient_names
  result
}

compile_contrast_vector <- function(spec, weights, metadata, design_formula,
                                    coefficient_names) {
  if (identical(spec$atom_kind, "coef")) {
    unknown <- setdiff(names(weights), coefficient_names)
    if (length(unknown)) {
      stop(
        "estimand ", spec$estimand_id, ": unknown coefficient(s) ",
        paste(unknown, collapse = ", "), "; available: ",
        paste(coefficient_names, collapse = ", "), call. = FALSE
      )
    }
    vector <- stats::setNames(rep(0, length(coefficient_names)), coefficient_names)
    vector[names(weights)] <- as.numeric(weights)
    return(vector)
  }
  if (abs(sum(weights)) > 1e-9) {
    stop("estimand ", spec$estimand_id, ": cell weights do not sum to zero", call. = FALSE)
  }
  newdata <- cell_newdata(spec$cells, names(weights), metadata, design_formula)
  model_matrix <- stats::model.matrix(design_formula, data = newdata)
  model_matrix <- align_to_coefficients(model_matrix, coefficient_names)
  vector <- as.numeric(t(model_matrix) %*% as.numeric(weights))
  vector[abs(vector) < 1e-12] <- 0
  stats::setNames(vector, coefficient_names)
}
