"""Design validation: cell membership, rank checks, confounding detection."""

from __future__ import annotations

import re
from collections.abc import Sequence

import numpy as np

from .parser import CELL_KEY_SEPARATOR
from .types import Family, TermTest


def _observed_levels(sample_rows: Sequence[dict[str, str]], column: str) -> list[str]:
    """Extract unique observed levels for a column, sorted."""
    return sorted({str(row[column]) for row in sample_rows if column in row})


def cell_membership(family: Family, sample_rows: Sequence[dict[str, str]]) -> dict[tuple[str, ...], list[str]]:
    """Map each observed design cell to the sample ids it contains."""
    membership: dict[tuple[str, ...], list[str]] = {}
    for row in sample_rows:
        if any(column not in row for column in family.cells):
            continue
        key = tuple(str(row[column]) for column in family.cells)
        membership.setdefault(key, []).append(str(row.get("sample_id", "")))
    return membership


def cell_mean_matrix(family: Family, sample_rows: Sequence[dict[str, str]]) -> tuple[np.ndarray, list[tuple[str, ...]]]:
    """Indicator matrix of samples against observed cells, plus the cell order."""
    membership = cell_membership(family, sample_rows)
    cells = sorted(membership)
    index = {cell: position for position, cell in enumerate(cells)}
    matrix = np.zeros((len(sample_rows), len(cells)), dtype=float)
    for position, row in enumerate(sample_rows):
        if any(column not in row for column in family.cells):
            continue
        key = tuple(str(row[column]) for column in family.cells)
        matrix[position, index[key]] = 1.0
    return matrix, cells


def validate_family_design(family: Family, sample_rows: Sequence[dict[str, str]]) -> list[str]:
    """Every hard fail from spec §6 that needs only the sample table."""
    # Import here to avoid circular dependency
    from .compiler import design_variables

    errors: list[str] = []
    columns = set().union(*(row.keys() for row in sample_rows)) if sample_rows else set()

    for column in family.cells:
        if column not in columns:
            errors.append(
                f"family {family.id}: cell column {column!r} is absent from samples.tsv"
            )
    for variable in design_variables(family.design):
        if variable not in columns:
            errors.append(
                f"family {family.id}: design variable {variable!r} is absent from samples.tsv"
            )
    if errors:
        return errors

    membership = cell_membership(family, sample_rows)
    for estimand in family.estimands:
        for cell in estimand.expression.cell_weights:
            # Check for CELL_KEY_SEPARATOR in any level before checking membership,
            # because a level containing the separator would be mis-parsed in R
            for level in cell:
                if CELL_KEY_SEPARATOR in level:
                    errors.append(
                        f"estimand {estimand.id}: level {level!r} in cell "
                        f"({CELL_KEY_SEPARATOR.join(cell)}) contains the reserved "
                        f"separator {CELL_KEY_SEPARATOR!r}; the level must be renamed"
                    )
            if cell not in membership:
                observed = ", ".join(
                    CELL_KEY_SEPARATOR.join(key) for key in sorted(membership)
                )
                errors.append(
                    f"estimand {estimand.id}: cell "
                    f"({CELL_KEY_SEPARATOR.join(cell)}) is empty or references an "
                    f"unknown level; observed cells: {observed}"
                )

    matrix, cells = cell_mean_matrix(family, sample_rows)
    if matrix.size:
        rank = int(np.linalg.matrix_rank(matrix))
        if rank < matrix.shape[1]:
            errors.append(
                f"family {family.id}: cell-mean design is rank-deficient "
                f"(rank {rank} < {matrix.shape[1]} cells)"
            )
        # Skip confounding validation for coefficient-only families (cells=()).
        # With a single degenerate cell containing all samples, both confounding
        # branches are vacuous: every nuisance level trivially occurs in the only
        # cell, so the check cannot say anything meaningful.
        if family.cells:
            nuisance = [
                variable
                for variable in design_variables(family.design)
                if variable not in set(family.cells)
            ]
            cell_index = {cell: position for position, cell in enumerate(cells)}
            for variable in nuisance:
                per_cell: dict[tuple[str, ...], set[str]] = {}
                per_value: dict[str, set[tuple[str, ...]]] = {}
                for row in sample_rows:
                    if variable not in row:
                        continue
                    key = tuple(str(row[column]) for column in family.cells)
                    if key not in cell_index:
                        continue
                    value = str(row[variable])
                    per_cell.setdefault(key, set()).add(value)
                    per_value.setdefault(value, set()).add(key)
                constant_within_cells = per_cell and all(
                    len(values) == 1 for values in per_cell.values()
                )
                determines_cell = per_value and all(
                    len(keys) == 1 for keys in per_value.values()
                )
                if constant_within_cells and len(per_cell) > 1:
                    errors.append(
                        f"family {family.id}: nuisance covariate {variable!r} is perfectly "
                        "confounded with cell membership (constant within every cell)"
                    )
                elif determines_cell and len(per_value) > 1:
                    errors.append(
                        f"family {family.id}: nuisance covariate {variable!r} is perfectly "
                        "confounded with cell membership (each level occurs in one cell)"
                    )

        # This parameter count is deliberately conservative: it counts cell-mean parameters
        # only (product of observed level counts per factor), ignoring nuisance covariates,
        # so it undercounts the true rank. This means a marginally non-estimable design
        # can pass this check and fail later at fit time, which is acceptable because the
        # real rank, condition number, and residual df are computed from the actual model
        # matrix by a later stage, which warns appropriately. Erring toward permissive
        # here means an analyst gets feedback at fit time rather than rejecting a design
        # that might be estimable despite marginal rank.
        parameters = 1
        for column in family.cells:
            parameters *= max(1, len(_observed_levels(sample_rows, column)))
        residual = len(sample_rows) - parameters
        if residual < 1:
            errors.append(
                f"family {family.id}: residual degrees of freedom is {residual} "
                f"({len(sample_rows)} samples, {parameters} model parameters); "
                "the model is not estimable"
            )

    if family.replicate_unit:
        if family.replicate_unit not in columns:
            errors.append(
                f"family {family.id}: replicate_unit {family.replicate_unit!r} is "
                "absent from samples.tsv"
            )
        else:
            per_cell_units: dict[tuple[str, ...], list[str]] = {}
            for row in sample_rows:
                key = tuple(str(row[column]) for column in family.cells)
                per_cell_units.setdefault(key, []).append(str(row[family.replicate_unit]))
            for key, units in per_cell_units.items():
                duplicates = sorted({unit for unit in units if units.count(unit) > 1})
                if duplicates:
                    errors.append(
                        f"family {family.id}: replicate_unit value(s) "
                        f"{', '.join(duplicates)} appear more than once inside cell "
                        f"({CELL_KEY_SEPARATOR.join(key)}). The engine has no "
                        "mixed-model path, so these repeated measures cannot be "
                        "treated as independent; aggregate them or remove "
                        "replicate_unit deliberately."
                    )
    return errors


def resolve_term_test_df(family: Family, sample_rows: Sequence[dict[str, str]]) -> tuple[TermTest, ...]:
    """Recompute each derived term test's df from observed level counts."""
    from .compiler import parse_formula_terms

    levels = {
        column: max(1, len(_observed_levels(sample_rows, column)))
        for column in family.cells
    }
    cell_set = set(family.cells)
    all_terms = parse_formula_terms(family.design)
    resolved: list[TermTest] = []
    for test in family.term_tests:
        kept = set(parse_formula_terms(test.reduced))
        dropped = [
            term for term in all_terms if term <= cell_set and term and term not in kept
        ]
        df = 0
        for term in dropped:
            contribution = 1
            for factor in term:
                contribution *= levels.get(factor, 2) - 1
            df += contribution
        resolved.append(
            TermTest(id=test.id, family_id=test.family_id, reduced=test.reduced, df=df)
        )
    return tuple(resolved)


def family_slug(design: str) -> str:
    """A deterministic, filesystem-safe family id derived from a design formula."""
    body = design.split("~", 1)[-1].lower()
    slug = re.sub(r"[^a-z0-9]+", "_", body).strip("_")
    return f"design_{slug}" if slug else "design_intercept"


def _quote_level(level: str) -> str:
    """Quote a level if it contains special characters that require escaping.

    Uses double quotes if the level contains a single quote, single quotes if it
    contains a double quote, and single quotes for any other special characters.
    """
    if not re.search(r"[,()'\"]", level):
        return level
    # Prefer single quotes, but use double if single quote is in the level
    if "'" not in level:
        return f"'{level}'"
    else:
        return f'"{level}"'
