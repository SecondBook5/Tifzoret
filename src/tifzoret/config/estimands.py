"""Declared estimands: expression parsing, compilation, and design validation.

An estimand is any linear combination of fitted model coefficients. Studies
declare one as a signed sum of design *cells* — ordered level-tuples over the
family's crossed factor columns — which is parameterization-free: unlike a
coefficient name, a cell tuple means the same thing however the model is
parameterized. This module turns that text into a cell-weight vector. Mapping
those weights onto actual model-matrix columns happens in R
(``workflow/scripts/estimands.R``), which is the only place ``model.matrix``
semantics are interpreted, so coefficient coding is never implemented twice.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from itertools import combinations

import numpy as np


class EstimandSyntaxError(ValueError):
    """Raised when an estimand expression cannot be parsed."""


# One atom plus its leading sign and optional numeric coefficient. Levels are
# either bare words or quoted (so a level may itself contain a comma or paren).
_LEVEL = r"""(?:'[^']*'|"[^"]*"|[^,()'"]+)"""
_CELL = rf"\(\s*{_LEVEL}(?:\s*,\s*{_LEVEL})*\s*\)"
_COEF = r"coef\(\s*[^)]+\s*\)"
_TERM = re.compile(
    rf"^\s*(?:(?P<number>[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*\*\s*)?"
    rf"(?P<atom>{_CELL}|{_COEF})\s*$"
)


@dataclass(frozen=True)
class ParsedExpression:
    """A parsed estimand: either cell weights or coefficient weights, never both."""

    source: str
    cell_weights: dict[tuple[str, ...], float] = field(default_factory=dict)
    coef_weights: dict[str, float] = field(default_factory=dict)

    @property
    def atom_kind(self) -> str:
        """``"coef"`` for a ``coef(...)`` expression, ``"cell"`` otherwise."""
        return "coef" if self.coef_weights else "cell"

    @property
    def weight_sum(self) -> float:
        """Sum of the signed cell weights; 0 for a well-formed contrast."""
        return float(sum(self.cell_weights.values()))


def _unquote(level: str) -> str:
    level = level.strip()
    if len(level) >= 2 and level[0] == level[-1] and level[0] in "'\"":
        return level[1:-1]
    return level


def _split_terms(text: str) -> list[tuple[float, str]]:
    """Split an expression into (sign, atom-text) pairs, preserving order."""
    padded = text.strip()
    if not padded:
        raise EstimandSyntaxError("estimand expression is empty")
    # Normalize a leading sign so the split below always yields sign/atom pairs.
    if padded[0] not in "+-":
        padded = "+" + padded
    pieces: list[tuple[float, str]] = []
    index = 0
    while index < len(padded):
        sign_char = padded[index]
        if sign_char not in "+-":
            raise EstimandSyntaxError(f"expected '+' or '-' at: {padded[index:]!r}")
        sign = 1.0 if sign_char == "+" else -1.0
        index += 1
        depth = 0
        start = index
        while index < len(padded):
            char = padded[index]
            if char in "'\"":
                closing = padded.find(char, index + 1)
                if closing == -1:
                    raise EstimandSyntaxError(f"unterminated quote in: {text!r}")
                index = closing + 1
                continue
            if char == "(":
                depth += 1
            elif char == ")":
                depth -= 1
                if depth < 0:
                    raise EstimandSyntaxError(f"unbalanced parentheses in: {text!r}")
            elif char in "+-" and depth == 0:
                break
            index += 1
        if depth != 0:
            raise EstimandSyntaxError(f"unbalanced parentheses in: {text!r}")
        atom = padded[start:index].strip()
        if not atom:
            raise EstimandSyntaxError(f"missing term after '{sign_char}' in: {text!r}")
        pieces.append((sign, atom))
    return pieces


def parse_expression(text: str, arity: int) -> ParsedExpression:
    """Parse an estimand expression into cell weights or coefficient weights.

    ``arity`` is ``len(family.cells)``; every cell tuple must match it.
    """
    source = text.strip()
    cell_weights: dict[tuple[str, ...], float] = {}
    coef_weights: dict[str, float] = {}
    for sign, atom_text in _split_terms(source):
        match = _TERM.match(atom_text)
        if match is None:
            raise EstimandSyntaxError(f"cannot parse estimand term: {atom_text!r}")
        magnitude = float(match.group("number")) if match.group("number") else 1.0
        weight = sign * magnitude
        atom = match.group("atom")
        if atom.startswith("coef"):
            name = atom[atom.index("(") + 1 : atom.rindex(")")].strip()
            if not name:
                raise EstimandSyntaxError("coef() requires a coefficient name")
            coef_weights[name] = coef_weights.get(name, 0.0) + weight
        else:
            inner = atom[atom.index("(") + 1 : atom.rindex(")")]
            levels = tuple(_unquote(part) for part in _split_levels(inner))
            if len(levels) != arity:
                raise EstimandSyntaxError(
                    f"cell {atom!r} has arity {len(levels)}, expected {arity}"
                )
            if any(not level for level in levels):
                raise EstimandSyntaxError(f"cell {atom!r} has an empty level")
            cell_weights[levels] = cell_weights.get(levels, 0.0) + weight
    if cell_weights and coef_weights:
        raise EstimandSyntaxError(
            "estimand expression mixed cell and coef() atoms; use one form"
        )
    if not cell_weights and not coef_weights:
        raise EstimandSyntaxError("estimand expression is empty")
    return ParsedExpression(
        source=source, cell_weights=cell_weights, coef_weights=coef_weights
    )


def _split_levels(inner: str) -> list[str]:
    """Split a cell's interior on commas that are not inside quotes."""
    levels: list[str] = []
    current: list[str] = []
    index = 0
    while index < len(inner):
        char = inner[index]
        if char in "'\"":
            closing = inner.find(char, index + 1)
            if closing == -1:
                raise EstimandSyntaxError(f"unterminated quote in cell: {inner!r}")
            current.append(inner[index : closing + 1])
            index = closing + 1
            continue
        if char == ",":
            levels.append("".join(current))
            current = []
        else:
            current.append(char)
        index += 1
    levels.append("".join(current))
    return levels


SHRINKAGE_CHOICES = ("apeglm", "ashr", "none")
FILTER_CHOICES = ("design_aware", "total_count", "none")
WEIGHT_SUM_TOLERANCE = 1e-9
CELL_KEY_SEPARATOR = "|"


@dataclass(frozen=True)
class Estimand:
    """One row of a family's contrast matrix."""

    id: str
    family_id: str
    label: str
    role: str
    expression: ParsedExpression


@dataclass(frozen=True)
class TermTest:
    """One full-vs-reduced likelihood-ratio comparison within a family."""

    id: str
    family_id: str
    reduced: str
    df: int


@dataclass(frozen=True)
class Family:
    """A design plus the crossed factor columns: the unit of fitting."""

    id: str
    design: str
    cells: tuple[str, ...]
    reference_levels: dict[str, str]
    replicate_unit: str | None
    shrinkage: str
    filter: str
    declared: bool
    estimands: tuple[Estimand, ...]
    term_tests: tuple[TermTest, ...] = ()


def design_variables(design: str) -> list[str]:
    """The variable names appearing in an R-style one-sided formula."""
    body = design.split("~", 1)[-1]
    tokens = re.findall(r"[A-Za-z_.][A-Za-z0-9_.]*", body)
    return [token for token in dict.fromkeys(tokens) if token not in {"I", "log", "log2"}]


def build_families(
    families_config: dict, *, known_ids: set[str]
) -> tuple[list[Family], list[str]]:
    """Build validated ``Family`` objects from ``analysis.families``.

    Performs only the checks that need no sample table: syntax, arity,
    sum-to-zero, knob choices, primary-role count, and id collisions. The
    sample-dependent gate lives in :func:`validate_family_design`.
    """
    errors: list[str] = []
    families: list[Family] = []
    seen: set[str] = set(known_ids)
    for family_id, raw in (families_config or {}).items():
        design = str(raw.get("design", "")).strip()
        cells = tuple(str(column) for column in (raw.get("cells") or ()))
        variables = set(design_variables(design))
        if cells:
            for column in cells:
                if column not in variables:
                    errors.append(
                        f"family {family_id}: cell column {column!r} is absent from design {design!r}"
                    )
        shrinkage = str(raw.get("shrinkage", SHRINKAGE_CHOICES[0]))
        if shrinkage not in SHRINKAGE_CHOICES:
            errors.append(
                f"family {family_id}: shrinkage {shrinkage!r} must be one of "
                f"{', '.join(SHRINKAGE_CHOICES)}"
            )
        gene_filter = str(raw.get("filter", FILTER_CHOICES[0]))
        if gene_filter not in FILTER_CHOICES:
            errors.append(
                f"family {family_id}: filter {gene_filter!r} must be one of "
                f"{', '.join(FILTER_CHOICES)}"
            )
        raw_estimands = raw.get("estimands") or []
        declared = bool(raw.get("_declared", True))
        if declared and not raw_estimands:
            errors.append(f"family {family_id}: requires at least one estimand")
        estimands: list[Estimand] = []
        primaries = 0
        for entry in raw_estimands:
            estimand_id = str(entry.get("id", "")).strip()
            if not estimand_id:
                errors.append(f"family {family_id}: an estimand is missing its id")
                continue
            if estimand_id in seen:
                errors.append(
                    f"family {family_id}: estimand id {estimand_id!r} collides with "
                    "another estimand or contrast id"
                )
            seen.add(estimand_id)
            role = str(entry.get("role", "")).strip()
            if role == "primary":
                primaries += 1
            elif role:
                errors.append(
                    f"estimand {estimand_id}: role {role!r} must be 'primary' or absent"
                )
            try:
                parsed = parse_expression(str(entry.get("expression", "")), len(cells))
            except EstimandSyntaxError as error:
                errors.append(f"estimand {estimand_id}: {error}")
                continue
            if (
                parsed.atom_kind == "cell"
                and abs(parsed.weight_sum) > WEIGHT_SUM_TOLERANCE
            ):
                errors.append(
                    f"estimand {estimand_id}: cell weights must sum to zero, got "
                    f"{parsed.weight_sum:g}"
                )
            estimands.append(
                Estimand(
                    id=estimand_id,
                    family_id=family_id,
                    label=str(entry.get("label", estimand_id)),
                    role=role,
                    expression=parsed,
                )
            )
        if declared and raw_estimands and primaries != 1:
            errors.append(
                f"family {family_id}: a declared family requires exactly one "
                f"role: primary estimand, found {primaries}"
            )
        reference_levels = {
            str(key): str(value)
            for key, value in (raw.get("reference_levels") or {}).items()
        }
        replicate_unit = raw.get("replicate_unit")
        families.append(
            Family(
                id=family_id,
                design=design,
                cells=cells,
                reference_levels=reference_levels,
                replicate_unit=str(replicate_unit) if replicate_unit else None,
                shrinkage=shrinkage,
                filter=gene_filter,
                declared=declared,
                estimands=tuple(estimands),
                term_tests=(),
            )
        )
    return families, errors


TERM_TEST_CEILING = 12


def parse_formula_terms(design: str) -> list[frozenset[str]]:
    """Expand a one-sided R formula into its model terms.

    ``a * b`` becomes ``a``, ``b``, ``a:b``; ``a:b`` stays a single term. Order
    of first appearance is preserved so a rendered formula reads like the
    declared one.
    """
    body = design.split("~", 1)[-1]
    terms: list[frozenset[str]] = []
    for piece in body.split("+"):
        piece = piece.strip()
        if not piece or piece == "1":
            continue
        if "*" in piece:
            factors = [part.strip() for part in piece.split("*") if part.strip()]
            for size in range(1, len(factors) + 1):
                for combo in combinations(factors, size):
                    candidate = frozenset(combo)
                    if candidate not in terms:
                        terms.append(candidate)
        else:
            factors = [part.strip() for part in piece.split(":") if part.strip()]
            candidate = frozenset(factors)
            if candidate and candidate not in terms:
                terms.append(candidate)
    return terms


def _term_label(term: frozenset[str], order: list[str]) -> str:
    return ":".join(sorted(term, key=order.index))


def format_formula(terms: Sequence[frozenset[str]]) -> str:
    """Render model terms back into a one-sided formula string."""
    materialized = list(terms)
    if not materialized:
        return "~ 1"
    order: list[str] = []
    for term in materialized:
        for factor in sorted(term):
            if factor not in order:
                order.append(factor)
    ordered = sorted(materialized, key=lambda term: (len(term), sorted(term)))
    return "~ " + " + ".join(_term_label(term, order) for term in ordered)


def derive_term_tests(design: str, cells: Sequence[str], family_id: str) -> list[TermTest]:
    """Derive the nested full-vs-reduced lattice over the family's cell factors.

    Two generators, per spec §7.4:

    * for each non-empty subset ``S`` of cell factors, drop every cell term that
      intersects ``S`` (so ``{factor_b}`` drops both ``factor_b`` and
      ``factor_a:factor_b`` — "does this gene respond to factor_b in *any*
      state?");
    * for each interaction term of order >= 2, drop that term and every
      higher-order term containing it, keeping the main effects.

    Nuisance terms — every term using no cell factor — stay in both models and
    are never tested.

    Note: ``df`` here counts dropped *terms*, which equals dropped coefficients
    only when every cell factor has two levels. Task 4 recomputes ``df`` from
    the actual level counts once ``samples.tsv`` is available; this value is
    a structural placeholder used only for the nesting check.
    """
    cell_set = set(cells)
    order = list(cells)
    all_terms = parse_formula_terms(design)
    cell_terms = [term for term in all_terms if term <= cell_set and term]
    nuisance_terms = [term for term in all_terms if not (term & cell_set)]
    tests: list[TermTest] = []
    seen: set[str] = set()

    def emit(test_id: str, dropped: list[frozenset[str]]) -> None:
        if test_id in seen or not dropped:
            return
        seen.add(test_id)
        kept = [term for term in cell_terms if term not in dropped]
        tests.append(
            TermTest(
                id=test_id,
                family_id=family_id,
                reduced=format_formula(kept + nuisance_terms),
                df=len(dropped),
            )
        )

    for size in range(1, len(order) + 1):
        for combo in combinations(order, size):
            subset = set(combo)
            dropped = [term for term in cell_terms if term & subset]
            test_id = (
                "any_effect"
                if subset == cell_set
                else "_".join(sorted(combo, key=order.index)) + "_any"
            )
            emit(test_id, dropped)

    for term in cell_terms:
        if len(term) < 2:
            continue
        dropped = [other for other in cell_terms if term <= other]
        emit(_term_label(term, order).replace(":", "_") + "_interaction", dropped)

    return tests


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
        # Skip confounding validation for coefficient-only families (empty cells).
        # These families are labels only and the coefficient is checked at DE time.
        if not family.cells:
            return errors
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


def desugar_contrast_rows(
    global_design: str, contrast_rows: Sequence[dict[str, str]]
) -> tuple[list[Family], list[str]]:
    """Translate every ``contrasts.tsv`` row into a family plus an estimand.

    Pairwise rows become a k=1 family whose estimand is ``(numerator) -
    (denominator)``; coefficient rows become a ``coef(...)`` estimand; omnibus
    rows become a term test. Rows sharing a design share one family — and
    therefore one fit — because a cell-mean expression does not depend on the
    reference level, so the per-contrast releveling that forced separate fits is
    no longer needed.
    """
    errors: list[str] = []
    # Group by exact design string, not by family_id, to prevent collisions
    grouped: dict[str, dict] = {}
    design_order: list[str] = []  # Track first appearance order for stable IDs

    for row in contrast_rows:
        contrast_id = str(row.get("contrast_id", "")).strip()
        row_type = (str(row.get("type", "")).strip() or "pairwise").lower()
        design = str(row.get("design", "")).strip() or global_design

        # Track first appearance of each design
        if design not in grouped:
            design_order.append(design)

        # Use exact design string as key
        bucket = grouped.setdefault(
            design,
            {
                "design": design,
                "cells": [],
                "estimands": [],
                "term_tests": [],
                "reference_levels": {},
                "_declared": False,
            },
        )
        factor = str(row.get("factor", "")).strip()
        references: dict[str, str] = {}
        raw_references = str(row.get("reference_levels", "")).strip()
        for piece in (part for part in raw_references.split(";") if part.strip()):
            if "=" not in piece:
                errors.append(
                    f"contrast {contrast_id}: reference_levels entry {piece!r} must be factor=level"
                )
                continue
            key, value = piece.split("=", 1)
            references[key.strip()] = value.strip()
        for key, value in references.items():
            existing = bucket["reference_levels"].get(key)
            if existing is not None and existing != value:
                errors.append(
                    f"family on design {design!r}: conflicting reference_levels for {key!r} "
                    f"({existing!r} vs {value!r}); give one of these contrasts its own design"
                )
            bucket["reference_levels"][key] = value

        if row_type == "omnibus":
            reduced = str(row.get("reduced", "")).strip()
            if not reduced:
                errors.append(f"contrast {contrast_id}: omnibus row requires a reduced formula")
                continue
            bucket["term_tests"].append(
                {"id": contrast_id, "reduced": reduced, "df": 0}
            )
            if factor and factor not in bucket["cells"]:
                bucket["cells"].append(factor)
            continue

        if row_type == "coefficient":
            coefficient = str(row.get("coefficient", "")).strip()
            if not coefficient:
                errors.append(
                    f"contrast {contrast_id}: coefficient row requires a coefficient"
                )
                continue
            bucket["estimands"].append(
                {"id": contrast_id, "label": contrast_id, "expression": f"coef({coefficient})"}
            )
            continue

        numerator = str(row.get("numerator", "")).strip()
        denominator = str(row.get("denominator", "")).strip()
        if not factor or not numerator or not denominator:
            errors.append(
                f"contrast {contrast_id}: pairwise row requires factor, numerator, denominator"
            )
            continue
        if factor not in bucket["cells"]:
            bucket["cells"].append(factor)
        if len(bucket["cells"]) > 1:
            errors.append(
                f"family on design {design!r}: contrasts on different factors ({', '.join(bucket['cells'])}) "
                "cannot share one desugared family; give one of them its own design"
            )
            continue
        bucket["estimands"].append(
            {
                "id": contrast_id,
                "label": f"{numerator} vs {denominator}",
                "expression": (
                    f"({_quote_level(numerator)}) - ({_quote_level(denominator)})"
                ),
            }
        )

    # Second pass: assign family IDs, handling slug collisions
    family_id_map: dict[str, str] = {}  # design -> family_id
    slug_designs: dict[str, list[str]] = {}  # slug -> list of designs that map to it

    # First, collect which designs map to which slugs
    for design in design_order:
        if design != global_design:
            slug = family_slug(design)
            slug_designs.setdefault(slug, []).append(design)

    # Then assign IDs based on order of first appearance
    for design in design_order:
        if design == global_design:
            family_id = "main"
        else:
            slug = family_slug(design)
            designs_with_this_slug = slug_designs[slug]

            if len(designs_with_this_slug) == 1:
                # No collision
                family_id = slug
            else:
                # Collision: append numeric suffix based on position
                position = designs_with_this_slug.index(design) + 1
                if position == 1:
                    family_id = slug
                else:
                    family_id = f"{slug}_{position}"

        family_id_map[design] = family_id

    # Build config with assigned family IDs
    config = {}
    for design in design_order:
        family_id = family_id_map[design]
        bucket = grouped[design]
        config[family_id] = {
            key: value for key, value in bucket.items() if key != "term_tests"
        }

    families, build_errors = build_families(config, known_ids=set())
    errors.extend(build_errors)
    attached: list[Family] = []
    for family in families:
        # Map family.design back to grouped to get term_tests and reference_levels
        design = family.design
        bucket = grouped[design]
        raw_tests = bucket["term_tests"]
        tests = tuple(
            TermTest(
                id=entry["id"],
                family_id=family.id,
                reduced=entry["reduced"],
                df=entry["df"],
            )
            for entry in raw_tests
        )
        attached.append(
            Family(
                id=family.id,
                design=family.design,
                cells=family.cells,
                reference_levels=bucket["reference_levels"],
                replicate_unit=family.replicate_unit,
                shrinkage=family.shrinkage,
                filter=family.filter,
                declared=False,
                estimands=family.estimands,
                term_tests=tests,
            )
        )
    return attached, errors


FAMILY_FIELDS = (
    "family_id", "design", "cells", "reference_levels",
    "replicate_unit", "shrinkage", "filter", "declared",
)
ESTIMAND_FIELDS = ("family_id", "estimand_id", "label", "role", "expression", "atom_kind")
CELL_WEIGHT_FIELDS = ("family_id", "estimand_id", "cell_key", "weight")
TERM_TEST_FIELDS = ("family_id", "term_test_id", "reduced", "df")


def compile_project_families(
    config: dict, sample_rows: Sequence[dict[str, str]], contrast_rows: Sequence[dict[str, str]]
) -> tuple[list[Family], list[str]]:
    """Compile declared families plus the desugared ``contrasts.tsv`` families.

    Declared families are validated first so their ids reserve the shared
    ``contrasts/<id>/`` namespace; the desugared families reuse the contrast ids
    that are already unique by construction.
    """
    errors: list[str] = []
    analysis = config.get("analysis", {}) or {}
    global_design = str(analysis.get("design", "")).strip()
    contrast_ids = {str(row.get("contrast_id", "")).strip() for row in contrast_rows}

    declared, declared_errors = build_families(
        analysis.get("families") or {}, known_ids=contrast_ids
    )
    errors.extend(declared_errors)
    desugared, desugar_errors = desugar_contrast_rows(global_design, contrast_rows)
    errors.extend(desugar_errors)

    resolved: list[Family] = []
    raw_families = analysis.get("families") or {}
    for family in declared + desugared:
        requested = raw_families.get(family.id, {}).get("term_tests", "auto")
        if requested is False:
            tests: tuple[TermTest, ...] = family.term_tests
        elif isinstance(requested, list):
            tests = tuple(
                TermTest(
                    id=str(entry["id"]),
                    family_id=family.id,
                    reduced=str(entry["reduced"]),
                    df=0,
                )
                for entry in requested
            )
        elif family.declared:
            tests = tuple(derive_term_tests(family.design, family.cells, family.id))
        else:
            tests = family.term_tests
        if len(tests) > TERM_TEST_CEILING:
            errors.append(
                f"family {family.id}: term_tests derived {len(tests)} tests, above the "
                f"ceiling of {TERM_TEST_CEILING}; declare them explicitly or set "
                "term_tests: false"
            )
        candidate = Family(
            id=family.id,
            design=family.design,
            cells=family.cells,
            reference_levels=family.reference_levels,
            replicate_unit=family.replicate_unit,
            shrinkage=family.shrinkage,
            filter=family.filter,
            declared=family.declared,
            estimands=family.estimands,
            term_tests=tests,
        )
        errors.extend(validate_family_design(candidate, sample_rows))
        resolved.append(
            Family(
                id=candidate.id,
                design=candidate.design,
                cells=candidate.cells,
                reference_levels=candidate.reference_levels,
                replicate_unit=candidate.replicate_unit,
                shrinkage=candidate.shrinkage,
                filter=candidate.filter,
                declared=candidate.declared,
                estimands=candidate.estimands,
                term_tests=resolve_term_test_df(candidate, sample_rows),
            )
        )
    return resolved, errors


def family_rows(families: Sequence[Family]) -> list[dict[str, str]]:
    """Flat ``inputs/families.tsv`` projection."""
    return [
        {
            "family_id": family.id,
            "design": family.design,
            "cells": CELL_KEY_SEPARATOR.join(family.cells),
            "reference_levels": ";".join(
                f"{key}={value}" for key, value in sorted(family.reference_levels.items())
            ),
            "replicate_unit": family.replicate_unit or "",
            "shrinkage": family.shrinkage,
            "filter": family.filter,
            "declared": "true" if family.declared else "false",
        }
        for family in families
    ]


def estimand_rows(families: Sequence[Family]) -> list[dict[str, str]]:
    """Flat ``inputs/estimands.tsv`` projection."""
    return [
        {
            "family_id": family.id,
            "estimand_id": estimand.id,
            "label": estimand.label,
            "role": estimand.role,
            "expression": estimand.expression.source,
            "atom_kind": estimand.expression.atom_kind,
        }
        for family in families
        for estimand in family.estimands
    ]


def cell_weight_rows(families: Sequence[Family]) -> list[dict[str, str]]:
    """Long-form ``inputs/estimand_cell_weights.tsv`` projection.

    ``coef(...)`` estimands emit their coefficient name as the ``cell_key``, so
    one table carries both atom kinds and R distinguishes them via
    ``estimands.tsv``'s ``atom_kind``.
    """
    rows: list[dict[str, str]] = []
    for family in families:
        for estimand in family.estimands:
            weights = estimand.expression.cell_weights or estimand.expression.coef_weights
            for key, weight in weights.items():
                cell_key = (
                    CELL_KEY_SEPARATOR.join(key) if isinstance(key, tuple) else str(key)
                )
                rows.append(
                    {
                        "family_id": family.id,
                        "estimand_id": estimand.id,
                        "cell_key": cell_key,
                        "weight": repr(float(weight)),
                    }
                )
    return rows


def term_test_rows(families: Sequence[Family]) -> list[dict[str, str]]:
    """Flat ``inputs/term_tests.tsv`` projection."""
    return [
        {
            "family_id": family.id,
            "term_test_id": test.id,
            "reduced": test.reduced,
            "df": str(test.df),
        }
        for family in families
        for test in family.term_tests
    ]


def synthesized_contrast_rows(families: Sequence[Family], existing_ids: set[str]) -> list[dict[str, str]]:
    """One ``contrast_rows`` entry per declared estimand (spec §8).

    This is the whole compatibility seam: figures/resolve.py, report.py,
    figures/gallery.py, hypothesis validation and expected_effects all iterate
    ``project.contrast_rows``, so synthesizing here makes estimands first-class
    to every one of them without editing any of those files.
    """
    rows: list[dict[str, str]] = []
    for family in families:
        if not family.declared:
            continue
        for estimand in family.estimands:
            if estimand.id in existing_ids:
                continue
            rows.append(
                {
                    "contrast_id": estimand.id,
                    "factor": "",
                    "numerator": estimand.label,
                    "denominator": "",
                    "type": "estimand",
                    "family": family.id,
                    "role": estimand.role,
                }
            )
    return rows
