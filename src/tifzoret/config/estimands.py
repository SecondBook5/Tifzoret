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
from dataclasses import dataclass, field
from itertools import combinations


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


def format_formula(terms) -> str:
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


def derive_term_tests(design: str, cells, family_id: str) -> list[TermTest]:
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
