"""Family compilation, desugaring, and row generation.

Takes declared families or legacy contrasts.tsv rows and produces Family objects
with their estimands and term tests. Also generates TSV rows for staging.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from itertools import combinations

from .parser import (
    CELL_KEY_SEPARATOR,
    FILTER_CHOICES,
    SHRINKAGE_CHOICES,
    WEIGHT_SUM_TOLERANCE,
    EstimandSyntaxError,
    parse_expression,
)
from .types import Estimand, Family, TermTest
from .validator import _quote_level, family_slug


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
    from .validator import resolve_term_test_df, validate_family_design

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
