"""The config-time design-validity gate (spec §6)."""

from __future__ import annotations

import dataclasses

from tifzoret.config.estimands import (
    build_families,
    cell_membership,
    derive_term_tests,
    resolve_term_test_df,
    validate_family_design,
)


def _family(**overrides):
    config = {
        "fam_a": {
            "design": "~ factor_a * factor_b",
            "cells": ["factor_a", "factor_b"],
            "estimands": [
                {
                    "id": "est_interaction",
                    "role": "primary",
                    "expression": "(a2,b2) - (a2,b1) - (a1,b2) + (a1,b1)",
                }
            ],
            **overrides,
        }
    }
    families, errors = build_families(config, known_ids=set())
    assert errors == [], errors
    return families[0]


def _with_lattice(family):
    """Attach the derived term tests lattice to a family."""
    return dataclasses.replace(
        family,
        term_tests=tuple(derive_term_tests(family.design, family.cells, family.id)),
    )


def _balanced(per_cell=3, **extra):
    rows = []
    index = 0
    for level_a in ("a1", "a2"):
        for level_b in ("b1", "b2"):
            for _ in range(per_cell):
                index += 1
                row = {
                    "sample_id": f"s{index}",
                    "factor_a": level_a,
                    "factor_b": level_b,
                }
                row.update(extra)
                rows.append(row)
    return rows


def test_a_balanced_two_by_two_passes():
    assert validate_family_design(_family(), _balanced()) == []


def test_cell_membership_counts_every_cell():
    membership = cell_membership(_family(), _balanced())
    assert set(membership) == {("a1", "b1"), ("a1", "b2"), ("a2", "b1"), ("a2", "b2")}
    assert all(len(ids) == 3 for ids in membership.values())


def test_an_empty_referenced_cell_is_a_hard_fail():
    rows = [row for row in _balanced() if not (row["factor_a"] == "a2" and row["factor_b"] == "b2")]
    errors = validate_family_design(_family(), rows)
    assert any("empty" in error and "a2" in error for error in errors)


def test_an_unknown_level_in_an_expression_is_a_hard_fail():
    family = _family(
        estimands=[{"id": "est_bad", "role": "primary", "expression": "(a3,b1) - (a1,b1)"}]
    )
    errors = validate_family_design(family, _balanced())
    assert any("a3" in error for error in errors)


def test_a_level_containing_separator_is_a_hard_fail():
    from tifzoret.config.estimands import CELL_KEY_SEPARATOR
    # A level containing | would mis-parse silently in R when split
    family = _family(
        estimands=[{"id": "est_bad", "role": "primary", "expression": "(a2,b|special) - (a1,b1)"}]
    )
    errors = validate_family_design(family, _balanced())
    assert any(
        "b|special" in error and "separator" in error and CELL_KEY_SEPARATOR in error
        for error in errors
    )


def test_a_cell_column_absent_from_samples_is_a_hard_fail():
    rows = [{"sample_id": row["sample_id"], "factor_a": row["factor_a"]} for row in _balanced()]
    errors = validate_family_design(_family(), rows)
    assert any("factor_b" in error and "samples" in error for error in errors)


def test_perfect_confounding_with_a_nuisance_covariate_is_a_hard_fail():
    rows = _balanced()
    for row in rows:
        row["nuisance_z"] = "z1" if row["factor_a"] == "a1" else "z2"
    family = _family(design="~ factor_a * factor_b + nuisance_z")
    errors = validate_family_design(family, rows)
    assert any("confound" in error for error in errors)


def test_a_nuisance_covariate_varying_within_cells_is_accepted():
    rows = _balanced()
    for index, row in enumerate(rows):
        row["nuisance_z"] = "z1" if index % 3 == 0 else "z2"
    family = _family(design="~ factor_a * factor_b + nuisance_z")
    assert validate_family_design(family, rows) == []


def test_too_few_samples_for_the_model_is_a_hard_fail():
    errors = validate_family_design(_family(), _balanced(per_cell=1))
    assert any("residual" in error for error in errors)


def test_a_repeated_replicate_unit_within_a_cell_is_a_hard_fail():
    rows = _balanced()
    rows[0]["unit_id"] = "u1"
    rows[1]["unit_id"] = "u1"
    for index, row in enumerate(rows[2:], start=2):
        row["unit_id"] = f"u{index}"
    family = _family(replicate_unit="unit_id")
    errors = validate_family_design(family, rows)
    assert any("u1" in error and "mixed-model" in error for error in errors)


def test_distinct_replicate_units_pass():
    rows = _balanced()
    for index, row in enumerate(rows):
        row["unit_id"] = f"u{index}"
    family = _family(replicate_unit="unit_id")
    assert validate_family_design(family, rows) == []


def test_an_absent_replicate_unit_column_is_a_hard_fail():
    family = _family(replicate_unit="unit_id")
    errors = validate_family_design(family, _balanced())
    assert any("unit_id" in error for error in errors)


def test_term_test_df_uses_real_level_counts_for_a_two_by_two():
    tests = {test.id: test for test in resolve_term_test_df(_with_lattice(_family()), _balanced())}
    assert tests["factor_a_factor_b_interaction"].df == 1
    assert tests["any_effect"].df == 3


def test_term_test_df_grows_with_a_three_level_factor():
    rows = []
    index = 0
    for level_a in ("a1", "a2"):
        for level_b in ("b1", "b2", "b3"):
            for _ in range(3):
                index += 1
                rows.append(
                    {"sample_id": f"s{index}", "factor_a": level_a, "factor_b": level_b}
                )
    family = _family(
        estimands=[{"id": "est_x", "role": "primary", "expression": "(a2,b3) - (a1,b1)"}]
    )
    tests = {test.id: test for test in resolve_term_test_df(_with_lattice(family), rows)}
    # factor_b contributes 2 df; the interaction contributes 1*2 = 2 df.
    assert tests["factor_a_factor_b_interaction"].df == 2
    assert tests["factor_b_any"].df == 4
    assert tests["any_effect"].df == 5


def test_coefficient_only_family_does_not_produce_spurious_confounding_error():
    """Coefficient-only families (cells=()) have one degenerate cell containing
    all samples, so confounding checks are vacuous and must be skipped."""
    base = _family()
    # Build a coefficient-only family by clearing cells.
    coef_family = dataclasses.replace(base, cells=())
    rows = _balanced(batch="b1")
    # Add a nuisance covariate that would trigger confounding if checked.
    for index, row in enumerate(rows):
        row["batch"] = "b1" if index < 6 else "b2"
    errors = validate_family_design(coef_family, rows)
    # Should NOT produce confounding error.
    assert not any("confounded" in error for error in errors)


def test_coefficient_only_family_still_enforces_replicate_unit_constraint():
    """Coefficient-only families must still reject repeated replicate_unit values
    because the engine has no mixed-model path."""
    base = _family()
    coef_family = dataclasses.replace(base, cells=(), replicate_unit="unit_id")
    rows = _balanced()
    # Repeat unit_id "u1" twice (in the degenerate cell).
    rows[0]["unit_id"] = "u1"
    rows[1]["unit_id"] = "u1"
    for index, row in enumerate(rows[2:], start=2):
        row["unit_id"] = f"u{index}"
    errors = validate_family_design(coef_family, rows)
    # Should produce the mixed-model hard fail.
    assert any("u1" in error and "mixed-model" in error for error in errors)
