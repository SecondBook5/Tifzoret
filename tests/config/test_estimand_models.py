"""Family/estimand model construction and structural validation (spec §3, §4, §6)."""

from __future__ import annotations

from tifzoret.config.estimands import (
    FILTER_CHOICES,
    SHRINKAGE_CHOICES,
    build_families,
)


def _config(**overrides):
    family = {
        "design": "~ factor_a * factor_b",
        "cells": ["factor_a", "factor_b"],
        "estimands": [
            {
                "id": "est_interaction",
                "label": "Interaction",
                "role": "primary",
                "expression": "(a2,b2) - (a2,b1) - (a1,b2) + (a1,b1)",
            }
        ],
    }
    family.update(overrides)
    return {"fam_a": family}


def test_builds_a_family_with_defaults_applied():
    families, errors = build_families(_config(), known_ids=set())
    assert errors == []
    assert len(families) == 1
    family = families[0]
    assert family.id == "fam_a"
    assert family.cells == ("factor_a", "factor_b")
    assert family.shrinkage == "apeglm"
    assert family.filter == "design_aware"
    assert family.declared is True
    assert family.replicate_unit is None
    assert family.estimands[0].id == "est_interaction"
    assert family.estimands[0].family_id == "fam_a"
    assert family.estimands[0].role == "primary"


def test_default_choices_are_the_documented_best_practice():
    assert SHRINKAGE_CHOICES[0] == "apeglm"
    assert FILTER_CHOICES[0] == "design_aware"


def test_explicit_knobs_override_defaults():
    families, errors = build_families(
        _config(shrinkage="ashr", filter="total_count", replicate_unit="unit_id"),
        known_ids=set(),
    )
    assert errors == []
    assert families[0].shrinkage == "ashr"
    assert families[0].filter == "total_count"
    assert families[0].replicate_unit == "unit_id"


def test_rejects_an_unknown_shrinkage_choice():
    _, errors = build_families(_config(shrinkage="magic"), known_ids=set())
    assert any("shrinkage" in error and "magic" in error for error in errors)


def test_rejects_an_unknown_filter_choice():
    _, errors = build_families(_config(filter="magic"), known_ids=set())
    assert any("filter" in error and "magic" in error for error in errors)


def test_rejects_cell_arity_mismatch_with_a_useful_message():
    _, errors = build_families(
        _config(
            estimands=[
                {"id": "est_bad", "role": "primary", "expression": "(a2,b2,c2) - (a1,b1,c1)"}
            ]
        ),
        known_ids=set(),
    )
    assert any("arity" in error for error in errors)


def test_rejects_cell_weights_that_do_not_sum_to_zero():
    _, errors = build_families(
        _config(
            estimands=[{"id": "est_bad", "role": "primary", "expression": "(a2,b2) + (a1,b1)"}]
        ),
        known_ids=set(),
    )
    assert any("sum to zero" in error for error in errors)


def test_coef_expressions_are_exempt_from_sum_to_zero():
    families, errors = build_families(
        _config(
            estimands=[
                {"id": "est_slope", "role": "primary", "expression": "coef(covariate_x)"}
            ]
        ),
        known_ids=set(),
    )
    assert errors == []
    assert families[0].estimands[0].expression.atom_kind == "coef"


def test_requires_exactly_one_primary_in_a_declared_family():
    _, errors = build_families(
        _config(
            estimands=[
                {"id": "est_one", "expression": "(a2,b1) - (a1,b1)"},
                {"id": "est_two", "expression": "(a1,b2) - (a1,b1)"},
            ]
        ),
        known_ids=set(),
    )
    assert any("role: primary" in error for error in errors)


def test_rejects_two_primaries():
    _, errors = build_families(
        _config(
            estimands=[
                {"id": "est_one", "role": "primary", "expression": "(a2,b1) - (a1,b1)"},
                {"id": "est_two", "role": "primary", "expression": "(a1,b2) - (a1,b1)"},
            ]
        ),
        known_ids=set(),
    )
    assert any("role: primary" in error for error in errors)


def test_rejects_an_estimand_id_colliding_with_a_known_id():
    _, errors = build_families(_config(), known_ids={"est_interaction"})
    assert any("collides" in error for error in errors)


def test_rejects_duplicate_estimand_ids_across_families():
    config = _config()
    config["fam_b"] = dict(config["fam_a"])
    _, errors = build_families(config, known_ids=set())
    assert any("collides" in error for error in errors)


def test_rejects_a_cell_column_absent_from_the_design():
    _, errors = build_families(
        _config(design="~ factor_a", cells=["factor_a", "factor_b"]), known_ids=set()
    )
    assert any("factor_b" in error and "design" in error for error in errors)


def test_rejects_an_empty_estimand_list():
    _, errors = build_families(_config(estimands=[]), known_ids=set())
    assert any("at least one estimand" in error for error in errors)


def test_accumulates_every_error_rather_than_raising_on_the_first():
    _, errors = build_families(
        _config(shrinkage="magic", filter="magic", estimands=[]), known_ids=set()
    )
    assert len(errors) >= 3


def test_supports_a_single_factor_family():
    families, errors = build_families(
        {
            "fam_k1": {
                "design": "~ factor_a",
                "cells": ["factor_a"],
                "estimands": [
                    {"id": "est_pair", "role": "primary", "expression": "(a2) - (a1)"}
                ],
            }
        },
        known_ids=set(),
    )
    assert errors == []
    assert families[0].cells == ("factor_a",)
