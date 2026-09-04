"""contrasts.tsv desugars into families with one shared fit (spec §9)."""

from __future__ import annotations

from tifzoret.config.estimands import desugar_contrast_rows, family_slug


def _row(**overrides):
    row = {
        "contrast_id": "treated_vs_ctrl",
        "factor": "condition",
        "numerator": "treated",
        "denominator": "ctrl",
        "type": "",
        "design": "",
        "coefficient": "",
        "reduced": "",
        "reference_levels": "",
    }
    row.update(overrides)
    return row


def test_a_pairwise_row_becomes_a_single_factor_family():
    families, errors = desugar_contrast_rows("~ condition", [_row()])
    assert errors == []
    assert len(families) == 1
    family = families[0]
    assert family.id == "main"
    assert family.cells == ("condition",)
    assert family.declared is False
    estimand = family.estimands[0]
    assert estimand.id == "treated_vs_ctrl"
    assert estimand.expression.cell_weights == {("treated",): 1.0, ("ctrl",): -1.0}


def test_a_desugared_family_carries_no_primary_role():
    families, _ = desugar_contrast_rows("~ condition", [_row()])
    assert families[0].estimands[0].role == ""


def test_rows_sharing_a_design_share_one_family():
    rows = [
        _row(contrast_id="a_vs_ctrl", numerator="a"),
        _row(contrast_id="b_vs_ctrl", numerator="b"),
        _row(contrast_id="c_vs_ctrl", numerator="c"),
    ]
    families, errors = desugar_contrast_rows("~ condition", rows)
    assert errors == []
    assert len(families) == 1
    assert len(families[0].estimands) == 3


def test_a_per_row_design_forms_its_own_family():
    rows = [
        _row(contrast_id="a_vs_ctrl", numerator="a"),
        _row(contrast_id="b_vs_ctrl", numerator="b", design="~ condition + nuisance_z"),
    ]
    families, errors = desugar_contrast_rows("~ condition", rows)
    assert errors == []
    assert len(families) == 2
    assert {family.id for family in families} == {"main", "design_condition_nuisance_z"}


def test_a_coefficient_row_becomes_a_coef_estimand():
    rows = [
        _row(
            contrast_id="interaction_row",
            type="coefficient",
            coefficient="factor_aa2.factor_bb2",
            design="~ factor_a * factor_b",
        )
    ]
    families, errors = desugar_contrast_rows("~ factor_a", rows)
    assert errors == []
    estimand = families[0].estimands[0]
    assert estimand.expression.atom_kind == "coef"
    assert estimand.expression.coef_weights == {"factor_aa2.factor_bb2": 1.0}


def test_an_omnibus_row_becomes_a_term_test_not_an_estimand():
    rows = [
        _row(
            contrast_id="condition_omnibus",
            type="omnibus",
            reduced="~ 1",
            numerator="",
            denominator="",
        )
    ]
    families, errors = desugar_contrast_rows("~ condition", rows)
    assert errors == []
    family = families[0]
    assert family.estimands == ()
    assert len(family.term_tests) == 1
    assert family.term_tests[0].id == "condition_omnibus"
    assert family.term_tests[0].reduced == "~ 1"


def test_reference_levels_carry_onto_the_family():
    rows = [_row(reference_levels="condition=ctrl")]
    families, _ = desugar_contrast_rows("~ condition", rows)
    assert families[0].reference_levels == {"condition": "ctrl"}


def test_conflicting_reference_levels_in_one_family_is_an_error():
    rows = [
        _row(contrast_id="a_vs_ctrl", numerator="a", reference_levels="condition=ctrl"),
        _row(contrast_id="b_vs_a", numerator="b", denominator="a", reference_levels="condition=a"),
    ]
    _, errors = desugar_contrast_rows("~ condition", rows)
    assert any("reference_levels" in error for error in errors)


def test_desugared_families_default_to_best_practice_knobs():
    families, _ = desugar_contrast_rows("~ condition", [_row()])
    assert families[0].shrinkage == "apeglm"
    assert families[0].filter == "design_aware"


def test_family_slug_is_deterministic_and_filesystem_safe():
    assert family_slug("~ condition + nuisance_z") == "design_condition_nuisance_z"
    assert family_slug("~ factor_a * factor_b") == "design_factor_a_factor_b"
    assert family_slug("~ condition") == family_slug("~  condition  ")


def test_two_designs_that_slug_identically_produce_separate_families():
    """Two designs that slug to the same ID produce separate families with numeric suffixes."""
    rows = [
        _row(contrast_id="one", numerator="a", design="~ condition + batch"),
        _row(contrast_id="two", numerator="b", design="~ condition * batch"),
    ]
    families, errors = desugar_contrast_rows("~ condition", rows)
    assert errors == []
    assert len(families) == 2

    # Find families by their designs
    family_by_design = {f.design: f for f in families}

    assert "~ condition + batch" in family_by_design
    assert "~ condition * batch" in family_by_design

    # Check that the estimands are assigned to the right families
    additive_family = family_by_design["~ condition + batch"]
    multiplicative_family = family_by_design["~ condition * batch"]

    assert len(additive_family.estimands) == 1
    assert additive_family.estimands[0].id == "one"

    assert len(multiplicative_family.estimands) == 1
    assert multiplicative_family.estimands[0].id == "two"

    # Check the IDs have numeric suffixes for collision
    assert additive_family.id == "design_condition_batch"
    assert multiplicative_family.id == "design_condition_batch_2"


def test_a_level_containing_a_single_quote_roundtrips():
    """A level containing a single quote is correctly quoted and parsed back."""
    rows = [
        _row(
            contrast_id="quote_level",
            numerator="it's",
            denominator="plain",
        )
    ]
    families, errors = desugar_contrast_rows("~ condition", rows)
    assert errors == []
    estimand = families[0].estimands[0]

    # The expression should have the correct cell weights with the original level string
    assert estimand.expression.cell_weights == {("it's",): 1.0, ("plain",): -1.0}
