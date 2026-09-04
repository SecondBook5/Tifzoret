"""Auto-derived full-vs-reduced term tests (spec §7.4)."""

from __future__ import annotations

from tifzoret.config.estimands import (
    TERM_TEST_CEILING,
    derive_term_tests,
    format_formula,
    parse_formula_terms,
)


def test_expands_a_crossed_two_factor_formula():
    terms = parse_formula_terms("~ factor_a * factor_b")
    assert terms == [
        frozenset({"factor_a"}),
        frozenset({"factor_b"}),
        frozenset({"factor_a", "factor_b"}),
    ]


def test_keeps_additive_terms_separate():
    terms = parse_formula_terms("~ factor_a + factor_b")
    assert frozenset({"factor_a", "factor_b"}) not in terms


def test_parses_an_explicit_colon_interaction():
    terms = parse_formula_terms("~ factor_a + factor_b + factor_a:factor_b")
    assert frozenset({"factor_a", "factor_b"}) in terms


def test_expands_a_three_way_crossing():
    terms = parse_formula_terms("~ factor_a * factor_b * factor_c")
    assert len(terms) == 7
    assert frozenset({"factor_a", "factor_b", "factor_c"}) in terms


def test_formats_terms_back_to_a_formula():
    assert format_formula(parse_formula_terms("~ factor_a * factor_b")) == (
        "~ factor_a + factor_b + factor_a:factor_b"
    )


def test_formats_an_empty_term_list_as_intercept_only():
    assert format_formula([]) == "~ 1"


def test_derives_exactly_four_tests_for_a_two_by_two():
    tests = derive_term_tests("~ factor_a * factor_b", ["factor_a", "factor_b"], "fam_a")
    assert len(tests) == 4
    by_id = {test.id: test for test in tests}
    assert by_id["factor_a_factor_b_interaction"].df == 1
    assert by_id["factor_a_any"].df == 2
    assert by_id["factor_b_any"].df == 2
    assert by_id["any_effect"].df == 3


def test_the_interaction_test_keeps_both_main_effects():
    tests = derive_term_tests("~ factor_a * factor_b", ["factor_a", "factor_b"], "fam_a")
    interaction = next(t for t in tests if t.id == "factor_a_factor_b_interaction")
    assert interaction.reduced == "~ factor_a + factor_b"


def test_the_any_effect_test_drops_every_cell_term():
    tests = derive_term_tests("~ factor_a * factor_b", ["factor_a", "factor_b"], "fam_a")
    any_effect = next(t for t in tests if t.id == "any_effect")
    assert any_effect.reduced == "~ 1"


def test_nuisance_terms_appear_in_every_reduced_model():
    tests = derive_term_tests(
        "~ factor_a * factor_b + nuisance_z", ["factor_a", "factor_b"], "fam_a"
    )
    assert len(tests) == 4
    for test in tests:
        assert "nuisance_z" in test.reduced


def test_nuisance_only_reduced_model_is_not_intercept_only():
    tests = derive_term_tests(
        "~ factor_a * factor_b + nuisance_z", ["factor_a", "factor_b"], "fam_a"
    )
    any_effect = next(t for t in tests if t.id == "any_effect")
    assert any_effect.reduced == "~ nuisance_z"


def test_derives_eleven_tests_for_a_fully_crossed_three_factor_design():
    tests = derive_term_tests(
        "~ factor_a * factor_b * factor_c",
        ["factor_a", "factor_b", "factor_c"],
        "fam_a",
    )
    assert len(tests) == 11


def test_a_single_factor_family_yields_one_test():
    tests = derive_term_tests("~ factor_a", ["factor_a"], "fam_k1")
    assert len(tests) == 1
    assert tests[0].id == "any_effect"
    assert tests[0].df == 1


def test_every_derived_test_carries_its_family_id():
    tests = derive_term_tests("~ factor_a * factor_b", ["factor_a", "factor_b"], "fam_a")
    assert {test.family_id for test in tests} == {"fam_a"}


def test_derived_ids_are_unique():
    tests = derive_term_tests(
        "~ factor_a * factor_b * factor_c",
        ["factor_a", "factor_b", "factor_c"],
        "fam_a",
    )
    assert len({test.id for test in tests}) == len(tests)


def test_ceiling_is_the_documented_value():
    assert TERM_TEST_CEILING == 12
