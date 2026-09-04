"""The estimand expression grammar (spec §4)."""

from __future__ import annotations

import pytest

from tifzoret.config.estimands import (
    EstimandSyntaxError,
    ParsedExpression,
    parse_expression,
)


def test_parses_a_two_factor_interaction():
    parsed = parse_expression("(a2,b2) - (a2,b1) - (a1,b2) + (a1,b1)", arity=2)
    assert parsed.atom_kind == "cell"
    assert parsed.cell_weights == {
        ("a2", "b2"): 1.0,
        ("a2", "b1"): -1.0,
        ("a1", "b2"): -1.0,
        ("a1", "b1"): 1.0,
    }
    assert parsed.source == "(a2,b2) - (a2,b1) - (a1,b2) + (a1,b1)"


def test_parses_a_simple_difference():
    parsed = parse_expression("(a2,b1) - (a1,b1)", arity=2)
    assert parsed.cell_weights == {("a2", "b1"): 1.0, ("a1", "b1"): -1.0}


def test_parses_numeric_coefficients_for_a_marginal_effect():
    parsed = parse_expression(
        "0.5*(a2,b1) + 0.5*(a2,b2) - 0.5*(a1,b1) - 0.5*(a1,b2)", arity=2
    )
    assert parsed.cell_weights == {
        ("a2", "b1"): 0.5,
        ("a2", "b2"): 0.5,
        ("a1", "b1"): -0.5,
        ("a1", "b2"): -0.5,
    }


def test_accumulates_repeated_cells_rather_than_overwriting():
    parsed = parse_expression("(a1,b1) + (a1,b1) - (a2,b1)", arity=2)
    assert parsed.cell_weights[("a1", "b1")] == 2.0


def test_parses_a_coefficient_escape_hatch():
    parsed = parse_expression("coef(covariate_x)", arity=2)
    assert parsed.atom_kind == "coef"
    assert parsed.coef_weights == {"covariate_x": 1.0}
    assert parsed.cell_weights == {}


def test_parses_arity_one_cells():
    parsed = parse_expression("(a2) - (a1)", arity=1)
    assert parsed.cell_weights == {("a2",): 1.0, ("a1",): -1.0}


def test_parses_arity_three_cells():
    parsed = parse_expression("(a2,b2,c2) - (a1,b1,c1)", arity=3)
    assert parsed.cell_weights == {("a2", "b2", "c2"): 1.0, ("a1", "b1", "c1"): -1.0}


def test_quoted_levels_may_contain_separators():
    parsed = parse_expression("('a,1',b1) - (\"a)2\",b1)", arity=2)
    assert parsed.cell_weights == {("a,1", "b1"): 1.0, ("a)2", "b1"): -1.0}


def test_whitespace_is_irrelevant():
    a = parse_expression("(a2,b2)-(a1,b1)", arity=2)
    b = parse_expression("  ( a2 , b2 )  -  ( a1 , b1 )  ", arity=2)
    assert a.cell_weights == b.cell_weights


def test_rejects_arity_mismatch():
    with pytest.raises(EstimandSyntaxError, match="arity"):
        parse_expression("(a2,b2,c2) - (a1,b1,c1)", arity=2)


def test_rejects_mixed_cell_and_coefficient_atoms():
    with pytest.raises(EstimandSyntaxError, match="mixed"):
        parse_expression("(a2,b1) - coef(covariate_x)", arity=2)


def test_rejects_an_empty_expression():
    with pytest.raises(EstimandSyntaxError, match="empty"):
        parse_expression("   ", arity=2)


def test_rejects_unbalanced_parentheses():
    with pytest.raises(EstimandSyntaxError):
        parse_expression("(a2,b1 - (a1,b1)", arity=2)


def test_rejects_a_trailing_operator():
    with pytest.raises(EstimandSyntaxError):
        parse_expression("(a2,b1) -", arity=2)


def test_rejects_a_non_numeric_coefficient():
    with pytest.raises(EstimandSyntaxError):
        parse_expression("x*(a2,b1)", arity=2)


def test_parsed_expression_is_frozen():
    parsed = parse_expression("(a2,b1) - (a1,b1)", arity=2)
    with pytest.raises(Exception):
        parsed.source = "other"  # type: ignore[misc]
    assert isinstance(parsed, ParsedExpression)
