"""Estimand subsystem: expression parsing, family compilation, design validation.

This subpackage provides the core logic for turning declared estimands (cell-means
arithmetic) into validated Family objects ready for DESeq2 fitting.
"""

from .compiler import (
    TERM_TEST_CEILING,
    build_families,
    cell_weight_rows,
    compile_project_families,
    derive_term_tests,
    desugar_contrast_rows,
    estimand_rows,
    family_rows,
    format_formula,
    parse_expression,
    parse_formula_terms,
    synthesized_contrast_rows,
    term_test_rows,
)
from .parser import (
    CELL_KEY_SEPARATOR,
    FILTER_CHOICES,
    SHRINKAGE_CHOICES,
    EstimandSyntaxError,
    ParsedExpression,
)
from .types import Estimand, Family, TermTest
from .validator import cell_membership, family_slug, resolve_term_test_df, validate_family_design

__all__ = [
    "CELL_KEY_SEPARATOR",
    "Estimand",
    "EstimandSyntaxError",
    "FILTER_CHOICES",
    "Family",
    "ParsedExpression",
    "SHRINKAGE_CHOICES",
    "TERM_TEST_CEILING",
    "TermTest",
    "build_families",
    "cell_membership",
    "cell_weight_rows",
    "compile_project_families",
    "derive_term_tests",
    "desugar_contrast_rows",
    "estimand_rows",
    "family_rows",
    "family_slug",
    "format_formula",
    "parse_expression",
    "parse_formula_terms",
    "resolve_term_test_df",
    "synthesized_contrast_rows",
    "term_test_rows",
    "validate_family_design",
]
