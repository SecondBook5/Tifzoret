"""Estimand data types: Family, Estimand, TermTest."""

from __future__ import annotations

from dataclasses import dataclass

from .parser import ParsedExpression


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
