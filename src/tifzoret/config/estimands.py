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
