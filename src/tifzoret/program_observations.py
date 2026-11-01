"""Observation-only program scores from explicitly normalized expression.

This boundary does not materialize integer counts, normalize libraries, perform DE,
or estimate cell fractions. It permits reuse of published normalized abundance
with exact gene coverage and explicit observation-scale provenance.
"""

from __future__ import annotations
from collections.abc import Mapping, Sequence
import numpy as np
import pandas as pd


def normalized_program_observations(
    expression: pd.DataFrame, annotation: pd.DataFrame, programs: Mapping[str, Sequence[str]], *, scale: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return equal-gene mean log1p-abundance observations and strict coverage.

    Args:
        expression: Nonnegative genes-by-library normalized abundances, unique IDs.
        annotation: Unique gene_id rows and gene_symbol column; ambiguous mappings fail.
        programs: Externally specified nonempty, nonoverlapping symbol lists.
        scale: Explicit supported abundance declaration; raw counts are not accepted.
    Returns:
        Long program/library observations and gene mapping/coverage rows.
    Raises:
        ValueError: Invalid, ambiguous or incomplete input. No gene is silently dropped.
    """
    if scale not in {"CPM", "TPM"}:
        raise ValueError("This observation-only boundary requires declared CPM or TPM")
    if not expression.index.is_unique or not expression.columns.is_unique:
        raise ValueError("Expression gene and library identifiers must be unique")
    values = expression.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Normalized abundances must be finite and nonnegative")
    ann = annotation[["gene_id", "gene_symbol"]].copy()
    if not ann.gene_id.is_unique:
        raise ValueError("Annotation gene IDs must be unique")
    coverage = []
    scores = []
    seen = set()
    for name, genes in programs.items():
        if not genes or len(set(genes)) != len(genes) or seen.intersection(genes):
            raise ValueError("Programs must be nonempty with distinct nonoverlapping genes")
        seen.update(genes)
        ids = []
        for gene in genes:
            selected = ann[ann.gene_symbol.eq(gene) & ann.gene_id.isin(expression.index)]
            if len(selected) != 1:
                raise ValueError(f"{name}/{gene}: expected exactly one measured mapping, found {len(selected)}")
            gid = selected.gene_id.iloc[0]
            ids.append(gid)
            coverage.append(
                dict(
                    program=name,
                    gene_symbol=gene,
                    gene_id=gid,
                    nonzero_libraries=int(expression.loc[gid].gt(0).sum()),
                    libraries=len(expression.columns),
                )
            )
        coordinate = np.log1p(expression.loc[ids].to_numpy(float)).mean(axis=0)
        for library, value in zip(expression.columns, coordinate):
            scores.append(
                dict(
                    program=name,
                    library=library,
                    score=float(value),
                    input_scale=scale,
                    definition="equal-gene mean log1p abundance; no physical calibration",
                )
            )
    return pd.DataFrame(scores), pd.DataFrame(coverage)
