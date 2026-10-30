"""A shrunken effect that the unshrunk data rule out is an optimiser failure, not shrinkage.

apeglm maximises a posterior that is bimodal for a gene whose counts are zero in one
group (the likelihood in the log fold change is flat, so a Cauchy prior leaves a
second mode near zero). On a real 3 v 3 study 54 of 417 strongly significant genes
with a raw |log2FC| above 4 were reported at under a quarter of their raw value --
one gene with 16/262/793 reads against 0/0/0 was displayed at -0.15 -- and which genes
collapse changed with the optimiser and with the coefficient parameterisation. Genuine
shrinkage moves an estimate toward zero but leaves it compatible with the unshrunk
Wald interval; a collapse does not. The helper marks the latter so the table and the
summary can say how many displayed effects are not to be trusted.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
UTILS_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "utils.R"


def _flags(shrunk, raw, se, eligible=None):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")

    def vec(values):
        return "c(" + ",".join("NA" if v is None else repr(float(v)) for v in values) + ")"

    code = (
        "suppressPackageStartupMessages({source(commandArgs(TRUE)[1])});"
        f"cat(jsonlite::toJSON(shrinkage_outside_ci({vec(shrunk)}, {vec(raw)}, {vec(se)}"
        + ("" if eligible is None else f", eligible = c({','.join('TRUE' if e else 'FALSE' for e in eligible)})")
        + ")))"
    )
    result = subprocess.run(["Rscript", "--vanilla", "-e", code, str(UTILS_R)], capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(result.stderr)
    return json.loads(result.stdout)


def test_a_collapsed_effect_is_flagged():
    # raw -11.2 +- 2.3 (CI -15.7 .. -6.7) displayed at -0.15
    assert _flags([-0.149], [-11.18], [2.27]) == [True]


def test_ordinary_shrinkage_stays_inside_the_interval():
    # raw 3 +- 3 pulled to 0.5; raw 2.42 +- 0.39 barely moved
    assert _flags([0.5, 2.30], [3.0, 2.42], [3.0, 0.39]) == [False, False]


def test_an_unshrunken_effect_is_never_flagged():
    assert _flags([1.7, -4.0], [1.7, -4.0], [0.3, 0.9]) == [False, False]


def test_missing_or_non_finite_inputs_are_not_flagged():
    assert _flags([None, 0.2, 0.2], [1.0, None, 1.0], [0.3, 0.3, None]) == [False, False, False]


def test_only_the_genes_under_test_are_judged():
    # a non-significant gene is meant to be shrunk hard, so the same pull is not a failure
    assert _flags([-0.149, -0.149], [-11.18, -11.18], [2.27, 2.27], eligible=[True, False]) == [True, False]
