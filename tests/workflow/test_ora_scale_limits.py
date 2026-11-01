"""The GO over-representation bubble plot must draw every displayed term.

``ontology.R`` fixed the colour scale at ``limits = c(4, max)`` and the size scale at
``limits = c(0, 80)``. A ggplot continuous scale turns every value outside its
limits into ``NA`` and removes the point, so a displayed term with
-log10(FDR) below 4 (an adjusted P above 1e-4) or with more than 80 genes was
listed in ``ontology_displayed.tsv`` and absent from the figure. The limits are
therefore derived from the displayed terms.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
UTILS_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "utils.R"


def _limits(scores: list[float], counts: list[int]) -> dict:
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    code = (
        "suppressPackageStartupMessages({source(commandArgs(TRUE)[1])});"
        f"s <- c({','.join(map(str, scores))}); n <- c({','.join(map(str, counts))});"
        "cat(jsonlite::toJSON(ora_scale_limits(s, n), auto_unbox = TRUE))"
    )
    result = subprocess.run(["Rscript", "--vanilla", "-e", code, str(UTILS_R)], capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(result.stderr)
    return json.loads(result.stdout)


def test_every_displayed_term_falls_inside_both_scales():
    scores = [1.4, 2.5, 3.1, 9.7, 16.2]
    counts = [5, 12, 40, 95, 130]
    limits = _limits(scores, counts)
    assert limits["score"][0] <= min(scores) and limits["score"][1] >= max(scores)
    assert limits["count"][0] <= min(counts) and limits["count"][1] >= max(counts)
    assert all(limits["score"][0] <= b <= limits["score"][1] for b in limits["score_breaks"])
    assert all(0 < b <= limits["count"][1] for b in limits["count_breaks"])


def test_the_reference_range_keeps_its_familiar_breaks():
    # scores 4.3-16.1 and counts up to 72 are the range the manuscript panel spans
    limits = _limits([4.3, 8.0, 12.5, 16.1], [10, 30, 55, 72])
    assert limits["score"] == [4, 17]
    assert limits["score_breaks"] == [4, 8, 12, 16]
    assert limits["count"] == [0, 80]
    assert limits["count_breaks"] == [10, 30, 50, 70]


def test_a_single_term_does_not_collapse_the_scale():
    limits = _limits([6.2], [18])
    assert limits["score"][0] < limits["score"][1]
