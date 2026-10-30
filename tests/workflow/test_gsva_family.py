"""A curated GSVA panel must be tested against the whole collection it was curated from.

``pathways.R`` used to score and test only the sets named in ``gsva_sets``. limma's
Benjamini-Hochberg adjustment then ran over that hand-picked list, so every adjusted
P was computed in a family that had already been filtered by effect size. On the
50 Hallmark sets of a real 3 v 3 study that turned 13 significant sets into 18.
``gsva_sets`` now chooses what is *drawn*; scoring and multiplicity cover every
measured set of the collection(s) the listed sets belong to.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
UTILS_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "utils.R"


def _family(gene_sets: dict[str, list[str]], requested: list[str]) -> dict:
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    sets = ",".join(f'"{name}" = c("a")' for name in gene_sets)
    ids = ",".join(f'"{name}"' for name in requested)
    code = (
        "suppressPackageStartupMessages({source(commandArgs(TRUE)[1])});"
        f"gs <- list({sets}); r <- tryCatch(gsva_family_sets(gs, c({ids})), error = function(e) conditionMessage(e));"
        "cat(jsonlite::toJSON(list(result = r), auto_unbox = FALSE))"
    )
    result = subprocess.run(["Rscript", "--vanilla", "-e", code, str(UTILS_R)], capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(result.stderr)
    return json.loads(result.stdout)["result"]


SETS = {
    "HALLMARK_A": [], "HALLMARK_B": [], "HALLMARK_C": [],
    "REACTOME_X": [], "REACTOME_Y": [], "GO_BP_Z": [],
}


def test_the_family_is_the_whole_collection_not_the_listed_sets():
    family = _family(SETS, ["HALLMARK_B"])
    assert sorted(family) == ["HALLMARK_A", "HALLMARK_B", "HALLMARK_C"]


def test_sets_from_two_collections_pull_in_both_collections():
    family = _family(SETS, ["HALLMARK_A", "REACTOME_Y"])
    assert sorted(family) == ["HALLMARK_A", "HALLMARK_B", "HALLMARK_C", "REACTOME_X", "REACTOME_Y"]


def test_a_listed_set_that_is_not_measured_is_an_error_not_a_silent_skip():
    message = _family(SETS, ["HALLMARK_A", "HALLMARK_NOT_THERE"])
    assert "HALLMARK_NOT_THERE" in " ".join(map(str, message))
    assert "not measured" in " ".join(map(str, message))
