"""DE figure titles and axes name an estimand by its label alone (no dangling 'versus' or '−')."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from _factorial_support import require_r

ROOT = Path(__file__).resolve().parents[2]
DE_RENDER = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "de_render.R"


def _labels(tmp_path: Path) -> dict:
    script = tmp_path / "probe.R"
    script.write_text(
        "suppressPackageStartupMessages({library(jsonlite)})\n"
        f'source("{DE_RENDER.as_posix()}")\n'
        'out <- list(pair = comparison_text("A", "B"),\n'
        '            joined = comparison_text("A", "B", " versus "),\n'
        '            empty = comparison_text("Interaction", ""),\n'
        '            na = comparison_text("Interaction", NA_character_),\n'
        '            null = comparison_text("Interaction", NULL),\n'
        '            axis = comparison_text("Interaction", "", " − "))\n'
        'jsonlite::write_json(out, "out.json", auto_unbox = TRUE)\n',
        encoding="utf-8",
    )
    result = subprocess.run(
        ["Rscript", "--vanilla", str(script)], capture_output=True, text=True, cwd=tmp_path
    )
    assert result.returncode == 0, result.stderr
    return json.loads((tmp_path / "out.json").read_text(encoding="utf-8"))


def test_pairwise_contrasts_keep_their_comparison_text(tmp_path):
    require_r("jsonlite")
    out = _labels(tmp_path)
    assert out["pair"] == "A vs B"
    assert out["joined"] == "A versus B"


def test_estimands_without_a_denominator_have_no_dangling_connector(tmp_path):
    require_r("jsonlite")
    out = _labels(tmp_path)
    for key in ("empty", "na", "null", "axis"):
        assert out[key] == "Interaction", (key, out[key])
