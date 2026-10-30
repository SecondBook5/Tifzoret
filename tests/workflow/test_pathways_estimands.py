"""Pathways accept weighted estimands, not just pairwise contrasts."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from _factorial_support import require_r, stage_family

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "src" / "tifzoret" / "workflow" / "scripts"
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"


def _run_r(body: str, tmp_path: Path) -> dict:
    script = tmp_path / "probe.R"
    script.write_text(
        "suppressPackageStartupMessages({library(jsonlite); library(DESeq2); library(limma)})\n"
        f'source("{(SCRIPTS / "utils.R").as_posix()}")\n'
        f'source("{(SCRIPTS / "estimands.R").as_posix()}")\n' + body,
        encoding="utf-8",
    )
    result = subprocess.run(
        ["Rscript", "--vanilla", str(script)], capture_output=True, text=True, cwd=tmp_path
    )
    assert result.returncode == 0, result.stderr
    return json.loads((tmp_path / "out.json").read_text(encoding="utf-8"))


def test_estimand_rows_are_rejected_unless_the_caller_opts_in(tmp_path):
    require_r("jsonlite")
    out = _run_r(
        'row <- data.frame(contrast_id = "e", type = "estimand", stringsAsFactors = FALSE)\n'
        'strict <- tryCatch(resolve_contrast(row, "~ g"), error = function(e) conditionMessage(e))\n'
        'ok <- resolve_contrast(row, "~ g", allow_estimand = TRUE)\n'
        'jsonlite::write_json(list(strict = strict, type = ok$type), "out.json", auto_unbox = TRUE)\n',
        tmp_path,
    )
    allowed = out["strict"].split("got:")[0]
    assert "must be" in allowed and "estimand" not in allowed
    assert out["type"] == "estimand"


def test_score_contrast_recovers_the_weighted_sum_of_cell_means(tmp_path):
    """A saturated 2x2 design: limma's contrast estimate on per-sample scores must
    equal sum(weight * cell mean) exactly, for the interaction, a simple arm and a composite."""
    require_r("DESeq2", "apeglm", "limma", "jsonlite")
    family_dir = stage_family(tmp_path, SCRIPTS / "03_differential" / "family_fit.R", TEMPLATE)
    out = _run_r(
        f'dds <- readRDS("{(family_dir / "objects" / "deseq2.rds").as_posix()}")\n'
        'family <- read_family_spec("families.tsv", "fam_a")\n'
        'result <- list()\n'
        'for (id in c("est_interaction", "est_arm_a1", "est_composite")) {\n'
        '  spec <- read_estimand_spec("estimands.tsv", id)\n'
        '  weights <- read_cell_weights("estimand_cell_weights.tsv", id)\n'
        '  sc <- estimand_score_contrast(dds, family, spec, weights)\n'
        '  meta <- as.data.frame(SummarizedExperiment::colData(dds))\n'
        '  key <- do.call(paste, c(meta[family$cells], sep = CELL_KEY_SEPARATOR))\n'
        '  set.seed(3); score <- matrix(rnorm(2 * nrow(meta)) + as.integer(factor(key)), nrow = 2,\n'
        '                              dimnames = list(c("s1", "s2"), rownames(meta)))\n'
        '  fit <- limma::contrasts.fit(limma::lmFit(score, sc$design), matrix(sc$contrast, ncol = 1))\n'
        '  expected <- vapply(1:2, function(i) sum(weights * tapply(score[i, ], key, mean)[names(weights)]), 0)\n'
        '  result[[id]] <- list(estimate = as.numeric(fit$coefficients[, 1]), expected = expected)\n'
        '}\n'
        'jsonlite::write_json(result, "out.json", auto_unbox = FALSE, digits = NA)\n',
        tmp_path,
    )
    for id_ in ("est_interaction", "est_arm_a1", "est_composite"):
        for got, want in zip(out[id_]["estimate"], out[id_]["expected"]):
            assert abs(got - want) < 1e-9, (id_, got, want)


def test_condition_palette_falls_back_for_levels_the_palette_does_not_name(tmp_path):
    """A pairwise contrast on `treatment` with a palette keyed on another grouping used
    to return NA colours (ComplexHeatmap then refused them); now it draws."""
    require_r("jsonlite")
    out = _run_r(
        'cfg <- list(figures = list(palette = list(a_x = "#112233", b_x = "#445566")))\n'
        'known <- condition_palette(cfg, c("b_x", "a_x"))\n'
        'unknown <- suppressWarnings(condition_palette(cfg, c("control", "treated")))\n'
        'warned <- tryCatch({condition_palette(cfg, "control"); FALSE}, warning = function(w) TRUE)\n'
        'jsonlite::write_json(list(known = known, known_names = names(known), unknown = unknown,\n'
        '                          unknown_names = names(unknown), warned = warned),\n'
        '                     "out.json", auto_unbox = FALSE)\n',
        tmp_path,
    )
    assert out["known"] == ["#445566", "#112233"] and out["known_names"] == ["b_x", "a_x"]
    assert out["unknown_names"] == ["control", "treated"]
    assert all(c.startswith("#") and len(c) == 7 for c in out["unknown"])
    assert len(set(out["unknown"])) == 2
    assert out["warned"] is True or out["warned"] == [True]
