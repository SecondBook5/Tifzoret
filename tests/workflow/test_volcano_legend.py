"""The clipped volcano lists all four gene classes in its legend with a key each.

A class whose only genes sit off-scale (drawn as legend-less triangles), or that has no genes at all, used
to get a blank legend key or none. An empty layer now carries the key for such classes; when every class
has a visible point nothing is added, so a volcano that was already complete is unchanged.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from _factorial_support import require_r

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "src" / "tifzoret" / "workflow" / "scripts"
UTILS = SCRIPTS / "utils.R"
DE_RENDER = SCRIPTS / "03_differential" / "de_render.R"

PROBE = """
suppressPackageStartupMessages({{library(dplyr); library(ggplot2); library(ggrepel); library(jsonlite)}})
source("{utils}")
source("{de_render}")
make_rows <- function(n, id, lfc, nlp, padj, class) tibble::tibble(
  gene_id = paste0(id, seq_len(n)), gene_symbol = paste0(id, seq_len(n)),
  log2_fold_change = lfc + seq_len(n) * 0.01, negative_log10_p = nlp + seq_len(n) * 0.05,
  adjusted_p_value = padj, significance_class = class, base_mean = 100)
rows <- bind_rows(
  make_rows(40, "ns", 0.1, 0.5, 0.9, "ns"),
  make_rows(4, "padj", 0.5, 4, 0.01, "padj_only"),
  make_rows(4, "lfc", 1.5, 1.5, 0.3, "lfc_only"),
  {up_rows}
)
cfg <- list(figures = list(de = list(fdr = 0.05, abs_log2fc = 1, top_labels = 6)))
built <- build_clipped_volcano(rows, cfg, "Treated", "Control")
legend_layers <- Filter(function(l) is.data.frame(l$data) && nrow(l$data) > 0L &&
                          all(is.na(l$data$log2_fold_change)), built$plot$layers)
keyed <- as.character(unlist(lapply(legend_layers, function(l) as.character(l$data$significance_class))))
scale <- built$plot$scales$get_scales("colour")
out <- list(legend_only_classes = sort(keyed), breaks = scale$breaks, labels = unname(scale$labels))
jsonlite::write_json(out, "out.json", auto_unbox = FALSE)
"""

VISIBLE_UP = 'make_rows(3, "up", 2, 4.5, 0.001, "significant_up")'
OFFSCALE_UP = 'make_rows(1, "up", 2, 1e6, 0.001, "significant_up")'


def _probe(tmp_path: Path, up_rows: str) -> dict:
    script = tmp_path / "probe.R"
    script.write_text(
        PROBE.format(utils=UTILS.as_posix(), de_render=DE_RENDER.as_posix(), up_rows=up_rows), encoding="utf-8"
    )
    result = subprocess.run(["Rscript", "--vanilla", str(script)], capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    return json.loads((tmp_path / "out.json").read_text(encoding="utf-8"))


def test_classes_without_a_visible_point_get_a_legend_key(tmp_path):
    require_r("jsonlite", "ggrepel")
    out = _probe(tmp_path, OFFSCALE_UP)
    # the only significant-up gene is off-scale and there is no significant-down gene
    assert out["legend_only_classes"] == ["significant_down", "significant_up"], out
    assert len(out["labels"]) == 4, out


def test_a_volcano_with_every_class_visible_is_unchanged(tmp_path):
    require_r("jsonlite", "ggrepel")
    out = _probe(tmp_path, VISIBLE_UP + ",\n  " + 'make_rows(3, "down", -2, 4.5, 0.001, "significant_down")')
    assert out["legend_only_classes"] == [], out
    assert len(out["labels"]) == 4, out
