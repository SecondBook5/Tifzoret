"""Every sample-labelled panel uses the same display labels ("Control 1", "Treated 2"), not a mix with raw ids.

The PCA already used `sample_display_labels`; the detected-genes axis and the DE heatmap columns drew the raw
library ids ("Ctr_1"), so one figure carried two naming schemes. Tables keep the ids; only axis text changes.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from _factorial_support import require_r

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "src" / "tifzoret" / "workflow" / "scripts"


def test_display_labels_follow_condition_and_replicate_number(tmp_path):
    require_r("jsonlite", "stringr")
    script = tmp_path / "probe.R"
    script.write_text(
        f'source("{(SCRIPTS / "utils.R").as_posix()}")\n'
        'ids <- c("Ctr_1", "Ctr_5", "Trt_2", "control3", "Ko4")\n'
        'conditions <- c("control", "control", "Treated", "control", "KO")\n'
        'jsonlite::write_json(sample_display_labels(ids, conditions), "out.json")\n',
        encoding="utf-8",
    )
    result = subprocess.run(["Rscript", "--vanilla", str(script)], capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert json.loads((tmp_path / "out.json").read_text(encoding="utf-8")) == [
        "Control 1", "Control 5", "Treated 2", "Control 3", "KO 4"]


def test_sample_axes_do_not_bypass_the_display_labels():
    qc = (SCRIPTS / "02_qc" / "qc.R").read_text(encoding="utf-8")
    de_render = (SCRIPTS / "03_differential" / "de_render.R").read_text(encoding="utf-8")
    assert "scale_x_discrete(labels = detected_genes_labels)" in qc
    assert "sample_display_labels(as.character(detected_genes_table$sample_id)" in qc
    assert "scale_x_discrete(labels = stats::setNames(" in de_render
    assert "sample_display_labels(column_order, metadata[column_order, display_group_col])" in de_render
    assert "aes(label = sample_display_labels(sample_id, .data[[display_group_col]]))" in de_render
