from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from PIL import Image
from reportlab.pdfgen.canvas import Canvas


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "bulk_rna_frame" / "templates" / "minimal"
SCRIPT = ROOT / "src" / "bulk_rna_frame" / "workflow" / "scripts" / "front_door.py"


def test_front_door_promotes_figure_pair_and_metadata(tmp_path: Path):
    results = tmp_path / "results" / "synthetic_demo" / "all"
    stem = results / "qc" / "figures" / "pca_correlation"
    stem.parent.mkdir(parents=True)
    Image.new("RGB", (240, 120), "white").save(stem.with_suffix(".png"))
    canvas = Canvas(str(stem.with_suffix(".pdf")), pagesize=(240, 120))
    canvas.drawString(20, 60, "PCA and correlation")
    canvas.save()
    figures_index = results / "figures" / "index.json"
    tables_index = results / "tables" / "index.json"
    subprocess.run([
        sys.executable, str(SCRIPT),
        "--project-config", str(TEMPLATE / "project.yaml"),
        "--results", str(results),
        "--figures-index", str(figures_index),
        "--tables-index", str(tables_index),
    ], check=True)
    index = json.loads(figures_index.read_text())
    assert index["figures"][0]["id"] == "qc__pca_correlation"
    metadata = results / index["figures"][0]["metadata"]
    assert metadata.is_file()
    assert (results / "figures" / "qc__pca_correlation.pdf").is_file()
    assert (results / "figures" / "qc__pca_correlation.png").is_file()
