"""The report falls back to the per-analysis figures when no publication figures exist."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "10_figures" / "report.py"
PNG = b"\x89PNG\r\n\x1a\n" + b"0" * 16


def _load_report():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("report_default_figures", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _entry(figure_id: str, png: bool = True) -> dict:
    artifacts = [{"path": f"figures/{figure_id}.pdf"}]
    if png:
        artifacts.append({"path": f"figures/{figure_id}.png"})
    return {"id": figure_id, "artifacts": artifacts}


def _write_run(root: Path, entries: list[dict]) -> None:
    (root / "figures").mkdir(parents=True)
    (root / "figures" / "index.json").write_text(json.dumps({"figures": entries}), encoding="utf-8")
    for entry in entries:
        for artifact in entry["artifacts"]:
            if artifact["path"].endswith(".png"):
                (root / artifact["path"]).write_bytes(PNG)


def test_default_figures_fill_an_otherwise_empty_report(tmp_path):
    entries = [
        _entry("qc__pca_correlation"),
        _entry("treated_vs_control__de_overview"),
        _entry("treated_vs_control__gsva_heatmap"),
    ]
    _write_run(tmp_path, entries)
    figures = _load_report()._collect_figures(tmp_path, None)
    assert figures["assembled"] == []
    panels = figures["panels"]
    assert [p["label"] for p in panels] == [
        "PCA and sample correlation",
        "DE overview",
        "GSVA pathway heatmap",
    ]
    assert [p["contrast"] for p in panels] == [None, "treated_vs_control", "treated_vs_control"]
    # "heatmap" must not be mistaken for a DE "ma" plot.
    assert [p["group"] for p in panels] == ["qc", "de", "enrichment"]
    assert all(p["png"].startswith("data:image/png;base64,") for p in panels)
    assert all(p["selected"] is False for p in panels)


def test_entries_without_a_png_are_skipped(tmp_path):
    _write_run(tmp_path, [_entry("qc__pca_correlation", png=False), _entry("c__de_overview")])
    panels = _load_report()._collect_figures(tmp_path, None)["panels"]
    assert [p["constructor"] for p in panels] == ["de_overview"]


def test_assembled_figures_suppress_the_fallback(tmp_path):
    _write_run(tmp_path, [_entry("fig1__assembled"), _entry("c__de_overview")])
    figures = _load_report()._collect_figures(tmp_path, None)
    assert [f["id"] for f in figures["assembled"]] == ["fig1__assembled"]
    assert figures["panels"] == []
