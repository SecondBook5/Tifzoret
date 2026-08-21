from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

# curvature.py solves the 1-Wasserstein transport with scipy.optimize.linprog.
# scipy is a conda-only runtime dependency (workflow/envs/network.yaml and the
# --no-conda environment.yaml), not a pip dependency of the package, so skip this
# execution test cleanly where the minimal install lacks it.
pytest.importorskip("scipy")


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "curvature.py"
PROJECT = ROOT / "src" / "tifzoret" / "templates" / "minimal" / "project.yaml"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_curvature_flags_the_bridge_edge_as_negative(tmp_path: Path) -> None:
    """Ollivier-Ricci curvature is the whole point of the module: an edge inside a
    dense (redundant) neighbourhood is positively curved, while a bridge joining
    two otherwise-separate communities is negatively curved. Two triangles joined
    by a single C--D edge is the canonical case — C--D must come out negative and
    be reported as an inter-module bridge; a triangle edge must come out positive."""
    nodes = tmp_path / "nodes.tsv"
    edges = tmp_path / "edges.tsv"
    outdir = tmp_path / "out"

    write_tsv(
        nodes,
        ["gene_id", "gene_symbol", "module"],
        [
            {"gene_id": "A", "gene_symbol": "Aaa", "module": "M1"},
            {"gene_id": "B", "gene_symbol": "Bbb", "module": "M1"},
            {"gene_id": "C", "gene_symbol": "Ccc", "module": "M1"},
            {"gene_id": "D", "gene_symbol": "Ddd", "module": "M2"},
            {"gene_id": "E", "gene_symbol": "Eee", "module": "M2"},
            {"gene_id": "F", "gene_symbol": "Fff", "module": "M2"},
        ],
    )
    write_tsv(
        edges,
        ["source", "target", "weight"],
        [
            {"source": "A", "target": "B", "weight": 1.0},
            {"source": "A", "target": "C", "weight": 1.0},
            {"source": "B", "target": "C", "weight": 1.0},
            {"source": "D", "target": "E", "weight": 1.0},
            {"source": "D", "target": "F", "weight": 1.0},
            {"source": "E", "target": "F", "weight": 1.0},
            {"source": "C", "target": "D", "weight": 1.0},  # the bridge
        ],
    )

    subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--project-config", str(PROJECT),
            "--edges", str(edges),
            "--nodes", str(nodes),
            "--contrast-id", "treated_vs_control",
            "--outdir", str(outdir),
        ],
        check=True,
    )

    summary = json.loads((outdir / "curvature_summary.json").read_text(encoding="utf-8"))
    assert summary["method"] == "ollivier_ricci"
    assert summary["alpha"] == 0.5
    assert (summary["nodes"], summary["edges"]) == (6, 7)
    assert summary["modules"] == 2
    assert summary["curvature_failures"] == 0
    assert summary["negative_edges"] >= 1
    assert summary["bridges_reported"] >= 1

    curvature = {
        frozenset({row["source"], row["target"]}): float(row["curvature"])
        for row in _read_tsv(outdir / "tables" / "edge_curvature.tsv")
    }
    assert curvature[frozenset({"C", "D"})] < 0  # bridge: bottleneck
    assert curvature[frozenset({"A", "B"})] > 0  # triangle: locally redundant

    bridges = _read_tsv(outdir / "tables" / "curvature_bridges.tsv")
    assert any(frozenset({row["source"], row["target"]}) == frozenset({"C", "D"}) for row in bridges)
    assert all(row["inter_module"] == "True" and float(row["curvature"]) < 0 for row in bridges)

    for stem in ("curvature_histogram", "curvature_network"):
        for suffix in ("pdf", "png"):
            assert (outdir / "figures" / f"{stem}.{suffix}").stat().st_size > 0
