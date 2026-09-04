"""The four estimand artifacts reach inputs/ for the R stages (spec §5)."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"
SCRIPT = (
    ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "01_inputs" / "materialize_inputs.py"
)


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _run(tmp_path: Path) -> Path:
    project_dir = tmp_path / "project"
    shutil.copytree(TEMPLATE, project_dir)
    config_path = project_dir / "project.yaml"
    outputs = tmp_path / "inputs"
    outputs.mkdir()
    subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--project-config", str(config_path),
            "--counts", str(outputs / "counts.tsv"),
            "--samples", str(outputs / "samples.tsv"),
            "--annotation", str(outputs / "annotation.tsv"),
            "--contrasts", str(outputs / "contrasts.tsv"),
            "--families", str(outputs / "families.tsv"),
            "--estimands", str(outputs / "estimands.tsv"),
            "--cell-weights", str(outputs / "estimand_cell_weights.tsv"),
            "--term-tests", str(outputs / "term_tests.tsv"),
            "--manifest", str(outputs / "input_manifest.json"),
        ],
        check=True,
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )
    return outputs


def test_all_four_artifacts_are_written(tmp_path):
    outputs = _run(tmp_path)
    for name in (
        "families.tsv", "estimands.tsv", "estimand_cell_weights.tsv", "term_tests.tsv"
    ):
        assert (outputs / name).is_file(), name


def test_families_table_has_the_documented_header(tmp_path):
    outputs = _run(tmp_path)
    with (outputs / "families.tsv").open(encoding="utf-8") as handle:
        header = handle.readline().strip().split("\t")
    assert header == [
        "family_id", "design", "cells", "reference_levels",
        "replicate_unit", "shrinkage", "filter", "declared",
    ]


def test_cell_weights_table_has_the_documented_header(tmp_path):
    outputs = _run(tmp_path)
    with (outputs / "estimand_cell_weights.tsv").open(encoding="utf-8") as handle:
        header = handle.readline().strip().split("\t")
    assert header == ["family_id", "estimand_id", "cell_key", "weight"]


def test_every_estimand_has_at_least_two_cell_weights(tmp_path):
    outputs = _run(tmp_path)
    estimands = _read(outputs / "estimands.tsv")
    weights = _read(outputs / "estimand_cell_weights.tsv")
    assert estimands
    by_estimand: dict[str, int] = {}
    for row in weights:
        by_estimand[row["estimand_id"]] = by_estimand.get(row["estimand_id"], 0) + 1
    for row in estimands:
        if row["atom_kind"] == "cell":
            assert by_estimand.get(row["estimand_id"], 0) >= 2


def test_cell_weights_sum_to_zero_per_estimand(tmp_path):
    outputs = _run(tmp_path)
    estimands = {r["estimand_id"]: r for r in _read(outputs / "estimands.tsv")}
    totals: dict[str, float] = {}
    for row in _read(outputs / "estimand_cell_weights.tsv"):
        totals[row["estimand_id"]] = totals.get(row["estimand_id"], 0.0) + float(row["weight"])
    for estimand_id, total in totals.items():
        if estimands[estimand_id]["atom_kind"] == "cell":
            assert abs(total) < 1e-9


def test_every_estimand_references_a_declared_family(tmp_path):
    outputs = _run(tmp_path)
    family_ids = {row["family_id"] for row in _read(outputs / "families.tsv")}
    for row in _read(outputs / "estimands.tsv"):
        assert row["family_id"] in family_ids
