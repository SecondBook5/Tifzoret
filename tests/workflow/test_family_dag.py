"""The family rules replace contrast_de/contrast_omnibus (spec §7.1, §7.2)."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / "src" / "tifzoret" / "workflow"
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"


def test_de_r_and_omnibus_r_are_gone():
    assert not (WORKFLOW / "scripts" / "03_differential" / "de.R").exists()
    assert not (WORKFLOW / "scripts" / "03_differential" / "omnibus.R").exists()


def test_core_smk_declares_the_three_family_rules():
    text = (WORKFLOW / "rules" / "core.smk").read_text(encoding="utf-8")
    for rule in ("rule family_fit:", "rule family_estimand:", "rule family_term_test:"):
        assert rule in text, rule
    assert "rule contrast_de:" not in text
    assert "rule contrast_omnibus:" not in text


def test_estimand_outputs_keep_the_existing_contrast_paths():
    text = (WORKFLOW / "rules" / "core.smk").read_text(encoding="utf-8")
    assert 'analysis("de", "tables/de_results.tsv")' in text


def test_snakefile_exposes_the_family_id_lists():
    text = (WORKFLOW / "Snakefile").read_text(encoding="utf-8")
    for name in ("FAMILY_IDS", "ESTIMAND_IDS", "FAMILY_ESTIMANDS", "FAMILY_TERM_TESTS"):
        assert name in text, name


@pytest.mark.parametrize("target", ["-n"])
def test_the_dag_builds_for_the_minimal_template(tmp_path, target):
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake not available")
    project = tmp_path / "project"
    shutil.copytree(TEMPLATE, project)
    result = subprocess.run(
        [sys.executable, "-m", "snakemake", "--snakefile", str(WORKFLOW / "Snakefile"),
         "--configfile", str(project / "project.yaml"), target, "--quiet"],
        capture_output=True, text=True, cwd=tmp_path,
        env={"PYTHONPATH": str(ROOT / "src"), "PATH": "/usr/bin:/bin"})
    assert result.returncode == 0, result.stdout + result.stderr


def test_the_dag_builds_with_a_declared_family(tmp_path):
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake not available")
    project = tmp_path / "project"
    shutil.copytree(TEMPLATE, project)
    config_path = project / "project.yaml"
    data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    samples = (project / "samples.tsv").read_text(encoding="utf-8").splitlines()
    header = samples[0].split("\t")
    factor = "condition" if "condition" in header else header[1]
    levels = sorted({line.split("\t")[header.index(factor)] for line in samples[1:] if line.strip()})
    data["analysis"]["families"] = {
        "fam_a": {
            "design": f"~ {factor}",
            "cells": [factor],
            "estimands": [
                {"id": "est_declared", "role": "primary",
                 "expression": f"({levels[1]}) - ({levels[0]})"}
            ],
        }
    }
    data["analysis"].setdefault("modules", {})["factorial"] = False
    config_path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, "-m", "tifzoret", "dry-run", str(config_path),
         "--no-conda", "--cores", "1"],
        capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "family_fit" in result.stdout or "family_fit" in result.stderr


def test_readme_stage_count_was_updated():
    text = (ROOT / "docs" / "README.md").read_text(encoding="utf-8")
    assert "35 stages" not in text


def _de_patterns_block() -> str:
    """The literal ``DE_PATTERNS = [...]`` list body from the Snakefile."""
    text = (WORKFLOW / "Snakefile").read_text(encoding="utf-8")
    return text.split("DE_PATTERNS = [", 1)[1].split("]", 1)[0]


def test_de_patterns_declare_exactly_what_estimand_r_writes():
    """``family_estimand``'s declared outputs must match what ``estimand.R`` writes.

    ``de.R`` used to ``saveRDS`` a per-contrast DESeqDataSet, so ``DE_PATTERNS``
    declared ``objects/deseq2.rds``. When ``de.R`` was retired the declaration
    outlived its writer: ``estimand.R`` reads the *family* fit and persists no
    object of its own, so every ``family_estimand`` job exited 0 and then died
    with ``MissingOutputException`` on an empty ``objects/`` directory.

    Asserting the two sides move together catches the omission in either
    direction -- declaring an object nothing writes, or writing one nothing
    declares (which Snakemake would not own, and so would not regenerate).
    """
    declares_object = "objects/" in _de_patterns_block()
    writes_object = "saveRDS" in (
        WORKFLOW / "scripts" / "03_differential" / "estimand.R"
    ).read_text(encoding="utf-8")
    assert declares_object == writes_object, (
        f"DE_PATTERNS declares an objects/ output: {declares_object}; "
        f"estimand.R calls saveRDS: {writes_object}"
    )


def test_the_canonical_fit_is_the_family_fit():
    """Only ``family_fit`` persists a DESeqDataSet; consumers read that one.

    One fit per family is the point of the family architecture -- estimands
    extract from a shared dispersion and covariance -- so a per-estimand copy of
    the same object would be both redundant and a second source of truth.
    """
    assert "saveRDS(dds, file.path(dirs$objects, \"deseq2.rds\"))" in (
        WORKFLOW / "scripts" / "03_differential" / "family_fit.R"
    ).read_text(encoding="utf-8")
    advanced = (WORKFLOW / "rules" / "advanced.smk").read_text(encoding="utf-8")
    assert 'analysis("de", "objects/deseq2.rds")' not in advanced
