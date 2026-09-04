"""analysis.families flows through load_project (spec §4, §8, §9)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
import yaml

from tifzoret.config import ProjectValidationError, load_project

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"


def _project(tmp_path: Path) -> Path:
    destination = tmp_path / "project"
    shutil.copytree(TEMPLATE, destination)
    return destination / "project.yaml"


def _write(config_path: Path, mutate) -> Path:
    data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    mutate(data)
    config_path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return config_path


def _samples_columns(config_path: Path) -> list[str]:
    samples = config_path.parent / "samples.tsv"
    return samples.read_text(encoding="utf-8").splitlines()[0].split("\t")


def test_a_project_without_families_still_loads(tmp_path):
    project = load_project(_project(tmp_path))
    assert project.families  # desugared from contrasts.tsv
    assert all(family.declared is False for family in project.families)


def test_pairwise_contrasts_desugar_into_a_family(tmp_path):
    project = load_project(_project(tmp_path))
    ids = {estimand.id for family in project.families for estimand in family.estimands}
    contrast_ids = {row["contrast_id"] for row in project.contrast_rows}
    assert ids <= contrast_ids


def test_contrast_rows_are_not_duplicated_by_desugaring(tmp_path):
    project = load_project(_project(tmp_path))
    ids = [row["contrast_id"] for row in project.contrast_rows]
    assert len(ids) == len(set(ids))


def test_a_declared_family_appears_in_contrast_rows_as_an_estimand(tmp_path):
    config_path = _project(tmp_path)
    columns = _samples_columns(config_path)
    factor = "condition" if "condition" in columns else columns[1]
    levels = sorted(
        {
            line.split("\t")[columns.index(factor)]
            for line in (config_path.parent / "samples.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[1:]
            if line.strip()
        }
    )
    assert len(levels) >= 2

    def mutate(data):
        data["analysis"]["families"] = {
            "fam_a": {
                "design": f"~ {factor}",
                "cells": [factor],
                "estimands": [
                    {
                        "id": "est_declared",
                        "label": "Declared estimand",
                        "role": "primary",
                        "expression": f"({levels[1]}) - ({levels[0]})",
                    }
                ],
            }
        }

    project = load_project(_write(config_path, mutate))
    row = next(r for r in project.contrast_rows if r["contrast_id"] == "est_declared")
    assert row["type"] == "estimand"
    assert row["family"] == "fam_a"
    assert row["role"] == "primary"
    assert row["numerator"] == "Declared estimand"


def test_estimand_rows_and_cell_weights_are_exported(tmp_path):
    config_path = _project(tmp_path)
    columns = _samples_columns(config_path)
    factor = "condition" if "condition" in columns else columns[1]
    levels = sorted(
        {
            line.split("\t")[columns.index(factor)]
            for line in (config_path.parent / "samples.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[1:]
            if line.strip()
        }
    )

    def mutate(data):
        data["analysis"]["families"] = {
            "fam_a": {
                "design": f"~ {factor}",
                "cells": [factor],
                "estimands": [
                    {
                        "id": "est_declared",
                        "role": "primary",
                        "expression": f"({levels[1]}) - ({levels[0]})",
                    }
                ],
            }
        }

    project = load_project(_write(config_path, mutate))
    weights = [r for r in project.cell_weight_rows if r["estimand_id"] == "est_declared"]
    assert {r["cell_key"]: float(r["weight"]) for r in weights} == {
        levels[1]: 1.0,
        levels[0]: -1.0,
    }
    families = {r["family_id"]: r for r in project.family_rows}
    assert families["fam_a"]["shrinkage"] == "apeglm"
    assert families["fam_a"]["filter"] == "design_aware"
    assert families["fam_a"]["declared"] == "true"


def test_term_tests_are_derived_with_real_degrees_of_freedom(tmp_path):
    config_path = _project(tmp_path)
    columns = _samples_columns(config_path)
    factor = "condition" if "condition" in columns else columns[1]
    levels = sorted(
        {
            line.split("\t")[columns.index(factor)]
            for line in (config_path.parent / "samples.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[1:]
            if line.strip()
        }
    )

    def mutate(data):
        data["analysis"]["families"] = {
            "fam_a": {
                "design": f"~ {factor}",
                "cells": [factor],
                "term_tests": "auto",
                "estimands": [
                    {
                        "id": "est_declared",
                        "role": "primary",
                        "expression": f"({levels[1]}) - ({levels[0]})",
                    }
                ],
            }
        }

    project = load_project(_write(config_path, mutate))
    rows = [r for r in project.term_test_rows if r["family_id"] == "fam_a"]
    assert rows
    assert all(int(r["df"]) >= 1 for r in rows)


def test_an_invalid_expression_is_a_load_error(tmp_path):
    config_path = _project(tmp_path)
    columns = _samples_columns(config_path)
    factor = "condition" if "condition" in columns else columns[1]

    def mutate(data):
        data["analysis"]["families"] = {
            "fam_a": {
                "design": f"~ {factor}",
                "cells": [factor],
                "estimands": [
                    {"id": "est_bad", "role": "primary", "expression": "(nonexistent_level)"}
                ],
            }
        }

    with pytest.raises(ProjectValidationError):
        load_project(_write(config_path, mutate))


def test_an_estimand_id_colliding_with_a_contrast_id_is_a_load_error(tmp_path):
    config_path = _project(tmp_path)
    columns = _samples_columns(config_path)
    factor = "condition" if "condition" in columns else columns[1]
    levels = sorted(
        {
            line.split("\t")[columns.index(factor)]
            for line in (config_path.parent / "samples.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[1:]
            if line.strip()
        }
    )
    existing = load_project(config_path).contrast_rows[0]["contrast_id"]

    def mutate(data):
        data["analysis"]["families"] = {
            "fam_a": {
                "design": f"~ {factor}",
                "cells": [factor],
                "estimands": [
                    {
                        "id": existing,
                        "role": "primary",
                        "expression": f"({levels[1]}) - ({levels[0]})",
                    }
                ],
            }
        }

    with pytest.raises(ProjectValidationError, match="collides"):
        load_project(_write(config_path, mutate))


def test_an_unknown_analysis_key_is_still_rejected(tmp_path):
    config_path = _project(tmp_path)

    def mutate(data):
        data["analysis"]["not_a_real_key"] = True

    with pytest.raises(ProjectValidationError):
        load_project(_write(config_path, mutate))
