"""Shared helpers for factorial fixture tests."""

from __future__ import annotations

import csv
import importlib.util
import shutil
import subprocess
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT / "tests" / "fixtures" / "factorial" / "generate.py"

# Cache R package availability checks at module level to avoid repeated expensive probes
_R_PACKAGE_CACHE: dict[str, str | None] = {}  # package -> skip_reason (None if available)


def require_r(*packages: str) -> None:
    """pytest.skip unless Rscript and every named R package are available."""
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")

    for package in packages:
        # Check cache first
        if package not in _R_PACKAGE_CACHE:
            try:
                subprocess.run(
                    ["Rscript", "-e", f"library({package})"],
                    check=True,
                    capture_output=True,
                    timeout=120,
                )
                _R_PACKAGE_CACHE[package] = None  # Available
            except subprocess.TimeoutExpired:
                _R_PACKAGE_CACHE[package] = f"{package} probe timed out"
                pytest.skip(_R_PACKAGE_CACHE[package])
            except subprocess.CalledProcessError:
                _R_PACKAGE_CACHE[package] = f"{package} not available"
                pytest.skip(_R_PACKAGE_CACHE[package])

        # Use cached result
        skip_reason = _R_PACKAGE_CACHE[package]
        if skip_reason is not None:
            pytest.skip(skip_reason)


def fixtures() -> ModuleType:
    """Import and return the fixture generator module (tests/fixtures/factorial/generate.py)."""
    spec = importlib.util.spec_from_file_location("factorial_fixtures", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def build_fixture(name: str, destination: Path) -> Path:
    """Build one scenario; returns the destination directory."""
    module = fixtures()
    return module.build(name, destination)


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    """Read a TSV into a list of row dicts."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# Shared estimand definitions for factorial fixture tests
ESTIMANDS = {
    "est_interaction": {"a2|b2": 1.0, "a2|b1": -1.0, "a1|b2": -1.0, "a1|b1": 1.0},
    "est_arm_a1": {"a1|b2": 1.0, "a1|b1": -1.0},
    "est_arm_a2": {"a2|b2": 1.0, "a2|b1": -1.0},
    "est_composite": {"a2|b2": 1.0, "a1|b2": -1.0},
}


def stage_family(
    tmp_path: Path,
    script: Path,
    template: Path,
    scenario: str = "positive_interaction",
    design: str = "~ factor_a * factor_b",
    cells: str = "factor_a|factor_b",
    gene_filter: str = "design_aware",
) -> Path:
    """
    Stage a family_fit.R invocation with the specified scenario and design.

    Returns the output directory.
    """
    import shutil
    import subprocess

    import yaml

    fixture = build_fixture(scenario, tmp_path / "fx")
    project = tmp_path / "project"
    shutil.copytree(template, project)
    config = yaml.safe_load((project / "project.yaml").read_text(encoding="utf-8"))
    (project / "project.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
    )
    (tmp_path / "families.tsv").write_text(
        "family_id\tdesign\tcells\treference_levels\treplicate_unit\tshrinkage\tfilter\tdeclared\n"
        f"fam_a\t{design}\t{cells}\t\t\tapeglm\t{gene_filter}\ttrue\n",
        encoding="utf-8",
    )
    lines = ["family_id\testimand_id\tlabel\trole\texpression\tatom_kind"]
    weights = ["family_id\testimand_id\tcell_key\tweight"]
    for index, (name, cell_weights) in enumerate(ESTIMANDS.items()):
        role = "primary" if index == 0 else ""
        lines.append(f"fam_a\t{name}\t{name}\t{role}\tsource\tcell")
        for key, weight in cell_weights.items():
            weights.append(f"fam_a\t{name}\t{key}\t{weight!r}")
    (tmp_path / "estimands.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (tmp_path / "estimand_cell_weights.tsv").write_text(
        "\n".join(weights) + "\n", encoding="utf-8"
    )
    outdir = tmp_path / "out"
    result = subprocess.run(
        [
            "Rscript",
            "--vanilla",
            str(script),
            "--project-config",
            str(project / "project.yaml"),
            "--counts",
            str(fixture / "counts.tsv"),
            "--samples",
            str(fixture / "samples.tsv"),
            "--families",
            str(tmp_path / "families.tsv"),
            "--estimands",
            str(tmp_path / "estimands.tsv"),
            "--cell-weights",
            str(tmp_path / "estimand_cell_weights.tsv"),
            "--family-id",
            "fam_a",
            "--outdir",
            str(outdir),
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    return outdir
