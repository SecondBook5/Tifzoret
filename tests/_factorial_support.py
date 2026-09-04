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


def require_r(*packages: str) -> None:
    """pytest.skip unless Rscript and every named R package are available."""
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")

    for package in packages:
        try:
            subprocess.run(
                ["Rscript", "-e", f"library({package})"],
                check=True,
                capture_output=True,
                timeout=5,
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            pytest.skip(f"{package} not available")


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
