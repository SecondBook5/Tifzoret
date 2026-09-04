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
