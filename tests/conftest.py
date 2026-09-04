"""Shared pytest fixtures for Tifzoret tests."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from collections.abc import Iterator


ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ROOT / "src" / "tifzoret" / "templates"


@pytest.fixture
def template_dir() -> Path:
    """Return the path to Tifzoret's template directory."""
    return TEMPLATES


@pytest.fixture
def minimal_template() -> Path:
    """Return path to the minimal template (smallest end-to-end example)."""
    return TEMPLATES / "minimal"


@pytest.fixture
def factorial_template() -> Path:
    """Return path to the factorial template (2×2 interaction design)."""
    return TEMPLATES / "factorial"


@pytest.fixture
def publication_template() -> Path:
    """Return path to the publication template (flagship scaffold)."""
    return TEMPLATES / "publication"


@pytest.fixture
def temp_project(tmp_path: Path, minimal_template: Path) -> Iterator[Path]:
    """Stage a complete minimal project in a temporary directory.

    Returns the project directory path with project.yaml ready to use.
    Automatically cleaned up after the test completes.
    """
    project = tmp_path / "project"
    shutil.copytree(minimal_template, project)
    yield project


@pytest.fixture
def source_root() -> Path:
    """Return the Tifzoret repository root."""
    return ROOT


@pytest.fixture
def workflow_scripts() -> Path:
    """Return path to workflow scripts directory."""
    return ROOT / "src" / "tifzoret" / "workflow" / "scripts"
