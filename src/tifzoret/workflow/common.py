"""Shared utilities for workflow stage scripts.

This module consolidates functions duplicated across 12+ workflow scripts.
Import these instead of defining them locally in each script.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def read_tsv(path: Path, *, comments: bool = False) -> tuple[list[str], list[dict[str, str]]]:
    """Read a tab-separated file into header list and row dicts.

    Args:
        path: Path to TSV file
        comments: If True, skip lines starting with '#'

    Returns:
        Tuple of (fieldnames, rows as list of dicts)
    """
    with path.open(newline="", encoding="utf-8") as handle:
        lines = (line for line in handle if not comments or not line.startswith("#"))
        reader = csv.DictReader(lines, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    """Write rows to a tab-separated file.

    Args:
        path: Output path
        fieldnames: Column headers
        rows: List of dicts (one per row)
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    """Return hex SHA-256 digest of file at path.

    Reads in 1MB blocks to handle large files efficiently.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_workflow_args(*required_args: str) -> dict[str, str]:
    """Standard argparse setup for workflow scripts.

    Args:
        *required_args: Names of required arguments (without -- prefix)

    Returns:
        Dict mapping argument names to their values

    Example:
        args = parse_workflow_args("project-config", "outdir")
        config = load_project(args["project-config"])
    """
    parser = argparse.ArgumentParser()
    for arg in required_args:
        parser.add_argument(f"--{arg}", required=True)
    parsed = parser.parse_args()
    return {arg: getattr(parsed, arg.replace("-", "_")) for arg in required_args}


def ensure_output_dirs(outdir: Path) -> dict[str, Path]:
    """Create standard output subdirectories.

    Args:
        outdir: Base output directory

    Returns:
        Dict with keys 'tables', 'figures', 'objects' pointing to created subdirs
    """
    outdir = Path(outdir)
    subdirs = {
        "tables": outdir / "tables",
        "figures": outdir / "figures",
        "objects": outdir / "objects",
    }
    for subdir in subdirs.values():
        subdir.mkdir(parents=True, exist_ok=True)
    return subdirs
