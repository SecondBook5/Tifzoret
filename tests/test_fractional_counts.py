"""Test that read_counts_contract() handles fractional counts correctly.

The guard detects non-integer counts, warns, and rounds them (matching
DESeqDataSetFromTximport behavior). Integer inputs are byte-identical to
the current behavior (no warning, no change).
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_integer_counts_are_unchanged():
    """All-integer counts pass through unchanged with no warning."""
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")

    with tempfile.TemporaryDirectory() as tmpdir:
        counts_file = Path(tmpdir) / "counts.tsv"
        counts_file.write_text(
            "gene_id\tsample1\tsample2\n"
            "g001\t100\t200\n"
            "g002\t50\t75\n"
        )

        script = f"""
        source("{ROOT / 'src' / 'tifzoret' / 'workflow' / 'scripts' / 'utils.R'}")
        mat <- read_counts_contract("{counts_file}")
        stopifnot(mat["g001", "sample1"] == 100)
        stopifnot(mat["g001", "sample2"] == 200)
        stopifnot(mat["g002", "sample1"] == 50)
        stopifnot(mat["g002", "sample2"] == 75)
        """

        result = subprocess.run(
            ["Rscript", "--vanilla", "-e", script],
            capture_output=True,
            text=True,
            timeout=10,
        )

        # Should succeed with no warnings.
        assert result.returncode == 0, f"stdout: {result.stdout}\nstderr: {result.stderr}"
        assert "warning" not in result.stderr.lower(), "Integer counts should not warn"


def test_fractional_counts_are_rounded_with_warning():
    """Fractional counts trigger a warning and are rounded to nearest integer."""
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")

    with tempfile.TemporaryDirectory() as tmpdir:
        counts_file = Path(tmpdir) / "counts.tsv"
        counts_file.write_text(
            "gene_id\tsample1\tsample2\n"
            "g001\t100.7\t200.3\n"
            "g002\t50.2\t75.8\n"
        )

        script = f"""
        source("{ROOT / 'src' / 'tifzoret' / 'workflow' / 'scripts' / 'utils.R'}")
        mat <- read_counts_contract("{counts_file}")
        stopifnot(mat["g001", "sample1"] == 101)  # 100.7 rounds to 101
        stopifnot(mat["g001", "sample2"] == 200)  # 200.3 rounds to 200
        stopifnot(mat["g002", "sample1"] == 50)   # 50.2 rounds to 50
        stopifnot(mat["g002", "sample2"] == 76)   # 75.8 rounds to 76
        """

        result = subprocess.run(
            ["Rscript", "--vanilla", "-e", script],
            capture_output=True,
            text=True,
            timeout=10,
        )

        # Should succeed with a warning about non-integer counts.
        assert result.returncode == 0, f"stdout: {result.stdout}\nstderr: {result.stderr}"
        assert "non-integer counts detected" in result.stderr, (
            "Fractional counts should trigger a warning"
        )
        assert "rounding to nearest integer" in result.stderr
