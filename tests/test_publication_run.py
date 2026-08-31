"""End-to-end execution test for the publication profile.

The publication profile is the flagship scaffold: it exercises composition.R,
regulators.R, and publication.R (plus qc, de, pathways, ontology, hypotheses,
report) in one integrated run. This test runs the scaffold to completion and
asserts that key artifacts from composition, regulators, and publication exist.

This is the executing counterpart to test_workflow_dag.py's publication dry-run
test — it verifies the R stages actually complete successfully, not just that
the DAG assembles.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "publication"


def test_publication_profile_executes_to_completion(tmp_path):
    """The publication scaffold runs end-to-end offline, producing the report,
    manifest, and representative artifacts from composition, regulators, and
    publication modules.

    This executing test closes much of the "15 R stages have no executing test"
    gap: it transitively exercises composition.R, regulators.R, and publication.R
    plus all their dependencies (qc, de, pathways, ontology, hypotheses, report).
    """
    # Skip if R/conda is unavailable (same guard as other executing tests).
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")

    # Check if required R packages are available. The publication profile needs
    # DESeq2 (for de.R), viper (for regulators.R), and fgsea (for pathways.R).
    for package in ("DESeq2", "viper", "fgsea"):
        try:
            subprocess.run(
                ["Rscript", "-e", f"library({package})"],
                check=True, capture_output=True, timeout=5,
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            pytest.skip(f"{package} not available")

    # Copy the publication template to a temp directory.
    project_dir = tmp_path / "project"
    shutil.copytree(TEMPLATE, project_dir)

    # Run the publication profile end-to-end with --no-conda.
    # Use --cores 2 to allow some parallelism without overwhelming CI.
    result = subprocess.run(
        [
            "tifzoret", "run", str(project_dir / "project.yaml"),
            "--cores", "2", "--no-conda",
        ],
        capture_output=True,
        text=True,
    )

    # Assert the run succeeded.
    assert result.returncode == 0, (
        f"tifzoret run failed:\nSTDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
    )

    # Assert the core front-door artifacts exist.
    results = project_dir / "results" / "publication_demo" / "all"
    assert (results / "REPORT.html").exists(), "REPORT.html missing"
    assert (results / "manifest.json").exists(), "manifest.json missing"

    # Assert representative artifacts from the three key R stages.
    contrast_dir = results / "contrasts" / "treated_vs_control" / "analyses"

    # Composition: cell-state signature scores.
    composition_table = contrast_dir / "composition" / "tables" / "cell_state_scores.tsv"
    assert composition_table.exists(), f"composition artifact missing: {composition_table}"

    # Regulators: summary artifact (always produced).
    regulators_summary = contrast_dir / "regulators" / "regulators_summary.json"
    assert regulators_summary.exists(), f"regulators artifact missing: {regulators_summary}"

    # Publication: program definitions table.
    publication_table = contrast_dir / "publication" / "tables" / "program_definitions.tsv"
    assert publication_table.exists(), f"publication artifact missing: {publication_table}"
