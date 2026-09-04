"""Regression test for the numerator-minus-denominator sign convention.

The engine's single most load-bearing scientific invariant is that positive
log2_fold_change = numerator - denominator, single-sourced through
resolve_contrast() in utils.R. This executing test runs the family path
(materialize_inputs → family_fit → estimand) on a fixture with unambiguous
up- and down-in-numerator genes and asserts the signs + direction labels are correct.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"
MATERIALIZE_INPUTS = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "01_inputs" / "materialize_inputs.py"
FAMILY_FIT_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "family_fit.R"
ESTIMAND_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "estimand.R"


def test_de_sign_convention_numerator_minus_denominator(tmp_path):
    """A gene higher in the numerator gets log2_fold_change > 0 and
    direction == "up_in_numerator"; a gene lower in the numerator gets
    log2_fold_change < 0 and direction == "down_in_numerator".

    The minimal template already encodes this: treatment_a has ~6x counts vs
    control for genes g001-g016 (clear UP) and ~0.2x counts for genes g017-g024
    (clear DOWN). This test freezes that invariant with an executing check.
    """
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")

    # Check if DESeq2 is available (same pattern as other R-executing tests).
    try:
        subprocess.run(
            ["Rscript", "-e", "library(DESeq2)"],
            check=True, capture_output=True, timeout=120,
        )
    except subprocess.TimeoutExpired:
        pytest.skip("DESeq2 probe timed out")
    except subprocess.CalledProcessError:
        pytest.skip("DESeq2 not available")

    # Copy the minimal template to a temp directory.
    project_dir = tmp_path / "project"
    shutil.copytree(TEMPLATE, project_dir)

    # Run materialize_inputs to generate families/estimands/cell_weights.
    inputs_dir = tmp_path / "inputs"
    subprocess.run(
        [
            "python", str(MATERIALIZE_INPUTS),
            "--project-config", str(project_dir / "project.yaml"),
            "--counts", str(inputs_dir / "counts.tsv"),
            "--samples", str(inputs_dir / "samples.tsv"),
            "--annotation", str(inputs_dir / "annotation.tsv"),
            "--contrasts", str(inputs_dir / "contrasts.tsv"),
            "--families", str(inputs_dir / "families.tsv"),
            "--estimands", str(inputs_dir / "estimands.tsv"),
            "--cell-weights", str(inputs_dir / "estimand_cell_weights.tsv"),
            "--term-tests", str(inputs_dir / "term_tests.tsv"),
            "--manifest", str(inputs_dir / "input_manifest.json"),
            "--threads", "1",
        ],
        check=True, capture_output=True, text=True,
    )

    # Run family_fit.R for the treatment_a_vs_control family.
    family_dir = tmp_path / "families" / "treatment_a_vs_control"
    subprocess.run(
        [
            "Rscript", "--vanilla", str(FAMILY_FIT_R),
            "--project-config", str(project_dir / "project.yaml"),
            "--counts", str(inputs_dir / "counts.tsv"),
            "--samples", str(inputs_dir / "samples.tsv"),
            "--families", str(inputs_dir / "families.tsv"),
            "--estimands", str(inputs_dir / "estimands.tsv"),
            "--cell-weights", str(inputs_dir / "estimand_cell_weights.tsv"),
            "--family-id", "treatment_a_vs_control",
            "--outdir", str(family_dir),
        ],
        check=True, capture_output=True, text=True,
    )

    # Run estimand.R to extract the estimand.
    outdir = tmp_path / "out"
    subprocess.run(
        [
            "Rscript", "--vanilla", str(ESTIMAND_R),
            "--project-config", str(project_dir / "project.yaml"),
            "--counts", str(inputs_dir / "counts.tsv"),
            "--samples", str(inputs_dir / "samples.tsv"),
            "--annotation", str(inputs_dir / "annotation.tsv"),
            "--families", str(inputs_dir / "families.tsv"),
            "--estimands", str(inputs_dir / "estimands.tsv"),
            "--cell-weights", str(inputs_dir / "estimand_cell_weights.tsv"),
            "--family-dir", str(family_dir),
            "--estimand-id", "treatment_a_vs_control",
            "--outdir", str(outdir),
        ],
        check=True, capture_output=True, text=True,
    )

    # Read the DE results table.
    de_results = outdir / "tables" / "de_results.tsv"
    assert de_results.exists(), "de_results.tsv not produced"

    with de_results.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        results_by_gene = {row["gene_id"]: row for row in reader}

    # Assert sign convention for a gene clearly HIGHER in numerator (treatment_a).
    # g001 has ~550 counts in treatment_a vs ~87 in control (6x up).
    up_gene = results_by_gene["g001"]
    assert float(up_gene["log2_fold_change"]) > 0, (
        f"Gene g001 is higher in numerator (treatment_a) but got "
        f"log2_fold_change={up_gene['log2_fold_change']} <= 0"
    )
    assert up_gene["direction"] == "up_in_numerator", (
        f"Gene g001 is higher in numerator but got direction={up_gene['direction']}"
    )

    # Assert sign convention for a gene clearly LOWER in numerator (treatment_a).
    # g017 has ~119 counts in treatment_a vs ~520 in control (0.23x down).
    down_gene = results_by_gene["g017"]
    assert float(down_gene["log2_fold_change"]) < 0, (
        f"Gene g017 is lower in numerator (treatment_a) but got "
        f"log2_fold_change={down_gene['log2_fold_change']} >= 0"
    )
    assert down_gene["direction"] == "down_in_numerator", (
        f"Gene g017 is lower in numerator but got direction={down_gene['direction']}"
    )
