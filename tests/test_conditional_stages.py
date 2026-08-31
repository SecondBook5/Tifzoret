"""Executing tests for conditional analysis stages.

These stages are opt-in or conditional on specific config flags/columns:
- de_confirm.R: edgeR confirmatory DE (requires de_confirm flag or module enabled)
- spia.R: SPIA pathway topology (requires spia flag/module and KEGG provider)
- batch.R: batch-corrected ordination (requires batch column in samples.tsv)

Each test runs the stage via subprocess with the smallest fixture that triggers it,
then asserts a real output property to verify the stage executed correctly.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"
DE_CONFIRM_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "de_confirm.R"
SPIA_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "04_enrichment" / "spia.R"
BATCH_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "02_qc" / "batch.R"


def _has_rscript_package(package: str) -> bool:
    """Check if an R package is available to system Rscript."""
    if shutil.which("Rscript") is None:
        return False
    try:
        subprocess.run(
            ["Rscript", "-e", f"library({package})"],
            check=True, capture_output=True, timeout=5,
        )
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return False


def test_de_confirm_concordance_with_edger(tmp_path):
    """de_confirm.R runs edgeR quasi-likelihood DE on the same contrast as DESeq2,
    producing a concordance table that classifies genes by whether both engines,
    one only, or neither called them significant.

    This test runs de.R first to produce the primary DESeq2 results, then runs
    de_confirm.R and asserts that:
    - The edger_results.tsv is produced with log2 fold-change and FDR columns
    - The concordance.tsv is produced with concordance classes
    - At least one gene is classified as "both" (significant in both engines)
    """
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")

    # de_confirm requires both DESeq2 (for de.R) and edgeR.
    for package in ("DESeq2", "edgeR"):
        if not _has_rscript_package(package):
            pytest.skip(f"{package} not available")

    # Copy the minimal template.
    project_dir = tmp_path / "project"
    shutil.copytree(TEMPLATE, project_dir)

    # Run de.R first to produce the primary DESeq2 results that de_confirm reads.
    de_r = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "de.R"
    de_outdir = tmp_path / "de"
    subprocess.run(
        [
            "Rscript", "--vanilla", str(de_r),
            "--project-config", str(project_dir / "project.yaml"),
            "--counts", str(project_dir / "counts.tsv"),
            "--samples", str(project_dir / "samples.tsv"),
            "--annotation", str(project_dir / "annotation.tsv"),
            "--contrasts", str(project_dir / "contrasts.tsv"),
            "--contrast-id", "treatment_a_vs_control",
            "--outdir", str(de_outdir),
        ],
        check=True, capture_output=True, text=True,
    )

    # Run de_confirm.R on the same contrast.
    confirm_outdir = tmp_path / "de_confirm"
    subprocess.run(
        [
            "Rscript", "--vanilla", str(DE_CONFIRM_R),
            "--project-config", str(project_dir / "project.yaml"),
            "--counts", str(project_dir / "counts.tsv"),
            "--samples", str(project_dir / "samples.tsv"),
            "--annotation", str(project_dir / "annotation.tsv"),
            "--contrasts", str(project_dir / "contrasts.tsv"),
            "--contrast-id", "treatment_a_vs_control",
            "--de", str(de_outdir / "tables" / "de_results.tsv"),
            "--outdir", str(confirm_outdir),
        ],
        check=True, capture_output=True, text=True,
    )

    # Assert edgeR results table exists with expected columns.
    edger_results = confirm_outdir / "tables" / "edger_results.tsv"
    assert edger_results.exists(), "edger_results.tsv not produced"
    with edger_results.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        edger_rows = list(reader)
    assert len(edger_rows) > 0, "edger_results.tsv is empty"
    # Check that edgeR produces log2 fold-change and FDR columns.
    first_row = edger_rows[0]
    assert "gene_id" in first_row
    assert "edger_log2_fold_change" in first_row
    assert "edger_fdr" in first_row
    assert float(first_row["edger_log2_fold_change"]) != 0.0, "edgeR log2FC should be non-zero for DE genes"

    # Assert concordance table exists with concordance classes.
    concordance = confirm_outdir / "tables" / "de_concordance_displayed.tsv"
    assert concordance.exists(), "de_concordance_displayed.tsv not produced"
    with concordance.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        concordance_rows = list(reader)
    assert len(concordance_rows) > 0, "de_concordance_displayed.tsv is empty"
    # Check that concordance classification was performed.
    concordance_classes = {row["concordance_class"] for row in concordance_rows}
    # The minimal template has clear up/down genes, so at least some should be
    # called by both engines (concordance_class == "both").
    assert "both" in concordance_classes, (
        "Expected at least some genes to be significant in both engines "
        f"(concordance classes: {concordance_classes})"
    )


def test_spia_pathway_topology_analysis(tmp_path):
    """spia.R performs topology-aware pathway perturbation analysis using SPIA.

    SPIA is opt-in and requires SPIA/graphite packages plus species-matched org.db.
    It may gracefully skip if dependencies are unavailable. This test runs de.R
    first, then runs spia.R and asserts that either:
    - Real SPIA results are produced (spia_pathways.tsv with pathways, summary.json
      with method="SPIA"), OR
    - The stage gracefully skipped (empty tables, summary.json with skipped=true).

    The test does NOT fail if SPIA is unavailable (consistent with the stage's
    graceful degradation), but it DOES assert that the stage ran to completion
    and produced well-formed outputs in either case.
    """
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")

    # SPIA requires DESeq2 for de.R upstream.
    if not _has_rscript_package("DESeq2"):
        pytest.skip("DESeq2 not available")

    # Copy the minimal template and configure it for mouse (SPIA needs a real species).
    project_dir = tmp_path / "project"
    shutil.copytree(TEMPLATE, project_dir)
    config_path = project_dir / "project.yaml"
    config = yaml.safe_load(config_path.read_text())
    # Set species to mouse so SPIA can attempt to use org.Mm.eg.db.
    config["species"] = {
        "provider": "mouse",
        "scientific_name": "Mus musculus",
        "taxonomy_id": 10090,
    }
    config["reference"] = {
        "genome_build": "GRCm39",
        "annotation_release": 107,
    }
    # Enable SPIA module and KEGG provider.
    config["analysis"]["modules"] = {"spia": True}
    config["resources"] = {
        "cache": "~/.cache/tifzoret/resources",
        "offline": False,
        "refresh": False,
        "gene_sets": {"gmt": "gene_sets.gmt"},
        "providers": {"kegg": True},
    }
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    # Run de.R first.
    de_r = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "de.R"
    de_outdir = tmp_path / "de"
    subprocess.run(
        [
            "Rscript", "--vanilla", str(de_r),
            "--project-config", str(config_path),
            "--counts", str(project_dir / "counts.tsv"),
            "--samples", str(project_dir / "samples.tsv"),
            "--annotation", str(project_dir / "annotation.tsv"),
            "--contrasts", str(project_dir / "contrasts.tsv"),
            "--contrast-id", "treatment_a_vs_control",
            "--outdir", str(de_outdir),
        ],
        check=True, capture_output=True, text=True,
    )

    # Run spia.R.
    spia_outdir = tmp_path / "spia"
    subprocess.run(
        [
            "Rscript", "--vanilla", str(SPIA_R),
            "--project-config", str(config_path),
            "--de", str(de_outdir / "tables" / "de_results.tsv"),
            "--contrast-id", "treatment_a_vs_control",
            "--outdir", str(spia_outdir),
        ],
        check=True, capture_output=True, text=True,
    )

    # Assert that the stage produced well-formed outputs (either real results or
    # a graceful skip).
    pathways_table = spia_outdir / "tables" / "spia_pathways.tsv"
    summary_json = spia_outdir / "spia_summary.json"
    assert pathways_table.exists(), "spia_pathways.tsv not produced"
    assert summary_json.exists(), "spia_summary.json not produced"

    # Check the summary to see if SPIA actually ran or gracefully skipped.
    summary = json.loads(summary_json.read_text())
    assert "method" in summary
    assert "skipped" in summary
    # The summary must record either real results (skipped=false, method="SPIA")
    # or a graceful skip (skipped=true, reason explaining why).
    if summary["skipped"]:
        assert "reason" in summary, "Skipped SPIA must record a reason"
        # Pathways table should be empty but well-formed (header only).
        with pathways_table.open() as handle:
            lines = handle.readlines()
        assert len(lines) == 1, "Skipped SPIA should emit header-only pathways table"
    else:
        # Real SPIA run: method should be "SPIA", and pathways count should be >= 0.
        assert "SPIA" in summary["method"]
        assert "pathways" in summary
        # If pathways were found, the table should have rows.
        with pathways_table.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            pathway_rows = list(reader)
        # Can't assert pathway_rows is non-empty because the synthetic fixture's
        # genes may not map to real KEGG pathways, but we can assert that if
        # the summary says N pathways, the table has N rows.
        assert len(pathway_rows) == summary["pathways"]


def test_batch_corrected_pca_and_distance(tmp_path):
    """batch.R produces batch-corrected PCA and sample-distance views when a batch
    column is configured in analysis.batch.

    This test creates a fixture with a batch column in samples.tsv, runs qc.R to
    produce the VST expression (which batch.R requires), then runs batch.R and
    asserts that:
    - The batch_corrected_expression.tsv is produced
    - The PCA comparison figure (before/after batch correction) is produced
    - The batch_summary.json records that correction was attempted
    """
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")

    # batch.R requires DESeq2 (for qc.R upstream) and limma (for removeBatchEffect).
    for package in ("DESeq2", "limma"):
        if not _has_rscript_package(package):
            pytest.skip(f"{package} not available")

    # Copy the minimal template and add a batch column to samples.tsv.
    project_dir = tmp_path / "project"
    shutil.copytree(TEMPLATE, project_dir)

    samples_path = project_dir / "samples.tsv"
    samples = []
    with samples_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            # Assign batch labels: alternate between "batch1" and "batch2".
            sample_idx = len(samples)
            row["batch"] = "batch1" if sample_idx % 2 == 0 else "batch2"
            samples.append(row)

    # Write the updated samples.tsv with the batch column.
    with samples_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(samples[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(samples)

    # Update project.yaml to configure analysis.batch.
    config_path = project_dir / "project.yaml"
    config = yaml.safe_load(config_path.read_text())
    config["analysis"]["batch"] = "batch"
    config["analysis"]["modules"] = {"batch": True}
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    # Run qc.R first to produce the VST expression that batch.R requires.
    qc_r = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "02_qc" / "qc.R"
    qc_outdir = tmp_path / "qc"
    subprocess.run(
        [
            "Rscript", "--vanilla", str(qc_r),
            "--project-config", str(config_path),
            "--counts", str(project_dir / "counts.tsv"),
            "--samples", str(samples_path),
            "--annotation", str(project_dir / "annotation.tsv"),
            "--outdir", str(qc_outdir),
        ],
        check=True, capture_output=True, text=True,
    )

    # Run batch.R.
    batch_outdir = tmp_path / "batch"
    subprocess.run(
        [
            "Rscript", "--vanilla", str(BATCH_R),
            "--project-config", str(config_path),
            "--vst", str(qc_outdir / "objects" / "vst.rds"),
            "--samples", str(samples_path),
            "--outdir", str(batch_outdir),
        ],
        check=True, capture_output=True, text=True,
    )

    # Assert batch-corrected expression table is produced.
    corrected_expr = batch_outdir / "tables" / "batch_corrected_expression.tsv"
    assert corrected_expr.exists(), "batch_corrected_expression.tsv not produced"
    with corrected_expr.open() as handle:
        lines = handle.readlines()
    # Should have header + gene rows.
    assert len(lines) > 1, "batch_corrected_expression.tsv should have gene rows"

    # Assert PCA comparison figure is produced (before/after batch correction).
    pca_pdf = batch_outdir / "figures" / "batch_pca.pdf"
    pca_png = batch_outdir / "figures" / "batch_pca.png"
    assert pca_pdf.exists(), "batch_pca.pdf not produced"
    assert pca_png.exists(), "batch_pca.png not produced"
    assert pca_pdf.stat().st_size > 0, "batch_pca.pdf is empty"

    # Assert batch summary records that correction was applied.
    summary_json = batch_outdir / "batch_summary.json"
    assert summary_json.exists(), "batch_summary.json not produced"
    summary = json.loads(summary_json.read_text())
    assert "batch_variable" in summary
    assert summary["batch_variable"] == "batch"
    # The summary should record whether correction was applied (can be false if
    # batch is confounded with the biological group, but for this fixture with
    # alternating batches it should succeed).
    assert "correction_applied" in summary
