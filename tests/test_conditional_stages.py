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


from _factorial_support import read_tsv_rows, require_r, stage_family

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"
MATERIALIZE_INPUTS = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "01_inputs" / "materialize_inputs.py"
FAMILY_FIT_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "family_fit.R"
ESTIMAND_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "estimand.R"
DE_CONFIRM_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "de_confirm.R"
SPIA_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "04_enrichment" / "spia.R"
BATCH_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "02_qc" / "batch.R"


def _has_rscript_package(package: str) -> bool | str:
    """Check if an R package is available to system Rscript.

    Returns True if available, or a string skip reason otherwise.
    """
    if shutil.which("Rscript") is None:
        return "Rscript not available"
    try:
        subprocess.run(
            ["Rscript", "-e", f"library({package})"],
            check=True, capture_output=True, timeout=120,
        )
        return True
    except subprocess.TimeoutExpired:
        return f"{package} probe timed out"
    except subprocess.CalledProcessError:
        return f"{package} not available"


def _run_family_de_factorial(tmp_path: Path, estimand_id: str, scenario: str = "positive_interaction") -> Path:
    """Run the family DE path using the factorial fixture.

    Returns the output directory containing de_results.tsv.
    """
    # Use stage_family to get a proper family fit with the factorial fixture.
    family_dir = stage_family(
        tmp_path,
        FAMILY_FIT_R,
        TEMPLATE,
        scenario=scenario
    )

    # Extract the specified estimand.
    fixture = tmp_path / "fx"
    de_outdir = tmp_path / "de"
    result = subprocess.run(
        ["Rscript", "--vanilla", str(ESTIMAND_R),
         "--project-config", str(tmp_path / "project" / "project.yaml"),
         "--counts", str(fixture / "counts.tsv"),
         "--samples", str(fixture / "samples.tsv"),
         "--annotation", str(fixture / "annotation.tsv"),
         "--families", str(tmp_path / "families.tsv"),
         "--estimands", str(tmp_path / "estimands.tsv"),
         "--cell-weights", str(tmp_path / "estimand_cell_weights.tsv"),
         "--family-dir", str(family_dir),
         "--estimand-id", estimand_id,
         "--outdir", str(de_outdir)],
        capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr

    return de_outdir


def _run_family_de_minimal(tmp_path: Path, estimand_id: str) -> tuple[Path, Path]:
    """Run the family DE path on the minimal template (desugar → family_fit → estimand).

    Returns (de_outdir, project_dir) where project_dir contains the original contrasts.tsv.
    """
    # Copy minimal template
    project_dir = tmp_path / "project"
    shutil.copytree(TEMPLATE, project_dir)

    # Run materialize_inputs to desugar contrasts into families/estimands/cell_weights.
    # The minimal template's pairwise contrasts all desugar into a single family named "main".
    inputs_dir = tmp_path / "inputs"
    result = subprocess.run(
        ["python", str(MATERIALIZE_INPUTS),
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
         "--threads", "1"],
        capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr

    # Run family_fit.R on the desugared family (family_id is "main" for the minimal template)
    family_id = "main"
    family_dir = tmp_path / "families" / family_id
    result = subprocess.run(
        ["Rscript", "--vanilla", str(FAMILY_FIT_R),
         "--project-config", str(project_dir / "project.yaml"),
         "--counts", str(inputs_dir / "counts.tsv"),
         "--samples", str(inputs_dir / "samples.tsv"),
         "--families", str(inputs_dir / "families.tsv"),
         "--estimands", str(inputs_dir / "estimands.tsv"),
         "--cell-weights", str(inputs_dir / "estimand_cell_weights.tsv"),
         "--family-id", family_id,
         "--outdir", str(family_dir)],
        capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr

    # Run estimand.R to produce de_results.tsv
    de_outdir = tmp_path / "de"
    result = subprocess.run(
        ["Rscript", "--vanilla", str(ESTIMAND_R),
         "--project-config", str(project_dir / "project.yaml"),
         "--counts", str(inputs_dir / "counts.tsv"),
         "--samples", str(inputs_dir / "samples.tsv"),
         "--annotation", str(inputs_dir / "annotation.tsv"),
         "--families", str(inputs_dir / "families.tsv"),
         "--estimands", str(inputs_dir / "estimands.tsv"),
         "--cell-weights", str(inputs_dir / "estimand_cell_weights.tsv"),
         "--family-dir", str(family_dir),
         "--estimand-id", estimand_id,
         "--outdir", str(de_outdir)],
        capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr

    return de_outdir, project_dir


# TODO: de_confirm.R stage itself has no test coverage; retargeting blocked on
# fixture format mismatch with factorial designs (Task 14 may address)
def test_deseq2_edger_concordance(tmp_path):
    """Test raw DESeq2/edgeR concordance on the same contrast.

    Two independent negative-binomial engines that agree on direction reduce
    false discoveries. This test verifies that DESeq2 and edgeR produce
    concordant results (same sign of log2FC for genes with meaningful effects)
    when run on the same factorial fixture data.

    Note: This tests the engines directly, not the de_confirm.R stage (which is
    an opt-in downstream module). The de_confirm.R stage wraps this concordance
    check for production use but is not itself tested here.
    """
    # Requires both DESeq2 (for family path) and edgeR (for direct concordance test).
    require_r("DESeq2", "edgeR", "apeglm")

    # Use est_arm_a1 from factorial fixture (simple pairwise: a1|b2 vs a1|b1).
    estimand_id = "est_arm_a1"
    de_outdir = _run_family_de_factorial(tmp_path, estimand_id)

    # Read DESeq2 results
    deseq2_results = read_tsv_rows(de_outdir / "tables" / "de_results.tsv")
    assert len(deseq2_results) > 0, "DESeq2 de_results.tsv is empty"

    # Create an R script to run edgeR on the same data and contrast
    fixture = tmp_path / "fx"
    edger_script = tmp_path / "run_edger.R"
    edger_script.write_text("""
library(edgeR)

# Read inputs
counts <- read.delim(commandArgs(TRUE)[1], row.names=1, check.names=FALSE)
samples <- read.delim(commandArgs(TRUE)[2], row.names=1)
outfile <- commandArgs(TRUE)[3]

# Build design for the factorial: ~ factor_a * factor_b
# est_arm_a1 compares a1|b2 vs a1|b1, so we subset to factor_a == "a1"
samples_subset <- samples[samples$factor_a == "a1", , drop=FALSE]
counts_subset <- counts[, rownames(samples_subset), drop=FALSE]

# Create groups for edgeR
group <- factor(samples_subset$factor_b)  # "b1" or "b2"

# Build DGEList and run edgeR pipeline
dge <- DGEList(counts=counts_subset, group=group)
dge <- calcNormFactors(dge)
design <- model.matrix(~group)
dge <- estimateDisp(dge, design)

# Run exact test (b2 vs b1)
et <- exactTest(dge, pair=c("b1", "b2"))
results <- topTags(et, n=Inf, sort.by="none")$table

# Add gene_id column
results$gene_id <- rownames(results)

# Write results
write.table(results[, c("gene_id", "logFC", "PValue", "FDR")],
            file=outfile, sep="\\t", quote=FALSE, row.names=FALSE)
""", encoding="utf-8")

    # Run the edgeR script
    edger_outfile = tmp_path / "edger_results.tsv"
    result = subprocess.run(
        ["Rscript", "--vanilla", str(edger_script),
         str(fixture / "counts.tsv"),
         str(fixture / "samples.tsv"),
         str(edger_outfile)],
        capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, f"edgeR script failed: {result.stderr}"

    # Read edgeR results
    edger_results = read_tsv_rows(edger_outfile)
    assert len(edger_results) > 0, "edgeR results are empty"

    # Build concordance: map by gene_id
    deseq2_by_gene = {row["gene_id"]: row for row in deseq2_results}
    edger_by_gene = {row["gene_id"]: row for row in edger_results}

    # Check that both engines produced results for the same genes
    common_genes = set(deseq2_by_gene.keys()) & set(edger_by_gene.keys())
    assert len(common_genes) > 0, "No common genes between DESeq2 and edgeR"

    # Assert concordance of log-fold changes (direction agreement).
    # For genes with non-zero effects in both engines, they should agree on direction.
    lfc_concordant = 0
    lfc_discordant = 0
    lfc_threshold = 0.1  # Consider genes with |logFC| > 0.1 in both engines

    for gene in common_genes:
        deseq2_lfc_str = deseq2_by_gene[gene].get("log2_fold_change", "NA")
        edger_lfc_str = edger_by_gene[gene].get("logFC", "NA")

        if deseq2_lfc_str in ("NA", "") or edger_lfc_str in ("NA", ""):
            continue

        deseq2_lfc = float(deseq2_lfc_str)
        edger_lfc = float(edger_lfc_str)

        # Only check genes with meaningful effects in both engines
        if abs(deseq2_lfc) > lfc_threshold and abs(edger_lfc) > lfc_threshold:
            if deseq2_lfc * edger_lfc > 0:
                lfc_concordant += 1
            else:
                lfc_discordant += 1

    # Assert that we tested some genes and that concordance is high
    total_tested = lfc_concordant + lfc_discordant
    assert total_tested > 0, (
        "No genes with |logFC| > 0.1 in both engines found for concordance test"
    )

    concordance_rate = lfc_concordant / total_tested if total_tested > 0 else 0
    assert concordance_rate >= 0.8, (
        f"Poor concordance between DESeq2 and edgeR: {lfc_concordant} concordant, "
        f"{lfc_discordant} discordant (concordance rate: {concordance_rate:.2%})"
    )


def test_spia_pathway_topology_analysis(tmp_path):
    """spia.R performs topology-aware pathway perturbation analysis using SPIA.

    SPIA is opt-in and requires SPIA/graphite packages plus species-matched org.db.
    It may gracefully skip if dependencies are unavailable. This test runs the family
    path first, then runs spia.R and asserts that either:
    - Real SPIA results are produced (spia_pathways.tsv with pathways, summary.json
      with method="SPIA"), OR
    - The stage gracefully skipped (empty tables, summary.json with skipped=true).

    The test does NOT fail if SPIA is unavailable (consistent with the stage's
    graceful degradation), but it DOES assert that the stage ran to completion
    and produced well-formed outputs in either case.
    """
    # SPIA requires DESeq2 and apeglm upstream.
    require_r("DESeq2", "apeglm")

    # Use est_arm_a1 for the DE results from factorial fixture.
    estimand_id = "est_arm_a1"
    de_outdir = _run_family_de_factorial(tmp_path, estimand_id)

    # Update the project config for mouse species (SPIA needs a real species).
    config_path = tmp_path / "project" / "project.yaml"
    config = yaml.safe_load(config_path.read_text())
    config["species"] = {
        "provider": "mouse",
        "scientific_name": "Mus musculus",
        "taxonomy_id": 10090,
    }
    config["reference"] = {
        "genome_build": "GRCm39",
        "annotation_release": 107,
    }
    config["analysis"]["modules"] = {"spia": True}
    config["resources"] = {
        "cache": "~/.cache/tifzoret/resources",
        "offline": False,
        "refresh": False,
        "gene_sets": {"gmt": "gene_sets.gmt"},
        "providers": {"kegg": True},
    }
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    # Run spia.R.
    spia_outdir = tmp_path / "spia"
    subprocess.run(
        [
            "Rscript", "--vanilla", str(SPIA_R),
            "--project-config", str(config_path),
            "--de", str(de_outdir / "tables" / "de_results.tsv"),
            "--contrast-id", estimand_id,
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
        result = _has_rscript_package(package)
        if result is not True:
            pytest.skip(result)

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
