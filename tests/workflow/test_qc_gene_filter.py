"""The QC variance-stabilising transform must not see genes that carry no information.

``qc.R`` used to build its DESeqDataSet from every row of the count matrix. On a
whole-transcriptome matrix a third to a half of the rows are zero or near-zero in
every sample (a GRCm39 matrix has 56,748 rows of which about 31,800 sum to fewer
than ten reads). Those rows are constant after the transform, they tie in the
ranks that ssGSEA and GSVA take over the same matrix, and they dilute every
principal component, so PC1 reported 55.5 % where the expressed genes alone give
58.7 %. The study-wide VST feeds PCA, sample distance, composition scores and
pathway scores, so the unfiltered transform skewed all of them.

Library metrics are a statement about the raw sequencing, so they must still be
computed on the unfiltered matrix.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"
QC_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "02_qc" / "qc.R"

ZERO_GENES = 12
LOW_GENES = 7  # summed count 9 across the nine samples: below the threshold of 10
MIN_TOTAL_COUNT = 10


def _rscript_can_load(package: str) -> bool:
    probe = subprocess.run(
        ["Rscript", "--vanilla", "-e", f"quit(status = !requireNamespace('{package}', quietly = TRUE))"],
        capture_output=True,
        text=True,
    )
    return probe.returncode == 0


def _append_uninformative_genes(project_dir: Path) -> None:
    counts_path = project_dir / "counts.tsv"
    annotation_path = project_dir / "annotation.tsv"
    with counts_path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    width = len(rows[0]) - 1
    for index in range(ZERO_GENES):
        rows.append([f"zero{index:02d}"] + ["0"] * width)
    for index in range(LOW_GENES):
        # nine reads in total: one per sample, so detected in every sample yet below 10
        rows.append([f"low{index:02d}"] + ["1"] * width)
    with counts_path.open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)

    with annotation_path.open(newline="") as handle:
        annotation = list(csv.reader(handle, delimiter="\t"))
    header = annotation[0]
    for prefix, count in (("zero", ZERO_GENES), ("low", LOW_GENES)):
        for index in range(count):
            gene = f"{prefix}{index:02d}"
            row = [gene if column == "gene_id" else (gene if column == "gene_symbol" else "") for column in header]
            annotation.append(row)
    with annotation_path.open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(annotation)


@pytest.fixture(scope="module")
def qc_outputs(tmp_path_factory) -> Path:
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    if not _rscript_can_load("DESeq2"):
        pytest.skip("DESeq2 not available")
    root = tmp_path_factory.mktemp("qc_filter")
    project_dir = root / "project"
    shutil.copytree(TEMPLATE, project_dir)
    _append_uninformative_genes(project_dir)
    config_path = project_dir / "project.yaml"
    outdir = root / "qc"
    subprocess.run(
        [
            "Rscript",
            "--vanilla",
            str(QC_R),
            "--project-config",
            str(config_path),
            "--counts",
            str(project_dir / "counts.tsv"),
            "--samples",
            str(project_dir / "samples.tsv"),
            "--annotation",
            str(project_dir / "annotation.tsv"),
            "--outdir",
            str(outdir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return outdir


def _expected_universe(qc_dir: Path) -> tuple[list[str], list[str]]:
    """Independent oracle: split the written count matrix at the documented threshold."""
    with (qc_dir.parent / "project" / "counts.tsv").open(newline="") as handle:
        table = list(csv.DictReader(handle, delimiter="\t"))
    samples = [name for name in table[0] if name != "gene_id"]
    kept, dropped = [], []
    for row in table:
        total = sum(float(row[sample]) for sample in samples)
        (kept if total >= MIN_TOTAL_COUNT else dropped).append(row["gene_id"])
    return kept, dropped


def _vst_gene_ids(qc_dir: Path) -> list[str]:
    command = (
        "suppressPackageStartupMessages(library(SummarizedExperiment));"
        f"v <- readRDS('{qc_dir / 'objects' / 'vst.rds'}');"
        "writeLines(rownames(v))"
    )
    result = subprocess.run(["Rscript", "--vanilla", "-e", command], check=True, capture_output=True, text=True)
    return result.stdout.split()


def test_the_vst_object_holds_only_genes_with_enough_reads(qc_outputs):
    kept, dropped = _expected_universe(qc_outputs)
    assert dropped, "the fixture must contain genes below the threshold"
    assert any(gene.startswith("zero") for gene in dropped)
    assert any(gene.startswith("low") for gene in dropped)
    assert sorted(_vst_gene_ids(qc_outputs)) == sorted(kept)


def test_library_metrics_still_describe_the_unfiltered_matrix(qc_outputs):
    with (qc_outputs / "tables" / "library_metrics.tsv").open(newline="") as handle:
        reported = {
            row["sample_id"]: int(float(row["detected_genes"])) for row in csv.DictReader(handle, delimiter="\t")
        }
    with (qc_outputs.parent / "project" / "counts.tsv").open(newline="") as handle:
        table = list(csv.DictReader(handle, delimiter="\t"))
    expected = {sample: sum(1 for row in table if float(row[sample]) > 0) for sample in table[0] if sample != "gene_id"}
    assert reported == expected
    # the 7 near-empty genes carry one read in every sample, so they are *detected*
    # even though the transform excludes them
    assert all(value >= LOW_GENES for value in reported.values())


def test_the_summary_records_how_many_genes_the_transform_excluded(qc_outputs):
    summary = json.loads((qc_outputs / "qc_summary.json").read_text())
    kept, dropped = _expected_universe(qc_outputs)
    assert summary["genes"] == len(kept) + len(dropped)
    assert summary["genes_in_vst"] == len(kept)
    assert summary["genes_excluded_low_count"] == len(dropped)
    assert summary["qc_min_total_count"] == MIN_TOTAL_COUNT


def test_the_reported_pca_is_the_pca_of_the_stored_vst(qc_outputs):
    command = (
        "suppressPackageStartupMessages(library(SummarizedExperiment));"
        f"v <- readRDS('{qc_outputs / 'objects' / 'vst.rds'}');"
        "p <- prcomp(t(assay(v)), center = TRUE, scale. = FALSE);"
        "cat(100 * p$sdev[1]^2 / sum(p$sdev^2))"
    )
    result = subprocess.run(["Rscript", "--vanilla", "-e", command], check=True, capture_output=True, text=True)
    independent = float(result.stdout)
    summary = json.loads((qc_outputs / "qc_summary.json").read_text())
    assert summary["pca_variance_pc1"] == pytest.approx(independent, abs=1e-3)
