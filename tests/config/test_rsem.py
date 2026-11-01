"""Exercise the RSEM boundary, normalization preservation, and failure gates."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from tifzoret.config import ProjectValidationError, load_project

ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "src/tifzoret/templates/rsem"
MATERIALIZE = ROOT / "src/tifzoret/workflow/scripts/01_inputs/materialize_inputs.py"
UTILS = ROOT / "src/tifzoret/workflow/scripts/utils.R"


@pytest.fixture
def study(tmp_path: Path) -> Path:
    """Copy the portable RSEM example so tests cannot mutate package templates."""
    destination = tmp_path / "study"
    shutil.copytree(TEMPLATE, destination)
    return destination / "project.yaml"


def prepare_boundary(study: Path) -> None:
    """Execute only the quantification-import boundary for a validated study."""
    project = load_project(study)
    out = project.result_root / "inputs"
    args = [sys.executable, str(MATERIALIZE), "--project-config", str(study)]
    for argument, filename in [
        ("counts", "counts.tsv"),
        ("samples", "samples.tsv"),
        ("annotation", "annotation.tsv"),
        ("contrasts", "contrasts.tsv"),
        ("families", "families.tsv"),
        ("estimands", "estimands.tsv"),
        ("cell-weights", "estimand_cell_weights.tsv"),
        ("term-tests", "term_tests.tsv"),
        ("manifest", "input_manifest.json"),
        ("tximport", "tximport.rds"),
    ]:
        args.extend(["--" + argument, str(out / filename)])
    subprocess.run(args, check=True, capture_output=True, text=True)


def test_rsem_validation_preserves_order_and_sources(study: Path) -> None:
    """Require a unique quantification source per canonical biological sample."""
    project = load_project(study)
    assert project.source_kind == "rsem"
    assert len(project.rsem_paths) == len(project.sample_rows)
    assert project.source_files == (*project.rsem_paths, project.annotation)
    assert [path.stem for path in project.rsem_paths] == [row["sample_id"] + ".genes" for row in project.sample_rows]


@pytest.mark.parametrize(
    "corruption", ["duplicate_path", "negative_count", "missing_column", "gene_order", "missing_file"]
)
def test_rsem_rejects_malformed_sources(study: Path, corruption: str) -> None:
    """Fail on source identity and numeric defects before any downstream work."""
    project = load_project(study)
    path = project.rsem_paths[0]
    if corruption == "duplicate_path":
        text = project.samples.read_text()
        text = text.replace(project.rsem_paths[1].name, path.name)
        project.samples.write_text(text)
    elif corruption == "missing_file":
        path.unlink()
    else:
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields, rows = reader.fieldnames, list(reader)
        if corruption == "negative_count":
            rows[0]["expected_count"] = "-1"
        elif corruption == "missing_column":
            fields.remove("effective_length")
        else:
            rows = rows[::-1]
        with path.open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    with pytest.raises(ProjectValidationError):
        load_project(study)


@pytest.mark.parametrize("module", ["de_confirm", "deconvolution"])
def test_count_only_optional_modules_fail_closed(study: Path, module: str) -> None:
    """Do not silently lose effective lengths through optional count-only paths."""
    cfg = yaml.safe_load(study.read_text())
    cfg["analysis"]["modules"] = {module: True}
    study.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ProjectValidationError, match="length-aware adapter"):
        load_project(study)


@pytest.mark.parametrize("mode", ["no"])
def test_tximport_normalization_matches_direct_package(study: Path, mode: str) -> None:
    """Compare shared QC/family constructors against direct package import.

    Deliberately varying effective lengths distinguish a genuine length-aware
    boundary from rounding counts. Estimate normalization only, without DE.
    """
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript unavailable")
    cfg = yaml.safe_load(study.read_text())
    cfg["inputs"]["counts_from_abundance"] = mode
    study.write_text(yaml.safe_dump(cfg))
    project = load_project(study)
    out = project.result_root / "inputs"
    prepare_boundary(study)
    check = f"""
    source({str(UTILS)!r})
    cfg <- read_project({str(study)!r})
    path <- {str(out / 'counts.tsv')!r}
    counts <- read_counts_contract(path)
    meta <- as.data.frame(readr::read_tsv({str(out / 'samples.tsv')!r}, show_col_types=FALSE))
    rownames(meta) <- meta$sample_id
    direct_samples <- read.delim({str(project.samples)!r})
    files <- setNames(file.path({str(study.parent)!r}, direct_samples$rsem), direct_samples$sample_id)
    direct <- tximport::tximport(files, type='rsem', txIn=FALSE, txOut=FALSE, countsFromAbundance={mode!r})
    actual <- make_deseq_dataset(cfg, path, counts, meta, ~ batch + condition)
    expected <- DESeq2::DESeqDataSetFromTximport(direct, meta, ~ batch + condition)
    actual <- DESeq2::estimateSizeFactors(actual)
    expected <- DESeq2::estimateSizeFactors(expected)
    stopifnot(isTRUE(all.equal(DESeq2::counts(actual, normalized=TRUE), DESeq2::counts(expected, normalized=TRUE))))
    if ({mode!r} == 'no') {{
      stopifnot('avgTxLength' %in% SummarizedExperiment::assayNames(actual))
      stopifnot(!is.null(DESeq2::normalizationFactors(actual)))
      plain <- DESeq2::estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(counts, meta, ~ batch + condition))
      stopifnot(!isTRUE(all.equal(DESeq2::counts(actual, normalized=TRUE), DESeq2::counts(plain, normalized=TRUE))))
    }} else {{
      stopifnot(!'avgTxLength' %in% SummarizedExperiment::assayNames(actual))
    }}
    subset <- make_deseq_dataset(cfg, path, counts[1:6,,drop=FALSE], meta, ~ batch + condition)
    stopifnot(identical(rownames(subset), rownames(counts)[1:6]))
    # A changed count artifact must not silently pass the length-aware guard.
    counts[1,1] <- counts[1,1] + 1L
    failed <- try(make_deseq_dataset(cfg, path, counts, meta, ~ batch + condition), silent=TRUE)
    stopifnot(inherits(failed, 'try-error'))
    """
    result = subprocess.run(["Rscript", "--vanilla", "-e", check], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


def test_rsem_dag_declares_bundle_dependency(study: Path) -> None:
    """QC and family fits must depend on the materialized quantification bundle."""
    if shutil.which("snakemake") is None:
        pytest.skip("Snakemake unavailable")
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
    result = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(ROOT / "src/tifzoret/workflow/Snakefile"),
            "--configfile",
            str(study),
            "--cores",
            "1",
            "--dry-run",
            str(load_project(study).result_root / "manifest.json"),
        ],
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "tximport.rds" in result.stdout
    assert "rsem_effective_lengths.tsv" in result.stdout


@pytest.mark.parametrize("mode", ["scaledTPM", "lengthScaledTPM"])
def test_gene_level_rsem_rejects_transcript_only_modes(study: Path, mode: str) -> None:
    """Gene-result files cannot provide tximport transcript-abundance scaling."""
    cfg = yaml.safe_load(study.read_text())
    cfg["inputs"]["counts_from_abundance"] = mode
    study.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ProjectValidationError):
        load_project(study)


@pytest.mark.parametrize("policy", ["error", "one"])
def test_zero_lengths_have_explicit_audited_policy(study: Path, policy: str) -> None:
    """Preserve original lengths and audit corrections rather than deleting genes."""
    cfg = yaml.safe_load(study.read_text())
    cfg["inputs"]["zero_length_policy"] = policy
    study.write_text(yaml.safe_dump(cfg))
    project = load_project(study)
    path = project.rsem_paths[0]
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields, rows = reader.fieldnames, list(reader)
    rows[0]["effective_length"] = "0"
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if policy == "error":
        with pytest.raises(subprocess.CalledProcessError):
            prepare_boundary(study)
        return
    prepare_boundary(study)
    out = project.result_root / "inputs"
    with (out / "rsem_length_adjustments.tsv").open() as handle:
        adjustments = list(csv.DictReader(handle, delimiter="\t"))
    assert len(adjustments) == 1
    assert adjustments[0]["original_length"] == "0"
    assert adjustments[0]["model_length"] == "1"
    with (out / "rsem_effective_lengths.tsv").open() as handle:
        raw = next(csv.DictReader(handle, delimiter="\t"))
    with (out / "rsem_model_lengths.tsv").open() as handle:
        model = next(csv.DictReader(handle, delimiter="\t"))
    sample = project.sample_rows[0]["sample_id"]
    assert raw[sample] == "0"
    assert model[sample] == "1"


@pytest.mark.parametrize("operator,retained", [("gt", False), ("ge", True)])
def test_count_threshold_comparator_is_exact(study: Path, operator: str, retained: bool) -> None:
    """A gene on the threshold must follow the declared strict/inclusive rule."""
    project = load_project(study)
    for path in project.rsem_paths:
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields, rows = reader.fieldnames, list(reader)
        rows[0]["expected_count"] = "1"
        with path.open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    cfg = yaml.safe_load(study.read_text())
    cfg["inputs"].update(
        min_total_count=len(project.sample_rows), count_threshold_operator=operator, exclude_gene_regex="^g002$"
    )
    study.write_text(yaml.safe_dump(cfg))
    prepare_boundary(study)
    with (project.result_root / "inputs/counts.tsv").open() as handle:
        ids = {row["gene_id"] for row in csv.DictReader(handle, delimiter="\t")}
    assert ("g001" in ids) == retained
    assert "g002" not in ids


def test_boundary_only_resources_and_explicit_families(study: Path) -> None:
    """Preparation needs no invented pathways or duplicate legacy contrasts."""
    cfg = yaml.safe_load(study.read_text())
    cfg["analysis"]["modules"] = {"pathways": False, "ontology": False, "de": False, "report": False}
    cfg["analysis"]["families"] = {
        "response": {
            "design": "~ condition",
            "cells": ["condition"],
            "filter": "none",
            "shrinkage": "none",
            "estimands": [{"id": "response", "role": "primary", "expression": "(treatment_a) - (control)"}],
        }
    }
    study.write_text(yaml.safe_dump(cfg))
    (study.parent / "gene_sets.gmt").write_text("")
    (study.parent / "contrasts.tsv").write_text("contrast_id\tfactor\tnumerator\tdenominator\n")
    project = load_project(study)
    assert project.contrast_rows[0]["contrast_id"] == "response"
    assert "de" not in project.modules


def test_provider_backed_study_allows_empty_custom_gmt(study: Path) -> None:
    """A real GO provider supplies resources without invented placeholder sets."""
    cfg = yaml.safe_load(study.read_text())
    cfg["resources"]["providers"] = {"go": True}
    (study.parent / "gene_sets.gmt").write_text("")
    study.write_text(yaml.safe_dump(cfg))
    assert load_project(study).config["resources"]["providers"]["go"] is True


def test_declared_pairwise_labels_survive_materialization(study: Path) -> None:
    """The enrichment seam must retain levels, type and family design."""
    cfg = yaml.safe_load(study.read_text())
    cfg["analysis"]["families"] = {
        "response": {
            "design": "~ batch + condition",
            "cells": ["condition"],
            "filter": "none",
            "shrinkage": "none",
            "estimands": [{"id": "response", "role": "primary", "expression": "(treatment_a) - (control)"}],
        }
    }
    (study.parent / "contrasts.tsv").write_text("contrast_id\tfactor\tnumerator\tdenominator\n")
    study.write_text(yaml.safe_dump(cfg))
    prepare_boundary(study)
    out = load_project(study).result_root / "inputs/contrasts.tsv"
    with out.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["type"] == "pairwise"
    assert row["factor"] == "condition"
    assert row["numerator"] == "treatment_a"
    assert row["denominator"] == "control"
