from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "10_figures" / "manifest.py"
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"


def _load_manifest_module():
    spec = importlib.util.spec_from_file_location("_manifest_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_analysis_r_version_prefers_recorded_over_path(tmp_path, monkeypatch):
    """The R that ran the analysis is recorded in the QC summary; the manifest
    must report it rather than whatever ``R`` happens to be on the system PATH
    (the manifest step runs in the core env, which has no R)."""
    module = _load_manifest_module()
    monkeypatch.setattr(module, "command_output", lambda command: "R version 0.0.0 (path fallback)")
    qc = tmp_path / "qc"
    qc.mkdir()
    (qc / "qc_summary.json").write_text(
        json.dumps({"r_version": "R version 4.5.3 (2026-03-11)"}), encoding="utf-8"
    )
    assert module.analysis_r_version(tmp_path) == "R version 4.5.3 (2026-03-11)"


def test_analysis_r_version_falls_back_when_unrecorded(tmp_path, monkeypatch):
    module = _load_manifest_module()
    monkeypatch.setattr(module, "command_output", lambda command: "R version 0.0.0 (path fallback)")
    # No qc_summary.json at all.
    assert module.analysis_r_version(tmp_path) == "R version 0.0.0 (path fallback)"
    # Present but missing/blank the field → also falls back.
    qc = tmp_path / "qc"
    qc.mkdir()
    (qc / "qc_summary.json").write_text(json.dumps({"r_version": "  "}), encoding="utf-8")
    assert module.analysis_r_version(tmp_path) == "R version 0.0.0 (path fallback)"


def test_manifest_expands_environment_variables_in_input_paths(tmp_path):
    (tmp_path / "samples.tsv").write_text(
        "sample_id\tbam\tcondition\n"
        "control_1\tcontrol_1.bam\tcontrol\n"
        "control_2\tcontrol_2.bam\tcontrol\n"
        "treated_1\ttreated_1.bam\ttreated\n"
        "treated_2\ttreated_2.bam\ttreated\n",
        encoding="utf-8",
    )
    (tmp_path / "contrasts.tsv").write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\n"
        "treated_vs_control\tcondition\ttreated\tcontrol\n",
        encoding="utf-8",
    )
    (tmp_path / "sets.gmt").write_text("fixture\tfixture\tGene1\tGene2\n", encoding="utf-8")
    (tmp_path / "genes.gtf").write_text(
        'chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "gene1"; gene_name "Gene1";\n',
        encoding="utf-8",
    )
    (tmp_path / "control_1.bam").touch()
    (tmp_path / "control_2.bam").touch()
    (tmp_path / "treated_1.bam").touch()
    (tmp_path / "treated_2.bam").touch()

    config = {
        "version": 2,
        "project": {"id": "manifest_fixture", "title": "Manifest fixture"},
        "species": {"provider": "mouse", "scientific_name": "Mus musculus", "taxonomy_id": 10090},
        "reference": {"genome_build": "GRCm39", "annotation_release": 107},
        "inputs": {
            "kind": "bam",
            "bam_root": ".",
            "samples": "samples.tsv",
            "gtf": "${BULK_RNA_FRAME_TEST_GTF}",
        },
        "counting": {
            "threads": 1, "feature_type": "exon", "attribute": "gene_id",
            "paired_end": False, "count_read_pairs": False,
            "require_both_ends_aligned": False, "exclude_chimeric_fragments": False,
            "strandedness": "unstranded", "strand_test_modes": [0],
            "strand_min_dominance": 0.8,
        },
        "analysis": {
            "design": "~ condition", "contrasts": "contrasts.tsv", "profile": "standard",
            "modules": {"de": False, "pathways": False, "ontology": False, "report": False},
        },
        "resources": {
            "gene_sets": {"gmt": "sets.gmt", "min_size": 2, "max_size": 100},
        },
        "figures": {
            "group": "condition", "palette": {"control": "#A6CEE3", "treated": "#F4A6A6"},
            "pca": {"ellipse_level": 0.8},
            "de": {"fdr": 0.05, "abs_log2fc": 1, "top_labels": 5, "top_heatmap_genes": 10, "z_limit": 1.5},
            "pathways": {"top_ora_terms": 5, "top_gsva_terms": 5, "gsea_curves_per_direction": 1, "seed": 1},
        },
        "output": {"root": "results"},
    }
    config_path = tmp_path / "project.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    results = tmp_path / "results"
    results.mkdir()
    (results / "inputs").mkdir()
    prepared_gtf_checksum = "a" * 64
    (results / "inputs" / "input_manifest.json").write_text(
        json.dumps(
            {
                "source": {
                    "gtf": {
                        "path": str(tmp_path / "genes.gtf"),
                        "bytes": (tmp_path / "genes.gtf").stat().st_size,
                        "sha256": prepared_gtf_checksum,
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    (results / "result.tsv").write_text("value\n1\n", encoding="utf-8")
    output = results / "manifest.json"

    environment = os.environ.copy()
    environment["BULK_RNA_FRAME_TEST_GTF"] = str(tmp_path / "genes.gtf")
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--project-config",
            str(config_path),
            "--results",
            str(results),
            "--output",
            str(output),
        ],
        check=True,
        env=environment,
    )

    manifest = json.loads(output.read_text(encoding="utf-8"))
    assert manifest["schema_version"] == 2
    assert manifest["contrast_semantics"] == "all signed effects are numerator minus denominator"
    input_paths = {record["path"] for record in manifest["inputs"]}
    assert str(tmp_path / "genes.gtf") in input_paths
    gtf_record = next(record for record in manifest["inputs"] if record["path"] == str(tmp_path / "genes.gtf"))
    assert gtf_record["sha256"] == prepared_gtf_checksum


def test_manifest_checksums_binding_prior_edges(tmp_path):
    """The unsigned binding-prior regulon is a declared analysis input and MUST be
    checksummed in the manifest -- previously it was silently omitted, leaving the
    second regulator view's source data unpinned. Uses the minimal template with
    the regulators module and a DoRothEA primary regulon so the binding prior is
    the only edge file on disk, and runs the manifest step directly (pure Python;
    no R or Snakemake needed)."""
    destination = tmp_path / "project"
    shutil.copytree(TEMPLATE, destination)
    config_path = destination / "project.yaml"
    data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    data["species"] = {"provider": "mouse", "scientific_name": "Mus musculus", "taxonomy_id": 10090}
    data["reference"] = {"genome_build": "GRCm39", "annotation_release": 107}
    data["analysis"]["modules"] = {"regulators": True}
    data["resources"]["providers"] = {"dorothea": True, "gtrd": True}
    data["resources"]["binding_prior_edges"] = "binding.tsv"
    (destination / "binding.tsv").write_text("source\ttarget\nTF1\tGene1\n", encoding="utf-8")
    config_path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")

    results = destination / "results"
    (results / "regulators").mkdir(parents=True)
    (results / "regulators" / "result.tsv").write_text("value\n1\n", encoding="utf-8")
    output = results / "manifest.json"

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--project-config",
            str(config_path),
            "--results",
            str(results),
            "--output",
            str(output),
        ],
        check=True,
    )

    manifest = json.loads(output.read_text(encoding="utf-8"))
    input_paths = {record["path"] for record in manifest["inputs"]}
    binding_path = str((destination / "binding.tsv").resolve())
    assert binding_path in input_paths, sorted(input_paths)
    binding_record = next(record for record in manifest["inputs"] if record["path"] == binding_path)
    assert len(binding_record["sha256"]) == 64


def _run_manifest(config_path: Path, results: Path) -> subprocess.CompletedProcess:
    output = results / "manifest.json"
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--project-config",
            str(config_path),
            "--results",
            str(results),
            "--output",
            str(output),
        ],
        check=False,
        capture_output=True,
        text=True,
    )


def _strict_manifest_project(tmp_path: Path, strict: bool, *, degradation: bool = True) -> tuple[Path, Path]:
    """Minimal-template project (optionally strict) with one stage summary that
    carries a warning, ready for the terminal manifest gate. ``degradation`` picks
    the severity: a degradation-tagged object (which blocks a strict build) versus
    a bare-string standing caveat (which is recorded but never blocks)."""
    destination = tmp_path / "project"
    shutil.copytree(TEMPLATE, destination)
    config_path = destination / "project.yaml"
    data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if strict:
        data["execution"] = {"strict": True}
    config_path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    results = destination / "results"
    (results / "qc").mkdir(parents=True)
    warning = (
        {"message": "deterministic proxy used -- NOT canonical scoring", "severity": "degradation"}
        if degradation
        else "STRING edges are association evidence, not causal edges."
    )
    (results / "qc" / "qc_summary.json").write_text(
        json.dumps({"warnings": [warning]}),
        encoding="utf-8",
    )
    return config_path, results


def test_strict_mode_fails_on_a_degradation_warning(tmp_path):
    """In strict (publication) mode the terminal manifest step fails if any stage
    recorded a DEGRADATION-severity warning, so a publication run cannot report
    success on a degraded result. The manifest is still written (for inspection)
    before failing."""
    config_path, results = _strict_manifest_project(tmp_path, strict=True, degradation=True)
    proc = _run_manifest(config_path, results)
    assert proc.returncode != 0, proc.stdout
    assert "strict mode" in proc.stderr
    assert "degradation" in proc.stderr
    manifest = json.loads((results / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["strict"] is True
    assert manifest["warnings"]  # recorded on disk despite the failure
    assert manifest["warnings"][0]["severity"] == "degradation"


def test_strict_mode_tolerates_standing_caveats(tmp_path):
    """A standing scientific CAVEAT (a bare-string interpretation note a clean run
    always emits, e.g. the STRING association disclaimer) is recorded but must NOT
    fail a strict build -- otherwise strict mode would be unusable for any real
    study, all of which emit such caveats. When a run's only warnings are caveats,
    strict passes."""
    config_path, results = _strict_manifest_project(tmp_path, strict=True, degradation=False)
    proc = _run_manifest(config_path, results)
    assert proc.returncode == 0, proc.stderr
    manifest = json.loads((results / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["strict"] is True
    assert manifest["warnings"]  # the caveat is recorded as provenance
    assert manifest["warnings"][0]["severity"] == "caveat"


def test_lenient_mode_tolerates_stage_warnings(tmp_path):
    """Absent execution.strict, even a degradation warning is recorded but the run
    succeeds -- byte-for-byte the historical behavior."""
    config_path, results = _strict_manifest_project(tmp_path, strict=False, degradation=True)
    proc = _run_manifest(config_path, results)
    assert proc.returncode == 0, proc.stderr
    manifest = json.loads((results / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["strict"] is False
    assert manifest["warnings"]


def test_strict_mode_succeeds_without_warnings(tmp_path):
    """Strict mode is not punitive on a clean run: with no recorded warnings the
    terminal gate passes."""
    destination = tmp_path / "project"
    shutil.copytree(TEMPLATE, destination)
    config_path = destination / "project.yaml"
    data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    data["execution"] = {"strict": True}
    config_path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    results = destination / "results"
    (results / "qc").mkdir(parents=True)
    (results / "qc" / "qc_summary.json").write_text(json.dumps({"warnings": []}), encoding="utf-8")
    proc = _run_manifest(config_path, results)
    assert proc.returncode == 0, proc.stderr
