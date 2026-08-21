from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "manifest.py"


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
        "sample_id\tbam\tcondition\ncontrol_1\tcontrol.bam\tcontrol\n"
        "treated_1\ttreated.bam\ttreated\n",
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
    (tmp_path / "control.bam").touch()
    (tmp_path / "treated.bam").touch()

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
