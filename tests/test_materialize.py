from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "bulk_rna_frame" / "workflow" / "scripts" / "materialize_inputs.py"
FIXTURES = ROOT / "tests" / "fixtures"
MINIMAL = ROOT / "src" / "bulk_rna_frame" / "templates" / "minimal"


def test_count_adapter_materializes_the_canonical_contract(tmp_path):
    output = tmp_path / "canonical"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--project-config",
            str(MINIMAL / "project.yaml"),
            "--counts",
            str(output / "counts.tsv"),
            "--samples",
            str(output / "samples.tsv"),
            "--annotation",
            str(output / "annotation.tsv"),
            "--manifest",
            str(output / "input_manifest.json"),
        ],
        check=True,
    )
    manifest = json.loads((output / "input_manifest.json").read_text())
    assert manifest["source"]["kind"] == "counts"
    assert manifest["source"]["selected_samples"] == [
        "control_1",
        "control_2",
        "control_3",
        "treatment_a_1",
        "treatment_a_2",
        "treatment_a_3",
        "treatment_b_1",
        "treatment_b_2",
        "treatment_b_3",
    ]
    assert (output / "counts.tsv").is_file()
    assert (output / "samples.tsv").is_file()
    assert (output / "annotation.tsv").is_file()


@pytest.mark.skipif(
    shutil.which("samtools") is None or shutil.which("featureCounts") is None,
    reason="BAM acceptance test requires samtools and featureCounts",
)
@pytest.mark.parametrize("kind", ["bam", "nfcore_rnaseq"])
def test_bam_adapters_run_featurecounts(tmp_path, kind):
    bam_root = tmp_path / "upstream"
    bam_root.mkdir()
    for sample_id in ("control_1", "treated_1"):
        subprocess.run(
            [
                "samtools",
                "view",
                "-b",
                "-o",
                str(bam_root / f"{sample_id}.sorted.bam"),
                str(FIXTURES / "tiny.sam"),
            ],
            check=True,
        )

    samples = tmp_path / "samples.tsv"
    with samples.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "bam", "condition"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            [
                {"sample_id": "control_1", "bam": "control_1.sorted.bam", "condition": "control"},
                {"sample_id": "treated_1", "bam": "treated_1.sorted.bam", "condition": "treated"},
            ]
        )
    (tmp_path / "contrasts.tsv").write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\n"
        "treated_vs_control\tcondition\ttreated\tcontrol\n"
    )
    (tmp_path / "sets.gmt").write_text("fixture\tfixture\tGeneOne\tUnusedGene\n")
    inputs = {
        "kind": kind,
        "samples": "samples.tsv",
        "gtf": str(FIXTURES / "genes.gtf"),
    }
    if kind == "bam":
        inputs["bam_root"] = "upstream"
    else:
        inputs["root"] = "upstream"
        inputs["bam_pattern"] = "{sample_id}.sorted.bam"
    config = {
        "version": 1,
        "project": {"id": "bam_fixture", "title": "BAM fixture"},
        "inputs": inputs,
        "counting": {
            "threads": 1,
            "feature_type": "exon",
            "attribute": "gene_id",
            "paired_end": False,
            "count_read_pairs": False,
            "require_both_ends_aligned": False,
            "exclude_chimeric_fragments": False,
            "strandedness": "infer",
            "strand_test_modes": [0, 1, 2],
            "strand_min_dominance": 0.8,
        },
        "design": {"formula": "~ condition"},
        "contrasts": "contrasts.tsv",
        "gene_sets": {"gmt": "sets.gmt", "min_size": 2, "max_size": 100},
        "modules": {"qc": True, "de": False, "pathways": False},
        "figures": {
            "group": "condition",
            "palette": {"control": "#A6CEE3", "treated": "#F4A6A6"},
            "pca": {"ellipse_level": 0.8},
            "de": {
                "fdr": 0.05,
                "abs_log2fc": 1.0,
                "top_labels": 10,
                "top_heatmap_genes": 20,
                "z_limit": 1.5,
            },
            "pathways": {
                "top_ora_terms": 8,
                "top_gsva_terms": 12,
                "gsea_curves_per_direction": 2,
                "seed": 1,
            },
        },
        "output": {"root": "results"},
    }
    config_path = tmp_path / "project.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    output = tmp_path / "canonical"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--project-config",
            str(config_path),
            "--counts",
            str(output / "counts.tsv"),
            "--samples",
            str(output / "samples.tsv"),
            "--annotation",
            str(output / "annotation.tsv"),
            "--manifest",
            str(output / "input_manifest.json"),
        ],
        check=True,
    )
    counts = list(csv.DictReader((output / "counts.tsv").open(), delimiter="\t"))
    manifest = json.loads((output / "input_manifest.json").read_text())
    assert counts == [{"gene_id": "gene1", "control_1": "1", "treated_1": "1"}]
    assert manifest["source"]["kind"] == kind
    assert manifest["source"]["counting"]["resolved_strand_mode"] in {0, 1, 2}
    assert manifest["source"]["bams"][0]["quickcheck"] is True
