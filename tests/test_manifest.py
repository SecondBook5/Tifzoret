from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "bulk_rna_frame" / "workflow" / "scripts" / "manifest.py"


def test_manifest_expands_environment_variables_in_input_paths(tmp_path):
    for name in ("samples.tsv", "contrasts.tsv", "sets.gmt", "genes.gtf"):
        (tmp_path / name).write_text(f"fixture {name}\n", encoding="utf-8")

    config = {
        "project": {"id": "manifest_fixture", "title": "Manifest fixture"},
        "inputs": {
            "kind": "bam",
            "samples": "samples.tsv",
            "gtf": "${BULK_RNA_FRAME_TEST_GTF}",
        },
        "contrasts": "contrasts.tsv",
        "gene_sets": {"gmt": "sets.gmt"},
    }
    config_path = tmp_path / "project.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    results = tmp_path / "results"
    results.mkdir()
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
    input_paths = {record["path"] for record in manifest["inputs"]}
    assert str(tmp_path / "genes.gtf") in input_paths
