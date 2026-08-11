from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from bulk_rna_frame.cli import main
import yaml


def test_init_creates_valid_project(tmp_path):
    destination = tmp_path / "new_project"
    assert main(["init", str(destination)]) == 0
    assert (destination / "project.yaml").is_file()
    assert main(["validate", str(destination / "project.yaml")]) == 0


def test_init_refuses_nonempty_destination(tmp_path):
    destination = tmp_path / "existing"
    destination.mkdir()
    (destination / "keep.txt").write_text("user data")
    assert main(["init", str(destination)]) == 2
    assert (destination / "keep.txt").read_text() == "user data"


def test_init_can_scaffold_each_bam_boundary(tmp_path):
    bam = tmp_path / "bam_project"
    nfcore = tmp_path / "nfcore_project"
    assert main(["init", str(bam), "--input", "bam"]) == 0
    assert main(["init", str(nfcore), "--input", "nfcore-rnaseq"]) == 0
    assert "kind: bam" in (bam / "project.yaml").read_text()
    assert "kind: nfcore_rnaseq" in (nfcore / "project.yaml").read_text()
    assert "bam_pattern:" in (nfcore / "project.yaml").read_text()


def test_prepare_targets_only_the_canonical_input_manifest(tmp_path, monkeypatch):
    destination = tmp_path / "project"
    assert main(["init", str(destination)]) == 0
    observed = {}

    def fake_run(command, **kwargs):
        observed["command"] = command
        observed["cwd"] = kwargs["cwd"]
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("bulk_rna_frame.cli.subprocess.run", fake_run)
    assert main(["prepare", str(destination / "project.yaml"), "--no-conda"]) == 0
    assert observed["command"][-1].endswith(
        "results/synthetic_demo/all/inputs/input_manifest.json"
    )
    assert "--dry-run" not in observed["command"]


def test_migrate_config_writes_valid_v2(tmp_path):
    destination = tmp_path / "project"
    assert main(["init", str(destination)]) == 0
    v2 = yaml.safe_load((destination / "project.yaml").read_text())
    v1 = {
        "version": 1, "project": v2["project"], "inputs": v2["inputs"],
        "design": {"formula": v2["analysis"]["design"]}, "contrasts": v2["analysis"]["contrasts"],
        "gene_sets": v2["resources"]["gene_sets"],
        "modules": {"qc": True, "de": True, "pathways": True},
        "figures": v2["figures"], "output": v2["output"],
    }
    source = destination / "v1.yaml"
    source.write_text(yaml.safe_dump(v1, sort_keys=False))
    output = destination / "v2.yaml"
    assert main(["migrate-config", str(source), "--output", str(output), "--species", "mouse", "--genome-build", "GRCm39"]) == 0
    assert yaml.safe_load(output.read_text())["version"] == 2


def test_verify_uses_project_result_as_candidate(tmp_path):
    destination = tmp_path / "project"
    assert main(["init", str(destination)]) == 0
    project = yaml.safe_load((destination / "project.yaml").read_text())
    candidate = destination / "results" / project["project"]["id"] / "all"
    reference = tmp_path / "reference"
    candidate.mkdir(parents=True)
    reference.mkdir()
    content = "gene_id\ts1\ng1\t10\n"
    (candidate / "counts.tsv").write_text(content)
    (reference / "counts.tsv").write_text(content)
    output = tmp_path / "verification.json"
    assert main([
        "verify", str(destination / "project.yaml"), "--reference", str(reference),
        "--output", str(output),
    ]) == 0
    assert output.is_file()
