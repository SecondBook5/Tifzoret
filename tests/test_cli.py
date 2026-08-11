from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from bulk_rna_frame.cli import main


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
        "results/synthetic_demo/inputs/input_manifest.json"
    )
    assert "--dry-run" not in observed["command"]
