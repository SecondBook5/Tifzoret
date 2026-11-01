"""The raw featureCounts result is cached, so editing unrelated config never recounts the BAMs."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "01_inputs" / "materialize_inputs.py"


@pytest.fixture()
def mat(monkeypatch):
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("materialize_cache_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "command_version", lambda command: "featureCounts v-test")
    monkeypatch.delenv("TIFZORET_RECOUNT", raising=False)
    return module


def make_project(tmp_path: Path, **counting) -> SimpleNamespace:
    bam, gtf = tmp_path / "a.bam", tmp_path / "genes.gtf"
    bam.write_bytes(b"bam")
    gtf.write_text("gtf")
    settings = {"strandedness": "reverse", "threads": 8, "paired_end": True, **counting}
    return SimpleNamespace(
        bam_paths=(bam,),
        gtf=gtf,
        sample_rows=[{"sample_id": "a"}],
        source_kind="bam",
        config={"counting": settings, "analysis": {"design": "~ condition"}},
    )


def test_fingerprint_ignores_edits_outside_the_counting_inputs(tmp_path, mat):
    project = make_project(tmp_path)
    before = mat.counting_fingerprint(project)
    project.config["analysis"]["design"] = "~ batch + condition"
    project.config["figures"] = {"group": "condition"}
    assert mat.counting_fingerprint(project) == before


def test_fingerprint_changes_with_counting_settings_bams_gtf_and_tool_version(tmp_path, mat, monkeypatch):
    project = make_project(tmp_path)
    before = mat.counting_fingerprint(project)
    project.config["counting"]["strandedness"] = "forward"
    assert mat.counting_fingerprint(project) != before
    project.config["counting"]["strandedness"] = "reverse"
    assert mat.counting_fingerprint(project) == before
    project.bam_paths[0].write_bytes(b"a longer bam")
    assert mat.counting_fingerprint(project) != before
    project.bam_paths[0].write_bytes(b"bam")
    project.gtf.write_text("a different gtf")
    assert mat.counting_fingerprint(project) != before
    project.gtf.write_text("gtf")
    monkeypatch.setattr(mat, "command_version", lambda command: "featureCounts v-other")
    assert mat.counting_fingerprint(project) != before


def test_cache_round_trip_restores_the_strand_decision(tmp_path, mat):
    cache = tmp_path / "results" / ".cache" / "counting" / "abc"
    raw = tmp_path / "raw.tsv"
    raw.write_text("Geneid\ta\ng1\t5\n")
    mat.save_counting_cache(cache, 2, {1: 3.9, 2: 58.3}, raw, "log text", 0.937)
    work = tmp_path / "work"
    work.mkdir()
    mode, rates, restored, log, dominance = mat.load_counting_cache(cache, work)
    assert (mode, rates, log, dominance) == (2, {1: 3.9, 2: 58.3}, "log text", 0.937)
    assert restored.read_text() == raw.read_text()
    assert not cache.with_name("abc.partial").exists()


def test_a_missing_cache_or_a_forced_recount_is_a_miss(tmp_path, mat, monkeypatch):
    work = tmp_path / "work"
    work.mkdir()
    cache = tmp_path / "c"
    assert mat.load_counting_cache(cache, work) is None
    raw = tmp_path / "raw.tsv"
    raw.write_text("x")
    mat.save_counting_cache(cache, 0, {}, raw, "log", None)
    assert mat.load_counting_cache(cache, work) is not None
    monkeypatch.setenv("TIFZORET_RECOUNT", "1")
    assert mat.load_counting_cache(cache, work) is None


def test_the_cache_lives_under_the_runs_cache_folder(tmp_path, mat):
    project = make_project(tmp_path)
    counts = tmp_path / "results" / "inputs" / "counts.tsv"
    cache = mat.counting_cache_dir(counts, project)
    assert cache.parent == tmp_path / "results" / ".cache" / "counting"
    assert len(cache.name) == 24


def test_a_second_run_does_not_count_again(tmp_path, mat, monkeypatch):
    project = make_project(tmp_path)
    counts, annotation = (
        tmp_path / "results" / "inputs" / "counts.tsv",
        tmp_path / "results" / "inputs" / "annotation.tsv",
    )
    calls = {"count": 0}

    def fake_featurecounts(proj, output, mode):
        calls["count"] += 1
        Path(output).write_text("Geneid\ta\ng1\t5\n")
        return "Successfully assigned fragments : 5 (50.0%)"

    monkeypatch.setattr(mat, "run_featurecounts", fake_featurecounts)
    monkeypatch.setattr(mat, "parse_featurecounts", lambda raw, proj, out: (["g1"], {}, None))
    monkeypatch.setattr(mat, "write_annotation", lambda *args, **kwargs: None)
    monkeypatch.setattr(mat, "bam_record", lambda sample_id, bam: {"sample_id": sample_id})
    monkeypatch.setattr(mat, "_symbol_resolution_settings", lambda proj: None)
    first = mat.materialize_bams(project, counts, annotation)
    second = mat.materialize_bams(project, counts, annotation)
    assert calls["count"] == 1
    assert first["counting"]["resolved_strand_mode"] == second["counting"]["resolved_strand_mode"] == 2
    assert second["counting"]["final_assignment_log"] == first["counting"]["final_assignment_log"]
    monkeypatch.setenv("TIFZORET_RECOUNT", "1")
    mat.materialize_bams(project, counts, annotation)
    assert calls["count"] == 2
