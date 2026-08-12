from tools.parity.merge_fragments import merge_fragments, coverage_report

def test_merge_sorts_by_severity(tmp_path):
    (tmp_path / "A1.yaml").write_text(
        "- {id: A1-1, track: A, area: de, verdict: OK, severity: S3, "
        "summary: s, detail: d, reference: [{file: r, line: 1}], "
        "frame: [{file: f, line: 1}], evidence: run, target: figures, "
        "status: open, commit: null}\n")
    (tmp_path / "A2.yaml").write_text(
        "- {id: A2-1, track: A, area: pathways, verdict: REGRESSION, severity: S1, "
        "summary: s, detail: d, reference: [{file: r, line: 1}], "
        "frame: [{file: f, line: 1}], evidence: run, target: stats, "
        "status: open, commit: null}\n")
    merged = merge_fragments(str(tmp_path))
    assert [m["id"] for m in merged] == ["A2-1", "A1-1"]  # S1 before S3

def test_coverage_report_flags_missing(tmp_path):
    (tmp_path / "A1.yaml").write_text(
        "- {id: A1-1, track: A, area: de, verdict: OK, severity: S3, "
        "summary: s, detail: d, reference: [{file: r, line: 1}], "
        "frame: [{file: f, line: 1}], evidence: run, target: figures, "
        "status: open, commit: null}\n")
    cov = coverage_report(str(tmp_path))
    assert "de" in cov["stages_seen"]
    assert cov["panels_missing"]  # 16 panels not all covered

def test_merge_detects_duplicate_ids(tmp_path):
    (tmp_path / "A1.yaml").write_text(
        "- {id: dup-1, track: A, area: qc, verdict: OK, severity: S3, "
        "summary: s, detail: d, reference: [{file: r, line: 1}], "
        "frame: [{file: f, line: 1}], evidence: run, target: figures, "
        "status: open, commit: null}\n")
    (tmp_path / "A2.yaml").write_text(
        "- {id: dup-1, track: A, area: de, verdict: OK, severity: S3, "
        "summary: s, detail: d, reference: [{file: r, line: 1}], "
        "frame: [{file: f, line: 1}], evidence: run, target: figures, "
        "status: open, commit: null}\n")
    try:
        merge_fragments(str(tmp_path))
        assert False, "expected ValueError for duplicate id"
    except ValueError as e:
        assert "duplicate" in str(e).lower()

def test_merge_sorts_by_track_within_severity(tmp_path):
    (tmp_path / "A1.yaml").write_text(
        "- {id: A1-1, track: A, area: qc, verdict: OK, severity: S2, "
        "summary: s, detail: d, reference: [{file: r, line: 1}], "
        "frame: [{file: f, line: 1}], evidence: run, target: figures, "
        "status: open, commit: null}\n")
    (tmp_path / "B1.yaml").write_text(
        "- {id: B1-1, track: B, area: Set1-A, verdict: OK, severity: S2, "
        "summary: s, detail: d, reference: [{file: r, line: 1}], "
        "frame: [{file: f, line: 1}], evidence: run, target: figures, "
        "status: open, commit: null}\n")
    (tmp_path / "C1.yaml").write_text(
        "- {id: C1-1, track: C, area: config, verdict: OK, severity: S2, "
        "summary: s, detail: d, reference: [{file: r, line: 1}], "
        "frame: [{file: f, line: 1}], evidence: run, target: figures, "
        "status: open, commit: null}\n")
    merged = merge_fragments(str(tmp_path))
    # Same severity S2, so sort by track A < B < C
    assert [m["id"] for m in merged] == ["A1-1", "B1-1", "C1-1"]
