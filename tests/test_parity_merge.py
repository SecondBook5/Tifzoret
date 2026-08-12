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

def _entry(id_, area, track="A"):
    return (
        f"- {{id: {id_}, track: {track}, area: {area}, verdict: OK, severity: S3, "
        "summary: s, detail: d, reference: [{file: r, line: 1}], "
        "frame: [{file: f, line: 1}], evidence: run, target: figures, "
        "status: open, commit: null}\n"
    )

def test_gsva_area_does_not_satisfy_sva_stage(tmp_path):
    # Regression: substring matching wrongly counted any 'gsva' finding as
    # covering the 'sva' stage. Token matching must not.
    (tmp_path / "A1.yaml").write_text(_entry("A1-1", "'Set2-F Hallmark GSVA (scoring method)'"))
    cov = coverage_report(str(tmp_path))
    assert "sva" not in cov["stages_seen"]
    assert "sva" in cov["stages_missing"]

def test_de_stage_not_satisfied_by_lookalike_tokens(tmp_path):
    # Regression: the 2-letter 'de' was matched by 'node', 'dendrograms',
    # 'config-defaults', 'seed-determinism'. None of these are the DE stage.
    for i, area in enumerate([
        "'Set1-B network node fill'",
        "'Set2-A PCA/correlation dendrograms'",
        "config-defaults",
        "pathways/seed-determinism",
    ]):
        (tmp_path / f"A{i}.yaml").write_text(_entry(f"A{i}-1", area))
    cov = coverage_report(str(tmp_path))
    assert "de" not in cov["stages_seen"]
    assert "de" in cov["stages_missing"]
    # But a genuine DE area still registers.
    (tmp_path / "Z.yaml").write_text(_entry("Z-1", "de/annotation"))
    assert "de" in coverage_report(str(tmp_path))["stages_seen"]

def test_plural_stem_keeps_real_stage_labels(tmp_path):
    # 'networks', 'pathways', 'regulators' must still cover their stages.
    (tmp_path / "A1.yaml").write_text(_entry("A1-1", "'networks (STRING score threshold)'"))
    (tmp_path / "A2.yaml").write_text(_entry("A2-1", "pathways/ora-universe"))
    (tmp_path / "A3.yaml").write_text(_entry("A3-1", "'regulators (VIPER aREA signature scaling)'"))
    seen = coverage_report(str(tmp_path))["stages_seen"]
    assert {"network", "pathways", "regulators"} <= seen

def test_missing_stage_detected_when_removed(tmp_path):
    # A single unrelated-area entry leaves 12 of 13 stages missing.
    (tmp_path / "A1.yaml").write_text(_entry("A1-1", "wgcna"))
    cov = coverage_report(str(tmp_path))
    assert cov["stages_seen"] == {"wgcna"}
    assert "de" in cov["stages_missing"] and "pathways" in cov["stages_missing"]

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
