"""term_test.R: full-vs-reduced LRTs on shared dispersions (spec §7.4)."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from _factorial_support import ESTIMANDS, read_tsv_rows, require_r, stage_family

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "term_test.R"
FAMILY_FIT_SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "family_fit.R"
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"

LATTICE = [
    ("factor_a_factor_b_interaction", "~ factor_a + factor_b", 1),
    ("factor_a_any", "~ factor_b", 2),
    ("factor_b_any", "~ factor_a", 2),
    ("any_effect", "~ 1", 3),
]


def _run(tmp_path: Path, term_test_id: str, scenario="positive_interaction", family_dir=None):
    if family_dir is None:
        family_dir = stage_family(tmp_path, FAMILY_FIT_SCRIPT, TEMPLATE, scenario=scenario)
    lines = ["family_id\tterm_test_id\treduced\tdf"]
    for name, reduced, df in LATTICE:
        lines.append(f"fam_a\t{name}\t{reduced}\t{df}")
    (tmp_path / "term_tests.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    outdir = tmp_path / "tt" / term_test_id
    result = subprocess.run(
        ["Rscript", "--vanilla", str(SCRIPT),
         "--project-config", str(tmp_path / "project" / "project.yaml"),
         "--annotation", str(tmp_path / "fx" / "annotation.tsv"),
         "--families", str(tmp_path / "families.tsv"),
         "--term-tests", str(tmp_path / "term_tests.tsv"),
         "--family-dir", str(family_dir),
         "--family-id", "fam_a",
         "--term-test-id", term_test_id,
         "--outdir", str(outdir)],
        capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    return outdir


def test_produces_the_lrt_contract(tmp_path):
    require_r("DESeq2")
    outdir = _run(tmp_path, "factor_a_factor_b_interaction")
    rows = read_tsv_rows(outdir / "tables" / "lrt_results.tsv")
    assert rows
    for column in (
        "gene_id", "gene_symbol", "base_mean", "lrt_statistic", "df",
        "p_value", "adjusted_p_value", "significant", "family_id", "term_test_id",
    ):
        assert column in rows[0], column
    assert (outdir / "figures" / "pvalue_distribution.png").is_file()


def test_records_the_reduced_model_and_df(tmp_path):
    require_r("DESeq2")
    summary = json.loads(
        (_run(tmp_path, "factor_a_factor_b_interaction") / "term_test_summary.json")
        .read_text(encoding="utf-8")
    )
    assert summary["reduced"] == "~ factor_a + factor_b"
    assert summary["df"] == 1
    assert summary["test"] == "LRT"


def test_uses_the_families_shared_gene_universe(tmp_path):
    require_r("DESeq2")
    family_dir = stage_family(tmp_path, FAMILY_FIT_SCRIPT, TEMPLATE)
    universe = {
        row["gene_id"] for row in read_tsv_rows(family_dir / "tables" / "tested_gene_universe.tsv")
        if row["retained"] == "TRUE"
    }
    outdir = _run(tmp_path, "any_effect", family_dir=family_dir)
    tested = {row["gene_id"] for row in read_tsv_rows(outdir / "tables" / "lrt_results.tsv")}
    assert tested <= universe


def test_the_interaction_lrt_detects_the_injected_interaction(tmp_path):
    require_r("DESeq2")
    outdir = _run(tmp_path, "factor_a_factor_b_interaction", scenario="positive_interaction")
    rows = {r["gene_id"]: r for r in read_tsv_rows(outdir / "tables" / "lrt_results.tsv")}
    signal = [g for g in rows if g.startswith("gene_signal_")]
    hits = [
        g for g in signal
        if rows[g]["adjusted_p_value"] not in ("NA", "")
        and float(rows[g]["adjusted_p_value"]) < 0.05
    ]
    assert len(hits) > 0.5 * len(signal)


def test_the_interaction_lrt_is_null_when_effects_are_additive(tmp_path):
    require_r("DESeq2")
    outdir = _run(
        tmp_path, "factor_a_factor_b_interaction", scenario="baseline_shift_no_interaction"
    )
    rows = {r["gene_id"]: r for r in read_tsv_rows(outdir / "tables" / "lrt_results.tsv")}
    signal = [g for g in rows if g.startswith("gene_signal_")]
    hits = [
        g for g in signal
        if rows[g]["adjusted_p_value"] not in ("NA", "")
        and float(rows[g]["adjusted_p_value"]) < 0.05
    ]
    assert len(hits) < 0.1 * len(signal)


def test_df_consistency_between_python_and_r(tmp_path):
    """Cross-language df consistency: declared df must match the model matrix rank difference."""
    require_r("DESeq2")
    # Stage a fixture to get metadata
    from _factorial_support import build_fixture
    fixture = build_fixture("positive_interaction", tmp_path / "fx")
    samples = read_tsv_rows(fixture / "samples.tsv")

    # Write metadata to a temp file for R to read
    metadata_path = tmp_path / "metadata_for_r.tsv"
    with metadata_path.open("w", encoding="utf-8") as f:
        f.write("sample_id\tfactor_a\tfactor_b\n")
        for row in samples:
            f.write(f"{row['sample_id']}\t{row['factor_a']}\t{row['factor_b']}\n")

    full_design = "~ factor_a * factor_b"

    # Compute the df for each test using R's model.matrix
    for term_test_id, reduced_design, expected_df in LATTICE:
        r_script = f"""
        metadata <- read.delim("{metadata_path}", stringsAsFactors = TRUE)
        full <- model.matrix({full_design}, data = metadata)
        reduced <- model.matrix({reduced_design}, data = metadata)
        cat(ncol(full) - ncol(reduced))
        """
        result = subprocess.run(
            ["Rscript", "-e", r_script],
            capture_output=True,
            text=True,
            cwd=tmp_path,
        )
        assert result.returncode == 0, f"R script failed for {term_test_id}: {result.stderr}"
        computed_df = int(result.stdout.strip())
        assert computed_df == expected_df, (
            f"df mismatch for {term_test_id}: "
            f"declared={expected_df}, computed={computed_df} "
            f"(full: {full_design}, reduced: {reduced_design})"
        )


def test_term_test_r_is_named_in_the_engine_map():
    assert "term_test.R" in (ROOT / "docs" / "engine-map.md").read_text(encoding="utf-8")
