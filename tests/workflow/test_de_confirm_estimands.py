"""edgeR confirmation accepts arbitrary estimands (spec §10)."""

from __future__ import annotations

import subprocess
from pathlib import Path

from _factorial_support import read_tsv_rows, require_r

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "de_confirm.R"


def test_pairwise_only_stop_is_gone():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "de_confirm supports pairwise contrasts only" not in text


def test_confirms_an_interaction_estimand(tmp_path):
    require_r("DESeq2", "edgeR", "apeglm")
    from test_estimand_extraction import _extract

    de_dir, family_dir = _extract(tmp_path, "est_interaction")
    outdir = tmp_path / "confirm"
    result = subprocess.run(
        ["Rscript", "--vanilla", str(SCRIPT),
         "--project-config", str(tmp_path / "project" / "project.yaml"),
         "--counts", str(tmp_path / "fx" / "counts.tsv"),
         "--samples", str(tmp_path / "fx" / "samples.tsv"),
         "--annotation", str(tmp_path / "fx" / "annotation.tsv"),
         "--families", str(tmp_path / "families.tsv"),
         "--estimands", str(tmp_path / "estimands.tsv"),
         "--cell-weights", str(tmp_path / "estimand_cell_weights.tsv"),
         "--family-dir", str(family_dir),
         "--estimand-id", "est_interaction",
         "--de", str(de_dir / "tables" / "de_results.tsv"),
         "--outdir", str(outdir)],
        capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    rows = read_tsv_rows(outdir / "tables" / "edger_results.tsv")
    assert rows
    assert (outdir / "figures" / "de_concordance.png").is_file()


def test_both_engines_agree_in_sign_on_the_injected_interaction(tmp_path):
    """Spec §11 test 5: one shared contrast-vector function, so signs must agree."""
    require_r("DESeq2", "edgeR", "apeglm")
    from test_estimand_extraction import _extract

    de_dir, family_dir = _extract(tmp_path, "est_interaction")
    outdir = tmp_path / "confirm"
    subprocess.run(
        ["Rscript", "--vanilla", str(SCRIPT),
         "--project-config", str(tmp_path / "project" / "project.yaml"),
         "--counts", str(tmp_path / "fx" / "counts.tsv"),
         "--samples", str(tmp_path / "fx" / "samples.tsv"),
         "--annotation", str(tmp_path / "fx" / "annotation.tsv"),
         "--families", str(tmp_path / "families.tsv"),
         "--estimands", str(tmp_path / "estimands.tsv"),
         "--cell-weights", str(tmp_path / "estimand_cell_weights.tsv"),
         "--family-dir", str(family_dir),
         "--estimand-id", "est_interaction",
         "--de", str(de_dir / "tables" / "de_results.tsv"),
         "--outdir", str(outdir)],
        check=True, capture_output=True, text=True, cwd=tmp_path)
    rows = read_tsv_rows(outdir / "tables" / "de_concordance_displayed.tsv")
    signal = [
        row for row in rows
        if row["gene_id"].startswith("gene_signal_")
        and row["deseq2_log2_fold_change"] not in ("NA", "")
        and row["edger_log2_fold_change"] not in ("NA", "")
    ]
    assert signal
    agreeing = [
        row for row in signal
        if (float(row["deseq2_log2_fold_change"]) > 0) == (float(row["edger_log2_fold_change"]) > 0)
    ]
    assert len(agreeing) > 0.9 * len(signal)


def test_concordance_uses_the_shared_gene_universe(tmp_path):
    require_r("DESeq2", "edgeR", "apeglm")
    from test_estimand_extraction import _extract

    de_dir, family_dir = _extract(tmp_path, "est_interaction")
    universe = {
        row["gene_id"]
        for row in read_tsv_rows(family_dir / "tables" / "tested_gene_universe.tsv")
        if row["retained"] == "TRUE"
    }
    outdir = tmp_path / "confirm"
    subprocess.run(
        ["Rscript", "--vanilla", str(SCRIPT),
         "--project-config", str(tmp_path / "project" / "project.yaml"),
         "--counts", str(tmp_path / "fx" / "counts.tsv"),
         "--samples", str(tmp_path / "fx" / "samples.tsv"),
         "--annotation", str(tmp_path / "fx" / "annotation.tsv"),
         "--families", str(tmp_path / "families.tsv"),
         "--estimands", str(tmp_path / "estimands.tsv"),
         "--cell-weights", str(tmp_path / "estimand_cell_weights.tsv"),
         "--family-dir", str(family_dir),
         "--estimand-id", "est_interaction",
         "--de", str(de_dir / "tables" / "de_results.tsv"),
         "--outdir", str(outdir)],
        check=True, capture_output=True, text=True, cwd=tmp_path)
    tested = {row["gene_id"] for row in read_tsv_rows(outdir / "tables" / "edger_results.tsv")}
    assert tested <= universe
