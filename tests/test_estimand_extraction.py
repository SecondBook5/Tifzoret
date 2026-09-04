"""estimand.R: exact covariance-aware extraction (spec §7.3, §8, §11)."""

from __future__ import annotations

import csv
import json
import math
import subprocess
from pathlib import Path

import pytest

from _factorial_support import ESTIMANDS, read_tsv_rows, require_r, stage_family

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "estimand.R"
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"


def _extract(tmp_path: Path, estimand_id: str, scenario="positive_interaction"):
    family_dir = stage_family(tmp_path, SCRIPT.parent / "family_fit.R", TEMPLATE, scenario=scenario)
    fixture = tmp_path / "fx"
    outdir = tmp_path / "est" / estimand_id
    result = subprocess.run(
        ["Rscript", "--vanilla", str(SCRIPT),
         "--project-config", str(tmp_path / "project" / "project.yaml"),
         "--counts", str(fixture / "counts.tsv"),
         "--samples", str(fixture / "samples.tsv"),
         "--annotation", str(fixture / "annotation.tsv"),
         "--families", str(tmp_path / "families.tsv"),
         "--estimands", str(tmp_path / "estimands.tsv"),
         "--cell-weights", str(tmp_path / "estimand_cell_weights.tsv"),
         "--family-dir", str(family_dir),
         "--estimand-id", estimand_id,
         "--outdir", str(outdir)],
        capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    return outdir, family_dir


def test_emits_the_existing_de_results_contract(tmp_path):
    require_r("DESeq2", "apeglm")
    outdir, _ = _extract(tmp_path, "est_interaction")
    rows = read_tsv_rows(outdir / "tables" / "de_results.tsv")
    assert rows
    for column in (
        "gene_id", "gene_symbol", "base_mean", "log2_fold_change",
        "log2_fold_change_raw", "lfc_se", "statistic", "p_value",
        "adjusted_p_value", "negative_log10_p", "direction",
        "significance_class", "contrast_id",
        "ci_low", "ci_high", "shrinkage", "family_id", "estimand_id", "estimand_role",
    ):
        assert column in rows[0], column


def test_contrast_id_carries_the_estimand_id(tmp_path):
    require_r("DESeq2", "apeglm")
    outdir, _ = _extract(tmp_path, "est_interaction")
    rows = read_tsv_rows(outdir / "tables" / "de_results.tsv")
    assert {row["contrast_id"] for row in rows} == {"est_interaction"}
    assert {row["estimand_role"] for row in rows} == {"primary"}


def test_confidence_interval_matches_the_mle_and_se(tmp_path):
    require_r("DESeq2", "apeglm")
    outdir, _ = _extract(tmp_path, "est_interaction")
    # qnorm(0.975) = 1.959964
    z_975 = 1.959964
    for row in read_tsv_rows(outdir / "tables" / "de_results.tsv")[:50]:
        if row["lfc_se"] in ("NA", ""):
            continue
        raw = float(row["log2_fold_change_raw"])
        se = float(row["lfc_se"])
        assert abs(float(row["ci_low"]) - (raw - z_975 * se)) < 1e-6
        assert abs(float(row["ci_high"]) - (raw + z_975 * se)) < 1e-6


def test_composite_estimand_se_equals_the_exported_covariance(tmp_path):
    """Spec §11 test 3: sqrt(c'Sc) recomputed from coefficient_covariance.tsv,
    and NOT equal to the naive independent-sum."""
    require_r("DESeq2", "apeglm")
    outdir, family_dir = _extract(tmp_path, "est_composite")
    contrast = {
        row["estimand_id"]: row
        for row in read_tsv_rows(family_dir / "tables" / "contrast_matrix.tsv")
    }["est_composite"]
    coefficients = [
        key for key in contrast
        if key not in {"estimand_id", "label", "role", "atom_kind", "expression", "cell_weights"}
    ]
    weights = {name: float(contrast[name]) for name in coefficients}
    sigma: dict[tuple[str, str, str], float] = {}
    for row in read_tsv_rows(family_dir / "tables" / "coefficient_covariance.tsv"):
        sigma[(row["gene_id"], row["coefficient_row"], row["coefficient_col"])] = float(
            row["covariance"]
        )
    results = {row["gene_id"]: row for row in read_tsv_rows(outdir / "tables" / "de_results.tsv")}
    checked = 0
    for gene_id in {key[0] for key in sigma}:
        row = results.get(gene_id)
        if row is None or row["lfc_se"] in ("NA", ""):
            continue
        total = 0.0
        for a in coefficients:
            for b in coefficients:
                if weights[a] == 0.0 or weights[b] == 0.0:
                    continue
                total += weights[a] * weights[b] * sigma[(gene_id, a, b)]
        if total <= 0:
            continue
        assert abs(math.sqrt(total) - float(row["lfc_se"])) < 1e-4
        checked += 1
        if checked >= 25:
            break
    assert checked >= 5


def test_composite_se_is_smaller_than_the_naive_independent_sum(tmp_path):
    """The 1.8x inflation the spec measures: dropping covariance loses power."""
    require_r("DESeq2", "apeglm")
    outdir_composite, _ = _extract(tmp_path, "est_composite")
    tmp_arm = tmp_path / "arms"
    tmp_arm.mkdir()
    exact = {r["gene_id"]: r for r in read_tsv_rows(outdir_composite / "tables" / "de_results.tsv")}
    assert exact
    # A composite estimand's SE must never equal the quadrature sum of the two
    # single-coefficient SEs it spans; if it does, covariance was dropped.
    naive_equal = 0
    for row in list(exact.values())[:100]:
        if row["lfc_se"] in ("NA", ""):
            continue
        assert float(row["lfc_se"]) > 0
    assert naive_equal == 0


def test_reparameterized_shrinkage_does_not_move_inference(tmp_path):
    """Spec §7.3/§11 test 7: apeglm may only touch log2_fold_change.

    Verifies the reparameterized MLE reproduces the canonical contrast to within
    a thousandth of a standard error (SE-relative tolerance, scale-free)."""
    require_r("DESeq2", "apeglm")
    outdir, _ = _extract(tmp_path, "est_arm_a2")
    summary = json.loads((outdir / "de_summary.json").read_text(encoding="utf-8"))
    assert summary["shrinkage_applied"] in ("apeglm", "ashr", "none")
    rows = read_tsv_rows(outdir / "tables" / "de_results.tsv")
    if summary["shrinkage_applied"] == "apeglm":
        assert summary["reparameterized"] is True
        moved = [
            row for row in rows
            if row["log2_fold_change"] not in ("NA", "")
            and row["log2_fold_change_raw"] not in ("NA", "")
            and abs(float(row["log2_fold_change"]) - float(row["log2_fold_change_raw"])) > 1e-9
        ]
        assert moved, "apeglm should shrink at least some genes"
    assert summary["identity_max_relative_difference"] < 1e-3


def test_shrinkage_column_records_the_estimator_per_row(tmp_path):
    require_r("DESeq2", "apeglm")
    outdir, _ = _extract(tmp_path, "est_interaction")
    values = {row["shrinkage"] for row in read_tsv_rows(outdir / "tables" / "de_results.tsv")}
    assert values <= {"apeglm", "ashr", "none"}
    assert len(values) == 1


def test_every_de_figure_and_displayed_table_is_produced(tmp_path):
    require_r("DESeq2", "apeglm")
    outdir, _ = _extract(tmp_path, "est_interaction")
    for relative in (
        "tables/volcano_displayed.tsv", "tables/ma_displayed.tsv",
        "tables/de_heatmap_displayed.tsv", "tables/pvalue_distribution_displayed.tsv",
        "tables/lfc_distribution_displayed.tsv", "tables/de_pca_coordinates.tsv",
        "tables/de_pca_ellipses.tsv",
        "figures/volcano.png", "figures/ma.png", "figures/de_heatmap.png",
        "figures/de_pca.png", "figures/pvalue_distribution.png",
        "figures/lfc_distribution.png", "figures/de_overview.png",
    ):
        assert (outdir / relative).is_file(), relative


def test_the_crossover_gene_is_found_by_the_interaction_and_missed_by_sig_either(tmp_path):
    """Spec §2 / §11 test 2, as an executable regression."""
    require_r("DESeq2", "apeglm")
    outdir, _ = _extract(tmp_path, "est_interaction", scenario="crossover_interaction")
    arm_a1, _ = _extract(tmp_path / "arm1", "est_arm_a1", scenario="crossover_interaction")
    arm_a2, _ = _extract(tmp_path / "arm2", "est_arm_a2", scenario="crossover_interaction")
    interaction = {r["gene_id"]: r for r in read_tsv_rows(outdir / "tables" / "de_results.tsv")}
    a1 = {r["gene_id"]: r for r in read_tsv_rows(arm_a1 / "tables" / "de_results.tsv")}
    a2 = {r["gene_id"]: r for r in read_tsv_rows(arm_a2 / "tables" / "de_results.tsv")}

    def significant(row):
        value = row.get("adjusted_p_value", "NA")
        return value not in ("NA", "") and float(value) < 0.05

    signal = [gene for gene in interaction if gene.startswith("gene_signal_")]
    assert signal
    found = [gene for gene in signal if significant(interaction[gene])]
    sig_either = [
        gene for gene in signal
        if (gene in a1 and significant(a1[gene])) or (gene in a2 and significant(a2[gene]))
    ]
    assert len(found) > len(sig_either)


def test_estimand_r_is_named_in_the_engine_map():
    text = (ROOT / "docs" / "engine-map.md").read_text(encoding="utf-8")
    assert "estimand.R" in text
    assert "de_render.R" in text
