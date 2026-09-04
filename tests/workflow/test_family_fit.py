"""family_fit.R: one Wald fit per family (spec §7.1, §6)."""

from __future__ import annotations

import json
from pathlib import Path

from _factorial_support import ESTIMANDS, read_tsv_rows, require_r, stage_family

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "family_fit.R"
)
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"


def test_produces_every_declared_output(tmp_path):
    require_r("DESeq2")
    outdir = stage_family(tmp_path, SCRIPT, TEMPLATE)
    for relative in (
        "objects/deseq2.rds",
        "tables/coefficient_covariance.tsv",
        "tables/contrast_matrix.tsv",
        "tables/tested_gene_universe.tsv",
        "tables/design_diagnostics.tsv",
        "figures/design_diagnostics.png",
        "figures/design_diagnostics.pdf",
        "family_summary.json",
    ):
        assert (outdir / relative).is_file(), relative


def test_contrast_matrix_has_one_row_per_estimand(tmp_path):
    require_r("DESeq2")
    rows = read_tsv_rows(stage_family(tmp_path, SCRIPT, TEMPLATE) / "tables" / "contrast_matrix.tsv")
    assert {row["estimand_id"] for row in rows} == set(ESTIMANDS)


def test_contrast_matrix_records_the_interaction_vector(tmp_path):
    require_r("DESeq2")
    rows = {
        r["estimand_id"]: r
        for r in read_tsv_rows(
            stage_family(tmp_path, SCRIPT, TEMPLATE) / "tables" / "contrast_matrix.tsv"
        )
    }
    interaction = rows["est_interaction"]
    assert float(interaction["factor_aa2.factor_bb2"]) == 1.0
    assert float(interaction["Intercept"]) == 0.0
    assert float(interaction["factor_a_a2_vs_a1"]) == 0.0


def test_gene_universe_records_the_filter_used(tmp_path):
    require_r("DESeq2", "edgeR")
    rows = read_tsv_rows(
        stage_family(tmp_path, SCRIPT, TEMPLATE) / "tables" / "tested_gene_universe.tsv"
    )
    assert rows
    assert {row["filter"] for row in rows} == {"design_aware"}
    assert any(row["retained"] == "TRUE" for row in rows)


def test_total_count_filter_is_selectable(tmp_path):
    require_r("DESeq2")
    rows = read_tsv_rows(
        stage_family(tmp_path, SCRIPT, TEMPLATE, gene_filter="total_count")
        / "tables"
        / "tested_gene_universe.tsv"
    )
    assert {row["filter"] for row in rows} == {"total_count"}


def test_summary_reports_rank_and_condition_number(tmp_path):
    require_r("DESeq2")
    summary = json.loads(
        (stage_family(tmp_path, SCRIPT, TEMPLATE) / "family_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["family_id"] == "fam_a"
    assert summary["rank"] == 4
    assert summary["condition_number"] > 0
    assert summary["residual_df"] > 0
    assert summary["estimands"] == list(ESTIMANDS)


def test_diagnostics_include_per_sample_leverage(tmp_path):
    require_r("DESeq2")
    rows = read_tsv_rows(
        stage_family(tmp_path, SCRIPT, TEMPLATE) / "tables" / "design_diagnostics.tsv"
    )
    leverage = [row for row in rows if row["metric"] == "leverage"]
    assert leverage
    assert all(0.0 <= float(row["value"]) <= 1.0 + 1e-9 for row in leverage)


def test_cell_sizes_are_reported(tmp_path):
    require_r("DESeq2")
    rows = read_tsv_rows(
        stage_family(tmp_path, SCRIPT, TEMPLATE) / "tables" / "design_diagnostics.tsv"
    )
    cells = [row for row in rows if row["metric"] == "cell_size"]
    assert len(cells) == 4


def test_covariance_export_is_symmetric(tmp_path):
    require_r("DESeq2")
    rows = read_tsv_rows(
        stage_family(tmp_path, SCRIPT, TEMPLATE) / "tables" / "coefficient_covariance.tsv"
    )
    assert rows
    lookup = {
        (row["gene_id"], row["coefficient_row"], row["coefficient_col"]): float(
            row["covariance"]
        )
        for row in rows
    }
    gene = rows[0]["gene_id"]
    pairs = [(r, c) for (g, r, c) in lookup if g == gene]
    for row_name, col_name in pairs[:20]:
        mirrored = lookup.get((gene, col_name, row_name))
        assert mirrored is not None
        assert abs(mirrored - lookup[(gene, row_name, col_name)]) < 1e-9


def test_covariance_matches_deseq2_lfcse(tmp_path):
    """
    Validate that sqrt(c'Σc) from our exported coefficient_covariance.tsv matches
    DESeq2's lfcSE for single-coefficient estimands.
    """
    require_r("DESeq2")
    import subprocess

    outdir = stage_family(tmp_path, SCRIPT, TEMPLATE)

    # Load coefficient covariance
    cov_rows = read_tsv_rows(outdir / "tables" / "coefficient_covariance.tsv")
    # Load contrast matrix
    contrast_rows = {
        r["estimand_id"]: r
        for r in read_tsv_rows(outdir / "tables" / "contrast_matrix.tsv")
    }

    # Find a single-coefficient estimand (one with exactly one ±1 in the contrast vector)
    # est_arm_a1 should be a simple contrast
    estimand_id = "est_arm_a1"
    contrast = contrast_rows[estimand_id]

    # Extract contrast vector (all coefficient columns)
    coefficient_names = [
        k for k in contrast.keys() if k not in {"estimand_id", "label", "role", "atom_kind", "expression", "cell_weights"}
    ]
    contrast_vector = {name: float(contrast[name]) for name in coefficient_names}

    # Verify it's a single-coefficient contrast (only one nonzero entry with value ±1)
    nonzero = {k: v for k, v in contrast_vector.items() if abs(v) > 1e-9}
    assert len(nonzero) == 1, f"Expected single-coefficient contrast, got {nonzero}"
    coef_name = list(nonzero.keys())[0]

    # Write an R script to extract DESeq2's lfcSE for this coefficient and a few genes
    r_script = tmp_path / "validate_cov.R"
    r_script.write_text(
        f"""
suppressPackageStartupMessages(library(DESeq2))
dds <- readRDS("{outdir / 'objects' / 'deseq2.rds'}")
res <- results(dds, name = "{coef_name}")
# Pick first 5 genes that have finite lfcSE
valid <- which(is.finite(res$lfcSE))
genes <- rownames(res)[valid[1:min(5, length(valid))]]
for (gene_id in genes) {{
  cat(gene_id, res[gene_id, "lfcSE"], "\\n", sep="\\t")
}}
""",
        encoding="utf-8",
    )

    result = subprocess.run(
        ["Rscript", "--vanilla", str(r_script)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr

    # Parse DESeq2's lfcSE values
    deseq2_se = {}
    for line in result.stdout.strip().split("\n"):
        if line:
            parts = line.split("\t")
            deseq2_se[parts[0]] = float(parts[1])

    # Build covariance matrix per gene and compute sqrt(c'Σc)
    cov_lookup = {
        (row["gene_id"], row["coefficient_row"], row["coefficient_col"]): float(
            row["covariance"]
        )
        for row in cov_rows
    }

    discrepancies = []
    for gene_id, expected_se in deseq2_se.items():
        # Extract Σ for this gene
        variance = cov_lookup.get((gene_id, coef_name, coef_name))
        if variance is None:
            continue  # Gene not in top N exported
        computed_se = variance**0.5
        rel_error = abs(computed_se - expected_se) / (expected_se + 1e-12)
        discrepancies.append(
            {
                "gene_id": gene_id,
                "computed_se": computed_se,
                "deseq2_se": expected_se,
                "rel_error": rel_error,
            }
        )

    # Assert tight agreement (tolerance 1e-4 relative)
    for d in discrepancies:
        assert (
            d["rel_error"] < 1e-4
        ), f"Gene {d['gene_id']}: computed SE {d['computed_se']:.6f} vs DESeq2 {d['deseq2_se']:.6f} (rel error {d['rel_error']:.6e})"


def test_an_unbalanced_nuisance_design_still_fits(tmp_path):
    require_r("DESeq2")
    outdir = stage_family(
        tmp_path,
        SCRIPT,
        TEMPLATE,
        scenario="unbalanced_nuisance",
        design="~ factor_a * factor_b + nuisance_z",
    )
    rows = {
        r["estimand_id"]: r
        for r in read_tsv_rows(outdir / "tables" / "contrast_matrix.tsv")
    }
    nuisance = {
        key: float(value)
        for key, value in rows["est_interaction"].items()
        if key.startswith("nuisance_z")
    }
    assert nuisance
    assert all(value == 0.0 for value in nuisance.values())


def test_family_fit_is_named_in_the_engine_map():
    assert "family_fit.R" in (ROOT / "docs" / "engine-map.md").read_text(encoding="utf-8")
