"""estimands.R maps cell weights onto model-matrix columns (spec §5)."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from _factorial_support import build_fixture, require_r

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "src" / "tifzoret" / "workflow" / "scripts"
ESTIMANDS_R = SCRIPTS / "estimands.R"


def _run_r(body: str, tmp_path: Path) -> dict:
    script = tmp_path / "probe.R"
    script.write_text(
        "suppressPackageStartupMessages({library(jsonlite)})\n"
        f'source("{ESTIMANDS_R.as_posix()}")\n' + body,
        encoding="utf-8",
    )
    result = subprocess.run(
        ["Rscript", "--vanilla", str(script)], capture_output=True, text=True, cwd=tmp_path
    )
    assert result.returncode == 0, result.stderr
    return json.loads((tmp_path / "out.json").read_text(encoding="utf-8"))


def _write_inputs(tmp_path: Path, weights: dict[str, float], cells="factor_a|factor_b",
                  design="~ factor_a * factor_b", atom_kind="cell") -> None:
    (tmp_path / "families.tsv").write_text(
        "family_id\tdesign\tcells\treference_levels\treplicate_unit\tshrinkage\tfilter\tdeclared\n"
        f"fam_a\t{design}\t{cells}\t\t\tapeglm\tdesign_aware\ttrue\n",
        encoding="utf-8",
    )
    (tmp_path / "estimands.tsv").write_text(
        "family_id\testimand_id\tlabel\trole\texpression\tatom_kind\n"
        f"fam_a\test_x\tLabel\tprimary\tsource-text\t{atom_kind}\n",
        encoding="utf-8",
    )
    lines = ["family_id\testimand_id\tcell_key\tweight"]
    for key, weight in weights.items():
        lines.append(f"fam_a\test_x\t{key}\t{weight!r}")
    (tmp_path / "estimand_cell_weights.tsv").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def test_reads_a_family_spec(tmp_path):
    require_r()
    _write_inputs(tmp_path, {"a2|b2": 1.0, "a1|b1": -1.0})
    out = _run_r(
        'spec <- read_family_spec("families.tsv", "fam_a")\n'
        'jsonlite::write_json(spec, "out.json", auto_unbox = TRUE)\n',
        tmp_path,
    )
    assert out["design"] == "~ factor_a * factor_b"
    assert out["cells"] == ["factor_a", "factor_b"]
    assert out["shrinkage"] == "apeglm"
    assert out["filter"] == "design_aware"


def test_reads_cell_weights_as_a_named_vector(tmp_path):
    require_r()
    _write_inputs(tmp_path, {"a2|b2": 1.0, "a2|b1": -1.0, "a1|b2": -1.0, "a1|b1": 1.0})
    out = _run_r(
        'w <- read_cell_weights("estimand_cell_weights.tsv", "est_x")\n'
        'jsonlite::write_json(as.list(w), "out.json", auto_unbox = TRUE)\n',
        tmp_path,
    )
    assert out == {"a2|b2": 1.0, "a2|b1": -1.0, "a1|b2": -1.0, "a1|b1": 1.0}


def test_compiles_the_interaction_contrast_vector(tmp_path):
    require_r("DESeq2")
    _write_inputs(tmp_path, {"a2|b2": 1.0, "a2|b1": -1.0, "a1|b2": -1.0, "a1|b1": 1.0})
    build_fixture("two_by_two", tmp_path / "fx")
    out = _run_r(
        """
        suppressPackageStartupMessages(library(DESeq2))
        md <- read.delim("fx/samples.tsv", stringsAsFactors = FALSE)
        rownames(md) <- md$sample_id
        md$factor_a <- factor(md$factor_a); md$factor_b <- factor(md$factor_b)
        cts <- as.matrix(read.delim("fx/counts.tsv", row.names = 1))
        dds <- DESeqDataSetFromMatrix(cts[, rownames(md)], md, ~ factor_a * factor_b)
        dds <- DESeq(dds, quiet = TRUE)
        spec <- read_family_spec("families.tsv", "fam_a")
        est  <- read_estimand_spec("estimands.tsv", "est_x")
        w    <- read_cell_weights("estimand_cell_weights.tsv", "est_x")
        cvec <- compile_contrast_vector(
          c(spec, est), w, md, stats::as.formula(spec$design), resultsNames(dds))
        jsonlite::write_json(as.list(cvec), "out.json", auto_unbox = TRUE)
        """,
        tmp_path,
    )
    # ~ factor_a * factor_b with refs a1/b1 -> the interaction IS the last coefficient.
    assert out["Intercept"] == 0.0
    assert out["factor_a_a2_vs_a1"] == 0.0
    assert out["factor_b_b2_vs_b1"] == 0.0
    assert out["factor_aa2.factor_bb2"] == 1.0


def test_compiles_a_simple_effect_contrast_vector(tmp_path):
    require_r("DESeq2")
    _write_inputs(tmp_path, {"a2|b1": 1.0, "a1|b1": -1.0})
    build_fixture("two_by_two", tmp_path / "fx")
    out = _run_r(
        """
        suppressPackageStartupMessages(library(DESeq2))
        md <- read.delim("fx/samples.tsv", stringsAsFactors = FALSE)
        rownames(md) <- md$sample_id
        md$factor_a <- factor(md$factor_a); md$factor_b <- factor(md$factor_b)
        cts <- as.matrix(read.delim("fx/counts.tsv", row.names = 1))
        dds <- DESeqDataSetFromMatrix(cts[, rownames(md)], md, ~ factor_a * factor_b)
        dds <- DESeq(dds, quiet = TRUE)
        spec <- read_family_spec("families.tsv", "fam_a")
        est  <- read_estimand_spec("estimands.tsv", "est_x")
        w    <- read_cell_weights("estimand_cell_weights.tsv", "est_x")
        cvec <- compile_contrast_vector(
          c(spec, est), w, md, stats::as.formula(spec$design), resultsNames(dds))
        jsonlite::write_json(as.list(cvec), "out.json", auto_unbox = TRUE)
        """,
        tmp_path,
    )
    assert out["factor_a_a2_vs_a1"] == 1.0
    assert out["factor_aa2.factor_bb2"] == 0.0


def test_nuisance_coefficients_receive_exactly_zero_weight(tmp_path):
    """Spec §5: sum-to-zero plus holding nuisance fixed must cancel it exactly,
    even when the covariate is unbalanced across cells."""
    require_r("DESeq2")
    _write_inputs(
        tmp_path,
        {"a2|b2": 1.0, "a2|b1": -1.0, "a1|b2": -1.0, "a1|b1": 1.0},
        design="~ factor_a * factor_b + nuisance_z",
    )
    build_fixture("unbalanced_nuisance", tmp_path / "fx")
    out = _run_r(
        """
        suppressPackageStartupMessages(library(DESeq2))
        md <- read.delim("fx/samples.tsv", stringsAsFactors = FALSE)
        rownames(md) <- md$sample_id
        for (v in c("factor_a", "factor_b", "nuisance_z")) md[[v]] <- factor(md[[v]])
        cts <- as.matrix(read.delim("fx/counts.tsv", row.names = 1))
        dds <- DESeqDataSetFromMatrix(
          cts[, rownames(md)], md, ~ factor_a * factor_b + nuisance_z)
        dds <- DESeq(dds, quiet = TRUE)
        spec <- read_family_spec("families.tsv", "fam_a")
        est  <- read_estimand_spec("estimands.tsv", "est_x")
        w    <- read_cell_weights("estimand_cell_weights.tsv", "est_x")
        cvec <- compile_contrast_vector(
          c(spec, est), w, md, stats::as.formula(spec$design), resultsNames(dds))
        jsonlite::write_json(as.list(cvec), "out.json", auto_unbox = TRUE)
        """,
        tmp_path,
    )
    nuisance = {name: value for name, value in out.items() if name.startswith("nuisance_z")}
    assert nuisance
    assert all(value == 0.0 for value in nuisance.values())


def test_a_coef_estimand_selects_one_coefficient(tmp_path):
    require_r("DESeq2")
    _write_inputs(tmp_path, {"factor_aa2.factor_bb2": 1.0}, atom_kind="coef")
    build_fixture("two_by_two", tmp_path / "fx")
    out = _run_r(
        """
        suppressPackageStartupMessages(library(DESeq2))
        md <- read.delim("fx/samples.tsv", stringsAsFactors = FALSE)
        rownames(md) <- md$sample_id
        md$factor_a <- factor(md$factor_a); md$factor_b <- factor(md$factor_b)
        cts <- as.matrix(read.delim("fx/counts.tsv", row.names = 1))
        dds <- DESeqDataSetFromMatrix(cts[, rownames(md)], md, ~ factor_a * factor_b)
        dds <- DESeq(dds, quiet = TRUE)
        spec <- read_family_spec("families.tsv", "fam_a")
        est  <- read_estimand_spec("estimands.tsv", "est_x")
        w    <- read_cell_weights("estimand_cell_weights.tsv", "est_x")
        cvec <- compile_contrast_vector(
          c(spec, est), w, md, stats::as.formula(spec$design), resultsNames(dds))
        jsonlite::write_json(as.list(cvec), "out.json", auto_unbox = TRUE)
        """,
        tmp_path,
    )
    assert out["factor_aa2.factor_bb2"] == 1.0
    assert sum(abs(value) for value in out.values()) == 1.0


def test_an_unknown_coefficient_name_is_a_hard_error(tmp_path):
    require_r("DESeq2")
    _write_inputs(tmp_path, {"not_a_coefficient": 1.0}, atom_kind="coef")
    script = tmp_path / "probe.R"
    script.write_text(
        f'source("{ESTIMANDS_R.as_posix()}")\n'
        'spec <- read_family_spec("families.tsv", "fam_a")\n'
        'est  <- read_estimand_spec("estimands.tsv", "est_x")\n'
        'w    <- read_cell_weights("estimand_cell_weights.tsv", "est_x")\n'
        'md <- data.frame(sample_id="s1", factor_a=factor("a1"), factor_b=factor("b1"))\n'
        'compile_contrast_vector(c(spec, est), w, md, ~ factor_a * factor_b, c("Intercept"))\n',
        encoding="utf-8",
    )
    result = subprocess.run(
        ["Rscript", "--vanilla", str(script)], capture_output=True, text=True, cwd=tmp_path
    )
    assert result.returncode != 0
    assert "not_a_coefficient" in result.stderr


def test_malformed_cell_key_arity_is_a_hard_error(tmp_path):
    """A cell key with wrong arity (e.g., level containing | that got split) must fail loudly."""
    require_r("DESeq2")
    # Manually write a malformed weights table with 3 parts instead of 2
    (tmp_path / "families.tsv").write_text(
        "family_id\tdesign\tcells\treference_levels\treplicate_unit\tshrinkage\tfilter\tdeclared\n"
        "fam_a\t~ factor_a * factor_b\tfactor_a|factor_b\t\t\tapeglm\tdesign_aware\ttrue\n",
        encoding="utf-8",
    )
    (tmp_path / "estimands.tsv").write_text(
        "family_id\testimand_id\tlabel\trole\texpression\tatom_kind\n"
        "fam_a\test_x\tLabel\tprimary\tsource-text\tcell\n",
        encoding="utf-8",
    )
    # Cell key with 3 parts (a2|extra|b2) instead of 2
    (tmp_path / "estimand_cell_weights.tsv").write_text(
        "family_id\testimand_id\tcell_key\tweight\n"
        "fam_a\test_x\ta2|extra|b2\t1.0\n"
        "fam_a\test_x\ta1|b1\t-1.0\n",
        encoding="utf-8",
    )
    script = tmp_path / "probe.R"
    script.write_text(
        f'source("{ESTIMANDS_R.as_posix()}")\n'
        'spec <- read_family_spec("families.tsv", "fam_a")\n'
        'est  <- read_estimand_spec("estimands.tsv", "est_x")\n'
        'w    <- read_cell_weights("estimand_cell_weights.tsv", "est_x")\n'
        'md <- data.frame(sample_id="s1", factor_a=factor("a1"), factor_b=factor("b1"))\n'
        'compile_contrast_vector(c(spec, est), w, md, ~ factor_a * factor_b, c("Intercept"))\n',
        encoding="utf-8",
    )
    result = subprocess.run(
        ["Rscript", "--vanilla", str(script)], capture_output=True, text=True, cwd=tmp_path
    )
    assert result.returncode != 0
    assert "arity" in result.stderr
    assert "a2|extra|b2" in result.stderr


def test_contrast_alignment_matches_deseq2_results_on_saturated_model(tmp_path):
    """Independent cross-check: contrast-derived LFC equals cell-mean LFC for a
    saturated model (~ factor_a * factor_b over 4 cells). If the contrast vector
    were mis-aligned to coefficients, this test would diverge wildly. Saturation
    means fitted cell means equal observed cell means, which licenses the comparison."""
    require_r("DESeq2")
    _write_inputs(tmp_path, {"a2|b2": 1.0, "a2|b1": -1.0, "a1|b2": -1.0, "a1|b1": 1.0})
    build_fixture("positive_interaction", tmp_path / "fx")
    script = tmp_path / "check.R"
    script.write_text(
        f'source("{ESTIMANDS_R.as_posix()}")\n'
        'suppressPackageStartupMessages(library(DESeq2))\n'
        'md <- read.delim("fx/samples.tsv", stringsAsFactors = FALSE)\n'
        'rownames(md) <- md$sample_id\n'
        'md$factor_a <- factor(md$factor_a); md$factor_b <- factor(md$factor_b)\n'
        'cts <- as.matrix(read.delim("fx/counts.tsv", row.names = 1))\n'
        'dds <- DESeqDataSetFromMatrix(cts[, rownames(md)], md, ~ factor_a * factor_b)\n'
        'dds <- DESeq(dds, quiet = TRUE)\n'
        'spec <- read_family_spec("families.tsv", "fam_a")\n'
        'est  <- read_estimand_spec("estimands.tsv", "est_x")\n'
        'w    <- read_cell_weights("estimand_cell_weights.tsv", "est_x")\n'
        'cvec <- compile_contrast_vector(c(spec, est), w, md, '
        '  stats::as.formula(spec$design), resultsNames(dds))\n'
        'res <- results(dds, contrast = cvec, independentFiltering = FALSE)\n'
        'signal_genes <- grep("^gene_signal_", rownames(res), value = TRUE)\n'
        'contrast_lfc <- median(res[signal_genes, "log2FoldChange"], na.rm = TRUE)\n'
        '# Independent from-counts calculation\n'
        'norm_cts <- counts(dds, normalized = TRUE)\n'
        'cell_means <- function(a_level, b_level) {\n'
        '  samples <- md$sample_id[md$factor_a == a_level & md$factor_b == b_level]\n'
        '  rowMeans(norm_cts[, samples, drop = FALSE])\n'
        '}\n'
        'mean_a2b2 <- cell_means("a2", "b2")\n'
        'mean_a2b1 <- cell_means("a2", "b1")\n'
        'mean_a1b2 <- cell_means("a1", "b2")\n'
        'mean_a1b1 <- cell_means("a1", "b1")\n'
        'observed_lfc <- (log2(mean_a2b2) - log2(mean_a2b1)) - '
        '                (log2(mean_a1b2) - log2(mean_a1b1))\n'
        'observed_median <- median(observed_lfc[signal_genes], na.rm = TRUE)\n'
        'cat(sprintf("contrast_lfc=%.6f observed_lfc=%.6f diff=%.6f\\n", '
        '  contrast_lfc, observed_median, abs(contrast_lfc - observed_median)))\n'
        'out <- list(\n'
        '  contrast_lfc = contrast_lfc,\n'
        '  observed_lfc = observed_median,\n'
        '  diff = abs(contrast_lfc - observed_median),\n'
        '  declared_truth = 2.0\n'
        ')\n'
        'jsonlite::write_json(out, "out.json", auto_unbox = TRUE)\n',
        encoding="utf-8",
    )
    result = subprocess.run(
        ["Rscript", "--vanilla", str(script)], capture_output=True, text=True, cwd=tmp_path
    )
    assert result.returncode == 0, result.stderr
    out = json.loads((tmp_path / "out.json").read_text(encoding="utf-8"))
    # The two computations must match tightly (saturated model)
    assert out["diff"] < 1e-2, \
        f"contrast vs observed LFC diverged: {out['contrast_lfc']} vs {out['observed_lfc']}"
    # And both should be close to the declared truth of 2.0 (looser tolerance for NB noise at n=4/cell)
    assert abs(out["contrast_lfc"] - out["declared_truth"]) < 0.35, \
        f"recovered LFC {out['contrast_lfc']} far from declared truth 2.0"
    assert abs(out["observed_lfc"] - out["declared_truth"]) < 0.35, \
        f"observed LFC {out['observed_lfc']} far from declared truth 2.0"


def test_estimands_r_is_named_in_the_engine_map():
    text = (ROOT / "docs" / "engine-map.md").read_text(encoding="utf-8")
    assert "estimands.R" in text
