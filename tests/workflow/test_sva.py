"""The SVA sensitivity stage must not assume how its input fit is parameterized.

``de.R`` used to persist a per-contrast DESeqDataSet releveled so the contrast
*denominator* was the reference (``utils.R`` ``resolve_contrast`` defaults
``reference_levels[[factor]]`` to the denominator). That guaranteed a coefficient
named literally ``<factor>_<numerator>_vs_<denominator>``, and ``sva.R`` resolved
the comparison by grepping for exactly that name.

The canonical fit is now the family fit, which is deliberately NOT releveled per
contrast -- a cell-means estimand is parameterization-free, so the family keeps
one reference for all of its estimands. Whenever a contrast's denominator is not
that reference, the old coefficient name does not exist and the grep yields
``NA``. The fixture below is that case: ``control`` is the family reference
(alphabetically first) while the contrast runs ``control`` vs ``treated``.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[2]
SVA_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "02_qc" / "sva.R"


# A condition effect on genes 1-60 and an independent hidden effect on genes
# 150-300, crossed with condition so sva has latent structure to find (n_sv > 0)
# without confounding the contrast.
FIXTURE_BUILDER_R = r"""
suppressPackageStartupMessages(library(DESeq2))
root <- commandArgs(trailingOnly = TRUE)[[1]]
dir.create(root, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

n_genes <- 500L
sample_ids <- sprintf("s%02d", 1:12)
condition <- rep(c("control", "treated"), each = 6L)
hidden <- rep(c("h1", "h2"), times = 6L)

samples <- data.frame(sample_id = sample_ids, condition = condition,
                      stringsAsFactors = FALSE)
write.table(samples, file.path(root, "samples.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

gene_ids <- sprintf("G%03d", 1:n_genes)
mu <- matrix(200, nrow = n_genes, ncol = length(sample_ids),
             dimnames = list(gene_ids, sample_ids))
mu[1:60, condition == "treated"] <- mu[1:60, condition == "treated"] * 4
mu[150:300, hidden == "h2"] <- mu[150:300, hidden == "h2"] * 5
counts <- matrix(stats::rnbinom(length(mu), mu = as.vector(mu), size = 50),
                 nrow = n_genes, dimnames = dimnames(mu))

rownames(samples) <- samples$sample_id
dds <- DESeq2::DESeqDataSetFromMatrix(counts, samples, ~ condition)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
# Guard the premise: the family reference is "control", so the coefficient the
# retired de.R path relied on ("condition_control_vs_treated") does NOT exist.
stopifnot(!("condition_control_vs_treated" %in% DESeq2::resultsNames(dds)))
saveRDS(dds, file.path(root, "dds.rds"))
saveRDS(DESeq2::varianceStabilizingTransformation(dds), file.path(root, "vst.rds"))

contrasts <- data.frame(contrast_id = "control_vs_treated", factor = "condition",
                        numerator = "control", denominator = "treated",
                        stringsAsFactors = FALSE)
write.table(contrasts, file.path(root, "contrasts.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
writeLines("DUMMY\tna\tG001\tG002", file.path(root, "sets.gmt"))
"""

_PROJECT_YAML = {
    "version": 2,
    "project": {"id": "sva_fixture", "title": "sva fixture", "description": "fixture"},
    "species": {"provider": "custom", "scientific_name": "synthetic", "taxonomy_id": None},
    "reference": {"genome_build": "synthetic", "annotation_release": None},
    "inputs": {"kind": "counts", "samples": "samples.tsv", "counts": "counts.tsv",
               "annotation": "annotation.tsv"},
    "analysis": {"design": "~ condition", "contrasts": "contrasts.tsv",
                 "profile": "full", "random_seed": 1},
    "resources": {"gene_sets": {"gmt": "sets.gmt"}},
    "figures": {"group": "condition"},
    "output": {"root": "out"},
}


def _build_fixture(tmp_path: Path) -> Path:
    fixture = tmp_path / "fixture"
    builder = tmp_path / "build_sva.R"
    builder.write_text(FIXTURE_BUILDER_R, encoding="utf-8")
    subprocess.run(["Rscript", "--vanilla", str(builder), str(fixture)],
                   check=True, capture_output=True, text=True)
    (fixture / "project.yaml").write_text(
        yaml.safe_dump(_PROJECT_YAML, sort_keys=False), encoding="utf-8")
    return fixture


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _rscript_has_sva() -> bool:
    if shutil.which("Rscript") is None:
        return False
    return subprocess.run(
        ["Rscript", "--vanilla", "-e",
         'quit(status = as.integer(!requireNamespace("sva", quietly = TRUE)))'],
        capture_output=True,
    ).returncode == 0


@pytest.mark.skipif(not _rscript_has_sva(), reason="sva not installed for this Rscript")
def test_sva_resolves_the_contrast_without_a_matching_coefficient_name(tmp_path):
    """sva.R runs against a family fit whose reference is not the denominator.

    Three things are asserted, and the second is what keeps the test honest: if
    ``sva`` finds no surrogate variables the sensitivity comparison is skipped
    entirely and the stage would pass without ever touching the code under test.
    """
    fixture = _build_fixture(tmp_path)
    outdir = tmp_path / "out"
    proc = subprocess.run(
        [
            "Rscript", "--vanilla", str(SVA_R),
            "--project-config", str(fixture / "project.yaml"),
            "--dds", str(fixture / "dds.rds"),
            "--vst", str(fixture / "vst.rds"),
            "--contrasts", str(fixture / "contrasts.tsv"),
            "--contrast-id", "control_vs_treated",
            "--outdir", str(outdir),
        ],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr

    # (1) The surrogate variables were estimated and written.
    surrogates = _read_tsv(outdir / "tables" / "surrogate_variables.tsv")
    assert len(surrogates) == 12, surrogates

    # (2) Non-vacuity gate: at least one SV, so the refit-and-compare branch ran.
    sv_columns = [name for name in surrogates[0] if name.startswith("SV")]
    assert sv_columns, f"no surrogate variables found -- gate is inert: {surrogates[0]}"

    # (3) Direction is right, not merely non-crashing. The fixture raises genes
    # 1-60 in `treated`, and the contrast is control vs treated, so those genes
    # must come back NEGATIVE. A fix that resolved the comparison backwards would
    # still run clean, so the sign is the assertion that matters.
    sensitivity = _read_tsv(outdir / "tables" / "sva_de_sensitivity.tsv")
    assert sensitivity, "sensitivity comparison table is empty"
    effects = {
        row["gene_id"]: float(row["log2_fold_change_original"])
        for row in sensitivity
        if row["log2_fold_change_original"] not in ("", "NA")
    }
    raised_in_treated = [effects[f"G{index:03d}"] for index in range(1, 61)
                         if f"G{index:03d}" in effects]
    assert len(raised_in_treated) >= 55, len(raised_in_treated)
    assert all(value < 0 for value in raised_in_treated), sorted(raised_in_treated)[-5:]
