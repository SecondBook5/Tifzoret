"""``utils.R``'s shared DESeq2 fit must degrade instead of aborting the workflow.

``DESeq()``'s parametric dispersion-trend fit does not warn when it cannot be
fitted -- it raises:

    all gene-wise dispersion estimates are within 2 orders of magnitude
    from the minimum value, and so the standard curve fitting techniques
    will not work.

DESeq2 raises this only when every gene-wise estimate falls below
``100 * minDisp`` (1e-6), i.e. when there is essentially no overdispersion
anywhere. Random counts never reach that floor -- even pure Poisson sampling
leaves the maximum gene-wise estimate around 9e-3 -- but *deterministic* counts
do. That is why the shipped synthetic templates hit this in a real run while
every randomly-generated fixture in the suite sailed past it.

``family_fit.R`` had guarded this since the beginning. ``sva.R``, which refits the
design augmented with surrogate variables, called ``DESeq()`` bare, so an
exploratory sensitivity module could take the entire run down on data the main DE
path handled fine. Both now share one helper, and this test pins its contract
directly: the fallback engages when it must, stays out of the way when it must
not, and never masks an unrelated failure.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
UTILS_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "utils.R"


# Three cases through the one helper, reported as JSON for the assertions below.
#
#   deterministic -> exact counts, zero overdispersion. The bare call must fail
#                    (asserted here, so the fixture cannot silently stop being the
#                    hard case) and the helper must succeed via "gene-wise".
#   dispersed     -> ordinary negative-binomial counts. The helper must NOT
#                    degrade; a fallback that always fires would silently drop the
#                    dispersion trend for every real study.
#   unrelated     -> a design whose model matrix is rank-deficient. The helper
#                    must still raise, because it catches one specific DESeq2
#                    message and must never swallow anything else.
PROBE_R = r"""
suppressPackageStartupMessages({library(DESeq2)})
args <- commandArgs(trailingOnly = TRUE)
source(args[[1]], local = FALSE)

sample_ids <- sprintf("s%02d", 1:6)
condition <- rep(c("control", "treated"), each = 3L)
n_genes <- 200L
gene_ids <- sprintf("G%03d", 1:n_genes)
mu <- matrix(400, nrow = n_genes, ncol = 6L, dimnames = list(gene_ids, sample_ids))
mu[1:40, condition == "treated"] <- mu[1:40, condition == "treated"] * 3

samples <- data.frame(sample_id = sample_ids, condition = condition,
                      stringsAsFactors = FALSE)
rownames(samples) <- samples$sample_id

build <- function(counts, design) {
  DESeq2::DESeqDataSetFromMatrix(counts, samples, design)
}

set.seed(11)
deterministic <- matrix(as.integer(round(as.vector(mu))), nrow = n_genes,
                        dimnames = dimnames(mu))
dispersed <- matrix(stats::rnbinom(length(mu), mu = as.vector(mu), size = 12),
                    nrow = n_genes, dimnames = dimnames(mu))

result <- list()

# --- deterministic: the bare call must fail, the helper must recover ----------
bare <- try(DESeq2::DESeq(build(deterministic, ~ condition),
                          fitType = "parametric", quiet = TRUE), silent = TRUE)
result$bare_call_fails <- inherits(bare, "try-error")
result$bare_message <- if (inherits(bare, "try-error")) conditionMessage(attr(bare, "condition")) else ""
helped <- fit_deseq_with_dispersion_fallback(build(deterministic, ~ condition))
result$deterministic_fit <- helped$dispersion_fit
result$deterministic_has_results <- nrow(DESeq2::results(helped$dds)) == n_genes

# --- dispersed: the helper must NOT degrade ----------------------------------
ordinary <- fit_deseq_with_dispersion_fallback(build(dispersed, ~ condition))
result$dispersed_fit <- ordinary$dispersion_fit

# --- unrelated failure: must still propagate ---------------------------------
samples$constant <- "same_for_every_sample"
rank_deficient <- try(
  fit_deseq_with_dispersion_fallback(
    DESeq2::DESeqDataSetFromMatrix(dispersed, samples, ~ condition + constant)),
  silent = TRUE)
result$unrelated_still_raises <- inherits(rank_deficient, "try-error")

cat(jsonlite::toJSON(result, auto_unbox = TRUE))
"""


def _rscript_has_deseq2() -> bool:
    if shutil.which("Rscript") is None:
        return False
    return subprocess.run(
        ["Rscript", "--vanilla", "-e",
         'quit(status = as.integer(!requireNamespace("DESeq2", quietly = TRUE)))'],
        capture_output=True,
    ).returncode == 0


@pytest.mark.skipif(not _rscript_has_deseq2(), reason="DESeq2 not available for this Rscript")
def test_dispersion_fallback_engages_only_when_the_trend_cannot_be_fitted(tmp_path):
    probe = tmp_path / "probe.R"
    probe.write_text(PROBE_R, encoding="utf-8")
    proc = subprocess.run(
        ["Rscript", "--vanilla", str(probe), str(UTILS_R)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    report = json.loads(proc.stdout[proc.stdout.index("{"):])

    # Premise: deterministic counts really are the hard case. Without this the
    # "gene-wise" assertion below could pass on data that never needed a fallback.
    assert report["bare_call_fails"], (
        "deterministic counts no longer break the parametric fit -- this fixture "
        "has stopped exercising the failure it exists to cover"
    )
    assert "all gene-wise dispersion estimates are within" in report["bare_message"], (
        report["bare_message"])

    # The helper recovers, reports which path ran, and returns a usable fit.
    assert report["deterministic_fit"] == "gene-wise", report
    assert report["deterministic_has_results"], report

    # ...and does not degrade ordinary data. A fallback that always fired would
    # quietly discard the dispersion trend for every real study, which is a
    # statistical regression no output would advertise.
    assert report["dispersed_fit"] == "parametric", report

    # Only the one DESeq2 message is caught; anything else must still surface.
    assert report["unrelated_still_raises"], (
        "an unrelated fit failure was swallowed by the dispersion fallback")
