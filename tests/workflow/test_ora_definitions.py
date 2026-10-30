"""Over-representation must test genes that can be annotated, against effects that were measured.

Two defects in ``pathways.R``:

* the background was every tested symbol -- including the ~700 genes that carry only an
  Ensembl accession and, more generally, the third of the transcriptome that no GO term
  contains -- so a hypergeometric test counted genes that cannot be in any set as
  non-members, inflating every enrichment (FDR 1e-19 where the standard annotated
  background gives 1e-16). The background is now, per collection, the tested genes that
  belong to at least one set of that collection, which is clusterProfiler's definition.
* the foreground took its up/down classes from the apeglm-*shrunken* fold change, so
  whether a gene entered the enrichment depended on whether an optimiser converged
  (Npy, 0/0/0 reads against 16/262/793, was displayed at -0.15 and left out). The
  foreground now uses the unshrunk estimate, which is the quantity the test inferred.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
UTILS_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "utils.R"


def _r(code: str):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    wrapped = (
        "suppressPackageStartupMessages({source(commandArgs(TRUE)[1]); library(dplyr)});"
        + code
    )
    result = subprocess.run(["Rscript", "--vanilla", "-e", wrapped, str(UTILS_R)], capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(result.stderr)
    return json.loads(result.stdout)


DE = (
    "de <- data.frame(gene_symbol = c('a','b','c','d','e','f'),"
    " adjusted_p_value = c(0.001, 0.001, 0.001, 0.001, NA, 0.5),"
    " log2_fold_change_raw = c(3, -3, 0.4, 2, 5, 5),"
    " log2_fold_change = c(2.5, -0.1, 0.3, 0.2, 4, 4),"  # b's shrunken value collapsed
    " stringsAsFactors = FALSE);"
)


def test_the_foreground_uses_the_unshrunk_effect_and_the_adjusted_p():
    out = _r(DE + "cat(jsonlite::toJSON(ora_foreground(de, fdr = 0.05, min_abs_log2fc = 1), auto_unbox = FALSE))")
    assert sorted(out["up"]) == ["a", "d"]
    assert out["down"] == ["b"]  # kept although its displayed (shrunken) effect collapsed to -0.1


def test_a_zero_threshold_classifies_by_sign_alone():
    out = _r(DE + "cat(jsonlite::toJSON(ora_foreground(de, 0.05, 0), auto_unbox = FALSE))")
    assert sorted(out["up"]) == ["a", "c", "d"]
    assert out["down"] == ["b"]


def test_the_background_is_the_annotated_genes_of_each_collection():
    code = (
        "sets <- list(go1 = c('a','b','x1','x2'), go2 = c('a','c','x3'), cu1 = c('d','e'));"
        "source_of <- c(go1 = 'go', go2 = 'go', cu1 = 'custom');"
        "tested <- c('a','b','c','d','e','f','g','h','ENSMUSG1');"
        "res <- run_ora(sets, selected = c('a','b'), universe = tested, direction = 'up', set_source = source_of);"
        "cat(jsonlite::toJSON(res[, c('pathway','count','selected_genes','pathway_genes','universe_genes')], dataframe = 'rows'))"
    )
    rows = {row["pathway"]: row for row in _r(code)}
    # go universe = tested genes in go1 or go2 = {a, b, c}; custom universe = {d, e}
    assert rows["go1"]["universe_genes"] == 3 and rows["go2"]["universe_genes"] == 3
    assert rows["cu1"]["universe_genes"] == 2
    # the foreground is restricted to the same background: both a and b are annotated
    assert rows["go1"]["selected_genes"] == 2 and rows["go1"]["count"] == 2
    # members count only tested genes: x1/x2 were not tested so go1 has 2 members
    assert rows["go1"]["pathway_genes"] == 2


def test_the_hypergeometric_p_is_computed_on_that_background():
    code = (
        "sets <- list(s = c('a','b','c'), t = c('d','e','f','g'));"
        "tested <- c('a','b','c','d','e','f','g','z1','z2','z3','z4','z5');"
        "res <- run_ora(sets, c('a','b','d'), tested, 'up', set_source = c(s = 'go', t = 'go'));"
        "cat(jsonlite::toJSON(list(p = res$p_value[res$pathway == 's'], N = res$universe_genes[1]), auto_unbox = TRUE, digits = NA))"
    )
    out = _r(code)
    from math import comb

    N, K, k, x = 7, 3, 3, 2  # 7 annotated genes, 3 in the set, 3 selected, 2 overlap
    expected = sum(comb(K, i) * comb(N - K, k - i) for i in range(x, min(K, k) + 1)) / comb(N, k)
    assert out["N"] == 7
    assert out["p"] == pytest.approx(expected)


def test_a_single_factor_difference_names_its_two_levels():
    # The volcano title and axis read "<numerator> vs <denominator>"; handing them the
    # estimand label printed "<group> vs <reference> vs )".
    out = _r(
        "cat(jsonlite::toJSON(list("
        "pair = pairwise_levels('condition', c(treated = 1, reference = -1)),"
        "composite = is.null(pairwise_levels(c('a','b'), c('x|y' = 1, 'x|z' = -1))),"
        "weighted = is.null(pairwise_levels('condition', c(p = 0.5, q = 0.5, r = -1)))"
        "), auto_unbox = TRUE))"
    )
    assert out["pair"] == {"numerator": "treated", "denominator": "reference"}
    assert out["composite"] is True and out["weighted"] is True


def test_the_background_comes_from_every_annotated_gene_not_only_those_in_testable_sets():
    # A set-size window (here: 2..3 genes) decides which terms are *tested*; it must not
    # shrink the set of genes the collection can annotate. 'big' has 5 genes (outside the
    # window) and contributes y1, y2 to the background; 'tiny' has 1 gene (outside) -> t1.
    code = (
        "annotation <- list(s = c('a','b','c'), big = c('a','y1','y2','y3','y4'), tiny = c('t1'));"
        "testable <- annotation['s'];"
        "src <- c(s = 'go', big = 'go', tiny = 'go');"
        "tested <- c('a','b','c','y1','y2','y3','y4','t1','z');"
        "res <- run_ora(testable, c('a','b'), tested, 'up', set_source = src, annotation_sets = annotation);"
        "cat(jsonlite::toJSON(list(N = res$universe_genes, selected = res$selected_genes), auto_unbox = TRUE))"
    )
    out = _r(code)
    assert out["N"] == 8  # a b c y1 y2 y3 y4 t1 -- not just {a, b, c}
    assert out["selected"] == 2
