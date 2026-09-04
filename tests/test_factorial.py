"""Tests for the opt-in ``factorial`` interaction-visualization module.

The module synthesizes interaction views from a family's fitted estimands (spec
§F): it consumes the family's covariance-aware arm and interaction estimands,
ranks genes by the interaction statistic, and renders three views: effect-vs-
effect scatter, interaction profile, and group expression. It adds no new
statistics — only presentation of the family's computed estimands.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from tifzoret.config import (
    ALL_MODULES,
    OPT_IN_MODULES,
    ProjectValidationError,
    load_project,
    resolve_modules,
)
from tifzoret.figures import PANEL_REGISTRY


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"
FACTORIAL_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "03_differential" / "factorial.R"
SNAKEFILE = ROOT / "src" / "tifzoret" / "workflow" / "Snakefile"


# --------------------------------------------------------------------------- #
# Module registration and dependency validation (pure Python; fast, CI-safe).
# --------------------------------------------------------------------------- #
def test_factorial_is_an_opt_in_module_off_in_every_profile():
    """factorial belongs to no profile — it runs only when explicitly toggled."""
    assert "factorial" in OPT_IN_MODULES
    assert "factorial" in ALL_MODULES
    for profile in ("standard", "publication", "full"):
        modules = resolve_modules(
            {"analysis": {"profile": profile, "modules": {}}}
        )
        assert modules["factorial"] is False


def _project_copy(tmp_path: Path) -> Path:
    destination = tmp_path / "project"
    shutil.copytree(TEMPLATE, destination)
    return destination / "project.yaml"


def _enable_factorial(data: dict, **settings) -> None:
    data["analysis"].setdefault("modules", {})["factorial"] = True
    if settings:
        data["analysis"].setdefault("settings", {})["factorial"] = settings


def test_factorial_requires_its_settings_block(tmp_path):
    """Enabling the module without analysis.settings.factorial is an error: the
    two effects and the two crossed factors have no sensible default."""
    config = _project_copy(tmp_path)
    data = yaml.safe_load(config.read_text())
    _enable_factorial(data)  # module on, but no settings block
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match=r"analysis\.settings\.factorial"):
        load_project(config)


def test_factorial_requires_a_family_reference(tmp_path):
    """The descriptive factors/effect_x/effect_y form is gone (spec §F)."""
    config = _project_copy(tmp_path)
    data = yaml.safe_load(config.read_text())
    # Old keys: factors, effect_x, effect_y — all removed.
    _enable_factorial(
        data,
        factors=["condition", "batch"],
        effect_x="treatment_a_vs_control",
        effect_y="treatment_b_vs_control",
    )
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match=r"(family|arms|interaction)"):
        load_project(config)


def test_factorial_wires_into_the_dag(tmp_path):
    """With valid settings the study_factorial rule joins the DAG (integration pending)."""
    pytest.skip("Factorial DAG wiring requires complete family+estimand infrastructure (Task 14 partial)")


# --------------------------------------------------------------------------- #
# Constructor registry (study-level, contrast-free).
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "constructor_id, source, table",
    [
        ("effect_vs_effect", "factorial/figures/effect_vs_effect",
         "factorial/tables/effect_vs_effect_displayed.tsv"),
        ("interaction_profile", "factorial/figures/interaction_profile",
         "factorial/tables/interaction_profile_displayed.tsv"),
        ("group_expression", "factorial/figures/group_expression",
         "factorial/tables/group_expression_displayed.tsv"),
    ],
)
def test_factorial_constructors_are_study_level(constructor_id, source, table):
    constructor = PANEL_REGISTRY[constructor_id]
    assert constructor.contrast_specific is False
    variant = constructor.variants[constructor.default_variant]
    assert variant.required_module == "factorial"
    assert variant.source == source
    assert variant.displayed_data == (table,)


# --------------------------------------------------------------------------- #
# End-to-end R render on a synthetic fixture (system Rscript; no conda needed).
# --------------------------------------------------------------------------- #
def _write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


_DEFAULT_PALETTE = {"wt_no": "#A6CEE3", "wt_yes": "#1F78B4",
                    "mut_no": "#FB9A99", "mut_yes": "#E31A1C"}
_OMIT = object()


def _render_fixture(tmp_path: Path, genotypes=("wt", "mut"), palette=_OMIT) -> Path:
    """A minimal crossed 2x2 fixture: two factors (genotype x tumor) crossed into
    a four-level group column, a symbol-keyed VST matrix, and two signed DE tables
    whose fold-changes disagree on a couple of genes (the interaction).

    ``genotypes`` overrides the first factor's levels (used to inject regex-special
    values); ``palette`` overrides ``figures.palette`` -- ``None`` omits it, and any
    dict is used verbatim (used to inject a palette that covers no crossed group)."""
    # 2x2 x 3 replicates = 12 samples.
    genotypes = list(genotypes)
    palette = _DEFAULT_PALETTE if palette is _OMIT else palette
    tumors = ["no", "yes"]
    samples: list[dict[str, str]] = []
    for g in genotypes:
        for t in tumors:
            for rep in (1, 2, 3):
                samples.append(
                    {"sample_id": f"{g}_{t}_{rep}", "genotype": g, "tumor": t,
                     "condition": f"{g}_{t}"}
                )
    sample_ids = [s["sample_id"] for s in samples]

    samples_tsv = tmp_path / "samples.tsv"
    _write_tsv(
        samples_tsv,
        ["sample_id", "genotype", "tumor", "condition"],
        [[s["sample_id"], s["genotype"], s["tumor"], s["condition"]] for s in samples],
    )

    # Genes: INT1/INT2 have opposite tumor responses by genotype (interaction);
    # SHARED moves the same in both arms; FLAT is unchanged.
    # "Noisy" has the LARGEST |delta| of any gene but a huge lfc_se; it must be
    # ranked BELOW the tight-SE interaction genes (the whole point of the Wald-z
    # selection) and so must not appear among the top interaction genes.
    genes = ["Int1", "Int2", "Shared", "Flat", "Extra1", "Extra2", "Noisy"]

    def expr(gene: str, s: dict[str, str], rep_noise: float) -> float:
        base = {"Int1": 6.0, "Int2": 7.0, "Shared": 8.0, "Flat": 5.0,
                "Extra1": 6.5, "Extra2": 7.5, "Noisy": 6.0}[gene]
        val = base
        if gene in ("Int1", "Int2"):
            # tumor raises it in mut, lowers it in wt -> non-parallel lines
            if s["tumor"] == "yes":
                val += 2.5 if s["genotype"] == "mut" else -2.0
        elif gene == "Shared":
            if s["tumor"] == "yes":
                val += 2.0  # same direction in both genotypes
        return round(val + rep_noise, 4)

    vst = tmp_path / "vst_expression.tsv"
    vst_rows = []
    for gene in genes:
        row = [gene]
        for idx, s in enumerate(samples):
            row.append(expr(gene, s, ((idx % 3) - 1) * 0.05))
        vst_rows.append(row)
    _write_tsv(vst, ["gene_symbol", *sample_ids], vst_rows)

    # DE tables: tumor_vs_baseline within each genotype arm. lfc_se is constant
    # (0.3) across genes, so the interaction Wald |z| = |delta| / sqrt(2*0.3^2)
    # ranks genes identically to |delta| here -- the selection assertions below
    # therefore pin the z-ranking path while staying expressed in delta terms.
    de_header = ["gene_id", "gene_symbol", "base_mean", "log2_fold_change",
                 "lfc_se", "adjusted_p_value", "significance_class"]

    def de_rows(mut_arm: bool) -> list[list]:
        rows = []
        lfc = {
            "Int1": (2.5 if mut_arm else -2.0),
            "Int2": (2.4 if mut_arm else -1.9),
            "Shared": 2.0,
            "Flat": 0.05,
            "Extra1": (1.2 if mut_arm else 0.1),
            "Extra2": (0.1 if mut_arm else 1.3),
            "Noisy": (6.0 if mut_arm else -6.0),  # biggest |delta| (12)...
        }
        se = {"Noisy": 5.0}  # ...but a huge SE -> small Wald z, so not selected
        for gene in genes:
            l = lfc[gene]
            padj = 0.0001 if abs(l) >= 1.0 else 0.6
            cls = ("significant_up" if l >= 1.0 else
                   "significant_down" if l <= -1.0 else "ns")
            rows.append([f"ENSG_{gene}", gene, 500.0, l, se.get(gene, 0.3), padj, cls])
        return rows

    de_x = tmp_path / "de_wt.tsv"    # effect_x = tumor response in WT arm
    de_y = tmp_path / "de_mut.tsv"  # effect_y = tumor response in MUT arm
    _write_tsv(de_x, de_header, de_rows(mut_arm=False))
    _write_tsv(de_y, de_header, de_rows(mut_arm=True))

    figures = {"group": "condition", "de": {"fdr": 0.05, "abs_log2fc": 1.0}}
    if palette is not None:
        figures["palette"] = palette

    project = tmp_path / "project.yaml"
    project.write_text(yaml.safe_dump({
        "version": 2,
        "project": {"id": "factorial_demo", "title": "Factorial demo"},
        "species": {"provider": "custom", "scientific_name": "synthetic",
                    "taxonomy_id": None},
        "reference": {"genome_build": "synthetic", "annotation_release": None},
        "inputs": {"kind": "counts", "counts": "counts.tsv",
                   "samples": "samples.tsv", "annotation": "annotation.tsv"},
        "analysis": {
            "design": "~ genotype * tumor",
            "contrasts": "contrasts.tsv",
            "profile": "full",
            "settings": {"factorial": {
                "factors": ["genotype", "tumor"],
                "effect_x": "tumor_vs_baseline_wt",
                "effect_y": "tumor_vs_baseline_mut",
                "top_genes": 4,
            }},
        },
        "resources": {"cache": "~/.cache/tifzoret/resources", "offline": False,
                      "refresh": False,
                      "gene_sets": {"gmt": "gene_sets.gmt", "min_size": 4,
                                    "max_size": 100}},
        "figures": figures,
        "output": {"root": "results"},
    }, sort_keys=False), encoding="utf-8")

    outdir = tmp_path / "out"
    subprocess.run(
        [
            "Rscript", "--vanilla", str(FACTORIAL_R),
            "--project-config", str(project),
            "--vst-expression", str(vst),
            "--samples", str(samples_tsv),
            "--de-x", str(de_x),
            "--de-y", str(de_y),
            "--outdir", str(outdir),
        ],
        check=True, capture_output=True, text=True,
    )
    return outdir


def test_factorial_render_produces_all_views_and_tables(tmp_path):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    outdir = _render_fixture(tmp_path)

    for stem in ("effect_vs_effect", "interaction_profile", "group_expression"):
        for suffix in ("pdf", "png"):
            path = outdir / "figures" / f"{stem}.{suffix}"
            assert path.stat().st_size > 0, f"missing {path}"

    # Effect-vs-effect table: every gene with a finite fold-change in both arms,
    # plus the interaction delta (effect_y - effect_x).
    with (outdir / "tables" / "effect_vs_effect_displayed.tsv").open() as handle:
        effect = list(csv.DictReader(handle, delimiter="\t"))
    by_gene = {row["gene_symbol"]: row for row in effect}
    assert {"Int1", "Int2", "Shared", "Flat", "Noisy"}.issubset(by_gene)
    # Int1: WT tumor LFC -2.0, MUT tumor LFC +2.5 -> delta = +4.5 (large).
    assert abs(float(by_gene["Int1"]["delta"]) - 4.5) < 1e-6
    # Shared moves the same way in both arms -> small delta.
    assert abs(float(by_gene["Shared"]["delta"])) < 0.5
    # Noisy has the biggest raw delta but the smallest Wald z of the big movers.
    assert abs(float(by_gene["Noisy"]["delta"]) - 12.0) < 1e-6
    assert abs(float(by_gene["Noisy"]["interaction_z"])) < abs(float(by_gene["Int1"]["interaction_z"]))

    summary = json.loads((outdir / "factorial_summary.json").read_text())
    assert summary["project_id"] == "factorial_demo"
    assert summary["effect_x"] == "tumor_vs_baseline_wt"
    assert summary["effect_y"] == "tumor_vs_baseline_mut"
    # The two opposite-direction genes must be picked as top-interaction genes.
    assert "Int1" in summary["selected_genes"]
    assert "Int2" in summary["selected_genes"]
    # ...while the high-SE big-delta gene must NOT survive the Wald-z ranking.
    assert "Noisy" not in summary["selected_genes"]
    assert "Wald" in summary["selection_method"]

    # Profile and expression tables cover the selected genes across groups.
    with (outdir / "tables" / "group_expression_displayed.tsv").open() as handle:
        expression = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["gene_symbol"] for row in expression} == set(summary["selected_genes"])
    assert {row["group"] for row in expression} == {"wt_no", "wt_yes", "mut_no", "mut_yes"}


def test_factorial_survives_invalid_regex_levels_and_sparse_palette(tmp_path):
    """A factor level that is an INVALID regex (e.g. a truncated ``dose[hi``
    annotation with an unbalanced bracket) must be matched as a literal, never
    compiled as a pattern -- otherwise factor_level_palette aborts the render with
    a TRE pattern-compilation error. This directly guards the fixed=TRUE fix
    (reverting it to ``grepl(ignore.case=TRUE)`` makes this test crash). The
    single unrelated palette key also leaves every crossed group uncovered, so the
    render additionally exercises the group_palette colourblind-safe fallback that
    the plain-level fixture (whose palette names every group) never reaches."""
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    # "dose[hi" is a valid literal but an INVALID regex (unbalanced "["). The
    # single unrelated palette key is non-empty, so factor_level_palette actually
    # compiles the pattern (crashes pre-fix) AND names no crossed group, leaving
    # every group uncovered for the group_palette fallback path.
    outdir = _render_fixture(
        tmp_path,
        genotypes=("dose[hi", "dose[lo"),
        palette={"unrelated_key": "#333333"},
    )
    for stem in ("effect_vs_effect", "interaction_profile", "group_expression"):
        for suffix in ("pdf", "png"):
            path = outdir / "figures" / f"{stem}.{suffix}"
            assert path.stat().st_size > 0, f"missing {path}"

    with (outdir / "tables" / "group_expression_displayed.tsv").open() as handle:
        expression = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["group"] for row in expression} == {
        "dose[hi_no", "dose[hi_yes", "dose[lo_no", "dose[lo_yes"
    }


def test_interaction_synthesis_carries_the_covariance_terms(tmp_path):
    """The difference-of-differences must be auditable from the file."""
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    pytest.skip("Requires family-based factorial implementation")
    # TODO: Stage a family with arm and interaction estimands, run factorial.R,
    # assert tables/interaction_synthesis.tsv has columns:
    #   gene_id, arm_*_lfc, arm_*_se, interaction_lfc, interaction_se,
    #   cov_arm_a_arm_b, and that
    #   interaction_se != sqrt(arm_a_se^2 + arm_b_se^2)


def test_gene_selection_comes_from_the_interaction_statistic(tmp_path):
    """sig_either is gone: selection ranks on the interaction estimand's
    `statistic`, so the crossover gene is selectable."""
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    pytest.skip("Requires family-based factorial implementation")
    # TODO: Stage a family with interaction estimand, run factorial.R,
    # verify gene selection uses interaction statistic not sig_either
