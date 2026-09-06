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


ROOT = Path(__file__).resolve().parents[2]
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


FACTORIAL_TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "factorial"


def test_factorial_wires_into_the_dag(tmp_path):
    """With valid settings the ``study_factorial`` rule joins the DAG.

    This is the test the three layers used to fail independently. ``load_project``
    validates ``settings.factorial.{family,arms,interaction}``; the Snakefile used
    to read ``effect_x``/``effect_y`` instead, so a *valid* config yielded the
    ``__factorial_disabled__`` sentinel and the rule's inputs resolved to paths no
    rule produces -- a "Missing input files" DAG error. A dry-run is the cheapest
    proof that config, Snakefile and rule now agree on the same keys, and that
    every input the rule names is producible (the three estimand DE tables and the
    family's contrast matrix and covariance).
    """
    destination = tmp_path / "factorial"
    shutil.copytree(FACTORIAL_TEMPLATE, destination)
    config_path = destination / "project.yaml"

    data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    data["analysis"].setdefault("modules", {})["factorial"] = True
    # The two arms are the family's two simple effects; their difference is the
    # family's interaction estimand, so all three come from one fit.
    data["analysis"].setdefault("settings", {})["factorial"] = {
        "family": "fam_factorial",
        "arms": ["simple_control", "simple_treated"],
        "interaction": "interaction",
        "factors": ["treatment", "genotype"],
        "top_genes": 4,
    }
    config_path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")

    # The settings must satisfy the validator, not merely the Snakefile.
    load_project(config_path)

    proc = subprocess.run(
        [sys.executable, "-m", "tifzoret", "dry-run", str(config_path),
         "--no-conda", "--cores", "1"],
        check=False, capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr or proc.stdout
    planned = {
        line.strip()[len("rule "):-1]
        for line in proc.stdout.splitlines()
        if line.strip().startswith("rule ") and line.strip().endswith(":")
    }
    assert "study_factorial" in planned, sorted(planned)
    # The interaction synthesis table is a declared output, so deleting it
    # reschedules the rule rather than silently never coming back.
    assert "factorial/tables/interaction_synthesis.tsv" in proc.stdout


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


# The 2x2 model's coefficients, in the order DESeq2 names them for
# ``~ genotype * tumor``. These are opaque labels to factorial.R -- it only joins
# contrast_matrix columns to coefficient_covariance keys -- so they stay fixed
# even when the fixture renames factor levels.
_COEFFICIENTS = ("Intercept", "genotype_mut_vs_wt", "tumor_yes_vs_no",
                 "genotypemut.tumoryes")

# One family, three estimands, expressed as contrast vectors over _COEFFICIENTS:
#   arm A = tumor effect in the WT (reference) genotype  -> tumor_yes_vs_no
#   arm B = tumor effect in the MUT genotype             -> tumor_yes_vs_no + interaction
#   interaction = arm B - arm A                          -> interaction alone
_CONTRAST_VECTORS = {
    "tumor_in_wt": (0.0, 0.0, 1.0, 0.0),
    "tumor_in_mut": (0.0, 0.0, 1.0, 1.0),
    "tumor_by_genotype": (0.0, 0.0, 0.0, 1.0),
}

# A per-gene coefficient covariance S with a deliberately NON-ZERO cross term
# between the two coefficients the arms share. Only the tumor/interaction block
# matters (both arms carry weight 0 on Intercept and genotype), but the fixture
# writes a full 4x4 the way family_fit.R does.
#
# With s33 = Var(tumor), s44 = Var(interaction), s34 = Cov(tumor, interaction):
#   Var(arm A)      = s33
#   Var(arm B)      = s33 + s44 + 2*s34
#   Cov(arm A,arm B)= s33 + s34
#   Var(interaction)= s44
# so Var(A) + Var(B) - 2*Cov(A,B) collapses to exactly s44 -- the identity the
# synthesis table must satisfy. Choosing s34 < 0 makes the naive independent-sum
# SE strictly LARGER than the truth, i.e. the naive statistic understates a real
# interaction, which is the failure mode worth pinning.
_S33, _S44, _S34 = 0.09, 0.18, -0.045
# "Noisy" gets the same shape scaled up: a huge interaction SE, so its very large
# fold-change still yields a small |z| and it must lose the ranking. The scale is
# chosen so Noisy's |z| (12 / sqrt(0.18 * scale) = 1.41) falls strictly BELOW
# every other mover's, with no tie -- at 100x it landed exactly on Extra2's 2.83
# and won the tie-break, which made the exclusion assertion depend on sort
# stability rather than on the statistic.
_NOISY_SCALE = 400.0


def _sigma(scale: float) -> dict[tuple[str, str], float]:
    """A 4x4 covariance keyed by coefficient-name pairs, scaled by ``scale``."""
    intercept, genotype, tumor, interaction = _COEFFICIENTS
    base = {
        (intercept, intercept): 0.04, (genotype, genotype): 0.05,
        (tumor, tumor): _S33, (interaction, interaction): _S44,
        (tumor, interaction): _S34, (interaction, tumor): _S34,
    }
    full = {}
    for row in _COEFFICIENTS:
        for col in _COEFFICIENTS:
            full[(row, col)] = base.get((row, col), 0.0) * scale
    return full


def _quadratic(vector_a, vector_b, sigma) -> float:
    """c_a' S c_b -- the covariance of two contrasts under one fit."""
    return sum(
        vector_a[i] * sigma[(_COEFFICIENTS[i], _COEFFICIENTS[j])] * vector_b[j]
        for i in range(len(_COEFFICIENTS))
        for j in range(len(_COEFFICIENTS))
    )


def _render_fixture(tmp_path: Path, genotypes=("wt", "mut"), palette=_OMIT) -> Path:
    """A minimal crossed 2x2 fixture staged as ONE family's outputs.

    Two factors (genotype x tumor) crossed into a four-level group column, a
    symbol-keyed VST matrix, and the three per-estimand ``de_results.tsv`` tables
    that ``family_estimand`` produces -- both arms and the interaction -- plus the
    family's ``contrast_matrix.tsv`` and ``coefficient_covariance.tsv``. Every SE
    is derived from a single shared covariance S, so the fixture is internally
    consistent the way a real family fit is: the interaction SE genuinely carries
    the cross term, and is NOT the independent sum of the arm SEs.

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

    # ----------------------------------------------------------------------- #
    # The family's three estimand tables. Fold-changes are chosen per gene; every
    # SE is DERIVED from that gene's shared covariance S, so arm and interaction
    # SEs are mutually consistent rather than independently invented.
    # ----------------------------------------------------------------------- #
    de_header = ["gene_id", "gene_symbol", "base_mean", "log2_fold_change",
                 "lfc_se", "statistic", "adjusted_p_value", "significance_class"]

    # Tumor response per gene in the WT arm and the MUT arm. The interaction
    # fold-change is exactly the difference, because the interaction estimand IS
    # arm B - arm A as a contrast vector.
    arm_lfc = {
        "Int1": (-2.0, 2.5),   # opposite directions -> delta +4.5
        "Int2": (-1.9, 2.4),   # opposite directions -> delta +4.3
        "Shared": (2.0, 2.0),  # same in both arms   -> delta 0
        "Flat": (0.05, 0.05),
        "Extra1": (0.1, 1.2),
        "Extra2": (1.3, 0.1),
        "Noisy": (-6.0, 6.0),  # biggest |delta| (12) but the largest SE
    }

    def gene_sigma(gene: str) -> dict[tuple[str, str], float]:
        return _sigma(_NOISY_SCALE if gene == "Noisy" else 1.0)

    def de_rows(estimand: str) -> list[list]:
        vector = _CONTRAST_VECTORS[estimand]
        rows = []
        for gene in genes:
            lfc_wt, lfc_mut = arm_lfc[gene]
            value = {"tumor_in_wt": lfc_wt, "tumor_in_mut": lfc_mut,
                     "tumor_by_genotype": lfc_mut - lfc_wt}[estimand]
            sigma = gene_sigma(gene)
            se = _quadratic(vector, vector, sigma) ** 0.5
            statistic = value / se if se > 0 else 0.0
            padj = 0.0001 if abs(statistic) >= 3.0 else 0.6
            cls = ("significant_up" if value >= 1.0 else
                   "significant_down" if value <= -1.0 else "ns")
            rows.append([f"ENSG_{gene}", gene, 500.0, value, se, statistic, padj, cls])
        return rows

    de_arm_a = tmp_path / "de_tumor_in_wt.tsv"
    de_arm_b = tmp_path / "de_tumor_in_mut.tsv"
    de_interaction = tmp_path / "de_tumor_by_genotype.tsv"
    _write_tsv(de_arm_a, de_header, de_rows("tumor_in_wt"))
    _write_tsv(de_arm_b, de_header, de_rows("tumor_in_mut"))
    _write_tsv(de_interaction, de_header, de_rows("tumor_by_genotype"))

    # contrast_matrix.tsv: one row per estimand, provenance columns then one
    # column per coefficient -- the layout family_fit.R writes.
    provenance = ["estimand_id", "label", "role", "atom_kind", "expression",
                  "cell_weights"]
    _write_tsv(
        tmp_path / "contrast_matrix.tsv",
        [*provenance, *_COEFFICIENTS],
        [
            [estimand, estimand, "primary", "difference", f"<{estimand}>", "",
             *_CONTRAST_VECTORS[estimand]]
            for estimand in ("tumor_in_wt", "tumor_in_mut", "tumor_by_genotype")
        ],
    )

    # coefficient_covariance.tsv: long-form S per gene.
    covariance_rows = []
    for gene in genes:
        sigma = gene_sigma(gene)
        for row_name in _COEFFICIENTS:
            for col_name in _COEFFICIENTS:
                covariance_rows.append(
                    [f"ENSG_{gene}", row_name, col_name, sigma[(row_name, col_name)]]
                )
    _write_tsv(
        tmp_path / "coefficient_covariance.tsv",
        ["gene_id", "coefficient_row", "coefficient_col", "covariance"],
        covariance_rows,
    )

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
                "family": "fam_tumor_genotype",
                "arms": ["tumor_in_wt", "tumor_in_mut"],
                "interaction": "tumor_by_genotype",
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
            "--de-arm-a", str(de_arm_a),
            "--de-arm-b", str(de_arm_b),
            "--de-interaction", str(de_interaction),
            "--contrast-matrix", str(tmp_path / "contrast_matrix.tsv"),
            "--covariance", str(tmp_path / "coefficient_covariance.tsv"),
            "--arms", "tumor_in_wt,tumor_in_mut",
            "--interaction", "tumor_by_genotype",
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
    assert summary["arm_a"] == "tumor_in_wt"
    assert summary["arm_b"] == "tumor_in_mut"
    assert summary["interaction"] == "tumor_by_genotype"
    # The covariance term was actually resolved for every gene, not silently
    # dropped -- otherwise the audit column would be all-NA and the identity test
    # below would be vacuous.
    assert summary["covariance_genes"] == len(
        {"Int1", "Int2", "Shared", "Flat", "Extra1", "Extra2", "Noisy"}
    ), summary["covariance_genes"]
    # The two opposite-direction genes must be picked as top-interaction genes.
    assert "Int1" in summary["selected_genes"]
    assert "Int2" in summary["selected_genes"]
    # ...while the high-SE big-delta gene must NOT survive the Wald-z ranking.
    assert "Noisy" not in summary["selected_genes"]
    assert "statistic" in summary["selection_method"]  # ranks by interaction |statistic|

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
    """The difference-of-differences must be auditable from the file.

    This is the test behind the engine's headline methodological claim. Two
    contrasts fitted independently would be differenced by adding variances:
    ``sqrt(arm_a_se^2 + arm_b_se^2)``. Arms extracted from ONE family fit are
    correlated, so the truth carries a cross term:

        Var(c_b'B - c_a'B) = c_a'S c_a + c_b'S c_b - 2 c_a'S c_b

    The assertions below pin all three properties that matter: the cross term is
    published, it is non-zero, and ``interaction_se`` actually reflects it (i.e.
    it differs from the naive sum by exactly the cross term). If the module ever
    regresses to differencing two independent tables, the reconstruction identity
    breaks and so does this test.
    """
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    outdir = _render_fixture(tmp_path)

    with (outdir / "tables" / "interaction_synthesis.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows, "interaction_synthesis.tsv is empty"

    required = {"gene_id", "gene_symbol", "arm_a_lfc", "arm_a_se", "arm_b_lfc",
                "arm_b_se", "interaction_lfc", "interaction_se", "se_naive",
                "cov_arm_a_arm_b"}
    assert required <= set(rows[0]), sorted(required - set(rows[0]))

    by_gene = {row["gene_symbol"]: row for row in rows}
    for gene in ("Int1", "Shared", "Noisy"):
        row = by_gene[gene]
        arm_a_se = float(row["arm_a_se"])
        arm_b_se = float(row["arm_b_se"])
        interaction_se = float(row["interaction_se"])
        covariance = float(row["cov_arm_a_arm_b"])
        naive = (arm_a_se ** 2 + arm_b_se ** 2) ** 0.5

        # (1) The cross term is present and genuinely non-zero -- a zero here
        # would make the remaining assertions pass trivially.
        assert abs(covariance) > 1e-9, f"{gene}: cross term is zero"

        # (2) The covariance-aware SE is NOT the independent sum. With the
        # fixture's negative cross term the naive value is strictly larger, so a
        # naive statistic would understate this interaction.
        assert abs(interaction_se - naive) > 1e-6, (
            f"{gene}: interaction_se {interaction_se} == naive {naive}")
        assert interaction_se < naive, f"{gene}: {interaction_se} !< {naive}"

        # (3) The exact identity, recomputable by a reader from these columns.
        reconstructed = (naive ** 2 - 2 * covariance) ** 0.5
        assert abs(reconstructed - interaction_se) < 1e-6, (
            f"{gene}: sqrt(naive^2 - 2*cov) = {reconstructed} != "
            f"interaction_se = {interaction_se}")
        # ...and the script publishes that reconstruction itself.
        assert abs(float(row["se_reconstructed"]) - interaction_se) < 1e-6


def test_gene_selection_comes_from_the_interaction_statistic(tmp_path):
    """sig_either is gone: selection ranks on the interaction estimand's
    `statistic`, so the crossover gene is selectable.

    ``Noisy`` is the discriminating case. It has by far the largest raw
    difference-of-differences (delta = 12, three times Int1's 4.5) but a
    proportionally huge covariance, so its covariance-aware |z| is the *smallest*
    of the big movers. Ranking on delta -- or on either arm's significance --
    selects it; ranking on the interaction statistic does not.
    """
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    outdir = _render_fixture(tmp_path)

    summary = json.loads((outdir / "factorial_summary.json").read_text())
    selected = summary["selected_genes"]
    assert summary["selection_method"] == (
        "top interaction |statistic| from the covariance-aware family fit"
    ), summary["selection_method"]

    with (outdir / "tables" / "effect_vs_effect_displayed.tsv").open() as handle:
        effect = list(csv.DictReader(handle, delimiter="\t"))
    by_gene = {row["gene_symbol"]: row for row in effect}

    # Premise: Noisy really does dominate on the descriptive delta, so excluding
    # it can only be the statistic's doing.
    assert abs(float(by_gene["Noisy"]["delta"])) > abs(float(by_gene["Int1"]["delta"]))
    # ...and really does lose on the covariance-aware statistic.
    assert abs(float(by_gene["Noisy"]["interaction_z"])) < abs(
        float(by_gene["Int1"]["interaction_z"]))

    assert "Int1" in selected and "Int2" in selected, selected
    assert "Noisy" not in selected, selected

    # The published interaction_z must BE the covariance-aware statistic, not the
    # naive one. They differ for every gene here (non-zero cross term), so this
    # distinguishes the two -- it was reporting the naive value while sorting by
    # the covariance-aware one.
    row = by_gene["Int1"]
    assert abs(float(row["interaction_z"])) > abs(float(row["interaction_z_naive"]))
    assert abs(float(row["interaction_z"])
               - float(row["interaction_lfc"]) / float(row["interaction_se"])) < 1e-6

    # Rows are ordered by the same statistic the selection uses.
    magnitudes = [abs(float(r["interaction_z"])) for r in effect
                  if r["interaction_z"] not in ("", "NA")]
    assert magnitudes == sorted(magnitudes, reverse=True), magnitudes
