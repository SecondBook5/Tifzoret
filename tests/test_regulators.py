"""Tests for the regulators stage, including the opt-in dual-view wiring.

The stage always renders a signed regulator-activity view (VIPER on a signed
regulon; DoRothEA or a custom ``regulon_edges`` table). When
``resources.binding_prior_edges`` is configured it renders a SECOND, unsigned
view in the same run -- a binding prior (e.g. a GTRD occupancy snapshot)
establishes occupancy, not direction, so its target modes are collapsed to +1
and its method is tagged with the provider (VIPER -> VIPER_unsigned_GTRD). One
``tifzoret run`` then regenerates both regulator panels from a single source of
truth.

Two guarantees are tested end-to-end on a synthetic fixture (system Rscript, no
conda; viper/dorothea absent -> the deterministic signed-weighted fallback):

1. Adding a binding prior does not perturb the primary view -- its
   ``regulator_differential.tsv`` is byte-identical with and without the prior.
2. The binding view is produced, is unsigned (every edge ``mor`` == 1), and its
   method string carries the provider tag.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from tifzoret.config import ProjectValidationError, load_project


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"
REGULATORS_R = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "06_regulators" / "regulators.R"


# --------------------------------------------------------------------------- #
# Config validation (pure Python; fast, CI-safe).
# --------------------------------------------------------------------------- #
def _regulators_config(tmp_path: Path) -> tuple[Path, dict]:
    destination = tmp_path / "project"
    shutil.copytree(TEMPLATE, destination)
    config = destination / "project.yaml"
    data = yaml.safe_load(config.read_text())
    data["species"] = {"provider": "mouse", "scientific_name": "Mus musculus", "taxonomy_id": 10090}
    data["reference"] = {"genome_build": "GRCm39", "annotation_release": 107}
    data["analysis"]["modules"] = {"regulators": True}
    # Primary regulon comes from DoRothEA; the binding prior is the GTRD snapshot.
    data["resources"]["providers"] = {"dorothea": True, "gtrd": True}
    return config, data


def test_binding_prior_requires_source_and_target(tmp_path):
    config, data = _regulators_config(tmp_path)
    (config.parent / "binding.tsv").write_text("regulator\tgene\nTF1\tGene1\n")
    data["resources"]["binding_prior_edges"] = "binding.tsv"
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="binding_prior_edges requires source and target"):
        load_project(config)


def test_binding_prior_satisfies_gtrd_provider_without_primary_regulon(tmp_path):
    """A GTRD-derived binding prior satisfies providers.gtrd; the primary regulon
    is DoRothEA, so no regulon_edges is required and the config loads."""
    config, data = _regulators_config(tmp_path)
    (config.parent / "binding.tsv").write_text("source\ttarget\nTF1\tGene1\n")
    data["resources"]["binding_prior_edges"] = "binding.tsv"
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    project = load_project(config)
    assert project.binding_prior_edges == config.parent / "binding.tsv"
    assert project.regulon_edges is None


def test_binding_prior_missing_file_is_reported(tmp_path):
    config, data = _regulators_config(tmp_path)
    data["resources"]["binding_prior_edges"] = "nope.tsv"
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="binding prior edge file does not exist"):
        load_project(config)


# --------------------------------------------------------------------------- #
# End-to-end R render on a synthetic fixture (system Rscript; no conda needed).
# --------------------------------------------------------------------------- #
FIXTURE_BUILDER_R = r"""
suppressPackageStartupMessages(library(SummarizedExperiment))
root <- commandArgs(trailingOnly = TRUE)[[1]]
dir.create(root, recursive = TRUE, showWarnings = FALSE)
conditions <- rep(c("control", "treated"), each = 3)
sample_ids <- paste0(rep(c("ctrl", "trt"), each = 3), rep(1:3, 2))
samples <- data.frame(sample_id = sample_ids, condition = conditions, stringsAsFactors = FALSE)
write.table(samples, file.path(root, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
n_genes <- 40L
gene_ids <- sprintf("ENSG%05d", 1:n_genes)
gene_symbols <- sprintf("Gene%02d", 1:n_genes)
annotation <- data.frame(gene_id = gene_ids, gene_symbol = gene_symbols, seqname = "1",
  start = 1:n_genes, end = (1:n_genes) + 100L, strand = "+",
  gene_biotype = "protein_coding", stringsAsFactors = FALSE)
write.table(annotation, file.path(root, "annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
mat <- matrix(0.0, nrow = n_genes, ncol = length(sample_ids), dimnames = list(gene_ids, sample_ids))
for (g in 1:n_genes) for (s in seq_along(sample_ids)) {
  treated <- conditions[s] == "treated"
  base <- 6 + (g %% 5) * 0.4 + (s %% 3) * 0.15
  shift <- if (g %% 2 == 0) 0.8 else -0.6
  mat[g, s] <- base + (if (treated) shift else 0) + sin(g + s) * 0.05
}
saveRDS(SummarizedExperiment(assays = list(vst = mat)), file.path(root, "vst.rds"))
reg_rows <- list(); regs <- c("REGA", "REGB", "REGC")
for (i in seq_along(regs)) {
  targets <- gene_symbols[((i - 1) * 8 + 1):((i - 1) * 8 + 8)]
  reg_rows[[i]] <- data.frame(source = regs[i], target = targets,
    mor = rep(c(1, -1), length.out = length(targets)), stringsAsFactors = FALSE)
}
regulon <- do.call(rbind, reg_rows)
write.table(regulon, file.path(root, "regulon_edges_input.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# Binding prior carries a SIGNED mor here on purpose: the engine must COLLAPSE it
# to +1 (occupancy, not direction). Writing signed input makes the unsigned-view
# assertion below non-tautological -- mor==1 in the output proves the collapse
# fired, not merely that a mor column was absent and defaulted to 1.
write.table(regulon[, c("source", "target", "mor")], file.path(root, "binding_prior_input.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
contrasts <- data.frame(contrast_id = "treated_vs_control", factor = "condition",
  numerator = "treated", denominator = "control", stringsAsFactors = FALSE)
write.table(contrasts, file.path(root, "contrasts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines("DUMMY\tna\tGene01\tGene02", file.path(root, "sets.gmt"))
writeLines("gene_id\tctrl1", file.path(root, "counts.tsv"))
"""

_BASE_YAML = {
    "version": 2,
    "project": {"id": "fixture", "title": "fixture", "description": "fixture"},
    "species": {"provider": "custom", "scientific_name": "synthetic", "taxonomy_id": None},
    "reference": {"genome_build": "synthetic", "annotation_release": None},
    "inputs": {"kind": "counts", "samples": "samples.tsv", "counts": "counts.tsv", "annotation": "annotation.tsv"},
    "analysis": {
        "design": "~condition",
        "contrasts": "contrasts.tsv",
        "profile": "standard",
        "settings": {"regulators": {"confidence": ["A", "B", "C"], "min_targets": 5, "max_targets": 1000, "top_regulators": 15}},
    },
    "resources": {
        "gene_sets": {"gmt": "sets.gmt"},
        "regulon_edges": "regulon_edges_input.tsv",
        "providers": {"dorothea": False, "gtrd": False},
    },
    "figures": {"group": "condition"},
    "output": {"root": "out"},
}


def _build_fixture(tmp_path: Path) -> Path:
    fixture = tmp_path / "fixture"
    builder = tmp_path / "build.R"
    builder.write_text(FIXTURE_BUILDER_R)
    subprocess.run(["Rscript", "--vanilla", str(builder), str(fixture)], check=True, capture_output=True, text=True)

    import copy

    primary = copy.deepcopy(_BASE_YAML)
    (fixture / "project_primary.yaml").write_text(yaml.safe_dump(primary, sort_keys=False))

    dual = copy.deepcopy(_BASE_YAML)
    dual["resources"]["binding_prior_edges"] = "binding_prior_input.tsv"
    dual["resources"]["providers"] = {"dorothea": False, "gtrd": True}
    (fixture / "project_dual.yaml").write_text(yaml.safe_dump(dual, sort_keys=False))
    return fixture


def _run_regulators(fixture: Path, config_name: str, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [
            "Rscript", "--vanilla", str(REGULATORS_R),
            "--project-config", str(fixture / config_name),
            "--samples", str(fixture / "samples.tsv"),
            "--annotation", str(fixture / "annotation.tsv"),
            "--contrasts", str(fixture / "contrasts.tsv"),
            "--contrast-id", "treated_vs_control",
            "--vst", str(fixture / "vst.rds"),
            "--outdir", str(outdir),
        ],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout + "\n" + proc.stderr)
    assert proc.returncode == 0, proc.stderr


def _rscript_has_viper() -> bool:
    """True iff the Rscript on PATH can load viper -- the CI standard profile uses
    a system Rscript without it (fallback path), while the pinned r.yaml conda env
    ships viper (canonical aREA path). The VIPER-backend test below is gated on
    this so it runs wherever viper is installed and skips cleanly where it is not."""
    if shutil.which("Rscript") is None:
        return False
    proc = subprocess.run(
        ["Rscript", "--vanilla", "-e",
         'quit(status = as.integer(!requireNamespace("viper", quietly = TRUE)))'],
        capture_output=True,
    )
    return proc.returncode == 0


@pytest.mark.skipif(not _rscript_has_viper(), reason="viper not installed for this Rscript")
def test_viper_backend_is_used_when_installed(tmp_path):
    """When viper IS importable, the stage MUST score with VIPER's aREA -- never
    silently downgrade to the deterministic signed-weighted proxy. The fixtures
    above run under a viper-absent Rscript and only exercise the fallback, so
    without this test the canonical (and scientifically load-bearing) scoring path
    is covered by nothing. Complements regulators.R's fail-loud handling: a viper
    error must stop the run, not quietly fall back. The method string is the tell
    -- the proxy path never writes "VIPER"."""
    fixture = _build_fixture(tmp_path)
    dual_dir = tmp_path / "dual_viper"
    _run_regulators(fixture, "project_dual.yaml", dual_dir)

    import json

    primary = json.loads((dual_dir / "regulators_summary.json").read_text())
    binding = json.loads((dual_dir / "regulators_binding_summary.json").read_text())
    assert primary["method"] == "VIPER", primary["method"]
    assert binding["method"] == "VIPER_unsigned_GTRD", binding["method"]


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")
def test_strict_mode_fails_on_non_viper_fallback(tmp_path):
    """Under a viper-absent Rscript the stage would silently score regulator
    activity with the deterministic signed-weighted PROXY. With execution.strict
    the run MUST stop instead -- a publication regulator panel may not be built on
    non-VIPER scores. Skipped where viper IS installed, because then the canonical
    aREA path runs and there is no degradation for strict mode to trip on."""
    if _rscript_has_viper():
        pytest.skip("viper installed: canonical VIPER path runs, no fallback to trip strict mode")
    fixture = _build_fixture(tmp_path)

    import copy

    strict = copy.deepcopy(_BASE_YAML)
    strict["execution"] = {"strict": True}
    (fixture / "project_strict.yaml").write_text(yaml.safe_dump(strict, sort_keys=False))

    outdir = tmp_path / "strict"
    outdir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [
            "Rscript", "--vanilla", str(REGULATORS_R),
            "--project-config", str(fixture / "project_strict.yaml"),
            "--samples", str(fixture / "samples.tsv"),
            "--annotation", str(fixture / "annotation.tsv"),
            "--contrasts", str(fixture / "contrasts.tsv"),
            "--contrast-id", "treated_vs_control",
            "--vst", str(fixture / "vst.rds"),
            "--outdir", str(outdir),
        ],
        capture_output=True, text=True,
    )
    assert proc.returncode != 0, "strict mode must fail on the non-VIPER fallback"
    assert "strict mode" in proc.stderr, proc.stderr
    # The stage fails BEFORE emitting a summary, so no degraded result is left behind.
    assert not (outdir / "regulators_summary.json").exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")
def test_dual_view_leaves_primary_identical_and_adds_unsigned_view(tmp_path):
    fixture = _build_fixture(tmp_path)
    primary_dir = tmp_path / "primary_only"
    dual_dir = tmp_path / "dual"
    _run_regulators(fixture, "project_primary.yaml", primary_dir)
    _run_regulators(fixture, "project_dual.yaml", dual_dir)

    # (1) Primary view is untouched by the presence of a binding prior.
    primary_diff = (primary_dir / "tables" / "regulator_differential.tsv").read_bytes()
    dual_primary_diff = (dual_dir / "tables" / "regulator_differential.tsv").read_bytes()
    assert primary_diff == dual_primary_diff

    # (2) Primary-only run emits NO binding artifacts.
    assert not (primary_dir / "tables" / "regulator_differential_binding.tsv").exists()

    # (3) Dual run emits the full binding view.
    for name in (
        "tables/regulon_edges_binding.tsv",
        "tables/regulator_differential_binding.tsv",
        "tables/regulator_activity_scores_binding.tsv",
        "figures/regulator_activity_binding.pdf",
        "figures/regulator_activity_binding.png",
        "regulators_binding_summary.json",
    ):
        assert (dual_dir / name).exists(), name

    # (4) Binding view is unsigned: every configured edge mor collapsed to +1.
    binding_edges = (dual_dir / "tables" / "regulon_edges_binding.tsv").read_text().splitlines()
    header = binding_edges[0].split("\t")
    mor_idx = header.index("mor")
    provider_idx = header.index("provider")
    mors = {row.split("\t")[mor_idx] for row in binding_edges[1:]}
    providers = {row.split("\t")[provider_idx] for row in binding_edges[1:]}
    assert mors == {"1"}, mors
    assert providers == {"gtrd"}, providers

    # (4b) The collapse is specific to the binding view: the primary (signed) view
    # PRESERVES the fixture's signed mor. With (4), this proves force_unsigned
    # actually fired rather than the input merely lacking a mor column.
    primary_edges = (dual_dir / "tables" / "regulon_edges.tsv").read_text().splitlines()
    p_header = primary_edges[0].split("\t")
    p_mor_idx = p_header.index("mor")
    primary_mors = {row.split("\t")[p_mor_idx] for row in primary_edges[1:]}
    assert primary_mors == {"1", "-1"}, primary_mors

    # (5) Method string carries the unsigned provider tag; primary does not.
    import json

    binding_summary = json.loads((dual_dir / "regulators_binding_summary.json").read_text())
    primary_summary = json.loads((dual_dir / "regulators_summary.json").read_text())
    assert binding_summary["method"].endswith("_unsigned_GTRD"), binding_summary["method"]
    assert "_unsigned_" not in primary_summary["method"], primary_summary["method"]


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")
def test_named_views_render_primary_canonical_and_suffixed_secondary(tmp_path):
    """The generalized analysis.settings.regulators.views list scores an arbitrary
    named set of views in one run. The FIRST view is the primary and keeps the
    canonical filenames (so the GRN/hypothesis stages read it unchanged); every
    later view writes `_<name>` variants. A signed first view + unsigned second
    view is the general form of the built-in dual view, so its primary output is
    byte-identical to the legacy single-view run and its second view is unsigned."""
    fixture = _build_fixture(tmp_path)

    import copy

    views = copy.deepcopy(_BASE_YAML)
    # views mode is mutually exclusive with resources.regulon_edges.
    del views["resources"]["regulon_edges"]
    views["analysis"]["settings"]["regulators"]["views"] = [
        {"name": "primary", "edges": "regulon_edges_input.tsv", "signed": True},
        {"name": "occupancy", "edges": "binding_prior_input.tsv", "signed": False, "provider_label": "gtrd"},
    ]
    (fixture / "project_views.yaml").write_text(yaml.safe_dump(views, sort_keys=False))

    primary_dir = tmp_path / "primary_only"
    views_dir = tmp_path / "views"
    _run_regulators(fixture, "project_primary.yaml", primary_dir)
    _run_regulators(fixture, "project_views.yaml", views_dir)

    # (1) The primary view writes the canonical filenames.
    for name in (
        "tables/regulon_edges.tsv",
        "tables/regulator_differential.tsv",
        "figures/regulator_activity.pdf",
        "regulators_summary.json",
    ):
        assert (views_dir / name).exists(), name

    # (2) Non-primary views write `_<name>` variants (never the canonical names).
    for name in (
        "tables/regulon_edges_occupancy.tsv",
        "tables/regulator_differential_occupancy.tsv",
        "tables/regulator_activity_scores_occupancy.tsv",
        "figures/regulator_activity_occupancy.pdf",
        "figures/regulator_activity_occupancy.png",
        "regulators_occupancy_summary.json",
    ):
        assert (views_dir / name).exists(), name

    # (3) The primary view is byte-identical to the legacy single-view run: the
    # first views entry drives the exact same code path over the same edge table.
    assert (views_dir / "tables" / "regulator_differential.tsv").read_bytes() == (
        primary_dir / "tables" / "regulator_differential.tsv"
    ).read_bytes()

    # (4) The second view is unsigned (every mor collapsed to +1) while the primary
    # preserves the fixture's signed mor.
    occ_edges = (views_dir / "tables" / "regulon_edges_occupancy.tsv").read_text().splitlines()
    occ_header = occ_edges[0].split("\t")
    occ_mor = {row.split("\t")[occ_header.index("mor")] for row in occ_edges[1:]}
    assert occ_mor == {"1"}, occ_mor
    primary_edges = (views_dir / "tables" / "regulon_edges.tsv").read_text().splitlines()
    p_header = primary_edges[0].split("\t")
    primary_mor = {row.split("\t")[p_header.index("mor")] for row in primary_edges[1:]}
    assert primary_mor == {"1", "-1"}, primary_mor

    # (5) Only the unsigned view's method carries the provider tag.
    import json

    occ_summary = json.loads((views_dir / "regulators_occupancy_summary.json").read_text())
    primary_summary = json.loads((views_dir / "regulators_summary.json").read_text())
    assert occ_summary["method"].endswith("_unsigned_GTRD"), occ_summary["method"]
    assert "_unsigned_" not in primary_summary["method"], primary_summary["method"]
