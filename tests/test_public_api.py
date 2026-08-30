"""Freeze Tifzoret's public import surface.

Every name here is imported by the engine, its tests, or a workflow script
today (verified 2026-08-29). The config.py -> config/ and figures.py -> figures/
package splits, and the top-level namespace re-export, must keep all of them
importable from the same dotted paths. If this test fails, a split moved or
dropped a public name — fix the __init__ re-exports, do not edit this test.
"""
import importlib


def test_config_public_surface():
    config = importlib.import_module("tifzoret.config")
    for name in (
        "load_project",
        "ProjectValidationError",
        "validation_report",
        "ResolvedProject",
        "_resolve",              # de-facto public: imported by materialize_inputs.py
    ):
        assert hasattr(config, name), f"tifzoret.config.{name} missing"


def test_figures_public_surface():
    figures = importlib.import_module("tifzoret.figures")
    for name in (
        "resolve_panel",
        "PANEL_REGISTRY",
        "build_gallery",
        "normalized_gene_panels",
    ):
        assert hasattr(figures, name), f"tifzoret.figures.{name} missing"


def test_submodule_import_forms_still_resolve():
    # 'from tifzoret import config as cfg' and dotted access must both work.
    from tifzoret import config as cfg
    from tifzoret import figures as figs
    assert cfg.load_project is not None
    assert figs.resolve_panel is not None


def test_top_level_namespace():
    import tifzoret as tfz
    for name in (
        "load_project", "ProjectValidationError", "validation_report", "ResolvedProject",
        "resolve_panel", "PANEL_REGISTRY", "build_gallery", "normalized_gene_panels",
    ):
        assert hasattr(tfz, name), f"tifzoret.{name} missing from top-level namespace"
