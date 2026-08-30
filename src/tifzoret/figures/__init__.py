"""Registered publication panels and hypothesis-driven review tooling."""

from .gallery import build_gallery
from .recipe import (
    initialize_figure_workflow,
    normalized_gene_panels,
    validate_recipe_contract,
)
from .registry import (
    PANEL_REGISTRY,
    PanelConstructor,
    PanelVariant,
    ResolvedPanel,
    constructor_catalog,
)
from .resolve import iter_review_panels, resolve_panel

__all__ = [
    "PANEL_REGISTRY",
    "PanelConstructor",
    "PanelVariant",
    "ResolvedPanel",
    "build_gallery",
    "constructor_catalog",
    "initialize_figure_workflow",
    "iter_review_panels",
    "normalized_gene_panels",
    "resolve_panel",
    "validate_recipe_contract",
]
