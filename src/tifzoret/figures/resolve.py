"""Panel resolution and review iteration."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable

from .registry import PANEL_REGISTRY, ResolvedPanel

if TYPE_CHECKING:
    from ..config import ResolvedProject


def resolve_panel(project: "ResolvedProject", panel: dict[str, Any]) -> ResolvedPanel:
    """Resolve a recipe panel to figure and displayed-data artifact paths."""
    constructor_id = panel.get("constructor")
    if not constructor_id:
        source = Path(panel["source"]).expanduser()
        source = source if source.is_absolute() else project.result_root / source
        source = source.with_suffix("") if source.suffix.lower() in {".pdf", ".png"} else source
        displayed = []
        for value in panel.get("displayed_data", []):
            path = Path(value).expanduser()
            displayed.append(path if path.is_absolute() else project.result_root / path)
        return ResolvedPanel(None, panel.get("variant"), panel.get("contrast"), panel.get("title", "Legacy source panel"), source, tuple(displayed), None)

    constructor = PANEL_REGISTRY[constructor_id]
    variant_id = panel.get("variant", constructor.default_variant)
    variant = constructor.variants[variant_id]
    contrast: str | None = panel.get("contrast")
    if constructor.contrast_specific and contrast is None:
        contrast = project.contrast_rows[0]["contrast_id"]
    values = {"contrast": contrast or ""}
    return ResolvedPanel(
        constructor_id,
        variant_id,
        contrast,
        panel.get("title", variant.label),
        project.result_root / variant.source.format_map(values),
        tuple(project.result_root / value.format_map(values) for value in variant.displayed_data),
        variant.required_module,
    )


def _selected_keys(project: "ResolvedProject") -> set[tuple[str, str, str]]:
    selected: set[tuple[str, str, str]] = set()
    for recipe in (project.recipe_config or {}).get("figure_sets", {}).values():
        for panel in recipe.get("panels", []):
            if panel.get("constructor"):
                resolved = resolve_panel(project, panel)
                selected.add((resolved.constructor or "", resolved.variant or "", resolved.contrast or ""))
    return selected


def iter_review_panels(project: "ResolvedProject") -> Iterable[dict[str, Any]]:
    """Yield every registered variant supported by the enabled project modules."""
    selected = _selected_keys(project)
    contrast_ids = [row["contrast_id"] for row in project.contrast_rows]
    for constructor in PANEL_REGISTRY.values():
        targets: list[str | None] = contrast_ids if constructor.contrast_specific else [None]
        for contrast in targets:
            for variant_id, variant in constructor.variants.items():
                if variant.required_module not in project.modules:
                    continue
                values = {"contrast": contrast or ""}
                source = project.result_root / variant.source.format_map(values)
                displayed = [project.result_root / value.format_map(values) for value in variant.displayed_data]
                yield {
                    "constructor": constructor.id,
                    "constructor_label": constructor.label,
                    "description": constructor.description,
                    "variant": variant_id,
                    "variant_label": variant.label,
                    "contrast": contrast,
                    "selected": (constructor.id, variant_id, contrast or "") in selected,
                    "source_pdf": str(source.with_suffix(".pdf")),
                    "source_png": str(source.with_suffix(".png")),
                    "displayed_data": [str(path) for path in displayed],
                    "available": source.with_suffix(".pdf").is_file() and source.with_suffix(".png").is_file(),
                }
