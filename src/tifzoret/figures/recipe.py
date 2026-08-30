"""Recipe normalization, validation, and initialization."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Iterable

import yaml

from .registry import PANEL_REGISTRY


def normalized_gene_panels(panel_config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Expose first-class programs through the legacy grouped-panel shape."""
    panels = {key: dict(value) for key, value in panel_config.get("gene_panels", {}).items()}
    for program_id, program in panel_config.get("programs", {}).items():
        label = program.get("label", program_id.replace("_", " ").title())
        panels[program_id] = {
            "description": program.get("description", label),
            "color": program.get("color"),
            "groups": {label: list(program["genes"])},
        }
    return panels


def validate_recipe_contract(
    recipe_config: dict[str, Any], modules: Iterable[str], contrast_ids: Iterable[str]
) -> list[str]:
    """Validate registered constructors beyond the structural JSON schema."""
    enabled = set(modules)
    contrasts = set(contrast_ids)
    errors: list[str] = []
    for figure_set, recipe in recipe_config.get("figure_sets", {}).items():
        for panel in recipe.get("panels", []):
            constructor_id = panel.get("constructor")
            if not constructor_id:
                continue
            constructor = PANEL_REGISTRY.get(constructor_id)
            prefix = f"figure set {figure_set!r} panel {panel.get('id')!r}"
            if constructor is None:
                errors.append(f"{prefix}: unknown constructor {constructor_id!r}")
                continue
            variant_id = panel.get("variant", constructor.default_variant)
            variant = constructor.variants.get(variant_id)
            if variant is None:
                errors.append(
                    f"{prefix}: constructor {constructor_id!r} has no variant {variant_id!r}; "
                    f"choose from {sorted(constructor.variants)}"
                )
                continue
            if variant.required_module not in enabled:
                errors.append(
                    f"{prefix}: constructor {constructor_id!r}/{variant_id!r} requires "
                    f"analysis module {variant.required_module!r}"
                )
            contrast = panel.get("contrast")
            if constructor.contrast_specific:
                if contrast is None and len(contrasts) != 1:
                    errors.append(f"{prefix}: contrast is required when the project has multiple contrasts")
                elif contrast is not None and contrast not in contrasts:
                    errors.append(f"{prefix}: unknown contrast {contrast!r}")
            elif contrast is not None:
                errors.append(f"{prefix}: constructor {constructor_id!r} is study-level and does not accept contrast")
    return errors


def _write_yaml(path: Path, value: dict[str, Any], force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"refusing to overwrite existing publication file: {path}")
    path.write_text(yaml.safe_dump(value, sort_keys=False), encoding="utf-8")


def _read_contrasts(project_path: Path, config: dict[str, Any]) -> list[dict[str, str]]:
    value = config.get("analysis", {}).get("contrasts")
    if not value:
        raise ValueError("project configuration has no analysis.contrasts path")
    path = Path(str(value)).expanduser()
    path = path if path.is_absolute() else project_path.parent / path
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"contrast file contains no rows: {path}")
    return rows


def initialize_figure_workflow(project_path: str | Path, *, force: bool = False) -> tuple[Path, ...]:
    """Scaffold story files and enable the generic publication constructor module."""
    path = Path(project_path).expanduser().resolve()
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or raw.get("version") != 2:
        raise ValueError("figures init requires a version 2 project.yaml")
    contrasts = _read_contrasts(path, raw)
    contrast_ids = [row["contrast_id"] for row in contrasts]
    claims_path = path.parent / "hypotheses.yaml"
    panels_path = path.parent / "hypothesis_panels.yaml"
    recipe_path = path.parent / "figure_recipe.yaml"
    if not force:
        existing = [candidate for candidate in (claims_path, panels_path, recipe_path) if candidate.exists()]
        if existing:
            raise FileExistsError(
                "refusing to overwrite existing publication files: "
                + ", ".join(str(candidate) for candidate in existing)
            )

    claims = {
        "study": raw.get("project", {}).get("title", raw.get("project", {}).get("id", "study")),
        "hypotheses": [
            {
                "id": f"{contrast_id}_working_hypothesis",
                "statement": "Replace with the biological hypothesis this contrast tests.",
                "contrast": contrast_id,
                "expected_direction": "context_dependent",
                "gene_panels": ["program_1"],
                "pathway_panels": ["pathways_1"],
            }
            for contrast_id in contrast_ids
        ],
    }
    panels = {
        "programs": {
            "program_1": {
                "label": "Biological program 1",
                "description": "Replace with a curated, hypothesis-relevant program.",
                "color": "#D97706",
                "genes": ["REPLACE_WITH_GENE_SYMBOLS"],
                "expected_direction": "not_specified",
            }
        },
        "pathway_panels": {
            "pathways_1": {
                "description": "Replace with pathways central to the hypothesis.",
                "pathways": [{"collection": "custom", "pathway": "REPLACE_WITH_PATHWAY_ID"}],
            }
        },
        "gsea_programs": ["program_1"],
        "program_order": ["Biological program 1"],
    }
    first_contrast = contrast_ids[0]
    recipe = {
        "figure_sets": {
            "primary": {
                "title": "Primary publication figure",
                "description": "Edit constructors, variants, and placement; run tifzoret figures gallery to review alternatives.",
                "width": 12,
                "height": 10,
                "units": "in",
                "columns": 2,
                "shared_legends": True,
                "panels": [
                    {"id": "A", "constructor": "pca_correlation", "row": 1, "column": 1, "column_span": 2},
                    {"id": "B", "constructor": "volcano", "contrast": first_contrast, "row": 2, "column": 1},
                    {"id": "C", "constructor": "de_heatmap", "variant": "global_clustered", "contrast": first_contrast, "row": 2, "column": 2},
                ],
            }
        }
    }
    _write_yaml(claims_path, claims, force)
    _write_yaml(panels_path, panels, force)
    _write_yaml(recipe_path, recipe, force)

    raw["hypotheses"] = {"claims": claims_path.name, "panels": panels_path.name}
    raw["publication"] = {"recipe": recipe_path.name}
    raw.setdefault("analysis", {}).setdefault("modules", {})["publication"] = True
    path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    return claims_path, panels_path, recipe_path
