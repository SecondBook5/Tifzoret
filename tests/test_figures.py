from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml
from PIL import Image
from reportlab.pdfgen import canvas

from tifzoret.cli import main
from tifzoret.config import ProjectValidationError, load_project
from tifzoret.figures import PANEL_REGISTRY, build_gallery, resolve_panel


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"


def _project(tmp_path: Path) -> Path:
    destination = tmp_path / "study"
    shutil.copytree(TEMPLATE, destination)
    return destination / "project.yaml"


def _figure_pair(stem: Path, label: str) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    drawing = canvas.Canvas(str(stem.with_suffix(".pdf")), pagesize=(360, 240))
    drawing.drawString(30, 120, label)
    drawing.save()
    Image.new("RGB", (600, 400), "white").save(stem.with_suffix(".png"))


def test_figures_init_scaffolds_programs_and_constructor_recipe(tmp_path):
    project_path = _project(tmp_path)
    assert main(["figures", "init", str(project_path)]) == 0
    project = load_project(project_path)
    assert "publication" in project.modules
    assert project.panel_config["programs"]["program_1"]["genes"] == [
        "REPLACE_WITH_GENE_SYMBOLS"
    ]
    panel = project.recipe_config["figure_sets"]["primary"]["panels"][2]
    assert panel["constructor"] == "de_heatmap"
    assert panel["variant"] == "global_clustered"
    before = (project_path.parent / "hypothesis_panels.yaml").read_text()
    assert main(["figures", "init", str(project_path)]) == 2
    assert (project_path.parent / "hypothesis_panels.yaml").read_text() == before


def test_constructor_catalog_is_public_and_contains_derived_variants(capsys):
    assert main(["figures", "catalog", "--json"]) == 0
    catalog = json.loads(capsys.readouterr().out)["constructors"]
    indexed = {item["id"]: item for item in catalog}
    assert set(indexed["de_heatmap"]["variants"]) == {
        "default",
        "global_clustered",
        "program_grouped",
        "direct_program_labels",
    }
    assert set(indexed["dorothea_grn"]["variants"]) == {"rectangular", "radial", "radial_legacy"}
    assert "program_violins" in PANEL_REGISTRY


def test_constructor_recipe_rejects_unknown_variant_and_missing_contrast(tmp_path):
    project_path = _project(tmp_path)
    assert main(["figures", "init", str(project_path)]) == 0
    recipe_path = project_path.parent / "figure_recipe.yaml"
    recipe = yaml.safe_load(recipe_path.read_text())
    recipe["figure_sets"]["primary"]["panels"][2]["variant"] = "invented"
    recipe_path.write_text(yaml.safe_dump(recipe, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="has no variant"):
        load_project(project_path)

    recipe["figure_sets"]["primary"]["panels"][2]["variant"] = "global_clustered"
    del recipe["figure_sets"]["primary"]["panels"][1]["contrast"]
    recipe_path.write_text(yaml.safe_dump(recipe, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="contrast is required"):
        load_project(project_path)


def test_gallery_and_assembly_stage_auditable_panel_artifacts(tmp_path):
    project_path = _project(tmp_path)
    assert main(["figures", "init", str(project_path)]) == 0
    project = load_project(project_path)
    recipe = project.recipe_config["figure_sets"]["primary"]
    for panel in recipe["panels"]:
        resolved = resolve_panel(project, panel)
        _figure_pair(resolved.source, f"{resolved.constructor}/{resolved.variant}")
        for displayed in resolved.displayed_data:
            displayed.parent.mkdir(parents=True, exist_ok=True)
            if displayed.suffix == ".json":
                displayed.write_text('{"test": true}\n')
            else:
                displayed.write_text("feature\tvalue\nexample\t1\n")

    gallery = build_gallery(project)
    assert gallery.is_file()
    gallery_data = json.loads((gallery.parent / "gallery.json").read_text())
    assert any(panel["selected"] and panel["available"] for panel in gallery_data["panels"])
    assert (gallery.parent / "contact_sheet.png").is_file()

    assembled = project.result_root / "publication" / "primary"
    command = [
        sys.executable,
        str(ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "assemble.py"),
        "--project-config", str(project_path),
        "--results", str(project.result_root),
        "--figure-set", "primary",
        "--pdf", str(assembled / "assembled" / "primary.pdf"),
        "--png", str(assembled / "assembled" / "primary.png"),
        "--metadata", str(assembled / "assembled" / "assembly.json"),
        "--panel-index", str(assembled / "panels" / "index.json"),
    ]
    subprocess.run(command, check=True)
    index = json.loads((assembled / "panels" / "index.json").read_text())
    assert len(index["panels"]) == 3
    assert index["panels"][2]["constructor"] == "de_heatmap"
    assert index["panels"][2]["displayed_data"]
    for panel in recipe["panels"]:
        panel_dir = assembled / "panels" / panel["id"]
        assert (panel_dir / "panel.pdf").is_file()
        assert (panel_dir / "panel.png").is_file()
        assert (panel_dir / "panel.json").is_file()
