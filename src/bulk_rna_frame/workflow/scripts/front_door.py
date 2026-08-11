#!/usr/bin/env python3
"""Promote review-facing artifacts and describe them with stable metadata."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

from PIL import Image
from pypdf import PdfReader

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from bulk_rna_frame.config import load_project  # noqa: E402


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_name(path: Path, root: Path) -> str:
    relative = path.relative_to(root).with_suffix("")
    return "__".join(relative.parts)


def dimensions(path: Path) -> dict[str, object]:
    if path.suffix.lower() == ".png":
        with Image.open(path) as image:
            return {"width_px": image.width, "height_px": image.height}
    reader = PdfReader(str(path))
    page = reader.pages[0]
    return {
        "pages": len(reader.pages),
        "width_points": float(page.mediabox.width),
        "height_points": float(page.mediabox.height),
    }


def resolve_panel_sources(project, root: Path) -> list[tuple[str, Path, dict]]:
    selected: list[tuple[str, Path, dict]] = []
    recipe = project.recipe_config or {}
    for figure_set, spec in recipe.get("figure_sets", {}).items():
        for panel in spec.get("panels", []):
            source = Path(panel["source"]).expanduser()
            source = source if source.is_absolute() else root / source
            stem = source.with_suffix("") if source.suffix.lower() in {".pdf", ".png"} else source
            selected.append((f"{figure_set}__panel_{panel['id']}", stem, panel))
        assembled = root / "publication" / figure_set / "assembled" / figure_set
        selected.append((f"{figure_set}__assembled", assembled, {"figure_set": figure_set}))
    if selected:
        return selected

    defaults = [
        ("qc__pca_correlation", root / "qc" / "figures" / "pca_correlation", {}),
    ]
    for row in project.contrast_rows:
        contrast = row["contrast_id"]
        base = root / "contrasts" / contrast / "analyses"
        for module, stem in (
            ("de", "de_overview"),
            ("pathways", "ora_bidirectional"),
            ("pathways", "gsva_heatmap"),
            ("pathways", "gsea_curves"),
        ):
            defaults.append((f"{contrast}__{stem}", base / module / "figures" / stem, {}))
    return defaults


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--results", required=True)
    parser.add_argument("--figures-index", required=True)
    parser.add_argument("--tables-index", required=True)
    args = parser.parse_args()

    project = load_project(args.project_config)
    root = Path(args.results).resolve()
    figures_dir = Path(args.figures_index).resolve().parent
    tables_dir = Path(args.tables_index).resolve().parent
    metadata_dir = figures_dir / "metadata"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    figure_records = []
    seen = set()
    for label, stem, selection in resolve_panel_sources(project, root):
        if label in seen:
            continue
        available = [stem.with_suffix(ext) for ext in (".pdf", ".png")]
        if not all(path.is_file() for path in available):
            continue
        seen.add(label)
        artifacts = []
        for source in available:
            destination = figures_dir / f"{label}{source.suffix.lower()}"
            shutil.copy2(source, destination)
            artifacts.append({
                "path": destination.relative_to(root).as_posix(),
                "source": source.relative_to(root).as_posix(),
                "sha256": checksum(destination),
                **dimensions(destination),
            })
        record = {
            "id": label,
            "selection": selection,
            "artifacts": artifacts,
        }
        metadata_path = metadata_dir / f"{label}.json"
        metadata_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        record["metadata"] = metadata_path.relative_to(root).as_posix()
        figure_records.append(record)

    table_records = []
    for source in sorted(root.rglob("*.tsv")):
        if tables_dir in source.parents or not source.is_file():
            continue
        if "displayed" not in source.stem and source.name not in {
            "de_results.tsv", "fgsea.tsv", "ora.tsv", "gsva_differential.tsv",
            "hypothesis_summary.tsv", "regulator_differential.tsv",
        }:
            continue
        name = f"{safe_name(source, root)}.tsv"
        destination = tables_dir / name
        shutil.copy2(source, destination)
        table_records.append({
            "path": destination.relative_to(root).as_posix(),
            "source": source.relative_to(root).as_posix(),
            "sha256": checksum(destination),
        })

    figures_index = {
        "schema_version": 1,
        "project_id": project.project_id,
        "analysis_set": project.analysis_set,
        "figures": figure_records,
    }
    tables_index = {
        "schema_version": 1,
        "project_id": project.project_id,
        "analysis_set": project.analysis_set,
        "tables": table_records,
    }
    Path(args.figures_index).write_text(json.dumps(figures_index, indent=2) + "\n", encoding="utf-8")
    Path(args.tables_index).write_text(json.dumps(tables_index, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
