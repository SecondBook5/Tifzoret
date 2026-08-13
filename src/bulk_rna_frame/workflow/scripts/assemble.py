#!/usr/bin/env python3
"""Assemble recipe-selected vector panels and a raster review image."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import shutil
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from pypdf import PageObject, PdfReader, PdfWriter, Transformation
from reportlab.pdfgen import canvas

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from bulk_rna_frame.config import load_project  # noqa: E402
from bulk_rna_frame.figures import resolve_panel  # noqa: E402


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _positions(recipe: dict[str, object]) -> tuple[int, int, list[tuple[int, int]]]:
    panels = recipe["panels"]
    columns = int(recipe.get("columns", 1))
    positions: list[tuple[int, int]] = []
    for index, panel in enumerate(panels):
        row = int(panel.get("row", index // columns + 1))
        column = int(panel.get("column", index % columns + 1))
        positions.append((row, column))
    rows = max(row + int(panel.get("row_span", 1)) - 1 for (row, _), panel in zip(positions, panels, strict=True))
    columns = max(columns, max(column + int(panel.get("column_span", 1)) - 1 for (_, column), panel in zip(positions, panels, strict=True)))
    return rows, columns, positions


def _label_overlay(
    width: float, height: float, label: str, x: float, y_top: float
) -> PageObject:
    buffer = io.BytesIO()
    drawing = canvas.Canvas(buffer, pagesize=(width, height))
    drawing.setFont("Helvetica-Bold", 14)
    drawing.drawString(x + 5, y_top - 17, label)
    drawing.save()
    buffer.seek(0)
    return PdfReader(buffer).pages[0]


def assemble(
    project_config: Path,
    results: Path,
    figure_set: str,
    pdf: Path,
    png: Path,
    metadata: Path,
    panel_index: Path,
) -> None:
    project = load_project(project_config)
    if results.resolve() != project.result_root.resolve():
        raise ValueError(
            f"assembly results root {results} does not match resolved project root {project.result_root}"
        )
    if project.recipe_config is None:
        raise RuntimeError("project has no validated figure recipe")
    recipe = project.recipe_config["figure_sets"][figure_set]
    units = recipe.get("units", "in")
    unit_scale = {"in": 1.0, "mm": 1 / 25.4, "cm": 1 / 2.54}[units]
    width_in = float(recipe["width"]) * unit_scale
    height_in = float(recipe["height"]) * unit_scale
    page_width, page_height = width_in * 72, height_in * 72
    rows, columns, positions = _positions(recipe)
    cell_width, cell_height = page_width / columns, page_height / rows
    output_page = PageObject.create_blank_page(width=page_width, height=page_height)
    review_dpi = int(recipe.get("review_dpi", 150))
    background = recipe.get("background", "#FFFFFF")
    review = Image.new("RGB", (round(width_in * review_dpi), round(height_in * review_dpi)), background)
    review_draw = ImageDraw.Draw(review)
    font = ImageFont.load_default(size=24)
    placements = []
    panel_records = []
    for panel, (row, column) in zip(recipe["panels"], positions, strict=True):
        resolved = resolve_panel(project, panel)
        constructor_defaults = (
            (project.panel_config or {}).get("constructor_defaults", {}).get(
                resolved.constructor, {}
            )
            if resolved.constructor
            else {}
        )
        panel_options = {**constructor_defaults, **panel.get("options", {})}
        source_pdf, source_png = resolved.source.with_suffix(".pdf"), resolved.source.with_suffix(".png")
        if not source_pdf.is_file() or not source_png.is_file():
            raise FileNotFoundError(f"panel {panel['id']}: expected {source_pdf} and {source_png}")

        staged_dir = panel_index.parent / str(panel["id"])
        staged_dir.mkdir(parents=True, exist_ok=True)
        staged_pdf, staged_png = staged_dir / "panel.pdf", staged_dir / "panel.png"
        shutil.copy2(source_pdf, staged_pdf)
        shutil.copy2(source_png, staged_png)
        displayed_records = []
        displayed_dir = staged_dir / "displayed_data"
        for data_index, data_path in enumerate(resolved.displayed_data, start=1):
            if not data_path.is_file():
                raise FileNotFoundError(
                    f"panel {panel['id']}: registered displayed-data artifact is missing: {data_path}"
                )
            displayed_dir.mkdir(parents=True, exist_ok=True)
            staged_data = displayed_dir / f"{data_index:02d}_{data_path.name}"
            shutil.copy2(data_path, staged_data)
            displayed_records.append({
                "source": str(data_path),
                "staged": str(staged_data),
                "sha256": _sha256(data_path),
            })
        panel_record = {
            "id": panel["id"],
            "constructor": resolved.constructor,
            "variant": resolved.variant,
            "contrast": resolved.contrast,
            "title": resolved.label,
            "caption": panel.get("caption"),
            "options": panel_options,
            "show_legend": panel.get("show_legend", True),
            "source_pdf": str(source_pdf),
            "source_png": str(source_png),
            "staged_pdf": str(staged_pdf),
            "staged_png": str(staged_png),
            "pdf_sha256": _sha256(source_pdf),
            "png_sha256": _sha256(source_png),
            "displayed_data": displayed_records,
            "warnings": [] if displayed_records else ["No displayed-data path was declared for this legacy source panel."],
        }
        (staged_dir / "panel.json").write_text(
            json.dumps({"schema_version": 1, **panel_record}, indent=2) + "\n",
            encoding="utf-8",
        )
        panel_records.append(panel_record)
        row_span = int(panel.get("row_span", 1))
        column_span = int(panel.get("column_span", 1))
        box_width, box_height = cell_width * column_span, cell_height * row_span
        source_page = PdfReader(str(source_pdf)).pages[0]
        source_width, source_height = float(source_page.mediabox.width), float(source_page.mediabox.height)
        fit_fraction = float(panel_options.get("scale", 0.96))
        if not 0.1 <= fit_fraction <= 1:
            raise ValueError(f"panel {panel['id']}: options.scale must be between 0.1 and 1")
        scale = min(box_width / source_width, box_height / source_height) * fit_fraction
        x = (column - 1) * cell_width + (box_width - source_width * scale) / 2
        y_top = page_height - (row - 1) * cell_height
        y = y_top - box_height + (box_height - source_height * scale) / 2
        source_page.add_transformation(Transformation().scale(scale).translate(x, y))
        output_page.merge_page(source_page)
        if panel_options.get("show_panel_label", True):
            output_page.merge_page(
                _label_overlay(
                    page_width,
                    page_height,
                    str(panel_options.get("panel_label", panel["id"])),
                    (column - 1) * cell_width,
                    y_top,
                )
            )

        image = Image.open(source_png).convert("RGB")
        px_box = (
            round((column - 1) * review.width / columns),
            round((row - 1) * review.height / rows),
            round((column - 1 + column_span) * review.width / columns),
            round((row - 1 + row_span) * review.height / rows),
        )
        available = (px_box[2] - px_box[0], px_box[3] - px_box[1])
        image.thumbnail(
            (round(available[0] * fit_fraction), round(available[1] * fit_fraction)),
            Image.Resampling.LANCZOS,
        )
        image_x = px_box[0] + (available[0] - image.width) // 2
        image_y = px_box[1] + (available[1] - image.height) // 2
        review.paste(image, (image_x, image_y))
        if panel_options.get("show_panel_label", True):
            review_draw.text(
                (px_box[0] + 5, px_box[1] + 4),
                str(panel_options.get("panel_label", panel["id"])),
                fill="black",
                font=font,
            )
        placements.append({
            "id": panel["id"], "constructor": resolved.constructor, "variant": resolved.variant,
            "contrast": resolved.contrast, "source_pdf": str(source_pdf), "source_png": str(source_png),
            "row": row, "column": column, "row_span": row_span, "column_span": column_span,
            "pdf_box_points": [x, y, source_width * scale, source_height * scale],
        })
    writer = PdfWriter()
    writer.add_page(output_page)
    pdf.parent.mkdir(parents=True, exist_ok=True)
    with pdf.open("wb") as handle:
        writer.write(handle)
    png.parent.mkdir(parents=True, exist_ok=True)
    review.save(png, dpi=(review_dpi, review_dpi))
    metadata.parent.mkdir(parents=True, exist_ok=True)
    metadata.write_text(json.dumps({
        "schema_version": 1, "figure_set": figure_set,
        "title": recipe.get("title"), "description": recipe.get("description"),
        "dimensions": {"width": recipe["width"], "height": recipe["height"], "units": units},
        "shared_legends": recipe.get("shared_legends", False), "review_dpi": review_dpi,
        "panels": placements,
    }, indent=2) + "\n", encoding="utf-8")
    panel_index.parent.mkdir(parents=True, exist_ok=True)
    panel_index.write_text(json.dumps({
        "schema_version": 1,
        "figure_set": figure_set,
        "panels": panel_records,
    }, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--figure-set", required=True)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--png", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--panel-index", type=Path, required=True)
    args = parser.parse_args()
    assemble(
        args.project_config.resolve(), args.results.resolve(), args.figure_set,
        args.pdf.resolve(), args.png.resolve(), args.metadata.resolve(), args.panel_index.resolve(),
    )


if __name__ == "__main__":
    main()
