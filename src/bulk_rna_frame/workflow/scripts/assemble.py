#!/usr/bin/env python3
"""Assemble recipe-selected vector panels and a raster review image."""

from __future__ import annotations

import argparse
import io
import json
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from pypdf import PageObject, PdfReader, PdfWriter, Transformation
from reportlab.pdfgen import canvas

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from bulk_rna_frame.config import load_project  # noqa: E402


def _panel_paths(results: Path, source: str) -> tuple[Path, Path]:
    path = Path(source).expanduser()
    path = path if path.is_absolute() else results / path
    if path.suffix.lower() == ".pdf":
        return path, path.with_suffix(".png")
    if path.suffix.lower() == ".png":
        return path.with_suffix(".pdf"), path
    return path.with_suffix(".pdf"), path.with_suffix(".png")


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


def assemble(project_config: Path, results: Path, figure_set: str, pdf: Path, png: Path, metadata: Path) -> None:
    project = load_project(project_config)
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
    review = Image.new("RGB", (round(width_in * 150), round(height_in * 150)), "white")
    review_draw = ImageDraw.Draw(review)
    font = ImageFont.load_default(size=24)
    placements = []
    for panel, (row, column) in zip(recipe["panels"], positions, strict=True):
        source_pdf, source_png = _panel_paths(results, panel["source"])
        if not source_pdf.is_file() or not source_png.is_file():
            raise FileNotFoundError(f"panel {panel['id']}: expected {source_pdf} and {source_png}")
        row_span = int(panel.get("row_span", 1))
        column_span = int(panel.get("column_span", 1))
        box_width, box_height = cell_width * column_span, cell_height * row_span
        source_page = PdfReader(str(source_pdf)).pages[0]
        source_width, source_height = float(source_page.mediabox.width), float(source_page.mediabox.height)
        scale = min(box_width / source_width, box_height / source_height) * 0.96
        x = (column - 1) * cell_width + (box_width - source_width * scale) / 2
        y_top = page_height - (row - 1) * cell_height
        y = y_top - box_height + (box_height - source_height * scale) / 2
        source_page.add_transformation(Transformation().scale(scale).translate(x, y))
        output_page.merge_page(source_page)
        output_page.merge_page(
            _label_overlay(
                page_width,
                page_height,
                str(panel["id"]),
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
        image.thumbnail((round(available[0] * 0.96), round(available[1] * 0.96)), Image.Resampling.LANCZOS)
        image_x = px_box[0] + (available[0] - image.width) // 2
        image_y = px_box[1] + (available[1] - image.height) // 2
        review.paste(image, (image_x, image_y))
        review_draw.text((px_box[0] + 5, px_box[1] + 4), str(panel["id"]), fill="black", font=font)
        placements.append({
            "id": panel["id"], "source_pdf": str(source_pdf), "source_png": str(source_png),
            "row": row, "column": column, "row_span": row_span, "column_span": column_span,
            "pdf_box_points": [x, y, source_width * scale, source_height * scale],
        })
    writer = PdfWriter()
    writer.add_page(output_page)
    pdf.parent.mkdir(parents=True, exist_ok=True)
    with pdf.open("wb") as handle:
        writer.write(handle)
    png.parent.mkdir(parents=True, exist_ok=True)
    review.save(png, dpi=(150, 150))
    metadata.parent.mkdir(parents=True, exist_ok=True)
    metadata.write_text(json.dumps({
        "schema_version": 1, "figure_set": figure_set,
        "dimensions": {"width": recipe["width"], "height": recipe["height"], "units": units},
        "panels": placements,
    }, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--figure-set", required=True)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--png", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    args = parser.parse_args()
    assemble(args.project_config.resolve(), args.results.resolve(), args.figure_set, args.pdf.resolve(), args.png.resolve(), args.metadata.resolve())


if __name__ == "__main__":
    main()
