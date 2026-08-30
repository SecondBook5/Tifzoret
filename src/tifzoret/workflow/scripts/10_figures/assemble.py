#!/usr/bin/env python3
# ┌─ TIFZORET STAGE ────────────────────────────────────────────────
# │ STAGE:     10_figures / assemble.py
# │ WHAT:      Multi-panel figure assembly (PDF + PNG raster preview)
# │ WHY:       Combines recipe-selected panel PDFs into publication figures with
# │            grid layout, labels, and SHA-256 provenance
# │ HOW:       pypdf grid assembly from figure_recipes.yaml; PIL raster preview
# │ INPUTS:    config (figure_recipes.yaml), panel PDFs from publication/modules
# │ PRODUCES:  figures/{recipe_id}.{pdf,png} + assembly_summary.json
# │ CALLED BY: rule assemble_figure (workflow/rules/publication.smk)
# │ ENV:       workflow/envs/core.yaml
# └─────────────────────────────────────────────────────────────────
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

SOURCE_ROOT = Path(__file__).resolve().parents[4]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from tifzoret.config import load_project  # noqa: E402
from tifzoret.figures import resolve_panel  # noqa: E402


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


# Outer page margin and inter-panel gutter for the assembled figure, in inches.
MARGIN_IN = 0.2
GUTTER_IN = 0.15


def _axis_sizes(panel_meta, count, pos_key, span_key, size_key):
    """Return the natural size of each grid track (column widths or row heights) in
    source points. A single-track panel forces its track to at least its own size; a
    panel spanning several tracks that still does not fit spreads the deficit evenly
    across the tracks it covers. Widest spans are resolved last so the result is
    stable regardless of panel order."""
    sizes = [0.0] * count
    for meta in panel_meta:
        if meta[span_key] == 1:
            index = meta[pos_key] - 1
            sizes[index] = max(sizes[index], meta[size_key])
    for meta in sorted((m for m in panel_meta if m[span_key] > 1), key=lambda m: m[span_key]):
        start, span = meta[pos_key] - 1, meta[span_key]
        covered = range(start, start + span)
        deficit = meta[size_key] - sum(sizes[i] for i in covered)
        if deficit > 0:
            for index in covered:
                sizes[index] += deficit / span
    return sizes


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
    """Lay out a figure set's recipe-selected panels onto one vector PDF page and a
    raster review PNG. Stage each source panel's PDF/PNG and its displayed-data
    artifacts with checksums, then write the assembly metadata and the panel
    index. Raise if the results root disagrees with the project or a panel source
    is missing."""
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
    # The recipe's declared size is a BOUNDING BOX, not a fixed frame: the grid is
    # scaled to the largest size that fits inside it, then the page is trimmed to the
    # packed content so no empty band survives when the content aspect differs from
    # the declared aspect. bound_* are that maximum extent, in points.
    bound_width = float(recipe["width"]) * unit_scale * 72
    bound_height = float(recipe["height"]) * unit_scale * 72
    rows, columns, positions = _positions(recipe)
    review_dpi = int(recipe.get("review_dpi", 150))
    background = recipe.get("background", "#FFFFFF")
    font = ImageFont.load_default(size=24)

    # Pre-pass: resolve every panel and read its native page size. Both the vector
    # page and its dimensions feed the shared-scale grid computed below.
    panel_meta = []
    for panel, (row, column) in zip(recipe["panels"], positions, strict=True):
        resolved = resolve_panel(project, panel)
        source_pdf = resolved.source.with_suffix(".pdf")
        source_png = resolved.source.with_suffix(".png")
        if not source_pdf.is_file() or not source_png.is_file():
            raise FileNotFoundError(f"panel {panel['id']}: expected {source_pdf} and {source_png}")
        source_page = PdfReader(str(source_pdf)).pages[0]
        panel_meta.append({
            "panel": panel, "resolved": resolved,
            "row": row, "column": column,
            "row_span": int(panel.get("row_span", 1)),
            "column_span": int(panel.get("column_span", 1)),
            "source_pdf": source_pdf, "source_png": source_png,
            "source_page": source_page,
            "sw": float(source_page.mediabox.width),
            "sh": float(source_page.mediabox.height),
        })

    # Content-proportional grid drawn at ONE shared scale. Each column takes the
    # width of its widest panel and each row the height of its tallest (deficits
    # from spanning panels spread across the tracks they cover); a single scale then
    # fits the whole grid inside the bounding box. Because every panel is placed at
    # that same scale, their absolute type sizes and line weights match and the sheet
    # reads as one figure — instead of each panel being shrunk independently into an
    # equal cell, which is what made dense and sparse panels clash.
    col_native = _axis_sizes(panel_meta, columns, "column", "column_span", "sw")
    row_native = _axis_sizes(panel_meta, rows, "row", "row_span", "sh")
    margin, gutter = MARGIN_IN * 72.0, GUTTER_IN * 72.0
    content_w = sum(col_native) or 1.0
    content_h = sum(row_native) or 1.0
    avail_w = max(bound_width - 2 * margin - (columns - 1) * gutter, 1.0)
    avail_h = max(bound_height - 2 * margin - (rows - 1) * gutter, 1.0)
    scale = min(avail_w / content_w, avail_h / content_h)
    # Trim the page to the packed grid: content at `scale`, plus interior gutters and
    # the outer margin. The declared bounding box caps the size; the actual page is
    # exactly the content, so the figure fills its frame in both dimensions.
    page_width = content_w * scale + (columns - 1) * gutter + 2 * margin
    page_height = content_h * scale + (rows - 1) * gutter + 2 * margin
    width_in, height_in = page_width / 72, page_height / 72
    output_page = PageObject.create_blank_page(width=page_width, height=page_height)
    review = Image.new("RGB", (round(width_in * review_dpi), round(height_in * review_dpi)), background)
    review_draw = ImageDraw.Draw(review)
    # Left edge (points from the left) of each column and top edge (points from the
    # bottom, PDF origin) of each row, including the outer margin and gutters.
    col_left = [margin]
    for index in range(1, columns):
        col_left.append(col_left[-1] + col_native[index - 1] * scale + gutter)
    row_top = [page_height - margin]
    for index in range(1, rows):
        row_top.append(row_top[-1] - row_native[index - 1] * scale - gutter)
    px_per_pt = review_dpi / 72.0

    placements = []
    panel_records = []
    for meta in panel_meta:
        panel, resolved = meta["panel"], meta["resolved"]
        row, column = meta["row"], meta["column"]
        row_span, column_span = meta["row_span"], meta["column_span"]
        source_pdf, source_png = meta["source_pdf"], meta["source_png"]
        source_page = meta["source_page"]
        source_width, source_height = meta["sw"], meta["sh"]
        constructor_defaults = (
            (project.panel_config or {}).get("constructor_defaults", {}).get(
                resolved.constructor, {}
            )
            if resolved.constructor
            else {}
        )
        panel_options = {**constructor_defaults, **panel.get("options", {})}

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
            # A panel with no declared displayed-data source cannot be covered by
            # the Source Data bundle, so it is a publication-blocking degradation
            # (not a standing caveat): strict mode must fail rather than ship it.
            "warnings": [] if displayed_records else [
                {"message": "No displayed-data path was declared for this legacy source panel.",
                 "severity": "degradation"}
            ],
        }
        (staged_dir / "panel.json").write_text(
            json.dumps({"schema_version": 1, **panel_record}, indent=2) + "\n",
            encoding="utf-8",
        )
        panel_records.append(panel_record)

        # Optional per-panel fine-tune multiplier. Default 1.0 keeps the shared
        # scale; a recipe may shrink one panel below the grid without disturbing
        # the rest, but panels are never enlarged past the common scale.
        fit_fraction = float(panel_options.get("scale", 1.0))
        if not 0.1 <= fit_fraction <= 1:
            raise ValueError(f"panel {panel['id']}: options.scale must be between 0.1 and 1")
        panel_scale = scale * fit_fraction

        # The panel's block spans col..col+column_span-1 and row..row+row_span-1,
        # covering any interior gutters. Panels are centered within their block so
        # narrower/shorter panels sit balanced against the shared gridlines.
        last_col = column - 1 + column_span - 1
        last_row = row - 1 + row_span - 1
        block_left = col_left[column - 1]
        block_right = col_left[last_col] + col_native[last_col] * scale
        block_top = row_top[row - 1]
        block_bottom = row_top[last_row] - row_native[last_row] * scale
        block_width = block_right - block_left
        block_height = block_top - block_bottom
        panel_w = source_width * panel_scale
        panel_h = source_height * panel_scale
        x = block_left + (block_width - panel_w) / 2
        y = block_bottom + (block_height - panel_h) / 2
        source_page.add_transformation(Transformation().scale(panel_scale).translate(x, y))
        output_page.merge_page(source_page)
        if panel_options.get("show_panel_label", True):
            output_page.merge_page(
                _label_overlay(
                    page_width,
                    page_height,
                    str(panel_options.get("panel_label", panel["id"])),
                    block_left,
                    block_top,
                )
            )

        # Mirror the exact same geometry into the raster review image (origin top-left).
        image = Image.open(source_png).convert("RGB")
        target = (max(1, round(panel_w * px_per_pt)), max(1, round(panel_h * px_per_pt)))
        image = image.resize(target, Image.Resampling.LANCZOS)
        block_left_px = round(block_left * px_per_pt)
        block_top_px = round((page_height - block_top) * px_per_pt)
        block_w_px = round(block_width * px_per_pt)
        block_h_px = round(block_height * px_per_pt)
        image_x = block_left_px + (block_w_px - image.width) // 2
        image_y = block_top_px + (block_h_px - image.height) // 2
        review.paste(image, (image_x, image_y))
        if panel_options.get("show_panel_label", True):
            review_draw.text(
                (block_left_px + 5, block_top_px + 4),
                str(panel_options.get("panel_label", panel["id"])),
                fill="black",
                font=font,
            )
        placements.append({
            "id": panel["id"], "constructor": resolved.constructor, "variant": resolved.variant,
            "contrast": resolved.contrast, "source_pdf": str(source_pdf), "source_png": str(source_png),
            "row": row, "column": column, "row_span": row_span, "column_span": column_span,
            "pdf_box_points": [x, y, panel_w, panel_h],
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
        "dimensions": {
            "width": round(page_width / 72 / unit_scale, 3),
            "height": round(page_height / 72 / unit_scale, 3),
            "units": units,
            "declared_width": recipe["width"],
            "declared_height": recipe["height"],
        },
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
    """Parse arguments and assemble the requested figure set into its PDF, review PNG, metadata, and panel index."""
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
