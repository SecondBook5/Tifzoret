"""HTML and PNG gallery generation."""

from __future__ import annotations

import hashlib
import html
import json
import shutil
from pathlib import Path
from typing import TYPE_CHECKING

from PIL import Image, ImageDraw, ImageFont

from .resolve import iter_review_panels

if TYPE_CHECKING:
    from ..config import ResolvedProject


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_gallery(project: "ResolvedProject", output: str | Path | None = None) -> Path:
    """Create an HTML/PNG review gallery of all available registered variants."""
    outdir = Path(output).expanduser().resolve() if output else project.result_root / "publication" / "gallery"
    assets = outdir / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    records = list(iter_review_panels(project))
    available = [record for record in records if record["available"]]
    for index, record in enumerate(available, start=1):
        name = f"{index:03d}_{record['constructor']}_{record['variant']}_{record['contrast'] or 'study'}.png"
        destination = assets / name
        shutil.copy2(record["source_png"], destination)
        record["review_image"] = f"assets/{name}"
        record["sha256"] = _sha256(Path(record["source_png"]))

    width, tile_width, tile_height, columns = 1600, 380, 300, 4
    rows = max(1, (len(available) + columns - 1) // columns)
    sheet = Image.new("RGB", (width, rows * tile_height), "white")
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.load_default(size=16)
    for index, record in enumerate(available):
        row, column = divmod(index, columns)
        x, y = column * tile_width, row * tile_height
        image = Image.open(outdir / record["review_image"]).convert("RGB")
        image.thumbnail((tile_width - 20, tile_height - 55), Image.Resampling.LANCZOS)
        sheet.paste(image, (x + (tile_width - image.width) // 2, y + 35))
        marker = "SELECTED · " if record["selected"] else ""
        label = f"{marker}{record['constructor']} / {record['variant']}"
        draw.text((x + 8, y + 8), label, fill="#14324A", font=font)
    sheet_path = outdir / "contact_sheet.png"
    sheet.save(sheet_path, dpi=(150, 150))

    cards = []
    for record in records:
        status = "selected" if record["selected"] else "available" if record["available"] else "not built"
        image = f'<img src="{html.escape(record.get("review_image", ""))}" alt="review image">' if record["available"] else '<div class="missing">Not built</div>'
        displayed = "".join(f"<li>{html.escape(value)}</li>" for value in record["displayed_data"])
        cards.append(
            f'<article class="card {status.replace(" ", "-")}">{image}'
            f'<h2>{html.escape(record["constructor_label"])} · {html.escape(record["variant_label"])}</h2>'
            f'<p><strong>{html.escape(status.upper())}</strong> · contrast: {html.escape(record["contrast"] or "study-level")}</p>'
            f'<p>{html.escape(record["description"])}</p><details><summary>Displayed data</summary><ul>{displayed}</ul></details></article>'
        )
    page = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Tifzoret figure gallery</title>
<style>body{{font:15px system-ui;margin:2rem;color:#14324A;background:#F7F9FA}}header{{max-width:70rem}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:1rem}}.card{{background:white;border:1px solid #DCE4E8;border-radius:10px;padding:1rem}}.card.selected{{border:3px solid #0F9D78}}img{{width:100%;height:240px;object-fit:contain;background:white}}.missing{{height:240px;display:grid;place-items:center;background:#EEF2F4;color:#6B7C87}}h2{{font-size:1.05rem}}details{{font-size:.8rem;overflow-wrap:anywhere}}</style></head>
<body><header><h1>Tifzoret publication figure gallery</h1><p>{html.escape(project.project_id)} · generated from registered constructors. Statistical results are unchanged; this page reviews presentation variants.</p><p><a href="contact_sheet.png">Open contact sheet</a></p></header><main class="grid">{''.join(cards)}</main></body></html>"""
    index = outdir / "index.html"
    index.write_text(page, encoding="utf-8")
    (outdir / "gallery.json").write_text(json.dumps({"schema_version": 1, "project": project.project_id, "panels": records}, indent=2) + "\n", encoding="utf-8")
    return index
