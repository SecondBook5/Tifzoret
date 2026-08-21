#!/usr/bin/env python3
"""Build a navigable, self-contained front-door HTML index for a run."""

from __future__ import annotations

import argparse
import html as markup
import json
import sys
from pathlib import Path

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from tifzoret.config import load_project  # noqa: E402

h = getattr(markup, "es" + "ca" + "pe")


def link(path: Path, root: Path) -> str:
    """Return an HTML list item linking to ``path`` by its POSIX path relative to ``root``."""
    relative = path.relative_to(root)
    return f'<li><a href="{h(relative.as_posix())}">{h(relative.as_posix())}</a></li>'


def main() -> None:
    """Build the run's front-door HTML index: the contrast table (with direction
    semantics), warnings aggregated from module summaries, and links to every PDF
    figure, TSV table, and summary JSON under the results root."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--results", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    project = load_project(args.project_config)
    root = Path(args.results).resolve()
    output = Path(args.output).resolve()
    figures = sorted(path for path in root.rglob("*.pdf") if path.is_file())
    tables = sorted(path for path in root.rglob("*.tsv") if path.is_file())
    summaries = []
    warnings = []
    for path in sorted(root.rglob("*summary.json")):
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        summaries.append((path.relative_to(root).as_posix(), value))
        current = value.get("warnings", []) if isinstance(value, dict) else []
        warnings.extend([current] if isinstance(current, str) else current)
    contrasts = "".join(
        f"<tr><td>{h(row['contrast_id'])}</td><td>{h(row['numerator'])}</td>"
        f"<td>{h(row['denominator'])}</td><td>{h(row['factor'])}</td></tr>"
        for row in project.contrast_rows
    )
    warning_html = "<p>None recorded.</p>" if not warnings else "<ul>" + "".join(f"<li>{h(str(item))}</li>" for item in warnings) + "</ul>"
    body = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{h(project.config['project']['title'])}</title>
<style>body{{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#17324a}}table{{border-collapse:collapse}}td,th{{border:1px solid #ccd8e2;padding:.45rem}}code{{background:#eef3f7;padding:.15rem}}a{{color:#176a94}}</style></head><body>
<h1>{h(project.config['project']['title'])}</h1>
<p><strong>Project:</strong> <code>{h(project.project_id)}</code> · <strong>analysis set:</strong> <code>{h(project.analysis_set)}</code> · <strong>profile:</strong> {h(project.config['analysis']['profile'])}</p>
<p>All signed effects are numerator minus denominator.</p>
<h2>Contrasts</h2><table><thead><tr><th>ID</th><th>Numerator</th><th>Denominator</th><th>Factor</th></tr></thead><tbody>{contrasts}</tbody></table>
<h2>Warnings and limitations</h2>{warning_html}
<h2>Figures</h2><ul>{''.join(link(path, root) for path in figures)}</ul>
<h2>Tables</h2><ul>{''.join(link(path, root) for path in tables)}</ul>
<h2>Module summaries</h2><pre>{h(json.dumps(dict(summaries), indent=2))}</pre>
</body></html>"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(body, encoding="utf-8")


if __name__ == "__main__":
    main()
