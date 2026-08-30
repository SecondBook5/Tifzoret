"""Guards for the self-contained HTML report builder (scripts/report.py).

These pin the properties a reviewer relies on when they open
the shared report offline and when they rerun the pipeline: no network dependency,
no HTML injection, no template-token corruption, and no leak of the author's
absolute filesystem layout into the shareable deliverable.
"""

from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "10_figures" / "report.py"


def _load_module():
    """Import the report script as a module (it adds src/ to sys.path itself)."""
    spec = importlib.util.spec_from_file_location("report_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


report = _load_module()


def _minimal_payload(**meta):
    base_meta = {
        "title": "Study title",
        "project_id": "proj",
        "analysis_set": "set",
        "profile": "standard",
        "contrast_semantics": "numerator minus denominator",
        "engine": {"name": "Tifzoret", "version": None},
        "counts": {"figures": 0, "assembled": 0, "panels": 0, "rows": 0, "contrasts": 0},
    }
    base_meta.update(meta)
    return {
        "meta": base_meta,
        "contrasts": [],
        "figures": {"assembled": [], "panels": []},
        "analyses": [],
        "provenance": {},
    }


def _embedded_json(html_doc: str) -> dict:
    """Extract and parse the JSON the front-end reads from the rendered report."""
    match = re.search(
        r'<script id="tifzoret-data" type="application/json">(.*?)</script>',
        html_doc,
        flags=re.DOTALL,
    )
    assert match, "payload <script> block not found"
    return json.loads(match.group(1))


def test_report_is_self_contained_and_offline():
    """No network dependency: nothing is fetched at view time, and fonts are inlined.

    The only permitted ``http`` references are XML/SVG *namespace* URIs (identifiers
    the browser never requests), so we forbid fetchable resource references instead of
    the literal substring ``http``.
    """
    doc = report.render(_minimal_payload(), "Study title")
    for needle in ("googleapis", "unpkg", "jsdelivr", "cdn.", "@import", "//fonts."):
        assert needle not in doc, f"external dependency leaked: {needle!r}"
    for pattern in (r'src="https?://', r'href="https?://', r"url\(\s*['\"]?https?://"):
        assert not re.search(pattern, doc), f"external resource load leaked: {pattern!r}"
    # Every http(s) reference that survives must be a non-fetched namespace URI.
    for url in re.findall(r"https?://[^\s\"')]+", doc):
        assert "www.w3.org" in url, f"unexpected external URL: {url!r}"
    assert 'id="tifzoret-data"' in doc
    # Fonts vendored inline keep it offline; ship at least one @font-face data URI.
    assert "data:font/woff2;base64," in doc


def test_report_escapes_html_in_title():
    """A title with angle brackets/quotes is escaped, never injected as raw markup."""
    doc = report.render(_minimal_payload(), '</title><script>alert(1)</script>')
    assert "<script>alert(1)</script>" not in doc
    assert "&lt;script&gt;alert(1)" in doc


def test_render_survives_template_tokens_inside_payload():
    """A result value containing the literal substitution tokens must not corrupt the
    file: the payload stays valid JSON and the JS asset is inlined exactly once."""
    hostile = "__JS__ __CSS__ __PAYLOAD__ __TITLE__"
    payload = _minimal_payload()
    payload["analyses"] = [{"key": "de", "label": "DE", "note": hostile}]
    doc = report.render(payload, "Study title")
    parsed = _embedded_json(doc)
    # The hostile value round-trips verbatim (no asset injected into the JSON).
    assert parsed["analyses"][0]["note"] == hostile
    # app.js is present and not duplicated by an accidental second substitution.
    app_js = (report.ASSETS / "app.js").read_text(encoding="utf-8")
    signature = 'JSON.parse(document.getElementById("tifzoret-data").textContent)'
    assert app_js.count(signature) == 1
    assert doc.count(signature) == 1


def test_render_is_deterministic():
    """Identical payload in -> byte-identical document out."""
    payload = _minimal_payload()
    assert report.render(payload, "Study title") == report.render(payload, "Study title")


def test_provenance_scrubs_absolute_local_paths(tmp_path):
    """The shareable report must not embed the author's absolute filesystem layout."""
    secret = "/home/secretuser/miniconda3/envs/study-env"
    image = "/home/secretuser/images/study.sif"
    manifest = {
        "generated_utc": "2026-08-24T00:00:00+00:00",
        "environment": {"conda_prefix": secret, "container": image},
        "repository": {"revision": "abc123", "dirty": False},
        "inputs": [],
        "results": [],
    }
    (tmp_path / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    prov = report._provenance(tmp_path)
    assert prov["environment"]["conda_prefix"] == "study-env"
    assert prov["environment"]["container"] == "study.sif"
    assert "/home/secretuser" not in json.dumps(prov)


def test_provenance_scrub_preserves_none():
    """A missing conda_prefix/container stays None rather than becoming a bogus name."""
    assert report._scrub_environment({"conda_prefix": None, "container": None}) == {
        "conda_prefix": None,
        "container": None,
    }
    assert report._scrub_environment(None) is None


def test_provenance_absent_manifest_is_empty(tmp_path):
    """No manifest -> empty provenance, never a fabricated block."""
    assert report._provenance(tmp_path) == {}
