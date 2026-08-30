"""Keep the navigation docs honest: they exist and cover every stage."""
from pathlib import Path

DOCS = Path(__file__).resolve().parents[1] / "docs"
SCRIPTS = Path(__file__).resolve().parents[1] / "src" / "tifzoret" / "workflow" / "scripts"


def test_nav_docs_exist():
    for name in ("engine-map.md", "walkthrough.md", "using-the-engine.md"):
        assert (DOCS / name).is_file(), f"missing docs/{name}"


def test_engine_map_lists_every_stage():
    text = (DOCS / "engine-map.md").read_text(encoding="utf-8")
    stems = [p.name for p in list(SCRIPTS.rglob("*.R")) + list(SCRIPTS.rglob("*.py"))
             if p.name != "utils.R" and "report_assets" not in p.parts and "__pycache__" not in p.parts]
    missing = [s for s in stems if s not in text]
    assert not missing, f"engine-map.md does not mention: {missing}"


def test_using_the_engine_covers_operator_questions():
    text = (DOCS / "using-the-engine.md").read_text(encoding="utf-8").lower()
    for cue in ("input", "tifzoret init", "profile", "skip", "output"):
        assert cue in text, f"using-the-engine.md missing coverage of: {cue}"
