"""Every pipeline stage script must carry the fixed self-describing header.

Fields make each file self-locating: its phase, what it does, WHY (the science),
how, its data in/out, its rule, and its conda env — without reading the body.
Uses rglob so it is agnostic to whether scripts are flat or already foldered.
"""
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "src" / "tifzoret" / "workflow" / "scripts"
REQUIRED = ("STAGE:", "WHAT:", "WHY:", "HOW:", "INPUTS:", "PRODUCES:", "CALLED BY:", "ENV:")
# utils.R, estimands.R, and de_render.R are shared helpers, not stages; report_assets holds no scripts.
EXCLUDE = {"utils.R", "estimands.R", "de_render.R"}


def _stage_files():
    for path in sorted(SCRIPTS.rglob("*.R")) + sorted(SCRIPTS.rglob("*.py")):
        if path.name in EXCLUDE or "report_assets" in path.parts or "__pycache__" in path.parts:
            continue
        yield path


def test_every_stage_has_all_header_fields():
    missing = {}
    for path in _stage_files():
        head = "\n".join(path.read_text(encoding="utf-8").splitlines()[:15])
        absent = [field for field in REQUIRED if field not in head]
        if absent:
            missing[path.name] = absent
    assert not missing, f"stage scripts missing header fields: {missing}"


def test_found_expected_stage_count():
    # Guards against the glob silently matching nothing (e.g. wrong path after a move).
    assert len(list(_stage_files())) == 36
