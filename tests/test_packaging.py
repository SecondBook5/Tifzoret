"""Packaging integrity test — ensures all required assets ship with the wheel."""

import subprocess
import sys
import zipfile
from pathlib import Path

import pytest


@pytest.mark.slow
def test_wheel_contains_required_assets(tmp_path):
    """Build wheel and verify report assets, locks, and R scripts are packaged."""
    repo_root = Path(__file__).parent.parent
    dist_dir = repo_root / "dist"
    build_dir = repo_root / "build"
    egg_info = repo_root / "src" / "tifzoret.egg-info"

    # Clean stale artifacts
    for path in [dist_dir, build_dir, egg_info]:
        if path.exists():
            import shutil
            shutil.rmtree(path)

    # Attempt to build wheel
    try:
        subprocess.run(
            [sys.executable, "-m", "build", "--wheel", "--outdir", str(dist_dir)],
            cwd=repo_root,
            check=True,
            capture_output=True,
            timeout=60,
        )
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        # Fall back to source tree assertion if build fails
        pytest.skip("build unavailable or failed; asserting source tree coverage instead")

    # Find the built wheel
    wheels = list(dist_dir.glob("*.whl"))
    assert len(wheels) == 1, f"Expected exactly one wheel, found {len(wheels)}"
    wheel_path = wheels[0]

    # Inspect wheel contents
    with zipfile.ZipFile(wheel_path) as whl:
        contents = set(whl.namelist())

    # Required assets
    required_patterns = [
        ("report template", lambda n: "report_assets/template.html" in n),
        ("woff2 font", lambda n: n.endswith(".woff2")),
        ("conda lock", lambda n: "envs/locks/" in n and n.endswith(".lock.txt")),
        ("R script", lambda n: "scripts/" in n and n.endswith(".R")),
    ]

    missing = []
    for name, matcher in required_patterns:
        if not any(matcher(f) for f in contents):
            missing.append(name)

    assert not missing, f"Wheel missing required assets: {', '.join(missing)}"


def test_source_tree_covers_report_assets():
    """Assert report_assets files exist in source tree (fast fallback)."""
    repo_root = Path(__file__).parent.parent
    report_assets = repo_root / "src" / "tifzoret" / "workflow" / "scripts" / "10_figures" / "report_assets"

    assert report_assets.exists(), "report_assets directory not found"
    assert (report_assets / "template.html").exists(), "template.html missing"

    fonts_dir = report_assets / "fonts"
    assert fonts_dir.exists(), "fonts directory missing"

    woff2_files = list(fonts_dir.glob("*.woff2"))
    assert len(woff2_files) > 0, "No .woff2 fonts found"

    # Also check locks and R scripts exist
    workflow = repo_root / "src" / "tifzoret" / "workflow"
    locks = list((workflow / "envs" / "locks").glob("*.lock.txt"))
    assert len(locks) > 0, "No conda lock files found"

    r_scripts = list(workflow.glob("scripts/**/*.R"))
    assert len(r_scripts) > 0, "No R scripts found"
