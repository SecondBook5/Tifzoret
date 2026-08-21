import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SELF = Path(__file__).resolve()

# Directories/artifacts that are never part of the engine's source tree.
_SKIP_DIRS = {"__pycache__", ".pytest_cache"}
_SCAN_DIRS = ("src", "docs", "tests")
# Root files that must also stay project-agnostic: the two docs plus the
# packaging/citation metadata, where a study name is just as damaging as in
# source (a leaked token in pyproject/MANIFEST/environment/CITATION ships too).
_ROOT_FILES = (
    "README.md",
    "CHANGELOG.md",
    "pyproject.toml",
    "MANIFEST.in",
    "environment.yaml",
    "CITATION.cff",
)
_SCANNED_SUFFIXES = {".py", ".R", ".smk", ".yaml", ".yml", ".md", ".sh"}


def _engine_source_files():
    """Yield every tracked source file that must stay project-agnostic.

    Walks the engine's source directories plus its root docs, skipping this
    guard (which necessarily contains the forbidden tokens as literals) and
    build/cache artifacts (``__pycache__``, ``*.egg-info``). Deliberately does
    not descend into run outputs (``results/``, ``.snakemake/``) so an
    incidental local run cannot make the check pass or fail spuriously.
    """
    for name in _SCAN_DIRS:
        for path in (ROOT / name).rglob("*"):
            if not path.is_file() or path.suffix not in _SCANNED_SUFFIXES:
                continue
            if path.resolve() == SELF:
                continue
            if any(part in _SKIP_DIRS or part.endswith(".egg-info") for part in path.parts):
                continue
            yield path
    for name in _ROOT_FILES:
        path = ROOT / name
        if path.is_file():
            yield path


def test_engine_contains_no_reference_project_names():
    # The engine is project-agnostic: reference-study names must never leak into
    # ANY source file (code, docs, tests, config) — not just workflow/. Each
    # token is matched only when not preceded by a letter, so genuine names
    # ("cape", "cape_thoracic_duct", "thoracicduct_cape_...") are caught while
    # innocent superstrings ("escape", "obligation", "intrathoracic") are not.
    # The list spans every planned reference cohort so a new study's name cannot
    # slip in via config or docs. Tokens are chosen to be fingerprint-specific:
    # "rela_ko" (not bare "rela", which would flag "related"/"relative") and
    # "nfkb"/"xizhao"/"pten"/"taxol" have no innocent-superstring collisions.
    forbidden = (
        "cape",
        "thoracic",
        "ligation",
        "lymphatic-flow-homeostasis",
        "xizhao",
        "pten",
        "taxol",
        "nfkb",
        "rela_ko",
    )
    patterns = [re.compile(rf"(?<![a-z]){re.escape(term)}") for term in forbidden]
    for path in _engine_source_files():
        text = path.read_text(errors="ignore").lower()
        hits = sorted({term for term, pattern in zip(forbidden, patterns) if pattern.search(text)})
        assert not hits, f"{path.relative_to(ROOT)}: {hits}"


def test_all_figure_contracts_declare_pdf_and_png():
    snakefile = (ROOT / "src" / "tifzoret" / "workflow" / "Snakefile").read_text()
    for stem in (
        "pca",
        "sample_correlation",
        "volcano",
        "de_heatmap",
        "ora_bidirectional",
        "gsva_heatmap",
        "gsea_curves",
    ):
        assert f'{stem}.pdf' in snakefile
        assert f'{stem}.png' in snakefile


def test_r_scripts_parse():
    # Full execution is covered by the synthetic acceptance workflow. This test
    # ensures every shipped R entry point remains present for that check.
    scripts = ROOT / "src" / "tifzoret" / "workflow" / "scripts"
    for name in ("utils.R", "qc.R", "de.R", "pathways.R"):
        assert (scripts / name).is_file()
