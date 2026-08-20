import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_engine_contains_no_reference_project_names():
    # The engine is project-agnostic: reference-study names must never leak into
    # any shipped workflow file. Each token is matched only when it is not
    # preceded by a letter, so genuine names ("cape", "cape_thoracic_duct",
    # "thoracicduct_cape_...") are caught while innocent superstrings
    # ("escape", "obligation", "intrathoracic") are not.
    forbidden = ("cape", "thoracic", "ligation", "lymphatic-flow-homeostasis")
    patterns = [re.compile(rf"(?<![a-z]){re.escape(term)}") for term in forbidden]
    for path in (ROOT / "src" / "bulk_rna_frame" / "workflow").rglob("*"):
        if path.is_file() and path.suffix in {".py", ".R", ".smk", ".yaml"}:
            text = path.read_text(errors="ignore").lower()
            hits = sorted({term for term, pattern in zip(forbidden, patterns) if pattern.search(text)})
            assert not hits, f"{path}: {hits}"


def test_all_figure_contracts_declare_pdf_and_png():
    snakefile = (ROOT / "src" / "bulk_rna_frame" / "workflow" / "Snakefile").read_text()
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
    scripts = ROOT / "src" / "bulk_rna_frame" / "workflow" / "scripts"
    for name in ("utils.R", "qc.R", "de.R", "pathways.R"):
        assert (scripts / name).is_file()
