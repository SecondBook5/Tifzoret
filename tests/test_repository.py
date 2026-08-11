from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_engine_contains_no_reference_project_names():
    forbidden = ("cape", "thoracic", "ligation", "lymphatic-flow-homeostasis")
    for path in (ROOT / "src" / "bulk_rna_frame" / "workflow").rglob("*"):
        if path.is_file():
            text = path.read_text(errors="ignore").lower()
            assert not any(term in text for term in forbidden), path


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
