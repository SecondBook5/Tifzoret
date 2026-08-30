"""Every stage script lives in exactly one phase folder — no orphans, no flat leftovers."""
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "src" / "tifzoret" / "workflow" / "scripts"

EXPECTED = {
    "01_inputs": {"materialize_inputs.py", "resources.R", "export_symbol_map.R"},
    "02_qc": {"qc.R", "batch.R", "sva.R", "variancepartition.R"},
    "03_differential": {"de.R", "de_confirm.R", "omnibus.R", "factorial.R"},
    "04_enrichment": {"pathways.R", "ontology.R", "spia.R", "enrichment_map.py"},
    "05_composition": {"composition.R", "deconvolution.R"},
    "06_regulators": {"regulators.R", "grn.py", "grn_radial.R"},
    "07_networks": {"networks.py", "string_network.R", "string_figures.R",
                    "wgcna.R", "curvature.py", "multilayer.py"},
    "08_causal": {"mediation.R", "power.py"},
    "09_synthesis": {"consensus.py", "hypotheses.py"},
    "10_figures": {"publication.R", "assemble.py", "report.py", "manifest.py", "front_door.py"},
}


def test_utils_stays_flat():
    assert (SCRIPTS / "utils.R").is_file()


def test_each_folder_has_exactly_its_scripts():
    for folder, names in EXPECTED.items():
        found = {p.name for p in (SCRIPTS / folder).glob("*") if p.suffix in {".R", ".py"}}
        assert found == names, f"{folder}: expected {names}, found {found}"


def test_no_stage_script_left_flat():
    flat = {p.name for p in SCRIPTS.glob("*") if p.suffix in {".R", ".py"}}
    assert flat == {"utils.R"}, f"unexpected flat scripts: {flat - {'utils.R'}}"


def test_report_assets_moved_with_report():
    assert (SCRIPTS / "10_figures" / "report_assets" / "template.html").is_file()
