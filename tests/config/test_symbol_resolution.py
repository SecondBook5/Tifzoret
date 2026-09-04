from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "01_inputs" / "materialize_inputs.py"


def _load_module():
    """Import the materialize_inputs script as a module (it adds src/ to sys.path itself)."""
    spec = importlib.util.spec_from_file_location("materialize_inputs_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


mi = _load_module()
DEFAULTS = mi.DEFAULT_PROVISIONAL_SYMBOL_PATTERNS


# --- _is_provisional -------------------------------------------------------

@pytest.mark.parametrize(
    "symbol, expected",
    [
        ("Gm12345", True),          # Ensembl/MGI predicted gene
        ("LOC108167921", True),     # NCBI uncharacterised locus
        ("2610305D13Rik", True),    # RIKEN cDNA clone
        ("Rik", True),              # bare Rik suffix
        ("Prox1", False),           # curated symbol
        ("Foxc2", False),
        ("ENSMUSG00000000001", False),
    ],
)
def test_is_provisional(symbol, expected):
    assert mi._is_provisional(symbol, DEFAULTS) is expected


# --- _resolve_symbol -------------------------------------------------------

def test_resolve_gtf_name_wins():
    assert mi._resolve_symbol("ENSMUSG1", "Prox1", "IgnoredSecondary", DEFAULTS) == (
        "Prox1",
        "gtf",
        "",
        False,
    )


def test_resolve_secondary_real_name_adopted():
    assert mi._resolve_symbol("ENSMUSG3", None, "Foxc2", DEFAULTS) == (
        "Foxc2",
        "secondary",
        "",
        False,
    )


def test_resolve_secondary_placeholder_keeps_accession():
    # A recovered placeholder is recorded but NEVER promoted to the displayed symbol.
    assert mi._resolve_symbol("ENSMUSG4", None, "Gm12345", DEFAULTS) == (
        "ENSMUSG4",
        "ensembl_id",
        "Gm12345",
        True,
    )


def test_resolve_nothing_recovered_retains_accession():
    assert mi._resolve_symbol("ENSMUSG2", None, None, DEFAULTS) == (
        "ENSMUSG2",
        "ensembl_id",
        "",
        False,
    )


# --- _load_symbol_map ------------------------------------------------------

def test_load_symbol_map_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        mi._load_symbol_map(tmp_path / "does_not_exist.tsv")


def test_load_symbol_map_dedup_lexically_smallest(tmp_path):
    path = tmp_path / "map.tsv"
    path.write_text(
        "# org.Mm.eg.db export, pinned\n"
        "gene_id\tsymbol\n"
        "ENSMUSG3\tFoxc2\n"
        "ENSMUSG4\t\n"          # blank symbol skipped
        "\tOrphan\n"            # blank id skipped
        "ENSMUSG5\tZeb1\n"
        "ENSMUSG5\tAaa1\n"      # duplicate id -> lexically smallest (Aaa1) wins
    )
    mapping = mi._load_symbol_map(path)
    assert mapping == {"ENSMUSG3": "Foxc2", "ENSMUSG5": "Aaa1"}


def test_load_symbol_map_requires_two_columns(tmp_path):
    path = tmp_path / "onecol.tsv"
    path.write_text("gene_id\nENSMUSG3\n")
    with pytest.raises(ValueError):
        mi._load_symbol_map(path)


# --- write_annotation ------------------------------------------------------

_GTF = (
    'chr1\tensembl\tgene\t100\t200\t.\t+\t.\tgene_id "ENSMUSG001"; gene_name "Prox1"; gene_biotype "protein_coding";\n'
    'chr1\tensembl\tgene\t300\t400\t.\t-\t.\tgene_id "ENSMUSG002"; gene_biotype "lncRNA";\n'
    'chr1\tensembl\tgene\t500\t600\t.\t+\t.\tgene_id "ENSMUSG003"; gene_biotype "protein_coding";\n'
    'chr1\tensembl\tgene\t700\t800\t.\t+\t.\tgene_id "ENSMUSG004"; gene_biotype "protein_coding";\n'
)
_GENE_IDS = ["ENSMUSG001", "ENSMUSG002", "ENSMUSG003", "ENSMUSG004"]


def _read(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    return header, rows


def test_write_annotation_disabled_is_legacy_seven_columns(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(_GTF)
    out = tmp_path / "annotation.tsv"
    mi.write_annotation(gtf, _GENE_IDS, out)  # resolution defaults to None (off)
    header, rows = _read(out)
    assert header == mi.ANNOTATION_BASE_FIELDS
    symbols = {r[0]: r[1] for r in rows}
    assert symbols["ENSMUSG001"] == "Prox1"          # gtf gene_name
    assert symbols["ENSMUSG002"] == "ENSMUSG002"     # no name -> accession fallback


def test_write_annotation_disabled_none_and_empty_dict_are_identical(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(_GTF)
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    mi.write_annotation(gtf, _GENE_IDS, a, resolution=None)
    mi.write_annotation(gtf, _GENE_IDS, b, resolution={"enabled": False})
    assert a.read_bytes() == b.read_bytes()


def test_write_annotation_enabled_resolves_and_adds_audit_columns(tmp_path, capsys):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(_GTF)
    secondary = tmp_path / "secondary.tsv"
    secondary.write_text(
        "gene_id\tsymbol\n"
        "ENSMUSG003\tFoxc2\n"       # real curated name -> adopted
        "ENSMUSG004\tGm12345\n"     # placeholder -> recorded, not promoted
    )
    out = tmp_path / "annotation.tsv"
    mi.write_annotation(
        gtf,
        _GENE_IDS,
        out,
        resolution={"enabled": True, "secondary_map": secondary},
    )
    header, rows = _read(out)
    assert header == [*mi.ANNOTATION_BASE_FIELDS, *mi.ANNOTATION_RESOLUTION_FIELDS]
    by_id = {r[0]: r for r in rows}

    # gene_symbol, symbol_source, provisional_name, provisional live at these indices:
    sym, src, prov_name, prov = 1, 7, 8, 9

    assert by_id["ENSMUSG001"][sym] == "Prox1"
    assert by_id["ENSMUSG001"][src] == "gtf"
    assert by_id["ENSMUSG001"][prov] == "FALSE"

    assert by_id["ENSMUSG002"][sym] == "ENSMUSG002"
    assert by_id["ENSMUSG002"][src] == "ensembl_id"
    assert by_id["ENSMUSG002"][prov] == "FALSE"

    assert by_id["ENSMUSG003"][sym] == "Foxc2"
    assert by_id["ENSMUSG003"][src] == "secondary"

    # placeholder is recorded but the stable accession stays as the symbol
    assert by_id["ENSMUSG004"][sym] == "ENSMUSG004"
    assert by_id["ENSMUSG004"][src] == "ensembl_id"
    assert by_id["ENSMUSG004"][prov_name] == "Gm12345"
    assert by_id["ENSMUSG004"][prov] == "TRUE"

    log = capsys.readouterr().out
    assert "[symbol-resolution] 4 genes" in log
    assert "gtf=1" in log
    assert "secondary_recovered=1" in log
    assert "provisional_flagged=1" in log
    assert "accession_retained=1" in log


def test_write_annotation_enabled_without_secondary_map(tmp_path):
    # Enabled but no secondary map: GTF names apply, everything else retains the accession.
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(_GTF)
    out = tmp_path / "annotation.tsv"
    mi.write_annotation(gtf, _GENE_IDS, out, resolution={"enabled": True})
    header, rows = _read(out)
    assert header == [*mi.ANNOTATION_BASE_FIELDS, *mi.ANNOTATION_RESOLUTION_FIELDS]
    by_id = {r[0]: r for r in rows}
    assert by_id["ENSMUSG001"][1] == "Prox1"
    assert by_id["ENSMUSG003"][1] == "ENSMUSG003"   # no secondary map -> accession retained
    assert by_id["ENSMUSG003"][7] == "ensembl_id"


def test_write_annotation_enabled_missing_secondary_map_raises(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(_GTF)
    out = tmp_path / "annotation.tsv"
    with pytest.raises(FileNotFoundError):
        mi.write_annotation(
            gtf,
            _GENE_IDS,
            out,
            resolution={"enabled": True, "secondary_map": tmp_path / "missing.tsv"},
        )


# --- _symbol_resolution_settings -------------------------------------------

def _fake_project(reference, config_dir):
    return SimpleNamespace(config={"reference": reference}, config_path=config_dir / "project.yaml")


def test_settings_none_when_block_absent(tmp_path):
    project = _fake_project({"genome_build": "GRCm39"}, tmp_path)
    assert mi._symbol_resolution_settings(project) is None


def test_settings_none_when_disabled(tmp_path):
    project = _fake_project({"symbol_resolution": {"enabled": False}}, tmp_path)
    assert mi._symbol_resolution_settings(project) is None


def test_settings_none_when_reference_missing(tmp_path):
    project = SimpleNamespace(config={}, config_path=tmp_path / "project.yaml")
    assert mi._symbol_resolution_settings(project) is None


def test_settings_resolves_secondary_map_relative_to_config(tmp_path):
    project = _fake_project(
        {"symbol_resolution": {"enabled": True, "secondary_map": "resources/map.tsv"}},
        tmp_path,
    )
    settings = mi._symbol_resolution_settings(project)
    assert settings["enabled"] is True
    # relative path is resolved against the config directory into an absolute path
    assert settings["secondary_map"] == (tmp_path / "resources" / "map.tsv").resolve()


def test_settings_secondary_map_none_stays_none(tmp_path):
    project = _fake_project({"symbol_resolution": {"enabled": True}}, tmp_path)
    settings = mi._symbol_resolution_settings(project)
    assert settings["enabled"] is True
    assert settings["secondary_map"] is None


# --- schema acceptance -----------------------------------------------------

def test_schema_accepts_and_constrains_symbol_resolution():
    from jsonschema import Draft202012Validator

    from tifzoret import config as cfg

    reference_schema = cfg._schema()["properties"]["reference"]
    validator = Draft202012Validator(reference_schema)

    # A full, valid block validates cleanly.
    valid = {
        "genome_build": "GRCm39",
        "annotation_release": 107,
        "symbol_resolution": {
            "enabled": True,
            "secondary_map": "resources/org_mm_eg_db_symbols.tsv",
            "provisional_patterns": [r"^Gm\d+$"],
        },
    }
    assert list(validator.iter_errors(valid)) == []

    # 'enabled' is required inside the block.
    assert list(validator.iter_errors({"genome_build": "GRCm39", "symbol_resolution": {}}))

    # Unknown keys inside the block are rejected (additionalProperties: false).
    assert list(
        validator.iter_errors(
            {"genome_build": "GRCm39", "symbol_resolution": {"enabled": True, "bogus": 1}}
        )
    )
