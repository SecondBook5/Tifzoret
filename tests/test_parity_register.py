from pathlib import Path
import pytest
from tools.parity.register_schema import validate_fragment, load_and_validate

FIX = Path(__file__).parent / "fixtures" / "parity"

def _entry(**over):
    base = dict(
        id="A1-001", track="A", area="de",
        verdict="REGRESSION", severity="S1",
        summary="apeglm shrinkage not applied",
        detail="frame uses normal shrinkage; reference uses apeglm",
        reference=[{"file": "workflow/stages/de/stage.R", "line": 42}],
        frame=[{"file": "src/bulk_rna_frame/workflow/scripts/de.R", "line": 88}],
        evidence="run", target="stats", status="open", commit=None,
    )
    base.update(over)
    return base

def test_valid_entry_has_no_errors():
    assert validate_fragment([_entry()]) == []

def test_missing_key_reported():
    bad = _entry(); del bad["severity"]
    errs = validate_fragment([bad])
    assert any("severity" in e for e in errs)

def test_bad_enum_reported():
    errs = validate_fragment([_entry(verdict="MAYBE")])
    assert any("verdict" in e and "MAYBE" in e for e in errs)

def test_fixed_inline_requires_commit():
    errs = validate_fragment([_entry(status="fixed-inline", commit=None)])
    assert any("commit" in e for e in errs)

def test_duplicate_ids_reported():
    errs = validate_fragment([_entry(id="A1-001"), _entry(id="A1-001")])
    assert any("duplicate" in e.lower() for e in errs)

def test_load_and_validate_raises_on_invalid_file():
    with pytest.raises(ValueError):
        load_and_validate(str(FIX / "invalid_fragment.yaml"))

def test_duplicate_ids_across_invalid_and_valid_entries():
    # Invalid entry (missing severity) with id="A1-001"
    invalid = _entry(id="A1-001"); del invalid["severity"]
    # Valid entry with same id
    valid = _entry(id="A1-001")
    errs = validate_fragment([invalid, valid])
    # Should report both the missing key AND the duplicate
    assert any("severity" in e for e in errs)
    assert any("duplicate" in e.lower() for e in errs)

def test_reference_file_and_line_value_types():
    # file must be a non-empty string
    errs = validate_fragment([_entry(reference=[{"file": None, "line": 42}])])
    assert any("file must be a non-empty string" in e for e in errs)

    errs = validate_fragment([_entry(reference=[{"file": "", "line": 42}])])
    assert any("file must be a non-empty string" in e for e in errs)

    # line must be an int
    errs = validate_fragment([_entry(reference=[{"file": "foo.R", "line": "not-a-number"}])])
    assert any("line must be an int" in e for e in errs)

    errs = validate_fragment([_entry(reference=[{"file": "foo.R", "line": None}])])
    assert any("line must be an int" in e for e in errs)
