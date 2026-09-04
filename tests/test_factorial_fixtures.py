"""The seeded factorial fixture generator (spec §11)."""

from __future__ import annotations

from pathlib import Path

import pytest

from _factorial_support import build_fixture, fixtures, read_tsv_rows

ROOT = Path(__file__).resolve().parents[1]


def test_every_declared_scenario_builds(tmp_path):
    module = fixtures()
    for name in module.STRUCTURAL + module.STATISTICAL:
        destination = module.build(name, tmp_path / name)
        assert (destination / "counts.tsv").is_file(), name
        assert (destination / "samples.tsv").is_file(), name
        assert (destination / "annotation.tsv").is_file(), name


def test_generation_is_deterministic(tmp_path):
    module = fixtures()
    first = module.build("positive_interaction", tmp_path / "a")
    second = module.build("positive_interaction", tmp_path / "b")
    assert (first / "counts.tsv").read_bytes() == (second / "counts.tsv").read_bytes()


def test_two_by_two_has_four_balanced_cells(tmp_path):
    module = fixtures()
    rows = read_tsv_rows(module.build("two_by_two", tmp_path / "x") / "samples.tsv")
    cells: dict[tuple[str, str], int] = {}
    for row in rows:
        key = (row["factor_a"], row["factor_b"])
        cells[key] = cells.get(key, 0) + 1
    assert len(cells) == 4
    assert len(set(cells.values())) == 1


def test_two_by_three_has_six_cells(tmp_path):
    module = fixtures()
    rows = read_tsv_rows(module.build("two_by_three", tmp_path / "x") / "samples.tsv")
    assert len({(r["factor_a"], r["factor_b"]) for r in rows}) == 6


def test_three_factor_fixture_has_eight_cells(tmp_path):
    module = fixtures()
    rows = read_tsv_rows(module.build("two_by_two_by_two", tmp_path / "x") / "samples.tsv")
    assert len({(r["factor_a"], r["factor_b"], r["factor_c"]) for r in rows}) == 8


def test_single_factor_fixture_has_one_cell_column(tmp_path):
    module = fixtures()
    rows = read_tsv_rows(module.build("single_factor", tmp_path / "x") / "samples.tsv")
    assert "factor_b" not in rows[0]
    assert len({r["factor_a"] for r in rows}) == 2


def test_unbalanced_nuisance_varies_within_cells(tmp_path):
    module = fixtures()
    rows = read_tsv_rows(module.build("unbalanced_nuisance", tmp_path / "x") / "samples.tsv")
    per_cell: dict[tuple[str, str], list[str]] = {}
    for row in rows:
        per_cell.setdefault((row["factor_a"], row["factor_b"]), []).append(row["nuisance_z"])
    # At least one cell varies within
    assert any(len(set(values)) > 1 for values in per_cell.values())
    # Mix differs between cells (a1 cells have different z1:z2 ratio than a2 cells)
    ratios = {cell: values.count("z1") / len(values) for cell, values in per_cell.items()}
    a1_ratios = {ratio for cell, ratio in ratios.items() if cell[0] == "a1"}
    a2_ratios = {ratio for cell, ratio in ratios.items() if cell[0] == "a2"}
    assert a1_ratios != a2_ratios


def test_confounded_batch_is_constant_within_every_cell(tmp_path):
    module = fixtures()
    rows = read_tsv_rows(module.build("confounded_batch", tmp_path / "x") / "samples.tsv")
    per_cell: dict[tuple[str, str], set[str]] = {}
    for row in rows:
        per_cell.setdefault((row["factor_a"], row["factor_b"]), set()).add(row["nuisance_z"])
    assert all(len(values) == 1 for values in per_cell.values())


def test_scenarios_record_their_truth_for_assertions(tmp_path):
    module = fixtures()
    truth = module.SCENARIOS["crossover_interaction"]["truth"]
    assert truth["n_signal_genes"] > 0
    assert truth["signal_prefix"] == "gene_signal_"


def test_counts_are_non_negative_integers(tmp_path):
    module = fixtures()
    rows = read_tsv_rows(module.build("high_dispersion", tmp_path / "x") / "counts.tsv")
    for row in rows:
        for key, value in row.items():
            if key == "gene_id":
                continue
            assert int(value) >= 0


def test_annotation_covers_every_gene(tmp_path):
    module = fixtures()
    destination = module.build("null", tmp_path / "x")
    counts = {row["gene_id"] for row in read_tsv_rows(destination / "counts.tsv")}
    annotation = {row["gene_id"] for row in read_tsv_rows(destination / "annotation.tsv")}
    assert counts == annotation


def test_rejects_an_unknown_scenario(tmp_path):
    module = fixtures()
    with pytest.raises(KeyError):
        module.build("not_a_scenario", tmp_path / "x")


def test_sample_ids_are_unique_and_match_count_columns(tmp_path):
    module = fixtures()
    destination = module.build("two_by_two", tmp_path / "x")
    samples = [row["sample_id"] for row in read_tsv_rows(destination / "samples.tsv")]
    with (destination / "counts.tsv").open(encoding="utf-8") as handle:
        header = handle.readline().strip().split("\t")
    assert len(samples) == len(set(samples))
    assert header[1:] == samples
