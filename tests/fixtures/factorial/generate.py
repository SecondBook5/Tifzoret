"""Deterministic factorial fixtures with recorded ground truth (spec §11).

Every scenario is generated from a fixed seed, so a committed matrix is never
needed and a reviewer can read the truth rather than infer it. Factor and level
names are deliberately neutral (``factor_a``, ``a1``) because
``tests/test_repository.py`` guards the source tree against study-specific
vocabulary.

``truth`` records what each scenario injected: signal genes are named with the
``gene_signal_`` prefix and null genes ``gene_null_``, so a test asserts on gene
identity rather than on a p-value threshold.
"""

from __future__ import annotations

import csv
import hashlib
from collections.abc import Callable
from pathlib import Path

import numpy as np

BASE_MEAN = 200.0
DISPERSION_SIZE = 8.0
N_NULL = 240
N_SIGNAL = 60
SIGNAL_PREFIX = "gene_signal_"
NULL_PREFIX = "gene_null_"


def _seed(name: str) -> int:
    """Stable seed from scenario name for deterministic generation."""
    return int.from_bytes(hashlib.sha256(name.encode()).digest()[:4], "big")


def _cells(levels: dict[str, tuple[str, ...]]) -> list[dict[str, str]]:
    columns = list(levels)
    combos: list[dict[str, str]] = [{}]
    for column in columns:
        combos = [
            {**combo, column: level} for combo in combos for level in levels[column]
        ]
    return combos


def _samples(
    levels: dict[str, tuple[str, ...]],
    per_cell: int | Callable[[dict[str, str]], int],
    nuisance: Callable[[dict[str, str], int, int], str] | None = None,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    index = 0
    for cell in _cells(levels):
        count = per_cell(cell) if callable(per_cell) else per_cell
        for position in range(count):
            index += 1
            row = {"sample_id": f"s{index:03d}", **cell}
            if nuisance is not None:
                row["nuisance_z"] = nuisance(cell, position, index)
            rows.append(row)
    return rows


# Each scenario: levels, samples-per-cell, optional nuisance assignment, and a
# log2 effect function giving the signal genes' log2 shift for a given cell.
SCENARIOS: dict[str, dict] = {
    "two_by_two": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 4,
        "effect": lambda cell: 0.0,
        "truth": {"n_signal_genes": 0, "signal_prefix": SIGNAL_PREFIX, "note": "structural only"},
    },
    "two_by_three": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2", "b3")},
        "per_cell": 4,
        "effect": lambda cell: 0.0,
        "truth": {"n_signal_genes": 0, "signal_prefix": SIGNAL_PREFIX, "note": "structural only"},
    },
    "two_by_two_by_two": {
        "levels": {
            "factor_a": ("a1", "a2"),
            "factor_b": ("b1", "b2"),
            "factor_c": ("c1", "c2"),
        },
        "per_cell": 3,
        "effect": lambda cell: 0.0,
        "truth": {"n_signal_genes": 0, "signal_prefix": SIGNAL_PREFIX, "note": "structural only"},
    },
    "single_factor": {
        "levels": {"factor_a": ("a1", "a2")},
        "per_cell": 6,
        "effect": lambda cell: 1.5 if cell["factor_a"] == "a2" else 0.0,
        "truth": {
            "n_signal_genes": N_SIGNAL,
            "signal_prefix": SIGNAL_PREFIX,
            "interaction_log2": 0.0,
        },
    },
    "unbalanced_nuisance": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 4,
        # Varies within cells, and the mix differs between cells, so a naive
        # per-cell average of real model-matrix rows would leak a nuisance
        # difference into the estimand (spec §5).
        "nuisance": lambda cell, position, index: (
            "z1" if position < (3 if cell["factor_a"] == "a1" else 1) else "z2"
        ),
        "effect": lambda cell: 0.0,
        "truth": {"n_signal_genes": 0, "signal_prefix": SIGNAL_PREFIX},
    },
    "confounded_batch": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 4,
        "nuisance": lambda cell, position, index: (
            f"z_{cell['factor_a']}_{cell['factor_b']}"
        ),
        "effect": lambda cell: 0.0,
        "truth": {"n_signal_genes": 0, "signal_prefix": SIGNAL_PREFIX},
    },
    "strongly_unbalanced": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": lambda cell: 8 if cell == {"factor_a": "a1", "factor_b": "b1"} else 2,
        "effect": lambda cell: 0.0,
        "truth": {"n_signal_genes": 0, "signal_prefix": SIGNAL_PREFIX},
    },
    "null": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 4,
        "effect": lambda cell: 0.0,
        "truth": {"n_signal_genes": 0, "signal_prefix": SIGNAL_PREFIX, "interaction_log2": 0.0},
    },
    "factor_a_only": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 4,
        "effect": lambda cell: 1.5 if cell["factor_a"] == "a2" else 0.0,
        "truth": {"n_signal_genes": N_SIGNAL, "signal_prefix": SIGNAL_PREFIX, "interaction_log2": 0.0},
    },
    "factor_b_only": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 4,
        "effect": lambda cell: 1.5 if cell["factor_b"] == "b2" else 0.0,
        "truth": {"n_signal_genes": N_SIGNAL, "signal_prefix": SIGNAL_PREFIX, "interaction_log2": 0.0},
    },
    "baseline_shift_no_interaction": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 4,
        "effect": lambda cell: (
            (1.2 if cell["factor_a"] == "a2" else 0.0)
            + (0.9 if cell["factor_b"] == "b2" else 0.0)
        ),
        "truth": {"n_signal_genes": N_SIGNAL, "signal_prefix": SIGNAL_PREFIX, "interaction_log2": 0.0},
    },
    "positive_interaction": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 4,
        "effect": lambda cell: (
            2.0 if cell["factor_a"] == "a2" and cell["factor_b"] == "b2" else 0.0
        ),
        "truth": {"n_signal_genes": N_SIGNAL, "signal_prefix": SIGNAL_PREFIX, "interaction_log2": 2.0},
    },
    "negative_interaction": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 4,
        "effect": lambda cell: (
            (1.5 if cell["factor_a"] == "a2" else 0.0)
            + (1.5 if cell["factor_b"] == "b2" else 0.0)
            + (-2.0 if cell["factor_a"] == "a2" and cell["factor_b"] == "b2" else 0.0)
        ),
        "truth": {"n_signal_genes": N_SIGNAL, "signal_prefix": SIGNAL_PREFIX, "interaction_log2": -2.0},
    },
    "crossover_interaction": {
        # The spec §2 case: each arm moves 0.8 in opposite directions, so neither
        # simple effect is large, but the interaction is 1.6.
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 5,
        "effect": lambda cell: (
            0.8
            if cell["factor_a"] == "a1" and cell["factor_b"] == "b2"
            else (-0.8 if cell["factor_a"] == "a2" and cell["factor_b"] == "b2" else 0.0)
        ),
        "truth": {
            "n_signal_genes": N_SIGNAL,
            "signal_prefix": SIGNAL_PREFIX,
            "arm_a1_log2": 0.8,
            "arm_a2_log2": -0.8,
            "interaction_log2": -1.6,
        },
    },
    "interaction_without_arm_significance": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 6,
        "effect": lambda cell: (
            0.6
            if cell["factor_a"] == "a1" and cell["factor_b"] == "b2"
            else (-0.6 if cell["factor_a"] == "a2" and cell["factor_b"] == "b2" else 0.0)
        ),
        "truth": {
            "n_signal_genes": N_SIGNAL,
            "signal_prefix": SIGNAL_PREFIX,
            "interaction_log2": -1.2,
        },
    },
    "influential_outlier": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 4,
        "effect": lambda cell: 0.0,
        "outlier": True,
        "truth": {"n_signal_genes": 0, "signal_prefix": SIGNAL_PREFIX, "outlier_sample": "s001"},
    },
    "low_counts": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 4,
        "base_mean": 3.0,
        "effect": lambda cell: 1.5 if cell["factor_a"] == "a2" else 0.0,
        "truth": {"n_signal_genes": N_SIGNAL, "signal_prefix": SIGNAL_PREFIX},
    },
    "high_dispersion": {
        "levels": {"factor_a": ("a1", "a2"), "factor_b": ("b1", "b2")},
        "per_cell": 4,
        "size": 0.6,
        "effect": lambda cell: 1.5 if cell["factor_a"] == "a2" else 0.0,
        "truth": {"n_signal_genes": N_SIGNAL, "signal_prefix": SIGNAL_PREFIX},
    },
}

STRUCTURAL = (
    "two_by_two",
    "two_by_three",
    "two_by_two_by_two",
    "single_factor",
    "unbalanced_nuisance",
)
STATISTICAL = (
    "null",
    "factor_a_only",
    "factor_b_only",
    "positive_interaction",
    "negative_interaction",
    "crossover_interaction",
    "baseline_shift_no_interaction",
    "interaction_without_arm_significance",
    "strongly_unbalanced",
    "confounded_batch",
    "influential_outlier",
    "low_counts",
    "high_dispersion",
)


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def build(name: str, destination: Path) -> Path:
    """Write ``counts.tsv``, ``samples.tsv`` and ``annotation.tsv`` for ``name``."""
    scenario = SCENARIOS[name]
    rng = np.random.default_rng(_seed(name))
    levels = scenario["levels"]
    samples = _samples(levels, scenario["per_cell"], scenario.get("nuisance"))
    columns = ["sample_id", *levels]
    if scenario.get("nuisance") is not None:
        columns.append("nuisance_z")
    destination.mkdir(parents=True, exist_ok=True)
    _write_tsv(destination / "samples.tsv", columns, samples)

    base_mean = float(scenario.get("base_mean", BASE_MEAN))
    size = float(scenario.get("size", DISPERSION_SIZE))
    genes = [f"{NULL_PREFIX}{index:04d}" for index in range(N_NULL)] + [
        f"{SIGNAL_PREFIX}{index:04d}" for index in range(N_SIGNAL)
    ]
    count_rows: list[dict[str, object]] = []
    for gene in genes:
        row: dict[str, object] = {"gene_id": gene}
        is_signal = gene.startswith(SIGNAL_PREFIX)
        for sample in samples:
            cell = {column: sample[column] for column in levels}
            shift = scenario["effect"](cell) if is_signal else 0.0
            mean = base_mean * (2.0**shift)
            row[sample["sample_id"]] = int(
                rng.negative_binomial(size, size / (size + mean))
            )
        count_rows.append(row)
    if scenario.get("outlier"):
        target = scenario["truth"]["outlier_sample"]
        for row in count_rows[:20]:
            row[target] = int(base_mean * 60)
    _write_tsv(
        destination / "counts.tsv",
        ["gene_id", *[sample["sample_id"] for sample in samples]],
        count_rows,
    )
    _write_tsv(
        destination / "annotation.tsv",
        ["gene_id", "gene_symbol"],
        [{"gene_id": gene, "gene_symbol": gene.upper()} for gene in genes],
    )
    return destination
