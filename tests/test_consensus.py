from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "consensus.py"
PROJECT = ROOT / "src" / "tifzoret" / "templates" / "minimal" / "project.yaml"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _de_results(path: Path, contrast_id: str, rows: list[dict[str, object]]) -> None:
    write_tsv(
        path,
        ["contrast_id", "gene_id", "gene_symbol", "direction", "log2_fold_change"],
        [{"contrast_id": contrast_id, **row} for row in rows],
    )


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_consensus_scores_agreement_and_flags_discordance(tmp_path: Path) -> None:
    """Cross-contrast consensus takes direction verbatim from each contrast's
    significance call: a gene agreeing in >= min_contrasts is a consensus gene,
    a gene up in one contrast and down in another is discordant (score = the
    larger of its up/down counts, so it never reaches consensus on a split)."""
    de_a = tmp_path / "de_a.tsv"
    de_b = tmp_path / "de_b.tsv"
    outdir = tmp_path / "out"
    # G1 up in both, G3 down in both (both consensus); G2 up then down (discordant).
    _de_results(de_a, "cond_a", [
        {"gene_id": "G1", "gene_symbol": "Aaa", "direction": "up_in_numerator", "log2_fold_change": 2.0},
        {"gene_id": "G2", "gene_symbol": "Bbb", "direction": "up_in_numerator", "log2_fold_change": 1.5},
        {"gene_id": "G3", "gene_symbol": "Ccc", "direction": "down_in_numerator", "log2_fold_change": -2.0},
    ])
    _de_results(de_b, "cond_b", [
        {"gene_id": "G1", "gene_symbol": "Aaa", "direction": "up_in_numerator", "log2_fold_change": 1.8},
        {"gene_id": "G2", "gene_symbol": "Bbb", "direction": "down_in_numerator", "log2_fold_change": -1.2},
        {"gene_id": "G3", "gene_symbol": "Ccc", "direction": "down_in_numerator", "log2_fold_change": -1.9},
    ])

    subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--project-config", str(PROJECT),
            "--de", str(de_b), str(de_a),  # deliberately reversed: output order must be deterministic
            "--outdir", str(outdir),
        ],
        check=True,
    )

    summary = json.loads((outdir / "consensus_summary.json").read_text(encoding="utf-8"))
    assert summary["contrasts"] == ["cond_a", "cond_b"]  # sorted regardless of CLI order
    assert summary["n_contrasts"] == 2
    assert summary["genes_considered"] == 3
    assert summary["consensus_genes"] == 2
    assert summary["discordant_genes"] == 1
    assert summary["intersections"] == 1  # every gene significant in exactly {cond_a, cond_b}

    consensus = {row["gene_symbol"]: row for row in _read_tsv(outdir / "tables" / "consensus_genes.tsv")}
    assert set(consensus) == {"Aaa", "Ccc"}
    assert consensus["Aaa"]["sign_consistent"] == "True"
    assert consensus["Ccc"]["n_down"] == "2"

    membership = {row["gene_symbol"]: row for row in _read_tsv(outdir / "tables" / "consensus_membership.tsv")}
    assert membership["Bbb"]["sign_consistent"] == "False"
    assert (membership["Bbb"]["n_up"], membership["Bbb"]["n_down"]) == ("1", "1")

    overlap = _read_tsv(outdir / "tables" / "contrast_overlap.tsv")
    assert len(overlap) == 1
    assert float(overlap[0]["jaccard"]) == 1.0
    assert abs(float(overlap[0]["sign_agreement"]) - 2 / 3) < 1e-9  # G1,G3 agree; G2 disagrees

    for stem in ("consensus_upset", "consensus_sign_heatmap"):
        for suffix in ("pdf", "png"):
            assert (outdir / "figures" / f"{stem}.{suffix}").stat().st_size > 0
