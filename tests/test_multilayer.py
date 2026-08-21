from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


SCRIPT = (
    Path(__file__).parents[1]
    / "src"
    / "tifzoret"
    / "workflow"
    / "scripts"
    / "multilayer.py"
)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_multilayer_uses_string_edge_endpoints(tmp_path: Path) -> None:
    grn = tmp_path / "grn.tsv"
    hubs = tmp_path / "hubs.tsv"
    string_up = tmp_path / "string_up.tsv"
    string_down = tmp_path / "string_down.tsv"
    outdir = tmp_path / "out"

    write_tsv(grn, ["source", "target"], [{"source": "Tf1", "target": "GeneA"}])
    write_tsv(hubs, ["gene_symbol"], [{"gene_symbol": "GeneA"}, {"gene_symbol": "GeneB"}])
    write_tsv(
        string_up,
        ["source", "target", "combined_score"],
        [{"source": "GeneA", "target": "GeneC", "combined_score": 0.9}],
    )
    write_tsv(
        string_down,
        ["source", "target", "combined_score"],
        [{"source": "GeneD", "target": "GeneB", "combined_score": 0.8}],
    )

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--grn-edges",
            str(grn),
            "--wgcna-hubs",
            str(hubs),
            "--string-up",
            str(string_up),
            "--string-down",
            str(string_down),
            "--contrast-id",
            "treated_vs_control",
            "--outdir",
            str(outdir),
        ],
        check=True,
    )

    with (outdir / "tables" / "multilayer_nodes.tsv").open(encoding="utf-8") as handle:
        nodes = {row["gene_symbol"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert nodes["GENEA"]["layers"] == "regulatory;coexpression;string_association"
    assert nodes["GENEA"]["triangulated"] == "True"
    assert nodes["GENEB"]["layers"] == "coexpression;string_association"
    assert {"GENEC", "GENED", "TF1"}.issubset(nodes)

    summary = json.loads((outdir / "multilayer_summary.json").read_text(encoding="utf-8"))
    assert summary["contrast_id"] == "treated_vs_control"
    assert summary["triangulated"] == 2
    assert (outdir / "figures" / "multilayer_network.pdf").stat().st_size > 0
    assert (outdir / "figures" / "multilayer_network.png").stat().st_size > 0
