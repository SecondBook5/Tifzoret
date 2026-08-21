from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "enrichment_map.py"
PROJECT = ROOT / "src" / "tifzoret" / "templates" / "minimal" / "project.yaml"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_enrichment_map_links_shared_gene_terms_and_clusters(tmp_path: Path) -> None:
    """Enriched terms become nodes; an edge is kept when their driver-gene
    signatures overlap by >= min_similarity (Jaccard). Two terms sharing 3 of 5
    genes (J=0.6) are linked into one community; a disjoint significant term is a
    separate community. A sub-threshold-FDR term contributes no node."""
    fgsea = tmp_path / "fgsea.tsv"
    ora = tmp_path / "ora.tsv"
    outdir = tmp_path / "out"

    write_tsv(
        fgsea,
        ["pathway", "pathway_label", "gene_set_source", "padj", "NES", "direction", "leadingEdge"],
        [
            {"pathway": "TERM_A", "pathway_label": "Term A", "gene_set_source": "custom",
             "padj": 0.001, "NES": 2.0, "direction": "up_in_numerator", "leadingEdge": "G1;G2;G3;G4"},
            {"pathway": "TERM_B", "pathway_label": "Term B", "gene_set_source": "custom",
             "padj": 0.010, "NES": 1.8, "direction": "up_in_numerator", "leadingEdge": "G1;G2;G3;G5"},
        ],
    )
    write_tsv(
        ora,
        ["pathway", "pathway_label", "gene_set_source", "adjusted_p_value", "count", "direction", "overlap_genes"],
        [
            {"pathway": "TERM_C", "pathway_label": "Term C", "gene_set_source": "custom",
             "adjusted_p_value": 0.020, "count": 5, "direction": "down_in_numerator", "overlap_genes": "G10;G11;G12"},
            {"pathway": "TERM_D", "pathway_label": "Term D", "gene_set_source": "custom",
             "adjusted_p_value": 0.500, "count": 4, "direction": "up_in_numerator", "overlap_genes": "G20;G21"},
        ],
    )

    subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--project-config", str(PROJECT),
            "--fgsea", str(fgsea),
            "--ora", str(ora),
            "--contrast-id", "treated_vs_control",
            "--outdir", str(outdir),
        ],
        check=True,
    )

    summary = json.loads((outdir / "enrichment_map_summary.json").read_text(encoding="utf-8"))
    assert summary["contrast_id"] == "treated_vs_control"
    assert summary["terms"] == 3  # TERM_D dropped (adjusted_p_value >= FDR)
    assert summary["edges"] == 1
    assert summary["clusters"] == 2

    edges = _read_tsv(outdir / "tables" / "enrichment_map_edges.tsv")
    assert len(edges) == 1
    assert {edges[0]["term_a"], edges[0]["term_b"]} == {"TERM_A", "TERM_B"}
    assert edges[0]["n_shared"] == "3"
    assert abs(float(edges[0]["jaccard"]) - 0.6) < 1e-9
    assert edges[0]["shared_genes"] == "G1;G2;G3"

    nodes = {row["term"]: row for row in _read_tsv(outdir / "tables" / "enrichment_map_nodes.tsv")}
    assert set(nodes) == {"TERM_A", "TERM_B", "TERM_C"}
    assert nodes["TERM_A"]["cluster"] == nodes["TERM_B"]["cluster"]  # linked -> same community
    assert nodes["TERM_C"]["cluster"] != nodes["TERM_A"]["cluster"]  # disjoint -> its own
    assert nodes["TERM_C"]["direction"] == "down_in_numerator"

    assert len(_read_tsv(outdir / "tables" / "enrichment_map_clusters.tsv")) == 2
    for suffix in ("pdf", "png"):
        assert (outdir / "figures" / f"enrichment_map.{suffix}").stat().st_size > 0
