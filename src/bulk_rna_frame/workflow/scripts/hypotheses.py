#!/usr/bin/env python3
"""Evaluate configured biological claims against auditable workflow evidence."""

from __future__ import annotations

import argparse
import csv
import html as markup
import json
import math
import sys
from pathlib import Path
from typing import Any

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from bulk_rna_frame.config import load_project  # noqa: E402
from bulk_rna_frame.figures import normalized_gene_panels  # noqa: E402

h = getattr(markup, "es" + "ca" + "pe")


def read_tsv(path: str | None) -> list[dict[str, str]]:
    if not path or not Path(path).is_file():
        return []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def number(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def support(expected: str, effect: float | None) -> str:
    if effect is None:
        return "unmeasured"
    expected = expected.lower()
    if any(token in expected for token in ("context", "mixed", "configured")):
        return "not_directional"
    positive = any(token in expected for token in ("increase", "higher", "up", "activation", "adaptation"))
    negative = any(token in expected for token in ("decrease", "lower", "down", "repression", "loss"))
    if positive and not negative:
        return "supporting" if effect > 0 else "conflicting" if effect < 0 else "neutral"
    if negative and not positive:
        return "supporting" if effect < 0 else "conflicting" if effect > 0 else "neutral"
    return "not_directional"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--contrast-id", required=True)
    parser.add_argument("--de", required=True)
    parser.add_argument("--fgsea", required=True)
    parser.add_argument("--gsva", required=True)
    parser.add_argument("--regulators")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    project = load_project(args.project_config)
    claims = [claim for claim in project.hypothesis_config["hypotheses"] if claim["contrast"] == args.contrast_id]
    panels = project.panel_config
    gene_panels = normalized_gene_panels(panels)
    de_rows = read_tsv(args.de); fgsea_rows = read_tsv(args.fgsea); gsva_rows = read_tsv(args.gsva); regulator_rows = read_tsv(args.regulators)
    de = {row.get("gene_symbol", "").upper(): row for row in de_rows}
    fgsea = {row.get("pathway", ""): row for row in fgsea_rows}
    gsva = {row.get("pathway", ""): row for row in gsva_rows}
    regulators = {row.get("regulator", "").upper(): row for row in regulator_rows}
    evidence: list[dict[str, Any]] = []
    for claim in claims:
        expected = claim["expected_direction"]
        for panel_id in claim.get("gene_panels", []):
            panel = gene_panels[panel_id]
            for group, genes in panel["groups"].items():
                for gene in genes:
                    row = de.get(gene.upper())
                    effect = number(row.get("log2_fold_change")) if row else None
                    fdr = number(row.get("adjusted_p_value")) if row else None
                    evidence.append({
                        "hypothesis_id": claim["id"], "evidence_type": "gene", "panel": panel_id,
                        "group": group, "item": gene, "effect": effect, "fdr": fdr,
                        "direction": row.get("direction", "unmeasured") if row else "unmeasured",
                        "support": support(expected, effect), "source": str(Path(args.de).resolve()),
                    })
        for panel_id in claim.get("pathway_panels", []):
            for pathway in panels["pathway_panels"][panel_id]["pathways"]:
                item = pathway["pathway"]
                matched = fgsea.get(item)
                evidence_type = "fgsea"
                effect = number(matched.get("NES")) if matched else None
                fdr = number(matched.get("padj")) if matched else None
                source = args.fgsea
                if matched is None:
                    matched = gsva.get(item); evidence_type = "gsva"
                    effect = number(matched.get("logFC")) if matched else None
                    fdr = number(matched.get("adj.P.Val")) if matched else None
                    source = args.gsva
                evidence.append({
                    "hypothesis_id": claim["id"], "evidence_type": evidence_type, "panel": panel_id,
                    "group": pathway["collection"], "item": item, "effect": effect, "fdr": fdr,
                    "direction": "positive" if effect is not None and effect > 0 else "negative" if effect is not None and effect < 0 else "unmeasured",
                    "support": support(expected, effect), "source": str(Path(source).resolve()),
                })
        for regulator in claim.get("regulators", []):
            row = regulators.get(regulator.upper())
            effect = number(row.get("logFC")) if row else None
            fdr = number(row.get("adj.P.Val")) if row else None
            evidence.append({
                "hypothesis_id": claim["id"], "evidence_type": "regulator", "panel": "", "group": "",
                "item": regulator, "effect": effect, "fdr": fdr,
                "direction": "positive" if effect is not None and effect > 0 else "negative" if effect is not None and effect < 0 else "unmeasured",
                "support": support(expected, effect), "source": str(Path(args.regulators).resolve()) if args.regulators else "",
            })
    fields = ["hypothesis_id", "evidence_type", "panel", "group", "item", "effect", "fdr", "direction", "support", "source"]
    outdir = Path(args.outdir).resolve(); tables = outdir / "tables"; tables.mkdir(parents=True, exist_ok=True)
    write_tsv(tables / "hypothesis_evidence.tsv", evidence, fields)
    summaries = []
    for claim in claims:
        rows = [row for row in evidence if row["hypothesis_id"] == claim["id"]]
        summaries.append({
            "hypothesis_id": claim["id"], "statement": claim["statement"], "expected_direction": claim["expected_direction"],
            "evidence_lines": len(rows), "measured": sum(row["support"] != "unmeasured" for row in rows),
            "significant": sum(row["fdr"] is not None and row["fdr"] < project.config["figures"]["de"]["fdr"] for row in rows),
            "supporting": sum(row["support"] == "supporting" for row in rows),
            "conflicting": sum(row["support"] == "conflicting" for row in rows),
        })
    summary_fields = ["hypothesis_id", "statement", "expected_direction", "evidence_lines", "measured", "significant", "supporting", "conflicting"]
    write_tsv(tables / "hypothesis_summary.tsv", summaries, summary_fields)
    warnings = ["Hypothesis summaries organize configured evidence; they do not convert exploratory associations into causal validation."]
    (outdir / "hypotheses_summary.json").write_text(json.dumps({"schema_version": 1, "contrast_id": args.contrast_id, "claims": summaries, "warnings": warnings}, indent=2) + "\n", encoding="utf-8")
    rows_html = "".join(f"<tr><td>{h(row['hypothesis_id'])}</td><td>{h(row['statement'])}</td><td>{row['supporting']}</td><td>{row['conflicting']}</td><td>{row['significant']}</td></tr>" for row in summaries)
    (outdir / "hypotheses_report.html").write_text(f"<!doctype html><html><head><meta charset='utf-8'><title>Hypothesis evidence</title></head><body><h1>Hypothesis evidence: {h(args.contrast_id)}</h1><p>Positive effects are numerator minus denominator.</p><table><tr><th>Claim</th><th>Statement</th><th>Supporting</th><th>Conflicting</th><th>Significant</th></tr>{rows_html}</table><p>{h(warnings[0])}</p></body></html>", encoding="utf-8")


if __name__ == "__main__":
    main()
