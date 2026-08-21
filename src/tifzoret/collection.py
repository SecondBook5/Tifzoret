"""Collection-level signed Stouffer meta-analysis and leave-one-out checks."""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from statistics import NormalDist

from .config import ResolvedCollection


NORMAL = NormalDist()


def _read_de(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def _de_path(project_root: Path, contrast: str) -> Path:
    candidates = (
        project_root / "contrasts" / contrast / "analyses" / "de" / "tables" / "de_results.tsv",
        project_root / "contrasts" / contrast / "tables" / "de_results.tsv",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"DE results not found for contrast {contrast!r} beneath {project_root}")


def _signed_z(row: dict[str, str]) -> float | None:
    try:
        effect = float(row.get("log2_fold_change", row.get("log2FoldChange", "NA")))
        raw_p = float(row.get("p_value", row.get("pvalue", "NA")))
    except (KeyError, TypeError, ValueError):
        return None
    if math.isnan(effect) or math.isnan(raw_p):
        return None
    p = min(max(raw_p, 1e-300), 1.0)
    return math.copysign(NORMAL.inv_cdf(1 - p / 2), effect)


def _combine(values: list[tuple[float, float]]) -> tuple[float, float]:
    numerator = sum(weight * z for z, weight in values)
    denominator = math.sqrt(sum(weight * weight for _, weight in values))
    z = numerator / denominator
    p = 2 * (1 - NORMAL.cdf(abs(z)))
    return z, max(p, 0.0)


def _bh(rows: list[dict[str, object]]) -> None:
    ordered = sorted(enumerate(rows), key=lambda item: float(item[1]["meta_pvalue"]))
    total = len(rows)
    running = 1.0
    for reverse_rank, (index, row) in enumerate(reversed(ordered), start=1):
        rank = total - reverse_rank + 1
        running = min(running, float(row["meta_pvalue"]) * total / rank)
        rows[index]["meta_padj"] = running


def run_collection(collection: ResolvedCollection) -> Path:
    """Combine each study's per-gene DE evidence into one signed meta-analysis.

    Reads the chosen contrast's DE table from every member study, converts each
    gene's effect and p-value into a direction-signed z (positive = numerator
    minus denominator, preserved across studies), and combines the weighted
    z-scores with Stouffer's method into a two-sided meta p-value with BH
    correction. Genes seen in fewer than two studies are skipped; with three or
    more, an optional leave-one-out pass records the z-range and whether the
    direction stays stable. Writes ``meta_analysis.tsv`` plus a provenance
    ``manifest.json`` and returns the table path.
    """
    output = collection.result_root
    output.mkdir(parents=True, exist_ok=True)
    studies = []
    for spec, project in zip(collection.config["studies"], collection.projects, strict=True):
        path = _de_path(project.result_root, spec["contrast"])
        studies.append((spec, project, path, _read_de(path)))
    genes = sorted(set().union(*(set(data) for _, _, _, data in studies)))
    rows: list[dict[str, object]] = []
    leave_one_out = collection.config.get("methods", {}).get("leave_one_out", True)
    for gene_id in genes:
        values: list[tuple[str, float, float]] = []
        symbol = ""
        for spec, _, _, data in studies:
            row = data.get(gene_id)
            if row is None:
                continue
            z = _signed_z(row)
            if z is None:
                continue
            symbol = symbol or row.get("gene_symbol", "")
            values.append((spec["id"], z, float(spec.get("weight", 1.0))))
        if len(values) < 2:
            continue
        z, p = _combine([(value, weight) for _, value, weight in values])
        result: dict[str, object] = {
            "gene_id": gene_id,
            "gene_symbol": symbol,
            "studies": len(values),
            "meta_z": z,
            "meta_pvalue": p,
            "direction": "positive" if z > 0 else "negative" if z < 0 else "zero",
        }
        if leave_one_out and len(values) > 2:
            loo = {
                study_id: _combine([(other_z, other_weight) for other_id, other_z, other_weight in values if other_id != study_id])[0]
                for study_id, _, _ in values
            }
            result["leave_one_out_min_z"] = min(loo.values())
            result["leave_one_out_max_z"] = max(loo.values())
            result["leave_one_out_direction_stable"] = all(value * z > 0 for value in loo.values())
        rows.append(result)
    _bh(rows)
    fields = list(rows[0]) if rows else ["gene_id", "gene_symbol", "studies", "meta_z", "meta_pvalue", "meta_padj", "direction"]
    table = output / "meta_analysis.tsv"
    with table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "collection": collection.config["collection"],
        "method": "signed weighted Stouffer with two-sided p-values and BH correction",
        "contrast_semantics": "positive input effects are project numerator minus denominator",
        "studies": [
            {"id": spec["id"], "project": str(project.config_path), "contrast": spec["contrast"], "de": str(path)}
            for spec, project, path, _ in studies
        ],
        "genes_tested": len(rows),
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return table
