"""Numerical and artifact verification for migrated BulkRNAFrame runs."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from PIL import Image
from pypdf import PdfReader

from .config import ResolvedProject


FEATURECOUNTS_FIELDS = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
DE_FIELDS = {
    "base_mean": ("raw", "baseMean", 1e-10),
    "log2_fold_change_raw": ("raw", "log2FoldChange", 1e-5),
    "log2_fold_change": ("shrunken", "log2FoldChange", 1e-5),
    "statistic": ("raw", "stat", 1e-5),
    "p_value": ("raw", "pvalue", 1e-6),
    "adjusted_p_value": ("raw", "padj", 1e-6),
}


@dataclass(frozen=True)
class VerificationResult:
    passed: bool
    report: dict[str, Any]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def _read_tsv_comments(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def _number(value: str | None) -> float:
    if value is None or value in {"", "NA", "NaN", "nan"}:
        return math.nan
    return float(value)


def _decision(padj: float, lfc: float, fdr: float, abs_lfc: float) -> str:
    if math.isnan(padj) or math.isnan(lfc):
        return "not_testable"
    if padj >= fdr:
        return "not_significant"
    if lfc >= abs_lfc:
        return "up_in_numerator"
    if lfc <= -abs_lfc:
        return "down_in_numerator"
    return "significant_below_lfc_threshold"


def _numeric_equal(left: str, right: str, *, atol: float, rtol: float) -> bool:
    if left == right:
        return True
    try:
        a, b = float(left), float(right)
    except ValueError:
        return False
    if math.isnan(a) and math.isnan(b):
        return True
    return math.isclose(a, b, abs_tol=atol, rel_tol=rtol)


def _compare_tsv(reference: Path, candidate: Path, *, atol: float, rtol: float) -> dict[str, Any]:
    reference_header, reference_rows = _read_tsv(reference)
    candidate_header, candidate_rows = _read_tsv(candidate)
    differences: list[dict[str, Any]] = []
    if reference_header != candidate_header:
        differences.append({"kind": "header", "reference": reference_header, "candidate": candidate_header})
    key = "gene_id" if "gene_id" in reference_header else None
    if key and key in candidate_header:
        reference_map = {row[key]: row for row in reference_rows}
        candidate_map = {row[key]: row for row in candidate_rows}
        missing = sorted(set(reference_map) - set(candidate_map))
        extra = sorted(set(candidate_map) - set(reference_map))
        if missing:
            differences.append({"kind": "missing_rows", "count": len(missing), "first": missing[:10]})
        if extra:
            differences.append({"kind": "extra_rows", "count": len(extra), "first": extra[:10]})
        pairs = ((row_key, reference_map[row_key], candidate_map[row_key]) for row_key in sorted(set(reference_map) & set(candidate_map)))
    else:
        if len(reference_rows) != len(candidate_rows):
            differences.append({"kind": "row_count", "reference": len(reference_rows), "candidate": len(candidate_rows)})
        pairs = ((str(index + 2), a, b) for index, (a, b) in enumerate(zip(reference_rows, candidate_rows)))
    differing_cells = 0
    examples: list[dict[str, str]] = []
    shared_columns = [column for column in reference_header if column in candidate_header]
    for row_key, reference_row, candidate_row in pairs:
        for column in shared_columns:
            if not _numeric_equal(reference_row[column], candidate_row[column], atol=atol, rtol=rtol):
                differing_cells += 1
                if len(examples) < 20:
                    examples.append({
                        "row": row_key,
                        "column": column,
                        "reference": reference_row[column],
                        "candidate": candidate_row[column],
                    })
    if differing_cells:
        differences.append({"kind": "cells", "count": differing_cells, "examples": examples})
    return {"passed": not differences, "differences": differences}


def verify_runs(reference: str | Path, candidate: str | Path, *, atol: float = 1e-8, rtol: float = 1e-6) -> VerificationResult:
    reference_root = Path(reference).expanduser().resolve()
    candidate_root = Path(candidate).expanduser().resolve()
    if not reference_root.is_dir() or not candidate_root.is_dir():
        missing = [str(path) for path in (reference_root, candidate_root) if not path.is_dir()]
        return VerificationResult(False, {"passed": False, "error": "run directory missing", "paths": missing})
    reference_tables = {
        path.relative_to(reference_root): path
        for path in reference_root.rglob("*.tsv")
    }
    candidate_tables = {
        path.relative_to(candidate_root): path
        for path in candidate_root.rglob("*.tsv")
    }
    missing = sorted(str(path) for path in set(reference_tables) - set(candidate_tables))
    extra = sorted(str(path) for path in set(candidate_tables) - set(reference_tables))
    comparisons: dict[str, Any] = {}
    for relative in sorted(set(reference_tables) & set(candidate_tables)):
        comparisons[str(relative)] = _compare_tsv(
            reference_tables[relative], candidate_tables[relative], atol=atol, rtol=rtol
        )
    failed = [name for name, comparison in comparisons.items() if not comparison["passed"]]
    report = {
        "schema_version": 1,
        "passed": not missing and not failed,
        "reference": str(reference_root),
        "candidate": str(candidate_root),
        "tolerances": {"absolute": atol, "relative": rtol},
        "missing_tables": missing,
        "candidate_only_tables": extra,
        "failed_tables": failed,
        "comparisons": comparisons,
    }
    return VerificationResult(bool(report["passed"]), report)


def _verify_legacy_counts(project: ResolvedProject, reference: Path, candidate: Path) -> dict[str, Any]:
    old_path = reference / ".cache" / "counts" / "gene_counts.tsv"
    new_path = candidate / "inputs" / "counts.tsv"
    if not old_path.is_file() or not new_path.is_file():
        return {"passed": False, "error": "count table missing", "reference": str(old_path), "candidate": str(new_path)}
    old_header, old_rows = _read_tsv_comments(old_path)
    new_header, new_rows = _read_tsv(new_path)
    sample_ids = [row["sample_id"] for row in project.sample_rows]
    bam_to_sample = {
        Path(row.get("bam", "")).name: row["sample_id"]
        for row in project.sample_rows
        if row.get("bam")
    }
    old_columns = [column for column in old_header if column not in FEATURECOUNTS_FIELDS]
    mapping = {
        column: bam_to_sample[Path(column).name]
        for column in old_columns
        if Path(column).name in bam_to_sample
    }
    if set(mapping.values()) != set(sample_ids):
        return {
            "passed": False,
            "error": "could not map selected samples to legacy featureCounts columns",
            "expected_samples": sample_ids,
            "mapped_samples": sorted(mapping.values()),
        }
    header_ok = new_header == ["gene_id", *sample_ids]
    old = {row["Geneid"]: row for row in old_rows}
    new = {row["gene_id"]: row for row in new_rows}
    missing = sorted(set(old) - set(new))
    extra = sorted(set(new) - set(old))
    summaries = {}
    passed = header_ok and not missing and not extra
    for sample_id in sample_ids:
        old_column = next(column for column, mapped in mapping.items() if mapped == sample_id)
        mismatches = 0
        maximum = 0
        for gene_id in set(old) & set(new):
            difference = abs(int(old[gene_id][old_column]) - int(new[gene_id][sample_id]))
            mismatches += difference != 0
            maximum = max(maximum, difference)
        summaries[sample_id] = {"mismatched_genes": mismatches, "maximum_absolute_difference": maximum}
        passed = passed and mismatches == 0
    return {
        "passed": passed,
        "canonical_header": header_ok,
        "reference_genes": len(old),
        "candidate_genes": len(new),
        "missing_candidate_genes": missing,
        "extra_candidate_genes": extra,
        "per_sample": summaries,
    }


def _verify_legacy_de(
    project: ResolvedProject,
    reference: Path,
    candidate: Path,
    *,
    rtol: float,
) -> dict[str, Any]:
    raw_path = reference / ".cache" / "de" / "de_raw.tsv"
    shrunken_path = reference / ".cache" / "de" / "de_shrunken.tsv"
    if not raw_path.is_file() or not shrunken_path.is_file():
        return {"passed": False, "error": "legacy DE tables missing"}
    _, raw_rows = _read_tsv(raw_path)
    _, shrunken_rows = _read_tsv(shrunken_path)
    raw = {row["gene_id"]: row for row in raw_rows}
    shrunken = {row["gene_id"]: row for row in shrunken_rows}
    reference_genes = set(raw) & set(shrunken)
    fdr = float(project.config["figures"]["de"]["fdr"])
    abs_lfc = float(project.config["figures"]["de"]["abs_log2fc"])
    contrasts = {}
    all_passed = True
    for contrast in project.contrast_rows:
        contrast_id = contrast["contrast_id"]
        candidate_path = candidate / "contrasts" / contrast_id / "analyses" / "de" / "tables" / "de_results.tsv"
        if not candidate_path.is_file():
            contrasts[contrast_id] = {"passed": False, "error": f"candidate DE table missing: {candidate_path}"}
            all_passed = False
            continue
        _, candidate_rows = _read_tsv(candidate_path)
        new = {row["gene_id"]: row for row in candidate_rows}
        missing = sorted(reference_genes - set(new))
        extra = sorted(set(new) - reference_genes)
        shared = reference_genes & set(new)
        field_results = {}
        passed = not missing and not extra
        sources = {"raw": raw, "shrunken": shrunken}
        for candidate_field, (source, legacy_field, atol) in DE_FIELDS.items():
            mismatches = 0
            maximum = 0.0
            for gene_id in shared:
                old_value = _number(sources[source][gene_id].get(legacy_field))
                new_value = _number(new[gene_id].get(candidate_field))
                equal = (math.isnan(old_value) and math.isnan(new_value)) or math.isclose(
                    old_value, new_value, rel_tol=rtol, abs_tol=atol
                )
                mismatches += not equal
                if math.isfinite(old_value) and math.isfinite(new_value):
                    maximum = max(maximum, abs(old_value - new_value))
            field_results[candidate_field] = {
                "absolute_tolerance": atol,
                "relative_tolerance": rtol,
                "mismatched_genes": mismatches,
                "maximum_absolute_difference": maximum,
            }
            passed = passed and mismatches == 0
        decisions = []
        for gene_id in sorted(shared):
            old_decision = _decision(
                _number(raw[gene_id].get("padj")),
                _number(shrunken[gene_id].get("log2FoldChange")),
                fdr,
                abs_lfc,
            )
            new_decision = _decision(
                _number(new[gene_id].get("adjusted_p_value")),
                _number(new[gene_id].get("log2_fold_change")),
                fdr,
                abs_lfc,
            )
            if old_decision != new_decision:
                decisions.append({"gene_id": gene_id, "reference": old_decision, "candidate": new_decision})
        passed = passed and not decisions
        contrasts[contrast_id] = {
            "passed": passed,
            "reference_genes": len(reference_genes),
            "candidate_genes": len(new),
            "missing_candidate_genes": missing,
            "extra_candidate_genes": extra,
            "fields": field_results,
            "decision_comparison": {
                "fdr": fdr,
                "absolute_log2fc": abs_lfc,
                "mismatched_genes": len(decisions),
                "examples": decisions[:20],
            },
        }
        all_passed = all_passed and passed
    return {"passed": all_passed, "contrasts": contrasts}


def _verify_figure_contract(candidate: Path) -> dict[str, Any]:
    invalid = []
    missing_pairs = []
    for pdf in sorted(candidate.rglob("*.pdf")):
        if not pdf.with_suffix(".png").is_file():
            missing_pairs.append(str(pdf.relative_to(candidate)))
        try:
            reader = PdfReader(str(pdf))
            if not reader.pages or float(reader.pages[0].mediabox.width) <= 0 or float(reader.pages[0].mediabox.height) <= 0:
                invalid.append(str(pdf.relative_to(candidate)))
        except Exception as error:  # malformed PDFs must be reported, not abort verification
            invalid.append(f"{pdf.relative_to(candidate)}: {error}")
    for png in sorted(candidate.rglob("*.png")):
        if not png.with_suffix(".pdf").is_file():
            missing_pairs.append(str(png.relative_to(candidate)))
        try:
            with Image.open(png) as image:
                image.verify()
        except Exception as error:
            invalid.append(f"{png.relative_to(candidate)}: {error}")
    front_door = [candidate / "figures" / "index.json", candidate / "tables" / "index.json"]
    missing_indexes = [str(path.relative_to(candidate)) for path in front_door if not path.is_file()]
    return {
        "passed": not invalid and not missing_pairs and not missing_indexes,
        "invalid_artifacts": invalid,
        "missing_pdf_png_pairs": missing_pairs,
        "missing_front_door_indexes": missing_indexes,
    }


def verify_project(
    project: ResolvedProject,
    reference: str | Path,
    candidate: str | Path | None = None,
    *,
    atol: float = 1e-8,
    rtol: float = 1e-6,
    scope: str = "all",
) -> VerificationResult:
    """Compare v2 results with either a v2 run or the established legacy layout."""
    reference_root = Path(reference).expanduser().resolve()
    candidate_root = Path(candidate).expanduser().resolve() if candidate else project.result_root
    if (reference_root / ".cache" / "counts" / "gene_counts.tsv").is_file():
        counts = _verify_legacy_counts(project, reference_root, candidate_root)
        report = {
            "schema_version": 2,
            "mode": "legacy_migration",
            "scope": scope,
            "project": project.project_id,
            "analysis_set": project.analysis_set,
            "reference": str(reference_root),
            "candidate": str(candidate_root),
            "counts": counts,
        }
        gates = [counts["passed"]]
        if scope in {"core", "all"}:
            report["de"] = _verify_legacy_de(project, reference_root, candidate_root, rtol=rtol)
            gates.append(report["de"]["passed"])
        if scope == "all":
            report["figure_contract"] = _verify_figure_contract(candidate_root)
            gates.append(report["figure_contract"]["passed"])
        report["passed"] = all(gates)
        return VerificationResult(bool(report["passed"]), report)
    return verify_runs(reference_root, candidate_root, atol=atol, rtol=rtol)


def write_verification(result: VerificationResult, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result.report, indent=2) + "\n", encoding="utf-8")
