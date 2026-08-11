"""Strict project loading and cross-file validation for BulkRNAFrame."""

from __future__ import annotations

import csv
import json
import os
import re
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Any

import yaml
from jsonschema import Draft202012Validator


class ProjectValidationError(ValueError):
    """Raised when a project cannot satisfy the workflow contract."""


@dataclass(frozen=True)
class ResolvedProject:
    config_path: Path
    config: dict[str, Any]
    source_kind: str
    samples: Path
    contrasts: Path
    gmt: Path
    output_root: Path
    sample_rows: tuple[dict[str, str], ...]
    contrast_rows: tuple[dict[str, str], ...]
    counts: Path | None = None
    annotation: Path | None = None
    gtf: Path | None = None
    source_root: Path | None = None
    bam_paths: tuple[Path, ...] = ()

    @property
    def project_id(self) -> str:
        return str(self.config["project"]["id"])

    @property
    def source_files(self) -> tuple[Path, ...]:
        if self.source_kind == "counts":
            return tuple(path for path in (self.counts, self.annotation) if path is not None)
        return (*self.bam_paths, *((self.gtf,) if self.gtf is not None else ()))


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def _resolve(base: Path, value: str) -> Path:
    expanded = os.path.expandvars(os.path.expanduser(value))
    candidate = Path(expanded)
    return (candidate if candidate.is_absolute() else base / candidate).resolve()


def _schema() -> dict[str, Any]:
    path = resources.files("bulk_rna_frame").joinpath("schemas/project.schema.yaml")
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _selected_samples(
    samples: list[dict[str, str]], sample_header: list[str], analysis_set: str | None
) -> list[dict[str, str]]:
    if analysis_set is None:
        return samples
    if "analysis_set" not in sample_header:
        raise ProjectValidationError(
            "inputs.analysis_set requires an analysis_set column in samples.tsv"
        )
    selected = [
        row
        for row in samples
        if analysis_set in {value.strip() for value in row["analysis_set"].split(",") if value.strip()}
    ]
    if not selected:
        raise ProjectValidationError(f"analysis_set {analysis_set!r} selects no samples")
    return selected


def _resolve_bams(
    *,
    kind: str,
    inputs: dict[str, Any],
    base: Path,
    sample_header: list[str],
    samples: list[dict[str, str]],
    errors: list[str],
) -> tuple[Path | None, tuple[Path, ...]]:
    root_key = "bam_root" if kind == "bam" else "root"
    root = _resolve(base, inputs[root_key])
    if not root.is_dir():
        errors.append(f"{root_key} directory does not exist: {root}")

    resolved: list[Path] = []
    if kind == "bam" and "bam" not in sample_header:
        errors.append("BAM input requires a bam column in samples.tsv")
        return root, ()

    for row in samples:
        sample_id = row.get("sample_id", "<missing>")
        if kind == "bam":
            relative = row.get("bam", "").strip()
            if not relative:
                errors.append(f"sample {sample_id}: bam must be non-empty")
                continue
        else:
            try:
                relative = str(inputs["bam_pattern"]).format_map(row)
            except KeyError as error:
                errors.append(
                    f"sample {sample_id}: bam_pattern references missing metadata column {error.args[0]!r}"
                )
                continue
        path = _resolve(root, relative)
        try:
            path.relative_to(root)
        except ValueError:
            errors.append(f"sample {sample_id}: BAM resolves outside the configured root: {path}")
            continue
        if not path.is_file():
            errors.append(f"sample {sample_id}: BAM does not exist: {path}")
        resolved.append(path)
    return root, tuple(resolved)


def load_project(config_path: str | Path) -> ResolvedProject:
    path = Path(config_path).expanduser().resolve()
    if not path.is_file():
        raise ProjectValidationError(f"Project configuration does not exist: {path}")
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise ProjectValidationError("Project configuration must be a YAML mapping.")

    schema_errors = sorted(
        Draft202012Validator(_schema()).iter_errors(config),
        key=lambda error: tuple(error.absolute_path),
    )
    errors = [
        f"schema {'.'.join(map(str, error.absolute_path)) or '<root>'}: {error.message}"
        for error in schema_errors
    ]
    if errors:
        raise ProjectValidationError("\n".join(errors))

    base = path.parent
    inputs = config["inputs"]
    kind = str(inputs["kind"])
    samples_path = _resolve(base, inputs["samples"])
    contrasts_path = _resolve(base, config["contrasts"])
    gmt_path = _resolve(base, config["gene_sets"]["gmt"])
    for name, resolved in (
        ("samples", samples_path),
        ("contrasts", contrasts_path),
        ("gmt", gmt_path),
    ):
        if not resolved.is_file():
            errors.append(f"{name} file does not exist: {resolved}")
    if errors:
        raise ProjectValidationError("\n".join(errors))

    sample_header, all_samples = _read_tsv(samples_path)
    if "sample_id" not in sample_header:
        errors.append("samples.tsv requires a sample_id column")
    if not all_samples:
        errors.append("samples.tsv contains no samples")
    if errors:
        raise ProjectValidationError("\n".join(errors))
    try:
        samples = _selected_samples(all_samples, sample_header, inputs.get("analysis_set"))
    except ProjectValidationError as error:
        errors.append(str(error))
        samples = []

    sample_ids = [row.get("sample_id", "").strip() for row in samples]
    if any(not sample_id for sample_id in sample_ids):
        errors.append("sample_id values must be non-empty")
    if len(sample_ids) != len(set(sample_ids)):
        errors.append("selected sample_id values must be unique")

    counts_path: Path | None = None
    annotation_path: Path | None = None
    gtf_path: Path | None = None
    source_root: Path | None = None
    bam_paths: tuple[Path, ...] = ()
    gene_ids: list[str] = []

    if kind == "counts":
        counts_path = _resolve(base, inputs["counts"])
        annotation_path = _resolve(base, inputs["annotation"])
        for name, resolved in (("counts", counts_path), ("annotation", annotation_path)):
            if not resolved.is_file():
                errors.append(f"{name} file does not exist: {resolved}")
        if not errors:
            count_header, count_rows = _read_tsv(counts_path)
            if not count_header or count_header[0] != "gene_id":
                errors.append("counts.tsv must begin with a gene_id column")
            count_samples = count_header[1:]
            missing_selected = sorted(set(sample_ids) - set(count_samples))
            unknown_counts = sorted(set(count_samples) - {row.get("sample_id", "") for row in all_samples})
            if missing_selected or unknown_counts:
                errors.append(
                    "count-matrix columns must contain every selected sample and no undeclared samples; "
                    f"missing_selected={missing_selected}, unknown_counts={unknown_counts}"
                )
            if inputs.get("analysis_set") is None and set(count_samples) != set(sample_ids):
                errors.append("count-matrix sample columns must exactly match sample_id values")
            if not count_rows:
                errors.append("counts.tsv contains no genes")
            for row_number, row in enumerate(count_rows, start=2):
                gene_id = row.get("gene_id", "").strip()
                gene_ids.append(gene_id)
                for sample_id in count_samples:
                    raw = row.get(sample_id, "")
                    try:
                        value = int(raw)
                    except (TypeError, ValueError):
                        errors.append(f"counts.tsv row {row_number}, {sample_id}: expected an integer")
                        continue
                    if value < 0:
                        errors.append(f"counts.tsv row {row_number}, {sample_id}: counts cannot be negative")
            if any(not gene_id for gene_id in gene_ids):
                errors.append("gene_id values must be non-empty")
            if len(gene_ids) != len(set(gene_ids)):
                errors.append("gene_id values must be unique")

            annotation_header, annotations = _read_tsv(annotation_path)
            if not {"gene_id", "gene_symbol"}.issubset(annotation_header):
                errors.append("annotation.tsv requires gene_id and gene_symbol columns")
            annotation_ids = {row.get("gene_id", "") for row in annotations}
            missing_annotation = sorted(set(gene_ids) - annotation_ids)
            if missing_annotation:
                errors.append(
                    f"annotation.tsv is missing {len(missing_annotation)} count gene_id values "
                    f"(first: {missing_annotation[:5]})"
                )
    else:
        gtf_path = _resolve(base, inputs["gtf"])
        if not gtf_path.is_file():
            errors.append(f"gtf file does not exist: {gtf_path}")
        counting = config["counting"]
        if counting["strandedness"] == "infer" and not {1, 2}.issubset(
            counting["strand_test_modes"]
        ):
            errors.append(
                "counting.strand_test_modes must include 1 and 2 when strandedness is infer"
            )
        if not counting["paired_end"] and any(
            counting[key]
            for key in (
                "count_read_pairs",
                "require_both_ends_aligned",
                "exclude_chimeric_fragments",
            )
        ):
            errors.append(
                "count_read_pairs, require_both_ends_aligned, and "
                "exclude_chimeric_fragments require counting.paired_end: true"
            )
        source_root, bam_paths = _resolve_bams(
            kind=kind,
            inputs=inputs,
            base=base,
            sample_header=sample_header,
            samples=samples,
            errors=errors,
        )

    contrast_header, contrasts = _read_tsv(contrasts_path)
    required_contrast = {"contrast_id", "factor", "numerator", "denominator"}
    if not required_contrast.issubset(contrast_header):
        errors.append("contrasts.tsv requires contrast_id, factor, numerator, and denominator columns")
    if not contrasts:
        errors.append("contrasts.tsv contains no contrasts")
    contrast_ids = [row.get("contrast_id", "").strip() for row in contrasts]
    if len(contrast_ids) != len(set(contrast_ids)):
        errors.append("contrast_id values must be unique")
    for row in contrasts:
        contrast_id = row.get("contrast_id", "<missing>")
        if not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", contrast_id):
            errors.append(f"invalid contrast_id {contrast_id!r}")
        factor = row.get("factor", "")
        if factor not in sample_header:
            errors.append(f"contrast {contrast_id}: factor {factor!r} is absent from samples.tsv")
            continue
        levels = {sample.get(factor, "") for sample in samples}
        numerator = row.get("numerator", "")
        denominator = row.get("denominator", "")
        if numerator == denominator:
            errors.append(f"contrast {contrast_id}: numerator and denominator must differ")
        for level in (numerator, denominator):
            if level not in levels:
                errors.append(f"contrast {contrast_id}: level {level!r} is absent from {factor}")
        if not re.search(rf"\b{re.escape(factor)}\b", config["design"]["formula"]):
            errors.append(f"contrast {contrast_id}: factor {factor!r} is absent from design formula")

    group = config["figures"]["group"]
    if group not in sample_header:
        errors.append(f"figures.group {group!r} is absent from samples.tsv")
    else:
        groups = {row[group] for row in samples}
        palette_groups = set(config["figures"]["palette"])
        missing_colors = sorted(groups - palette_groups)
        if missing_colors:
            errors.append(f"figures.palette has no colors for: {missing_colors}")

    if config["gene_sets"]["min_size"] > config["gene_sets"]["max_size"]:
        errors.append("gene_sets.min_size cannot exceed gene_sets.max_size")
    modules = config["modules"]
    if not any(modules.values()):
        errors.append("at least one analysis module must be enabled")
    if modules["pathways"] and not modules["de"]:
        errors.append("modules.pathways requires modules.de")
    if modules["pathways"] and not modules["qc"]:
        errors.append("modules.pathways requires modules.qc because it consumes VST expression")
    valid_gmt_lines = 0
    with gmt_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                errors.append(
                    f"GMT line {line_number} must contain a name, description, and at least two genes"
                )
            else:
                valid_gmt_lines += 1
    if valid_gmt_lines == 0:
        errors.append("GMT file contains no valid gene sets")

    if errors:
        raise ProjectValidationError("\n".join(errors))

    return ResolvedProject(
        config_path=path,
        config=config,
        source_kind=kind,
        counts=counts_path,
        samples=samples_path,
        annotation=annotation_path,
        contrasts=contrasts_path,
        gmt=gmt_path,
        gtf=gtf_path,
        source_root=source_root,
        bam_paths=bam_paths,
        output_root=_resolve(base, config["output"]["root"]),
        sample_rows=tuple(samples),
        contrast_rows=tuple(contrasts),
    )


def validation_report(project: ResolvedProject) -> dict[str, Any]:
    report: dict[str, Any] = {
        "status": "ok",
        "project_id": project.project_id,
        "config": str(project.config_path),
        "samples": len(project.sample_rows),
        "contrasts": [row["contrast_id"] for row in project.contrast_rows],
        "input_kind": project.source_kind,
        "output": str(project.output_root / project.project_id),
    }
    if project.bam_paths:
        report["bams"] = len(project.bam_paths)
        report["source_root"] = str(project.source_root)
    return report


def report_json(project: ResolvedProject) -> str:
    return json.dumps(validation_report(project), indent=2) + "\n"
