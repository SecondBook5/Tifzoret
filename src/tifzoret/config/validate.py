"""Schema validation and reporting utilities."""

from __future__ import annotations

import json
from importlib import resources
from pathlib import Path
from typing import Any

import yaml
from jsonschema import Draft202012Validator

from .types import ResolvedProject


def _schema(name: str = "project") -> dict[str, Any]:
    path = resources.files("tifzoret").joinpath(f"schemas/{name}.schema.yaml")
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _validate_document(document: dict[str, Any], schema_name: str) -> list[str]:
    validation_errors = sorted(
        Draft202012Validator(_schema(schema_name)).iter_errors(document),
        key=lambda error: tuple(error.absolute_path),
    )
    return [
        f"{schema_name} schema {'.'.join(map(str, error.absolute_path)) or '<root>'}: {error.message}"
        for error in validation_errors
    ]


def _load_schema_document(path: Path, schema_name: str) -> tuple[dict[str, Any] | None, list[str]]:
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        return None, [f"{schema_name} could not be read: {error}"]
    if not isinstance(document, dict):
        return None, [f"{schema_name} must be a YAML mapping: {path}"]
    return document, _validate_document(document, schema_name)


def validation_report(project: ResolvedProject) -> dict[str, Any]:
    """Summarize a validated project as a JSON-serializable dictionary.

    Reports the identity, sample/contrast counts, input kind, active profile and
    modules, species, and resolved output root — the human-facing confirmation
    that ``validate`` prints.
    """
    report: dict[str, Any] = {
        "status": "ok",
        "project_id": project.project_id,
        "config": str(project.config_path),
        "samples": len(project.sample_rows),
        "contrasts": [row["contrast_id"] for row in project.contrast_rows],
        "input_kind": project.source_kind,
        "analysis_set": project.analysis_set,
        "profile": project.config["analysis"]["profile"],
        "strict": project.strict,
        "modules": list(project.modules),
        "species": project.config["species"],
        "output": str(project.result_root),
    }
    if project.bam_paths:
        report["bams"] = len(project.bam_paths)
        report["source_root"] = str(project.source_root)
    return report


def report_json(project: ResolvedProject) -> str:
    """Render :func:`validation_report` as an indented JSON string with a trailing newline."""
    return json.dumps(validation_report(project), indent=2) + "\n"
