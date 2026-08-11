#!/usr/bin/env python3
"""Write a checksummed, auditable release manifest for one completed run."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from bulk_rna_frame.config import load_project  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(
    path: Path,
    relative_to: Path | None = None,
    known: dict[str, object] | None = None,
) -> dict[str, object]:
    stat = path.stat()
    reusable = (
        known is not None
        and known.get("sha256")
        and known.get("bytes") == stat.st_size
        and (known.get("mtime_ns") is None or known.get("mtime_ns") == stat.st_mtime_ns)
    )
    return {
        "path": str(path.relative_to(relative_to) if relative_to else path),
        "bytes": stat.st_size,
        "sha256": str(known["sha256"]) if reusable else sha256(path),
    }


def prepared_checksums(results: Path) -> dict[str, dict[str, object]]:
    """Index checksums already established by canonical input preparation."""
    path = results / "inputs" / "input_manifest.json"
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}

    indexed: dict[str, dict[str, object]] = {}

    def visit(value: object) -> None:
        if isinstance(value, dict):
            candidate = value.get("path")
            if candidate and value.get("sha256") and value.get("bytes") is not None:
                candidate_path = Path(str(candidate))
                if candidate_path.is_absolute():
                    indexed[str(candidate_path.resolve())] = value
            for child in value.values():
                visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)

    visit(manifest)
    return indexed


def command_output(command: list[str]) -> str:
    try:
        result = subprocess.run(
            command,
            text=True,
            encoding="utf-8",
            errors="replace",
            capture_output=True,
            check=False,
        )
    except OSError:
        return "unavailable"
    lines = [line for line in (result.stdout + result.stderr).splitlines() if line.strip()]
    return lines[0] if lines else "unavailable"


def repository_revision(start: Path) -> dict[str, object]:
    revision = command_output(["git", "-C", str(start), "rev-parse", "HEAD"])
    status = command_output(["git", "-C", str(start), "status", "--porcelain"])
    return {
        "revision": revision,
        "dirty": status not in {"", "unavailable"},
    }


def collect_warnings(results: Path) -> list[dict[str, str]]:
    warnings: list[dict[str, str]] = []
    for path in sorted(results.rglob("*.json")):
        if path.name == "manifest.json":
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        values = data.get("warnings", []) if isinstance(data, dict) else []
        if isinstance(values, str):
            values = [values]
        for value in values:
            warnings.append({"source": str(path.relative_to(results)), "message": str(value)})
    return warnings


def collect_resource_receipts(results: Path) -> list[dict[str, Any]]:
    receipts: list[dict[str, Any]] = []
    cache = results / ".cache" / "resources"
    if not cache.is_dir():
        return receipts
    for path in sorted(cache.rglob("*.json")):
        try:
            receipt = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if isinstance(receipt, dict):
            receipt = {**receipt, "receipt": str(path.relative_to(results))}
            receipts.append(receipt)
    return receipts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--results", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    project = load_project(args.project_config)
    results = Path(args.results).resolve()
    output = Path(args.output).resolve()
    input_paths = [
        project.config_path,
        project.samples,
        project.contrasts,
        project.gmt,
        *project.source_files,
    ]
    input_paths.extend(
        path
        for path in (
            project.hypotheses,
            project.hypothesis_panels,
            project.figure_recipe,
            project.cell_state_signatures,
            project.regulon_edges,
        )
        if path is not None
    )
    unique_inputs = list(dict.fromkeys(path.resolve() for path in input_paths))
    known_checksums = prepared_checksums(results)
    result_paths = sorted(
        path
        for path in results.rglob("*")
        if path.is_file() and path != output and not path.name.endswith(".log")
    )
    manifest = {
        "schema_version": 2,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": project.config["project"],
        "analysis_set": project.analysis_set,
        "profile": project.config["analysis"]["profile"],
        "modules": list(project.modules),
        "species": project.config["species"],
        "reference": project.config["reference"],
        "configuration": project.config,
        "contrast_semantics": "all signed effects are numerator minus denominator",
        "contrasts": [
            {
                "contrast_id": row["contrast_id"],
                "factor": row["factor"],
                "numerator": row["numerator"],
                "denominator": row["denominator"],
            }
            for row in project.contrast_rows
        ],
        "random_seeds": {
            "global": project.config["analysis"].get("random_seed", 1),
            "pathways": project.config["figures"]["pathways"]["seed"],
        },
        "platform": platform.platform(),
        "repository": repository_revision(SOURCE_ROOT),
        "environment": {
            "conda_prefix": os.environ.get("CONDA_PREFIX"),
            "container": os.environ.get("APPTAINER_CONTAINER") or os.environ.get("SINGULARITY_CONTAINER"),
        },
        "tools": {
            "python": platform.python_version(),
            "snakemake": command_output(["snakemake", "--version"]),
            "R": command_output(["R", "--version"]),
            "samtools": command_output(["samtools", "--version"]),
            "featureCounts": command_output(["featureCounts", "-v"]),
        },
        "resource_snapshots": collect_resource_receipts(results),
        "warnings": collect_warnings(results),
        "inputs": [record(path, known=known_checksums.get(str(path))) for path in unique_inputs],
        "results": [record(path, results) for path in result_paths],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output} with {len(result_paths)} result files")


if __name__ == "__main__":
    main()
