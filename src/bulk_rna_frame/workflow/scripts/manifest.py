#!/usr/bin/env python3
"""Write a checksummed release manifest for one completed project."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(path: Path, relative_to: Path | None = None) -> dict[str, object]:
    return {
        "path": str(path.relative_to(relative_to) if relative_to else path),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def command_version(command: list[str]) -> str:
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


def resolve(base: Path, value: str) -> Path:
    path = Path(os.path.expandvars(value)).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--results", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config_path = Path(args.project_config).resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    base = config_path.parent
    results = Path(args.results).resolve()
    output = Path(args.output).resolve()
    input_paths = [
        config_path,
        resolve(base, config["inputs"]["samples"]),
        resolve(base, config["contrasts"]),
        resolve(base, config["gene_sets"]["gmt"]),
    ]
    if config["inputs"]["kind"] == "counts":
        input_paths.extend(
            [
                resolve(base, config["inputs"]["counts"]),
                resolve(base, config["inputs"]["annotation"]),
            ]
        )
    else:
        input_paths.append(resolve(base, config["inputs"]["gtf"]))
    result_paths = sorted(
        path for path in results.rglob("*")
        if path.is_file() and path != output and not path.name.endswith(".log")
    )
    manifest = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": config["project"],
        "configuration": config,
        "contrast_direction": "all signed effects are numerator minus denominator",
        "platform": platform.platform(),
        "tools": {
            "python": platform.python_version(),
            "snakemake": command_version(["snakemake", "--version"]),
            "R": command_version(["R", "--version"]),
            "samtools": command_version(["samtools", "--version"]),
            "featureCounts": command_version(["featureCounts", "-v"]),
        },
        "inputs": [record(path) for path in input_paths],
        "results": [record(path, results) for path in result_paths],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output} with {len(result_paths)} result files")


if __name__ == "__main__":
    main()
