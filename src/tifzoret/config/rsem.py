"""Validate gene-level RSEM inputs without transforming estimated counts."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from .types import ProjectValidationError, _resolve


def validate_rsem_files(base: Path, samples: list[dict[str, str]]) -> tuple[tuple[Path, ...], list[str]]:
    """Check one gene-result file per selected biological sample.

    Paths resolve relative to project.yaml. Require unique files, unique ordered
    gene identifiers shared across files, and finite nonnegative estimates.
    Zero effective lengths are allowed here because tximport owns their handling.
    Raise ProjectValidationError on malformed input instead of rounding values.
    """
    # Collect paths in the canonical sample order.
    paths: list[Path] = []
    # Preserve the first file's feature order for cross-sample validation.
    canonical: list[str] | None = None
    # Validate every selected sample rather than inferring filenames.
    for sample in samples:
        # Require an explicit gene-level RSEM path.
        value = sample.get("rsem", "").strip()
        # Fail clearly when a sample lacks its quantification file.
        if not value:
            raise ProjectValidationError(f"sample {sample['sample_id']}: requires an rsem path")
        # Resolve paths with the same rules as the other input adapters.
        path = _resolve(base, value)
        # Reject missing files and reuse of another sample's quantification.
        if not path.is_file() or path in paths:
            raise ProjectValidationError(f"missing or duplicate RSEM file: {path}")
        # Read one file without retaining all numeric rows in memory.
        with path.open(newline="", encoding="utf-8") as handle:
            # Parse RSEM's tab-delimited gene-result format.
            reader = csv.DictReader(handle, delimiter="\t")
            # Require the fields consumed by tximport's RSEM gene adapter.
            required = {"gene_id", "length", "effective_length", "expected_count", "TPM"}
            # Reject transcript-level or malformed tables explicitly.
            if not required.issubset(reader.fieldnames or []):
                raise ProjectValidationError(f"RSEM gene-result columns missing: {path}")
            # Accumulate only the ordered feature identifiers.
            genes: list[str] = []
            # Validate numeric estimates before handing the files to R.
            for number, row in enumerate(reader, start=2):
                # Preserve the gene identifier without guessing a mapping.
                genes.append(row.get("gene_id", "").strip())
                # Require finite nonnegative quantification fields.
                for field in ("length", "effective_length", "expected_count", "TPM"):
                    # Convert one numeric field with a contextual error.
                    try:
                        numeric = float(row[field])
                    except (ValueError, TypeError, KeyError) as error:
                        raise ProjectValidationError(f"{path}:{number}: invalid {field}") from error
                    # Reject impossible or non-finite estimates.
                    if not math.isfinite(numeric) or numeric < 0:
                        raise ProjectValidationError(f"{path}:{number}: invalid {field}")
        # Require a nonempty unique gene universe.
        if not genes or any(not gene for gene in genes) or len(genes) != len(set(genes)):
            raise ProjectValidationError(f"empty or duplicate RSEM gene identifiers: {path}")
        # Require all files to share the exact feature order.
        if canonical is not None and genes != canonical:
            raise ProjectValidationError(f"RSEM gene identity/order differs: {path}")
        # Retain the validated universe and source path.
        canonical = genes
        paths.append(path)
    # Return canonical paths and identifiers for provenance and annotation checks.
    return tuple(paths), canonical or []
