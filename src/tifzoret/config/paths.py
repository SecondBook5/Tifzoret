"""Path resolution and TSV reading utilities."""

from __future__ import annotations

import csv
import os
import re
from pathlib import Path

from .types import ProjectValidationError


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a tab-separated file into its header list and list of row dicts."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def _resolve(base: Path, value: str) -> Path:
    """Resolve a configured path against ``base``, expanding ``~`` and ``$VAR``.

    Relative paths resolve from the configuration directory. Any environment
    variable that stays unexpanded (i.e. is unset) raises
    :class:`ProjectValidationError` naming the missing variable, so a run never
    silently proceeds with a wrong or empty path.
    """
    expanded = os.path.expandvars(os.path.expanduser(value))
    unresolved = re.findall(r"\$\{([^}]+)\}|\$([A-Za-z_][A-Za-z0-9_]*)", expanded)
    if unresolved:
        names = sorted({left or right for left, right in unresolved})
        raise ProjectValidationError(
            f"required environment variable(s) are not set: {', '.join(names)}"
        )
    candidate = Path(expanded)
    return (candidate if candidate.is_absolute() else base / candidate).resolve()
