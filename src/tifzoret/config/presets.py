"""Deconvolution preset utilities."""

from __future__ import annotations

from importlib import resources
from pathlib import Path

from .types import ProjectValidationError


def deconvolution_presets() -> tuple[str, ...]:
    """Return the names of the cell-type reference matrices shipped with the package.

    A preset is any ``<name>.tsv`` under the packaged ``data/deconvolution``
    directory (documented in that directory's ``PROVENANCE.md``). Names are the
    file stems, sorted, so the set is discovered from what is installed rather
    than hard-coded here.
    """
    directory = resources.files("tifzoret").joinpath("data/deconvolution")
    if not directory.is_dir():
        return ()
    return tuple(sorted(entry.name[:-4] for entry in directory.iterdir() if entry.name.endswith(".tsv")))


def _deconvolution_preset_path(name: str) -> Path:
    """Resolve a named deconvolution preset to its packaged signature-matrix path.

    Raises :class:`ProjectValidationError` (naming the available presets) when the
    requested preset is not installed, so a mistyped name fails at validation
    rather than silently deconvolving against nothing.
    """
    entry = resources.files("tifzoret").joinpath(f"data/deconvolution/{name}.tsv")
    if not entry.is_file():
        available = ", ".join(deconvolution_presets()) or "(none installed)"
        raise ProjectValidationError(
            f"unknown resources.deconvolution_preset {name!r}; available presets: {available}"
        )
    return Path(str(entry))
