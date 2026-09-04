"""Core types for project configuration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    from .estimands import Family


class ProjectValidationError(ValueError):
    """Raised when a project cannot satisfy the workflow contract."""


@dataclass(frozen=True)
class ResolvedProject:
    """A fully validated project with every path resolved and input checked.

    Produced by :func:`load_project` once the configuration, the tabular inputs
    (samples, contrasts, counts/BAMs, GMT), and any companion documents all
    satisfy the workflow contract. Immutable, so it can be passed to the CLI,
    the figure layer, and the workflow rules as a single source of truth.
    """

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
    archive: Path | None = None
    analysis_set: str = "all"
    modules: tuple[str, ...] = ()
    hypotheses: Path | None = None
    hypothesis_panels: Path | None = None
    figure_recipe: Path | None = None
    hypothesis_config: dict[str, Any] | None = None
    panel_config: dict[str, Any] | None = None
    recipe_config: dict[str, Any] | None = None
    cell_state_signatures: Path | None = None
    regulon_edges: Path | None = None
    binding_prior_edges: Path | None = None
    # Generalized named regulator views (analysis.settings.regulators.views). Each
    # entry is {name, edges (resolved absolute path str), signed (bool),
    # provider_label (str | None)}; empty when the legacy regulon_edges(+binding)
    # path is in use. The first view is the primary (canonical filenames).
    regulator_views: tuple[dict[str, Any], ...] = ()
    deconvolution_signature: Path | None = None
    # Compiled estimand layer. `families` carries the validated Family objects;
    # the *_rows tuples are the flat, TSV-ready projections that
    # materialize_inputs.py stages into inputs/ for the R stages to read, so an
    # expression is never parsed twice (once in Python, once in R).
    families: tuple[Family, ...] = ()
    family_rows: tuple[dict[str, str], ...] = ()
    estimand_rows: tuple[dict[str, str], ...] = ()
    cell_weight_rows: tuple[dict[str, str], ...] = ()
    term_test_rows: tuple[dict[str, str], ...] = ()

    @property
    def project_id(self) -> str:
        """The filesystem-safe project identifier (``project.id``)."""
        return str(self.config["project"]["id"])

    @property
    def strict(self) -> bool:
        """Whether ``execution.strict`` selects publication-execution mode.

        In strict mode any degraded or fallback path (e.g. a non-VIPER
        regulator-activity proxy) is a hard failure rather than a recorded
        warning, and the terminal manifest step fails if any stage recorded a
        warning. Absent/false is the default and behaves identically to the
        historical engine.
        """
        return bool(self.config.get("execution", {}).get("strict", False))

    @property
    def source_files(self) -> tuple[Path, ...]:
        """The primary input files for this boundary (counts+annotation, archive+GTF, or BAMs+GTF)."""
        if self.source_kind == "counts":
            return tuple(path for path in (self.counts, self.annotation) if path is not None)
        if self.source_kind == "archive":
            return tuple(path for path in (self.archive, self.gtf) if path is not None)
        return (*self.bam_paths, *((self.gtf,) if self.gtf is not None else ()))

    @property
    def result_root(self) -> Path:
        """Where this run writes: ``<output.root>/<project.id>/<analysis_set>/``."""
        return self.output_root / self.project_id / self.analysis_set
