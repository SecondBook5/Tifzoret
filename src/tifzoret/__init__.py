"""Tifzoret: reproducible downstream bulk RNA-seq analysis."""

from tifzoret.config import (
    ProjectValidationError,
    ResolvedProject,
    load_project,
    validation_report,
)
from tifzoret.figures import (
    PANEL_REGISTRY,
    build_gallery,
    normalized_gene_panels,
    resolve_panel,
)

__version__ = "0.1.0.dev0"

__all__ = [
    "ProjectValidationError", "ResolvedProject", "load_project", "validation_report",
    "PANEL_REGISTRY", "build_gallery", "normalized_gene_panels", "resolve_panel",
    "__version__",
]
