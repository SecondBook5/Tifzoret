"""Strict project loading and cross-file validation for Tifzoret."""

# The underscore-prefixed names below are re-exported on purpose, so existing callers of `tifzoret.config` keep working.
from .types import ProjectValidationError, ResolvedProject, _read_tsv, _resolve  # noqa: F401
from .normalize import ALL_MODULES, OPT_IN_MODULES, PROFILE_MODULES  # noqa: F401
from .normalize import migrate_v1_mapping, normalize_config, resolve_modules
from .validate import _load_schema_document, _schema, _validate_document  # noqa: F401
from .validate import report_json, validation_report
from .load import _deconvolution_preset_path, deconvolution_presets  # noqa: F401
from .load import _load_companion, _resolve_bams, _selected_samples, load_project  # noqa: F401
from .collection import ResolvedCollection, collection_report, load_collection
from .estimands import (
    CELL_KEY_SEPARATOR,
    Estimand,
    EstimandSyntaxError,
    Family,
    TermTest,
    compile_project_families,
)

__all__ = [
    "ProjectValidationError", "ResolvedProject", "_resolve", "load_project",
    "validation_report", "report_json", "normalize_config", "resolve_modules",
    "deconvolution_presets", "ResolvedCollection", "load_collection", "collection_report",
    "migrate_v1_mapping", "ALL_MODULES", "OPT_IN_MODULES",
    "CELL_KEY_SEPARATOR", "Estimand", "EstimandSyntaxError", "Family", "TermTest",
    "compile_project_families",
]
