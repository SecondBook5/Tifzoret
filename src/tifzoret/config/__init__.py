"""Strict project loading and cross-file validation for Tifzoret."""

from .types import ProjectValidationError, ResolvedProject
from .paths import _read_tsv, _resolve
from .normalize import (
    ALL_MODULES, OPT_IN_MODULES, PROFILE_MODULES,
    migrate_v1_mapping, normalize_config, resolve_modules,
)
from .validate import (
    _load_schema_document, _schema, _validate_document,
    report_json, validation_report,
)
from .presets import _deconvolution_preset_path, deconvolution_presets
from .load import _load_companion, _resolve_bams, _selected_samples, load_project
from .collection import ResolvedCollection, collection_report, load_collection

__all__ = [
    "ProjectValidationError", "ResolvedProject", "_resolve", "load_project",
    "validation_report", "report_json", "normalize_config", "resolve_modules",
    "deconvolution_presets", "ResolvedCollection", "load_collection", "collection_report",
    "migrate_v1_mapping", "ALL_MODULES", "OPT_IN_MODULES",
]
