"""Configuration normalization and module resolution."""

from __future__ import annotations

import copy
from typing import Any

from .types import ProjectValidationError


PROFILE_MODULES: dict[str, frozenset[str]] = {
    "standard": frozenset({"qc", "de", "pathways", "ontology", "report"}),
    "publication": frozenset(
        {
            "qc", "de", "pathways", "ontology", "composition", "regulators",
            "networks", "hypotheses", "publication", "report",
        }
    ),
    "full": frozenset(
        {
            "qc", "de", "pathways", "ontology", "composition", "regulators",
            "networks", "hypotheses", "publication", "report", "sva", "wgcna",
            "mediation", "multilayer",
        }
    ),
}
# Opt-in modules belong to no profile — they run only when explicitly toggled on
# via analysis.modules — so union them in rather than deriving ALL_MODULES from
# the profiles alone. batch: corrected PCA/distance views; de_confirm: edgeR
# confirmatory DE; deconvolution: signature-matrix cell fractions; curvature:
# Ollivier-Ricci on the WGCNA co-expression graph; consensus: cross-contrast
# direction consensus; spia: signed pathway-topology impact; variance_partition:
# variance explained per design covariate; enrichment_map: term-similarity map;
# factorial: interaction views (effect-vs-effect, profile, four-group expression).
OPT_IN_MODULES = frozenset({
    "batch", "de_confirm", "deconvolution", "curvature",
    "consensus", "spia", "variance_partition", "enrichment_map", "factorial",
})
ALL_MODULES = frozenset().union(*PROFILE_MODULES.values()) | OPT_IN_MODULES


def migrate_v1_mapping(config: dict[str, Any]) -> dict[str, Any]:
    """Convert the development v1 mapping into the public v2 contract."""
    if config.get("version") != 1:
        raise ProjectValidationError("migrate_v1_mapping requires version: 1")
    inputs = copy.deepcopy(config["inputs"])
    gtf = inputs.get("gtf")
    legacy_modules = config.get("modules", {})
    module_overrides = {name: False for name in ALL_MODULES}
    module_overrides.update({name: bool(value) for name, value in legacy_modules.items()})
    migrated: dict[str, Any] = {
        "version": 2,
        "project": copy.deepcopy(config["project"]),
        "species": {
            "provider": "custom",
            "scientific_name": "unspecified",
            "taxonomy_id": None,
        },
        "reference": {
            "genome_build": "unspecified",
            "annotation_release": None,
        },
        "inputs": inputs,
        "analysis": {
            "design": config["design"]["formula"],
            "contrasts": config["contrasts"],
            "profile": "standard",
            "random_seed": config.get("figures", {}).get("pathways", {}).get("seed", 1),
            # v1 predates profiles; explicitly disable newly introduced
            # modules so migration preserves the workflow that actually ran.
            "modules": module_overrides,
        },
        "resources": {
            "cache": "~/.cache/tifzoret/resources",
            "offline": False,
            "refresh": False,
            "gene_sets": copy.deepcopy(config["gene_sets"]),
        },
        "figures": copy.deepcopy(config["figures"]),
        "output": copy.deepcopy(config["output"]),
    }
    if "counting" in config:
        migrated["counting"] = copy.deepcopy(config["counting"])
    if gtf is not None:
        migrated["inputs"]["gtf"] = gtf
    return migrated


def normalize_config(config: dict[str, Any]) -> dict[str, Any]:
    """Return one v2 configuration model for either supported file version."""
    version = config.get("version")
    if version == 1:
        return migrate_v1_mapping(config)
    if version == 2:
        return copy.deepcopy(config)
    raise ProjectValidationError(f"unsupported project configuration version: {version!r}")


def resolve_modules(config: dict[str, Any]) -> dict[str, bool]:
    """Return the enabled/disabled state of every module for this project.

    Starts from the profile's default module set (``standard``/``publication``/
    ``full``) and applies any explicit ``analysis.modules`` overrides on top, so
    the result is the exact set of stages the workflow will run.
    """
    analysis = config["analysis"]
    enabled = {name: name in PROFILE_MODULES[analysis["profile"]] for name in ALL_MODULES}
    enabled.update(analysis.get("modules", {}))
    return enabled
