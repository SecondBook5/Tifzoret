"""Strict project loading and cross-file validation for BulkRNAFrame."""

from __future__ import annotations

import csv
import copy
import json
import os
import re
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from pathlib import PurePosixPath
from typing import Any

import yaml
from jsonschema import Draft202012Validator


class ProjectValidationError(ValueError):
    """Raised when a project cannot satisfy the workflow contract."""


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
ALL_MODULES = frozenset().union(*PROFILE_MODULES.values())


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
            "cache": "~/.cache/bulk-rna-frame/resources",
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
    analysis = config["analysis"]
    enabled = {name: name in PROFILE_MODULES[analysis["profile"]] for name in ALL_MODULES}
    enabled.update(analysis.get("modules", {}))
    return enabled


@dataclass(frozen=True)
class ResolvedProject:
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

    @property
    def project_id(self) -> str:
        return str(self.config["project"]["id"])

    @property
    def source_files(self) -> tuple[Path, ...]:
        if self.source_kind == "counts":
            return tuple(path for path in (self.counts, self.annotation) if path is not None)
        if self.source_kind == "archive":
            return tuple(path for path in (self.archive, self.gtf) if path is not None)
        return (*self.bam_paths, *((self.gtf,) if self.gtf is not None else ()))

    @property
    def result_root(self) -> Path:
        return self.output_root / self.project_id / self.analysis_set


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def _resolve(base: Path, value: str) -> Path:
    expanded = os.path.expandvars(os.path.expanduser(value))
    unresolved = re.findall(r"\$\{([^}]+)\}|\$([A-Za-z_][A-Za-z0-9_]*)", expanded)
    if unresolved:
        names = sorted({left or right for left, right in unresolved})
        raise ProjectValidationError(
            f"required environment variable(s) are not set: {', '.join(names)}"
        )
    candidate = Path(expanded)
    return (candidate if candidate.is_absolute() else base / candidate).resolve()


def _schema(name: str = "project") -> dict[str, Any]:
    path = resources.files("bulk_rna_frame").joinpath(f"schemas/{name}.schema.yaml")
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _load_schema_document(path: Path, schema_name: str) -> tuple[dict[str, Any] | None, list[str]]:
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        return None, [f"{schema_name} could not be read: {error}"]
    if not isinstance(document, dict):
        return None, [f"{schema_name} must be a YAML mapping: {path}"]
    validation_errors = sorted(
        Draft202012Validator(_schema(schema_name)).iter_errors(document),
        key=lambda error: tuple(error.absolute_path),
    )
    return document, [
        f"{schema_name} schema {'.'.join(map(str, error.absolute_path)) or '<root>'}: {error.message}"
        for error in validation_errors
    ]


def _selected_samples(
    samples: list[dict[str, str]], sample_header: list[str], analysis_set: str | None
) -> list[dict[str, str]]:
    if analysis_set is None:
        return samples
    if "analysis_set" not in sample_header:
        raise ProjectValidationError(
            "inputs.analysis_set requires an analysis_set column in samples.tsv"
        )
    selected = [
        row
        for row in samples
        if analysis_set in {value.strip() for value in row["analysis_set"].split(",") if value.strip()}
    ]
    if not selected:
        raise ProjectValidationError(f"analysis_set {analysis_set!r} selects no samples")
    return selected


def _resolve_bams(
    *,
    kind: str,
    inputs: dict[str, Any],
    base: Path,
    sample_header: list[str],
    samples: list[dict[str, str]],
    errors: list[str],
) -> tuple[Path | None, tuple[Path, ...]]:
    root_key = "bam_root" if kind == "bam" else "root"
    root = _resolve(base, inputs[root_key])
    if not root.is_dir():
        errors.append(f"{root_key} directory does not exist: {root}")

    resolved: list[Path] = []
    if kind == "bam" and "bam" not in sample_header:
        errors.append("BAM input requires a bam column in samples.tsv")
        return root, ()

    for row in samples:
        sample_id = row.get("sample_id", "<missing>")
        if kind == "bam":
            relative = row.get("bam", "").strip()
            if not relative:
                errors.append(f"sample {sample_id}: bam must be non-empty")
                continue
        else:
            try:
                relative = str(inputs["bam_pattern"]).format_map(row)
            except KeyError as error:
                errors.append(
                    f"sample {sample_id}: bam_pattern references missing metadata column {error.args[0]!r}"
                )
                continue
        path = _resolve(root, relative)
        try:
            path.relative_to(root)
        except ValueError:
            errors.append(f"sample {sample_id}: BAM resolves outside the configured root: {path}")
            continue
        if not path.is_file():
            errors.append(f"sample {sample_id}: BAM does not exist: {path}")
        resolved.append(path)
    return root, tuple(resolved)


def load_project(config_path: str | Path) -> ResolvedProject:
    path = Path(config_path).expanduser().resolve()
    if not path.is_file():
        raise ProjectValidationError(f"Project configuration does not exist: {path}")
    raw_config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw_config, dict):
        raise ProjectValidationError("Project configuration must be a YAML mapping.")
    config = normalize_config(raw_config)

    schema_errors = sorted(
        Draft202012Validator(_schema()).iter_errors(config),
        key=lambda error: tuple(error.absolute_path),
    )
    errors = [
        f"schema {'.'.join(map(str, error.absolute_path)) or '<root>'}: {error.message}"
        for error in schema_errors
    ]
    if errors:
        raise ProjectValidationError("\n".join(errors))

    base = path.parent
    inputs = config["inputs"]
    kind = str(inputs["kind"])
    samples_path = _resolve(base, inputs["samples"])
    contrasts_path = _resolve(base, config["analysis"]["contrasts"])
    gmt_path = _resolve(base, config["resources"]["gene_sets"]["gmt"])
    for name, resolved in (
        ("samples", samples_path),
        ("contrasts", contrasts_path),
        ("gmt", gmt_path),
    ):
        if not resolved.is_file():
            errors.append(f"{name} file does not exist: {resolved}")
    if errors:
        raise ProjectValidationError("\n".join(errors))

    sample_header, all_samples = _read_tsv(samples_path)
    if "sample_id" not in sample_header:
        errors.append("samples.tsv requires a sample_id column")
    if not all_samples:
        errors.append("samples.tsv contains no samples")
    if errors:
        raise ProjectValidationError("\n".join(errors))
    try:
        samples = _selected_samples(all_samples, sample_header, inputs.get("analysis_set"))
    except ProjectValidationError as error:
        errors.append(str(error))
        samples = []

    sample_ids = [row.get("sample_id", "").strip() for row in samples]
    if any(not sample_id for sample_id in sample_ids):
        errors.append("sample_id values must be non-empty")
    if len(sample_ids) != len(set(sample_ids)):
        errors.append("selected sample_id values must be unique")

    counts_path: Path | None = None
    annotation_path: Path | None = None
    gtf_path: Path | None = None
    source_root: Path | None = None
    bam_paths: tuple[Path, ...] = ()
    archive_path: Path | None = None
    gene_ids: list[str] = []

    if kind == "counts":
        counts_path = _resolve(base, inputs["counts"])
        annotation_path = _resolve(base, inputs["annotation"])
        for name, resolved in (("counts", counts_path), ("annotation", annotation_path)):
            if not resolved.is_file():
                errors.append(f"{name} file does not exist: {resolved}")
        if not errors:
            count_header, count_rows = _read_tsv(counts_path)
            if not count_header or count_header[0] != "gene_id":
                errors.append("counts.tsv must begin with a gene_id column")
            count_samples = count_header[1:]
            missing_selected = sorted(set(sample_ids) - set(count_samples))
            unknown_counts = sorted(set(count_samples) - {row.get("sample_id", "") for row in all_samples})
            if missing_selected or unknown_counts:
                errors.append(
                    "count-matrix columns must contain every selected sample and no undeclared samples; "
                    f"missing_selected={missing_selected}, unknown_counts={unknown_counts}"
                )
            if inputs.get("analysis_set") is None and set(count_samples) != set(sample_ids):
                errors.append("count-matrix sample columns must exactly match sample_id values")
            if not count_rows:
                errors.append("counts.tsv contains no genes")
            for row_number, row in enumerate(count_rows, start=2):
                gene_id = row.get("gene_id", "").strip()
                gene_ids.append(gene_id)
                for sample_id in count_samples:
                    raw = row.get(sample_id, "")
                    try:
                        value = int(raw)
                    except (TypeError, ValueError):
                        errors.append(f"counts.tsv row {row_number}, {sample_id}: expected an integer")
                        continue
                    if value < 0:
                        errors.append(f"counts.tsv row {row_number}, {sample_id}: counts cannot be negative")
            if any(not gene_id for gene_id in gene_ids):
                errors.append("gene_id values must be non-empty")
            if len(gene_ids) != len(set(gene_ids)):
                errors.append("gene_id values must be unique")

            annotation_header, annotations = _read_tsv(annotation_path)
            if not {"gene_id", "gene_symbol"}.issubset(annotation_header):
                errors.append("annotation.tsv requires gene_id and gene_symbol columns")
            annotation_ids = {row.get("gene_id", "") for row in annotations}
            missing_annotation = sorted(set(gene_ids) - annotation_ids)
            if missing_annotation:
                errors.append(
                    f"annotation.tsv is missing {len(missing_annotation)} count gene_id values "
                    f"(first: {missing_annotation[:5]})"
                )
    elif kind in {"bam", "nfcore_rnaseq"}:
        gtf_path = _resolve(base, inputs["gtf"])
        if not gtf_path.is_file():
            errors.append(f"gtf file does not exist: {gtf_path}")
        counting = config["counting"]
        if counting["strandedness"] == "infer" and not {1, 2}.issubset(
            counting["strand_test_modes"]
        ):
            errors.append(
                "counting.strand_test_modes must include 1 and 2 when strandedness is infer"
            )
        if not counting["paired_end"] and any(
            counting[key]
            for key in (
                "count_read_pairs",
                "require_both_ends_aligned",
                "exclude_chimeric_fragments",
            )
        ):
            errors.append(
                "count_read_pairs, require_both_ends_aligned, and "
                "exclude_chimeric_fragments require counting.paired_end: true"
            )
        source_root, bam_paths = _resolve_bams(
            kind=kind,
            inputs=inputs,
            base=base,
            sample_header=sample_header,
            samples=samples,
            errors=errors,
        )
    else:
        archive_path = _resolve(base, inputs["archive"])
        gtf_path = _resolve(base, inputs["gtf"])
        if not archive_path.is_file():
            errors.append(f"archive file does not exist: {archive_path}")
        if not gtf_path.is_file():
            errors.append(f"gtf file does not exist: {gtf_path}")
        if "bam" not in sample_header:
            errors.append("archive input requires a bam column containing archive member paths")
        else:
            member_root = inputs.get("member_root", "").strip("/\\")
            for row in samples:
                relative = row.get("bam", "").strip().replace("\\", "/")
                combined = f"{member_root}/{relative}" if member_root else relative
                member = PurePosixPath(combined)
                if (
                    not relative
                    or member.is_absolute()
                    or any(part in {"", ".", ".."} for part in member.parts)
                ):
                    errors.append(
                        f"sample {row.get('sample_id', '<missing>')}: unsafe archive BAM member {combined!r}"
                    )

    contrast_header, contrasts = _read_tsv(contrasts_path)
    required_contrast = {"contrast_id", "factor", "numerator", "denominator"}
    if not required_contrast.issubset(contrast_header):
        errors.append("contrasts.tsv requires contrast_id, factor, numerator, and denominator columns")
    if not contrasts:
        errors.append("contrasts.tsv contains no contrasts")
    contrast_ids = [row.get("contrast_id", "").strip() for row in contrasts]
    if len(contrast_ids) != len(set(contrast_ids)):
        errors.append("contrast_id values must be unique")
    for row in contrasts:
        contrast_id = row.get("contrast_id", "<missing>")
        if not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", contrast_id):
            errors.append(f"invalid contrast_id {contrast_id!r}")
        factor = row.get("factor", "")
        if factor not in sample_header:
            errors.append(f"contrast {contrast_id}: factor {factor!r} is absent from samples.tsv")
            continue
        levels = {sample.get(factor, "") for sample in samples}
        numerator = row.get("numerator", "")
        denominator = row.get("denominator", "")
        if numerator == denominator:
            errors.append(f"contrast {contrast_id}: numerator and denominator must differ")
        for level in (numerator, denominator):
            if level not in levels:
                errors.append(f"contrast {contrast_id}: level {level!r} is absent from {factor}")
        if not re.search(rf"\b{re.escape(factor)}\b", config["analysis"]["design"]):
            errors.append(f"contrast {contrast_id}: factor {factor!r} is absent from design formula")

    group = config["figures"]["group"]
    if group not in sample_header:
        errors.append(f"figures.group {group!r} is absent from samples.tsv")
    else:
        groups = {row[group] for row in samples}
        palette_groups = set(config["figures"]["palette"])
        missing_colors = sorted(groups - palette_groups)
        if missing_colors:
            errors.append(f"figures.palette has no colors for: {missing_colors}")

    gene_sets = config["resources"]["gene_sets"]
    if gene_sets["min_size"] > gene_sets["max_size"]:
        errors.append("gene_sets.min_size cannot exceed gene_sets.max_size")
    modules = resolve_modules(config)
    if not any(modules.values()):
        errors.append("at least one analysis module must be enabled")
    if modules["pathways"] and not modules["de"]:
        errors.append("modules.pathways requires modules.de")
    if modules["pathways"] and not modules["qc"]:
        errors.append("modules.pathways requires modules.qc because it consumes VST expression")
    if modules["ontology"] and not modules["pathways"]:
        errors.append("modules.ontology requires modules.pathways")
    if any(modules[name] for name in ("composition", "regulators", "hypotheses")) and not modules["qc"]:
        errors.append("composition, regulators, and hypotheses require modules.qc")
    if modules["hypotheses"] and (not modules["de"] or not modules["pathways"]):
        errors.append("modules.hypotheses requires modules.de and modules.pathways")
    if modules["publication"] and (not modules["qc"] or not modules["de"]):
        errors.append("modules.publication requires modules.qc and modules.de")
    if any(modules[name] for name in ("regulators", "networks")) and not modules["de"]:
        errors.append("regulators and networks require modules.de")
    if modules["networks"] and not modules["ontology"]:
        errors.append("modules.networks requires modules.ontology")
    if modules["networks"]:
        if not config["resources"].get("providers", {}).get("string", False):
            errors.append("modules.networks requires resources.providers.string: true")
        if config["species"].get("taxonomy_id") is None:
            errors.append("modules.networks requires species.taxonomy_id")

    signature_path: Path | None = None
    regulon_path: Path | None = None
    if modules["composition"]:
        signature_value = config["resources"].get("cell_state_signatures")
        if not signature_value:
            errors.append("modules.composition requires resources.cell_state_signatures")
        else:
            signature_path = _resolve(base, signature_value)
            if not signature_path.is_file():
                errors.append(f"cell-state signature file does not exist: {signature_path}")
            else:
                signature_config, document_errors = _load_schema_document(signature_path, "signatures")
                errors.extend(document_errors)
                if signature_config is not None:
                    ids = [item["id"] for item in signature_config["signatures"]]
                    if len(ids) != len(set(ids)):
                        errors.append("cell-state signature ids must be unique")
    if modules["regulators"] and config["resources"].get("regulon_edges"):
        regulon_path = _resolve(base, config["resources"]["regulon_edges"])
        if not regulon_path.is_file():
            errors.append(f"regulon edge file does not exist: {regulon_path}")
        else:
            regulon_header, _ = _read_tsv(regulon_path)
            if not {"source", "target"}.issubset(regulon_header):
                errors.append("resources.regulon_edges requires source and target columns")
    if modules["regulators"] and regulon_path is None:
        if config["resources"].get("providers", {}).get("gtrd", False):
            errors.append(
                "resources.providers.gtrd requires resources.regulon_edges pointing to an "
                "exported GTRD-derived source/target snapshot; BulkRNAFrame does not silently "
                "redistribute GTRD data"
            )
        if not config["resources"].get("providers", {}).get("dorothea", False):
            errors.append("modules.regulators requires resources.regulon_edges or resources.providers.dorothea: true")

    hypothesis_path: Path | None = None
    panel_path: Path | None = None
    recipe_path: Path | None = None
    hypothesis_config: dict[str, Any] | None = None
    panel_config: dict[str, Any] | None = None
    recipe_config: dict[str, Any] | None = None
    if modules["hypotheses"] or modules["publication"]:
        if "hypotheses" not in config:
            errors.append(
                "hypothesis/publication modules require hypotheses.panels; the hypothesis "
                "module also requires hypotheses.claims"
            )
        else:
            panel_path = _resolve(base, config["hypotheses"]["panels"])
            if not panel_path.is_file():
                errors.append(f"hypothesis panels file does not exist: {panel_path}")
            else:
                panel_config, document_errors = _load_schema_document(
                    panel_path, "hypothesis_panels"
                )
                errors.extend(document_errors)
            if modules["hypotheses"]:
                hypothesis_path = _resolve(base, config["hypotheses"]["claims"])
                if not hypothesis_path.is_file():
                    errors.append(f"hypothesis claims file does not exist: {hypothesis_path}")
                else:
                    hypothesis_config, document_errors = _load_schema_document(
                        hypothesis_path, "hypotheses"
                    )
                    errors.extend(document_errors)
    if modules["publication"]:
        if "publication" not in config:
            errors.append("modules.publication requires publication.recipe")
        else:
            recipe_path = _resolve(base, config["publication"]["recipe"])
            if not recipe_path.is_file():
                errors.append(f"figure recipe file does not exist: {recipe_path}")
            else:
                recipe_config, document_errors = _load_schema_document(
                    recipe_path, "figure_recipe"
                )
                errors.extend(document_errors)
    if hypothesis_config is not None:
        known_contrasts = set(contrast_ids)
        known_gene_panels = set((panel_config or {}).get("gene_panels", {})) | set(
            (panel_config or {}).get("programs", {})
        )
        hypothesis_ids = [item["id"] for item in hypothesis_config["hypotheses"]]
        if len(hypothesis_ids) != len(set(hypothesis_ids)):
            errors.append("hypothesis ids must be unique")
        for item in hypothesis_config["hypotheses"]:
            if item["contrast"] not in known_contrasts:
                errors.append(
                    f"hypothesis {item['id']}: unknown contrast {item['contrast']!r}"
                )
            if panel_config is not None:
                for panel in item.get("gene_panels", []):
                    if panel not in known_gene_panels:
                        errors.append(f"hypothesis {item['id']}: unknown gene panel {panel!r}")
                for panel in item.get("pathway_panels", []):
                    if panel not in panel_config.get("pathway_panels", {}):
                        errors.append(f"hypothesis {item['id']}: unknown pathway panel {panel!r}")
    if panel_config is not None:
        known_gene_panels = set(panel_config.get("gene_panels", {})) | set(
            panel_config.get("programs", {})
        )
        known_pathway_panels = set(panel_config.get("pathway_panels", {}))
        order = panel_config.get("figure_order", {})
        for panel in order.get("gene_panels", []):
            if panel not in known_gene_panels:
                errors.append(f"figure_order references unknown gene panel {panel!r}")
        for panel in order.get("pathway_panels", []):
            if panel not in known_pathway_panels:
                errors.append(f"figure_order references unknown pathway panel {panel!r}")
        for panel in panel_config.get("gsea_programs", []):
            if panel not in known_gene_panels:
                errors.append(f"gsea_programs references unknown gene panel {panel!r}")
        for effect_id, effect in panel_config.get("expected_effects", {}).items():
            if effect["contrast"] not in set(contrast_ids):
                errors.append(f"expected_effects {effect_id!r}: unknown contrast {effect['contrast']!r}")
            if effect["target_type"] == "program" and effect["target"] not in known_gene_panels:
                errors.append(f"expected_effects {effect_id!r}: unknown program {effect['target']!r}")
            if effect["target_type"] == "pathway" and effect["target"] not in known_pathway_panels:
                errors.append(f"expected_effects {effect_id!r}: unknown pathway panel {effect['target']!r}")
    if recipe_config is not None:
        for figure_set, recipe in recipe_config["figure_sets"].items():
            panel_ids = [panel["id"] for panel in recipe["panels"]]
            if len(panel_ids) != len(set(panel_ids)):
                errors.append(f"figure set {figure_set!r} contains duplicate panel ids")
        from .figures import validate_recipe_contract

        errors.extend(validate_recipe_contract(recipe_config, (name for name, enabled in modules.items() if enabled), contrast_ids))
    valid_gmt_lines = 0
    with gmt_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                errors.append(
                    f"GMT line {line_number} must contain a name, description, and at least two genes"
                )
            else:
                valid_gmt_lines += 1
    if valid_gmt_lines == 0:
        errors.append("GMT file contains no valid gene sets")

    if errors:
        raise ProjectValidationError("\n".join(errors))

    return ResolvedProject(
        config_path=path,
        config=config,
        source_kind=kind,
        counts=counts_path,
        samples=samples_path,
        annotation=annotation_path,
        contrasts=contrasts_path,
        gmt=gmt_path,
        gtf=gtf_path,
        source_root=source_root,
        bam_paths=bam_paths,
        archive=archive_path,
        analysis_set=str(inputs.get("analysis_set", "all")),
        modules=tuple(sorted(name for name, enabled in modules.items() if enabled)),
        hypotheses=hypothesis_path,
        hypothesis_panels=panel_path,
        figure_recipe=recipe_path,
        hypothesis_config=hypothesis_config,
        panel_config=panel_config,
        recipe_config=recipe_config,
        cell_state_signatures=signature_path,
        regulon_edges=regulon_path,
        output_root=_resolve(base, config["output"]["root"]),
        sample_rows=tuple(samples),
        contrast_rows=tuple(contrasts),
    )


def validation_report(project: ResolvedProject) -> dict[str, Any]:
    report: dict[str, Any] = {
        "status": "ok",
        "project_id": project.project_id,
        "config": str(project.config_path),
        "samples": len(project.sample_rows),
        "contrasts": [row["contrast_id"] for row in project.contrast_rows],
        "input_kind": project.source_kind,
        "analysis_set": project.analysis_set,
        "profile": project.config["analysis"]["profile"],
        "modules": list(project.modules),
        "species": project.config["species"],
        "output": str(project.result_root),
    }
    if project.bam_paths:
        report["bams"] = len(project.bam_paths)
        report["source_root"] = str(project.source_root)
    return report


def report_json(project: ResolvedProject) -> str:
    return json.dumps(validation_report(project), indent=2) + "\n"


@dataclass(frozen=True)
class ResolvedCollection:
    config_path: Path
    config: dict[str, Any]
    projects: tuple[ResolvedProject, ...]
    output_root: Path

    @property
    def collection_id(self) -> str:
        return str(self.config["collection"]["id"])

    @property
    def result_root(self) -> Path:
        return self.output_root / self.collection_id


def load_collection(config_path: str | Path) -> ResolvedCollection:
    path = Path(config_path).expanduser().resolve()
    if not path.is_file():
        raise ProjectValidationError(f"Collection configuration does not exist: {path}")
    document, errors = _load_schema_document(path, "collection")
    if document is None:
        raise ProjectValidationError("\n".join(errors))
    studies = document["studies"]
    study_ids = [study["id"] for study in studies]
    if len(study_ids) != len(set(study_ids)):
        errors.append("collection study ids must be unique")
    projects: list[ResolvedProject] = []
    for study in studies:
        project_path = _resolve(path.parent, study["project"])
        try:
            project = load_project(project_path)
        except ProjectValidationError as error:
            errors.append(f"study {study['id']}: {error}")
            continue
        if study["contrast"] not in {row["contrast_id"] for row in project.contrast_rows}:
            errors.append(
                f"study {study['id']}: contrast {study['contrast']!r} is absent from project"
            )
        projects.append(project)
    if errors:
        raise ProjectValidationError("\n".join(errors))
    return ResolvedCollection(
        config_path=path,
        config=document,
        projects=tuple(projects),
        output_root=_resolve(path.parent, document["output"]["root"]),
    )


def collection_report(collection: ResolvedCollection) -> str:
    data = {
        "status": "ok",
        "collection_id": collection.collection_id,
        "studies": [
            {
                "id": study["id"],
                "project": str(project.config_path),
                "contrast": study["contrast"],
            }
            for study, project in zip(
                collection.config["studies"], collection.projects, strict=True
            )
        ],
        "output": str(collection.result_root),
    }
    return json.dumps(data, indent=2) + "\n"
