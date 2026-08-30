"""Multi-study collection loading and validation."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .types import ProjectValidationError, ResolvedProject
from .paths import _resolve
from .validate import _load_schema_document
from .load import load_project


@dataclass(frozen=True)
class ResolvedCollection:
    """A validated multi-study collection for cross-study meta-analysis.

    Holds each member study as an already-validated :class:`ResolvedProject`
    plus the chosen contrast per study, so the meta-analysis stage combines only
    comparisons that actually exist and share consistent direction semantics.
    """

    config_path: Path
    config: dict[str, Any]
    projects: tuple[ResolvedProject, ...]
    output_root: Path

    @property
    def collection_id(self) -> str:
        """The filesystem-safe collection identifier (``collection.id``)."""
        return str(self.config["collection"]["id"])

    @property
    def result_root(self) -> Path:
        """Where the collection writes: ``<output.root>/<collection.id>/``."""
        return self.output_root / self.collection_id


def load_collection(config_path: str | Path) -> ResolvedCollection:
    """Load and validate a collection: its schema, each member project, and its contrast.

    Every referenced study is loaded through :func:`load_project`, and each
    declared contrast must exist in its project. Study ids must be unique. All
    problems are collected and raised together as one
    :class:`ProjectValidationError`.
    """
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
    """Render a validated collection (its studies, projects, and contrasts) as JSON."""
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
