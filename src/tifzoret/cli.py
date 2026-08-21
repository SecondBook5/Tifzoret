"""Command-line interface for Tifzoret."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from importlib import resources
from pathlib import Path

import yaml

from .config import (
    ProjectValidationError,
    collection_report,
    load_collection,
    load_project,
    migrate_v1_mapping,
    report_json,
)
from .collection import run_collection
from .figures import build_gallery, constructor_catalog, initialize_figure_workflow
from .verification import verify_project, write_verification


def _copy_tree(source, destination: Path) -> None:
    """Recursively copy a packaged template tree into ``destination``."""
    destination.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        target = destination / item.name
        if item.is_dir():
            _copy_tree(item, target)
        else:
            target.write_bytes(item.read_bytes())


def command_init(args: argparse.Namespace) -> int:
    """Scaffold a new project directory from an input-boundary template.

    ``--input`` chooses the template (``counts``, ``bam``, ``nfcore-rnaseq``, or
    ``archive``). The destination must not already contain files. Returns 0.
    """
    destination = Path(args.directory).expanduser().resolve()
    if destination.exists() and any(destination.iterdir()):
        raise ProjectValidationError(f"Destination is not empty: {destination}")
    template_name = {
        "counts": "minimal",
        "bam": "bam",
        "nfcore-rnaseq": "nfcore_rnaseq",
        "archive": "archive",
    }[args.input]
    template = resources.files("tifzoret").joinpath(f"templates/{template_name}")
    _copy_tree(template, destination)
    print(f"Created {args.input} project template: {destination}")
    return 0


def command_validate(args: argparse.Namespace) -> int:
    """Validate the project and every tabular input, printing a JSON summary.

    Loading raises :class:`ProjectValidationError` (exit 2 via :func:`main`) on
    any contract violation; on success prints the resolved summary and returns 0.
    """
    project = load_project(args.project)
    print(report_json(project), end="")
    return 0


def _snakemake(
    args: argparse.Namespace, dry_run: bool, targets: tuple[str, ...] = ()
) -> int:
    """Invoke Snakemake against the packaged workflow for a validated project.

    The project is validated first, then Snakemake runs with the study's config
    file, from the config directory (so relative paths resolve identically to
    ``validate``). ``targets`` defaults to the release ``manifest.json``; conda
    provisioning is on unless ``--no-conda`` is given, and ``dry_run`` adds
    ``--dry-run``. Returns Snakemake's exit code unchanged.
    """
    project = load_project(args.project)
    if not targets:
        targets = (str(project.result_root / "manifest.json"),)
    snakefile = resources.files("tifzoret").joinpath("workflow/Snakefile")
    command = [
        args.snakemake,
        "--snakefile",
        str(snakefile),
        "--configfile",
        str(project.config_path),
        "--cores",
        str(args.cores),
        "--rerun-incomplete",
        "--printshellcmds",
    ]
    if not args.no_conda:
        command.append("--use-conda")
    if dry_run:
        command.append("--dry-run")
    command.extend(targets)
    return subprocess.run(command, cwd=project.config_path.parent, check=False).returncode


def command_prepare(args: argparse.Namespace) -> int:
    """Run only up to the canonical inputs (counting/annotation), then stop.

    Targets the input manifest so a study's data boundary can be materialized
    and inspected before committing to the full analysis. Returns Snakemake's code.
    """
    project = load_project(args.project)
    target = project.result_root / "inputs" / "input_manifest.json"
    return _snakemake(args, dry_run=False, targets=(str(target),))


def command_report(args: argparse.Namespace) -> int:
    """Build the navigable ``REPORT.html`` (running any missing upstream steps)."""
    project = load_project(args.project)
    target = project.result_root / "REPORT.html"
    return _snakemake(args, dry_run=False, targets=(str(target),))


def command_assemble(args: argparse.Namespace) -> int:
    """Assemble the configured multi-panel publication figure(s).

    Requires ``publication.recipe`` (a file or inline mapping); raises
    :class:`ProjectValidationError` otherwise. Returns Snakemake's exit code.
    """
    project = load_project(args.project)
    if project.figure_recipe is None:
        raise ProjectValidationError("assemble requires publication.recipe in project.yaml")
    return _snakemake(args, dry_run=False, targets=("assemble_publication",))


def command_figures_init(args: argparse.Namespace) -> int:
    """Scaffold the three publication companion files and validate the recipe.

    Writes ``hypotheses.yaml``, ``hypothesis_panels.yaml``, and
    ``figure_recipe.yaml`` beside the project (skipped unless ``--force`` when
    they already exist, which surfaces as exit 2). Returns 0 on success.
    """
    try:
        outputs = initialize_figure_workflow(args.project, force=args.force)
        project = load_project(args.project)
    except (OSError, ValueError, yaml.YAMLError) as error:
        raise ProjectValidationError(str(error)) from error
    print("Initialized hypothesis-driven publication files:")
    for output in outputs:
        print(f"  {output}")
    print(f"Validated constructor recipe for: {project.project_id}")
    return 0


def command_figures_catalog(args: argparse.Namespace) -> int:
    """List the registered figure constructors, their variants, and scope.

    ``--json`` emits the machine-readable catalog; otherwise prints an aligned
    table of ``id``, contrast-vs-study scope, and available variants. Returns 0.
    """
    catalog = constructor_catalog()
    if args.json:
        print(json.dumps({"constructors": catalog}, indent=2))
    else:
        for constructor in catalog:
            variants = ", ".join(constructor["variants"])
            scope = "contrast" if constructor["contrast_specific"] else "study"
            print(f"{constructor['id']:<28} {scope:<8} {variants}")
    return 0


def command_figures_build(args: argparse.Namespace) -> int:
    """Render the recipe's panels and assemble its figure set(s).

    With ``--figure-set`` only that set's panels are built (targeting its panel
    index); without it, every configured set is assembled. Requires
    ``publication.recipe``. Returns Snakemake's exit code.
    """
    project = load_project(args.project)
    if project.figure_recipe is None:
        raise ProjectValidationError("figures build requires publication.recipe in project.yaml")
    if args.figure_set:
        if args.figure_set not in project.recipe_config["figure_sets"]:
            raise ProjectValidationError(f"unknown figure set: {args.figure_set}")
        target = project.result_root / "publication" / args.figure_set / "panels" / "index.json"
        return _snakemake(args, dry_run=False, targets=(str(target),))
    return _snakemake(args, dry_run=False, targets=("assemble_publication",))


def command_figures_gallery(args: argparse.Namespace) -> int:
    """Build an HTML + contact-sheet review of every built panel variant.

    Lets a reviewer compare variants side by side before choosing the final
    panel. Requires ``publication.recipe``; ``--output`` overrides the path.
    Returns 0.
    """
    project = load_project(args.project)
    if project.figure_recipe is None:
        raise ProjectValidationError("figures gallery requires publication.recipe in project.yaml")
    output = build_gallery(project, args.output)
    print(f"Figure review gallery written: {output}")
    return 0


def command_verify(args: argparse.Namespace) -> int:
    """Compare a project's result against a reference run within tolerances.

    ``--scope`` restricts the comparison (``counts``, ``core``, or ``all``);
    ``--candidate`` overrides the directory checked (defaults to the resolved
    result). Writes a verification report and returns 0 when it passes, 1 otherwise.
    """
    project = load_project(args.project)
    candidate = Path(args.candidate).expanduser().resolve() if args.candidate else project.result_root
    result = verify_project(
        project, args.reference, candidate, atol=args.atol, rtol=args.rtol, scope=args.scope
    )
    output = Path(args.output).expanduser().resolve()
    write_verification(result, output)
    print(f"Verification {'passed' if result.passed else 'failed'}: {output}")
    return 0 if result.passed else 1


def command_collection_validate(args: argparse.Namespace) -> int:
    """Validate a multi-study collection and print its resolved study summary."""
    print(collection_report(load_collection(args.collection)), end="")
    return 0


def command_collection_run(args: argparse.Namespace) -> int:
    """Run the cross-study (meta-analysis) collection and report the output path."""
    output = run_collection(load_collection(args.collection))
    print(f"Collection analysis written: {output}")
    return 0


def command_migrate_config(args: argparse.Namespace) -> int:
    """Convert a development v1 configuration into the public v2 contract.

    Applies the requested species/reference metadata (v1 predates these),
    refuses to overwrite an existing ``--output`` unless ``--force``, and
    re-validates the written file in place so relative paths resolve. Returns 0.
    """
    source = Path(args.project).expanduser().resolve()
    raw = yaml.safe_load(source.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ProjectValidationError("Project configuration must be a YAML mapping.")
    migrated = migrate_v1_mapping(raw)
    species_defaults = {
        "mouse": ("Mus musculus", 10090),
        "human": ("Homo sapiens", 9606),
        "custom": ("unspecified", None),
    }
    scientific_name, taxonomy_id = species_defaults[args.species]
    migrated["species"] = {
        "provider": args.species,
        "scientific_name": args.scientific_name or scientific_name,
        "taxonomy_id": args.taxonomy_id if args.taxonomy_id is not None else taxonomy_id,
    }
    migrated["reference"]["genome_build"] = args.genome_build
    migrated["reference"]["annotation_release"] = args.annotation_release
    output = Path(args.output).expanduser().resolve()
    if output.exists() and not args.force:
        raise ProjectValidationError(f"Refusing to overwrite existing output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(yaml.safe_dump(migrated, sort_keys=False), encoding="utf-8")
    # Validate in its final directory because all relative paths are resolved
    # relative to the configuration file.
    load_project(output)
    print(f"Migrated configuration v1 -> v2: {output}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Construct the ``tifzoret`` argument parser with every subcommand wired to a handler."""
    parser = argparse.ArgumentParser(prog="tifzoret")
    commands = parser.add_subparsers(dest="command", required=True)

    init_parser = commands.add_parser("init", help="create a project scaffold")
    init_parser.add_argument("directory")
    init_parser.add_argument(
        "--input",
        choices=("counts", "bam", "nfcore-rnaseq", "archive"),
        default="counts",
        help="source boundary for the scaffold (default: counts)",
    )
    init_parser.set_defaults(handler=command_init)

    validate_parser = commands.add_parser("validate", help="validate a project and all tabular inputs")
    validate_parser.add_argument("project")
    validate_parser.set_defaults(handler=command_validate)

    migrate_parser = commands.add_parser(
        "migrate-config", help="convert a version 1 project configuration to version 2"
    )
    migrate_parser.add_argument("project")
    migrate_parser.add_argument("--output", required=True)
    migrate_parser.add_argument("--species", choices=("mouse", "human", "custom"), default="custom")
    migrate_parser.add_argument("--scientific-name")
    migrate_parser.add_argument("--taxonomy-id", type=int)
    migrate_parser.add_argument("--genome-build", default="unspecified")
    migrate_parser.add_argument("--annotation-release")
    migrate_parser.add_argument("--force", action="store_true")
    migrate_parser.set_defaults(handler=command_migrate_config)

    prepare_parser = commands.add_parser(
        "prepare", help="validate and materialize canonical inputs only"
    )
    prepare_parser.add_argument("project")
    prepare_parser.add_argument("--cores", type=int, default=1)
    prepare_parser.add_argument("--snakemake", default="snakemake")
    prepare_parser.add_argument("--no-conda", action="store_true")
    prepare_parser.set_defaults(handler=command_prepare)

    report_parser = commands.add_parser("report", help="build the navigable HTML report")
    report_parser.add_argument("project")
    report_parser.add_argument("--cores", type=int, default=1)
    report_parser.add_argument("--snakemake", default="snakemake")
    report_parser.add_argument("--no-conda", action="store_true")
    report_parser.set_defaults(handler=command_report)

    assemble_parser = commands.add_parser("assemble", help="assemble configured publication figures")
    assemble_parser.add_argument("project")
    assemble_parser.add_argument("--cores", type=int, default=1)
    assemble_parser.add_argument("--snakemake", default="snakemake")
    assemble_parser.add_argument("--no-conda", action="store_true")
    assemble_parser.set_defaults(handler=command_assemble)

    figures_parser = commands.add_parser(
        "figures", help="initialize, build, and review hypothesis-driven publication figures"
    )
    figure_commands = figures_parser.add_subparsers(dest="figures_command", required=True)
    figures_init = figure_commands.add_parser(
        "init", help="scaffold hypothesis, program, and constructor recipe files"
    )
    figures_init.add_argument("project")
    figures_init.add_argument("--force", action="store_true")
    figures_init.set_defaults(handler=command_figures_init)
    figures_catalog = figure_commands.add_parser(
        "catalog", help="list registered publication constructors and variants"
    )
    figures_catalog.add_argument("--json", action="store_true")
    figures_catalog.set_defaults(handler=command_figures_catalog)
    figures_build = figure_commands.add_parser(
        "build", help="render configured panels and assemble publication figure sets"
    )
    figures_build.add_argument("project")
    figures_build.add_argument("--figure-set")
    figures_build.add_argument("--cores", type=int, default=1)
    figures_build.add_argument("--snakemake", default="snakemake")
    figures_build.add_argument("--no-conda", action="store_true")
    figures_build.set_defaults(handler=command_figures_build)
    figures_gallery = figure_commands.add_parser(
        "gallery", help="create an HTML/contact-sheet review of all built variants"
    )
    figures_gallery.add_argument("project")
    figures_gallery.add_argument("--output")
    figures_gallery.set_defaults(handler=command_figures_gallery)

    verify_parser = commands.add_parser("verify", help="compare a project result with a reference run")
    verify_parser.add_argument("project", help="project.yaml whose resolved result is the candidate")
    verify_parser.add_argument("--reference", required=True)
    verify_parser.add_argument("--candidate", help="override the candidate result directory")
    verify_parser.add_argument("--scope", choices=("counts", "core", "all"), default="all")
    verify_parser.add_argument("--output", default="verification.json")
    verify_parser.add_argument("--atol", type=float, default=1e-8)
    verify_parser.add_argument("--rtol", type=float, default=1e-6)
    verify_parser.set_defaults(handler=command_verify)

    collection_parser = commands.add_parser("collection", help="validate or run a multi-study collection")
    collection_commands = collection_parser.add_subparsers(dest="collection_command", required=True)
    collection_validate = collection_commands.add_parser("validate")
    collection_validate.add_argument("collection")
    collection_validate.set_defaults(handler=command_collection_validate)
    collection_run = collection_commands.add_parser("run")
    collection_run.add_argument("collection")
    collection_run.set_defaults(handler=command_collection_run)

    for name, dry_run in (("dry-run", True), ("run", False)):
        run_parser = commands.add_parser(name)
        run_parser.add_argument("project")
        run_parser.add_argument("--cores", type=int, default=1)
        run_parser.add_argument("--snakemake", default="snakemake")
        run_parser.add_argument("--no-conda", action="store_true")
        run_parser.set_defaults(handler=lambda args, dry_run=dry_run: _snakemake(args, dry_run))
    return parser


def main(argv: list[str] | None = None) -> int:
    """Parse arguments, dispatch to the selected handler, and map errors to exit codes.

    Returns the handler's exit code, 2 for a :class:`ProjectValidationError`
    (with the message on stderr), or 127 when a required executable is missing.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except ProjectValidationError as error:
        print(f"validation error:\n{error}", file=sys.stderr)
        return 2
    except FileNotFoundError as error:
        print(f"command not found: {error.filename}", file=sys.stderr)
        return 127


if __name__ == "__main__":
    raise SystemExit(main())
