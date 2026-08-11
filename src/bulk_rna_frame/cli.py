"""Command-line interface for BulkRNAFrame."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from importlib import resources
from pathlib import Path

from .config import ProjectValidationError, load_project, report_json


def _copy_tree(source, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        target = destination / item.name
        if item.is_dir():
            _copy_tree(item, target)
        else:
            target.write_bytes(item.read_bytes())


def command_init(args: argparse.Namespace) -> int:
    destination = Path(args.directory).expanduser().resolve()
    if destination.exists() and any(destination.iterdir()):
        raise ProjectValidationError(f"Destination is not empty: {destination}")
    template_name = {
        "counts": "minimal",
        "bam": "bam",
        "nfcore-rnaseq": "nfcore_rnaseq",
    }[args.input]
    template = resources.files("bulk_rna_frame").joinpath(f"templates/{template_name}")
    _copy_tree(template, destination)
    print(f"Created {args.input} project template: {destination}")
    return 0


def command_validate(args: argparse.Namespace) -> int:
    project = load_project(args.project)
    print(report_json(project), end="")
    return 0


def _snakemake(
    args: argparse.Namespace, dry_run: bool, targets: tuple[str, ...] = ()
) -> int:
    project = load_project(args.project)
    snakefile = resources.files("bulk_rna_frame").joinpath("workflow/Snakefile")
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
    project = load_project(args.project)
    target = project.output_root / project.project_id / "inputs" / "input_manifest.json"
    return _snakemake(args, dry_run=False, targets=(str(target),))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bulk-rna")
    commands = parser.add_subparsers(dest="command", required=True)

    init_parser = commands.add_parser("init", help="create a project scaffold")
    init_parser.add_argument("directory")
    init_parser.add_argument(
        "--input",
        choices=("counts", "bam", "nfcore-rnaseq"),
        default="counts",
        help="source boundary for the scaffold (default: counts)",
    )
    init_parser.set_defaults(handler=command_init)

    validate_parser = commands.add_parser("validate", help="validate a project and all tabular inputs")
    validate_parser.add_argument("project")
    validate_parser.set_defaults(handler=command_validate)

    prepare_parser = commands.add_parser(
        "prepare", help="validate and materialize canonical inputs only"
    )
    prepare_parser.add_argument("project")
    prepare_parser.add_argument("--cores", type=int, default=1)
    prepare_parser.add_argument("--snakemake", default="snakemake")
    prepare_parser.add_argument("--no-conda", action="store_true")
    prepare_parser.set_defaults(handler=command_prepare)

    for name, dry_run in (("dry-run", True), ("run", False)):
        run_parser = commands.add_parser(name)
        run_parser.add_argument("project")
        run_parser.add_argument("--cores", type=int, default=1)
        run_parser.add_argument("--snakemake", default="snakemake")
        run_parser.add_argument("--no-conda", action="store_true")
        run_parser.set_defaults(handler=lambda args, dry_run=dry_run: _snakemake(args, dry_run))
    return parser


def main(argv: list[str] | None = None) -> int:
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
