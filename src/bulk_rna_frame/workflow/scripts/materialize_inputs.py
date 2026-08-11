#!/usr/bin/env python3
"""Materialize supported sources into BulkRNAFrame's canonical input contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from bulk_rna_frame.config import ResolvedProject, load_project  # noqa: E402


STRAND_MODES = {"unstranded": 0, "forward": 1, "reverse": 2}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path, *, comments: bool = False) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        lines = (line for line in handle if not comments or not line.startswith("#"))
        reader = csv.DictReader(lines, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def materialize_samples(project: ResolvedProject, output: Path) -> None:
    header, _ = read_tsv(project.samples)
    write_tsv(output, header, [dict(row) for row in project.sample_rows])


def materialize_count_matrix(
    project: ResolvedProject, counts_output: Path, annotation_output: Path
) -> dict[str, object]:
    assert project.counts is not None and project.annotation is not None
    header, rows = read_tsv(project.counts)
    selected = [row["sample_id"] for row in project.sample_rows]
    fields = ["gene_id", *selected]
    write_tsv(counts_output, fields, rows)

    annotation_header, annotations = read_tsv(project.annotation)
    write_tsv(annotation_output, annotation_header, annotations)
    return {
        "kind": "counts",
        "source_files": [
            {"path": str(project.counts), "bytes": project.counts.stat().st_size, "sha256": sha256(project.counts)},
            {
                "path": str(project.annotation),
                "bytes": project.annotation.stat().st_size,
                "sha256": sha256(project.annotation),
            },
        ],
        "source_samples": header[1:],
        "selected_samples": selected,
    }


def featurecounts_command(
    project: ResolvedProject, output: Path, mode: int
) -> list[str]:
    assert project.gtf is not None
    counting = project.config["counting"]
    command = [
        "featureCounts",
        "-T",
        str(counting["threads"]),
        "-F",
        "GTF",
        "-t",
        str(counting["feature_type"]),
        "-g",
        str(counting["attribute"]),
        "-a",
        str(project.gtf),
        "-o",
        str(output),
        "-s",
        str(mode),
    ]
    if counting["paired_end"]:
        command.append("-p")
    if counting["count_read_pairs"]:
        command.append("--countReadPairs")
    if counting["require_both_ends_aligned"]:
        command.append("-B")
    if counting["exclude_chimeric_fragments"]:
        command.append("-C")
    command.extend(str(path) for path in project.bam_paths)
    return command


def run_featurecounts(project: ResolvedProject, output: Path, mode: int) -> str:
    result = subprocess.run(
        featurecounts_command(project, output, mode),
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(
            f"featureCounts failed in strand mode {mode}:\n{result.stdout}\n{result.stderr}"
        )
    return result.stdout + result.stderr


def infer_strand_mode(
    project: ResolvedProject, temporary: Path
) -> tuple[int, dict[int, float], Path, str, float]:
    modes = project.config["counting"]["strand_test_modes"]
    means: dict[int, float] = {}
    outputs: dict[int, Path] = {}
    logs: dict[int, str] = {}
    for mode in modes:
        output = temporary / f"strand_{mode}.tsv"
        log = run_featurecounts(project, output, int(mode))
        outputs[int(mode)] = output
        logs[int(mode)] = log
        values = [
            float(value)
            for value in re.findall(
                r"Successfully assigned (?:alignments|fragments)\s*:\s*\d+\s*\(([\d.]+)%\)",
                log,
            )
        ]
        if values:
            means[int(mode)] = sum(values) / len(values)
    if not {1, 2}.issubset(means):
        raise RuntimeError(
            "strandedness inference requires featureCounts assignment rates for modes 1 and 2"
        )
    directional_total = means[1] + means[2]
    dominance = max(means[1], means[2]) / directional_total if directional_total else 0.0
    threshold = float(project.config["counting"]["strand_min_dominance"])
    mode = max((1, 2), key=means.get) if dominance >= threshold else 0
    if mode not in outputs:
        output = temporary / f"strand_{mode}.tsv"
        logs[mode] = run_featurecounts(project, output, mode)
        outputs[mode] = output
    return mode, means, outputs[mode], logs[mode], dominance


def parse_featurecounts(
    raw_path: Path, project: ResolvedProject, counts_output: Path
) -> list[str]:
    header, rows = read_tsv(raw_path, comments=True)
    annotation_fields = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
    raw_samples = [name for name in header if name not in annotation_fields]
    if len(raw_samples) != len(project.bam_paths):
        raise RuntimeError(
            f"featureCounts returned {len(raw_samples)} sample columns for {len(project.bam_paths)} BAMs"
        )

    sample_ids = [row["sample_id"] for row in project.sample_rows]
    basename_to_sample: dict[str, str] = {}
    for sample_id, bam in zip(sample_ids, project.bam_paths, strict=True):
        if bam.name in basename_to_sample:
            raise RuntimeError(f"BAM basenames are not unique: {bam.name}")
        basename_to_sample[bam.name] = sample_id
    raw_to_sample: dict[str, str] = {}
    for raw in raw_samples:
        basename = Path(raw).name
        if basename not in basename_to_sample:
            raise RuntimeError(f"Could not match featureCounts column to a declared BAM: {raw}")
        raw_to_sample[raw] = basename_to_sample[basename]

    canonical: list[dict[str, object]] = []
    for row in rows:
        output_row: dict[str, object] = {"gene_id": row["Geneid"]}
        for raw, sample_id in raw_to_sample.items():
            value = int(row[raw])
            if value < 0:
                raise RuntimeError(f"featureCounts emitted a negative count for {row['Geneid']}")
            output_row[sample_id] = value
        canonical.append(output_row)
    write_tsv(counts_output, ["gene_id", *sample_ids], canonical)
    return [str(row["gene_id"]) for row in canonical]


def gtf_symbols(gtf: Path) -> dict[str, str]:
    symbols: dict[str, str] = {}
    with gtf.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attributes = fields[8]
            gene_match = re.search(r'(?:^|;\s*)gene_id "([^"]+)"', attributes)
            if gene_match is None:
                continue
            gene_id = gene_match.group(1)
            name_match = re.search(r'(?:^|;\s*)gene_name "([^"]+)"', attributes)
            symbol = name_match.group(1) if name_match else re.sub(r"\.\d+$", "", gene_id)
            symbols.setdefault(gene_id, symbol)
            symbols.setdefault(re.sub(r"\.\d+$", "", gene_id), symbol)
    return symbols


def write_annotation(gtf: Path, gene_ids: list[str], output: Path) -> None:
    symbols = gtf_symbols(gtf)
    rows = [
        {
            "gene_id": gene_id,
            "gene_symbol": symbols.get(gene_id, symbols.get(re.sub(r"\.\d+$", "", gene_id), gene_id)),
        }
        for gene_id in gene_ids
    ]
    write_tsv(output, ["gene_id", "gene_symbol"], rows)


def bam_record(sample_id: str, path: Path) -> dict[str, object]:
    quickcheck = subprocess.run(
        ["samtools", "quickcheck", "-v", str(path)], text=True, capture_output=True, check=False
    )
    if quickcheck.returncode:
        raise RuntimeError(f"samtools quickcheck failed for {path}: {quickcheck.stderr}")
    header = subprocess.run(
        ["samtools", "view", "-H", str(path)], text=True, capture_output=True, check=True
    ).stdout
    sort_match = re.search(r"^@HD\t.*(?:^|\t)SO:([^\t\n]+)", header, flags=re.MULTILINE)
    sort_order = sort_match.group(1) if sort_match else None
    if sort_order not in {"coordinate", "unknown"}:
        raise RuntimeError(f"BAM is not coordinate sorted: {path} (SO={sort_order!r})")
    stat = path.stat()
    return {
        "sample_id": sample_id,
        "path": str(path),
        "bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sort_order": sort_order,
        "header_sha256": hashlib.sha256(header.encode()).hexdigest(),
        "quickcheck": True,
    }


def command_version(command: list[str]) -> str:
    result = subprocess.run(
        command,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    lines = [line for line in (result.stdout + result.stderr).splitlines() if line.strip()]
    return lines[0] if lines else "unavailable"


def materialize_bams(
    project: ResolvedProject, counts_output: Path, annotation_output: Path
) -> dict[str, object]:
    assert project.gtf is not None
    sample_ids = [row["sample_id"] for row in project.sample_rows]
    records = [
        bam_record(sample_id, bam)
        for sample_id, bam in zip(sample_ids, project.bam_paths, strict=True)
    ]
    with tempfile.TemporaryDirectory(prefix="bulk-rna-frame-") as directory:
        temporary = Path(directory)
        configured = project.config["counting"]["strandedness"]
        strand_rates: dict[int, float] = {}
        strand_dominance: float | None = None
        if configured == "infer":
            mode, strand_rates, raw_counts, final_log, strand_dominance = infer_strand_mode(
                project, temporary
            )
        else:
            mode = STRAND_MODES[configured]
            raw_counts = temporary / "featurecounts.tsv"
            final_log = run_featurecounts(project, raw_counts, mode)
        gene_ids = parse_featurecounts(raw_counts, project, counts_output)
    write_annotation(project.gtf, gene_ids, annotation_output)
    return {
        "kind": project.source_kind,
        "upstream": "nf-core/rnaseq" if project.source_kind == "nfcore_rnaseq" else "aligned BAM files",
        "selected_samples": sample_ids,
        "bams": records,
        "gtf": {"path": str(project.gtf), "bytes": project.gtf.stat().st_size, "sha256": sha256(project.gtf)},
        "counting": {
            **project.config["counting"],
            "resolved_strand_mode": mode,
            "strand_test_mean_assignment_percent": strand_rates,
            "strand_directional_dominance": strand_dominance,
            "final_assignment_log": final_log,
        },
        "tools": {
            "samtools": command_version(["samtools", "--version"]),
            "featureCounts": command_version(["featureCounts", "-v"]),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--threads", type=int)
    args = parser.parse_args()

    project = load_project(args.project_config)
    if args.threads is not None and project.source_kind != "counts":
        if args.threads < 1:
            parser.error("--threads must be at least 1")
        project.config["counting"]["threads"] = args.threads
    counts_output = Path(args.counts).resolve()
    samples_output = Path(args.samples).resolve()
    annotation_output = Path(args.annotation).resolve()
    manifest_output = Path(args.manifest).resolve()
    for output in (counts_output, samples_output, annotation_output, manifest_output):
        output.parent.mkdir(parents=True, exist_ok=True)

    materialize_samples(project, samples_output)
    if project.source_kind == "counts":
        provenance = materialize_count_matrix(project, counts_output, annotation_output)
    else:
        provenance = materialize_bams(project, counts_output, annotation_output)
    manifest = {
        "schema_version": 1,
        "project_id": project.project_id,
        "source": provenance,
        "canonical": {
            "counts": {"path": str(counts_output), "sha256": sha256(counts_output)},
            "samples": {"path": str(samples_output), "sha256": sha256(samples_output)},
            "annotation": {"path": str(annotation_output), "sha256": sha256(annotation_output)},
        },
    }
    manifest_output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Materialized {project.source_kind} input for {len(project.sample_rows)} samples")


if __name__ == "__main__":
    main()
