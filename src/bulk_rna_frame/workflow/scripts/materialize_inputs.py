#!/usr/bin/env python3
"""Materialize supported sources into BulkRNAFrame's canonical input contract."""

from __future__ import annotations

import argparse
import csv
from contextlib import contextmanager
from dataclasses import replace
import hashlib
import json
import re
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path
from pathlib import PurePosixPath
from typing import Iterator

SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from bulk_rna_frame.config import ResolvedProject, load_project  # noqa: E402


STRAND_MODES = {"unstranded": 0, "forward": 1, "reverse": 2}


def safe_archive_member(value: str) -> PurePosixPath:
    """Normalize one archive member and reject traversal or ambiguous paths."""
    normalized = value.replace("\\", "/")
    member = PurePosixPath(normalized)
    if (
        not normalized
        or member.is_absolute()
        or any(part in {"", ".", ".."} for part in member.parts)
        or normalized.endswith("/")
    ):
        raise RuntimeError(f"unsafe BAM archive member path: {value!r}")
    return member


@contextmanager
def extracted_archive_bams(project: ResolvedProject) -> Iterator[tuple[Path, ...]]:
    """Extract only declared BAM members into an isolated temporary directory."""
    assert project.archive is not None
    member_root = project.config["inputs"].get("member_root", "").strip("/\\")
    declared = [
        safe_archive_member(
            f"{member_root}/{row['bam'].strip()}" if member_root else row["bam"].strip()
        )
        for row in project.sample_rows
    ]
    if len(declared) != len(set(declared)):
        raise RuntimeError("archive BAM member paths must be unique")

    with tempfile.TemporaryDirectory(prefix="bulk-rna-frame-archive-") as directory:
        root = Path(directory)
        extracted: list[Path] = []
        if zipfile.is_zipfile(project.archive):
            with zipfile.ZipFile(project.archive) as archive:
                available = set(archive.namelist())
                for member in declared:
                    name = member.as_posix()
                    if name not in available:
                        raise RuntimeError(f"declared BAM is absent from archive: {name}")
                    info = archive.getinfo(name)
                    if info.is_dir():
                        raise RuntimeError(f"declared BAM archive member is a directory: {name}")
                    destination = root.joinpath(*member.parts)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    with archive.open(info) as source, destination.open("wb") as target:
                        for block in iter(lambda: source.read(1024 * 1024), b""):
                            target.write(block)
                    extracted.append(destination)
        elif tarfile.is_tarfile(project.archive):
            with tarfile.open(project.archive, mode="r:*") as archive:
                available = {item.name: item for item in archive.getmembers()}
                for member in declared:
                    name = member.as_posix()
                    info = available.get(name)
                    if info is None:
                        raise RuntimeError(f"declared BAM is absent from archive: {name}")
                    if not info.isfile() or info.issym() or info.islnk():
                        raise RuntimeError(f"declared BAM archive member is not a regular file: {name}")
                    source = archive.extractfile(info)
                    if source is None:
                        raise RuntimeError(f"could not read BAM archive member: {name}")
                    destination = root.joinpath(*member.parts)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    with source, destination.open("wb") as target:
                        for block in iter(lambda: source.read(1024 * 1024), b""):
                            target.write(block)
                    extracted.append(destination)
        else:
            raise RuntimeError(
                f"unsupported archive format for {project.archive}; expected ZIP or TAR"
            )
        yield tuple(extracted)


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


def materialize_contrasts(project: ResolvedProject, output: Path) -> None:
    header, _ = read_tsv(project.contrasts)
    write_tsv(output, header, [dict(row) for row in project.contrast_rows])


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


def gtf_annotations(gtf: Path) -> dict[str, dict[str, object]]:
    annotations: dict[str, dict[str, object]] = {}
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
            biotype_match = re.search(r'(?:^|;\s*)gene_(?:bio)?type "([^"]+)"', attributes)
            record = annotations.setdefault(
                gene_id,
                {
                    "gene_symbol": symbol,
                    "seqname": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "gene_biotype": biotype_match.group(1) if biotype_match else "",
                },
            )
            record["start"] = min(int(record["start"]), int(fields[3]))
            record["end"] = max(int(record["end"]), int(fields[4]))
    return annotations


def write_annotation(gtf: Path, gene_ids: list[str], output: Path) -> None:
    annotations = gtf_annotations(gtf)
    rows = []
    for gene_id in gene_ids:
        record = annotations.get(gene_id, annotations.get(re.sub(r"\.\d+$", "", gene_id), {}))
        rows.append({
            "gene_id": gene_id,
            "gene_symbol": record.get("gene_symbol", gene_id),
            "seqname": record.get("seqname", ""),
            "start": record.get("start", ""),
            "end": record.get("end", ""),
            "strand": record.get("strand", ""),
            "gene_biotype": record.get("gene_biotype", ""),
        })
    write_tsv(output, ["gene_id", "gene_symbol", "seqname", "start", "end", "strand", "gene_biotype"], rows)


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
        "sha256": sha256(path),
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


def materialize_archive(
    project: ResolvedProject, counts_output: Path, annotation_output: Path
) -> dict[str, object]:
    assert project.archive is not None
    member_root = project.config["inputs"].get("member_root", "").strip("/\\")
    member_names = [
        f"{member_root}/{row['bam'].strip()}" if member_root else row["bam"].strip()
        for row in project.sample_rows
    ]
    with extracted_archive_bams(project) as bams:
        extracted_project = replace(project, source_root=bams[0].parent, bam_paths=bams)
        provenance = materialize_bams(extracted_project, counts_output, annotation_output)
    provenance["kind"] = "archive"
    provenance["upstream"] = "archive containing aligned BAM files"
    provenance["archive"] = {
        "path": str(project.archive),
        "bytes": project.archive.stat().st_size,
        "sha256": sha256(project.archive),
        "members": member_names,
    }
    for record, member in zip(provenance["bams"], member_names, strict=True):
        record["path"] = f"archive://{project.archive.name}/{member}"
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--contrasts", required=True)
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
    contrasts_output = Path(args.contrasts).resolve()
    manifest_output = Path(args.manifest).resolve()
    for output in (counts_output, samples_output, annotation_output, contrasts_output, manifest_output):
        output.parent.mkdir(parents=True, exist_ok=True)

    materialize_samples(project, samples_output)
    materialize_contrasts(project, contrasts_output)
    if project.source_kind == "counts":
        provenance = materialize_count_matrix(project, counts_output, annotation_output)
    elif project.source_kind == "archive":
        provenance = materialize_archive(project, counts_output, annotation_output)
    else:
        provenance = materialize_bams(project, counts_output, annotation_output)
    manifest = {
        "schema_version": 2,
        "project_id": project.project_id,
        "analysis_set": project.analysis_set,
        "species": project.config["species"],
        "reference": project.config["reference"],
        "contrast_semantics": "positive effects are numerator minus denominator",
        "source": provenance,
        "canonical": {
            "counts": {"path": str(counts_output), "sha256": sha256(counts_output)},
            "samples": {"path": str(samples_output), "sha256": sha256(samples_output)},
            "annotation": {"path": str(annotation_output), "sha256": sha256(annotation_output)},
            "contrasts": {"path": str(contrasts_output), "sha256": sha256(contrasts_output)},
        },
    }
    manifest_output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Materialized {project.source_kind} input for {len(project.sample_rows)} samples")


if __name__ == "__main__":
    main()
