"""Early check that BAM and GTF name their chromosomes the same way.

A BAM aligned to UCSC-style names (``chr1``) counted against an Ensembl-style GTF (``1``) assigns
no reads, and the mismatch otherwise surfaces only after hours of counting. The check reads the
contig names from the first BAM header (no samtools needed) and samples sequence names from the
head and tail of the GTF, which is enough to tell the naming styles apart.
"""

from __future__ import annotations

import gzip
import re
import struct
from pathlib import Path

# Main chromosomes in either naming style (1, chr1, X, chrX, M, MT, chrM). Unplaced scaffolds such as
# JH584295.1 carry the same name in UCSC- and Ensembl-style files, so they must not decide the check.
_MAIN_CHROMOSOME = re.compile(r"^(?:chr)?(?:\d{1,2}|[XYWZ]|M|MT)$", re.IGNORECASE)
_GTF_HEAD_LINES = 20000
_GTF_TAIL_BYTES = 2_000_000


def read_bam_contigs(path: Path) -> list[str] | None:
    """Reference names from a BAM header, or None when the file cannot be read as a BAM."""
    try:
        with gzip.open(path, "rb") as handle:
            if handle.read(4) != b"BAM\x01":
                return None
            (text_length,) = struct.unpack("<i", handle.read(4))
            handle.read(text_length)
            (n_ref,) = struct.unpack("<i", handle.read(4))
            names: list[str] = []
            for _ in range(n_ref):
                (name_length,) = struct.unpack("<i", handle.read(4))
                names.append(handle.read(name_length).rstrip(b"\0").decode("ascii", "replace"))
                handle.read(4)
            return names
    except (OSError, EOFError, struct.error):
        return None


def sample_gtf_seqnames(path: Path) -> set[str]:
    """Sequence names from the first lines (and, for plain text, the last 2 MB) of a GTF."""
    names: set[str] = set()
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            seen = 0
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                names.add(line.split("\t", 1)[0])
                seen += 1
                if seen >= _GTF_HEAD_LINES:
                    break
        if path.suffix != ".gz":
            size = path.stat().st_size
            if size > _GTF_TAIL_BYTES:
                with open(path, "rb") as handle:
                    handle.seek(size - _GTF_TAIL_BYTES)
                    handle.readline()  # drop the partial first line
                    for raw in handle:
                        line = raw.decode("utf-8", "replace")
                        if line.strip() and not line.startswith("#"):
                            names.add(line.split("\t", 1)[0])
    except OSError:
        return set()
    return names


def contig_mismatch_message(bam: Path, gtf: Path) -> str | None:
    """An error message when the BAM and GTF share no sequence names, else None.

    Returns None whenever the check cannot be made (unreadable BAM or GTF), so it never blocks a
    project whose inputs are checked elsewhere."""
    bam_names = read_bam_contigs(bam)
    gtf_names = sample_gtf_seqnames(gtf)
    if not bam_names or not gtf_names:
        return None
    bam_main = {name for name in bam_names if _MAIN_CHROMOSOME.match(name)}
    gtf_main = {name for name in gtf_names if _MAIN_CHROMOSOME.match(name)}
    if bam_main and gtf_main:
        if bam_main & gtf_main:
            return None
        bam_names, gtf_names = sorted(bam_main), gtf_main
    elif set(bam_names) & gtf_names:
        return None

    def strip(values):
        return {value[3:] if value.lower().startswith("chr") else value for value in values}

    hint = (
        " They differ only in a 'chr' prefix: add or remove it in the GTF (or use an annotation "
        "built for the same reference as the BAMs)."
        if strip(bam_names) & strip(gtf_names)
        else " Use an annotation built for the same reference as the BAMs."
    )

    def example(values):
        return ", ".join(sorted(values, key=lambda v: (len(v), v))[:3])

    return (
        f"BAM and GTF share no chromosome names, so featureCounts would assign no reads: "
        f"{bam.name} uses {example(bam_names)} ... but {gtf.name} uses {example(gtf_names)} ...." + hint
    )
