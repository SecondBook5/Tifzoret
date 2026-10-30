"""validate fails early when BAMs and the GTF name their chromosomes differently."""

from __future__ import annotations

import gzip
import struct
from pathlib import Path

import pytest

from tifzoret.config import ProjectValidationError, load_project
from tifzoret.config.contigs import contig_mismatch_message, read_bam_contigs, sample_gtf_seqnames

from test_config import _replace_with_bam_input, project_copy


def write_bam_header(path: Path, names: list[str]) -> None:
    """A BAM with a header and no reads: enough for the contig list."""
    text = b"@HD\tVN:1.6\n" + b"".join(f"@SQ\tSN:{n}\tLN:1000\n".encode() for n in names)
    body = b"BAM\x01" + struct.pack("<i", len(text)) + text + struct.pack("<i", len(names))
    for name in names:
        encoded = name.encode() + b"\0"
        body += struct.pack("<i", len(encoded)) + encoded + struct.pack("<i", 1000)
    with gzip.open(path, "wb") as handle:
        handle.write(body)


def write_gtf(path: Path, seqname: str) -> None:
    path.write_text(f'{seqname}\ttest\texon\t1\t10\t.\t+\t.\tgene_id "g1"; gene_name "G1";\n')


def test_reads_contig_names_from_a_bam_header(tmp_path):
    bam = tmp_path / "a.bam"
    write_bam_header(bam, ["chr1", "chr2", "chrM"])
    assert read_bam_contigs(bam) == ["chr1", "chr2", "chrM"]


def test_an_unreadable_or_empty_bam_is_skipped_not_an_error(tmp_path):
    empty = tmp_path / "empty.bam"
    empty.touch()
    text = tmp_path / "text.bam"
    text.write_text("not a bam")
    gtf = tmp_path / "g.gtf"
    write_gtf(gtf, "1")
    assert read_bam_contigs(empty) is None and read_bam_contigs(text) is None
    assert contig_mismatch_message(empty, gtf) is None
    assert contig_mismatch_message(text, gtf) is None


def test_matching_names_pass_and_the_message_names_a_chr_prefix_mismatch(tmp_path):
    bam = tmp_path / "a.bam"
    write_bam_header(bam, ["chr1", "chr2"])
    ucsc, ensembl = tmp_path / "ucsc.gtf", tmp_path / "ensembl.gtf"
    write_gtf(ucsc, "chr1")
    write_gtf(ensembl, "1")
    assert contig_mismatch_message(bam, ucsc) is None
    message = contig_mismatch_message(bam, ensembl)
    assert message and "share no chromosome names" in message and "'chr' prefix" in message


def test_unrelated_references_get_the_generic_hint(tmp_path):
    bam = tmp_path / "a.bam"
    write_bam_header(bam, ["chr1"])
    gtf = tmp_path / "g.gtf"
    write_gtf(gtf, "scaffold_9")
    message = contig_mismatch_message(bam, gtf)
    assert message and "'chr' prefix" not in message and "same reference" in message


def test_scaffolds_with_identical_names_do_not_hide_a_chromosome_mismatch(tmp_path):
    """mm39 BAMs list unplaced scaffolds (JH584295.1) that Ensembl and UCSC GTFs name identically."""
    bam = tmp_path / "a.bam"
    write_bam_header(bam, ["chr1", "chr2", "chrM", "JH584295.1"])
    gtf = tmp_path / "ensembl.gtf"
    gtf.write_text(
        '1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "a"; gene_name "A";\n'
        'MT\ttest\texon\t1\t10\t.\t+\t.\tgene_id "b"; gene_name "B";\n'
        'JH584295.1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "c"; gene_name "C";\n'
    )
    message = contig_mismatch_message(bam, gtf)
    assert message and "share no chromosome names" in message and "'chr' prefix" in message


def test_gtf_sampling_reads_the_tail_of_a_large_file(tmp_path):
    gtf = tmp_path / "big.gtf"
    filler = 'chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "g"; gene_name "G";\n'
    gtf.write_text(filler * 30000 + 'chrX\ttest\texon\t1\t10\t.\t+\t.\tgene_id "x"; gene_name "X";\n' * 40000)
    assert sample_gtf_seqnames(gtf) == {"chr1", "chrX"}


@pytest.mark.parametrize("kind", ["bam", "nfcore_rnaseq"])
def test_validate_rejects_a_chr_prefix_mismatch_before_any_counting(tmp_path, kind):
    config = project_copy(tmp_path)
    _replace_with_bam_input(config, kind)
    for bam in (config.parent / "bams").iterdir():
        write_bam_header(bam, ["chr1", "chr2"])
    write_gtf(config.parent / "genes.gtf", "1")
    with pytest.raises(ProjectValidationError) as error:
        load_project(config)
    assert "share no chromosome names" in str(error.value)


def test_validate_accepts_matching_names(tmp_path):
    config = project_copy(tmp_path)
    _replace_with_bam_input(config, "bam")
    for bam in (config.parent / "bams").iterdir():
        write_bam_header(bam, ["chr1", "chr2"])
    write_gtf(config.parent / "genes.gtf", "chr1")
    assert len(load_project(config).bam_paths) == 9
