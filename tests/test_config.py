from __future__ import annotations

import csv
import shutil
from pathlib import Path

import pytest
import yaml

from bulk_rna_frame.config import ProjectValidationError, load_project, validation_report


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "bulk_rna_frame" / "templates" / "minimal"


def project_copy(tmp_path: Path) -> Path:
    destination = tmp_path / "project"
    shutil.copytree(TEMPLATE, destination)
    return destination / "project.yaml"


def test_template_is_valid_and_multi_contrast():
    project = load_project(TEMPLATE / "project.yaml")
    report = validation_report(project)
    assert report["status"] == "ok"
    assert report["samples"] == 9
    assert report["contrasts"] == ["treatment_a_vs_control", "treatment_b_vs_control"]


def test_count_samples_must_match_metadata(tmp_path):
    config = project_copy(tmp_path)
    counts = config.parent / "counts.tsv"
    lines = counts.read_text().splitlines()
    lines[0] = lines[0].replace("\ttreatment_b_3", "")
    lines[1:] = [line.rsplit("\t", 1)[0] for line in lines[1:]]
    counts.write_text("\n".join(lines) + "\n")
    with pytest.raises(ProjectValidationError, match="exactly match"):
        load_project(config)


def test_duplicate_contrast_ids_are_rejected(tmp_path):
    config = project_copy(tmp_path)
    path = config.parent / "contrasts.tsv"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    rows[1]["contrast_id"] = rows[0]["contrast_id"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ProjectValidationError, match="contrast_id values must be unique"):
        load_project(config)


def test_palette_must_cover_display_groups(tmp_path):
    config = project_copy(tmp_path)
    data = yaml.safe_load(config.read_text())
    del data["figures"]["palette"]["treatment_b"]
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="no colors"):
        load_project(config)


def test_pathways_requires_qc_and_de(tmp_path):
    config = project_copy(tmp_path)
    data = yaml.safe_load(config.read_text())
    data["modules"]["qc"] = False
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="requires modules.qc"):
        load_project(config)


def _counting_config(*, paired_end: bool = True) -> dict[str, object]:
    return {
        "threads": 2,
        "feature_type": "exon",
        "attribute": "gene_id",
        "paired_end": paired_end,
        "count_read_pairs": paired_end,
        "require_both_ends_aligned": paired_end,
        "exclude_chimeric_fragments": paired_end,
        "strandedness": "infer",
        "strand_test_modes": [0, 1, 2],
        "strand_min_dominance": 0.8,
    }


def _replace_with_bam_input(config: Path, kind: str) -> None:
    data = yaml.safe_load(config.read_text())
    sample_path = config.parent / "samples.tsv"
    rows = list(csv.DictReader(sample_path.open(), delimiter="\t"))
    fieldnames = [*rows[0], "bam"]
    bam_root = config.parent / "bams"
    bam_root.mkdir()
    for row in rows:
        row["bam"] = f"{row['sample_id']}.sorted.bam"
        (bam_root / row["bam"]).touch()
    with sample_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    gtf = config.parent / "genes.gtf"
    gtf.write_text('chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "gene1"; gene_name "Gene1";\n')
    if kind == "bam":
        data["inputs"] = {
            "kind": "bam",
            "bam_root": "bams",
            "samples": "samples.tsv",
            "gtf": "genes.gtf",
        }
    else:
        data["inputs"] = {
            "kind": "nfcore_rnaseq",
            "root": "bams",
            "samples": "samples.tsv",
            "gtf": "genes.gtf",
            "bam_pattern": "{sample_id}.sorted.bam",
        }
    data["counting"] = _counting_config()
    config.write_text(yaml.safe_dump(data, sort_keys=False))


@pytest.mark.parametrize("kind", ["bam", "nfcore_rnaseq"])
def test_bam_boundaries_resolve_declared_files(tmp_path, kind):
    config = project_copy(tmp_path)
    _replace_with_bam_input(config, kind)
    project = load_project(config)
    assert project.source_kind == kind
    assert len(project.bam_paths) == 9
    assert project.gtf == config.parent / "genes.gtf"


def test_nfcore_pattern_rejects_unknown_metadata_field(tmp_path):
    config = project_copy(tmp_path)
    _replace_with_bam_input(config, "nfcore_rnaseq")
    data = yaml.safe_load(config.read_text())
    data["inputs"]["bam_pattern"] = "{missing_field}.bam"
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="missing metadata column"):
        load_project(config)


def test_analysis_set_selects_a_subset_of_count_columns(tmp_path):
    config = project_copy(tmp_path)
    samples = config.parent / "samples.tsv"
    rows = list(csv.DictReader(samples.open(), delimiter="\t"))
    fieldnames = [*rows[0], "analysis_set"]
    for row in rows:
        row["analysis_set"] = "primary" if row["condition"] != "treatment_b" else "secondary"
    with samples.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    data = yaml.safe_load(config.read_text())
    data["inputs"]["analysis_set"] = "primary"
    data["contrasts"] = "primary_contrast.tsv"
    (config.parent / "primary_contrast.tsv").write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\n"
        "treatment_a_vs_control\tcondition\ttreatment_a\tcontrol\n"
    )
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    project = load_project(config)
    assert len(project.sample_rows) == 6
    assert {row["condition"] for row in project.sample_rows} == {"control", "treatment_a"}
