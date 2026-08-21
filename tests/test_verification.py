import csv
import shutil
from pathlib import Path

from tifzoret.config import load_project
from tifzoret.verification import verify_project, verify_runs


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"


def test_verification_detects_and_reports_table_differences(tmp_path: Path):
    reference = tmp_path / "reference"; candidate = tmp_path / "candidate"
    reference.mkdir(); candidate.mkdir()
    (reference / "counts.tsv").write_text("gene_id\ts1\ng1\t10\n")
    (candidate / "counts.tsv").write_text("gene_id\ts1\ng1\t11\n")
    result = verify_runs(reference, candidate)
    assert not result.passed
    assert result.report["failed_tables"] == ["counts.tsv"]


def test_legacy_verification_compares_counts_de_decisions_and_figures(tmp_path: Path):
    study = tmp_path / "study"
    shutil.copytree(TEMPLATE, study)
    samples_path = study / "samples.tsv"
    with samples_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = [*rows[0], "bam"]
    for row in rows:
        row["bam"] = f"{row['sample_id']}.bam"
    with samples_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    project = load_project(study / "project.yaml")
    candidate = project.result_root
    inputs = candidate / "inputs"
    inputs.mkdir(parents=True)
    shutil.copy2(study / "counts.tsv", inputs / "counts.tsv")
    (candidate / "figures").mkdir(); (candidate / "tables").mkdir()
    (candidate / "figures" / "index.json").write_text('{"figures": []}\n')
    (candidate / "tables" / "index.json").write_text('{"tables": []}\n')

    with (study / "counts.tsv").open(newline="") as handle:
        counts = list(csv.DictReader(handle, delimiter="\t"))
    reference = tmp_path / "reference"
    count_cache = reference / ".cache" / "counts"
    de_cache = reference / ".cache" / "de"
    count_cache.mkdir(parents=True); de_cache.mkdir(parents=True)
    sample_ids = [row["sample_id"] for row in rows]
    with (count_cache / "gene_counts.tsv").open("w", newline="") as handle:
        handle.write("# featureCounts fixture\n")
        writer = csv.DictWriter(handle, fieldnames=["Geneid", "Chr", "Start", "End", "Strand", "Length", *[f"/old/{sample}.bam" for sample in sample_ids]], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in counts:
            writer.writerow({
                "Geneid": row["gene_id"], "Chr": "1", "Start": "1", "End": "2", "Strand": "+", "Length": "2",
                **{f"/old/{sample}.bam": row[sample] for sample in sample_ids},
            })
    (de_cache / "de_raw.tsv").write_text(
        "gene_id\tbaseMean\tlog2FoldChange\tstat\tpvalue\tpadj\n"
        "g1\t10\t1.2\t3\t0.01\t0.02\n"
    )
    (de_cache / "de_shrunken.tsv").write_text(
        "gene_id\tbaseMean\tlog2FoldChange\tpvalue\tpadj\n"
        "g1\t10\t1.1\t0.01\t0.02\n"
    )
    for contrast in project.contrast_rows:
        table = candidate / "contrasts" / contrast["contrast_id"] / "analyses" / "de" / "tables"
        table.mkdir(parents=True)
        (table / "de_results.tsv").write_text(
            "gene_id\tbase_mean\tlog2_fold_change_raw\tlog2_fold_change\tstatistic\tp_value\tadjusted_p_value\n"
            "g1\t10\t1.2\t1.1\t3\t0.01\t0.02\n"
        )
    result = verify_project(project, reference)
    assert result.passed
    assert result.report["mode"] == "legacy_migration"
    assert result.report["counts"]["passed"]
    assert result.report["de"]["passed"]
