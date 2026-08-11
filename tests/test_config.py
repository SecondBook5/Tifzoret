from __future__ import annotations

import csv
import shutil
import subprocess
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
    data["analysis"].setdefault("modules", {})["qc"] = False
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
    data["analysis"]["contrasts"] = "primary_contrast.tsv"
    (config.parent / "primary_contrast.tsv").write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\n"
        "treatment_a_vs_control\tcondition\ttreatment_a\tcontrol\n"
    )
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    project = load_project(config)
    assert len(project.sample_rows) == 6
    assert {row["condition"] for row in project.sample_rows} == {"control", "treatment_a"}


def test_publication_profile_validates_cross_file_contracts(tmp_path):
    config = project_copy(tmp_path)
    data = yaml.safe_load(config.read_text())
    data["species"] = {"provider": "mouse", "scientific_name": "Mus musculus", "taxonomy_id": 10090}
    data["reference"] = {"genome_build": "GRCm39", "annotation_release": 107}
    data["analysis"]["profile"] = "publication"
    data["analysis"]["random_seed"] = 1
    data["analysis"]["settings"] = {
        "composition": {"min_genes": 2},
        "regulators": {"min_targets": 2, "top_regulators": 5},
        "networks": {"required_score": 700, "max_nodes": 40, "seed": 1},
    }
    data["resources"]["cell_state_signatures"] = "signatures.yaml"
    data["resources"]["providers"] = {"dorothea": True, "string": True}
    data["hypotheses"] = {"claims": "hypotheses.yaml", "panels": "panels.yaml"}
    data["publication"] = {"recipe": "figure_recipe.yaml"}
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    (config.parent / "signatures.yaml").write_text(
        "signatures:\n  - id: contractile\n    label: Contractile\n    category: mural\n"
        "    genes: [ACTA2, CNN1, MYH11]\n"
    )
    (config.parent / "hypotheses.yaml").write_text(
        "hypotheses:\n  - id: response\n    statement: Treatment changes contractile state.\n"
        "    contrast: treatment_a_vs_control\n    expected_direction: increased\n"
        "    gene_panels: [contractile]\n    pathway_panels: [response]\n"
    )
    (config.parent / "panels.yaml").write_text(
        "gene_panels:\n  contractile:\n    description: Contractile genes.\n"
        "    groups:\n      contractile: [ACTA2, CNN1, MYH11]\n"
        "program_annotations:\n  ACTA2: Contractile / cytoskeletal\n"
        "program_colors:\n  Contractile / cytoskeletal: '#D55E00'\n"
        "program_order: [Contractile / cytoskeletal]\n"
        "gsea_programs: [contractile]\n"
        "pathway_panels:\n  response:\n    description: Response pathways.\n"
        "    pathways:\n      - {collection: custom, pathway: CONTRACTILE_PROGRAM}\n"
    )
    (config.parent / "figure_recipe.yaml").write_text(
        "figure_sets:\n  primary:\n    width: 12\n    height: 6\n    columns: 2\n"
        "    panels:\n      - {id: A, source: qc/figures/pca_correlation}\n"
        "      - {id: B, source: contrasts/treatment_a_vs_control/analyses/composition/figures/cell_state_signatures}\n"
    )
    project = load_project(config)
    assert project.config["analysis"]["profile"] == "publication"
    assert {"composition", "regulators", "networks", "hypotheses", "publication"}.issubset(project.modules)
    assert project.recipe_config["figure_sets"]["primary"]["panels"][0]["id"] == "A"
    assert project.panel_config["program_colors"]["Contractile / cytoskeletal"] == "#D55E00"
    assert project.panel_config["gsea_programs"] == ["contractile"]
    publication_dry_run = subprocess.run(
        [
            "snakemake", "--snakefile", str(ROOT / "src" / "bulk_rna_frame" / "workflow" / "Snakefile"),
            "--configfile", str(config), "--cores", "1", "--dry-run",
            str(project.result_root / "manifest.json"),
        ],
        cwd=config.parent,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "contrast_publication" in publication_dry_run.stdout
    assert "assemble_figure" in publication_dry_run.stdout
    data["analysis"]["profile"] = "full"
    data["analysis"]["settings"].update({
        "sva": {"minimum_recommended_samples": 10},
        "wgcna": {
            "top_variable_genes": 100,
            "minimum_module_size": 5,
            "minimum_recommended_samples": 15,
        },
        "mediation": {
            "mediator_pathway": "CONTRACTILE_PROGRAM",
            "outcome_pathways": ["IMMUNE_PROGRAM"],
            "simulations": 100,
            "minimum_recommended_samples": 20,
        },
    })
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    full = load_project(config)
    assert {"sva", "wgcna", "mediation", "multilayer"}.issubset(full.modules)
    full_dry_run = subprocess.run(
        [
            "snakemake", "--snakefile", str(ROOT / "src" / "bulk_rna_frame" / "workflow" / "Snakefile"),
            "--configfile", str(config), "--cores", "1", "--dry-run",
            str(full.result_root / "manifest.json"),
        ],
        cwd=config.parent,
        check=True,
        capture_output=True,
        text=True,
    )
    for rule in ("contrast_sva", "contrast_wgcna", "contrast_mediation", "contrast_multilayer"):
        assert rule in full_dry_run.stdout


def test_human_provider_contract_uses_same_module_interfaces(tmp_path):
    config = project_copy(tmp_path)
    data = yaml.safe_load(config.read_text())
    data["species"] = {"provider": "human", "scientific_name": "Homo sapiens", "taxonomy_id": 9606}
    data["reference"] = {"genome_build": "GRCh38", "annotation_release": 110}
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    project = load_project(config)
    assert project.config["species"]["taxonomy_id"] == 9606
    assert project.modules == tuple(sorted(project.modules))


def test_gtrd_provider_requires_an_explicit_snapshot(tmp_path):
    config = project_copy(tmp_path)
    data = yaml.safe_load(config.read_text())
    data["species"] = {"provider": "mouse", "scientific_name": "Mus musculus", "taxonomy_id": 10090}
    data["reference"] = {"genome_build": "GRCm39"}
    data["analysis"]["modules"] = {"regulators": True}
    data["resources"]["providers"] = {"gtrd": True}
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="exported GTRD-derived"):
        load_project(config)

    (config.parent / "gtrd_edges.tsv").write_text("source\ttarget\nTF1\tGene1\n")
    data["resources"]["regulon_edges"] = "gtrd_edges.tsv"
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    project = load_project(config)
    assert project.regulon_edges == config.parent / "gtrd_edges.tsv"
