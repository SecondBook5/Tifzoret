from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from tifzoret.config import ProjectValidationError, load_project, validation_report


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"
MATERIALIZE = ROOT / "src" / "tifzoret" / "workflow" / "scripts" / "materialize_inputs.py"


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
            "snakemake", "--snakefile", str(ROOT / "src" / "tifzoret" / "workflow" / "Snakefile"),
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
            "snakemake", "--snakefile", str(ROOT / "src" / "tifzoret" / "workflow" / "Snakefile"),
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


# Companion documents (hypotheses claims/panels + figure recipe) may be authored
# either as sibling files or inlined in the main config, so a UI can emit one
# self-contained YAML. Both forms must resolve to identical loaded configuration,
# point at the canonical inputs directory, and materialize to faithful files.
CLAIMS_DOC = {
    "hypotheses": [
        {
            "id": "response",
            "statement": "Treatment changes contractile state.",
            "contrast": "treatment_a_vs_control",
            "expected_direction": "increased",
            "gene_panels": ["contractile"],
            "pathway_panels": ["response"],
        }
    ]
}
PANELS_DOC = {
    "gene_panels": {
        "contractile": {
            "description": "Contractile genes.",
            "groups": {"contractile": ["ACTA2", "CNN1", "MYH11"]},
        }
    },
    "program_annotations": {"ACTA2": "Contractile / cytoskeletal"},
    "program_colors": {"Contractile / cytoskeletal": "#D55E00"},
    "program_order": ["Contractile / cytoskeletal"],
    "gsea_programs": ["contractile"],
    "pathway_panels": {
        "response": {
            "description": "Response pathways.",
            "pathways": [{"collection": "custom", "pathway": "CONTRACTILE_PROGRAM"}],
        }
    },
}
RECIPE_DOC = {
    "figure_sets": {
        "primary": {
            "width": 12,
            "height": 6,
            "columns": 2,
            "panels": [
                {"id": "A", "source": "qc/figures/pca_correlation"},
                {
                    "id": "B",
                    "source": "contrasts/treatment_a_vs_control/analyses/composition/figures/cell_state_signatures",
                },
            ],
        }
    }
}


def _publication_project(root: Path, *, inline: bool) -> Path:
    shutil.copytree(TEMPLATE, root)
    config = root / "project.yaml"
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
    (root / "signatures.yaml").write_text(
        "signatures:\n  - id: contractile\n    label: Contractile\n    category: mural\n"
        "    genes: [ACTA2, CNN1, MYH11]\n"
    )
    if inline:
        data["hypotheses"] = {"claims": CLAIMS_DOC, "panels": PANELS_DOC}
        data["publication"] = {"recipe": RECIPE_DOC}
    else:
        data["hypotheses"] = {"claims": "hypotheses.yaml", "panels": "panels.yaml"}
        data["publication"] = {"recipe": "figure_recipe.yaml"}
        (root / "hypotheses.yaml").write_text(yaml.safe_dump(CLAIMS_DOC, sort_keys=False))
        (root / "panels.yaml").write_text(yaml.safe_dump(PANELS_DOC, sort_keys=False))
        (root / "figure_recipe.yaml").write_text(yaml.safe_dump(RECIPE_DOC, sort_keys=False))
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    return config


def test_inline_companion_documents_match_referenced_files(tmp_path):
    referenced = load_project(_publication_project(tmp_path / "referenced", inline=False))
    inline = load_project(_publication_project(tmp_path / "inline", inline=True))

    # The loaded companion documents are identical regardless of authoring form.
    assert inline.panel_config == referenced.panel_config == PANELS_DOC
    assert inline.hypothesis_config == referenced.hypothesis_config == CLAIMS_DOC
    assert inline.recipe_config == referenced.recipe_config == RECIPE_DOC

    # Both forms resolve companion paths to the canonical inputs directory.
    for project in (referenced, inline):
        inputs = project.result_root / "inputs"
        assert project.hypothesis_panels == inputs / "hypothesis_panels.yaml"
        assert project.hypotheses == inputs / "hypotheses.yaml"
        assert project.figure_recipe == inputs / "figure_recipe.yaml"


def test_materialize_stages_inline_companions(tmp_path):
    project = load_project(_publication_project(tmp_path / "inline", inline=True))
    inputs = project.result_root / "inputs"
    subprocess.run(
        [
            sys.executable, str(MATERIALIZE),
            "--project-config", str(project.config_path),
            "--counts", str(inputs / "counts.tsv"),
            "--samples", str(inputs / "samples.tsv"),
            "--annotation", str(inputs / "annotation.tsv"),
            "--contrasts", str(inputs / "contrasts.tsv"),
            "--manifest", str(inputs / "input_manifest.json"),
            "--panels", str(project.hypothesis_panels),
            "--claims", str(project.hypotheses),
            "--recipe", str(project.figure_recipe),
        ],
        check=True,
    )
    assert yaml.safe_load(project.hypothesis_panels.read_text()) == PANELS_DOC
    assert yaml.safe_load(project.hypotheses.read_text()) == CLAIMS_DOC
    assert yaml.safe_load(project.figure_recipe.read_text()) == RECIPE_DOC


def _dry_run_rules(config: Path) -> str:
    """Return `snakemake --dry-run` stdout for a project's manifest target."""
    project = load_project(config)
    proc = subprocess.run(
        [
            "snakemake", "--snakefile",
            str(ROOT / "src" / "tifzoret" / "workflow" / "Snakefile"),
            "--configfile", str(config), "--cores", "1", "--dry-run",
            str(project.result_root / "manifest.json"),
        ],
        cwd=config.parent, check=True, capture_output=True, text=True,
    )
    return proc.stdout


def test_opt_in_modules_wire_into_the_dag(tmp_path):
    """batch, de_confirm, deconvolution, and curvature belong to no profile and
    appear in the DAG only when toggled on under analysis.modules (curvature
    drags in its wgcna prerequisite)."""
    config = project_copy(tmp_path)
    data = yaml.safe_load(config.read_text())
    data["analysis"]["batch"] = "batch"
    data["analysis"]["modules"] = {
        "wgcna": True, "batch": True, "de_confirm": True,
        "deconvolution": True, "curvature": True,
    }
    data["analysis"]["settings"] = {
        "wgcna": {
            "top_variable_genes": 100, "minimum_module_size": 5,
            "minimum_recommended_samples": 5, "network_neighbors": 5,
        },
        "de": {"shrinkage": "ashr", "confirm_method": "edger"},
        "deconvolution": {"method": "nnls", "min_genes": 3},
        "curvature": {"alpha": 0.5, "top_bridges": 10},
    }
    data["resources"]["deconvolution_signature"] = "signature.tsv"
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    (config.parent / "signature.tsv").write_text(
        "gene\tMural\tEpithelial\tImmune\n"
        "Acta2\t8.0\t0.1\t0.2\nCnn1\t7.5\t0.1\t0.1\nMyh11\t9.0\t0.0\t0.1\n"
        "Epcam\t0.1\t8.0\t0.2\nKrt5\t0.1\t7.0\t0.1\nKrt18\t0.2\t7.5\t0.1\n"
        "Lst1\t0.1\t0.2\t8.0\nTyrobp\t0.0\t0.1\t7.5\nCtss\t0.1\t0.2\t7.0\n"
    )
    project = load_project(config)
    assert {"batch", "de_confirm", "deconvolution", "curvature", "wgcna"}.issubset(project.modules)
    stdout = _dry_run_rules(config)
    for rule in ("study_batch", "study_deconvolution", "contrast_de_confirm",
                 "contrast_curvature", "contrast_wgcna"):
        assert rule in stdout


def test_omnibus_contrast_type_is_accepted(tmp_path):
    """A `type: omnibus` contrast row carries no numerator/denominator and
    resolves to a distinct omnibus contrast id."""
    config = project_copy(tmp_path)
    contrasts = config.parent / "contrasts.tsv"
    contrasts.write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\ttype\treduced\n"
        "treatment_a_vs_control\tcondition\ttreatment_a\tcontrol\tpairwise\t\n"
        "condition_any\tcondition\t\t\tomnibus\t~ batch\n"
    )
    project = load_project(config)
    ids = [row["contrast_id"] for row in project.contrast_rows]
    assert "condition_any" in ids


def test_coefficient_contrast_type_is_accepted(tmp_path):
    """A `type: coefficient` row names a resultsNames() coefficient under a
    (possibly per-row) design; factor/numerator/denominator are labels only and
    the coefficient itself is checked at DE time against the fitted model."""
    config = project_copy(tmp_path)
    (config.parent / "contrasts.tsv").write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\ttype\tcoefficient\tdesign\n"
        "interaction_ab\tcondition\ttreatment_a\tcontrol\tcoefficient"
        "\tconditiontreatment_a.batchb2\t~ batch + condition\n"
    )
    project = load_project(config)
    assert "interaction_ab" in [row["contrast_id"] for row in project.contrast_rows]


def test_coefficient_contrast_requires_a_coefficient(tmp_path):
    """The named coefficient is the whole point of a coefficient contrast; an
    empty one must be rejected rather than silently degrade to a pairwise test."""
    config = project_copy(tmp_path)
    (config.parent / "contrasts.tsv").write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\ttype\tcoefficient\n"
        "interaction_ab\tcondition\ttreatment_a\tcontrol\tcoefficient\t\n"
    )
    with pytest.raises(ProjectValidationError, match="requires a non-empty coefficient"):
        load_project(config)


def test_coefficient_contrast_rejects_unknown_design_variable(tmp_path):
    """A per-row design must reference only columns present in samples.tsv."""
    config = project_copy(tmp_path)
    (config.parent / "contrasts.tsv").write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\ttype\tcoefficient\tdesign\n"
        "interaction_ab\tcondition\ttreatment_a\tcontrol\tcoefficient\tsomeCoef\t~ genotype\n"
    )
    with pytest.raises(ProjectValidationError, match="design variable 'genotype' is absent"):
        load_project(config)


def test_coefficient_contrast_validates_reference_levels(tmp_path):
    """`reference_levels` entries must be factor=level with a level that exists
    in the named design factor."""
    config = project_copy(tmp_path)
    (config.parent / "contrasts.tsv").write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\ttype\tcoefficient\tdesign\treference_levels\n"
        "interaction_ab\tcondition\ttreatment_a\tcontrol\tcoefficient\tsomeCoef"
        "\t~ batch + condition\tcondition=nonesuch\n"
    )
    with pytest.raises(ProjectValidationError, match="level 'nonesuch' is absent from condition"):
        load_project(config)


@pytest.mark.parametrize(
    "reduced,expected",
    [
        ("", "requires a non-empty reduced design formula"),
        ("~ batch + condition", "reduced formula must not contain the tested factor"),
    ],
)
def test_omnibus_reduced_formula_is_validated(tmp_path, reduced, expected):
    """An omnibus LRT needs a reduced design that drops the tested factor: a
    missing reduced formula, or one that still contains the factor, is rejected."""
    config = project_copy(tmp_path)
    (config.parent / "contrasts.tsv").write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\ttype\treduced\n"
        f"condition_any\tcondition\t\t\tomnibus\t{reduced}\n"
    )
    with pytest.raises(ProjectValidationError, match=expected):
        load_project(config)


def test_shrinkage_rejects_an_unknown_prior(tmp_path):
    """DE shrinkage is a fixed enum (apeglm/ashr/normal/none); a typo must fail
    schema validation rather than silently fall through to the default."""
    config = project_copy(tmp_path)
    data = yaml.safe_load(config.read_text())
    data["analysis"]["settings"] = {"de": {"shrinkage": "bogus"}}
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="shrinkage"):
        load_project(config)


def test_go_domains_are_opt_in_and_bounded(tmp_path):
    """CC/MF are strictly opt-in: with `resources.go_domains` unset the resolved
    config carries no domain list (the GO-BP default is applied downstream), so
    legacy GO-BP studies gain no cellular-component/molecular-function content
    unless explicitly requested. The domain vocabulary is also bounded — an
    unknown domain fails schema validation rather than resolving silently.
    (The faceted `ontology_domains` figure is emitted unconditionally as an
    additive artifact; the BP-only guarantee is about content, not its presence.)"""
    project = load_project(project_copy(tmp_path))  # standard profile, no go_domains set
    assert project.config["resources"].get("go_domains") is None

    config = project_copy(tmp_path / "unknown")
    data = yaml.safe_load(config.read_text())
    data.setdefault("resources", {})["go_domains"] = ["BP", "XYZ"]
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="go_domains"):
        load_project(config)


def test_counting_boundary_emits_tpm_fpkm(tmp_path):
    """BAM (counting) boundaries add TPM/FPKM and per-gene exon-length outputs
    to the materialize rule; count-matrix inputs do not."""
    config = project_copy(tmp_path)
    _replace_with_bam_input(config, "bam")
    stdout = _dry_run_rules(config)
    for token in ("tpm.tsv", "fpkm.tsv", "gene_lengths.tsv"):
        assert token in stdout


def test_wave2_opt_in_modules_wire_into_the_dag(tmp_path):
    """consensus, spia, variance_partition, and enrichment_map belong to no
    profile and appear in the DAG only when toggled on; enabling extra GO
    domains and Reactome still validates, and the additive ontology-domain
    outputs are produced (the front door requires them)."""
    config = project_copy(tmp_path)
    data = yaml.safe_load(config.read_text())
    data["analysis"]["batch"] = "batch"
    data["analysis"]["modules"] = {
        "consensus": True, "spia": True,
        "variance_partition": True, "enrichment_map": True,
    }
    data["analysis"]["settings"] = {
        "consensus": {"min_contrasts": 2},
        "spia": {"top_pathways": 15},
        "variance_partition": {"top_variable_genes": 200},
        "enrichment_map": {"min_similarity": 0.25, "top_terms": 40},
    }
    data["resources"]["go_domains"] = ["BP", "CC", "MF"]
    data["resources"]["providers"] = {"kegg": True, "reactome": True}
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    project = load_project(config)
    assert {"consensus", "spia", "variance_partition", "enrichment_map"}.issubset(project.modules)
    assert project.config["resources"]["go_domains"] == ["BP", "CC", "MF"]
    stdout = _dry_run_rules(config)
    for rule in ("study_consensus", "contrast_spia",
                 "study_variance_partition", "contrast_enrichment_map"):
        assert rule in stdout
    # The GO CC/MF + Reactome breadth is emitted as an additive faceted view
    # whose outputs the front-door index collects; a missing producer here is
    # exactly the class of DAG break this test guards.
    assert "ontology_domains.pdf" in stdout
    assert "ontology_domain_displayed.tsv" in stdout


def test_consensus_requires_two_pairwise_contrasts(tmp_path):
    """A single pairwise contrast cannot form a cross-contrast consensus."""
    config = project_copy(tmp_path)
    (config.parent / "contrasts.tsv").write_text(
        "contrast_id\tfactor\tnumerator\tdenominator\n"
        "treatment_a_vs_control\tcondition\ttreatment_a\tcontrol\n"
    )
    data = yaml.safe_load(config.read_text())
    data["analysis"]["modules"] = {"consensus": True}
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="at least two pairwise contrasts"):
        load_project(config)


def test_spia_requires_the_kegg_provider(tmp_path):
    """SPIA perturbs KEGG pathway topology, so it needs the KEGG provider."""
    config = project_copy(tmp_path)
    data = yaml.safe_load(config.read_text())
    data["analysis"]["modules"] = {"spia": True}
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="resources.providers.kegg"):
        load_project(config)


def test_deconvolution_preset_resolves_and_is_exclusive(tmp_path):
    """A named packaged preset resolves to its shipped signature matrix, is
    mutually exclusive with a supplied signature, and an unknown name fails
    with the installed presets listed."""
    config = project_copy(tmp_path)
    data = yaml.safe_load(config.read_text())
    data["analysis"]["modules"] = {"deconvolution": True}
    data["analysis"]["settings"] = {"deconvolution": {"method": "nnls", "min_genes": 3}}
    data["resources"]["deconvolution_preset"] = "mouse_immune"
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    project = load_project(config)
    assert "deconvolution" in project.modules
    assert project.deconvolution_signature.name == "mouse_immune.tsv"
    assert project.deconvolution_signature.is_file()

    data["resources"]["deconvolution_signature"] = "signature.tsv"
    (config.parent / "signature.tsv").write_text("gene\tA\tB\nActa2\t1\t0\nEpcam\t0\t1\n")
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="mutually exclusive"):
        load_project(config)

    del data["resources"]["deconvolution_signature"]
    data["resources"]["deconvolution_preset"] = "does_not_exist"
    config.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(ProjectValidationError, match="available presets"):
        load_project(config)
