from __future__ import annotations

import shutil
from pathlib import Path

import yaml

from tifzoret.cli import main


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"


def test_collection_validate_and_run_signed_meta_analysis(tmp_path: Path):
    projects = []
    for index, effect in ((1, 1.0), (2, 1.5), (3, 0.8)):
        destination = tmp_path / f"study_{index}"
        shutil.copytree(TEMPLATE, destination)
        config_path = destination / "project.yaml"
        config = yaml.safe_load(config_path.read_text())
        config["project"]["id"] = f"study_{index}"
        config_path.write_text(yaml.safe_dump(config, sort_keys=False))
        result = destination / "results" / f"study_{index}" / "all" / "contrasts" / "treatment_a_vs_control" / "analyses" / "de" / "tables"
        result.mkdir(parents=True)
        (result / "de_results.tsv").write_text(
            "gene_id\tgene_symbol\tlog2_fold_change\tp_value\n"
            f"g1\tGene1\t{effect}\t0.01\n"
        )
        projects.append(config_path)

    collection = tmp_path / "collection.yaml"
    collection.write_text(yaml.safe_dump({
        "version": 1,
        "collection": {"id": "three_studies", "title": "Three studies"},
        "studies": [
            {"id": f"study_{index}", "project": str(path), "contrast": "treatment_a_vs_control"}
            for index, path in enumerate(projects, start=1)
        ],
        "methods": {"stouffer": True, "leave_one_out": True, "bh_correction": True},
        "output": {"root": "collection_results"},
    }, sort_keys=False))
    assert main(["collection", "validate", str(collection)]) == 0
    assert main(["collection", "run", str(collection)]) == 0
    output = tmp_path / "collection_results" / "three_studies" / "meta_analysis.tsv"
    text = output.read_text()
    assert "leave_one_out_direction_stable" in text
    assert "positive" in text
