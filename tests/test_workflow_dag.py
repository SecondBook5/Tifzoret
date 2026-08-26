"""Snakemake dry-run integration tests: prove the rule graph assembles.

A dry-run (``snakemake -n``) builds the full job DAG without executing any rule,
so it needs neither R, conda, nor the network -- only the packaged Snakefile and
a valid project. It is the cheapest test that exercises rule *wiring*: every
output is reachable, no rule is ambiguous, and the graph is acyclic. The unit
tests elsewhere cover the individual scripts in isolation; these cover how the
Snakefile stitches them together.

The publication profile is the important case. Its rule graph -- composition,
the dual-view regulators stage (including the unsigned binding view), the GRN
and STRING figure seams, hypotheses, the publication panels, and the figure
assembler -- is otherwise unexercised by CI, which only runs the standard
profile end to end. A dropped or mis-wired publication rule would resolve to a
DAG error here (non-zero exit) before it ever reached a study run.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import yaml

from tifzoret.figures import PANEL_REGISTRY


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "src" / "tifzoret" / "templates" / "minimal"

# Rules that only the publication profile pulls in. Presence of each in the
# planned DAG is what distinguishes "publication graph wired" from "core graph
# happened to resolve". assemble_figure is the capstone: it depends on every
# recipe panel's source, so its presence also proves each panel -- including the
# binding-view regulator panel -- resolved to a producible file.
PUBLICATION_RULES = frozenset({
    "contrast_composition",
    "contrast_regulators",
    "contrast_grn",
    "contrast_grn_radial",
    "contrast_networks",
    "contrast_string_figures",
    "contrast_hypotheses",
    "contrast_publication",
    "assemble_figure",
})


def _dry_run(config_path: Path) -> subprocess.CompletedProcess:
    """Plan the DAG for ``config_path`` offline (no conda, single core)."""
    return subprocess.run(
        [
            sys.executable, "-m", "tifzoret", "dry-run",
            str(config_path), "--no-conda", "--cores", "1",
        ],
        check=False,
        capture_output=True,
        text=True,
    )


def _planned_rules(output: str) -> set[str]:
    """Rule names Snakemake reports it would run (one ``rule <name>:`` per job)."""
    rules = set()
    for line in output.splitlines():
        stripped = line.strip()
        if stripped.startswith("rule ") and stripped.endswith(":"):
            rules.add(stripped[len("rule "):-1])
    return rules


def _result_tail(rel: str) -> str:
    """A result-relative panel path, stripped of its ``contrasts/<id>/`` prefix.

    Panel paths in the figure catalog carry a ``{contrast}`` wildcard; a planned
    DAG names concrete contrasts. Dropping the ``contrasts/{contrast}/`` prefix
    yields a contrast-agnostic tail (e.g. ``analyses/publication/tables/x.tsv``)
    that is a substring of every concrete absolute path Snakemake prints for it,
    so a simple ``tail in dry_run_stdout`` membership test works regardless of
    which contrast produced it. QC panel paths carry no prefix and pass through.
    """
    prefix = "contrasts/{contrast}/"
    return rel[len(prefix):] if rel.startswith(prefix) else rel


def _publication_project(tmp_path: Path) -> Path:
    """A publication-profile project built on the minimal counts template.

    Adds the synthetic companion documents the publication modules require --
    cell-state signatures, hypothesis programs and claims, a figure recipe, and
    an unsigned binding-prior regulon -- plus the STRING/DoRothEA providers and a
    taxonomy id. Content is arbitrary but schema-valid; a dry-run never reads it,
    so only the wiring it induces is under test. The recipe deliberately requests
    both regulator-activity views (``default`` and ``binding``) so the dual-view
    seam is part of the assembled figure.
    """
    destination = tmp_path / "project"
    shutil.copytree(TEMPLATE, destination)

    yaml.safe_dump(
        {"signatures": [
            {"id": "smooth_muscle", "label": "Smooth muscle", "category": "structural",
             "description": "Synthetic smooth-muscle signature for the DAG fixture.",
             "genes": ["Acta2", "Cnn1", "Myh11", "Tagln", "Tpm2", "Mylk"]},
            {"id": "epithelial", "label": "Epithelial", "category": "structural",
             "description": "Synthetic epithelial signature for the DAG fixture.",
             "genes": ["Epcam", "Krt5", "Krt18", "Cldn7", "Gata3", "Foxq1"]},
        ]},
        (destination / "signatures.yaml").open("w"), sort_keys=False,
    )

    (destination / "binding.tsv").write_text(
        "source\ttarget\nGata3\tEpcam\nFoxq1\tKrt5\nGata3\tCldn7\n", encoding="utf-8"
    )

    yaml.safe_dump(
        {"programs": {
            "contractile": {
                "label": "Contractile program",
                "description": "Synthetic contractile program for the fixture.",
                "color": "#A6CEE3", "genes": ["Acta2", "Cnn1", "Myh11", "Tagln"],
                "expected_direction": "increased", "contrast": "treatment_a_vs_control"},
            "epithelial_id": {
                "label": "Epithelial identity",
                "description": "Synthetic epithelial-identity program for the fixture.",
                "color": "#F4A6A6", "genes": ["Epcam", "Krt5", "Krt18", "Cldn7"],
                "expected_direction": "decreased", "contrast": "treatment_a_vs_control"},
        }, "program_order": ["contractile", "epithelial_id"]},
        (destination / "panels.yaml").open("w"), sort_keys=False,
    )

    yaml.safe_dump(
        {"study": "synthetic_publication_fixture", "hypotheses": [
            {"id": "h1", "statement": "The contractile program increases under treatment A.",
             "contrast": "treatment_a_vs_control", "expected_direction": "increased",
             "gene_panels": ["contractile"]},
            {"id": "h2", "statement": "The epithelial-identity program decreases under treatment A.",
             "contrast": "treatment_a_vs_control", "expected_direction": "decreased",
             "gene_panels": ["epithelial_id"]},
        ]},
        (destination / "hypotheses.yaml").open("w"), sort_keys=False,
    )

    yaml.safe_dump(
        {"figure_sets": {"main": {
            "title": "Synthetic publication figure",
            "description": "Exercises the publication DAG across every publication module.",
            "width": 180, "height": 240, "units": "mm", "columns": 2,
            "panels": [
                {"id": "A", "constructor": "pca"},
                {"id": "B", "constructor": "volcano", "contrast": "treatment_a_vs_control"},
                {"id": "C", "constructor": "ora_bidirectional", "contrast": "treatment_a_vs_control"},
                {"id": "D", "constructor": "cell_state_effects", "contrast": "treatment_a_vs_control"},
                {"id": "E", "constructor": "regulator_activity", "variant": "default",
                 "contrast": "treatment_a_vs_control"},
                {"id": "F", "constructor": "regulator_activity", "variant": "binding",
                 "contrast": "treatment_a_vs_control"},
                {"id": "G", "constructor": "string_enrichment", "contrast": "treatment_a_vs_control"},
                {"id": "H", "constructor": "program_heatmap_effects", "contrast": "treatment_a_vs_control"},
            ],
        }}},
        (destination / "recipe.yaml").open("w"), sort_keys=False,
    )

    config_path = destination / "project.yaml"
    data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    data["analysis"]["profile"] = "publication"
    data["species"]["taxonomy_id"] = 10090
    data["resources"]["providers"] = {"string": True, "dorothea": True}
    data["resources"]["cell_state_signatures"] = "signatures.yaml"
    data["resources"]["binding_prior_edges"] = "binding.tsv"
    data["hypotheses"] = {"panels": "panels.yaml", "claims": "hypotheses.yaml"}
    data["publication"] = {"recipe": "recipe.yaml"}
    config_path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return config_path


def test_standard_profile_dag_resolves(tmp_path):
    """The stock minimal (standard-profile) project plans a clean core DAG and
    pulls in none of the publication-only rules -- proving profile gating keeps
    the extra stages off unless the profile asks for them."""
    config_path = tmp_path / "minimal" / "project.yaml"
    shutil.copytree(TEMPLATE, tmp_path / "minimal")
    proc = _dry_run(config_path)
    assert proc.returncode == 0, proc.stderr or proc.stdout
    planned = _planned_rules(proc.stdout)
    assert {"study_qc", "contrast_de", "contrast_pathways", "report_html"} <= planned, sorted(planned)
    assert not (PUBLICATION_RULES & planned), sorted(PUBLICATION_RULES & planned)


def test_publication_profile_dag_resolves(tmp_path):
    """The full publication DAG assembles offline: every publication module's
    rule is present and reachable, and the figure assembler resolves all recipe
    panels (including the dual-view/binding regulator panels). A non-zero exit
    would mean a mis-wired rule -- a missing input, an ambiguity, or a cycle."""
    config_path = _publication_project(tmp_path)
    proc = _dry_run(config_path)
    assert proc.returncode == 0, proc.stderr or proc.stdout
    planned = _planned_rules(proc.stdout)
    missing = PUBLICATION_RULES - planned
    assert not missing, f"publication rules absent from DAG: {sorted(missing)}"


def test_displayed_data_are_declared_outputs(tmp_path):
    """Every displayed-data table a *built* panel exposes must be a declared
    engine output, not an untracked side-effect.

    ``figures.py`` advertises, per panel variant, the ``displayed_data`` tables
    that back its figure (the Source Data a reviewer audits). A rule's script may
    *write* such a table, but unless the table is also named in that rule's
    ``output:`` list, Snakemake does not own it: it is an untracked side-effect,
    and deleting it while its declared siblings stay current would NOT reschedule
    the producing rule -- so the table silently never comes back. That
    reschedule-on-delete guarantee is exactly Snakemake's declared-output
    property, so asserting "the table is a planned output" is the offline,
    R-free, network-free way to assert "deleting it regenerates it" (a true
    delete-and-rerun test would need a built result tree).

    The contract is gated on the figure actually being built: for each variant
    whose source figure is a planned output of this publication run, every one of
    its ``displayed_data`` tables must also be planned. Variants whose figure this
    run does not build (opt-in modules, on-demand DE/pathway diagnostics no recipe
    surfaces) are skipped -- their tables are legitimately absent here. This is
    the regression guard for the ``program_effects_forest_displayed`` and
    ``string_{up,down}_network_displayed`` gaps, and a general guard for any
    future panel that writes an audit table without declaring it.
    """
    config_path = _publication_project(tmp_path)
    proc = _dry_run(config_path)
    assert proc.returncode == 0, proc.stderr or proc.stdout
    plan = proc.stdout

    def _figure_built(source_stem: str) -> bool:
        tail = _result_tail(source_stem)
        return (tail + ".pdf") in plan or (tail + ".png") in plan

    checked = 0
    undeclared: list[str] = []
    for constructor in PANEL_REGISTRY.values():
        for variant_name, variant in constructor.variants.items():
            if not variant.displayed_data or not _figure_built(variant.source):
                continue
            for rel in variant.displayed_data:
                checked += 1
                if _result_tail(rel) not in plan:
                    undeclared.append(f"{constructor.id}/{variant_name}: {rel}")

    # The fixture builds enough panels that the guard exercises real tables; a
    # zero count would mean the gate silently matched nothing (e.g. output-format
    # drift) rather than that the contract holds.
    assert checked, "no built-panel displayed tables were checked -- gate is inert"
    assert not undeclared, (
        "displayed-data tables exposed by a built panel but NOT declared as engine "
        "outputs (untracked side-effects -- deleting them would not reschedule "
        f"their rule):\n  " + "\n  ".join(sorted(undeclared))
    )
