# Parity & Best-Practice Audit — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a ranked, evidence-backed gap register + audit report comparing `bulk-rna-frame` against the `lymphatic-flow-homeostasis` reference pipeline (13 stages + 16 gold figure panels + config/reproducibility infra), applying a 3-way best-practice judgment.

**Architecture:** A small typed gap-register schema + validator is built first (TDD). A numeric-evidence step runs the existing bridge and captures per-stage diffs against the reference `results/`. Then N independent audit clusters (dispatched as parallel subagents) each read reference↔frame for their slice, emit a register **fragment** file (isolated, no write conflicts), and apply trivial inline fixes. A final synthesis step validates, dedupes, ranks, and merges fragments into `docs/parity-gaps.yaml` + `docs/parity-audit.md`.

**Tech Stack:** Python 3.13, PyYAML, pytest; R (reference stages); Snakemake workflow; the existing `analysis/bulk_rna_frame/` bridge in the lymphatic repo.

## Global Constraints

- Reference repo path: `/home/ajbook/projects/lymphatic-flow-homeostasis` (abbrev `REF`).
- Frame repo path: `/home/ajbook/projects/bulk-rna-frame` (abbrev `FRAME`, cwd).
- Gold figure PDFs: `/tmp/cape_ref/*.pdf` (Set 1, 5 panels) and `/tmp/cape_ref/set2/*.pdf` (Set 2, 11 panels). If `/tmp` was cleared, re-extract from `/mnt/c/Users/ajboo/Downloads/Re_ Bulk RNA-Seq analysis .zip` and `... 2.zip`.
- Numeric ground truth: `REF/results/{cape_thoracic_duct,ligation_vs_sham,old_vs_young,rela_ko_vs_wt,bulk_rna_frame,bulk_rna_frame_publication,bulk_rna_frame_v2}`.
- Verdicts (exactly one per finding): `REGRESSION` | `FIDELITY-GAP` | `BEST-PRACTICE` | `OK`.
- Severity (exactly one): `S1` stats-correctness | `S2` numeric-parity | `S3` figure-fidelity | `S4` generality | `S5` polish.
- Target routing (exactly one): `harness` (#2) | `figures` (#3) | `stats` (#4) | `inline` (fixed now).
- Evidence mode (exactly one): `run` (backed by an executed diff) | `code-only` (from reading).
- All fragment/register files are YAML and MUST validate against the Task 1 schema.
- Inline fixes: only trivial/obvious; each committed separately and referenced by `commit` sha in its register entry with `status: fixed-inline`. Non-trivial gaps stay `status: open` and route to #2/#3/#4.
- Contrast semantics are sacred: positive effect = `numerator − denominator`, never inferred. Any inference is at least `S1`.
- Every inline fix keeps `pytest` green.

---

### Task 1: Gap-register schema + validator

**Files:**
- Create: `tools/parity/register_schema.py`
- Create: `tools/parity/validate_register.py`
- Test: `tests/test_parity_register.py`
- Create (fixture): `tests/fixtures/parity/valid_fragment.yaml`, `tests/fixtures/parity/invalid_fragment.yaml`

**Interfaces:**
- Produces: `validate_fragment(data: list[dict]) -> list[str]` returning a list of human-readable error strings (empty = valid); `REQUIRED_KEYS`, `VERDICTS`, `SEVERITIES`, `TARGETS`, `EVIDENCE_MODES` constants; `load_and_validate(path: str) -> list[dict]` that raises `ValueError` on invalid.
- Consumes: nothing (first task).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_parity_register.py
from pathlib import Path
import pytest
from tools.parity.register_schema import validate_fragment, load_and_validate

FIX = Path(__file__).parent / "fixtures" / "parity"

def _entry(**over):
    base = dict(
        id="A1-001", track="A", area="de",
        verdict="REGRESSION", severity="S1",
        summary="apeglm shrinkage not applied",
        detail="frame uses normal shrinkage; reference uses apeglm",
        reference=[{"file": "workflow/stages/de/stage.R", "line": 42}],
        frame=[{"file": "src/bulk_rna_frame/workflow/scripts/de.R", "line": 88}],
        evidence="run", target="stats", status="open", commit=None,
    )
    base.update(over)
    return base

def test_valid_entry_has_no_errors():
    assert validate_fragment([_entry()]) == []

def test_missing_key_reported():
    bad = _entry(); del bad["severity"]
    errs = validate_fragment([bad])
    assert any("severity" in e for e in errs)

def test_bad_enum_reported():
    errs = validate_fragment([_entry(verdict="MAYBE")])
    assert any("verdict" in e and "MAYBE" in e for e in errs)

def test_fixed_inline_requires_commit():
    errs = validate_fragment([_entry(status="fixed-inline", commit=None)])
    assert any("commit" in e for e in errs)

def test_duplicate_ids_reported():
    errs = validate_fragment([_entry(id="A1-001"), _entry(id="A1-001")])
    assert any("duplicate" in e.lower() for e in errs)

def test_load_and_validate_raises_on_invalid_file():
    with pytest.raises(ValueError):
        load_and_validate(str(FIX / "invalid_fragment.yaml"))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_parity_register.py -v`
Expected: FAIL — `ModuleNotFoundError: tools.parity.register_schema`.

- [ ] **Step 3: Write minimal implementation**

```python
# tools/parity/register_schema.py
from __future__ import annotations
import yaml

REQUIRED_KEYS = {"id","track","area","verdict","severity","summary",
                 "detail","reference","frame","evidence","target","status","commit"}
VERDICTS = {"REGRESSION","FIDELITY-GAP","BEST-PRACTICE","OK"}
SEVERITIES = {"S1","S2","S3","S4","S5"}
TARGETS = {"harness","figures","stats","inline"}
EVIDENCE_MODES = {"run","code-only"}
STATUSES = {"open","fixed-inline"}
TRACKS = {"A","B","C"}

def _check_ref_list(name, value, out, idx):
    if not isinstance(value, list) or not value:
        out.append(f"entry {idx}: {name} must be a non-empty list")
        return
    for j, ref in enumerate(value):
        if not isinstance(ref, dict) or "file" not in ref or "line" not in ref:
            out.append(f"entry {idx}: {name}[{j}] needs file and line")

def validate_fragment(data):
    errors = []
    seen = set()
    if not isinstance(data, list):
        return ["fragment must be a list of entries"]
    for i, e in enumerate(data):
        if not isinstance(e, dict):
            errors.append(f"entry {i}: not a mapping"); continue
        missing = REQUIRED_KEYS - set(e)
        for m in sorted(missing):
            errors.append(f"entry {i}: missing key '{m}'")
        if missing:
            continue
        if e["track"] not in TRACKS: errors.append(f"entry {i}: bad track {e['track']!r}")
        if e["verdict"] not in VERDICTS: errors.append(f"entry {i}: bad verdict {e['verdict']!r}")
        if e["severity"] not in SEVERITIES: errors.append(f"entry {i}: bad severity {e['severity']!r}")
        if e["target"] not in TARGETS: errors.append(f"entry {i}: bad target {e['target']!r}")
        if e["evidence"] not in EVIDENCE_MODES: errors.append(f"entry {i}: bad evidence {e['evidence']!r}")
        if e["status"] not in STATUSES: errors.append(f"entry {i}: bad status {e['status']!r}")
        if e["status"] == "fixed-inline" and not e.get("commit"):
            errors.append(f"entry {i}: status fixed-inline requires a commit sha")
        _check_ref_list("reference", e.get("reference"), errors, i)
        _check_ref_list("frame", e.get("frame"), errors, i)
        if e["id"] in seen:
            errors.append(f"entry {i}: duplicate id {e['id']!r}")
        seen.add(e["id"])
    return errors

def load_and_validate(path):
    data = yaml.safe_load(open(path, encoding="utf-8")) or []
    errs = validate_fragment(data)
    if errs:
        raise ValueError("invalid parity register:\n" + "\n".join(errs))
    return data
```

Create `tools/parity/__init__.py` (empty) and `tools/__init__.py` (empty) so imports resolve. Create the two fixtures: `valid_fragment.yaml` with one well-formed entry; `invalid_fragment.yaml` with a missing `severity`.

```python
# tools/parity/validate_register.py
import sys
from tools.parity.register_schema import load_and_validate
if __name__ == "__main__":
    for p in sys.argv[1:]:
        load_and_validate(p); print(f"OK: {p}")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_parity_register.py -v`
Expected: PASS (6 tests).

- [ ] **Step 5: Commit**

```bash
git add tools/parity tools/__init__.py tests/test_parity_register.py tests/fixtures/parity
git commit -m "feat: parity gap-register schema and validator"
```

---

### Task 2: Numeric-evidence capture (bridge run)

**Files:**
- Create: `docs/parity/evidence/README.md` (records how evidence was produced)
- Create: `docs/parity/evidence/<cohort>/*.txt` (captured diff outputs)
- Create: `tools/parity/capture_evidence.sh`

**Interfaces:**
- Produces: for each available cohort, a text artifact summarizing counts/DE diffs between `FRAME` output and `REF/results/<cohort>`, plus an environment stamp. Track A clusters cite these files with `evidence: run`.
- Consumes: `REF/analysis/bulk_rna_frame/{runner.py,compare_counts.py,compare_de.py}`, `REF/results/`.

- [ ] **Step 1: Inventory what the bridge can compare.** Run and read:

```bash
sed -n '1,120p' /home/ajbook/projects/lymphatic-flow-homeostasis/analysis/bulk_rna_frame/compare_counts.py
sed -n '1,160p' /home/ajbook/projects/lymphatic-flow-homeostasis/analysis/bulk_rna_frame/compare_de.py
ls /home/ajbook/projects/lymphatic-flow-homeostasis/results/cape_thoracic_duct
```

Determine which cohorts have both a reference result tree and a runnable frame config. Record findings in `docs/parity/evidence/README.md`.

- [ ] **Step 2: Attempt a frame run + diff for `cape_thoracic_duct`.** Write `tools/parity/capture_evidence.sh` that, for a cohort, invokes the bridge (`python analysis/bulk_rna_frame/runner.py ...` from the REF repo, or the frame CLI directly on the cohort config) and pipes `compare_counts.py` / `compare_de.py` output into `docs/parity/evidence/<cohort>/`. Capture stdout+stderr; do not fail the script on diff mismatch (mismatches are the evidence).

- [ ] **Step 3: Run it for every cohort that can run.**

```bash
bash tools/parity/capture_evidence.sh cape_thoracic_duct
bash tools/parity/capture_evidence.sh ligation_vs_sham
bash tools/parity/capture_evidence.sh old_vs_young
bash tools/parity/capture_evidence.sh rela_ko_vs_wt
```

- [ ] **Step 4: Record the fallback.** If a cohort cannot run in this environment (missing data root, missing conda env, container unavailable), write that fact into `docs/parity/evidence/README.md` and note that Track A findings for the affected stages will be `evidence: code-only`. **An environment that cannot reproduce a published cohort is itself an `S2` finding** — add it to a fragment `docs/parity/fragments/evidence.yaml` (validate with Task 1).

- [ ] **Step 5: Commit**

```bash
git add tools/parity/capture_evidence.sh docs/parity/evidence docs/parity/fragments/evidence.yaml
git commit -m "chore: capture numeric parity evidence from bridge runs"
```

---

## Audit-cluster protocol (applies to Tasks 3–11)

Each cluster task is executed by **one subagent**. Tasks 3–11 are mutually
independent and dispatched in parallel. Every cluster subagent MUST:

1. Read its listed **reference** files (`REF/...`) and **frame** files (`FRAME/...`) in full.
2. For Track A: read the matching `docs/parity/evidence/<cohort>/` artifacts; cite them with `evidence: run` where a diff exists, else `evidence: code-only`.
3. For Track B: open the matching gold PDF(s) (Read tool renders PDFs) as the visual ground truth.
4. Apply the **3-way judgment** (reference vs frame vs best practice) and emit one register entry per meaningful comparison, using `id` prefix = cluster id (e.g. `A1-001`). Give `OK` entries too, so coverage is provable.
5. Apply trivial/obvious inline fixes only (typos, a wrong default that is unambiguous, a missing `set.seed`, a mislabeled axis); commit each and set `status: fixed-inline` + `commit`. Anything requiring judgment stays `open` and routes to `harness`/`figures`/`stats`.
6. Write results to `docs/parity/fragments/<cluster-id>.yaml` and validate:
   `python -m tools.parity.validate_register docs/parity/fragments/<cluster-id>.yaml`
7. Ensure `python -m pytest -q` is green if any code changed.
8. Commit the fragment.

Fragment `id` prefixes: A1,A2,A3,A4,A5 (Track A), B1,B2,B3 (Track B), C1 (Track C).

---

### Task 3: Cluster A1 — Counts, QC, DE

**Files (compare):**
- Reference: `REF/workflow/stages/core/{count_genes.py,deseq2.R,qc_gate.py,strand_test.py,preflight_bams.py,limma_sensitivity.R,validate_study.py}`, `REF/workflow/stages/qc/stage.R`, `REF/workflow/stages/de/{stage.R,de_annotation.R,de_io.R,de_labels.R,de_tables.R}`, `REF/workflow/lib/{contrast.R,annotate.R}`
- Frame: `FRAME/src/bulk_rna_frame/workflow/scripts/{de.R,qc.R,materialize_inputs.py,front_door.py}`
- Evidence: `docs/parity/evidence/*/compare_counts*, */compare_de*`
- Output: `docs/parity/fragments/A1.yaml`

- [ ] **Step 1:** Follow the audit-cluster protocol. Focus checks: featureCounts/count parameters (strandedness, multimappers), DESeq2 design/fitType/`betaPrior`, apeglm vs normal shrinkage, independent filtering, contrast direction (`numerator − denominator`), p-adjust method, gene annotation join. Cite the `compare_de` evidence for numeric verdicts.
- [ ] **Step 2:** Emit `docs/parity/fragments/A1.yaml`; validate; commit.

---

### Task 4: Cluster A2 — Pathways & Ontology

**Files (compare):**
- Reference: `REF/workflow/stages/pathways/{stage.R,pathway_fgsea.R,pathway_ora.R,pathway_gsva.R,pathway_gene_sets.R,pathway_rankings.R,pathway_tables.R}`, `REF/workflow/stages/ontology/{stage.R,ontology_go.R,ontology_kegg.R,ontology_string.R,ontology_tables.R}`
- Frame: `FRAME/src/bulk_rna_frame/workflow/scripts/{pathways.R,ontology.R}`
- Output: `docs/parity/fragments/A2.yaml`

- [ ] **Step 1:** Follow the protocol. Focus: fgsea ranking metric + ties + `eps`/`nPermSimple`, seed determinism, ORA universe definition + direction split, GSVA/ssGSEA method + kcdf, MSigDB/GO/KEGG collection versions and pinning, multiple-testing within vs across collections.
- [ ] **Step 2:** Emit `A2.yaml`; validate; commit.

---

### Task 5: Cluster A3 — Regulators, Networks, GRN

**Files (compare):**
- Reference: `REF/workflow/stages/regulators/{stage.R,viper.R,viper_dorothea.R,grn_dorothea.R,regulators_prior.R,regulators_tables.R}`, `REF/workflow/stages/network/{stage.R,network_integrate.R,network_tables.R}`
- Frame: `FRAME/src/bulk_rna_frame/workflow/scripts/{regulators.R,networks.py,grn.py}`
- Output: `docs/parity/fragments/A3.yaml`

- [ ] **Step 1:** Follow the protocol. Focus: VIPER method (`aREA`), regulon source (GTRD/DoRothEA confidence levels A–E), unsigned vs signed target sets, minsize, NES vs activity scaling, STRING score threshold + version + induced-network batching (see recent frame commit `4546a17`), DoRothEA edge sign handling, GRN edge-type distinctions (regulatory vs co-expression vs association — must stay explicit).
- [ ] **Step 2:** Emit `A3.yaml`; validate; commit.

---

### Task 6: Cluster A4 — Composition & Hypotheses

**Files (compare):**
- Reference: `REF/workflow/stages/composition/{stage.R,deconvolution_signatures.R,deconvolution_tables.R,deconvolution_io.R}`, `REF/workflow/stages/hypothesis/{stage.R,hypothesis_gene_panels.R,hypothesis_tables.R}`, `REF/workflow/lib/evidence.R`, `REF/config/deconvolution_signatures.yaml`
- Frame: `FRAME/src/bulk_rna_frame/workflow/scripts/{composition.R,hypotheses.py}`, `FRAME/src/bulk_rna_frame/schemas/{hypotheses.schema.yaml,hypothesis_panels.schema.yaml,signatures.schema.yaml}`
- Output: `docs/parity/fragments/A4.yaml`

- [ ] **Step 1:** Follow the protocol. Focus: cell-state = **relative signature scores, not cell fractions** (README claim — verify the frame never implies fractions), scoring method (GSVA/AUCell/mean-z), signature provenance, hypothesis evidence aggregation logic, direction handling in panels (see frame commit `6e5dac1` re: panels without color annotations).
- [ ] **Step 2:** Emit `A4.yaml`; validate; commit.

---

### Task 7: Cluster A5 — Advanced (frame-only) & shared libs

**Files (compare):**
- Reference (partial/none — flag as best-practice baseline): `REF/workflow/stages/core/limma_sensitivity.R`, `REF/workflow/envs/{sva_mediation.yaml,wgcna.yaml}`, `REF/workflow/lib/{analysis_set.py,configmerge.py,io.R,packages.R}`, `REF/workflow/stages/common/{adapter.py,runtime_utils.R}`
- Frame: `FRAME/src/bulk_rna_frame/workflow/scripts/{sva.R,wgcna.R,mediation.R,power.py,multilayer.py,utils.R,resources.R}`
- Output: `docs/parity/fragments/A5.yaml`

- [ ] **Step 1:** Follow the protocol. These advanced stages have **no reference stage script** — judge on best-practice alone (verdict mostly `BEST-PRACTICE` or `OK`, severity ≤ S4 unless a correctness bug). Focus: SVA n.sv estimation and whether surrogate variables feed DE correctly, WGCNA soft-power selection + determinism, mediation assumptions stated, power-analysis inputs, multilayer integration method. Flag small-sample guards (README says advanced analyses warn, not fail — verify).
- [ ] **Step 2:** Emit `A5.yaml`; validate; commit.

---

### Task 8: Cluster B1 — Figure Set 1 (STRING/Regulator panels A–E)

**Files (compare):**
- Reference: `REF/analysis/cape_xizhao_edits/{figure2_string_enrichment_panel_A.R,figure2_string_up_network_panel_B.R,figure2_string_down_network_panel_C.R,figure2_regulator_activity_panel_D.R,figure2_dorothea_grn_radial_option_E.R}` + the relevant sections of `REF/analysis/cape_xizhao_edits/publication_figure_functions.R`
- Frame: `FRAME/src/bulk_rna_frame/{figures.py,workflow/scripts/publication.R,workflow/scripts/assemble.py}`, `FRAME/src/bulk_rna_frame/schemas/figure_recipe.schema.yaml`
- Gold PDFs: `/tmp/cape_ref/Panel A - STRING Functional Enrichment.pdf`, `Panel B - Upregulated STRING Network.pdf`, `Panel C - Downregulated STRING Network.pdf`, `Panel D - Regulator Activity.pdf`, `Panel E - DoRothEA Radial GRN.pdf`
- Output: `docs/parity/fragments/B1.yaml`

- [ ] **Step 1:** Follow the protocol. For each of the 5 panels compare: data source, faceting (Down/Leading-edge/Up columns in A), shape=STRING category, size=mapped genes, color=−log10(FDR) scale, network node/edge encodings (B/C), heatmap+lollipop composition (D), radial layout + sector programs + edge sign colors (E). Record any panel with **no frame constructor** as `FIDELITY-GAP`, `S3`, target `figures`.
- [ ] **Step 2:** Emit `B1.yaml`; validate; commit.

---

### Task 9: Cluster B2 — Figure Set 2 core/DE panels (A–D)

**Files (compare):**
- Reference: `REF/analysis/cape_xizhao_edits/{pca_requested_edit.R,correlation_requested_edit.R,cell_state_requested_edit.R,de_heatmap_requested_edit.R,program_panels_requested_edit.R}` + `publication_figure_functions.R` sections
- Frame: `FRAME/src/bulk_rna_frame/{figures.py,workflow/scripts/publication.R}`
- Gold PDFs: `/tmp/cape_ref/set2/{Panel A - PCA and Sample Correlation.pdf,Panel B - Cell-State Signatures.pdf,Panel C - Original Volcano Plot.pdf,Panel D - Option 1 - Global Clustering.pdf,Panel D - Option 2 - Program Grouped.pdf,Panel D - Option 3 - Direct Program Labels.pdf}`
- Output: `docs/parity/fragments/B2.yaml`

- [ ] **Step 1:** Follow the protocol. Compare PCA variance-explained + correlation matrix, cell-state signature plot (relative scores), volcano thresholds/labels, and the **three DE-clustering options** (global / program-grouped / direct-labels) — note which options the frame supports.
- [ ] **Step 2:** Emit `B2.yaml`; validate; commit.

---

### Task 10: Cluster B3 — Figure Set 2 pathway panels (E–I)

**Files (compare):**
- Reference: `REF/analysis/cape_xizhao_edits/{ora_combined_requested_edit.R,gsva_heatmap_requested_edit.R,gsea_advanced_requested_edit.R}` + `de_heatmap_requested_edit.R` (integrated H) + `cell_state`/`program_panels` (violins I) + `publication_figure_functions.R`
- Frame: `FRAME/src/bulk_rna_frame/{figures.py,workflow/scripts/publication.R}`
- Gold PDFs: `/tmp/cape_ref/set2/{Panel E - GO Biological Process ORA.pdf,Panel F - Hallmark GSVA Heatmap.pdf,Panel G - Curated Program GSEA.pdf,Panel H - Integrated Heatmap and Effects.pdf,Panel I - Consolidated Violins.pdf}`
- Output: `docs/parity/fragments/B3.yaml`

- [ ] **Step 1:** Follow the protocol. Compare GO BP ORA dotplot (direction split, term selection), Hallmark GSVA heatmap (scaling, clustering), curated multi-track GSEA curves, integrated heatmap+effects composition, consolidated violins (grouping, stats annotation).
- [ ] **Step 2:** Emit `B3.yaml`; validate; commit.

---

### Task 11: Cluster C1 — Config generality & reproducibility infra

**Files (compare):**
- Reference: `REF/config/{config.schema.yaml,defaults.yaml,deconvolution_signatures.yaml}`, `REF/workflow/lib/configmerge.py`, `REF/results/` (7 cohorts as generality evidence)
- Frame: `FRAME/src/bulk_rna_frame/{config.py,cli.py,collection.py,verification.py}`, `FRAME/src/bulk_rna_frame/schemas/{project.schema.yaml,collection.schema.yaml}`, `FRAME/src/bulk_rna_frame/workflow/rules/{providers.smk,report.smk}`, `FRAME/src/bulk_rna_frame/workflow/scripts/manifest.py`, `FRAME/{Dockerfile,Apptainer.def,environment.yaml}`, `FRAME/src/bulk_rna_frame/workflow/envs/*.yaml`
- Output: `docs/parity/fragments/C1.yaml`

- [ ] **Step 1:** Follow the protocol. Focus: does `project.yaml`/config flex to a new cohort (not CAPE-hardcoded)? Are the 7 reference cohorts reproducible via frame config? Provider caching + receipts, `verify` tolerances, release manifest completeness (seeds, checksums, env, revision), container/env pinning, seed propagation across R and Python stages. Generality gaps route to `harness`; reproducibility gaps `S2`/`S4`.
- [ ] **Step 2:** Emit `C1.yaml`; validate; commit.

---

### Task 12: Synthesis — ranked register + audit report

**Files:**
- Create: `docs/parity-gaps.yaml` (merged, ranked register)
- Create: `docs/parity-audit.md` (human-readable report)
- Create: `tools/parity/merge_fragments.py`
- Test: `tests/test_parity_merge.py`

**Interfaces:**
- Consumes: all `docs/parity/fragments/*.yaml`.
- Produces: `merge_fragments(fragment_dir: str) -> list[dict]` sorted by (severity S1→S5, then track), with globally-unique ids; raises `ValueError` if any fragment is invalid or ids collide across fragments.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_parity_merge.py
from tools.parity.merge_fragments import merge_fragments, coverage_report

def test_merge_sorts_by_severity(tmp_path):
    (tmp_path / "A1.yaml").write_text(
        "- {id: A1-1, track: A, area: de, verdict: OK, severity: S3, "
        "summary: s, detail: d, reference: [{file: r, line: 1}], "
        "frame: [{file: f, line: 1}], evidence: run, target: figures, "
        "status: open, commit: null}\n")
    (tmp_path / "A2.yaml").write_text(
        "- {id: A2-1, track: A, area: pathways, verdict: REGRESSION, severity: S1, "
        "summary: s, detail: d, reference: [{file: r, line: 1}], "
        "frame: [{file: f, line: 1}], evidence: run, target: stats, "
        "status: open, commit: null}\n")
    merged = merge_fragments(str(tmp_path))
    assert [m["id"] for m in merged] == ["A2-1", "A1-1"]  # S1 before S3

def test_coverage_report_flags_missing(tmp_path):
    (tmp_path / "A1.yaml").write_text(
        "- {id: A1-1, track: A, area: de, verdict: OK, severity: S3, "
        "summary: s, detail: d, reference: [{file: r, line: 1}], "
        "frame: [{file: f, line: 1}], evidence: run, target: figures, "
        "status: open, commit: null}\n")
    cov = coverage_report(str(tmp_path))
    assert "de" in cov["stages_seen"]
    assert cov["panels_missing"]  # 16 panels not all covered
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_parity_merge.py -v`
Expected: FAIL — module missing.

- [ ] **Step 3: Implement `tools/parity/merge_fragments.py`.** Load every fragment via `load_and_validate` (Task 1), assert global id-uniqueness, sort by `(severity_rank, track, id)` where `severity_rank = {"S1":0,...,"S5":4}`. Add `coverage_report(dir)` returning `{"stages_seen": set(...), "panels_seen": set(...), "stages_missing": [...], "panels_missing": [...]}` checking the 13 expected stages (`qc,de,pathways,ontology,composition,regulators,network,hypothesis,sva,wgcna,mediation,power,multilayer`) and 16 expected panels (Set1 A–E, Set2 A,B,C,D1,D2,D3,E,F,G,H,I).

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_parity_merge.py -v`
Expected: PASS.

- [ ] **Step 5: Generate the register.** Run the merger over `docs/parity/fragments/`, write sorted result to `docs/parity-gaps.yaml`. Fail loudly (do not write) if `coverage_report` shows any of the 13 stages or 16 panels uncovered — that means an audit cluster missed scope; go back and fill it.

- [ ] **Step 6: Write `docs/parity-audit.md`.** Human-readable report generated/curated from the register: executive summary (counts by verdict×severity), a per-track section, the top S1/S2 findings with two-sided `file:line` evidence, the list of inline fixes applied (with commit shas), and an explicit "routed to sub-project #2/#3/#4" backlog. End with the acceptance-criteria checklist from the design spec, each checked.

- [ ] **Step 7: Commit**

```bash
git add tools/parity/merge_fragments.py tests/test_parity_merge.py docs/parity-gaps.yaml docs/parity-audit.md
git commit -m "feat: synthesize ranked parity gap register and audit report"
```

---

## Self-Review (completed during authoring)

- **Spec coverage:** Track A → Tasks 3–7 (13 stages: A1 qc/de, A2 pathways/ontology, A3 regulators/network/grn, A4 composition/hypothesis, A5 sva/wgcna/mediation/power/multilayer). Track B → Tasks 8–10 (all 16 panels). Track C → Task 11. 3-way judgment + verdicts/severities → Global Constraints + protocol. Evidence standard → Task 2 + protocol steps 2–3. Parallel execution → Tasks 3–11 independent. Deliverables (`parity-audit.md`, `parity-gaps.yaml`) → Task 12. Acceptance criteria → Task 12 Steps 5–6. Inline fixes → protocol step 5 + schema `fixed-inline`/`commit` rule (Task 1). ✔
- **Placeholder scan:** no TBD/TODO; all code steps carry real code; audit-cluster tasks carry exact file lists and focus checklists. ✔
- **Type consistency:** `validate_fragment`/`load_and_validate`/`REQUIRED_KEYS` defined in Task 1 and reused verbatim in Tasks 2–12; `merge_fragments`/`coverage_report` defined and tested in Task 12. Verdict/severity/target/evidence enums identical across schema, constraints, and cluster tasks. ✔
