# bulk-rna-frame Parity & Best-Practice Audit Report

**Date:** 2026-08-12  
**Auditor:** Parallel agent clusters (A1-A5, B1-B3, C1) + synthesis  
**Reference:** lymphatic-flow-homeostasis pipeline (2,937-line figure engine)  
**Frame:** bulk-rna-frame v0.x (generalized productization)  

---

## Executive Summary

This audit compares every stage (13) and figure panel (16) in `bulk-rna-frame` against the mature lymphatic reference implementation to assess statistical correctness, numeric parity, figure fidelity, and cohort generality. A total of **139 findings** across **3 tracks** were identified and ranked by severity (S1 correctness → S5 polish).

### Findings by Verdict × Severity

| Verdict         | S1 | S2 | S3 | S4 | S5 | Total |
|-----------------|----|----|----|----|----|----|
| **REGRESSION**  | 1 | 10 | 9 | 3 | 2 | 25 |
| **FIDELITY-GAP** | 0 | 2 | 25 | 10 | 2 | 39 |
| **BEST-PRACTICE** | 2 | 3 | 4 | 10 | 4 | 23 |
| **OK**          | 3 | 15 | 12 | 14 | 8 | 52 |
| **Total**       | **6** | **30** | **50** | **37** | **16** | **139** |

### Track Breakdown

- **Track A (Stages / Numeric fidelity):** 74 findings across 13 stages (qc, de, pathways, ontology, composition, regulators, network, hypothesis, sva, wgcna, mediation, power, multilayer)
- **Track B (Figure fidelity):** 53 findings across 16 gold panels (Set1 A-E, Set2 A-I)
- **Track C (Config / Infra):** 12 findings on generality, reproducibility, and tooling

---

## Track A: Stage / Numeric Fidelity

74 findings across 13 stages. The frame's core counting and DESeq2 layer is faithful (cape_thoracic_duct counts: 56748/56748 genes exact match; DE fold-changes and adjusted p-values numerically equivalent). Key issues:

### Highlights

- **Counts & QC:** featureCounts invocation byte-for-byte equivalent; FRAME's strand inference uses dominance ratio (correct) vs REF's biased argmax.
- **DE:** DESeq2 core (design, fit, shrinkage, BH correction) faithful; contrast direction explicit on both sides; FRAME survives failed parametric fits REF would crash.
- **Pathways/Ontology:** Multiple-testing scope regression (pooled BH across collections/directions instead of per collection×direction) — see systemic root cause below.
- **Regulators/Networks:** VIPER/DoRothEA/STRING methods faithful; FRAME's STRING score threshold differs (0.700 vs REF's 0.400 then clamp).
- **Composition/Hypothesis:** Cell-state signatures and hypothesis testing faithful; direction never inferred.
- **Advanced (SVA/WGCNA/mediation/power/multilayer):** Most stages are OK; SVA lacks determinism (no seed on permutation-based n.sv inference).

---

## Track B: Figure Fidelity

53 findings across 16 panels. The frame's figure architecture is modular (registry + per-panel constructors) but many panels lack directional faceting, advanced encodings, or gold-PDF-matching transforms.

### Highlights

**Set 1 (networks/regulators):**
- **A (STRING enrichment):** FRAME pools genes; REF facets Down / GSEA-LE / Up seed sets.
- **B/C (STRING networks):** FRAME lacks directional split; REF has separate up/down networks.
- **D (VIPER activity):** Shape encoding missing; FRAME omits gene-set size in labels.
- **E (DoRothEA radial):** Chord layout missing; FRAME uses force-directed.

**Set 2 (integrated analysis):**
- **A (PCA/correlation):** PCA biplot missing; correlation dendrogram differs.
- **B (cell-state):** Faithful.
- **C (volcano):** Significance line displays raw-p threshold on raw-p volcano (misleading FDR encoding); off-scale handling differs.
- **D (DE-clustering):** Clustering methods differ; program composition and scaling differ.
- **E-I (ORA/GSVA/GSEA/heatmap/violins):** Multiple-testing scope issues (E); GSVA z-score clamp differs (F); leading-edge labels missing (G); FDR encoding differs (H); bracket annotations use Wilcoxon recompute instead of modeled FDR (I).

---

## Track C: Config / Infra / Reproducibility

12 findings on cohort generality, reproducibility infra, and tooling. The frame has strong reproducibility foundations (manifest, seeds, provider caching, verification harness) but lacks config migration for legacy cohorts and REF's defaults inheritance mechanism.

### Highlights

- **Config migration:** 3 of 4 published cohorts cannot run without manual migration (evidence-001).
- **Reproducibility infra:** Manifest captures seeds/env/revisions; provider caching is content-addressable; verification harness has configurable tolerances.
- **Generality:** Hypothesis/figure schemas flex to new cohorts; network providers/caching support multiple sources.

---

## Critical S1/S2 Findings

The highest-severity findings requiring immediate attention:


### 1. A1-005 [S1 BEST-PRACTICE]

**FRAME survives a failed parametric dispersion fit that would crash REF, but its fallback uses un-shrunk gene-wise dispersions**

- **Reference:** workflow/stages/core/deseq2.R:29
- **Frame:** src/bulk_rna_frame/workflow/scripts/de.R:40

### 2. A1-007 [S1 OK]

**Contrast direction is derived explicitly from numerator/denominator on both sides; positive = up in numerator**

- **Reference:** workflow/stages/core/deseq2.R:26
- **Frame:** src/bulk_rna_frame/workflow/scripts/de.R:35

### 3. A3-006 [S1 OK]

**Regulator differential direction is read from the contrast, never inferred from group order**

- **Reference:** workflow/stages/regulators/regulators_prior.R:260
- **Frame:** src/bulk_rna_frame/workflow/scripts/regulators.R:101

### 4. A4-001 [S1 OK]

**Cell-state outputs are labeled relative signature/composition scores and never imply estimated cell fractions**

- **Reference:** workflow/stages/composition/deconvolution_signatures.R:198
- **Frame:** src/bulk_rna_frame/workflow/scripts/composition.R:85

### 5. A4-006 [S1 BEST-PRACTICE]

**REF hypothesis pathway/reference stats infer treatment-vs-control from sample-sheet appearance order, ignoring the contrast; FRAME consumes contrast-oriented effects**

- **Reference:** workflow/stages/hypothesis/hypothesis_gene_panels.R:376
- **Frame:** src/bulk_rna_frame/workflow/scripts/hypotheses.py:104

### 6. B3-017 [S1 REGRESSION]

**REF annotates brackets with the modeled DESeq2 DE FDR; FRAME recomputes a Wilcoxon test on VST expression that cannot reach the significance the panel shows at n=3/group**

- **Reference:** analysis/cape_xizhao_edits/program_panels_requested_edit.R:288
- **Frame:** src/bulk_rna_frame/workflow/scripts/publication.R:174

### 7. A1-001 [S2 OK]

**featureCounts command construction is byte-for-byte equivalent; counts exactly equal on cape**

- **Reference:** workflow/stages/core/count_genes.py:7
- **Frame:** src/bulk_rna_frame/workflow/scripts/materialize_inputs.py:170

### 8. A1-002 [S2 BEST-PRACTICE]

**REF strand auto-detection argmaxes assignment rate over modes {0,1,2} and is biased to unstranded; FRAME's dominance ratio is correct**

- **Reference:** workflow/stages/core/strand_test.py:48
- **Frame:** src/bulk_rna_frame/workflow/scripts/materialize_inputs.py:243

### 9. A1-004 [S2 OK]

**DESeq2 core (design, parametric fit, betaPrior=FALSE, independent filtering, BH) is equivalent; numerics within tol**

- **Reference:** workflow/stages/core/deseq2.R:29
- **Frame:** src/bulk_rna_frame/workflow/scripts/de.R:41

### 10. A1-008 [S2 OK]

**Significance classification uses raw padj plus shrunken |LFC| on both sides; NA-guarded**

- **Reference:** workflow/stages/de/de_tables.R:519
- **Frame:** src/bulk_rna_frame/workflow/scripts/de.R:87

### 11. A2-002 [S2 REGRESSION]

**FRAME passes tied ranks to fgsea with no tie-breaking and suppresses the warning; REF applies deterministic tie-breaking (an explicit audit fix)**

- **Reference:** workflow/stages/pathways/pathway_fgsea.R:135, workflow/stages/pathways/pathway_fgsea.R:424
- **Frame:** src/bulk_rna_frame/workflow/scripts/pathways.R:199, src/bulk_rna_frame/workflow/scripts/pathways.R:201

### 12. A2-005 [S2 REGRESSION]

**FRAME's ORA background is every gene with a symbol (including NA-padj / independent-filtered rows); REF restricts to genes with finite stat, padj and log2FC**

- **Reference:** workflow/stages/pathways/pathway_ora.R:173
- **Frame:** src/bulk_rna_frame/workflow/scripts/pathways.R:220

### 13. A2-009 [S2 REGRESSION]

**FRAME runs one fgseaMultilevel over all collections pooled (one BH family); REF runs fgsea per collection so BH is applied within each collection**

- **Reference:** workflow/stages/pathways/pathway_fgsea.R:602, workflow/stages/pathways/pathway_fgsea.R:438
- **Frame:** src/bulk_rna_frame/workflow/scripts/pathways.R:201, src/bulk_rna_frame/workflow/scripts/pathways.R:218

### 14. A2-010 [S2 REGRESSION]

**FRAME BH-adjusts ORA once across both directions and all collections pooled; REF adjusts within each collection x direction cell**

- **Reference:** workflow/stages/pathways/pathway_ora.R:646
- **Frame:** src/bulk_rna_frame/workflow/scripts/pathways.R:223, src/bulk_rna_frame/workflow/scripts/pathways.R:228

### 15. A2-011 [S2 REGRESSION]

**FRAME's ontology stage reuses the pooled ORA adjusted p-value; GO/KEGG FDR is not computed within the GO or KEGG family**

- **Reference:** workflow/stages/ontology/ontology_go.R:58, workflow/stages/ontology/ontology_kegg.R:56
- **Frame:** src/bulk_rna_frame/workflow/scripts/ontology.R:17, src/bulk_rna_frame/workflow/scripts/ontology.R:33


---

## Systemic Root Cause: Multiple-Testing Scope

**Affected findings:** A2-009 (fgsea), A2-010 (ORA), A2-011 (ontology FDR inheritance), B3-001 (Set2-E GO BP ORA FDR family)

All four findings share a single root cause: FRAME applies BH multiple-testing correction once across **pooled collections and both directions** (one large p-value family), while REF applies BH **separately within each collection×direction cell** (smaller, focused families). This inflates the effective family size in FRAME, reduces power, and violates the statistical best practice of controlling FDR within the natural structure of the data.

**Implication:** FRAME's pathway/ontology adjusted p-values are overly conservative and not directly comparable to REF's per-collection results.

**Fix route:** Sub-project #4 (statistical correctness audit + fixes) — refactor pathways/ontology stages to compute per-collection×direction BH adjustment.

---

## Inline Fixes

**Applied inline (status `fixed-inline`, commit `2e0c662`).** Four trivial/obvious defects were corrected during the serial inline pass; all tests remain green:

- **A5-002** (S2, BEST-PRACTICE): `sva.R` had no `set.seed` before `sva::num.sv(method="be")`, so the surrogate-variable count — and the whole SV-adjusted DE sensitivity result — was nondeterministic. Fixed by seeding from `cfg$analysis$random_seed` (matching `mediation.R`).
- **A5-014** (S4, REGRESSION): `utils.R::read_counts_contract` lacked REF's 2³¹ overflow guard, so counts ≥ 2³¹ would silently coerce to `NA`. Fixed by erroring on all-NA matrices and on values ≥ 2147483647 before integer coercion.
- **B2-012** (S3, REGRESSION): the volcano drew a horizontal `geom_hline` at `-log10(FDR)` on a raw-p axis, implying a raw-p cutoff that does not exist; REF deliberately omits it. Line removed.
- **B3-007** (S5, REGRESSION): GSVA row z-scores were clamped to ±2 vs the gold legend's ±1.5. Restored to ±1.5 (`row_zscore`'s own default).

**Carried `target: inline` but require no fix (verdict `OK`).** These are faithful/equivalent behaviors recorded for coverage, not gaps: A3-007, C1-005, C1-006, C1-007, C1-008, C1-009, C1-010.

**Re-routed off inline.** C1-004 (defaults-inheritance mechanism) is *not* a one-liner — it needs a schema field plus deep-merge logic — so it moves to the harness backlog rather than an inline edit.

---

## Routed Backlog by Sub-Project

Findings with `status: open` are routed to sub-projects #2–#4 for systematic remediation:

- **Figures:** 54 findings
- **Harness:** 15 findings
- **Stats:** 59 findings


**Next steps:**
- Sub-project #2 (Numerical-parity harness): Close 15 harness findings.
- Sub-project #3 (Figure engine port): Close 54 figure findings.
- Sub-project #4 (Statistical correctness): Close 59 stats findings + systemic multiple-testing fix.

---

## Acceptance Criteria

Per the design spec (`docs/superpowers/specs/2026-08-12-parity-audit-design.md`), the audit meets the following criteria:

- [x] **Every stage (13) and panel (16) has a verdict with two-sided file:line evidence.** All 13 stages (qc, de, pathways, ontology, composition, regulators, network, hypothesis, sva, wgcna, mediation, power, multilayer) and 16 panels (Set1 A-E, Set2 A-I) are covered with reference+frame citations.

- [x] **Numeric verdicts backed by real diff.** Findings marked `evidence: run` (vs `code-only`) are backed by actual re-runs through the diff bridge (`runner.py`, `compare_counts.py`, `compare_de.py`) on cape_thoracic_duct. Code-only findings are explicitly marked where the bridge lacks a comparator.

- [x] **Register ranked by severity and routed.** The machine-readable register (`docs/parity-gaps.yaml`) contains 139 entries sorted by (severity, track, id). Each entry specifies `target` (harness/figures/stats/inline) and `status` (open/fixed-inline). Inline fixes include commit references.

- [x] **Config/generality and reproducibility-infra findings captured (Track C).** 12 Track C findings assess cohort generality (hypothesis/figure schemas, defaults inheritance) and reproducibility infra (manifest, seeds, provider caching, verification harness, config migration).

---

## Appendix: Coverage Report

- **Stages seen:** 13 / 13 — composition, de, hypothesis, mediation, multilayer, network, ontology, pathways, power, qc, regulators, sva, wgcna
- **Panels seen:** 16 / 16 — Set1-A, Set1-B, Set1-C, Set1-D, Set1-E, Set2-A, Set2-B, Set2-C, Set2-D1, Set2-D2, Set2-D3, Set2-E, Set2-F, Set2-G, Set2-H, Set2-I
- **Stages missing:** 0 (gate: PASS)
- **Panels missing:** 0 (gate: PASS)

**Coverage gate:** ✓ PASS (all 13 stages and 16 panels covered)

---

**Register:** `docs/parity-gaps.yaml` (139 entries)  
**Evidence:** `analysis/bulk_rna_frame/` diff bridge + cape_thoracic_duct re-run  
**Next:** Inline fixes applied (commit `2e0c662`) → Sub-project #2 (harness) → Sub-project #3 (figures) → Sub-project #4 (stats)
