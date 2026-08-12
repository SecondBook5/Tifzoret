# Parity & Best-Practice Audit — Design Spec

**Date:** 2026-08-12
**Status:** Approved (design), pending implementation plan
**Sub-project:** 1 of 5 in the "impeccable bulk RNA-seq engine" program

## Program context

`bulk-rna-frame` is the productized generalization of the mature Snakemake
pipeline in `~/projects/lymphatic-flow-homeostasis`. The goal of the overall
program is a repeatable, cohort-configurable, hypothesis-driven bulk RNA-seq
engine (mouse-first, human later) that reliably produces publication-grade
figures matching the CAPE/lymphatic gold set — statistically correct,
reproducible, and visually equivalent.

The reference figure gold set (16 panels) lives in two archives and was produced
by `analysis/cape_xizhao_edits/` in the lymphatic repo:

- Set 1 (5 panels): STRING functional enrichment (A), up/down STRING networks
  (B/C), VIPER regulator activity (D), DoRothEA radial regulon map (E).
- Set 2 (11 panels): PCA + sample correlation (A), cell-state signatures (B),
  volcano (C), DE clustering options (D×3), GO BP ORA (E), Hallmark GSVA (F),
  curated program GSEA (G), integrated heatmap+effects (H), consolidated
  violins (I).

The whole engine `publication_figure_functions.R` reference is **2,937 lines**;
the frame currently reduces the figure layer to `figures.py` (~448) +
`publication.R` (~209). Each reference stage is decomposed into
`io/plots/tables/theme/stage.R` modules; the frame collapses each into a single
script. Both facts imply real divergence risk.

The program decomposes into five sub-projects. This spec covers **only #1**:

1. **Parity & best-practice audit** (this doc) — foundation.
2. Numerical-parity harness.
3. Port/generalize the figure engine.
4. Statistical-correctness audit + fixes.
5. Release hardening.

## Decisions locked with the user

- **Audit first.** Frame confirmed: `bulk-rna-frame` = productized
  generalization of the lymphatic Snakemake pipeline.
- **Best practice wins.** This is a 3-way judgment — reference behavior vs.
  frame behavior vs. statistical best practice. When the reference itself is
  questionable, the frame should do the correct thing, and the audit flags it
  rather than "matching" a flaw.
- **Parity acceptance bar** (for later sub-projects, recorded here): numeric +
  visual equivalence. Underlying numbers match the reference within declared
  tolerances AND panels are visually equivalent (same terms, genes, layout,
  encodings). Pixel-perfection is not required.
- **Parallel execution** of the audit is approved.
- **Inline fixes allowed.** The audit is read-first, but trivial/obvious gaps
  may be fixed as found; each such fix is recorded in the register with a
  `fixed-inline` status and a commit reference. Non-trivial fixes are deferred
  to sub-projects #2–#4.

## Purpose & output

Produce a single authoritative deliverable — an **audit report + a
machine-readable gap register** — that:

1. maps every `bulk-rna-frame` stage and figure to its lymphatic reference,
2. classifies each divergence,
3. ranks findings by severity,
4. routes each finding to the sub-project that will close it, and
5. crystallizes into the **contract spec** that #2–#3 implement against.

## Scope

**In scope:** reading and comparing both codebases; running the existing diff
bridge to produce evidence; writing the report + register; small extensions to
the diff bridge; inline fixes of trivial/obvious gaps.

**Out of scope (deferred to #2–#4):** porting the 2,937-line figure engine,
rewriting stages, building the full numeric-parity harness, non-trivial
statistical corrections. Human-species parity (mouse-first for now).

## Three audit tracks

### Track A — Stage / numeric fidelity
For each stage — `qc, de, pathways, ontology, composition, regulators, network,
hypothesis` plus advanced `sva, wgcna, mediation, power, multilayer` — compare
the reference modules (`workflow/stages/<stage>/{io,plots,tables,theme,stage}.R`)
against the frame's single script (`workflow/scripts/<stage>.R|.py`). Assess:

- statistical method and library (e.g. DESeq2 apeglm, fgsea, VIPER, DoRothEA);
- parameters, thresholds, and defaults;
- random seeds and determinism;
- edge-case handling (small-n, zero-variance genes, all-NA, ties, empty sets);
- contrast semantics (numerator − denominator must be explicit, never inferred).

### Track B — Figure fidelity
For each of the **16 gold panels**, compare the reference plotting code
(the relevant `publication_figure_functions.R` section + the panel wrapper in
`analysis/cape_xizhao_edits/`) against the frame constructor in `figures.py` /
`publication.R`. Assess: data source, transforms, geoms/encodings, faceting,
color/shape/size scales, legends, annotations, theme, page geometry. Flag panels
absent from the frame entirely.

### Track C — Config / hypothesis generality + reproducibility infra
- Can `config.py` + the hypothesis and figure-recipe schemas flex to a *new*
  cohort, or are defaults implicitly CAPE-shaped?
- Audit network providers/caching, `verify`, release manifest, seeds, and
  container definitions for genuine end-to-end reproducibility.

## Classification & severity

Every finding carries exactly one **verdict**:

- `REGRESSION` — frame diverges from the reference in a worse/incorrect way.
- `FIDELITY-GAP` — frame is missing a reference capability or panel.
- `BEST-PRACTICE` — the reference itself is questionable; the frame should
  improve on it (per "best practice wins").
- `OK` — faithful / statistically equivalent.

…and exactly one **severity**: `S1` statistical correctness → `S2` numeric
parity → `S3` figure fidelity → `S4` cohort generality → `S5` polish.

Each finding also records the **target sub-project** for its fix
(#2 harness / #3 figures / #4 stats / inline), and — if fixed inline — a
`fixed-inline` status plus commit ref.

## Evidence standard

Findings cite `file:line` on **both** sides. Numeric claims are backed by an
**actual re-run**, not code-reading: reuse
`analysis/bulk_rna_frame/{runner.py, compare_counts.py, compare_de.py}` to run
the lymphatic dataset through `bulk-rna-frame` and diff against the reference
`results/`. Where a stage lacks a comparator, note the gap and, if cheap, extend
the bridge; otherwise mark the finding `code-only` (verdict from reading, not a
run) so confidence is explicit.

## Execution method

Fan out **parallel audit agents**, one per stage cluster / figure cluster / infra
area, each returning a structured finding set (verdict, severity, evidence,
target sub-project). A synthesis step dedupes and ranks into the register. The
gold PDFs are the visual ground truth for Track B; the lymphatic `results/` are
the numeric ground truth for Track A.

## Deliverables

- `docs/superpowers/specs/2026-08-12-parity-audit-design.md` — this design.
- `docs/parity-audit.md` — human-readable audit report (build phase).
- `docs/parity-gaps.yaml` — machine-readable, ranked gap register (build phase).
- Any inline fixes, each committed and referenced from the register.

## Acceptance criteria

- Every stage (13) and every gold panel (16) has a verdict with two-sided
  `file:line` evidence.
- Numeric verdicts are backed by a real diff wherever the bridge allows;
  `code-only` findings are explicitly marked.
- The register is ranked by severity and each item is routed to a fix
  sub-project (or marked `fixed-inline`).
- Config/generality and reproducibility-infra findings are captured (Track C).

## Testing

The audit's numeric claims are validated by executing the diff bridge on the
lymphatic dataset — findings are evidence-based, not speculative. Any inline fix
must keep the existing `bulk-rna-frame` test suite green
(`pytest`) and, where relevant, add a regression test.

## Risks & controls

- **Scope creep into fixing.** Controlled by the inline-fix rule: only
  trivial/obvious fixes inline; everything else routed to #2–#4.
- **False parity.** "Matching" a reference flaw is prevented by the 3-way
  best-practice judgment and the `BEST-PRACTICE` verdict.
- **Environment drift** in re-runs. Controlled by using the frame's pinned
  envs/containers for the bridge run and recording the environment in evidence.
