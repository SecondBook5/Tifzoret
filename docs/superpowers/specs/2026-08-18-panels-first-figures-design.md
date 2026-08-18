# Panels-first, dataviz-validated figures — design

Date: 2026-08-18
Status: proposed (awaiting review)
Repo: `bulk-rna-frame` (engine). Palette configs in `lymphatic-flow-homeostasis`.
Supersedes in priority: the factorial-interaction work is paused until this ships.

## Problem

A visual audit of the assembled figures exposed systemic defects that the
end-to-end pipeline "pass" never checked:

1. **Composition collapse (every study).** The assembler fits each panel into a
   rigid 3-column grid *preserving aspect ratio, then centers it*. Any panel whose
   natural shape ≠ its cell collapses to a strip in whitespace. Measured
   cell-fill: cape D 44% W, H 29% W, I 28% W; A/C/E 64–67% H; and equivalents in
   all 5 studies (PTEN H 22% W, I 20% W).
2. **Detached legend (panel D).** The DE-heatmap program color strip is given a
   6%-wide patchwork column that also has to hold the strip's own "Program"
   legend, so the color bar detaches from the rows it annotates and the legend
   floats in the gap.
3. **Wrong form (panels H, I).** `program_violins` is a grid of ~40 micro-violins
   (unreadable at any placement, native 3900×12150 px); `program_integrated`
   crams a heatmap + forest into one cell.
4. **Palettes fail the dataviz validator (computed, not eyeballed):**
   - condition colors are ColorBrewer-*Paired* light tints used alone — too light,
     near-gray chroma, contrast < 2:1 on white; cape's pair is normal-vision
     ΔE 14.2 (below the 15 hard floor).
   - the program palette is 8–9 classes (> 7 ceiling); the two ambers are ΔE 11.1
     (normal vision); gray↔magenta is ΔE 2.2 (invisible to deutan viewers).

Root process failure: dataviz **step 7 — "render it and look at it"** was skipped.

## Decision (chosen by the user)

**Panels-first.** Every individual panel becomes a publication-grade, validated
standalone artifact. The auto-assembled composite is demoted to a **rough contact
sheet** (a labeled thumbnail index), and the final multi-panel figure is assembled
by hand per journal — matching the existing Xizhao-edits workflow. We do **not**
build an aspect-aware auto-composite.

## Goals

1. Each promoted panel is individually correct: right aspect, attached legend,
   validated palette, a defensible *form*.
2. Color correctness is **operationalized** — a test runs the validator so a
   washed-out or CVD-unsafe palette can never regress in silently.
3. The composite is honestly labeled as an index, not a figure.
4. No change to any numeric result — this is a presentation-layer change only.

## Non-goals (YAGNI)

- Aspect-aware auto-composite / grid solver (user chose panels-first).
- Dark mode, hover, tooltips (static PDF/PNG publication figures).
- Re-running featureCounts (counts are cached; only the figure/report tail reruns).

## Validated palettes (computed with the skill's validator)

All values below were run through `validate_palette.py` and **pass** on the light
surface `#fcfcfb`. They replace the failing palettes.

### Two-group condition (control/CAPE and every pairwise contrast)
`denominator → #2a78d6` (blue, slot 1), `numerator → #eb6834` (orange, slot 2).
ALL PASS; worst CVD ΔE 24.7. Set in each two-group study's `figures.palette`.

### PTEN 2×2 (genotype = hue, tumor = shade) — preserves the semantic
| group | hex | meaning |
|---|---|---|
| WT_baseline | `#6da7ec` | blue, light = baseline |
| WT_tumor | `#184f95` | blue, dark = tumor |
| PTEN_baseline | `#ef8683` | red, light = baseline |
| PTEN_tumor | `#b02724` | red, dark = tumor |

ALL PASS under `--pairs all`; worst CVD ΔE 15.2. (The naive "4 distinct hues"
alternative FAILS: green↔orange ΔE 3.2 protan — so the 2-hue×2-shade scheme is
the colorblind-safe choice, not merely the prettier one.) The two light members
carry a contrast WARN → satisfied by the relief rule (PCA points and heatmap
columns already carry direct group/sample labels).

### Program annotation — cap at 7 meaningful classes + gray "Other"
Replace `default_colors` with the seven validated slots:
`#2a78d6, #eb6834, #1baf7a, #eda100, #e87ba4, #008300, #4a3aa7`
and reserve a single gray for the catch-all (dedupe today's two near-greys
`#8C8C8C`/`#969696` into one `#8C8C8C` = "Other"). ALL PASS; worst adjacent
CVD ΔE 9.1. Programs beyond seven fold into "Other" (the dataviz-sanctioned tail
treatment). The magenta/yellow/aqua contrast WARN is satisfied by relief: the
strip always ships with its legend and text-labeled gene rows, so identity is
never color-alone.

## Per-panel changes (all in `publication.R` unless noted)

### A. Palette plumbing
- Study configs (`figures.palette`) updated to the validated values above
  (`lymphatic-flow-homeostasis/studies/*/bulk_rna_frame.yaml`).
- `publication.R`: `default_colors` → the 7 validated slots; single "Other" gray;
  cap distinct programs at 7 (overflow → "Other").

### B. de_heatmap legend detachment (panel D) — `publication.R:109–137`
- Collect guides across the composition: `plot_layout(guides = "collect") &
  theme(legend.position = "right")`.
- The program strip becomes a flush thin column with **no embedded legend**
  (remove its per-subplot `legend.position = "right"`); both legends (Program +
  Row z-score) gather on the right.
- Tighten the strip width so tiles sit adjacent to the heatmap rows.

### C. program_violins (panel I) — wrong form — `publication.R:170–190`
- Stop rendering all ~40 genes as micro-violins. Facet a **curated subset**: the
  genes of the configured `gsea_programs` (fallback: top DE genes per program),
  **capped at 16**, at a readable per-facet size (fixed facet dimensions, not a
  height that grows unbounded with gene count).
- The complete per-gene table stays in `program_violins_displayed.tsv` — nothing
  is gated, exactly as the dataviz "table view" rule requires. `log()` what was
  capped so the truncation is explicit.

### D. program_integrated (panel H) — crammed — `publication.R:150–168`
- Emit the heatmap and the effects/forest as **two properly sized standalone
  panels** (`program_heatmap` and `program_effects`) rather than one crushed
  composite, so each is legible on its own. The combined `program_integrated`
  may remain for the contact sheet but is no longer a promoted panel.

### E. Composite → contact sheet — assembly rule/script + `front_door.py`
- Replace the fit-and-center grid assembler with a **contact-sheet** generator:
  uniform-width thumbnails at each panel's **native aspect ratio**, labeled
  A/B/C…, generous surface gaps, titled "Contact sheet — assemble final figure by
  hand." No forced cell-fill, no centering-in-oversized-cells.
- `front_door.py` continues to promote the **individual** panels (the deliverable);
  the contact sheet is promoted as a single clearly-named index image.
- `figure_recipe.yaml` still declares which panels and in what order; only the
  rendering target changes (thumbnail index, not fixed-cell composite).

## Testing

1. **Palette-validation test (new, operationalizes dataviz step 3).** A test that
   shells the validator (or ports its checks) over: each study's configured
   `figures.palette`, and the program palette. Asserts PASS (WARN allowed only
   where the relief rule is documented). This is the guard that stops silent
   regression.
2. **Panel-geometry test.** Each promoted panel's saved PDF/PNG aspect ratio is
   within a sane band (no 0.32-style collapse); contact-sheet thumbnails preserve
   native aspect.
3. **Numeric-invariance regression.** DE/pathway/ontology tables are byte-identical
   before/after — proves this is presentation-only.
4. **Step 7 — manual visual sign-off.** After re-render, open every promoted panel
   and the contact sheet for all 5 studies and confirm against the dataviz
   anti-patterns list before calling it done. This step is non-optional and is the
   direct fix for how these defects shipped.

## Rollout

Engine branch `figure-engine`. TDD per the tests above. Then update the 5 study
palettes and re-run all 5 (`counting` cached → only the analysis→figures→report
tail re-executes). Rebuild the Xizhao bundle from the corrected panels.

## Risks

- **Curating violins may drop a gene a reader wanted.** Mitigated: full table
  retained + explicit `log()` of what was capped.
- **Palette change alters every figure's colors.** Intended; the old colors fail
  the validator. The numeric-invariance test proves only color/layout moved.
- **Contact sheet mistaken for a final figure.** Mitigated by the explicit title
  and by it not being a promoted "panel."
