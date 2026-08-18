# Factorial interaction contrasts — design

Date: 2026-08-18
Status: proposed (awaiting review)
Repo: `bulk-rna-frame` (engine). Driver-side change in `lymphatic-flow-homeostasis`.

## Problem

The engine models every contrast as a **pairwise group comparison**: a
`contrasts.tsv` row names a `factor`, a `numerator` level, and a `denominator`
level, and DE extracts the apeglm-shrunken coefficient
`factor_numerator_vs_denominator` from a single **global** design formula
(`analysis.design`, e.g. `~ condition`).

This cannot express a **difference-in-differences (interaction) contrast** —
the question "does deleting PTEN change the *tumor response* differently than it
changes in WT?" That is the interaction coefficient `genotypePTEN.tumoryes` from
the design `~ genotype * tumor`. It is not any single pairwise group difference;
it is `(PTEN_tumor − PTEN_baseline) − (WT_tumor − WT_baseline)`.

The flagship study (`lec_pten_tumor`, 2×2: genotype × tumor, 18 samples) needs
this contrast as the mechanistic centerpiece. The engine is meant to be the
universal bulk-RNA-seq tool, so factorial interaction must be a first-class,
config-driven capability — not a study-specific hack.

## Goals

1. Add **interaction / named-coefficient contrasts** as a configurable contrast
   type, expressible entirely in `contrasts.tsv`.
2. **Zero regression** for existing pairwise contrasts: their DE, pathways,
   figures, and tables must be byte-for-byte identical to today.
3. **Minimal blast radius** and a single source of truth for contrast
   resolution shared by DE and pathways.
4. Config validation that catches malformed interaction rows early with clear
   messages, without loosening pairwise validation.

## Non-goals (YAGNI)

- Arbitrary numeric contrast vectors (`c(0, 1, -1, ...)`).
- Likelihood-ratio-test contrasts / ANOVA-style multi-coefficient drops.
- Designs beyond what `resultsNames(dds)` already produces (no manual contrast
  algebra). If DESeq2 doesn't name the coefficient, we don't support it yet.
- More than one interaction term per contrast row.

## Schema: four optional `contrasts.tsv` columns

Existing columns are unchanged and remain required:
`contrast_id, factor, numerator, denominator, description`.

Four **optional** columns are added. Absent column or empty cell ⇒ pairwise
behavior, exactly as today.

| column | meaning | pairwise | coefficient |
|---|---|---|---|
| `type` | `pairwise` (default) or `coefficient` | `pairwise`/empty | `coefficient` |
| `design` | per-contrast design formula override | empty ⇒ global `analysis.design` | e.g. `~ genotype * tumor` |
| `reference_levels` | `;`-separated `factor=level` relevels applied before fit | empty | e.g. `genotype=WT;tumor=no` |
| `coefficient` | exact `resultsNames(dds)` entry to extract | empty | e.g. `genotypePTEN.tumoryes` |

For a `coefficient` contrast, `factor`/`numerator`/`denominator` become
**labels and metadata only** — they drive the direction enum naming
(`up_in_numerator` etc.), plot titles, and palette selection, but the effect
size itself comes from the named coefficient, not a level difference. Authors
set them to the biologically meaningful poles (numerator = the level whose
tumor-response is enhanced). This keeps every downstream consumer that reads
`factor/numerator/denominator` working unchanged.

Because `materialize_inputs.py::materialize_contrasts` copies **all** source
columns into the canonical `INPUTS/contrasts.tsv`, these columns reach `de.R`
and `pathways.R` with no rule/plumbing changes.

### Example rows (`lec_pten_tumor`)

```
contrast_id                factor    numerator  denominator  type         design              reference_levels     coefficient              description
pten_vs_wt_tumor           condition PTEN_tumor WT_tumor     pairwise                                                                        PTEN vs WT under tumor (primary pairwise)
pten_x_tumor_interaction   genotype  PTEN       WT           coefficient  ~ genotype * tumor   genotype=WT;tumor=no genotypePTEN.tumoryes    Difference-in-differences: does PTEN deletion change the tumor response vs WT?
```

(The four existing pairwise rows keep `type` empty; only the 5th row is new.)

## Architecture: one resolver, two branches

### New: `utils.R::resolve_contrast(contrast_row, global_design)`

Single source of truth for interpreting a contrast row. Returns a normalized
list consumed identically by `de.R` and `pathways.R`:

```
list(
  type              = "pairwise" | "coefficient",
  factor_name       = <chr>,      # always present (labels/direction/palette)
  numerator         = <chr>,
  denominator       = <chr>,
  design_formula    = <formula>,  # per-row `design` or global fallback
  reference_levels  = <named list factor -> level>,   # parsed from reference_levels
  coefficient_name  = <chr|NA>    # required when type == "coefficient"
)
```

Rules:
- `type`: `contrast_row$type` if the column exists and is non-empty, else
  `"pairwise"`. (A missing column reads as `NULL` in R, so guard with
  `is.null()`/`is.na()`/`nzchar()`.)
- `design_formula`: parse `contrast_row$design` when non-empty, else
  `as.formula(global_design)`.
- `reference_levels`: split on `;`, then `=`, into a named list; empty ⇒ empty
  list. Pairwise contrasts implicitly relevel `factor` to `denominator` (today's
  behavior) — the resolver encodes that so both branches share it.
- For `type == "coefficient"`, `coefficient_name` must be non-empty (validated
  upstream in `config.py`; resolver `stop()`s defensively if not).

### `de.R` — branch on `type`

- **pairwise (unchanged, must stay byte-identical):** current lines 33–35 relevel
  `factor` to `denominator` under the global design; coefficient pattern
  `^factor_numerator_vs_denominator$`; `results(dds, name=)` + `lfcShrink(coef=,
  type="apeglm")`. No behavioral change — ideally the diff to this branch is only
  that it now calls `resolve_contrast()` and reads fields off the returned list.
- **coefficient (new):**
  1. Factor every design variable; apply `reference_levels` relevels (so DESeq2
     names the coefficient the way the author wrote it).
  2. Build `dds` with the **per-contrast** `design_formula`.
  3. Verify `coefficient_name %in% resultsNames(dds)`; if not, `stop()` listing
     the available names (the single most likely author error — a mis-typed
     coefficient — gets an actionable message).
  4. `results(dds, name = coefficient_name)` + `lfcShrink(dds, coef =
     coefficient_name, type = "apeglm")`.
  5. Direction enum uses the same `up_in_numerator`/`down_in_numerator` semantics
     (positive LFC = numerator pole), so ORA/GSEA rankings are unaffected in shape.
  6. **Display set = all samples in the contrast's design**, not just
     `numerator`/`denominator` levels. Today `display_samples` (PCA at L182,
     heatmap at L198) filters to the two levels of `factor`; for an interaction
     across all four groups that would drop half the data. The resolver exposes
     which columns define the display set: pairwise ⇒ the two `factor` levels
     (unchanged); coefficient ⇒ all samples participating in `design_formula`.

### `pathways.R` — GSVA-differential branch (L284–313)

This is the only genuinely 2-group-specific pathways step (ORA/GSEA consume the
DE `direction`/ranking and need no change).

- **pairwise (unchanged):** relevel `contrast_group` to `denominator`, substitute
  `factor` → `contrast_group` in the design string, `model.matrix`, extract
  `^contrast_group{make.names(numerator)}$`, limma `lmFit`/`eBayes`.
- **coefficient (new):** apply `reference_levels`, build `model.matrix` from
  `design_formula`, extract the column matching `make.names(coefficient_name)`,
  same limma flow. Display set = all design samples (mirrors `de.R`).

### `config.py` — validation extension (current L424–450)

Keep all pairwise checks as-is. Branch on `type`:

- `type` (when present) must be `pairwise` or `coefficient`.
- **pairwise / absent:** unchanged — `factor` in samples header; numerator &
  denominator observed levels; `factor` appears in `analysis.design`.
- **coefficient:**
  - `coefficient` must be non-empty.
  - `design` (row override, or global if blank) must parse and every variable in
    it must be a column in `samples.tsv`.
  - `reference_levels`, if present, must parse as `factor=level` pairs where each
    factor is a design variable and each level is observed for that factor.
  - `factor`/`numerator`/`denominator` still required (labels); `factor` must be
    a samples column, but numerator/denominator are **not** required to be levels
    of `factor` and `factor` is **not** required to appear in the design (they're
    labels for an interaction whose poles may be a different variable).
  - We do **not** re-derive `resultsNames()` at config time (needs the fitted
    model); the mis-typed-coefficient case is caught at DE time with the
    available-names message. Document this explicitly.

## Driver-side change (`lymphatic-flow-homeostasis`)

- `studies/lec_pten_tumor/samples.tsv` already has `genotype` (WT/PTEN) and
  `tumor` (no/yes) — **no sample-sheet change needed.**
- Add the 5th row `pten_x_tumor_interaction` to
  `studies/lec_pten_tumor/contrasts.tsv` (see example above).
- Optionally add an interaction panel to `figure_recipe.yaml` later; not part of
  this spec's core.

## Testing

1. **`config.py` unit tests:** valid interaction row accepted; each malformed
   case rejected with its specific message (bad `type`; empty `coefficient`;
   `design` referencing a missing column; `reference_levels` with unobserved
   level / non-design factor). Existing pairwise validation tests still pass.
2. **`resolve_contrast()` unit test (R):** pairwise row → global design + relevel
   to denominator; coefficient row → per-row design, parsed reference levels,
   coefficient name.
3. **Integration test — tiny 2×2 synthetic fixture:** counts for genotype ×
   tumor with a known injected interaction. Assert the extracted
   `genotypePTEN.tumoryes` LFC matches the manual difference-in-differences of
   group means within tolerance, and that the display set includes all four
   groups.
4. **Pairwise regression / parity:** an existing pairwise study runs through the
   refactored `resolve_contrast()` path and produces byte-identical
   `de_results.tsv` / `gsva_differential.tsv` to the pre-change baseline.
5. **Constructor-recipe + front_door regression:** the `front_door.py` staged-
   panel path (recently fixed) stays green when an interaction contrast is
   present (guards the figure-promotion path against the new column set).

## Risks & mitigations

- **Silent wrong coefficient:** author typos `coefficient`. Mitigated by the
  DE-time `resultsNames()` membership check that prints the valid names.
- **Pairwise drift:** the refactor to route pairwise through `resolve_contrast()`
  could perturb output. Mitigated by test #4 (byte-identical parity gate) — this
  is the acceptance criterion for the refactor.
- **apeglm on interaction coef:** apeglm shrinks named coefficients fine; if a
  specific design makes the coefficient non-shrinkable, fall back is out of scope
  (fail loudly, don't silently switch estimators).

## Rollout

Engine change on branch `figure-engine` (current), TDD per the tests above,
then add the driver contrast row and re-run `lec_pten_tumor` (counting is cached;
only the ~15-job analysis→figures→report tail re-executes).
