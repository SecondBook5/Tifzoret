# Migration plan

## Phase 0: reference preservation

- Preserve the current lymphatic project and its successful result manifests.
- Record numerical and figure-level acceptance fixtures.
- Do not move project-specific biological annotations into the engine.

## Phase 1: neutral MVP

- Strict configuration and input validation.
- Delivered-BAM input with uniform featureCounts quantification.
- nf-core/rnaseq result-directory adapter at the BAM boundary.
- Count-matrix bypass and arbitrary covariates.
- Multiple named contrasts.
- QC, DE, custom-GMT ORA, fgsea, ssGSEA, and publication exports.
- Synthetic end-to-end test project.

## Phase 2: reusable manuscript figures

Promote the refined PCA/correlation, cell-state, DE heatmap variants, combined
ORA, GSVA, GSEA, and consolidated program panels behind generic interfaces.

## Phase 3: species adapters

Add explicit mouse and human providers for annotation, MSigDB, KEGG, STRING,
and regulator resources. Reject unsupported combinations instead of silently
falling back across species.

## Phase 4: networks and regulators

Promote connected STRING networks, regulator activity, rectangular GRNs, and
radial program-aware GRNs with complete node/edge audit tables.

## Phase 5: story recipes and assembly

Add declarative panel selection, dimensions, legend policy, alternative panel
variants, and vector multi-panel assembly.

## Phase 6: release

- Synthetic CI run.
- Mouse and human fixtures.
- Numerical regression against all configured lymphatic studies.
- Locked environments or containers.
- Versioned release and migration guide.
