# Parity Evidence Capture

This directory contains numeric evidence from bridge comparison runs that diff bulk-rna-frame CLI outputs against reference results from the lymphatic-flow-homeostasis repository.

## Capture Method

Evidence was captured using `tools/parity/capture_evidence.sh`, which:
1. Checks for pre-existing comparison outputs in `lymphatic-flow-homeostasis/results/bulk_rna_frame/migration/`
2. If not found, attempts to run the bridge (`lymphatic-flow-homeostasis/analysis/bulk_rna_frame/runner.py`) to:
   - Validate the study configuration
   - Run the bulk-rna-frame CLI to generate counts and DE results
   - Compare outputs using `compare_counts.py` and `compare_de.py`
3. Captures stdout/stderr; does not fail on numeric mismatches (mismatches ARE the evidence)

## Environment

- **bulk-rna CLI**: Installed and available (`bulk-rna --help` succeeds)
- **Data root**: `/mnt/e/rnaseq_projects` (exists)
- **GTF**: `/mnt/e/rnaseq_projects/Mus_musculus.GRCm39.107.gtf` (exists)
- **Reference repository**: `/home/ajbook/projects/lymphatic-flow-homeostasis`
- **Capture date**: 2026-08-12

## Results by Cohort

### cape_thoracic_duct ✓

**Status**: SUCCESS (pre-existing comparison outputs)

Pre-existing comparison outputs were found in the migration directory and copied:
- `counts_equivalence.json`: All samples exactly equal (0 mismatched genes, 0 maximum absolute difference)
- `de_equivalence.json`: All fields within tolerance (0 mismatched genes, 0 decision mismatches)

**Evidence quality**: `evidence: run` (actual numeric comparison from prior bridge run)

### ligation_vs_sham ✗

**Status**: BLOCKED (validation failed)

**Reason**: Unsupported project configuration version (None)

The study configuration at `lymphatic-flow-homeostasis/studies/ligation_vs_sham/config.yaml` lacks a `version` field and uses schema elements no longer supported by the current bulk-rna-frame CLI. Validation failed with exit code 2.

**Evidence quality**: `evidence: code-only` (fallback to static code inspection)

**Reproducibility impact**: This cohort was included in published analyses but cannot be re-run through the current bulk-rna-frame CLI without config migration. This represents an S2 reproducibility gap.

### old_vs_young ✗

**Status**: BLOCKED (validation failed)

**Reason**: Unsupported project configuration version (None)

The study configuration at `lymphatic-flow-homeostasis/studies/old_vs_young/config.yaml` lacks a `version` field and uses schema elements no longer supported by the current bulk-rna-frame CLI. Validation failed with exit code 2.

**Evidence quality**: `evidence: code-only` (fallback to static code inspection)

**Reproducibility impact**: This cohort was included in published analyses but cannot be re-run through the current bulk-rna-frame CLI without config migration. This represents an S2 reproducibility gap.

### rela_ko_vs_wt ✗

**Status**: BLOCKED (validation failed)

**Reason**: Unsupported project configuration version (None)

The study configuration at `lymphatic-flow-homeostasis/studies/rela_ko_vs_wt/config.yaml` lacks a `version` field and uses schema elements no longer supported by the current bulk-rna-frame CLI. Validation failed with exit code 2.

**Evidence quality**: `evidence: code-only` (fallback to static code inspection)

**Reproducibility impact**: This cohort was included in published analyses but cannot be re-run through the current bulk-rna-frame CLI without config migration. This represents an S2 reproducibility gap.

## Implications for Track A

Track A findings (numeric/functional parity) will use the following evidence modes:
- **cape_thoracic_duct stages**: `evidence: run` (numeric diffs available)
- **ligation_vs_sham, old_vs_young, rela_ko_vs_wt stages**: `evidence: code-only` (fallback to static code inspection since numeric runs are blocked)

## Reproducibility Finding

The inability to re-run 3 of 4 published cohorts through the current CLI represents a **reproducibility regression**. This has been documented as a parity register entry in `docs/parity/fragments/evidence.yaml`:
- **ID**: evidence-001
- **Verdict**: FIDELITY-GAP
- **Severity**: S2 (moderate: workaround exists via config migration, but not automated)
- **Track**: A (numeric/functional parity)
- **Area**: config-migration
- **Target**: harness (CLI validation layer)

## Re-running Capture

To re-run evidence capture for a cohort:

```bash
bash tools/parity/capture_evidence.sh <cohort>
```

Example:
```bash
bash tools/parity/capture_evidence.sh cape_thoracic_duct
```

The script will:
- Create `docs/parity/evidence/<cohort>/` if it doesn't exist
- Write a timestamped `capture.log`
- Copy or generate comparison outputs
- Exit 0 even if blocked (blocking is expected evidence)
