#!/usr/bin/env bash
# Capture numeric parity evidence by running bridge comparison scripts.
# Usage: bash tools/parity/capture_evidence.sh <cohort>
# Example: bash tools/parity/capture_evidence.sh cape_thoracic_duct

set -euo pipefail

COHORT="${1:-}"
if [[ -z "$COHORT" ]]; then
    echo "Usage: $0 <cohort>" >&2
    echo "Available cohorts: cape_thoracic_duct, ligation_vs_sham, old_vs_young, rela_ko_vs_wt" >&2
    exit 1
fi

# Paths
FRAME_ROOT="/home/ajbook/projects/bulk-rna-frame"
LYMPH_ROOT="/home/ajbook/projects/lymphatic-flow-homeostasis"
EVIDENCE_DIR="${FRAME_ROOT}/docs/parity/evidence/${COHORT}"
BRIDGE_DIR="${LYMPH_ROOT}/analysis/bulk_rna_frame"
STUDY_CONFIG="${LYMPH_ROOT}/studies/${COHORT}/config.yaml"
RUNNER="${BRIDGE_DIR}/runner.py"

# Validate cohort exists
if [[ ! -f "$STUDY_CONFIG" ]]; then
    echo "ERROR: Study config not found: ${STUDY_CONFIG}" >&2
    exit 1
fi

# Create evidence directory
mkdir -p "${EVIDENCE_DIR}"

echo "=== Capturing evidence for ${COHORT} ===" | tee "${EVIDENCE_DIR}/capture.log"
echo "Study config: ${STUDY_CONFIG}" | tee -a "${EVIDENCE_DIR}/capture.log"
echo "Timestamp: $(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "${EVIDENCE_DIR}/capture.log"
echo "" | tee -a "${EVIDENCE_DIR}/capture.log"

# Check if pre-existing comparison outputs exist in migration directory
MIGRATION_COUNTS="${LYMPH_ROOT}/results/bulk_rna_frame/migration/${COHORT}_counts_equivalence.json"
MIGRATION_DE="${LYMPH_ROOT}/results/bulk_rna_frame/migration/${COHORT}_de_equivalence.json"

if [[ -f "$MIGRATION_COUNTS" ]] && [[ -f "$MIGRATION_DE" ]]; then
    echo "Found pre-existing comparison outputs in migration directory" | tee -a "${EVIDENCE_DIR}/capture.log"
    echo "Copying: ${MIGRATION_COUNTS}" | tee -a "${EVIDENCE_DIR}/capture.log"
    echo "Copying: ${MIGRATION_DE}" | tee -a "${EVIDENCE_DIR}/capture.log"
    cp "$MIGRATION_COUNTS" "${EVIDENCE_DIR}/counts_equivalence.json"
    cp "$MIGRATION_DE" "${EVIDENCE_DIR}/de_equivalence.json"
    echo "SUCCESS: Pre-existing comparison outputs copied" | tee -a "${EVIDENCE_DIR}/capture.log"
    exit 0
fi

# Otherwise, attempt to run the bridge
echo "No pre-existing comparison outputs found; attempting fresh bridge run" | tee -a "${EVIDENCE_DIR}/capture.log"

# Step 1: Validate the study config
echo "" | tee -a "${EVIDENCE_DIR}/capture.log"
echo "--- Step 1: Validating study config ---" | tee -a "${EVIDENCE_DIR}/capture.log"
set +e
python3 "$RUNNER" validate "$STUDY_CONFIG" 2>&1 | tee "${EVIDENCE_DIR}/validate.txt"
VALIDATE_EXIT=$?
set -e
echo "Validation exit code: ${VALIDATE_EXIT}" | tee -a "${EVIDENCE_DIR}/capture.log"

if [[ $VALIDATE_EXIT -ne 0 ]]; then
    echo "BLOCKED: Validation failed for ${COHORT}" | tee -a "${EVIDENCE_DIR}/capture.log"
    exit 0  # Don't fail the script; this is expected evidence
fi

# Step 2: Run prepare to generate canonical counts
echo "" | tee -a "${EVIDENCE_DIR}/capture.log"
echo "--- Step 2: Running prepare (canonical counts) ---" | tee -a "${EVIDENCE_DIR}/capture.log"
set +e
python3 "$RUNNER" prepare "$STUDY_CONFIG" --cores 1 --no-conda 2>&1 | tee "${EVIDENCE_DIR}/prepare.txt"
PREPARE_EXIT=$?
set -e
echo "Prepare exit code: ${PREPARE_EXIT}" | tee -a "${EVIDENCE_DIR}/capture.log"

if [[ $PREPARE_EXIT -ne 0 ]]; then
    echo "BLOCKED: Prepare failed for ${COHORT}" | tee -a "${EVIDENCE_DIR}/capture.log"
    exit 0
fi

# Step 3: Run full pipeline to generate DE results
echo "" | tee -a "${EVIDENCE_DIR}/capture.log"
echo "--- Step 3: Running full pipeline ---" | tee -a "${EVIDENCE_DIR}/capture.log"
set +e
python3 "$RUNNER" run "$STUDY_CONFIG" --cores 1 --no-conda 2>&1 | tee "${EVIDENCE_DIR}/run.txt"
RUN_EXIT=$?
set -e
echo "Run exit code: ${RUN_EXIT}" | tee -a "${EVIDENCE_DIR}/capture.log"

if [[ $RUN_EXIT -ne 0 ]]; then
    echo "BLOCKED: Run failed for ${COHORT}" | tee -a "${EVIDENCE_DIR}/capture.log"
    exit 0
fi

# Step 4: Compare counts
echo "" | tee -a "${EVIDENCE_DIR}/capture.log"
echo "--- Step 4: Comparing counts ---" | tee -a "${EVIDENCE_DIR}/capture.log"

# Need to locate reference and candidate counts files
# Reference: typically in ${LYMPH_ROOT}/results/${COHORT}/{full,primary}/*.counts.txt
# Candidate: in ${LYMPH_ROOT}/studies/${COHORT}/results/tables/counts.tsv
# Samples: ${LYMPH_ROOT}/studies/${COHORT}/samples.tsv

SAMPLES_TSV="${LYMPH_ROOT}/studies/${COHORT}/samples.tsv"
CANDIDATE_COUNTS="${LYMPH_ROOT}/studies/${COHORT}/results/tables/counts.tsv"

# Find the reference counts file (look for featureCounts output)
REFERENCE_COUNTS=""
for subdir in full primary; do
    REF_DIR="${LYMPH_ROOT}/results/${COHORT}/${subdir}"
    if [[ -d "$REF_DIR" ]]; then
        # Look for .counts.txt files
        COUNTS_FILE=$(find "$REF_DIR" -name "*.counts.txt" -o -name "*featureCounts*" | head -1 || echo "")
        if [[ -n "$COUNTS_FILE" ]]; then
            REFERENCE_COUNTS="$COUNTS_FILE"
            break
        fi
    fi
done

if [[ -z "$REFERENCE_COUNTS" ]] || [[ ! -f "$REFERENCE_COUNTS" ]]; then
    echo "BLOCKED: Could not locate reference counts file for ${COHORT}" | tee -a "${EVIDENCE_DIR}/capture.log"
    exit 0
fi

if [[ ! -f "$SAMPLES_TSV" ]]; then
    echo "BLOCKED: samples.tsv not found: ${SAMPLES_TSV}" | tee -a "${EVIDENCE_DIR}/capture.log"
    exit 0
fi

if [[ ! -f "$CANDIDATE_COUNTS" ]]; then
    echo "BLOCKED: Candidate counts not found: ${CANDIDATE_COUNTS}" | tee -a "${EVIDENCE_DIR}/capture.log"
    exit 0
fi

# Determine analysis set (usually "primary" or the contrast name)
ANALYSIS_SET="primary"

set +e
python3 "${BRIDGE_DIR}/compare_counts.py" \
    --samples "$SAMPLES_TSV" \
    --analysis-set "$ANALYSIS_SET" \
    --reference "$REFERENCE_COUNTS" \
    --candidate "$CANDIDATE_COUNTS" \
    --output "${EVIDENCE_DIR}/counts_equivalence.json" \
    2>&1 | tee "${EVIDENCE_DIR}/compare_counts.txt"
COUNTS_EXIT=$?
set -e
echo "Counts comparison exit code: ${COUNTS_EXIT}" | tee -a "${EVIDENCE_DIR}/capture.log"

# Step 5: Compare DE results
echo "" | tee -a "${EVIDENCE_DIR}/capture.log"
echo "--- Step 5: Comparing DE results ---" | tee -a "${EVIDENCE_DIR}/capture.log"

# Find reference DE files and candidate DE file
# Reference raw: typically results/${COHORT}/{full,primary}/de_results_raw.tsv
# Reference shrunken: typically results/${COHORT}/{full,primary}/de_results_shrunken.tsv
# Candidate: studies/${COHORT}/results/analyses/<contrast>/de.tsv

# First, determine the contrast name
CONTRAST_DIR=$(find "${LYMPH_ROOT}/studies/${COHORT}/results/analyses" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | head -1 || echo "")
if [[ -z "$CONTRAST_DIR" ]] || [[ ! -d "$CONTRAST_DIR" ]]; then
    echo "BLOCKED: Could not locate contrast directory for ${COHORT}" | tee -a "${EVIDENCE_DIR}/capture.log"
    exit 0
fi

CANDIDATE_DE="${CONTRAST_DIR}/de.tsv"
if [[ ! -f "$CANDIDATE_DE" ]]; then
    echo "BLOCKED: Candidate DE file not found: ${CANDIDATE_DE}" | tee -a "${EVIDENCE_DIR}/capture.log"
    exit 0
fi

# Find reference DE files
REFERENCE_RAW=""
REFERENCE_SHRUNKEN=""
for subdir in full primary; do
    REF_DIR="${LYMPH_ROOT}/results/${COHORT}/${subdir}"
    if [[ -d "$REF_DIR" ]]; then
        RAW=$(find "$REF_DIR" -name "*raw*.tsv" -o -name "*de_raw*" | head -1 || echo "")
        SHRUNKEN=$(find "$REF_DIR" -name "*shrunken*.tsv" -o -name "*de_shrunken*" | head -1 || echo "")
        if [[ -n "$RAW" ]] && [[ -n "$SHRUNKEN" ]]; then
            REFERENCE_RAW="$RAW"
            REFERENCE_SHRUNKEN="$SHRUNKEN"
            break
        fi
    fi
done

if [[ -z "$REFERENCE_RAW" ]] || [[ ! -f "$REFERENCE_RAW" ]]; then
    echo "BLOCKED: Could not locate reference raw DE file for ${COHORT}" | tee -a "${EVIDENCE_DIR}/capture.log"
    exit 0
fi

if [[ -z "$REFERENCE_SHRUNKEN" ]] || [[ ! -f "$REFERENCE_SHRUNKEN" ]]; then
    echo "BLOCKED: Could not locate reference shrunken DE file for ${COHORT}" | tee -a "${EVIDENCE_DIR}/capture.log"
    exit 0
fi

set +e
python3 "${BRIDGE_DIR}/compare_de.py" \
    --reference-raw "$REFERENCE_RAW" \
    --reference-shrunken "$REFERENCE_SHRUNKEN" \
    --candidate "$CANDIDATE_DE" \
    --output "${EVIDENCE_DIR}/de_equivalence.json" \
    2>&1 | tee "${EVIDENCE_DIR}/compare_de.txt"
DE_EXIT=$?
set -e
echo "DE comparison exit code: ${DE_EXIT}" | tee -a "${EVIDENCE_DIR}/capture.log"

# Summary
echo "" | tee -a "${EVIDENCE_DIR}/capture.log"
echo "=== Capture complete ===" | tee -a "${EVIDENCE_DIR}/capture.log"
echo "Evidence directory: ${EVIDENCE_DIR}" | tee -a "${EVIDENCE_DIR}/capture.log"
if [[ -f "${EVIDENCE_DIR}/counts_equivalence.json" ]]; then
    echo "Counts comparison: COMPLETED" | tee -a "${EVIDENCE_DIR}/capture.log"
else
    echo "Counts comparison: FAILED or BLOCKED" | tee -a "${EVIDENCE_DIR}/capture.log"
fi
if [[ -f "${EVIDENCE_DIR}/de_equivalence.json" ]]; then
    echo "DE comparison: COMPLETED" | tee -a "${EVIDENCE_DIR}/capture.log"
else
    echo "DE comparison: FAILED or BLOCKED" | tee -a "${EVIDENCE_DIR}/capture.log"
fi

exit 0
