#!/bin/bash
#SBATCH --job-name=af_concordance
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --partition=norm
#
# Step 6 stage B, run standalone: re-derive the AF-concordance exclusion list at different knobs
# against an existing stage-A output. Step 6 already runs this in-job; after tuning here, re-run
# step 6 with SKIP_AF_BUILD=1 to apply the result. Method: af_concordance_build.py, METHODS.md §6.
#
# Input must be the UNFILTERED fileset: a list derived from already-filtered input would look
# clean whether or not it was.
#
# Knobs: THRESH (0.05) · ZMIN (5.0) · MIN_CELL (100) · DISC_RATE (0.50) · ANCS · DXS
#        HWE (1e-4, 0 disables) · MIN_HWE_CONTROLS (50) · HWE_BOTH_TAILS (unset)
#        HWE_REQUIRE_EXCESS (1; 0 = ungated union, which reproduces the 4,415-variant list)
#        MISHAP (0 = off; set 1e-4 to enable) · MIN_MISHAP (100)
#
# GUARDRAIL: run by the user; reads the ID-bearing annotation and genotypes. Writes the exclusion
# list; stdout is aggregate counts and variant IDs only.
set -o pipefail

# Resolve the bundle root and load the one file that knows where anything lives.
BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

QC_DIR=${QC_DIR:-${MERGED_DIR}/by_ancestry_qc/unfiltered}
WORK=${MERGED_DIR}/af_concordance
OUT=${AF_EXCLUDE:-${MERGED_DIR}/exclude_af_concordance.txt}
DISC=${DISC:-${MERGED_DIR}/concordance/per_variant_discordance.tsv}

THRESH=${THRESH:-0.05}
ZMIN=${ZMIN:-5.0}
MIN_CELL=${MIN_CELL:-100}
DISC_RATE=${DISC_RATE:-0.50}
# Stage 2 — per-callset HWE among controls (1e-4 is GenoTools' default; HWE=0 skips it).
HWE=${HWE:-1e-4}
MIN_HWE_CONTROLS=${MIN_HWE_CONTROLS:-50}
# Stage 3 — per-callset haplotype missingness (plink1.9 --test-mishap). Off by default: slow, and
# missingness is not the driver here.
MISHAP=${MISHAP:-0}
MIN_MISHAP=${MIN_MISHAP:-100}
# An HWE cell contributes only if its rejections exceed chance at $HWE; 0 = ungated union.
HWE_REQUIRE_EXCESS=${HWE_REQUIRE_EXCESS:-1}

# Both plink generations by absolute path: they are one module family here, so loading MOD_PLINK1
# unloads MOD_PLINK2 and a bare `plink2` then resolves to a non-executable file.
module load "${MOD_PLINK2}"; PLINK2_BIN=$(command -v plink2)
module load "${MOD_PLINK1}"; PLINK1_BIN=$(command -v plink)
[[ -x "$PLINK2_BIN" ]] || { echo "ERROR: plink2 not found via ${MOD_PLINK2}" >&2; exit 1; }
[[ -x "$PLINK1_BIN" ]] || { echo "ERROR: plink not found via ${MOD_PLINK1}" >&2; exit 1; }
module load "${MOD_PYTHON}"
source "${VENV}/bin/activate"

echo "=========================================="
echo "AF-concordance exclusion list (step 6 stage B, run standalone)"
echo "Job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}   start $(date)"
echo "QC dir: ${QC_DIR}   (must be the UNFILTERED generation)"
echo "Annot:  ${ANNOT}"
echo "Out:    ${OUT}"
echo "=========================================="

[[ -f "$ANNOT" ]] || { echo "ERROR: sample annot not found: $ANNOT — run clinical_core.py (§12a)" >&2; exit 1; }
[[ -d "$QC_DIR" ]] || { echo "ERROR: ${QC_DIR} not found — run step 6 (stage A builds it)" >&2; exit 1; }

python3 "${BUNDLE}/scripts/af_concordance_build.py" \
    --qc-dir "${QC_DIR}" \
    --annot "${ANNOT}" \
    --out "${OUT}" \
    --work "${WORK}" \
    --plink1 "${PLINK1_BIN}" \
    --plink2 "${PLINK2_BIN}" \
    --thresh "${THRESH}" \
    --zmin "${ZMIN}" \
    --min-cell "${MIN_CELL}" \
    --hwe "${HWE}" \
    --min-hwe-controls "${MIN_HWE_CONTROLS}" \
    ${HWE_BOTH_TAILS:+--hwe-both-tails} \
    $([[ "${HWE_REQUIRE_EXCESS}" == "0" ]] && echo --no-hwe-require-excess) \
    --mishap "${MISHAP}" \
    --min-mishap "${MIN_MISHAP}" \
    --discordance "${DISC}" \
    --disc-rate "${DISC_RATE}" \
    ${ANCS:+--ancs ${ANCS}} \
    ${DXS:+--dx ${DXS}}
RC=$?

echo ""
echo "=========================================="
echo "Done $(date)  (exit ${RC})"
echo "=========================================="
exit ${RC}
