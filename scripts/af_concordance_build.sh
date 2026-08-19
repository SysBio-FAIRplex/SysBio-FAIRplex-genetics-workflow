#!/bin/bash
#SBATCH --job-name=af_concordance
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --partition=norm
#
# STEP 6a — build the per-callset AF-concordance exclusion list that step 6 then applies.
#
# Runs on step 6's EXISTING output, so it needs no re-run to produce the list. Two stages, unioned:
#   1. AF concordance — callsets compared only WITHIN a (stratum x dx) cell, which holds disease
#      constant and is what keeps a genome-wide frequency filter from deleting APOE.
#   2. Per-callset HWE among controls — GenoTools' `hwe` step at its own threshold. Lives here
#      rather than in step 6 because it needs the grain for control labels, and step 6 cannot
#      depend on the grain without a circular ordering (grain <- §12 <- manifest <- step 6).
#   3. Per-callset haplotype missingness — GenoTools' `haplotype` step. OFF by default (MISHAP=0).
# See the docstring in af_concordance_build.py.
#
# Together these are a POST-MERGE equivalent of GenoTools --all_variant, minus its `case_control`
# step, which step 7 already covers per contrast with a stricter rule (chi-square AND a minimum
# |F_MISS| difference, vs the demo's single global p<0.01 pass).
#
# ORDER:
#   1. ./submit.sh scripts/af_concordance_build.sh        <- this, on the CURRENT step-6 output
#   2. read the BY CALLSET PAIR table + the sentinel-loci tripwire in the log
#   3. ./submit.sh scripts/06_ancestry_qc.sh            <- picks the list up automatically
#   4. notebook §12 (PCs changed -> the grain must be rebuilt)
#   5. ./submit.sh scripts/07_gwas.sh
#
# Knobs: THRESH (0.05) · ZMIN (5.0) · MIN_CELL (100) · DISC_RATE (0.50) · ANCS · DXS
#        HWE (1e-4, 0 disables) · MIN_HWE_CONTROLS (50) · HWE_BOTH_TAILS (unset)
#        MISHAP (0 = off; set 1e-4 to enable) · MIN_MISHAP (100)
#
# GUARDRAIL: sbatch script, run by the USER (it reads the id-bearing grain and genotypes).
# Writes the exclusion list; stdout is aggregate counts and variant IDs only.
set -o pipefail

# Resolve the bundle root and load the one file that knows where anything lives.
BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

QC_DIR=${MERGED_DIR}/by_ancestry_qc
WORK=${MERGED_DIR}/af_concordance
OUT=${AF_EXCLUDE:-${MERGED_DIR}/exclude_af_concordance.txt}
DISC=${DISC:-${MERGED_DIR}/concordance/per_variant_discordance.tsv}

THRESH=${THRESH:-0.05}
ZMIN=${ZMIN:-5.0}
MIN_CELL=${MIN_CELL:-100}
DISC_RATE=${DISC_RATE:-0.50}
# Stage 2 — per-callset HWE among controls, the GenoTools --all_variant equivalent. 1e-4 is
# GenoTools' default. Set HWE=0 to skip the stage. HWE_BOTH_TAILS=1 drops keep-fewhet.
HWE=${HWE:-1e-4}
MIN_HWE_CONTROLS=${MIN_HWE_CONTROLS:-50}
# Stage 3 — per-callset haplotype missingness (GenoTools' `haplotype`; plink1.9 --test-mishap).
# OFF by default: it is the slow one, and our own measurement says missingness is not the driver
# here. Run MISHAP=1e-4 once for completeness / demo equivalence, not on the critical path.
MISHAP=${MISHAP:-0}
MIN_MISHAP=${MIN_MISHAP:-100}

module load "${MOD_PLINK2}"
[[ "${MISHAP}" != "0" ]] && module load "${MOD_PLINK1}"   # --test-mishap is plink1.9 only
module load "${MOD_PYTHON}"
source "${VENV}/bin/activate"

echo "=========================================="
echo "Step 6a — per-callset AF-concordance exclusion list"
echo "Job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}   start $(date)"
echo "QC dir: ${QC_DIR}"
echo "Grain:  ${GRAIN}"
echo "Out:    ${OUT}"
echo "=========================================="

[[ -f "$GRAIN" ]] || { echo "ERROR: grain not found: $GRAIN (rsync it from local first)" >&2; exit 1; }
[[ -d "$QC_DIR" ]] || { echo "ERROR: no by_ancestry_qc — run step 6 at least once first" >&2; exit 1; }

python3 "${BUNDLE}/scripts/af_concordance_build.py" \
    --qc-dir "${QC_DIR}" \
    --grain "${GRAIN}" \
    --out "${OUT}" \
    --work "${WORK}" \
    --thresh "${THRESH}" \
    --zmin "${ZMIN}" \
    --min-cell "${MIN_CELL}" \
    --hwe "${HWE}" \
    --min-hwe-controls "${MIN_HWE_CONTROLS}" \
    ${HWE_BOTH_TAILS:+--hwe-both-tails} \
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
