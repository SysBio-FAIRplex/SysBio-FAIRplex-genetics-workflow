#!/bin/bash
#SBATCH --job-name=af_concordance
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --partition=norm
#
# Build the per-callset AF-concordance exclusion list.
#
# ─────────────────────────────────────────────────────────────────────────────
# YOU PROBABLY DO NOT NEED TO SUBMIT THIS. Step 6 runs it as stage B, in-job, between its
# unfiltered QC and the filtered association set. That is the supported path and it is a single
# submission. This wrapper survives for one purpose: re-deriving the list at DIFFERENT knobs
# (THRESH, ZMIN, MIN_CELL, HWE, MISHAP) against an existing stage-A output, without paying for
# another scan of cohort_merged. Having tuned them, re-run step 6 with SKIP_AF_BUILD=1 to apply
# the result.
#
# The file was formerly numbered 06a, and the ORDER block here used to read "run this, then run
# step 6 again". Both encoded a two-pass pipeline that no longer exists — see the header of
# 06_ancestry_qc.sh for why the ordering it was built around was never actually circular.
# ─────────────────────────────────────────────────────────────────────────────
#
# Three stages, unioned:
#   1. AF concordance — callsets compared only WITHIN a (stratum x dx) cell, which holds disease
#      constant and is what keeps a genome-wide frequency filter from deleting APOE.
#   2. Per-callset HWE among controls — GenoTools' `hwe` step at its own threshold. It needs the
#      per-sample control label, which is why this is a separate stage from step 6's pooled HWE.
#   3. Per-callset haplotype missingness — GenoTools' `haplotype` step. OFF by default (MISHAP=0).
# See the docstring in af_concordance_build.py.
#
# Together these are a POST-MERGE equivalent of GenoTools --all_variant, minus its `case_control`
# step, which step 7 already covers per contrast with a stricter rule (chi-square AND a minimum
# |F_MISS| difference, vs the demo's single global p<0.01 pass).
#
# INPUT IS THE UNFILTERED FILESET, and that is load-bearing: deriving the list from a directory
# with a previous list already applied would pre-filter this filter's own input, and the result
# would look clean whether or not it was.
#
# Knobs: THRESH (0.05) · ZMIN (5.0) · MIN_CELL (100) · DISC_RATE (0.50) · ANCS · DXS
#        HWE (1e-4, 0 disables) · MIN_HWE_CONTROLS (50) · HWE_BOTH_TAILS (unset)
#        MISHAP (0 = off; set 1e-4 to enable) · MIN_MISHAP (100)
#
# GUARDRAIL: sbatch script, run by the USER (it reads the id-bearing annotation and genotypes).
# Writes the exclusion list; stdout is aggregate counts and variant IDs only.
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
