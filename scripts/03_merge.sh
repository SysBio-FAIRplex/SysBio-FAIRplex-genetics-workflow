#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1500G
#SBATCH --partition=largemem
#
# STEP 3 — plink1.9 union merge of the three normalized bed filesets into one cohort.
#
# plink1.9 rather than plink2 because the non-concatenating --pmerge-list is unimplemented
# in this plink2 build. --allow-extra-chr covers the PAR codes plink1.9 does not recognize.
#
#   ./submit.sh scripts/03_merge.sh

set -o pipefail

# Resolve the bundle root and load the one file that knows where anything lives.
BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

module load "${MOD_PLINK1}"

# The three NORM_* bed stems come from config.sh.
OUT_DIR="${MERGED_DIR}"
MERGE_LIST="${OUT_DIR}/merge_list.txt"
OUT="${OUT_DIR}/cohort_merged"
TMP="${OUT_DIR}/tmp_merge"

mkdir -p "${OUT_DIR}" "${TMP}"

echo "=========================================="
echo "Cohort merge — union bed"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       ${SLURMD_NODENAME}"
echo "Start time: $(date)"
echo "=========================================="

# Build merge list (secondary filesets) — the step 2 NORMALIZED beds:
# chr-prefixed codes + reference-oriented REF/ALT + uniform chr:pos:REF:ALT IDs.
cat > "${MERGE_LIST}" << MERGEEOF
${NORM_WB}
${NORM_DC}
MERGEEOF

echo "Primary:   ${NORM_WGS}"
echo "Secondary: ${NORM_WB}, ${NORM_DC}"
echo "Mode:      union"

# ── Pass 1: attempt merge, catch strand flip conflicts ────────────────────────
echo "Pass 1: initial merge attempt..."
plink \
    --bfile "${NORM_WGS}" \
    --merge-list "${MERGE_LIST}" \
    --make-bed \
    --allow-no-sex \
    --allow-extra-chr \
    --out "${TMP}/pass1"

PASS1_EXIT=$?

# ── Pass 2: flip mismatching SNPs and retry if needed ────────────────────────
if [[ ${PASS1_EXIT} -ne 0 ]] && [[ -f "${TMP}/pass1-merge.missnp" ]]; then
    echo "Pass 1 found strand conflicts — flipping and retrying..."
    MISSNP="${TMP}/pass1-merge.missnp"
    NSNPS=$(wc -l < "${MISSNP}")
    echo "Flipping ${NSNPS} SNPs in ${NORM_WGS}..."

    plink \
        --bfile "${NORM_WGS}" \
        --flip "${MISSNP}" \
        --make-bed \
        --allow-no-sex \
        --allow-extra-chr \
        --out "${TMP}/wgs_harm_flipped"

    echo "Pass 2: re-attempting merge with flipped SNPs..."
    plink \
        --bfile "${TMP}/wgs_harm_flipped" \
        --merge-list "${MERGE_LIST}" \
        --make-bed \
        --allow-no-sex \
        --allow-extra-chr \
        --out "${OUT}"

    FINAL_EXIT=$?
else
    FINAL_EXIT=${PASS1_EXIT}
    if [[ ${FINAL_EXIT} -eq 0 ]]; then
        mv "${TMP}/pass1.bed" "${OUT}.bed"
        mv "${TMP}/pass1.bim" "${OUT}.bim"
        mv "${TMP}/pass1.fam" "${OUT}.fam"
        mv "${TMP}/pass1.log" "${OUT}.log"
    fi
fi

echo "=========================================="
echo "End time: $(date)"

if [[ ${FINAL_EXIT} -eq 0 ]]; then
    echo "Success."
    echo "Variants: $(wc -l < ${OUT}.bim)"
    echo "Samples:  $(wc -l < ${OUT}.fam)"
    du -sh "${OUT}.bed"
else
    echo "ERROR: merge failed with exit code ${FINAL_EXIT}" >&2
fi

exit ${FINAL_EXIT}