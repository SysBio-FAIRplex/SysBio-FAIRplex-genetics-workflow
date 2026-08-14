#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96g
#SBATCH --partition=norm
#
# STEP 2b — GENOME-WIDE variant normalization. The pass validated on chr22 by step 2a
# (3-way exact-ID 84k -> 330k). Run ONCE PER CALLSET; the three jobs are independent.
#
# Same pass as the chr22 check, minus --chr:
#   --fa <GRCh38 ref> --ref-from-fa force   -> REF/ALT set from the reference
#   --output-chr chrM                       -> chr-prefixed chromosome codes
#   --set-all-var-ids '@:#:$r:$a'           -> uniform IDs chr:pos:REF:ALT across all callsets
#   --new-id-max-allele-len 1000 missing    -> survive long WGS indels (fallback ID for >1000bp)
# Output is BED so the plink1.9 union merge in step 3 consumes it directly.
#
# Parameterized via --export (see launch lines in the handoff / below). Required vars:
#   PFILE = input --pfile stem (filtered pgen)
#   OUT   = output bed stem (normalized)
#   TAG   = short label for logging
#
# The three stems are already in config.sh, so the launch lines stay short (README §3 step 2):
#   ./submit.sh scripts/02b_normalize.sh --job-name=norm_wgs --export=PFILE=$PF_WGS,OUT=$NORM_WGS,TAG=wgs
#   ./submit.sh scripts/02b_normalize.sh --job-name=norm_wb  --export=PFILE=$PF_WB,OUT=$NORM_WB,TAG=wb
#   ./submit.sh scripts/02b_normalize.sh --job-name=norm_dc  --export=PFILE=$PF_DC,OUT=$NORM_DC,TAG=dc

set -o pipefail

# Resolve the bundle root and load the one file that knows where anything lives.
BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

module load "${MOD_PLINK2}"

REF="${REF_FASTA}"

: "${PFILE:?set PFILE via --export}"
: "${OUT:?set OUT via --export}"
: "${TAG:?set TAG via --export}"

echo "=========================================="
echo "genome-wide normalization — ${TAG}"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       ${SLURMD_NODENAME}"
echo "Start time: $(date)"
echo "Input:      ${PFILE}"
echo "Output:     ${OUT}"
echo "Reference:  ${REF}"
echo "=========================================="

if [[ ! -f "${REF}" ]]; then
    echo "ERROR: reference FASTA not found at ${REF}" >&2
    exit 1
fi

mkdir -p "$(dirname "${OUT}")"

plink2 \
    --pfile "${PFILE}" \
    --fa "${REF}" --ref-from-fa force \
    --output-chr chrM \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 1000 missing \
    --make-bed \
    --out "${OUT}"
RC=$?

echo "=========================================="
echo "End time: $(date)"
if [[ ${RC} -ne 0 ]]; then
    echo "ERROR: normalization failed for ${TAG} (exit ${RC})" >&2
    exit ${RC}
fi
echo "Success (${TAG})."
echo "Variants: $(wc -l < "${OUT}.bim")"
echo "Samples:  $(wc -l < "${OUT}.fam")"
du -sh "${OUT}.bed"
echo "First 3 IDs:"; head -3 "${OUT}.bim" | awk '{print "  "$2}'
