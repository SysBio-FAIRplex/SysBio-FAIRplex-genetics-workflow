#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96g
#SBATCH --partition=norm
#
# STEP 2 — genome-wide variant normalization. One job per callset; they are independent.
#
#   --fa <GRCh38> --ref-from-fa force     REF/ALT set from the reference (most cross-callset
#                                         disagreement is REF/ALT swaps)
#   --output-chr chrM                     chr-prefixed chromosome codes
#   --set-all-var-ids '@:#:$r:$a'         uniform chr:pos:REF:ALT IDs across all callsets
#   --new-id-max-allele-len 1000 missing  survive long indels (plink2's default aborts over 23 bp).
#                                         Keep identical across callsets or the IDs stop matching.
#
# Output is BED so the plink1.9 union merge in step 3 consumes it directly.
#
#   ./submit.sh scripts/02_normalize.sh --job-name=norm_<tag> \
#       --export=PFILE=$PF_<X>,OUT=$NORM_<X>,TAG=<tag>

set -eo pipefail

BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"
module load "${MOD_PLINK2}"

: "${PFILE:?set PFILE via --export}"
: "${OUT:?set OUT via --export}"
: "${TAG:?set TAG via --export}"

[[ -f "${REF_FASTA}" ]] || { echo "ERROR: no reference FASTA at ${REF_FASTA}" >&2; exit 1; }

echo "normalize ${TAG} — job ${SLURM_JOB_ID} on ${SLURMD_NODENAME} — $(date)"
echo "  in:  ${PFILE}"
echo "  out: ${OUT}"

mkdir -p "$(dirname "${OUT}")"

plink2 \
    --pfile "${PFILE}" \
    --fa "${REF_FASTA}" --ref-from-fa force \
    --output-chr chrM \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 1000 missing \
    --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --make-bed \
    --out "${OUT}"

echo "done ${TAG} — $(date)"
echo "  variants: $(wc -l < "${OUT}.bim")"
echo "  samples:  $(wc -l < "${OUT}.fam")"
echo "  first 3 IDs:"; head -3 "${OUT}.bim" | awk '{print "    "$2}'
