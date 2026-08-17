#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96g
#SBATCH --partition=norm
#
# STEP 2 — genome-wide variant normalization. One job per callset; they are independent.
#
#   --fa <GRCh38> --ref-from-fa force     REF/ALT set from the reference (~75% of the
#                                         cross-callset disagreement was REF/ALT swaps)
#   --output-chr chrM                     chr-prefixed chromosome codes
#   --set-all-var-ids '@:#:$r:$a'         uniform chr:pos:REF:ALT IDs across all callsets
#   --new-id-max-allele-len 1000 missing  survive long WGS indels. plink2 defaults to
#                                         `23 error`, which ABORTS on any allele over 23 bp.
#                                         Keep this identical across callsets or the IDs stop
#                                         being comparable.
#
# Output is BED so the plink1.9 union merge in step 3 consumes it directly.
#
# This pass was validated on chr22 before being run genome-wide: 3-way exact-ID overlap went
# from ~83,905 to ~331,015, which is the allele-concordant ceiling. The gate script that
# measured it (02a_normalize_check.sh) was retired once it passed — see PROJECT_LOG.md
# 2026-08-14 for the numbers and the reasoning.
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
