#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1500G
#SBATCH --partition=largemem
#
# STEP 3 — plink1.9 union merge of the normalized bed filesets into one cohort.
#
# plink1.9 rather than plink2 because the non-concatenating --pmerge-list is unimplemented
# in this plink2 build. --allow-extra-chr covers the PAR codes plink1.9 does not recognize.
#
# To add a callset: add its NORM_* stem to SECONDARY. Nothing else changes.
#
#   ./submit.sh scripts/03_merge.sh

set -o pipefail          # NOT set -e — pass 1 is allowed to fail; that is the retry trigger.

BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"
: "${NORM_WGS:?config.sh failed to load}"

module load "${MOD_PLINK1}"

PRIMARY="${NORM_WGS}"
SECONDARY=("${NORM_WB}" "${NORM_DC}" "${NORM_BR}")

OUT="${MERGED_DIR}/cohort_merged"
TMP="${MERGED_DIR}/tmp_merge"
MERGE_LIST="${MERGED_DIR}/merge_list.txt"
mkdir -p "${MERGED_DIR}" "${TMP}"

# plink1.9 sizes threads and memory from the NODE, not the allocation. 90% of the request
# leaves headroom for everything outside plink's own workspace.
OPTS=(--make-bed --allow-no-sex --allow-extra-chr --threads "${SLURM_CPUS_PER_TASK:-4}")
[[ -n "${SLURM_MEM_PER_NODE:-}" ]] && OPTS+=(--memory $(( SLURM_MEM_PER_NODE * 9 / 10 )))

echo "cohort merge — job ${SLURM_JOB_ID} on ${SLURMD_NODENAME} — $(date)"
echo "  primary:   ${PRIMARY}"
printf '  secondary: %s\n' "${SECONDARY[@]}"
printf '%s\n' "${SECONDARY[@]}" > "${MERGE_LIST}"

plink --bfile "${PRIMARY}" --merge-list "${MERGE_LIST}" "${OPTS[@]}" --out "${TMP}/merged"
RC=$?

# Retry once with the conflicting SNPs strand-flipped in the primary. Since step 2 encodes
# REF/ALT in every variant ID, a shared ID cannot carry mismatched alleles — so this should
# no longer be reachable. Kept as a cheap net.
if [[ ${RC} -ne 0 && -f "${TMP}/merged-merge.missnp" ]]; then
    echo "pass 1: $(wc -l < "${TMP}/merged-merge.missnp") allele conflicts — flipping, retrying..."
    plink --bfile "${PRIMARY}" --flip "${TMP}/merged-merge.missnp" \
          "${OPTS[@]}" --out "${TMP}/primary_flipped"
    plink --bfile "${TMP}/primary_flipped" --merge-list "${MERGE_LIST}" \
          "${OPTS[@]}" --out "${TMP}/merged"
    RC=$?
fi

[[ ${RC} -eq 0 ]] || { echo "ERROR: merge failed (exit ${RC})" >&2; exit ${RC}; }

# Staged in TMP throughout, so a failed re-run never destroys a good cohort_merged.
for ext in bed bim fam log; do mv "${TMP}/merged.${ext}" "${OUT}.${ext}"; done

echo "done — $(date)"
echo "  variants: $(wc -l < "${OUT}.bim")"
echo "  samples:  $(wc -l < "${OUT}.fam")"
