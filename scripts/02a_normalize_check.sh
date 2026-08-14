#!/bin/bash
#SBATCH --job-name=normalize_check
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128g
#SBATCH --partition=norm
#
# STEP 2a — chr22 NORMALIZATION CHECK. This is the gate for step 2b: it validates the
# variant-ID realignment on chr22 BEFORE committing to the genome-wide re-normalize + re-merge.
#
# For each of the 3 filtered pgens, extract chr22 and run the planned normalization pass:
#   --fa <GRCh38 ref> --ref-from-fa force   -> REF/ALT set from the reference (fixes ~75% ref/alt swaps)
#   --output-chr chrM                       -> chromosome codes come out chr-prefixed (chr22)
#   --set-all-var-ids '@:#:$r:$a'           -> uniform IDs chr22:pos:REF:ALT across ALL callsets
# Then recompute the 3-way EXACT variant-ID overlap on chr22.
#
# EXPECT (success criteria):
#   - IDs come out "chr22:..." for all three (DivCo prefix no longer special-cased)
#   - 3-way exact-ID overlap jumps from ~84k -> ~331k (~ the allele-concordant position count)
#   - few / no "no matching reference allele" variants in the plink2 logs
# If met -> run the identical pass genome-wide on all 3 (step 2b), then merge (step 3).
#
# NOTE (carries into the genome-wide run): WGS callsets contain long indels. plink2's
# --set-all-var-ids defaults to `--new-id-max-allele-len 23 error`, which ABORTS on any allele
# >23 bp. We set `--new-id-max-allele-len 1000 missing` so the run survives; the rare variant
# whose allele exceeds 1000 bp gets a fallback ID (chr22:pos) and simply won't align. Keep this
# flag identical in the genome-wide pass.

set -o pipefail

# Resolve the bundle root and load the one file that knows where anything lives.
BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

module load "${MOD_PLINK2}"

# REF_FASTA and the three PF_* filtered-pgen stems come from config.sh.
REF="${REF_FASTA}"

OUT_DIR="${MERGED_DIR}/norm_check"
mkdir -p "${OUT_DIR}"

echo "=========================================="
echo "Step 2a — chr22 normalization check (gate for step 2b)"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       ${SLURMD_NODENAME}"
echo "Start time: $(date)"
echo "Reference:  ${REF}"
echo "=========================================="

if [[ ! -f "${REF}" ]]; then
    echo "ERROR: reference FASTA not found at ${REF}" >&2
    echo "Place GRCh38_full_analysis_set_plus_decoy_hla.fa.zst under data/ref/ first." >&2
    exit 1
fi

# ── Normalize chr22 of one callset ────────────────────────────────────────────
#   $1 = --pfile stem   $2 = short tag (out name)
normalize () {
    local pf="$1" tag="$2"
    echo ""
    echo "---- normalizing ${tag} (chr22) ----"
    plink2 \
        --pfile "${pf}" \
        --chr 22 \
        --fa "${REF}" --ref-from-fa force \
        --output-chr chrM \
        --set-all-var-ids '@:#:$r:$a' \
        --new-id-max-allele-len 1000 missing \
        --make-bed \
        --out "${OUT_DIR}/${tag}_chr22_norm"
    local rc=$?
    if [[ ${rc} -ne 0 ]]; then
        echo "ERROR: normalization failed for ${tag} (exit ${rc})" >&2
        exit ${rc}
    fi
    # Surface ref-from-fa outcome from the log (mismatch / no-match reporting varies by build).
    echo "  ref-from-fa notes for ${tag}:"
    grep -iE "reference allele|ref-from-fa|no matching|mismatch|--fa" \
        "${OUT_DIR}/${tag}_chr22_norm.log" | sed 's/^/    /' || true
    echo "  variants (chr22, ${tag}): $(wc -l < "${OUT_DIR}/${tag}_chr22_norm.bim")"
    echo "  first 3 IDs:"
    head -3 "${OUT_DIR}/${tag}_chr22_norm.bim" | awk '{print "    "$2}'
}

normalize "${PF_WGS}" wgs
normalize "${PF_WB}"  wb
normalize "${PF_DC}"  dc

# ── 3-way EXACT variant-ID overlap after normalization ─────────────────────────
# IDs are now uniform (chr22:pos:REF:ALT) across all three — no chr-strip special-casing.
echo ""
echo "=== exact variant-ID overlap on chr22 AFTER normalization ==="
awk '{print $2}' "${OUT_DIR}/wgs_chr22_norm.bim" | sort -u > "${OUT_DIR}/wgs.id"
awk '{print $2}' "${OUT_DIR}/wb_chr22_norm.bim"  | sort -u > "${OUT_DIR}/wb.id"
awk '{print $2}' "${OUT_DIR}/dc_chr22_norm.bim"  | sort -u > "${OUT_DIR}/dc.id"

printf "WGS ids:          %s\n" "$(wc -l < "${OUT_DIR}/wgs.id")"
printf "WB ids:           %s\n" "$(wc -l < "${OUT_DIR}/wb.id")"
printf "DivCo ids:        %s\n" "$(wc -l < "${OUT_DIR}/dc.id")"
printf "WGS n WB:         %s\n" "$(comm -12 "${OUT_DIR}/wgs.id" "${OUT_DIR}/wb.id" | wc -l)"
printf "3-way exact-ID:   %s\n" "$(comm -12 "${OUT_DIR}/wgs.id" "${OUT_DIR}/wb.id" | comm -12 - "${OUT_DIR}/dc.id" | wc -l)"

echo ""
echo "INTERPRET: baseline 3-way exact-ID (pre-norm, chr-stripped) was ~83,905;"
echo "allele-concordant positions were ~331,015. A jump toward ~331k = fix works ->"
echo "proceed to step 2b (02b_normalize.sh) on all 3 callsets, then step 3 (03_merge.sh)."

echo "=========================================="
echo "End time: $(date)"
