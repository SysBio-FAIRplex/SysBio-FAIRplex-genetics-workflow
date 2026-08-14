#!/bin/bash
#SBATCH --job-name=relatedness
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --partition=norm
#
# STEP 4 — CROSS-DATASET relatedness on the merged cohort, stratified by reused per-callset
# ancestry. Runs after step 3. Report-only: it decides nothing, step 5 adjudicates.
#
# Approach:
#   1. Build the common-variant set on the full merged cohort (--geno 0.05 on the union =
#      variants genotyped across all callsets -> same call-rate denominator for every sample,
#      so cross-dataset KING is fair).
#   2. Reuse the per-callset genotools ancestry labels (SAME GP2 panel across all 3, so directly
#      comparable). Union them, restrict to the common set. NO re-derivation on the merge.
#   3. Split the common set by ancestry (--keep per stratum) and run genotools --related
#      REPORT-ONLY within each stratum (no --ancestry: each subset is already homogeneous).
#      Step 5 (05_excludelist.py) makes the dup/relative picks — not genotools.
#
# Why split instead of genotools --ancestry: reuses the per-callset labels as locked instead
# of re-projecting on the merged intersection (fewer variants -> worse
# calls + label reconciliation). KING-robust at close-kinship thresholds tolerates the residual
# within-stratum structure.
# KNOWN LIMITATION: a related/dup pair split across two ancestry labels is not compared. Same-
# person cross-dataset dups (Rush<->ROSMAP) get the same label -> same stratum -> caught. Only
# label noise near admixed groups is at risk; acceptable (report-only + human adjudication).
#
#   ./submit.sh scripts/04_relatedness.sh

set -o pipefail

# Resolve the bundle root and load the one file that knows where anything lives.
BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

module load "${MOD_PLINK2}"
module load "${MOD_PYTHON}"

MERGED="${MERGED_DIR}/cohort_merged"
OUT_DIR="${MERGED_DIR}/relatedness"
COMMON="${OUT_DIR}/cohort_common"
SPLIT_DIR="${OUT_DIR}/by_ancestry"

# LBL_WGS / LBL_WB / LBL_DC (the per-callset genotools ancestry predictions, FID<TAB>IID<TAB>label,
# all projected against the same GP2 panel) come from config.sh.

MERGED_LABELS="${OUT_DIR}/merged_ancestry_labels.txt"   # FID IID label, deduped
COMMON_LABELS="${OUT_DIR}/common_ancestry_labels.txt"    # restricted to the common set

RELATED_CUTOFF=0.0884   # up to 2nd degree
DUP_CUTOFF=0.354        # duplicate / MZ

mkdir -p "${OUT_DIR}" "${SPLIT_DIR}"
cd "${WGS_ROOT}"
source "${VENV}/bin/activate"

echo "=========================================="
echo "Step 4 — cross-dataset relatedness, split by ancestry (report only)"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       ${SLURMD_NODENAME}"
echo "Merged in:  ${MERGED}"
echo "Start time: $(date)"
echo "=========================================="

for f in "${MERGED}.bed" "${LBL_WGS}" "${LBL_WB}" "${LBL_DC}"; do
    [[ -f "$f" ]] || { echo "ERROR: missing input: $f" >&2; exit 1; }
done

# ── Step 1: common-variant subset (fair cross-dataset denominator) ────────────
echo "Step 1: extracting common (cross-callset genotyped) autosomal biallelic SNPs..."
plink2 \
    --bfile "${MERGED}" \
    --autosome --snps-only --max-alleles 2 \
    --geno 0.05 \
    --make-bed --threads 32 \
    --out "${COMMON}"
[[ $? -eq 0 ]] || { echo "ERROR: common-variant subset failed" >&2; exit 1; }
NCOMMON=$(wc -l < "${COMMON}.bim")
NSAMP=$(wc -l < "${COMMON}.fam")
echo "Common set: ${NCOMMON} variants, ${NSAMP} samples"
[[ "${NCOMMON}" -gt 0 ]] || { echo "ERROR: 0 common variants — merge still misaligned?" >&2; exit 1; }

# ── Step 2: fair per-sample call rate (same denominator for all samples) ──────
echo "Step 2: per-sample missingness on the common set..."
plink2 --bfile "${COMMON}" --missing --threads 32 --out "${COMMON}"
[[ $? -eq 0 ]] || { echo "ERROR: --missing failed" >&2; exit 1; }
echo "Call-rate file: ${COMMON}.smiss"

# ── Step 3: reuse + union per-callset ancestry labels ─────────────────────────
# Dedup on FID+IID (the 94 same-IID Mayo dups appear in 2 callsets — identical label expected;
# a conflict is surfaced, first label kept). No re-derivation.
echo "Step 3: building merged ancestry-label table (reused per-callset labels)..."
awk -F'\t' 'FNR>1 {
        key=$1 SUBSEP $2
        if (key in lab) { if (lab[key]!=$3) print "WARN: conflicting label "$1" "$2": "lab[key]" vs "$3 > "/dev/stderr" }
        else lab[key]=$3
     } END { for (k in lab){ split(k,a,SUBSEP); print a[1]"\t"a[2]"\t"lab[k] } }' \
     "${LBL_WGS}" "${LBL_WB}" "${LBL_DC}" | sort > "${MERGED_LABELS}"
echo "  merged label rows: $(wc -l < "${MERGED_LABELS}")"

# Restrict to samples actually in the common set, keyed FID+IID (robust to FID differences).
awk 'NR==FNR { k[$1 SUBSEP $2]=1; next }
     ($1 SUBSEP $2) in k { print $1"\t"$2"\t"$3 }' \
     "${COMMON}.fam" "${MERGED_LABELS}" > "${COMMON_LABELS}"
NLAB=$(wc -l < "${COMMON_LABELS}")
echo "  common-set samples with a label: ${NLAB} / ${NSAMP}"
if [[ "${NLAB}" -lt "${NSAMP}" ]]; then
    echo "  NOTE: $(( NSAMP - NLAB )) common-set samples have NO reused label (genotools skip_fails?)"
    echo "        -> they are omitted from the per-ancestry KING; listed below:"
    awk 'NR==FNR { lab[$1 SUBSEP $2]=1; next }
         !(($1 SUBSEP $2) in lab) { print "        UNLABELED "$1" "$2 }' \
         "${COMMON_LABELS}" "${COMMON}.fam"
fi

echo "  per-ancestry counts (common set):"
awk '{print $3}' "${COMMON_LABELS}" | sort | uniq -c | sort -rn | sed 's/^/    /'

# ── Step 4: split by ancestry + genotools --related report-only per stratum ───
echo "Step 4: per-ancestry relatedness (report only)..."
LABELS=$(awk '{print $3}' "${COMMON_LABELS}" | sort -u)
FAILED=""
for L in ${LABELS}; do
    KEEP="${SPLIT_DIR}/keep_${L}.txt"
    awk -v l="${L}" '$3==l { print $1"\t"$2 }' "${COMMON_LABELS}" > "${KEEP}"
    N=$(wc -l < "${KEEP}")
    echo ""
    echo "  ---- ancestry ${L}: ${N} samples ----"
    if [[ "${N}" -lt 2 ]]; then
        echo "    <2 samples — no pairs possible, skipping."
        continue
    fi

    # No per-stratum variant filtering here: per-ancestry call-rate/MAF/HWE/LD are ALL deferred to
    # step 6, applied together on the deduplicated grain. KING runs on the global
    # common set — negligible effect at close-kinship thresholds, report-only + human-adjudicated.
    SUB="${SPLIT_DIR}/cohort_common_${L}"
    plink2 --bfile "${COMMON}" --keep "${KEEP}" --make-bed --threads 32 --out "${SUB}"
    [[ $? -eq 0 ]] || { echo "    ERROR: subset for ${L} failed" >&2; FAILED="${FAILED} ${L}(subset)"; continue; }

    # genotools --related on an already-homogeneous stratum (no --ancestry).
    genotools \
        --bfile "${SUB}" \
        --out "${SPLIT_DIR}/relatedness_${L}" \
        --related \
        --related_cutoff "${RELATED_CUTOFF}" \
        --duplicated_cutoff "${DUP_CUTOFF}" \
        --prune_related false \
        --prune_duplicated false \
        --full_output \
        --warn
    RC=$?
    if [[ ${RC} -ne 0 ]]; then
        echo "    ERROR: genotools --related failed for ${L} (exit ${RC})" >&2
        FAILED="${FAILED} ${L}(related:${RC})"
        continue
    fi
    echo "    genotools outputs for ${L}:"
    ls -1 "${SPLIT_DIR}/relatedness_${L}"* 2>/dev/null | sed 's/^/      /'
done

echo ""
echo "=========================================="
echo "End time: $(date)"
if [[ -n "${FAILED}" ]]; then
    echo "COMPLETED WITH FAILURES:${FAILED}" >&2
    echo "(If the FIRST/smallest stratum failed, genotools likely needs a flag tweak — check its" >&2
    echo " stderr; --related may not run standalone in this build.)" >&2
    exit 1
fi
echo "Success. Per-ancestry relatedness reports in: ${SPLIT_DIR}"
echo "Next (step 5): 05_excludelist.py consumes ${SPLIT_DIR}/relatedness_*<pairs> + ${COMMON}.smiss"
