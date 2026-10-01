#!/bin/bash
#SBATCH --job-name=ancestry_qc
#SBATCH --time=10:00:00
# ~13-14 min per stratum x 11, dominated by scanning cohort_merged (~527 GiB .bed); 10 h also
# covers the AF build and the second prune+PCA.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --partition=norm
#
# STEP 6 — per-ancestry variant QC + within-ancestry PCA, in one submission.
# Rationale: METHODS.md §5-6.
#
# STAGES
#   A  per ANC: extract + variant QC from cohort_merged  -> unfiltered/cohort_<ANC>_qc
#      per ANC: LD-prune + PCA on that                   -> unfiltered/cohort_<ANC>_pca  (baseline)
#   B  build the AF-concordance exclusion list from the UNFILTERED filesets, never filtered input
#   C  per ANC: --exclude the list                       -> cohort_<ANC>_qc   (ASSOCIATION set)
#   D  per ANC: LD-prune + PCA on the filtered set       -> cohort_<ANC>_pca  (COVARIATES)
#   E  both retained_samples_manifest.csv files (review/plot_af_filter_effect.py reads the pair)
#
# INPUT: cohort_merged + step 5's retained_manifest.csv (only retained samples are kept, so no
# --remove) + clinical_core's sample_annot.csv.
#
# Long-range-LD regions are excluded at the PRUNE step only: they would otherwise dominate a top
# PC, but the MHC and 17q21.31 (MAPT) are real loci and must stay in the association set.
# The BED ships at ref/highld_exclude_hg38.bed (config.sh HIGHLD_BED); absent = warned, not fatal.
#
# Autosomes only (chrX needs --merge-par + sex-aware handling). All 11 strata run; the
# >=100-per-arm viability cut is applied later, at the phenotype boundary.
#
# KNOBS
#   AF_EXCLUDE=none       run the whole step unfiltered on purpose (stages B-D skip the exclusion).
#   SKIP_AF_BUILD=1       reuse the exclusion list on disk instead of rebuilding it.
#   FORCE_QC=1            rebuild stage A even when a valid unfiltered fileset is present.
#   HWE_REQUIRE_EXCESS=0  ungated HWE union (reproduces the 4,415-variant list).
#   plus every af_concordance_build knob (THRESH, ZMIN, MIN_CELL, HWE, MISHAP, ...) — passed through.
#
#   ./submit.sh scripts/06_ancestry_qc.sh
set -o pipefail

# Resolve the bundle root and load the one file that knows where anything lives.
BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

MERGED=${MERGED_DIR}/cohort_merged
MANIFEST=${MERGED_DIR}/relatedness/retained_manifest.csv
OUT_DIR=${MERGED_DIR}/by_ancestry_qc
UNF_DIR=${OUT_DIR}/unfiltered
KEEP_DIR=${OUT_DIR}/keep
LOG=${OUT_DIR}/step6_summary.txt
mkdir -p "$OUT_DIR" "$UNF_DIR" "$KEEP_DIR"

# ── QC thresholds (locked) ──
GENO=0.05
MAF=0.01
HWE_QC=1e-6
NPC=10
# LD pruning for the PCA input. The window is PHYSICAL (kb): at WGS density a 200-variant window
# spans only ~60 kb, far shorter than LD extends.
LD="1000kb 1 0.1"
THREADS=32

# Long-range-LD / inversion regions (hg38, 0-based BED, chr-prefixed). PCA input only.
EXCL=${EXCL:-${HIGHLD_BED}}
if [[ -f "$EXCL" ]]; then
    EXCL_ARG="--exclude bed0 $EXCL"
    EXCL_NOTE="$(wc -l < "$EXCL") regions"
else
    EXCL_ARG=""
    EXCL_NOTE="NONE FOUND at $EXCL — PCs may be driven by MHC/inversion structure"
fi

# Built at stage B, applied at stage C to the ASSOCIATION set as well as the PCA input: a variant
# mismapped enough to bend PC1 also gives a spurious association no PC adjustment reaches.
AF_EXCLUDE=${AF_EXCLUDE:-${MERGED_DIR}/exclude_af_concordance.txt}
AF_WORK=${MERGED_DIR}/af_concordance
DISC=${DISC:-${MERGED_DIR}/concordance/per_variant_discordance.tsv}

# Stage B rebuilds the list inside this job, so a stale list on disk cannot be applied. Only the
# SKIP_AF_BUILD path reuses one, and it keeps a date guard.

module load "${MOD_PLINK2}"

echo "==========================================" | tee "$LOG"
echo "Step 6 — per-ancestry variant QC + PCA (single pass)" | tee -a "$LOG"
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURMD_NODENAME}  Start: $(date)" | tee -a "$LOG"
echo "Merged:   ${MERGED}"                         | tee -a "$LOG"
echo "Manifest: ${MANIFEST}"                       | tee -a "$LOG"
echo "Annot:    ${ANNOT}"                          | tee -a "$LOG"
echo "geno<=${GENO} maf>=${MAF} hwe>=${HWE_QC} pcs=${NPC} ld='${LD}'" | tee -a "$LOG"
echo "PCA-input long-range-LD exclusion: ${EXCL_NOTE}" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"

if [[ ! -f "${MANIFEST}" ]]; then echo "ERROR: manifest not found: ${MANIFEST}" >&2; exit 1; fi
if [[ ! -f "${MERGED}.bed" ]]; then echo "ERROR: merged bed not found: ${MERGED}.bed" >&2; exit 1; fi

# ancestries present in the retained manifest (col 3; skip header + blanks)
ANCS=$(awk -F, 'NR>1 && $3!="" {print $3}' "$MANIFEST" | sort -u)
echo "Ancestries: $(echo $ANCS | tr '\n' ' ')" | tee -a "$LOG"

# ─────────────────────────────────────────────────────────────────────────────
# prune_and_pca <input-stem> <output-dir> <ANC> <n-samples>  -> echoes "<prune_in> <pcs>"
# Shared by stages A and D so the only difference between the two generations is the exclusion list.
# ─────────────────────────────────────────────────────────────────────────────
prune_and_pca() {
    local stem=$1 dir=$2 anc=$3 n=$4
    local prune=${dir}/${anc}_prune
    local pca=${dir}/cohort_${anc}_pca

    if ! plink2 --bfile "$stem" $EXCL_ARG --indep-pairwise $LD \
                --threads "$THREADS" --out "$prune" >/dev/null 2>&1; then
        echo "PRUNE_FAIL -"
        return
    fi
    local nin
    nin=$(wc -l < "${prune}.prune.in" 2>/dev/null || echo 0)

    if [[ "$n" -ge 3 && "$nin" -ge 2 ]]; then
        local npc_eff=$NPC
        [[ "$n" -le "$NPC" ]] && npc_eff=$((n-1))
        if plink2 --bfile "$stem" --extract "${prune}.prune.in" \
                  --pca "$npc_eff" --threads "$THREADS" --out "$pca" >/dev/null 2>&1; then
            echo "$nin $npc_eff"
        else
            echo "$nin PCA_FAIL"
        fi
    else
        echo "$nin too_small"
    fi
}

# ═════════════════════════════════════════════════════════════════════════════
# STAGE A — per-ancestry QC + PCA, UNFILTERED. The baseline, kept permanently.
# ═════════════════════════════════════════════════════════════════════════════
echo ""                                                        | tee -a "$LOG"
echo "--- STAGE A: unfiltered QC + PCA -> ${UNF_DIR} ---"       | tee -a "$LOG"
printf "%-5s %8s %12s %10s %5s  %s\n" "ANC" "samples" "qc_variants" "prune_in" "pcs" "note" | tee -a "$LOG"

declare -A ANC_N
for ANC in $ANCS; do
    KEEP=${KEEP_DIR}/keep_${ANC}.txt
    awk -F, -v a="$ANC" 'NR>1 && $3==a {print $1"\t"$2}' "$MANIFEST" > "$KEEP"
    N=$(wc -l < "$KEEP")
    ANC_N[$ANC]=$N
    if [[ "$N" -lt 2 ]]; then
        printf "%-5s %8s %12s %10s %5s  %s\n" "$ANC" "$N" "SKIP(<2)" "-" "-" "" | tee -a "$LOG"
        continue
    fi

    UQC=${UNF_DIR}/cohort_${ANC}_qc

    # Reuse a fileset that postdates the merge: stage A's inputs are fixed, so rebuilding it would
    # reproduce the same bytes at the cost of another full scan.
    NOTE=""
    if [[ -z "${FORCE_QC}" && -f "${UQC}.bed" && "${UQC}.bed" -nt "${MERGED}.bed" ]]; then
        NOTE="reused (postdates merge; FORCE_QC=1 to rebuild)"
    else
        plink2 --bfile "$MERGED" --keep "$KEEP" \
               --autosome \
               --geno "$GENO" --maf "$MAF" --hwe "$HWE_QC" keep-fewhet \
               --make-bed --threads "$THREADS" --out "$UQC" \
            || { printf "%-5s %8s %12s %10s %5s  %s\n" "$ANC" "$N" "QC_FAILED" "-" "-" "" | tee -a "$LOG"; continue; }
        NOTE="built"
    fi
    NVAR=$(wc -l < "${UQC}.bim")

    if [[ -z "${FORCE_QC}" && -f "${UNF_DIR}/cohort_${ANC}_pca.eigenvec" \
          && "${UNF_DIR}/cohort_${ANC}_pca.eigenvec" -nt "${UQC}.bed" ]]; then
        NIN=$(wc -l < "${UNF_DIR}/${ANC}_prune.prune.in" 2>/dev/null || echo "-")
        PCS=$(head -1 "${UNF_DIR}/cohort_${ANC}_pca.eigenvec" | awk '{print NF-2}')
        NOTE="${NOTE}; PCA reused"
    else
        read -r NIN PCS <<< "$(prune_and_pca "$UQC" "$UNF_DIR" "$ANC" "$N")"
    fi

    printf "%-5s %8s %12s %10s %5s  %s\n" "$ANC" "$N" "$NVAR" "$NIN" "$PCS" "$NOTE" | tee -a "$LOG"
done

# ═════════════════════════════════════════════════════════════════════════════
# STAGE B — build the AF-concordance exclusion list from the UNFILTERED filesets.
# ═════════════════════════════════════════════════════════════════════════════
echo ""                                                              | tee -a "$LOG"
echo "--- STAGE B: AF-concordance exclusion list ---"                | tee -a "$LOG"

if [[ "$AF_EXCLUDE" == "none" ]]; then
    echo "AF_EXCLUDE=none — running deliberately UNFILTERED. Stages C/D copy stage A's" | tee -a "$LOG"
    echo "variant set through unchanged; the before/after comparison will be a null one." | tee -a "$LOG"
    AFX_ARG=""
elif [[ -n "${SKIP_AF_BUILD}" && -f "$AF_EXCLUDE" ]]; then
    # The only path by which a list this job did not build reaches stage C, so it keeps a date guard.
    if [[ "$AF_EXCLUDE" -ot "${MERGED}.bed" ]]; then
        echo "REFUSING: SKIP_AF_BUILD was set but ${AF_EXCLUDE} predates ${MERGED}.bed," >&2
        echo "  so it was built from a different cohort. Unset SKIP_AF_BUILD to rebuild it." >&2
        exit 1
    fi
    echo "SKIP_AF_BUILD — reusing $(wc -l < "$AF_EXCLUDE") variants from ${AF_EXCLUDE}" | tee -a "$LOG"
    AFX_ARG="--exclude $AF_EXCLUDE"
else
    [[ -f "$ANNOT" ]] || { echo "ERROR: sample annot not found: ${ANNOT}" >&2
                           echo "  Run clinical_core.py (§12a writes it; no PCs needed)." >&2; exit 1; }
    module load "${MOD_PYTHON}"
    source "${VENV}/bin/activate"

    # Stage B needs both plink generations (plink1.9 for --assoc, plink2 for --hwe/--freq), and they
    # are one module family: resolve both to absolute paths, then reload MOD_PLINK2, because stages
    # C and D call bare `plink2`.
    PLINK2_BIN=$(command -v plink2)
    module load "${MOD_PLINK1}"; PLINK1_BIN=$(command -v plink)
    module load "${MOD_PLINK2}"          # restore plink2 for stages C/D — do not remove
    [[ -x "$PLINK1_BIN" ]] || { echo "ERROR: plink not found via ${MOD_PLINK1}" >&2; exit 1; }
    [[ -x "$PLINK2_BIN" ]] || { echo "ERROR: plink2 not found via ${MOD_PLINK2}" >&2; exit 1; }

    python3 "${BUNDLE}/scripts/af_concordance_build.py" \
        --qc-dir "${UNF_DIR}" \
        --annot "${ANNOT}" \
        --out "${AF_EXCLUDE}" \
        --work "${AF_WORK}" \
        --plink1 "${PLINK1_BIN}" \
        --plink2 "${PLINK2_BIN}" \
        --thresh "${THRESH:-0.05}" \
        --zmin "${ZMIN:-5.0}" \
        --min-cell "${MIN_CELL:-100}" \
        --hwe "${HWE:-1e-4}" \
        --min-hwe-controls "${MIN_HWE_CONTROLS:-50}" \
        ${HWE_BOTH_TAILS:+--hwe-both-tails} \
        $([[ "${HWE_REQUIRE_EXCESS:-1}" == "0" ]] && echo --no-hwe-require-excess) \
        --mishap "${MISHAP:-0}" \
        --min-mishap "${MIN_MISHAP:-100}" \
        --discordance "${DISC}" \
        --disc-rate "${DISC_RATE:-0.50}" \
        ${ANCS_AF:+--ancs ${ANCS_AF}} \
        ${DXS:+--dx ${DXS}} 2>&1 | tee -a "$LOG"
    RC=${PIPESTATUS[0]}
    if [[ "$RC" -ne 0 ]]; then
        echo "ERROR: af_concordance_build failed (exit ${RC}). Refusing to continue —" >&2
        echo "  proceeding would silently produce an UNFILTERED association set under the" >&2
        echo "  filtered names, which is the failure mode this step exists to prevent." >&2
        echo "  To run unfiltered on purpose: AF_EXCLUDE=none ./submit.sh scripts/06_ancestry_qc.sh" >&2
        exit "$RC"
    fi
    echo "built $(wc -l < "$AF_EXCLUDE") variants -> ${AF_EXCLUDE}" | tee -a "$LOG"
    AFX_ARG="--exclude $AF_EXCLUDE"
fi

# ═════════════════════════════════════════════════════════════════════════════
# STAGE C+D — apply the exclusion, then prune + PCA on the filtered set.
# Stage C reads a stratum-sized fileset, not cohort_merged: --geno/--maf/--hwe are per-variant
# on a fixed sample set, so they commute with --exclude and there is nothing to recompute.
# ═════════════════════════════════════════════════════════════════════════════
echo ""                                                        | tee -a "$LOG"
echo "--- STAGE C+D: filtered association set + PCA -> ${OUT_DIR} ---" | tee -a "$LOG"
printf "%-5s %8s %12s %12s %10s %5s\n" "ANC" "samples" "unfiltered" "qc_variants" "prune_in" "pcs" | tee -a "$LOG"

for ANC in $ANCS; do
    N=${ANC_N[$ANC]}
    UQC=${UNF_DIR}/cohort_${ANC}_qc
    [[ "$N" -ge 2 && -f "${UQC}.bed" ]] || continue

    NVAR_U=$(wc -l < "${UQC}.bim")
    QC=${OUT_DIR}/cohort_${ANC}_qc

    plink2 --bfile "$UQC" $AFX_ARG --make-bed --threads "$THREADS" --out "$QC" >/dev/null 2>&1 \
        || { printf "%-5s %8s %12s %12s %10s %5s\n" "$ANC" "$N" "$NVAR_U" "EXCLUDE_FAILED" "-" "-" | tee -a "$LOG"; continue; }
    NVAR=$(wc -l < "${QC}.bim")

    read -r NIN PCS <<< "$(prune_and_pca "$QC" "$OUT_DIR" "$ANC" "$N")"
    printf "%-5s %8s %12s %12s %10s %5s\n" "$ANC" "$N" "$NVAR_U" "$NVAR" "$NIN" "$PCS" | tee -a "$LOG"
done

# ═════════════════════════════════════════════════════════════════════════════
# STAGE E — assemble BOTH manifests, so the filter's effect is measured rather than asserted;
# review/plot_af_filter_effect.py reads the pair.
# ═════════════════════════════════════════════════════════════════════════════
echo ""                                                        | tee -a "$LOG"
echo "--- STAGE E: sample manifests (unfiltered + filtered) ---" | tee -a "$LOG"

module load "${MOD_PYTHON}"
source "${VENV}/bin/activate"

for GEN in unfiltered filtered; do
    if [[ "$GEN" == "unfiltered" ]]; then PCA_DIR=${UNF_DIR}; else PCA_DIR=${OUT_DIR}; fi
    echo "[${GEN}]" | tee -a "$LOG"
    python3 "${BUNDLE}/scripts/ancestry_qc_manifest.py" \
        --manifest "${MANIFEST}" \
        --pca-dir  "${PCA_DIR}" \
        --out      "${PCA_DIR}/retained_samples_manifest.csv" \
        --npc "${NPC}" 2>&1 | tee -a "$LOG"
done

echo "==========================================" | tee -a "$LOG"
echo "Step 6 complete: $(date)"                   | tee -a "$LOG"
echo "  ASSOCIATION set : ${OUT_DIR}/cohort_<ANC>_qc.{bed,bim,fam}   (read by step 7)"     | tee -a "$LOG"
echo "  COVARIATES      : ${OUT_DIR}/cohort_<ANC>_pca.eigenvec       (read by §12)"        | tee -a "$LOG"
echo "  BASELINE        : ${UNF_DIR}/  — unfiltered twin of both, for the before/after"    | tee -a "$LOG"
echo "  EXCLUSION LIST  : ${AF_EXCLUDE}"                                                   | tee -a "$LOG"
echo ""                                                                                    | tee -a "$LOG"
echo "NEXT: python3 analysis_grain.py  (PCs changed -> §12 rebuilds analysis_grain.csv)"   | tee -a "$LOG"
echo "      then ./submit.sh scripts/07_gwas.sh"                                           | tee -a "$LOG"
