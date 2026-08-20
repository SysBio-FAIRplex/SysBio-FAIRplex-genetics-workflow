#!/bin/bash
#SBATCH --job-name=ancestry_qc
#SBATCH --time=10:00:00
# MEASURED, not guessed: the 2026-07-24 run took ~2.5 h wall for all 11 strata — a flat ~13-14 min
# each, because cost is dominated by scanning cohort_merged (171M x 13,237 = ~527 GiB .bed), which
# every stratum pays regardless of how many samples it keeps. EUR (9.8k samples) cost barely more
# than CAH (71).
#
# 10 h rather than the old 6 h because this script now does in ONE job what used to be three
# submissions (step 6, af_concordance_build, step 6 again). The extra work is NOT another merged
# scan — stage A pays that once and stages C/D read stratum-sized filesets — so the increment is
# the AF build plus a second prune+PCA, not a second pass.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --partition=norm
#
# STEP 6 — per-ancestry variant QC + within-ancestry PCA, in a SINGLE pass.
# Produces the GWAS-ready grain: one QC'd genotype set + PCs per ancestry stratum.
#
# ─────────────────────────────────────────────────────────────────────────────
# WHY THIS IS ONE PASS NOW (it used to be two, with a hand-run script between them)
#
# The old shape was: run step 6 unfiltered -> run af_concordance_build against its output ->
# run step 6 AGAIN to apply the resulting exclusion list. That was justified by a claimed
# circular ordering: the AF build needs the analysis grain, the grain needs §12, §12 needs the
# PCs, and the PCs are step 6's output. Two facts dissolve it.
#
#   1. THE AF BUILD NEVER NEEDED THE GRAIN. It reads exactly IID -> (source_callset,
#      dx_detailed) and builds its comparison cells as (dx x callset); the stratum comes from
#      which fileset a sample is in, not from any column. Both fields are pure clinical output,
#      available from §7's crosswalk and §4's reconciliation before step 1 runs. They only ever
#      arrived via the grain because the grain is the file that happens to carry dx and the PCs
#      in one CSV. §12a now writes them separately as sample_annot.csv, and the cycle is gone.
#
#   2. THE QC PASS NEVER NEEDED RE-RUNNING. --geno, --maf and --hwe are per-variant statistics
#      computed on a fixed sample set (there is no --mind here), so they COMMUTE with --exclude:
#      QC-then-exclude and exclude-then-QC give the identical variant set. Pass 2's expensive
#      re-scan of cohort_merged was recomputing numbers it already had. Stage C therefore applies
#      the exclusion to stage A's output instead, reading a stratum-sized fileset rather than
#      527 GiB. This is the same identity the 2026-08-19 premise test relied on.
#
# So: QC once, build the list from that, apply it, prune and PCA on both. What used to be a
# two-submission sequence with a mandatory `mv` between the passes is one submission with no
# manual step and no one-way door.
#
# ─────────────────────────────────────────────────────────────────────────────
# STAGES
#
#   A  per ANC: extract + variant QC from cohort_merged   -> unfiltered/cohort_<ANC>_qc
#      per ANC: LD-prune + PCA on that                    -> unfiltered/cohort_<ANC>_pca
#              This is the BASELINE, and it is a permanent named output. It is not something
#              you have to remember to preserve before the filtered run overwrites it — that
#              requirement is exactly what made the old step 6 a one-way door.
#
#   B  build the AF-concordance exclusion list FROM THE UNFILTERED filesets.
#              Building it from unfiltered input is not incidental. The old ordering derived
#              each list from a directory that already had the previous list applied, so the
#              filter's input was pre-filtered by its own output and looked clean either way.
#              Here the input generation is fixed by construction.
#
#   C  per ANC: --exclude the list                        -> cohort_<ANC>_qc   (the ASSOCIATION set)
#   D  per ANC: LD-prune + PCA on the filtered set        -> cohort_<ANC>_pca  (the COVARIATES)
#   E  assemble both retained_samples_manifest.csv files (unfiltered + filtered) so the
#      before/after eta^2 comparison in review/plot_af_filter_effect.py has both generations.
#
# INPUT  : cohort_merged + step 5's retained_manifest.csv + clinical_core's sample_annot.csv
#          (the excludelist is applied simply by keeping only retained samples — dropped
#           dups/relatives/QC-fails are absent from the manifest, so no --remove is needed).
# PER ANC: 1. keep-list from the manifest (FID IID where ancestry==ANC)
#          2. extract + VARIANT QC in one plink2 pass (autosomes):
#               --geno 0.05  (variant call rate >=95%, computed WITHIN the stratum -> naturally
#                             drops variants not genotyped across the callsets present in ANC)
#               --maf 0.01   (common-variant GWAS; a rare/burden set would be a separate pass)
#               --hwe 1e-6 keep-fewhet  (remove excess-het genotyping artifacts; keep het-deficient
#                             real signal). HWE is ancestry-specific -> must run within stratum.
#          3. LD-prune (--indep-pairwise 1000kb 1 0.1) EXCLUDING long-range-LD/inversion
#             regions -> PCA input ONLY. The kb window is deliberate: a variant-count window
#             covers far too little physical distance at WGS marker density.
#          4. PCA on the pruned set (--pca 10; capped for tiny strata) -> eigenvec/eigenval
#     NOTE: the ASSOCIATION set is the full QC-passing variants; the LD-pruned set is used ONLY
#           for PCA (standard practice).
#     LONG-RANGE LD: pruning alone does not neutralize the MHC or the big inversions — enough
#           correlated structure survives --indep-pairwise for them to dominate a top PC, which
#           then encodes inversion/HLA haplotype instead of ancestry and propagates into every
#           GWAS as a covariate. The exclusion is applied at the PRUNE step only, so it shapes
#           the PCA input and never the association set. This matters concretely: 17q21.31 is
#           MAPT (top PSP locus, major PD locus) and the MHC is a real AD locus — masking them
#           from association would delete the signals we most expect to see.
#           The BED ships in this bundle at ref/highld_exclude_hg38.bed — rsync it to
#           $REF_DIR before running (README §2). If absent the run continues WITHOUT the exclusion and
#           says so loudly, rather than failing the whole grain.
#
# DECISIONS (locked 2026-07-23): run ALL 11 strata (no min-n gate — the >=100-cases GWAS-viability
#   selection is post-hoc at the phenotype boundary). Autosomes only (chrX is a separate later job
#   needing --merge-par + sex-aware handling). Relatives (<=2nd deg) already removed in step 5, so
#   PCA runs on an unrelated set -> clean PCs.
#
# KNOBS
#   AF_EXCLUDE=none     run the whole step unfiltered on purpose (stages B-D skip the exclusion).
#   SKIP_AF_BUILD=1     reuse the exclusion list already on disk instead of rebuilding it.
#   FORCE_QC=1          rebuild stage A even when a valid unfiltered fileset is already present.
#   plus every af_concordance_build knob (THRESH, ZMIN, MIN_CELL, HWE, MISHAP, ...) — passed through.
#
# GUARDRAIL: sbatch script, run by the user (reads genotypes = "the machine"). The AI writes it only.
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
# LD pruning for the PCA input. Window is PHYSICAL (kb), not variant-count: at WGS density
# (~250-350 bp/variant after QC) a 200-VARIANT window spans only ~60 kb, far shorter than the
# few-hundred-kb range over which LD actually extends — so the old "200 50 0.2" pruned almost
# nothing of consequence. 1 Mb window, step 1, r^2 0.1.
LD="1000kb 1 0.1"
THREADS=32

# Long-range-LD / inversion regions, hg38, 0-based BED, chr-prefixed.
# Ships at ref/highld_exclude_hg38.bed; rsync to $REF_DIR. PCA input only.
EXCL=${EXCL:-${HIGHLD_BED}}
if [[ -f "$EXCL" ]]; then
    EXCL_ARG="--exclude bed0 $EXCL"
    EXCL_NOTE="$(wc -l < "$EXCL") regions"
else
    EXCL_ARG=""
    EXCL_NOTE="NONE FOUND at $EXCL — PCs may be driven by MHC/inversion structure"
fi

# Per-callset AF-concordance exclusion. Built by stage B below, applied at stage C, so it shapes
# the association set as well as the PCA input. That is deliberate: a variant mismapped badly
# enough to bend PC1 produces a spurious association in the test itself, where no PC adjustment
# can reach it — filtering only the PCA input would leave it in the sumstats.
AF_EXCLUDE=${AF_EXCLUDE:-${MERGED_DIR}/exclude_af_concordance.txt}
AF_WORK=${MERGED_DIR}/af_concordance
DISC=${DISC:-${MERGED_DIR}/concordance/per_variant_discordance.tsv}

# The 2026-08-17 stale-list incident is now prevented by CONSTRUCTION rather than by a guard:
# stage B rebuilds the list from this merge inside this job, so an old file on disk is
# overwritten before it can be applied. The guard that used to live here (refuse a list older
# than cohort_merged.bed) only mattered because the list arrived from a separate submission.
# It survives in one place only — the SKIP_AF_BUILD path, which is the sole way a pre-existing
# list can still reach stage C.

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
# prune_and_pca <input-stem> <output-dir> <ANC> <n-samples>
# Shared by stage A (unfiltered) and stage D (filtered) so the two generations cannot drift
# apart in their pruning or PCA settings — the whole before/after comparison depends on the
# ONLY difference between them being the exclusion list.
# Echoes "<prune_in> <pcs>" for the summary table.
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

    # Reuse rule, and the reason it is a standing feature rather than a migration flag: stage A's
    # inputs are cohort_merged, the manifest and the locked thresholds. If a valid fileset already
    # postdates the merge, recomputing it reproduces the same bytes at the cost of another full
    # scan. This is also what lets an existing unfiltered pass be adopted without re-running it.
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
    # The only path by which a list this job did not build can reach stage C, so it keeps the
    # provenance guard the rest of the script no longer needs.
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
    [[ "${MISHAP:-0}" != "0" ]] && module load "${MOD_PLINK1}"   # --test-mishap is plink1.9 only

    python3 "${BUNDLE}/scripts/af_concordance_build.py" \
        --qc-dir "${UNF_DIR}" \
        --annot "${ANNOT}" \
        --out "${AF_EXCLUDE}" \
        --work "${AF_WORK}" \
        --thresh "${THRESH:-0.05}" \
        --zmin "${ZMIN:-5.0}" \
        --min-cell "${MIN_CELL:-100}" \
        --hwe "${HWE:-1e-4}" \
        --min-hwe-controls "${MIN_HWE_CONTROLS:-50}" \
        ${HWE_BOTH_TAILS:+--hwe-both-tails} \
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
# STAGE E — assemble BOTH manifests. The unfiltered one exists purely so the effect of the
# filter can be measured rather than asserted; review/plot_af_filter_effect.py reads the pair.
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
