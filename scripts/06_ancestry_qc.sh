#!/bin/bash
#SBATCH --job-name=ancestry_qc
#SBATCH --time=6:00:00
# MEASURED, not guessed: the 2026-07-24 run took ~2.5 h wall for all 11 strata — a flat ~13-14 min
# each, because cost is dominated by scanning cohort_merged (171M x 13,237 = ~527 GiB .bed), which
# every stratum pays regardless of how many samples it keeps. EUR (9.8k samples) cost barely more
# than CAH (71). 6 h is >2x margin; the old 48 h was padding and blocks scheduling.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --partition=norm
#
# STEP 6 — per-ancestry variant QC + within-ancestry PCA.
# Produces the GWAS-ready grain: one QC'd genotype set + PCs per ancestry stratum.
#
# INPUT  : cohort_merged (171M variants, 13,237 samples) + step 5's retained_manifest.csv
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
#     NOTE: the ASSOCIATION set is the full QC-passing variants (step 2 output); the LD-pruned set
#           is used ONLY for PCA (standard practice).
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
# GUARDRAIL: sbatch script, run by the user (reads genotypes = "the machine"). The AI writes it only.
#   ./submit.sh scripts/06_ancestry_qc.sh
set -o pipefail

# Resolve the bundle root and load the one file that knows where anything lives.
BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

MERGED=${MERGED_DIR}/cohort_merged
MANIFEST=${MERGED_DIR}/relatedness/retained_manifest.csv
OUT_DIR=${MERGED_DIR}/by_ancestry_qc
KEEP_DIR=${OUT_DIR}/keep
LOG=${OUT_DIR}/step6_summary.txt
mkdir -p "$OUT_DIR" "$KEEP_DIR"

# ── QC thresholds (locked) ──
GENO=0.05
MAF=0.01
HWE=1e-6
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

# Per-callset AF-concordance exclusion (built by scripts/af_concordance_build.sh). Applied at the
# VARIANT QC pass, so it shapes the association set as well as the PCA input ($PRUNE reads from
# $QC). That is deliberate: a
# variant mismapped badly enough to bend PC1 produces a spurious association in the test itself,
# where no PC adjustment can reach it — filtering only the PCA input would leave it in the sumstats.
# Absent => run unfiltered and say so, rather than failing the grain (same policy as $EXCL above).
AF_EXCLUDE=${AF_EXCLUDE:-${MERGED_DIR}/exclude_af_concordance.txt}
if [[ -f "$AF_EXCLUDE" ]]; then
    AFX_ARG="--exclude $AF_EXCLUDE"
    AFX_NOTE="$(wc -l < "$AF_EXCLUDE") variants from $AF_EXCLUDE"
else
    AFX_ARG=""
    AFX_NOTE="NONE at $AF_EXCLUDE — run scripts/af_concordance_build.sh first, or accept that the"
    AFX_NOTE="$AFX_NOTE callset AF divergence stays in both the PCs and the association set"
fi

module load "${MOD_PLINK2}"

echo "==========================================" | tee "$LOG"
echo "Step 6 — per-ancestry variant QC + PCA"     | tee -a "$LOG"
echo "Job: ${SLURM_JOB_ID}  Node: ${SLURMD_NODENAME}  Start: $(date)" | tee -a "$LOG"
echo "Merged: ${MERGED}   Manifest: ${MANIFEST}"   | tee -a "$LOG"
echo "geno<=${GENO} maf>=${MAF} hwe>=${HWE} pcs=${NPC} ld='${LD}'" | tee -a "$LOG"
echo "PCA-input long-range-LD exclusion: ${EXCL_NOTE}" | tee -a "$LOG"
echo "AF-concordance exclusion (assoc + PCA): ${AFX_NOTE}" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"

if [[ ! -f "${MANIFEST}" ]]; then echo "ERROR: manifest not found: ${MANIFEST}" >&2; exit 1; fi
if [[ ! -f "${MERGED}.bed" ]]; then echo "ERROR: merged bed not found: ${MERGED}.bed" >&2; exit 1; fi

# ancestries present in the retained manifest (col 3; skip header + blanks)
ANCS=$(awk -F, 'NR>1 && $3!="" {print $3}' "$MANIFEST" | sort -u)
echo "Ancestries: $(echo $ANCS | tr '\n' ' ')" | tee -a "$LOG"

printf "%-5s %8s %12s %10s %5s\n" "ANC" "samples" "qc_variants" "prune_in" "pcs" | tee -a "$LOG"

for ANC in $ANCS; do
    KEEP=${KEEP_DIR}/keep_${ANC}.txt
    awk -F, -v a="$ANC" 'NR>1 && $3==a {print $1"\t"$2}' "$MANIFEST" > "$KEEP"
    N=$(wc -l < "$KEEP")
    if [[ "$N" -lt 2 ]]; then
        printf "%-5s %8s %12s %10s %5s\n" "$ANC" "$N" "SKIP(<2)" "-" "-" | tee -a "$LOG"
        continue
    fi

    QC=${OUT_DIR}/cohort_${ANC}_qc
    PRUNE=${OUT_DIR}/${ANC}_prune
    PCA=${OUT_DIR}/cohort_${ANC}_pca

    # ── extract + variant QC (autosomes), one pass ──
    plink2 --bfile "$MERGED" --keep "$KEEP" \
           --autosome $AFX_ARG \
           --geno "$GENO" --maf "$MAF" --hwe "$HWE" keep-fewhet \
           --make-bed --threads "$THREADS" --out "$QC" \
        || { printf "%-5s %8s %12s %10s %5s\n" "$ANC" "$N" "QC_FAILED" "-" "-" | tee -a "$LOG"; continue; }
    NVAR=$(wc -l < "${QC}.bim")

    # ── LD-prune (PCA input only; long-range-LD regions dropped here, NOT from $QC) ──
    plink2 --bfile "$QC" $EXCL_ARG --indep-pairwise $LD --threads "$THREADS" --out "$PRUNE" \
        || { printf "%-5s %8s %12s %10s %5s\n" "$ANC" "$N" "$NVAR" "PRUNE_FAIL" "-" | tee -a "$LOG"; continue; }
    NIN=$(wc -l < "${PRUNE}.prune.in" 2>/dev/null || echo 0)

    # ── PCA on pruned set (cap #PCs for small strata; need N>npcs and enough variants) ──
    PCS="-"
    if [[ "$N" -ge 3 && "$NIN" -ge 2 ]]; then
        NPC_EFF=$NPC
        [[ "$N" -le "$NPC" ]] && NPC_EFF=$((N-1))
        if plink2 --bfile "$QC" --extract "${PRUNE}.prune.in" \
                  --pca "$NPC_EFF" --threads "$THREADS" --out "$PCA"; then
            PCS=$NPC_EFF
        else
            PCS="PCA_FAIL"
        fi
    else
        PCS="too_small"
    fi

    printf "%-5s %8s %12s %10s %5s\n" "$ANC" "$N" "$NVAR" "$NIN" "$PCS" | tee -a "$LOG"
done

echo "==========================================" | tee -a "$LOG"
echo "Step 6 complete: $(date)" | tee -a "$LOG"
echo "Outputs in ${OUT_DIR}: cohort_<ANC>_qc.{bed,bim,fam} (association set) + cohort_<ANC>_pca.eigenvec/.eigenval (covariates)" | tee -a "$LOG"
