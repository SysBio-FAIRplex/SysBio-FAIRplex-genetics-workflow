#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH --time=6:00:00
# MEASURED: the 2026-07-25 run did 44 contrasts in 100 min. The current script adds two --missing
# passes per contrast (differential missingness) plus awk passes over each ~700 MB output for
# BETA_MAX and hits, so budget ~3 h. 6 h is 2x margin.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=norm
#
# STEP 7 — per-ancestry, per-contrast GWAS on the QC'd grain.
#
# For each viable ancestry (EUR/AJ/AAC/AFR/AMR/CAH — the strata step 6 produced PCs for) and each
# case/control contrast on dx_detailed, this: builds a covar file (sex + PC1..PC10) and a case/control
# pheno file straight from analysis_grain.csv, then runs plink2 --glm on cohort_<ANC>_qc. Every
# attempt is logged to gwas_summary.csv with arm counts, the empirical cohort composition of each arm,
# a data-driven confound tag, and the >=100-cases/arm viability flag — so the post-hoc filter and the
# interpretation guide are computed, never hand-counted.
#
# CONTRASTS (case arm listed FIRST -> coded 2; control arm -> 1). An arm may be source-restricted as
#   "<dx>@amppd" (wb_dwgs) or "<dx>@ampad" (the AMP-AD callsets); bare dx = any source.
#   dx contrasts:  PD-vs-AD ; PD-vs-{DLB,MCI,PSP,control,other} ; AD-vs-{MCI,DLB,PSP,control,other}
#   control-vs-control (batch-artifact scan): control@amppd:control@ampad  — both arms disease-free, so
#     any hit is a cohort/batch frequency artifact; feed those variants to mask_cohort_artifacts.py to
#     clean the cross-cohort disease sumstats (post-hoc; no plink --exclude re-run).
#   within-cohort (confound-free — the trusted backbone): PD@amppd:control@amppd ; AD@ampad:control@ampad
#
# MODEL: plink2 --glm with firth-fallback (auto logistic for a 1/2 binary pheno; Firth on
#   quasi-separated SNPs -> robust for the imbalanced arms), covariates sex + PC1..PC10,
#   --covar-variance-standardize. NO mixed model: <=2nd-degree relatives were removed in step 5, so
#   samples are independent. Per-ancestry -> ancestry-stratified by construction.
#   --keep restricts to the contrast's arms, so --maf MIN_MAF (default 0.05) and the association are
#   computed on the tested samples (not the pooled stratum). MAF 0.05: at these arm sizes MAF<0.05 is
#   unstable + batch-artifact-prone (measured); rare/low-freq belongs in a separate burden pass.
#
# DIFFERENTIAL MISSINGNESS (pre-association filter, per contrast): the primary defence against the
#   cross-cohort confound. A variant called well in one sequencing program and poorly in the other
#   produces a clean, highly significant, entirely artifactual association — and because disease is
#   entangled with program here, that is the default failure mode, not an edge case. For each
#   contrast we compute per-arm variant missingness and drop variants failing BOTH tests: a 2x2
#   chi-square exceeding DIFFMISS_P (default 1e-4) AND an absolute arm difference of at least
#   DIFFMISS_MINDIFF (default 0.02). Requiring both is deliberate — a chi-square alone scales with
#   n, so at EUR's 10,059 samples a negligible difference clears any threshold while at CAH's 105 a
#   large one may not, and the filter would then mean something different in every stratum. 0.02 is
#   a meaningful share of the 0.05 missingness ceiling --geno already imposes. This catches the
#   mechanism directly and BEFORE the association, where the control-vs-control scan only catches
#   artifacts that happen to reach significance in that particular comparison.
#   Implemented via two --missing passes + an awk chi-square rather than --test-missing, which is a
#   plink1.9 feature and is not relied on being present in this plink2 build.
#
# AGE COVARIATE: included when the grain carries an age column, otherwise omitted with a warning.
#   Age is the dominant confounder for both AD and PD, so its absence is a real limitation — but the
#   cohorts do not currently supply a commensurable variable (AMP-AD gives age at death, AMP-PD gives
#   age at baseline/analysis), so this cannot simply be forced. Until the phenotype track emits a
#   harmonized age, the within-cohort contrasts are the only ones where age would be comparable
#   anyway. Detected by header name: age | age_analysis | agedeath | age_baseline | age_cov.
#
# CONFOUND (interpretation guide, NOT a filter): disease is confounded with study program
#   (AMP-PD = wb_dwgs ; AMP-AD = wgs_harm/divco_hs/fused). The summary reports each arm's %AMP-PD and
#   delta_amppd = |case%AMP-PD - ctrl%AMP-PD|. Heuristic tag: within_cohort (delta<=20, cleanest) /
#   partial / cross_cohort (delta>=70, interpret cautiously — e.g. PD-vs-AD and, per the data,
#   AD-vs-control since controls are ~88% AMP-PD).
#
# DUAL-SOURCE 87 (source_callset == "divco_hs|wgs_harm"): these are same-IID DivCo+WGS_Harm genomes
#   FUSED at merge time (one row, not a duplicate). Decision (2026-07-24): accept fused as-is; their
#   DivCo-only sites are dropped by the per-ancestry --geno anyway (they're a tiny minority of EUR).
#   Set EXCLUDE_DUAL=1 to drop the 87 from all pheno files for a cheap WGS_Harm-only sensitivity run
#   (no re-merge needed).
#
# RUN POLICY: attempt every ANC x contrast whose BOTH arms have >= MIN_ARM samples (default 20, for
#   model convergence); flag >=100-both as the GWAS-viability cut. MIN_ARM=0 forces literally all.
#
# GUARDRAIL: sbatch script, run by the USER (reads genotypes + the id-bearing grain = "the machine").
#   The AI writes it only. Sumstats + summary are the deliverables; stdout is aggregate counts.
#
# PREREQ: analysis_grain.csv must be at $GRAIN and must carry the PCs from the CURRENT step-6 run.
#   It is written by clinical_core.ipynb §12. Re-running step 6 invalidates it — regenerate and
#   re-copy, or this silently runs on stale covariates. See README §3 step 7.
#
#   ./submit.sh scripts/07_gwas.sh
#   (sensitivity: ./submit.sh scripts/07_gwas.sh --export=EXCLUDE_DUAL=1,OUT_SUBDIR=gwas_nodual)
set -o pipefail

# Resolve the bundle root and load the one file that knows where anything lives.
BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

QC_DIR=${MERGED_DIR}/by_ancestry_qc
OUT_DIR=${QC_DIR}/${OUT_SUBDIR:-gwas}
PHENO_DIR=${OUT_DIR}/pheno
COVAR_DIR=${OUT_DIR}/covar
SUMMARY=${OUT_DIR}/gwas_summary.csv
mkdir -p "$OUT_DIR" "$PHENO_DIR" "$COVAR_DIR"

# ── knobs (locked defaults; overridable via --export) ──
ANCS=${ANCS:-"EUR AJ AAC AFR AMR CAH"}                 # the 6 strata with PCs
CONTRASTS=${CONTRASTS:-"PD:AD PD:DLB PD:MCI PD:PSP PD:control PD:other AD:MCI AD:DLB AD:PSP AD:control AD:other control@amppd:control@ampad PD@amppd:control@amppd AD@ampad:control@ampad"}
MIN_ARM=${MIN_ARM:-20}                                  # both arms must reach this to attempt
MIN_MAF=${MIN_MAF:-0.05}                                 # common-variant floor. At these arm sizes MAF<0.05
                                                         # has too few minor alleles for a stable single-variant
                                                         # test AND is where cross-cohort batch artifacts
                                                         # concentrate (measured). Computed per-contrast (--keep
                                                         # the arms), so it's MAF in the tested samples, not the
                                                         # pooled stratum. Rare/low-freq -> separate burden pass.
EXCLUDE_DUAL=${EXCLUDE_DUAL:-0}                          # 1 -> drop the fused 87 (WGS_Harm-only sensitivity)
DIFFMISS_P=${DIFFMISS_P:-1e-4}                           # per-contrast differential-missingness cut
DIFFMISS_MINDIFF=${DIFFMISS_MINDIFF:-0.02}               # AND: minimum |arm missingness difference|.
                                                         # A p-value alone does not transfer across strata —
                                                         # at EUR's n=10,059 a trivial difference clears any
                                                         # threshold, while at CAH's n=105 a large one may not.
                                                         # Requiring practical as well as statistical
                                                         # significance keeps the filter comparable across
                                                         # strata. 0.02 is meaningful against the 0.05 ceiling
                                                         # --geno already imposes. Set 0 for p-value only.
GWSIG=${GWSIG:-5e-8}                                     # genome-wide significance for the hits file
BETA_MAX=${BETA_MAX:-5}                                  # |BETA| sanity bound (non-converged Firth)
THREADS=32
MISS_DIR=${OUT_DIR}/diffmiss
mkdir -p "$MISS_DIR"

# chi-square (1 df) critical value for DIFFMISS_P. Tabulated rather than computed so the filter
# has no scipy/python dependency inside the loop.
case "$DIFFMISS_P" in
    1e-3) DIFFMISS_CHI2=10.8276 ;;
    1e-4) DIFFMISS_CHI2=15.1367 ;;
    1e-5) DIFFMISS_CHI2=19.5114 ;;
    1e-6) DIFFMISS_CHI2=23.9281 ;;
    off)  DIFFMISS_CHI2="" ;;
    *) echo "ERROR: DIFFMISS_P must be one of 1e-3 1e-4 1e-5 1e-6 off (got '$DIFFMISS_P')" >&2; exit 2 ;;
esac

module load plink/6-alpha

# ── genomic inflation from a plink2 sumstats file (guardrail-safe: variant-level, no subject rows) ──
# lambda_GC = median(chi2_obs)/0.4549, with chi2 = Z_STAT^2 taken straight from plink2 (EXACT — no
# inverse-normal approximation). Rows with a non-numeric Z (e.g. an unconverged Firth SNP) are skipped;
# negligible against millions of variants. If the build emits no Z_STAT/T_STAT column at all, this
# prints "NA NA" and lambda is instead recomputed exactly (scipy) by the local plot_gwas.py.
# lambda_1000 rescales to a 1000/1000 study so lambdas are comparable across arm sizes; it is emitted
# only when the smaller arm >=200 (the rescale is unstable below that), else NA.
# Args: <sumstats_file> <n_case> <n_ctrl>  -> prints "lambda_gc lambda_1000" (or "NA NA").
lambda_from_sumstats() {
    local FILE=$1 NC=$2 NK=$3
    [[ -f "$FILE" ]] || { echo "NA NA"; return; }
    awk 'NR==1{ for(i=1;i<=NF;i++) h[$i]=i
                zc=(h["Z_STAT"]?h["Z_STAT"]:h["T_STAT"]); tc=h["TEST"]; next }
         zc { if(tc && $tc!="ADD") next; z=$zc; if(z=="NA"||z=="") next; print (z*z) }' "$FILE" \
    | sort -g \
    | awk -v nc="$NC" -v nk="$NK" '
        {v[NR]=$1}
        END{
            n=NR; if(n==0){print "NA NA"; exit}
            med=(n%2)?v[(n+1)/2]:(v[n/2]+v[n/2+1])/2
            lam=med/0.454936
            # lambda_1000 rescale factor ~500*(1/nc+1/nk) explodes for small arms -> only report it
            # when the smaller arm >=200 and lambda is sane; else NA (uninterpretable, not "broken").
            minarm=(nc<nk?nc:nk)
            if(nc>0 && nk>0 && minarm>=200 && lam>0){
                l=1+(lam-1)*((1/nc)+(1/nk))/(2/1000)
                l1000=(l<0)?"NA":sprintf("%.4f", l)
            } else l1000="NA"
            printf "%.4f %s\n", lam, l1000
        }'
}

echo "=========================================="
echo "Step 7 — per-ancestry per-contrast GWAS"
echo "Job: ${SLURM_JOB_ID:-NA}  Node: ${SLURMD_NODENAME:-NA}  Start: $(date)"
echo "Grain: ${GRAIN}"
echo "QC dir: ${QC_DIR}   Out: ${OUT_DIR}"
echo "Ancestries: ${ANCS}"
echo "MIN_ARM=${MIN_ARM}  EXCLUDE_DUAL=${EXCLUDE_DUAL}"
echo "=========================================="

[[ -f "$GRAIN" ]] || { echo "ERROR: grain not found: $GRAIN (rsync it from local first)" >&2; exit 1; }

# grain columns (1-based): 1 IID,2 individual_id,3 source_callset,4 ancestry,5 sex,6 pheno,
#   7 dx_detailed,8 pheno_conflict,9 dx_conflict,10 sex_conflict,11 call_rate,12 dup_cluster_id,
#   13..22 PC1..PC10
DUAL="divco_hs|wgs_harm"

# ── age covariate: present only if the grain carries a recognizable age column ──
AGE_IDX=$(head -1 "$GRAIN" | awk -F, '{for(i=1;i<=NF;i++){c=tolower($i); gsub(/^[ \t]+|[ \t]+$/,"",c);
    if(c=="age"||c=="age_analysis"||c=="agedeath"||c=="age_death"||c=="age_baseline"||c=="age_cov"){print i; exit}}}')
AGE_IDX=${AGE_IDX:-0}
if [[ "$AGE_IDX" -gt 0 ]]; then
    echo "age covariate: grain column ${AGE_IDX} — included"
else
    echo "age covariate: NOT FOUND in grain header — running WITHOUT age."
    echo "  Age is the dominant confounder for AD and PD. This is a known limitation, blocked on the"
    echo "  phenotype track emitting a harmonized age (AMP-AD age-at-death vs AMP-PD age-at-baseline"
    echo "  are not the same variable). Interpret accordingly."
fi

# summary header
echo "ancestry,contrast,case_arm,ctrl_arm,n_case,n_ctrl,case_pct_amppd,ctrl_pct_amppd,delta_amppd,confound_tag,viable_ge100,ran,lambda_gc,lambda_1000,n_diffmiss_excluded,n_beta_outlier,n_gwsig,sumstats_file" > "$SUMMARY"

printf "%-5s %-16s %7s %7s %8s %-13s %-7s %-12s %8s %8s %9s %6s\n" "ANC" "contrast" "n_case" "n_ctrl" "delta" "confound" "viable" "status" "lam_gc" "lam_1k" "diffmiss" "hits"

for ANC in $ANCS; do
    QC=${QC_DIR}/cohort_${ANC}_qc
    if [[ ! -f "${QC}.bed" ]]; then
        echo "  skip ${ANC}: no ${QC}.bed"
        continue
    fi

    # ── IID->FID map from the genotype fileset (source of truth for the FID plink matches on) ──
    # plink2 defaults a pheno/covar file's missing FID to "0" and matches on FID+IID, so the AMP-PD
    # (wb_dwgs) samples — which carry a non-zero FID — need their real FID or they never match.
    FAMMAP=${COVAR_DIR}/fid_${ANC}.map
    awk '{print $2"\t"$1}' "${QC}.fam" > "$FAMMAP"      # IID<tab>FID (tab-delim, no commas)

    # ── covar file for this ancestry (built once): #FID IID SEX PC1..PC10 ──
    COVAR=${COVAR_DIR}/covar_${ANC}.txt
    awk -F, -v A="$ANC" -v AGEI="$AGE_IDX" '
        NR==FNR { split($1,m,"\t"); fid[m[1]]=m[2]; next }   # map file (first arg): IID -> FID
        FNR==1 { h="#FID\tIID\tSEX"; if(AGEI>0) h=h "\tAGE";
                 for(p=1;p<=10;p++) h=h "\tPC" p; print h; next }
        $4==A {
            f=(($1 in fid)?fid[$1]:"0"); sex=($5==""?"NA":$5);
            line=f "\t" $1 "\t" sex;
            if(AGEI>0){ a=$AGEI; line=line "\t" (a==""?"NA":a) }
            for(i=13;i<=22;i++){ v=($i==""?"NA":$i); line=line "\t" v }
            print line
        }' "$FAMMAP" "$GRAIN" > "$COVAR"

    for C in $CONTRASTS; do
        CASE_ARM=${C%%:*}
        CTRL_ARM=${C##*:}
        # Optional per-arm source restriction "<dx>@amppd" | "<dx>@ampad" (amppd=wb_dwgs, ampad=rest).
        # Enables within-cohort scans and the control-vs-control batch-artifact scan; bare dx = any source.
        CASE_DX=${CASE_ARM%%@*}; CASE_SRC=${CASE_ARM#*@}; [ "$CASE_SRC" = "$CASE_ARM" ] && CASE_SRC=""
        CTRL_DX=${CTRL_ARM%%@*}; CTRL_SRC=${CTRL_ARM#*@}; [ "$CTRL_SRC" = "$CTRL_ARM" ] && CTRL_SRC=""
        TAG=$(printf '%s_vs_%s' "$CASE_ARM" "$CTRL_ARM" | tr '@' '_')   # filesystem-safe (@ -> _)
        PHENO=${PHENO_DIR}/pheno_${ANC}_${TAG}.txt

        # ── build pheno (#FID IID; case=2, ctrl=1) + count arms + per-arm AMP-PD share, one awk pass ──
        read -r NCASE NCTRL CASE_PD CTRL_PD < <(
            awk -F, -v A="$ANC" -v CDX="$CASE_DX" -v CSRC="$CASE_SRC" -v KDX="$CTRL_DX" -v KSRC="$CTRL_SRC" \
                    -v DUAL="$DUAL" -v EXD="$EXCLUDE_DUAL" -v PF="$PHENO" '
                function srcok(want, ispd){ return (want=="" || (want=="amppd"&&ispd) || (want=="ampad"&&!ispd)) }
                NR==FNR { split($1,m,"\t"); fid[m[1]]=m[2]; next }   # map file (first arg): IID -> FID
                FNR==1 { print "#FID\tIID\tpheno" > PF; next }
                $4==A {
                    if (EXD==1 && $3==DUAL) next
                    ispd = ($3=="wb_dwgs") ? 1 : 0
                    f=(($1 in fid)?fid[$1]:"0")
                    if ($7==CDX && srcok(CSRC,ispd)){ ncase++; case_pd+=ispd; print f"\t"$1"\t"2 > PF }
                    else if ($7==KDX && srcok(KSRC,ispd)){ nctrl++; ctrl_pd+=ispd; print f"\t"$1"\t"1 > PF }
                }
                END{ printf "%d %d %d %d\n", ncase, nctrl, case_pd, ctrl_pd }
            ' "$FAMMAP" "$GRAIN"
        )
        NCASE=${NCASE:-0}; NCTRL=${NCTRL:-0}; CASE_PD=${CASE_PD:-0}; CTRL_PD=${CTRL_PD:-0}

        # ── cohort-composition + confound tag (data-driven) ──
        CASE_PCT=$(awk -v n="$NCASE" -v p="$CASE_PD" 'BEGIN{printf (n>0)?"%.1f":"NA", (n>0)?100*p/n:0}')
        CTRL_PCT=$(awk -v n="$NCTRL" -v p="$CTRL_PD" 'BEGIN{printf (n>0)?"%.1f":"NA", (n>0)?100*p/n:0}')
        DELTA=$(awk -v a="$CASE_PCT" -v b="$CTRL_PCT" 'BEGIN{ if(a=="NA"||b=="NA"){print "NA"} else {d=a-b; print (d<0?-d:d)} }')
        TAGC=$(awk -v d="$DELTA" 'BEGIN{ if(d=="NA"){print "NA"} else if(d<=20){print "within_cohort"} else if(d>=70){print "cross_cohort"} else {print "partial"} }')

        VIABLE="no"; [[ "$NCASE" -ge 100 && "$NCTRL" -ge 100 ]] && VIABLE="yes"

        # ── run gate: both arms >= MIN_ARM ──
        SUMFILE=""; LAM_GC="NA"; LAM_1K="NA"; NDIFF="NA"; NBETA="NA"; NHIT="NA"
        if [[ "$NCASE" -ge "$MIN_ARM" && "$NCTRL" -ge "$MIN_ARM" ]]; then
            OUT=${OUT_DIR}/gwas_${ANC}_${TAG}

            # ── differential missingness between the arms (pre-association) ──
            # Two --missing passes, then a 2x2 chi-square on (missing, nonmissing) x (case, ctrl).
            # Variants above the critical value are excluded from THIS contrast only.
            DM_EXCL=""
            if [[ -n "$DIFFMISS_CHI2" ]]; then
                DMB=${MISS_DIR}/${ANC}_${TAG}
                awk 'NR>1 && $3==2 {print $1"\t"$2}' "$PHENO" > "${DMB}.case.keep"
                awk 'NR>1 && $3==1 {print $1"\t"$2}' "$PHENO" > "${DMB}.ctrl.keep"
                if plink2 --bfile "$QC" --keep "${DMB}.case.keep" --missing variant-only \
                          --threads "$THREADS" --out "${DMB}.case" >/dev/null 2>&1 \
                && plink2 --bfile "$QC" --keep "${DMB}.ctrl.keep" --missing variant-only \
                          --threads "$THREADS" --out "${DMB}.ctrl" >/dev/null 2>&1; then
                    # .vmiss: #CHROM ID MISSING_CT OBS_CT F_MISS   (F_MISS = MISSING_CT/OBS_CT)
                    # Drop a variant only if the arm difference is BOTH statistically detectable
                    # (chi-square) AND practically large (|F_MISS_case - F_MISS_ctrl|). Either test
                    # alone misbehaves: chi-square scales with n, the raw difference ignores it.
                    awk -v CRIT="$DIFFMISS_CHI2" -v MIND="$DIFFMISS_MINDIFF" '
                        FNR==1 { for(i=1;i<=NF;i++) h[$i]=i
                                 idc=h["ID"]; mc=h["MISSING_CT"]; oc=h["OBS_CT"]; next }
                        NR==FNR { miss[$idc]=$mc; obs[$idc]=$oc; next }
                        ($idc in miss) {
                            a=miss[$idc]; b=obs[$idc]-a          # case: missing, nonmissing
                            c=$mc;        d=$oc-c                # ctrl: missing, nonmissing
                            n=a+b+c+d
                            if(n==0) next
                            r1=a+b; r2=c+d; c1=a+c; c2=b+d
                            if(r1==0||r2==0||c1==0||c2==0) next  # no missingness anywhere -> nothing to test
                            num=(a*d-b*c); chi=n*num*num/(r1*r2*c1*c2)
                            fcase=(r1>0)?a/r1:0; fctrl=(r2>0)?c/r2:0
                            diff=fcase-fctrl; if(diff<0) diff=-diff
                            if(chi>CRIT && diff>=MIND) print $idc
                        }' "${DMB}.case.vmiss" "${DMB}.ctrl.vmiss" > "${DMB}.exclude"
                    NDIFF=$(wc -l < "${DMB}.exclude")
                    [[ "$NDIFF" -gt 0 ]] && DM_EXCL="--exclude ${DMB}.exclude"
                else
                    NDIFF="MISS_FAILED"
                fi
            fi

            if plink2 --bfile "$QC" \
                      --keep "$PHENO" \
                      $DM_EXCL \
                      --maf "$MIN_MAF" \
                      --pheno "$PHENO" --pheno-name pheno \
                      --covar "$COVAR" --covar-variance-standardize \
                      --glm hide-covar firth-fallback cols=+a1freq,+beta \
                      --threads "$THREADS" --out "$OUT" > "${OUT}.console" 2>&1; then
                # plink2 names it .pheno.glm.logistic.hybrid (firth-fallback) or .glm.logistic
                SUMFILE=$(ls "${OUT}".pheno.glm.logistic* 2>/dev/null | head -1)
                if [[ -n "$SUMFILE" ]]; then
                    RAN="yes"
                    read -r LAM_GC LAM_1K < <(lambda_from_sumstats "$SUMFILE" "$NCASE" "$NCTRL")
                    # ── |BETA| sanity pass: emit an ADD-only, outlier-free copy ──
                    # Firth fallback can return implausible effect sizes on unconverged fits; those
                    # are numerical failures, not findings. Original sumstats are left untouched.
                    NBETA=$(awk -v FILT="${OUT}.filtered.tsv" -v BMAX="$BETA_MAX" '
                        NR==1 { for(i=1;i<=NF;i++) h[$i]=i; bc=h["BETA"]; tc=h["TEST"]
                                print > FILT; next }
                        { if(tc && $tc!="ADD") next
                          v=$bc; if(v=="NA"||v=="") { print > FILT; next }
                          if(v>BMAX || v<-BMAX) { n++; next }
                          print > FILT }
                        END { print n+0 }' "$SUMFILE")

                    # ── genome-wide significant hits, from the |BETA|-filtered ADD set ──
                    # Drawn from the filtered copy on purpose: an unconverged Firth fit can carry
                    # a tiny p-value, and those are numerical failures, not associations.
                    NHIT=$(awk -v HITS="${OUT}.hits.tsv" -v SIG="$GWSIG" '
                        NR==1 { for(i=1;i<=NF;i++) h[$i]=i; pc=h["P"]; print > HITS; next }
                        { v=$pc; if(v=="NA"||v=="") next
                          if(v+0 < SIG+0) { n++; print > HITS } }
                        END { print n+0 }' "${OUT}.filtered.tsv")
                else
                    RAN="no_output"
                fi
            else
                RAN="FAILED"
            fi
        else
            RAN="skip_min_arm"
        fi

        printf "%-5s %-16s %7s %7s %8s %-13s %-7s %-12s %8s %8s %9s %6s\n" \
            "$ANC" "$TAG" "$NCASE" "$NCTRL" "$DELTA" "$TAGC" "$VIABLE" "$RAN" "$LAM_GC" "$LAM_1K" "$NDIFF" "$NHIT"
        echo "${ANC},${TAG},${CASE_ARM},${CTRL_ARM},${NCASE},${NCTRL},${CASE_PCT},${CTRL_PCT},${DELTA},${TAGC},${VIABLE},${RAN},${LAM_GC},${LAM_1K},${NDIFF},${NBETA},${NHIT},${SUMFILE##*/}" >> "$SUMMARY"
    done
done

echo "=========================================="
echo "GWAS complete: $(date)"
echo "Summary: ${SUMMARY}"
echo "Sumstats: ${OUT_DIR}/gwas_<ANC>_<CASE>_vs_<CTRL>.pheno.glm.logistic.hybrid"
echo "  + .filtered.tsv (ADD only, |BETA|<=${BETA_MAX})   + .hits.tsv (P<${GWSIG}, from the filtered set)"
echo "Viable (>=100/arm) rows  [ANC contrast (case/ctrl) confound  lambda_gc/lambda_1000]:"
awk -F, 'NR>1 && $11=="yes"{printf "  %-4s %-16s (%s/%s)  %-13s  lam=%s / l1000=%s\n",$1,$2,$5,$6,$10,$13,$14}' "$SUMMARY"
