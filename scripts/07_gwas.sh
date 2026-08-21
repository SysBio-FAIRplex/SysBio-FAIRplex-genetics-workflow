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
# STEP 7 — per-ancestry, per-contrast GWAS. Model, thresholds and rationale: METHODS.md §7.
#
# For each row of $CONTRASTS_CSV whose stratum has a QC'd fileset, runs plink2 --glm on
# cohort_<ANC>_qc using the pheno and covar files analysis_grain.py §13 wrote. Per-contrast results
# and arm composition go to gwas_summary.csv.
#
# READS §13'S FILES; DOES NOT REBUILD THEM. Three consequences that bite if forgotten:
#   * MIN_ARM can only RAISE the floor — §13 writes no file below 20 per arm.
#   * EXCLUDE_DUAL is not accepted here (§13 owns it). Passing it is a hard error, not ignored.
#   * The age covariate is decided in §13; this reports what covar_<ANC>.txt actually contains.
#
# PREREQ, CHECKED NOT WARNED: each covar file must be newer than its stratum's .eigenvec (exit 4),
# and every pheno IID must appear in the .fam being tested (exit 5) — which catches a mismatched
# generation that mtime alone would miss. contrasts.csv is read by header NAME, so a reordered
# column aborts naming the missing one (exit 3).
#
#   ./submit.sh scripts/07_gwas.sh
set -o pipefail

# Resolve the bundle root and load the one file that knows where anything lives.
BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

QC_DIR=${MERGED_DIR}/by_ancestry_qc
OUT_DIR=${QC_DIR}/${OUT_SUBDIR:-gwas}
SUMMARY=${OUT_DIR}/gwas_summary.csv
mkdir -p "$OUT_DIR"

# ── knobs (locked defaults; overridable via --export) ──
ANCS=${ANCS:-"EUR AJ AAC AFR AMR CAH"}   # the 6 strata with PCs. Filters contrasts.csv.
CONTRASTS=${CONTRASTS:-""}               # optional filter: contrast tags as they appear in
                                         # contrasts.csv. Empty = every row for the selected strata.
MIN_ARM=${MIN_ARM:-20}                   # RAISES the floor only — §13 wrote nothing below 20.
MIN_MAF=${MIN_MAF:-0.05}                 # computed per-contrast (--keep the arms), not pooled.
DIFFMISS_P=${DIFFMISS_P:-1e-4}           # differential missingness: BOTH this AND MINDIFF must fail
DIFFMISS_MINDIFF=${DIFFMISS_MINDIFF:-0.02}  # a variant. Requiring both is what makes the filter mean
                                         # the same thing at EUR's n and at CAH's. 0 = p-value only.
GWSIG=${GWSIG:-5e-8}
BETA_MAX=${BETA_MAX:-5}                  # |BETA| sanity bound (non-converged Firth)
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

# ── genomic inflation. Args: <sumstats_file> <n_case> <n_ctrl> -> "lambda_gc lambda_1000" ──
# lambda_GC = median(Z_STAT^2)/0.4549, exact — no inverse-normal approximation. A build emitting no
# Z_STAT/T_STAT prints "NA NA"; plot_gwas.py recomputes it with scipy in that case.
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
echo "Contrasts: ${CONTRASTS_CSV}"
echo "Pheno:     ${PHENO_SRC}"
echo "Covar:     ${COVAR_SRC}"
echo "QC dir: ${QC_DIR}   Out: ${OUT_DIR}"
echo "Ancestries: ${ANCS}   MIN_ARM=${MIN_ARM}"
echo "=========================================="

# ── EXCLUDE_DUAL moved to §13. REFUSE rather than ignore: a sensitivity run that silently came ──
# ── back as the primary is worse than a failed submission.                                     ──
if [[ -n "${EXCLUDE_DUAL:-}" && "${EXCLUDE_DUAL}" != "0" ]]; then
    echo "ERROR: EXCLUDE_DUAL is no longer handled here (removed 2026-08-20)." >&2
    echo "  Arms are defined once, in analysis_grain.py §13. Set EXCLUDE_DUAL=True there, rerun it" >&2
    echo "  to a separate CLINICAL_OUT, then point this script at that generation:" >&2
    echo "  --export=PHENO_SRC=<dir>/pheno,COVAR_SRC=<dir>/covar,CONTRASTS_CSV=<dir>/contrasts.csv,OUT_SUBDIR=gwas_nodual" >&2
    exit 2
fi

[[ -f "$CONTRASTS_CSV" ]] || { echo "ERROR: no contrasts.csv at $CONTRASTS_CSV — run analysis_grain.py (§13) after step 6" >&2; exit 1; }
[[ -d "$PHENO_SRC" ]]     || { echo "ERROR: no pheno dir at $PHENO_SRC — run analysis_grain.py (§13)" >&2; exit 1; }
[[ -d "$COVAR_SRC" ]]     || { echo "ERROR: no covar dir at $COVAR_SRC — run analysis_grain.py (§13)" >&2; exit 1; }

# ── read contrasts.csv BY HEADER NAME (rule 4) and emit a fixed-order tab stream ──
# Positional would mislabel a contrast the moment §13 reorders a column. Tab-separated because no
# field can contain a tab.
CROWS=${OUT_DIR}/.contrast_rows.tsv
awk -F, -v ANCS=" ${ANCS} " -v WANT=" ${CONTRASTS} " '
    NR==1 { for(i=1;i<=NF;i++){ k=$i; gsub(/^[ \t]+|[ \t]+$/,"",k); h[k]=i }
            need="ancestry contrast case_arm ctrl_arm n_case n_ctrl case_pct_amppd ctrl_pct_amppd delta_amppd confound_tag max_callset_delta worst_callset callset_one_sided viable_ge100"
            n=split(need,w," ")
            for(j=1;j<=n;j++) if(!(w[j] in h)){ printf "MISSING_COLUMN\t%s\n", w[j] > "/dev/stderr"; bad=1 }
            if(bad) exit 3
            next }
    { if (index(ANCS, " " $h["ancestry"] " ") == 0) next
      if (WANT !~ /^ *$/ && index(WANT, " " $h["contrast"] " ") == 0) next
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
        $h["ancestry"], $h["contrast"], $h["case_arm"], $h["ctrl_arm"],
        $h["n_case"], $h["n_ctrl"], $h["case_pct_amppd"], $h["ctrl_pct_amppd"],
        $h["delta_amppd"], $h["confound_tag"], $h["max_callset_delta"],
        $h["worst_callset"], $h["callset_one_sided"], $h["viable_ge100"] }
' "$CONTRASTS_CSV" > "$CROWS" || { echo "ERROR: $CONTRASTS_CSV is missing a required column (above). If it predates 2026-08-21 it has no callset-skew columns — rerun analysis_grain.py." >&2; exit 3; }

NROWS=$(wc -l < "$CROWS")
[[ "$NROWS" -gt 0 ]] || { echo "ERROR: no contrasts.csv rows matched ANCS='${ANCS}' CONTRASTS='${CONTRASTS}'" >&2; exit 1; }
echo "contrasts.csv rows selected: ${NROWS}"

# ── age: report what the covar files ACTUALLY carry, so the log cannot claim a term they lack ──
AGE_NOTE="absent"
for A in $ANCS; do
    if [[ -f "${COVAR_SRC}/covar_${A}.txt" ]]; then
        head -1 "${COVAR_SRC}/covar_${A}.txt" | grep -qw "AGE" && AGE_NOTE="present" || AGE_NOTE="absent"
        break
    fi
done
if [[ "$AGE_NOTE" == "present" ]]; then
    echo "age covariate: AGE column present in §13's covar files — included"
else
    echo "age covariate: no AGE column in §13's covar files — running WITHOUT age."
    echo "  Age is the dominant confounder for AD and PD. Known limitation (HANDOFF issue 9),"
    echo "  blocked on the phenotype track emitting a harmonized age: AMP-AD gives age-at-death,"
    echo "  AMP-PD age-at-baseline, and they are not the same variable. Interpret accordingly."
fi

# summary header
echo "ancestry,contrast,case_arm,ctrl_arm,n_case,n_ctrl,case_pct_amppd,ctrl_pct_amppd,delta_amppd,confound_tag,max_callset_delta,worst_callset,callset_one_sided,viable_ge100,ran,lambda_gc,lambda_1000,n_diffmiss_excluded,n_beta_outlier,n_gwsig,sumstats_file" > "$SUMMARY"

printf "%-5s %-16s %7s %7s %8s %-13s %-7s %-12s %8s %8s %9s %6s\n" "ANC" "contrast" "n_case" "n_ctrl" "delta" "confound" "viable" "status" "lam_gc" "lam_1k" "diffmiss" "hits"

# ── per-stratum staleness gate, once per ancestry ──
# §13's covar must post-date the PCs it claims to carry: an older covar means the grain was built on
# PCs that no longer exist.
declare -A ANC_OK=()
for ANC in $ANCS; do
    QC=${QC_DIR}/cohort_${ANC}_qc
    COVAR=${COVAR_SRC}/covar_${ANC}.txt
    EIG=${QC_DIR}/cohort_${ANC}_pca.eigenvec
    if [[ ! -f "${QC}.bed" ]]; then echo "  skip ${ANC}: no ${QC}.bed"; continue; fi
    if [[ ! -f "$COVAR" ]];     then echo "  skip ${ANC}: no ${COVAR} (§13 wrote no covar for it)"; continue; fi
    if [[ -f "$EIG" && "$COVAR" -ot "$EIG" ]]; then
        echo "ERROR: ${COVAR} is OLDER than ${EIG}." >&2
        echo "  Step 6 has run since the grain was built, so these covariates are the wrong PCs." >&2
        echo "  Rerun analysis_grain.py, then resubmit. (Refusing rather than warning: this is the" >&2
        echo "  failure mode that produced a GWAS on stale covariates before.)" >&2
        exit 4
    fi
    ANC_OK[$ANC]=1
done
[[ ${#ANC_OK[@]} -gt 0 ]] || { echo "ERROR: no stratum had both a fileset and a §13 covar file" >&2; exit 1; }

while IFS=$'\t' read -r ANC TAG CASE_ARM CTRL_ARM NCASE NCTRL CASE_PCT CTRL_PCT DELTA TAGC CSDELTA CSWORST CSONESIDED VIABLE_N; do
        [[ -n "${ANC_OK[$ANC]:-}" ]] || continue
        QC=${QC_DIR}/cohort_${ANC}_qc
        COVAR=${COVAR_SRC}/covar_${ANC}.txt
        PHENO=${PHENO_SRC}/pheno_${ANC}_${TAG}.txt

        # Arm sizes and tags come from contrasts.csv, off the same rows §13 wrote the pheno file
        # from. Re-deriving them here could only introduce disagreement.
        VIABLE="no"; [[ "$VIABLE_N" == "1" ]] && VIABLE="yes"

        # ── a callset present in one arm and absent from the other ──
        # ANNOTATION, NOT A WARNING. One-sidedness alone fires on 5 of 5 real contrasts and predicts
        # nothing; one-sidedness of a SPARSE callset does. n_diffmiss_excluded on this same line is
        # the signal — this string only names which callset explains it. within_cohort rows only.
        # METHODS.md §10 carries the measurement and the two rejected designs.
        SKEW=""
        if [[ "$TAGC" == "within_cohort" && -n "$CSONESIDED" \
              && "$CSONESIDED" != "none" && "$CSONESIDED" != "NA" ]]; then
            SKEW="  [one-sided: ${CSONESIDED} — read n_diffmiss beside it]"
        fi

        if [[ ! -f "$PHENO" ]]; then
            printf "%-5s %-16s %7s %7s %8s %-13s %-7s %-12s %8s %8s %9s %6s%s\n" \
                "$ANC" "$TAG" "$NCASE" "$NCTRL" "$DELTA" "$TAGC" "$VIABLE" "no_pheno_file" "NA" "NA" "NA" "NA" "$SKEW"
            echo "${ANC},${TAG},${CASE_ARM},${CTRL_ARM},${NCASE},${NCTRL},${CASE_PCT},${CTRL_PCT},${DELTA},${TAGC},${CSDELTA},${CSWORST},${CSONESIDED},${VIABLE},no_pheno_file,NA,NA,NA,NA,NA," >> "$SUMMARY"
            continue
        fi

        # ── the pheno file must describe THIS fileset ──
        # Catches a pheno directory from a different step-6 generation, which mtime alone can miss:
        # a covar can be newer than the eigenvec and still list samples this .fam does not contain.
        NOTINFAM=$(awk 'NR==FNR{f[$2]=1;next} FNR>1 && $2!="" && !($2 in f){n++} END{print n+0}' \
                       "${QC}.fam" "$PHENO")
        if [[ "$NOTINFAM" -gt 0 ]]; then
            echo "ERROR: ${PHENO} lists ${NOTINFAM} IIDs absent from ${QC}.fam." >&2
            echo "  That pheno file was built against a different sample set. Rerun analysis_grain.py." >&2
            exit 5
        fi

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
                    # BOTH tests must fail a variant: chi-square scales with n, the raw
                    # difference ignores it. Either alone means something different per stratum.
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
                    # ── |BETA| sanity pass: ADD-only, outlier-free COPY; originals untouched ──
                    # Unconverged Firth fits return implausible effect sizes — numerical failures,
                    # not findings.
                    NBETA=$(awk -v FILT="${OUT}.filtered.tsv" -v BMAX="$BETA_MAX" '
                        NR==1 { for(i=1;i<=NF;i++) h[$i]=i; bc=h["BETA"]; tc=h["TEST"]
                                print > FILT; next }
                        { if(tc && $tc!="ADD") next
                          v=$bc; if(v=="NA"||v=="") { print > FILT; next }
                          if(v>BMAX || v<-BMAX) { n++; next }
                          print > FILT }
                        END { print n+0 }' "$SUMFILE")

                    # ── hits from the |BETA|-filtered ADD set: an unconverged Firth fit can ──
                    # ── carry a tiny p-value, and that is not an association.                 ──
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

        printf "%-5s %-16s %7s %7s %8s %-13s %-7s %-12s %8s %8s %9s %6s%s\n" \
            "$ANC" "$TAG" "$NCASE" "$NCTRL" "$DELTA" "$TAGC" "$VIABLE" "$RAN" "$LAM_GC" "$LAM_1K" "$NDIFF" "$NHIT" "$SKEW"
        echo "${ANC},${TAG},${CASE_ARM},${CTRL_ARM},${NCASE},${NCTRL},${CASE_PCT},${CTRL_PCT},${DELTA},${TAGC},${CSDELTA},${CSWORST},${CSONESIDED},${VIABLE},${RAN},${LAM_GC},${LAM_1K},${NDIFF},${NBETA},${NHIT},${SUMFILE##*/}" >> "$SUMMARY"
done < "$CROWS"

echo "=========================================="
echo "GWAS complete: $(date)"
echo "Summary: ${SUMMARY}"
echo "Sumstats: ${OUT_DIR}/gwas_<ANC>_<CASE>_vs_<CTRL>.pheno.glm.logistic.hybrid"
echo "  + .filtered.tsv (ADD only, |BETA|<=${BETA_MAX})   + .hits.tsv (P<${GWSIG}, from the filtered set)"
echo "Viable (>=100/arm) rows  [ANC contrast (case/ctrl) confound  lambda_gc/lambda_1000]:"
# By header name, not position (rule 4) — this block pointed at the wrong columns the moment two
# were inserted before viable_ge100.
awk -F, '
    NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
    $h["viable_ge100"]=="yes" {
        skew=""
        if ($h["confound_tag"]=="within_cohort" && $h["max_callset_delta"]!="NA" \
            && $h["max_callset_delta"]+0 >= 20)
            skew=sprintf("  <<< callset-skewed: %s d%spp", $h["worst_callset"], $h["max_callset_delta"])
        printf "  %-4s %-30s (%s/%s)  %-13s  lam=%s / l1000=%s%s\n",
            $h["ancestry"], $h["contrast"], $h["n_case"], $h["n_ctrl"],
            $h["confound_tag"], $h["lambda_gc"], $h["lambda_1000"], skew }' "$SUMMARY"
echo
# ── within_cohort contrasts ranked by how much the diffmiss filter had to remove ──
# Threshold-free on purpose: a two-order-of-magnitude outlier is obvious side by side, and no
# percentage-point rule catches it (BR is 2.1pp of its arm). METHODS.md §10.
echo "within_cohort rows by variants removed for differential missingness (descending):"
echo "  a large value here means the arms differed TECHNICALLY despite the within_cohort tag."
awk -F, '
    NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
    $h["confound_tag"]=="within_cohort" && $h["viable_ge100"]=="yes" {
        printf "  %10s  %-4s %-30s one-sided=%s\n",
            $h["n_diffmiss_excluded"], $h["ancestry"], $h["contrast"], $h["callset_one_sided"] }' \
    "$SUMMARY" | sort -k1 -gr
echo "See HANDOFF known issue 10 — the tag measures PROGRAM, not callset."
