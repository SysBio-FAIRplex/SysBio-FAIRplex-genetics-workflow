#!/bin/bash
#
# DIAGNOSTIC — why did every sample come back CAH?  READ-ONLY. Runs in seconds on a login node.
#
#   bash scripts/diag_cah.sh [dataset] [control_dataset]
#   bash scripts/diag_cah.sh br_dsnwgs wgs_harm      # the default
#
# Nothing here submits a job, writes to a callset directory, or reruns genotools. It reads
# artifacts job 27211436 already left on disk.
#
# WHY THESE FIVE THINGS, from genotools/ancestry.py (v1.3.6 and current main are byte-identical
# for this file, so the source below IS what ran):
#
#   ancestry.py:357  mean_imp = SimpleImputer(strategy='mean').fit(X_train)   # fit on the REF panel
#   ancestry.py:421  geno_snps = mean_imp.transform(geno_snps)                # study NAs -> REF mean
#   ancestry.py:470  data = (data - train_mean) / sd                          # ...which is then ~0
#
# A study sample whose matrix is mostly NA is imputed to the reference mean, standardized to the
# origin, and the origin IS the 'ALL' centroid — which is the CAH rule. So "all CAH" is a SINK,
# not a diagnosis: genuinely missing calls, a collapsed SNP intersection, and an allele-orientation
# mismatch all land there. These tests separate them.
#
# NOTE ON A PREVIOUS MEASUREMENT: PROJECT_LOG's "99.8% position overlap" was taken on
# ${PGEN}_filtered_bed. GenoTools does not use that file — ancestry.py:92 first runs
# `--geno 0.1 --max-alleles 2` and intersects against the PRUNED output. Test 2 measures the
# overlap that genotools actually computed. The pruned bed/bim/fam are deleted by clean_up()
# (ancestry.py:255), so its log line is the only surviving record.

set -uo pipefail

BUNDLE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${BUNDLE}/config.sh"

DS="${1:-br_dsnwgs}"
CTRL="${2:-wgs_harm}"

gt_dir() {
    case "$1" in
        br_dsnwgs) echo "${DIR_BR}/genotools" ;;
        wgs_harm)  echo "${DIR_WGS}/genotools" ;;
        divco_hs)  echo "${DIR_DC}/genotools" ;;
        wb_dwgs)   echo "${DIR_WB}/genotools" ;;
        *) echo "" ;;
    esac
}

GT="$(gt_dir "${DS}")"
GTC="$(gt_dir "${CTRL}")"
OUT="${GT}/FILTERED.${DS}"
OUTC="${GTC}/FILTERED.${CTRL}"

[[ -d "${GT}" ]] || { echo "FATAL: no genotools dir for ${DS}: ${GT}"; exit 1; }

hdr() { echo; echo "=============================================================="; echo "$1"; echo "=============================================================="; }

# ── the ancestry step's out_path is NOT ${OUT} ───────────────────────────────
# pipeline.py:86-97 hands the AncestryPredictions object its own out_path, and every ancestry
# artifact is named from THAT, not from the pipeline prefix. ancestry.py:653 writes the labels as
# f'{out_path}_umap_linearsvc_predicted_labels.txt', so the labels file is the one artifact that
# reveals the prefix. Derive it by glob rather than hardcoding the infix — an earlier version of
# this script assumed ${OUT} and reported every ancestry file as MISSING.
#
# Note pipeline.py:89 only redirects out_path into a TemporaryDirectory when --full_output is
# ABSENT. Step 1 passes --full_output, so these files persist. If they are missing even under the
# globbed prefix, that is a real finding: ancestry did not write them, and any labels file present
# is left over from an earlier attempt. Hence the mtimes below.
anc_prefix() {
    local dir="$1" f
    for f in "${dir}"/*_umap_linearsvc_predicted_labels.txt; do
        [[ -e "${f}" ]] || return 1
        echo "${f%_umap_linearsvc_predicted_labels.txt}"
        return 0
    done
}
ANC="$(anc_prefix "${GT}" || true)"
ANCC="$(anc_prefix "${GTC}" || true)"
[[ -n "${ANC}" ]] || { echo "NOTE: no predicted_labels file in ${GT} — falling back to ${OUT}"; ANC="${OUT}"; }

stamp() {  # $1 = path, $2 = label — size + mtime, so stale artifacts are visible
    if [[ -f "$1" ]]; then
        printf '    %-52s %8s  %s\n' "$2" "$(du -h "$1" 2>/dev/null | cut -f1)" "$(date -r "$1" 2>/dev/null)"
    else
        printf '    %-52s %8s\n' "$2" "MISSING"
    fi
}

# ── TEST 0 — what is actually in the directory, and when was it written? ─────
hdr "TEST 0 — genotools output directory"
echo "  dir            : ${GT}"
echo "  ancestry prefix: ${ANC##*/}"
echo
ls -la "${GT}" 2>/dev/null | head -40 | sed 's/^/    /'
echo
echo "  READ: job 27211436 finished 2026-08-12. Any ancestry artifact whose mtime is NOT that"
echo "        date was written by an earlier attempt — including, possibly, the labels file"
echo "        that reports 97 CAH. Check that before trusting the CAH result at all."

# ── TEST 1 — did the plink2 module move under us? ────────────────────────────
# Every plink log carries its own version banner, so July's logs record July's plink2. This
# compares the binary that produced the control callset (worked) against the one that produced
# this run (all CAH), with no need to remember what was loaded at the time.
hdr "TEST 1 — plink2 version, this run vs the control callset"
for pair in "${DS}:${OUT}" "${CTRL}:${OUTC}"; do
    name="${pair%%:*}"; pre="${pair#*:}"
    log="${pre}_all_logs.log"
    if [[ -f "${log}" ]]; then
        ver="$(grep -m1 -E '^PLINK v' "${log}" 2>/dev/null)"
        printf '  %-12s %s\n' "${name}" "${ver:-<no PLINK banner in log>}"
        printf '  %-12s log mtime: %s\n' "" "$(date -r "${log}" 2>/dev/null)"
    else
        printf '  %-12s MISSING: %s\n' "${name}" "${log}"
    fi
done
echo
echo "  currently loaded:"
command -v plink2 >/dev/null 2>&1 && plink2 --version 2>&1 | head -1 | sed 's/^/    /' \
    || echo "    (plink2 not on PATH — 'module load plink/6-alpha' first)"
echo
# GenoTools does NOT use the module — dependencies.check_plink2() resolves a plink2 it ships
# with itself, so the banner in the logs above is genotools' own binary and the module version
# is irrelevant to anything inside ancestry.py. But steps 1-3 of 01_genotools.sh DO use the
# module to build ${PGEN}_filtered_bed, which is genotools' INPUT. That path is still exposed.
echo "  --- module-built inputs (01_genotools.sh steps 1-3, these DO use the module) ---"
for pair in "${DS}:${RAW_BR}" "${CTRL}:${RAW_WGS}"; do
    name="${pair%%:*}"; stem="${pair#*:}"
    for suf in _filtered _filtered_bed; do
        log="${stem}${suf}.log"
        if [[ -f "${log}" ]]; then
            printf '    %-12s %-16s %s | %s\n' "${name}" "${suf}" \
                "$(grep -m1 -E '^PLINK v' "${log}" 2>/dev/null || echo '<no banner>')" \
                "$(date -r "${log}" 2>/dev/null)"
        else
            printf '    %-12s %-16s MISSING %s\n' "${name}" "${suf}" "${log}"
        fi
    done
done
echo
echo "  READ: the two genotools banners matching means the ancestry path used the SAME binary"
echo "        in July and August — the module version above is genotools' business only if it"
echo "        differs there. If the genotools banners match but the module-built input banners"
echo "        do NOT, the update moved plink/6-alpha and changed the INPUT rather than the"
echo "        projection, and ${DS}'s pgen is no longer built the way wgs_harm's was."

# ── TEST 2 — what survived --geno 0.1, and what intersected the panel? ───────
hdr "TEST 2 — variant counts through the genotools pipeline"
LOG="${OUT}_all_logs.log"
if [[ -f "${LOG}" ]]; then
    echo "  --- lines recording variant counts (${LOG##*/}) ---"
    grep -nE 'variants loaded|variants remaining|--geno:|--max-alleles:|--extract:|pass filters|Error|Warning' \
        "${LOG}" | head -60 | sed 's/^/    /'
else
    echo "  MISSING: ${LOG}"
fi
echo
echo "  --- SNP set sizes ---"
REF_BIM="${REF_PANEL}.bim"
MODEL_SNPS="${ANC}_umap_linearsvc_ancestry_model.common_snps"   # ref side  (ancestry.py:101-104)
GENO_BIM="${ANC}_common_snps.bim"                               # study side (ancestry.py:197)
for f in "${REF_BIM}" "${MODEL_SNPS}" "${GENO_BIM}"; do
    if [[ -f "${f}" ]]; then printf '    %10s  %s\n' "$(wc -l < "${f}" | tr -d ' ')" "${f##*/}"
    else printf '    %10s  %s\n' "MISSING" "${f##*/}"; fi
done
echo
echo "  READ: the panel is ~209,517 variants. If the two common_snps counts are close to that,"
echo "        the intersection is healthy and the fault is downstream (go to TEST 4)."
echo "        If they are small, the intersection collapsed — and --geno 0.1 above will say"
echo "        whether pruning caused it. get_common_snps (utils.py:404) prints only"
echo "        'Getting Common SNPs' and never reports its overlap, so an empty intersection"
echo "        exits 0 with plausible-looking output. That is the silence being measured here."

# ── TEST 3 — the decisive one: NA fraction in the matrix that was projected ──
# clean_up() (ancestry.py:255) removes the pruned intermediates but NOT _common_snps.raw, so
# the exact matrix passed to mean_imp.transform() is still on disk. This measures it directly
# instead of inferring it from --missing on a different file.
hdr "TEST 3 — NA fraction in the matrix actually handed to the projection"
RAW="${ANC}_common_snps.raw"
if [[ -f "${RAW}" ]]; then
    echo "  ${RAW##*/}  ($(du -h "${RAW}" | cut -f1))"
    awk '
        NR==1 { nsnp = NF - 6; next }
        {
            na = 0
            for (i = 7; i <= NF; i++) if ($i == "NA" || $i == "nan" || $i == "") na++
            tot_na += na; tot_cell += (NF - 6); nsamp++
            f = (NF > 6) ? na / (NF - 6) : 0
            if (nsamp == 1 || f < min) min = f
            if (nsamp == 1 || f > max) max = f
        }
        END {
            printf "    SNP columns          : %d\n", nsnp
            printf "    samples              : %d\n", nsamp
            printf "    overall NA fraction  : %.4f\n", (tot_cell ? tot_na / tot_cell : 0)
            printf "    per-sample NA range  : %.4f - %.4f\n", min, max
        }' "${RAW}"
else
    echo "  MISSING: ${RAW}"
    echo "  (if genotools was run with a model_path or in a container this may not be written)"
fi
echo
echo "  READ: ~0.75 and up => hom-ref is encoded as missing, mean-imputation to the reference"
echo "        mean sends every sample to the origin, and CAH follows mechanically."
echo "        ~0.02 => the calls are fine. The matrix is dense, so the collapse (if TEST 5"
echo "        shows one) is orientation or scaling, not missingness — suspect the"
echo "        --alt1-allele switch at ancestry.py:216-220."

# ── TEST 4 — did the samples actually collapse to the origin? ────────────────
hdr "TEST 4 — projected study PCs vs reference PCs"
pcstat() {  # $1 = file, $2 = label.  Layout: FID IID PC1..PC50 [label]
    if [[ ! -f "$1" ]]; then printf '  %-24s MISSING: %s\n' "$2" "${1##*/}"; return; fi
    awk -v tag="$2" '
        NR == 1 { next }
        { n++; for (i = 3; i <= 12; i++) { s[i] += $i; ss[i] += $i * $i } }
        END {
            printf "  %s  (n=%d)\n", tag, n
            for (i = 3; i <= 12; i++) {
                m = s[i] / n; v = ss[i] / n - m * m; if (v < 0) v = 0
                printf "    PC%-2d  mean %11.5f   sd %10.5f\n", i - 2, m, sqrt(v)
            }
        }' "$1"
}
pcstat "${ANC}_projected_new_pca.txt" "projected (study samples)"
echo
pcstat "${ANC}_labeled_ref_pca.txt"   "reference panel"
echo
echo "  --- mtimes (is any of this from an earlier attempt?) ---"
stamp "${ANC}_projected_new_pca.txt"                "projected_new_pca.txt"
stamp "${ANC}_labeled_ref_pca.txt"                  "labeled_ref_pca.txt"
stamp "${ANC}_pca_eigenvalues.txt"                  "pca_eigenvalues.txt"
stamp "${ANC}_umap_linearsvc_predicted_labels.txt"  "predicted_labels.txt"
stamp "${RAW}"                                      "common_snps.raw"
echo
echo "  --- top PCs' explained variance (ancestry.py:406) ---"
if [[ -f "${ANC}_pca_eigenvalues.txt" ]]; then head -11 "${ANC}_pca_eigenvalues.txt" | sed 's/^/    /'
else echo "    MISSING ${ANC##*/}_pca_eigenvalues.txt"; fi
echo
echo "  READ: study sd near zero with means near zero, while the reference spreads normally,"
echo "        IS the collapse — every sample sitting on the 'ALL' centroid. Confirms the sink."
echo "        Study PCs with REAL spread that are still all labelled CAH means something else"
echo "        entirely: the samples are somewhere real but off the panel's manifold, and the"
echo "        suspect becomes scaling/orientation rather than missing data."

# ── TEST 5 — what the labels actually say ────────────────────────────────────
hdr "TEST 5 — label counts"
for pair in "${DS}:${ANC}" "${CTRL}:${ANCC}"; do
    name="${pair%%:*}"; pre="${pair#*:}"
    lbl="${pre}_umap_linearsvc_predicted_labels.txt"
    if [[ -n "${pre}" && -f "${lbl}" ]]; then
        echo "  ${name}:"
        awk 'NR>1 {print $NF}' "${lbl}" | sort | uniq -c | sort -rn | sed 's/^/    /'
    else
        echo "  ${name}: MISSING ${lbl##*/}"
    fi
done
echo
echo "  READ: the control callset's labels are from JULY, before the update. They are the"
echo "        record of what a working run produced, not a fresh result."

hdr "DONE — read-only, nothing was modified"
