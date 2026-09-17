#!/bin/bash
#SBATCH --job-name=amppd_release
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --partition=norm
#
# STEP 9 (follow-up, not tier 1) — subset the AMP-PD donors out of the final QC'd association
# set and publish them to a GCS bucket.
#
# This runs AFTER step 8 and changes nothing upstream. It is a release step: it reads step 6
# stage C's output in place and writes a new staging tree. No pipeline artifact is modified.
#
# WHAT "the final merged and QC'd dataset" MEANS HERE
#   data/merged/by_ancestry_qc/cohort_<ANC>_qc.{bed,bim,fam} — step 6 STAGE C. Post per-ancestry
#   QC (geno<=0.05, maf>=0.01, hwe>=1e-6, autosomes) AND post the AF-concordance exclusion list.
#   This is the exact variant set step 7 tested, which is why it is the thing worth shipping.
#   It is deliberately NOT cohort_merged (pre-QC) and NOT unfiltered/ (the baseline twin).
#
#   The release is therefore PER STRATUM, not one file — because the QC is per stratum. Each
#   cohort_<ANC>_qc has its own variant set; concatenating them would invent a fileset that no
#   GWAS ever ran on. See METHODS.md §5-6.
#
# WHAT IS NOT RE-APPLIED, AND WHY
#   Subsetting to AMP-PD drops ~2/3 of the samples, so some retained variants become monomorphic
#   or low-MAF within the subset. This step does NOT re-apply --maf/--geno/--hwe by default.
#   Re-filtering would produce a variant set that is not the one step 7 consumed, and the
#   released pgen would then silently disagree with the released sumstats. Set SUBSET_MAF /
#   SUBSET_GENO to opt in; the manifest and the release README record whichever was used.
#
# WHO COUNTS AS AMP-PD
#   source_callset in {wb_dwgs, br_dsnwgs}, read BY HEADER NAME from step 6 stage E's
#   retained_samples_manifest.csv (config.sh RETAINED_MANIFEST). That file is the single existing
#   definition of per-sample callset membership — this step does not re-derive it from the four
#   genotools label files, which would be a second implementation of the same concept.
#
# MODES — the build and the push run on DIFFERENT HOSTS.
#   MODE=build (default)  sbatch on biowulf. Subsets, hashes, writes MANIFEST.tsv + README.md.
#                           ./submit.sh scripts/09_amppd_release.sh
#   MODE=push             run DIRECTLY on helix.nih.gov, not under sbatch. Biowulf compute nodes
#                         have no general outbound network; helix is the transfer host. Refuses
#                         to run inside a compute-node allocation unless ALLOW_COMPUTE_PUSH=1.
#                           MODE=push GCS_DEST=gs://<bucket>/<prefix> bash scripts/09_amppd_release.sh
#
# THE PUSH IS OPT-IN AND GATED. MODE=build never contacts the network. MODE=push refuses to
# start unless GCS_DEST is set, and refuses to upload unless the destination bucket has uniform
# bucket-level access ON, public access prevention ENFORCED, and no allUsers /
# allAuthenticatedUsers IAM binding. This is individual-level genotype data from two
# controlled-access programs; redistribution is governed by the AMP-PD DUA. Widening bucket
# access is a separate decision that this script will not make for you.
#
# KNOBS
#   GCS_DEST            gs://bucket/prefix — required for MODE=push. No bucket is hardcoded.
#   STAGE_DIR           where the release tree is built (default data/merged/release_amppd)
#   AMPPD_CALLSETS      "wb_dwgs br_dsnwgs"
#   OUT_FMT             pgen (default) | bed | vcf
#   SUBSET_MAF          re-apply --maf inside the subset (default: unset, do not re-filter)
#   SUBSET_GENO         re-apply --geno inside the subset (default: unset)
#   OVERWRITE=1         allow MODE=build to replace a staging tree, or MODE=push to overwrite
#                       objects already under GCS_DEST. Both refuse without it.
#   ALLOW_NO_HASH=1     proceed when local CRC32C cannot be computed (see stage C).
#   MOD_GCLOUD          module providing gcloud (default google-cloud-sdk); skipped if gcloud
#                       is already on PATH.
#
# GUARDRAIL: this step's DELIVERABLE is sample IDs and genotypes. Run it yourself. The AI must
# not run it, and must not read the staging tree's .psam/.fam. All output below is aggregate.
set -o pipefail

BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

MODE=${MODE:-build}
QC_DIR=${MERGED_DIR}/by_ancestry_qc
STAGE_DIR=${STAGE_DIR:-${MERGED_DIR}/release_amppd}
AMPPD_CALLSETS=${AMPPD_CALLSETS:-"wb_dwgs br_dsnwgs"}
OUT_FMT=${OUT_FMT:-pgen}
THREADS=${SLURM_CPUS_PER_TASK:-16}
MANIFEST_TSV=${STAGE_DIR}/MANIFEST.tsv
RELEASE_README=${STAGE_DIR}/README.md
PROV=${STAGE_DIR}/.provenance

banner() { echo "=========================================="; echo "$@"; echo "=========================================="; }

# ═════════════════════════════════════════════════════════════════════════════
# MODE=build — stages A, B, C. No network.
# ═════════════════════════════════════════════════════════════════════════════
if [[ "$MODE" == "build" ]]; then

banner "Step 9 BUILD — AMP-PD subset of the QC'd association set"
echo "Job ${SLURM_JOB_ID:-local} on ${SLURMD_NODENAME:-$(hostname)}   start $(date)"
echo "Source (step 6 stage C): ${QC_DIR}/cohort_<ANC>_qc"
echo "Sample source:           ${RETAINED_MANIFEST}"
echo "AMP-PD callsets:         ${AMPPD_CALLSETS}"
echo "Staging:                 ${STAGE_DIR}"
echo "Format:                  ${OUT_FMT}"
echo "Re-filter inside subset: maf=${SUBSET_MAF:-<none>} geno=${SUBSET_GENO:-<none>}"
echo ""

# ── preconditions. Each one fails loudly; none of them is allowed to be silently skipped. ──
[[ -f "${RETAINED_MANIFEST}" ]] || {
    echo "ERROR: retained_samples_manifest.csv not found: ${RETAINED_MANIFEST}" >&2
    echo "  It is step 6 stage E's output. Run step 6, or point RETAINED_MANIFEST at it." >&2
    echo "  Do NOT substitute relatedness/retained_manifest.csv — that file has no" >&2
    echo "  source_callset column, so AMP-PD membership cannot be read from it." >&2
    exit 1; }

compgen -G "${QC_DIR}/cohort_*_qc.bed" >/dev/null || {
    echo "ERROR: no cohort_<ANC>_qc.bed under ${QC_DIR} — step 6 stage C has not run." >&2
    exit 1; }

if [[ -d "${STAGE_DIR}" && -z "${OVERWRITE}" ]]; then
    echo "ERROR: staging tree already exists: ${STAGE_DIR}" >&2
    echo "  Refusing to overwrite a release that may already have been pushed." >&2
    echo "  OVERWRITE=1 to rebuild it, or set STAGE_DIR elsewhere." >&2
    exit 1
fi
rm -rf "${STAGE_DIR}"
mkdir -p "${STAGE_DIR}/keep"

module load "${MOD_PLINK2}"

# ─────────────────────────────────────────────────────────────────────────────
# STAGE A — select the AMP-PD IIDs, by header name.
#
# Rule 4 (verify by name, not by position): retained_samples_manifest.csv is written by
# ancestry_qc_manifest.py as IID,source_callset,ancestry,call_rate,dup_cluster_id,PC1.. — but
# a PC column count change has already broken one positional reader in this project, so the
# column indices are resolved from the header row rather than assumed.
# ─────────────────────────────────────────────────────────────────────────────
echo "--- STAGE A: select AMP-PD samples ---"

AMPPD_IIDS=${STAGE_DIR}/keep/amppd_iids.txt
awk -F, -v want="${AMPPD_CALLSETS}" '
    NR==1 {
        for (i=1; i<=NF; i++) { h[$i]=i }
        if (!("IID" in h))            { print "ERROR: no IID column in manifest"            > "/dev/stderr"; exit 9 }
        if (!("source_callset" in h)) { print "ERROR: no source_callset column in manifest" > "/dev/stderr"; exit 9 }
        if (!("ancestry" in h))       { print "ERROR: no ancestry column in manifest"       > "/dev/stderr"; exit 9 }
        iid=h["IID"]; sc=h["source_callset"]; anc=h["ancestry"]
        n=split(want, w, /[ ,]+/); for (k=1; k<=n; k++) keep[w[k]]=1
        next
    }
    $iid=="" { next }
    ($sc in keep) { print $iid "\t" $anc "\t" $sc }
' "${RETAINED_MANIFEST}" > "${AMPPD_IIDS}"
RC=$?
[[ "$RC" -eq 0 ]] || { echo "ERROR: could not read ${RETAINED_MANIFEST} (exit ${RC})" >&2; exit "$RC"; }

N_AMPPD=$(wc -l < "${AMPPD_IIDS}" | tr -d ' ')
if [[ "$N_AMPPD" -eq 0 ]]; then
    echo "ERROR: zero samples matched source_callset in {${AMPPD_CALLSETS}}." >&2
    echo "  Observed source_callset values in the manifest:" >&2
    awk -F, 'NR==1 { for (i=1;i<=NF;i++) if ($i=="source_callset") c=i; next } { print "    " $c }' \
        "${RETAINED_MANIFEST}" | sort | uniq -c | sort -rn >&2
    echo "  This is a NAME mismatch, not an empty cohort. Fix AMPPD_CALLSETS; do not proceed." >&2
    exit 1
fi

echo "AMP-PD samples in the retained set: ${N_AMPPD}"
echo "  by callset:"
cut -f3 "${AMPPD_IIDS}" | sort | uniq -c | awk '{printf "    %-12s %6d\n", $2, $1}'
echo "  by ancestry:"
cut -f2 "${AMPPD_IIDS}" | sort | uniq -c | sort -k1,1nr | awk '{printf "    %-6s %6d\n", $2, $1}'
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STAGE B — one subset per stratum.
#
# The keep file is built by joining the AMP-PD IIDs against the stratum's OWN .fam rather than
# against the manifest's ancestry column. Two reasons: the .fam supplies the FID that plink2
# --keep wants, and the join makes the sample accounting exact — a sample can only be counted
# as shipped if it is actually in the fileset being subset. The .fam is read positionally
# because it is a fixed-format file with no header; every CSV above is read by name.
# ─────────────────────────────────────────────────────────────────────────────
echo "--- STAGE B: subset per stratum -> ${STAGE_DIR} ---"
printf "%-6s %9s %9s %12s  %s\n" "ANC" "in_strat" "amppd" "variants" "note"

case "${OUT_FMT}" in
    pgen) FMT_ARG="--make-pgen"; FMT_EXT="pgen pvar psam" ;;
    bed)  FMT_ARG="--make-bed";  FMT_EXT="bed bim fam"    ;;
    vcf)  FMT_ARG="--export vcf bgz"; FMT_EXT="vcf.gz"    ;;
    *)    echo "ERROR: OUT_FMT must be pgen, bed or vcf (got '${OUT_FMT}')" >&2; exit 2 ;;
esac

FILTER_ARGS=""
[[ -n "${SUBSET_MAF}"  ]] && FILTER_ARGS="${FILTER_ARGS} --maf ${SUBSET_MAF}"
[[ -n "${SUBSET_GENO}" ]] && FILTER_ARGS="${FILTER_ARGS} --geno ${SUBSET_GENO}"

N_SHIPPED=0
N_STRATA=0
declare -A SHIPPED_BY_ANC
FAILED=""

for BED in "${QC_DIR}"/cohort_*_qc.bed; do
    ANC=$(basename "$BED"); ANC=${ANC#cohort_}; ANC=${ANC%_qc.bed}
    STEM=${QC_DIR}/cohort_${ANC}_qc
    N_STRAT=$(wc -l < "${STEM}.fam" | tr -d ' ')

    KEEP=${STAGE_DIR}/keep/keep_amppd_${ANC}.txt
    awk 'NR==FNR { want[$1]=1; next } ($2 in want) { print $1 "\t" $2 }' \
        "${AMPPD_IIDS}" "${STEM}.fam" > "${KEEP}"
    N_KEEP=$(wc -l < "${KEEP}" | tr -d ' ')

    if [[ "$N_KEEP" -eq 0 ]]; then
        printf "%-6s %9s %9s %12s  %s\n" "$ANC" "$N_STRAT" "0" "-" "no AMP-PD samples — skipped"
        continue
    fi

    OUT=${STAGE_DIR}/amppd_${ANC}
    if ! plink2 --bfile "$STEM" --keep "$KEEP" ${FILTER_ARGS} \
                $FMT_ARG --threads "$THREADS" --out "$OUT" >"${OUT}.plink.log" 2>&1; then
        printf "%-6s %9s %9s %12s  %s\n" "$ANC" "$N_STRAT" "$N_KEEP" "FAILED" "see ${OUT}.plink.log"
        FAILED="${FAILED} ${ANC}"
        continue
    fi

    case "${OUT_FMT}" in
        pgen) NVAR=$(grep -vc '^#' "${OUT}.pvar") ;;
        bed)  NVAR=$(wc -l < "${OUT}.bim" | tr -d ' ') ;;
        vcf)  NVAR=$(zcat "${OUT}.vcf.gz" | grep -vc '^#') ;;
    esac

    NOTE=""
    if [[ -n "${FILTER_ARGS}" ]]; then
        NVAR_SRC=$(wc -l < "${STEM}.bim" | tr -d ' ')
        NOTE="re-filtered:${FILTER_ARGS# } (was ${NVAR_SRC})"
    fi
    printf "%-6s %9s %9s %12s  %s\n" "$ANC" "$N_STRAT" "$N_KEEP" "$NVAR" "$NOTE"
    N_SHIPPED=$((N_SHIPPED + N_KEEP))
    N_STRATA=$((N_STRATA + 1))
    SHIPPED_BY_ANC[$ANC]=$N_KEEP
done

echo ""
if [[ -n "$FAILED" ]]; then
    echo "ERROR: plink2 failed for stratum/strata:${FAILED}" >&2
    echo "  Refusing to publish a partial release. Fix and rerun with OVERWRITE=1." >&2
    exit 1
fi

# ── the accounting. Rule 3: a release that silently drops samples looks exactly like a
# release that had none to drop, so the reconciliation is mandatory and fatal, not advisory.
echo "--- sample accounting ---"
echo "AMP-PD in retained manifest : ${N_AMPPD}"
echo "AMP-PD written to release   : ${N_SHIPPED}  across ${N_STRATA} strata"
if [[ "$N_SHIPPED" -ne "$N_AMPPD" ]]; then
    MISSING=$((N_AMPPD - N_SHIPPED))
    echo ""
    echo "ERROR: ${MISSING} AMP-PD sample(s) are in the retained manifest but reached no" >&2
    echo "  release fileset. They are listed below by the ancestry the manifest assigns them;" >&2
    echo "  the usual cause is a stratum with <2 samples, for which step 6 stage A wrote no" >&2
    echo "  fileset at all. That is a real hole in the release, not a rounding difference." >&2
    for ANC in $(cut -f2 "${AMPPD_IIDS}" | sort -u); do
        WANT=$(awk -F"\t" -v a="$ANC" '$2==a' "${AMPPD_IIDS}" | wc -l | tr -d ' ')
        GOT=${SHIPPED_BY_ANC[$ANC]:-0}
        [[ "$WANT" -ne "$GOT" ]] && printf "    %-6s manifest %5d  shipped %5d  gap %5d\n" "$ANC" "$WANT" "$GOT" "$((WANT-GOT))" >&2
    done
    echo "  Resolve or accept each gap deliberately, then rerun with OVERWRITE=1." >&2
    exit 1
fi
echo "reconciled: every AMP-PD sample in the manifest is in exactly one released fileset."
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STAGE C — manifest + release README, both DERIVED from what stage B actually wrote.
#
# Neither is transcribed. The 2026-08-24 lesson (review/methods_numbers.py) was that every
# hand-typed number in a document drifts from its artifact invisibly; a release doc is the
# worst place for that, because the reader has no way to check it.
#
# CRC32C rather than MD5: GCS does not compute MD5 for composite (parallel-chunked) uploads,
# and 42 of the 46 objects in the previous sumstats release came back MD5-less for exactly
# that reason. CRC32C is present on every object and is what stage E compares.
# ─────────────────────────────────────────────────────────────────────────────
echo "--- STAGE C: manifest + release README ---"

# gcloud is only needed here for the local hash. If it is absent the manifest still gets
# written, but with no crc32c — and stage E then has nothing but size to verify against, so
# the build says so and exits non-zero unless that is explicitly accepted.
HAVE_HASH=1
if ! command -v gcloud >/dev/null 2>&1; then
    module load "${MOD_GCLOUD:-google-cloud-sdk}" 2>/dev/null || true
fi
command -v gcloud >/dev/null 2>&1 || HAVE_HASH=0

printf 'file\tbytes\tcrc32c\tsha256\n' > "${MANIFEST_TSV}"
for F in "${STAGE_DIR}"/amppd_*; do
    B=$(basename "$F")
    case "$B" in *.plink.log) continue ;; esac
    SZ=$(stat -c %s "$F" 2>/dev/null || stat -f %z "$F")
    CRC="—"
    if [[ "$HAVE_HASH" -eq 1 ]]; then
        CRC=$(gcloud storage hash --skip-md5 "$F" 2>/dev/null \
              | awk -F': *' 'tolower($1) ~ /crc32c/ {print $2; exit}')
        [[ -n "$CRC" ]] || CRC="—"
    fi
    SHA=$(sha256sum "$F" 2>/dev/null | cut -d' ' -f1)
    [[ -n "$SHA" ]] || SHA="—"
    printf '%s\t%s\t%s\t%s\n' "$B" "$SZ" "$CRC" "$SHA" >> "${MANIFEST_TSV}"
done

N_OBJ=$(($(wc -l < "${MANIFEST_TSV}" | tr -d ' ') - 1))
TOT_BYTES=$(awk -F'\t' 'NR>1 {s+=$2} END {print s+0}' "${MANIFEST_TSV}")
TOT_GIB=$(awk -v b="$TOT_BYTES" 'BEGIN {printf "%.2f", b/1073741824}')

if [[ "$HAVE_HASH" -eq 0 ]]; then
    echo ""
    echo "!! CRC32C COULD NOT BE COMPUTED — gcloud is not on PATH and MOD_GCLOUD did not load." >&2
    echo "   The manifest carries sizes and sha256 but no crc32c, so stage E will be able to" >&2
    echo "   verify that every object ARRIVED and is the right SIZE, and nothing about its" >&2
    echo "   CONTENT. That is a weaker check than this release is supposed to get." >&2
    if [[ -z "${ALLOW_NO_HASH}" ]]; then
        echo "   Set MOD_GCLOUD to the right module, or ALLOW_NO_HASH=1 to accept size-only." >&2
        exit 1
    fi
    echo "   ALLOW_NO_HASH=1 — continuing with size-only verification." >&2
fi

# ── provenance, machine-readable, consumed by MODE=push ──
{
    echo "built_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "build_job=${SLURM_JOB_ID:-local}"
    echo "build_host=${SLURMD_NODENAME:-$(hostname)}"
    echo "git_commit=$(cd "${BUNDLE}" && git rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "git_dirty=$(cd "${BUNDLE}" && [[ -n "$(git status --porcelain 2>/dev/null)" ]] && echo yes || echo no)"
    echo "source_dir=${QC_DIR}"
    echo "n_samples=${N_AMPPD}"
    echo "n_strata=${N_STRATA}"
    echo "n_objects=${N_OBJ}"
    echo "total_bytes=${TOT_BYTES}"
    echo "out_fmt=${OUT_FMT}"
    echo "subset_maf=${SUBSET_MAF:-none}"
    echo "subset_geno=${SUBSET_GENO:-none}"
} > "${PROV}"

# ── the release README, generated ──
{
    echo "# AMP-PD subset — AD-vs-PD WGS GWAS association set"
    echo ""
    echo "Generated by \`scripts/09_amppd_release.sh\` on $(date -u +%Y-%m-%d) from the project's"
    echo "final QC'd association set. Nothing in this file is typed by hand."
    echo ""
    echo "## What this is"
    echo ""
    echo "The **AMP-PD donors only** (\`source_callset\` in {$(echo ${AMPPD_CALLSETS} | tr ' ' ',')}),"
    echo "subset out of \`cohort_<ANC>_qc\` — step 6 stage C of the pipeline, which is the exact"
    echo "genotype set the released summary statistics were computed on."
    echo ""
    echo "QC already applied upstream, per ancestry stratum:"
    echo ""
    echo "- autosomes only"
    echo "- \`--geno 0.05\`, \`--maf 0.01\`, \`--hwe 1e-6 keep-fewhet\`"
    echo "- duplicates, relatives and sex-check failures removed (KING; step 5)"
    echo "- the AF-concordance exclusion list applied (step 6a) — cross-callset frequency"
    echo "  outliers, removed because disease is held constant within each (stratum x dx) cell"
    echo ""
    echo "## One fileset per ancestry stratum, on purpose"
    echo ""
    echo "QC is per stratum, so each stratum has its own variant set. These are **not**"
    echo "concatenated into one fileset: doing so would produce a variant set that no GWAS in"
    echo "this study ever ran on."
    echo ""
    echo "## Not re-filtered after subsetting"
    if [[ -z "${FILTER_ARGS}" ]]; then
        echo ""
        echo "Removing the AMP-AD donors leaves some variants monomorphic or low-MAF within this"
        echo "subset. They are **retained**. Re-applying \`--maf\`/\`--geno\` here would yield a"
        echo "variant set that disagrees with the published summary statistics. Apply your own"
        echo "thresholds downstream if you need them."
    else
        echo ""
        echo "**This release WAS re-filtered after subsetting** (\`${FILTER_ARGS# }\`), so its"
        echo "variant set does **not** match the published summary statistics one-for-one."
    fi
    echo ""
    echo "## Contents"
    echo ""
    echo "${N_OBJ} objects, ${TOT_GIB} GiB, ${N_AMPPD} samples across ${N_STRATA} strata."
    echo ""
    echo "| stratum | samples | variants |"
    echo "|---|---|---|"
    for ANC in $(printf '%s\n' "${!SHIPPED_BY_ANC[@]}" | sort); do
        case "${OUT_FMT}" in
            pgen) NV=$(grep -vc '^#' "${STAGE_DIR}/amppd_${ANC}.pvar") ;;
            bed)  NV=$(wc -l < "${STAGE_DIR}/amppd_${ANC}.bim" | tr -d ' ') ;;
            vcf)  NV=$(zcat "${STAGE_DIR}/amppd_${ANC}.vcf.gz" | grep -vc '^#') ;;
        esac
        echo "| \`amppd_${ANC}\` | ${SHIPPED_BY_ANC[$ANC]} | ${NV} |"
    done
    echo ""
    echo "\`MANIFEST.tsv\` carries per-object size, CRC32C and SHA-256. Verify with CRC32C, not"
    echo "MD5 — GCS does not compute MD5 for composite (parallel-chunked) uploads."
    echo ""
    echo "## Known limitations carried in from the study"
    echo ""
    echo "- **BR-DSNWGS contributes to roughly half the association set.** A 97-donor joint call"
    echo "  emits nothing at sites monomorphic in its own donors, so BR-DSNWGS genotypes are"
    echo "  missing at the other half by construction, not by QC failure."
    echo "- **No age covariate.** AMP-AD supplies age at death and AMP-PD age at baseline; the"
    echo "  two were not forced into one column. See \`METHODS.md\` §10."
    echo "- Strata below ~50 samples have no principal components — plink2 will not LD-prune at"
    echo "  those sizes. They are below the study's 100-per-arm viability cut regardless."
    echo ""
    echo "## Access"
    echo ""
    echo "Individual-level genotypes derived from AMP-PD controlled-access data. Redistribution"
    echo "is governed by the AMP-PD data use agreement. This bucket is private by construction:"
    echo "the publishing script refuses to upload unless uniform bucket-level access is on,"
    echo "public access prevention is enforced, and no \`allUsers\` or \`allAuthenticatedUsers\`"
    echo "binding exists. **Widening access is a DUA decision, not an operational one.**"
    echo ""
    echo "## Provenance"
    echo ""
    echo "\`\`\`"
    cat "${PROV}"
    echo "\`\`\`"
} > "${RELEASE_README}"

echo "manifest : ${MANIFEST_TSV}  (${N_OBJ} objects, ${TOT_GIB} GiB)"
echo "readme   : ${RELEASE_README}"
echo ""
banner "Step 9 BUILD complete $(date)"
echo "NOTHING HAS BEEN UPLOADED. The push runs on helix, not here:"
echo ""
echo "  ssh helix.nih.gov"
echo "  cd ${BUNDLE}"
echo "  MODE=push GCS_DEST=gs://<bucket>/<prefix> bash scripts/09_amppd_release.sh"
echo ""
exit 0
fi

# ═════════════════════════════════════════════════════════════════════════════
# MODE=push — stages D and E. Runs on helix.nih.gov.
# ═════════════════════════════════════════════════════════════════════════════
if [[ "$MODE" == "push" ]]; then

banner "Step 9 PUSH — publish the AMP-PD release to GCS"
echo "Host: $(hostname)   start $(date)"

[[ -n "${GCS_DEST}" ]] || {
    echo "ERROR: GCS_DEST is not set, and no bucket is hardcoded in this script." >&2
    echo "  MODE=push GCS_DEST=gs://<bucket>/<prefix> bash scripts/09_amppd_release.sh" >&2
    exit 2; }
[[ "${GCS_DEST}" == gs://* ]] || { echo "ERROR: GCS_DEST must start with gs:// (got '${GCS_DEST}')" >&2; exit 2; }
GCS_DEST=${GCS_DEST%/}
BUCKET=${GCS_DEST#gs://}; BUCKET=${BUCKET%%/*}

[[ -f "${PROV}" && -f "${MANIFEST_TSV}" ]] || {
    echo "ERROR: no completed build at ${STAGE_DIR} (missing MANIFEST.tsv or .provenance)." >&2
    echo "  Run MODE=build under sbatch on biowulf first." >&2
    exit 1; }

# Biowulf compute nodes have no general outbound network. A push attempted from one fails in a
# way that looks like a credentials or bucket problem, which is an expensive misdiagnosis.
if [[ -n "${SLURM_JOB_ID}" && -z "${ALLOW_COMPUTE_PUSH}" ]]; then
    echo "ERROR: MODE=push is running inside SLURM job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}." >&2
    echo "  Biowulf compute nodes have no general outbound network — this upload will fail," >&2
    echo "  and it will fail looking like an auth problem. Run it on helix.nih.gov instead." >&2
    echo "  ALLOW_COMPUTE_PUSH=1 if you know this node has egress." >&2
    exit 1
fi

if ! command -v gcloud >/dev/null 2>&1; then
    module load "${MOD_GCLOUD:-google-cloud-sdk}" 2>/dev/null || true
fi
command -v gcloud >/dev/null 2>&1 || {
    echo "ERROR: gcloud not found on PATH and MOD_GCLOUD='${MOD_GCLOUD:-google-cloud-sdk}' did not load." >&2
    echo "  Find it with: module spider google-cloud-sdk" >&2
    exit 1; }

gcloud auth print-access-token >/dev/null 2>&1 || {
    echo "ERROR: no active gcloud credentials. Run 'gcloud auth login' yourself, then rerun." >&2
    exit 1; }

# Both the bucket preflight and the verify read gcloud's --json output rather than its
# human-facing text. That needs a python3; helix has one, but say so rather than failing
# three stages later inside a pipe whose exit status is easy to lose.
command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 not found. It parses gcloud's --json output for the bucket preflight" >&2
    echo "  and the upload verification; without it neither check can run, and a check that" >&2
    echo "  cannot run must not be silently skipped." >&2
    exit 1; }

echo "gcloud   : $(command -v gcloud)"
echo "account  : $(gcloud config get-value account 2>/dev/null)"
echo "dest     : ${GCS_DEST}"
echo "staging  : ${STAGE_DIR}"
echo ""
sed 's/^/  /' "${PROV}"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STAGE D — destination preflight, then upload.
#
# This is individual-level controlled-access genotype data. The three checks below are the
# ones the 2026-09-15 audit of the sumstats bucket used, promoted from a post-hoc audit to a
# precondition: a bucket that fails any of them is not a place this data may go, and the
# script will not make that call for you by proceeding with a warning.
# ─────────────────────────────────────────────────────────────────────────────
echo "--- STAGE D: destination preflight ---"

DESC=$(gcloud storage buckets describe "gs://${BUCKET}" --format=json 2>/dev/null) || {
    echo "ERROR: cannot describe gs://${BUCKET} — it does not exist, or this account cannot see it." >&2
    exit 1; }

# gcloud has emitted these two fields under both snake_case and camelCase across releases, so
# they are looked up by name at any depth rather than at a fixed path. A key that is genuinely
# absent comes back as the empty string and FAILS the check below — it is never read as a pass.
read -r UBLA PAP <<< "$(printf '%s' "$DESC" | python3 -c '
import json, sys

def find(node, *names):
    """First value under any of `names`, searched breadth-first. Missing -> ""."""
    queue = [node]
    while queue:
        cur = queue.pop(0)
        if isinstance(cur, dict):
            for k, v in cur.items():
                if k in names:
                    return v
            queue.extend(cur.values())
        elif isinstance(cur, list):
            queue.extend(cur)
    return ""

try:
    d = json.load(sys.stdin)
except Exception:
    print("PARSE_ERROR PARSE_ERROR"); sys.exit(0)

ubla = find(d, "uniform_bucket_level_access", "uniformBucketLevelAccess")
if isinstance(ubla, dict):
    ubla = ubla.get("enabled", "")
pap = find(d, "public_access_prevention", "publicAccessPrevention")
# No f-string here: this runs inside a single-quoted shell argument, so it cannot contain a
# single quote, and an escaped double quote inside an f-string expression is a SyntaxError.
# A broken parser fails CLOSED (both values read as absent, both checks FAIL), which is the
# right direction but would also refuse a perfectly good bucket.
ubla = str(ubla) if ubla != "" else "<absent>"
pap = str(pap) if pap != "" else "<absent>"
print(ubla + " " + pap)
')"

IAM=$(gcloud storage buckets get-iam-policy "gs://${BUCKET}" --format=json 2>/dev/null) || {
    echo "ERROR: cannot read the IAM policy on gs://${BUCKET}." >&2
    echo "  The public-binding check CANNOT RUN, so there is no evidence this bucket is private." >&2
    echo "  Refusing to upload controlled-access genotypes to an unverified destination." >&2
    exit 1; }
PUBLIC=$(printf '%s' "$IAM" | grep -cE '"(allUsers|allAuthenticatedUsers)"')

FAIL=0
printf "  uniform bucket-level access : %s" "$UBLA"
if [[ "$UBLA" == "True" || "$UBLA" == "true" ]]; then echo "  OK"; else echo "  FAIL (must be enabled)"; FAIL=1; fi
printf "  public access prevention    : %s" "${PAP:-<unset>}"
if [[ "$PAP" == "enforced" ]]; then echo "  OK"; else echo "  FAIL (must be 'enforced')"; FAIL=1; fi
printf "  public IAM bindings         : %s" "$PUBLIC"
if [[ "$PUBLIC" -eq 0 ]]; then echo "  OK"; else echo "  FAIL (allUsers/allAuthenticatedUsers present)"; FAIL=1; fi

if [[ "$FAIL" -ne 0 ]]; then
    echo "" >&2
    echo "REFUSING TO UPLOAD. gs://${BUCKET} is not verifiably private, and this release is" >&2
    echo "  individual-level genotype data under the AMP-PD data use agreement. Fix the bucket" >&2
    echo "  configuration, or pick a different destination. There is no override flag here on" >&2
    echo "  purpose — the DUA decision does not belong to a shell variable." >&2
    exit 1
fi

EXISTING=$(gcloud storage ls "${GCS_DEST}/**" 2>/dev/null | grep -c . || true)
if [[ "${EXISTING:-0}" -gt 0 && -z "${OVERWRITE}" ]]; then
    echo "" >&2
    echo "ERROR: ${EXISTING} object(s) already exist under ${GCS_DEST}." >&2
    echo "  Refusing to overwrite a published release. Use a new prefix, or OVERWRITE=1." >&2
    exit 1
fi

echo ""
echo "--- uploading ---"

# The upload list comes from MANIFEST.tsv, not from a shell glob. A glob over the staging tree
# would also sweep up the per-stratum .plink.log files, and an unmatched brace pattern (say
# amppd_*.bed on a pgen build) expands to a literal nonexistent path that gcloud then errors on.
# The manifest already is the authoritative list of what this release contains.
UPLOAD_LIST=()
while IFS=$'\t' read -r F _ _ _; do
    [[ "$F" == "file" ]] && continue
    [[ -f "${STAGE_DIR}/${F}" ]] || { echo "ERROR: ${F} is in MANIFEST.tsv but not on disk" >&2; exit 1; }
    UPLOAD_LIST+=("${STAGE_DIR}/${F}")
done < "${MANIFEST_TSV}"
[[ ${#UPLOAD_LIST[@]} -gt 0 ]] || { echo "ERROR: MANIFEST.tsv lists no data objects" >&2; exit 1; }
UPLOAD_LIST+=("${MANIFEST_TSV}" "${RELEASE_README}")
echo "  ${#UPLOAD_LIST[@]} objects"

UPLOAD_START=$(date -u +%Y-%m-%dT%H:%M:%SZ)
gcloud storage cp "${UPLOAD_LIST[@]}" "${GCS_DEST}/" 2>&1 | grep -v '^$'
RC=${PIPESTATUS[0]}
UPLOAD_END=$(date -u +%Y-%m-%dT%H:%M:%SZ)
[[ "$RC" -eq 0 ]] || { echo "ERROR: gcloud storage cp failed (exit ${RC}). Release is PARTIAL — do not announce it." >&2; exit "$RC"; }

# ─────────────────────────────────────────────────────────────────────────────
# STAGE E — verify what landed against what was built.
#
# `cp` exiting 0 is not evidence that the release is complete: the failure this catches is a
# file that was never in the argument list at all, which no per-object transfer error would
# report. Presence and size are checked always; CRC32C whenever the build recorded it.
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "--- STAGE E: verify ---"

# Listed through a format projection rather than by parsing `ls -L` prose: that output is
# human-facing, has changed shape between gsutil and gcloud storage, and renders sizes as
# "12345 (12.05 KiB)" — a positional read of which yields "KiB)". Rule 4, applied to a listing.
REMOTE=$(mktemp)
gcloud storage ls --json --recursive "${GCS_DEST}/**" 2>/dev/null \
  | python3 -c '
import json, sys

try:
    items = json.load(sys.stdin)
except Exception:
    sys.exit(0)                     # empty output -> the guard below reports it
if isinstance(items, dict):
    items = [items]

for it in items:
    md = it.get("metadata", it)
    url = it.get("url") or md.get("name") or ""
    if not url or url.endswith("/"):
        continue
    name = url.rstrip("/").rsplit("/", 1)[-1]
    size = md.get("size", md.get("Content-Length", ""))
    crc  = md.get("crc32c", md.get("crc32c_hash", ""))
    size = str(size) if size != "" else "-"
    crc  = str(crc) if crc else "-"
    print(name + "\t" + size + "\t" + crc)
' > "$REMOTE"

if [[ ! -s "$REMOTE" ]]; then
    echo "ERROR: the bucket listing came back EMPTY after an upload that reported success." >&2
    echo "  Either nothing landed, or this gcloud build's --json shape is not the one parsed" >&2
    echo "  above. Either way the release is UNVERIFIED — check by hand before announcing it:" >&2
    echo "    gcloud storage ls -L '${GCS_DEST}/**'" >&2
    rm -f "$REMOTE"
    exit 1
fi

MISMATCH=0
NO_CRC=0
while IFS=$'\t' read -r F SZ CRC SHA; do
    [[ "$F" == "file" ]] && continue
    R=$(awk -F'\t' -v f="$F" '$1==f {print $2 "\t" $3; exit}' "$REMOTE")
    if [[ -z "$R" ]]; then
        echo "  MISSING in bucket: ${F}" >&2; MISMATCH=1; continue
    fi
    RSZ=${R%%$'\t'*}; RCRC=${R##*$'\t'}
    if [[ "$RSZ" != "$SZ" ]]; then
        echo "  SIZE MISMATCH ${F}: local ${SZ}, bucket ${RSZ}" >&2; MISMATCH=1; continue
    fi
    if [[ "$CRC" == "—" || "$RCRC" == "-" ]]; then
        NO_CRC=$((NO_CRC + 1))
    elif [[ "$CRC" != "$RCRC" ]]; then
        echo "  CRC32C MISMATCH ${F}: local ${CRC}, bucket ${RCRC}" >&2; MISMATCH=1
    fi
done < "${MANIFEST_TSV}"

# README.md and MANIFEST.tsv are uploaded but are not rows in the manifest (it describes the
# data objects), so they are checked for presence separately rather than assumed.
for EXTRA in README.md MANIFEST.tsv; do
    awk -F'\t' -v f="$EXTRA" '$1==f {found=1} END {exit !found}' "$REMOTE" \
        || { echo "  MISSING in bucket: ${EXTRA}" >&2; MISMATCH=1; }
done

N_LOCAL=$(($(wc -l < "${MANIFEST_TSV}" | tr -d ' ') - 1 + 2))
N_REMOTE=$(wc -l < "$REMOTE" | tr -d ' ')
echo "  objects: ${N_REMOTE} in bucket, ${N_LOCAL} expected"
[[ "$N_REMOTE" -eq "$N_LOCAL" ]] || { echo "  OBJECT COUNT MISMATCH" >&2; MISMATCH=1; }

if [[ "$NO_CRC" -gt 0 ]]; then
    echo "  !! ${NO_CRC} object(s) verified on SIZE ONLY — no CRC32C on one side." >&2
    echo "     Size agreement is not content agreement. Rebuild with gcloud available to" >&2
    echo "     get a real integrity check." >&2
fi
rm -f "$REMOTE"

if [[ "$MISMATCH" -ne 0 ]]; then
    echo "" >&2
    echo "VERIFY FAILED. The release at ${GCS_DEST} does not match the staging tree." >&2
    echo "  Do not announce or link it until this is resolved." >&2
    exit 1
fi

{
    echo "pushed_dest=${GCS_DEST}"
    echo "pushed_start_utc=${UPLOAD_START}"
    echo "pushed_end_utc=${UPLOAD_END}"
    echo "pushed_host=$(hostname)"
    echo "pushed_by=${USER:-unknown}"
    echo "pushed_account=$(gcloud config get-value account 2>/dev/null)"
} >> "${PROV}"

echo ""
banner "Step 9 PUSH complete $(date)"
echo "  ${GCS_DEST}  — ${N_REMOTE} objects, verified against MANIFEST.tsv"
echo ""
echo "RECORD THIS IN PROJECT_LOG.md: destination, object count, byte total, the upload window"
echo "above, and the git commit in ${PROV}. The previous release's exact cp command and host"
echo "were never written down and had to be recovered by audit five weeks later."
echo ""
exit 0
fi

echo "ERROR: MODE must be 'build' or 'push' (got '${MODE}')" >&2
exit 2
