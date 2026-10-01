#!/bin/bash
#SBATCH --job-name=amppd_subset
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --partition=norm
#
# STEP 9 (release, not tier 1) — subset the AMP-PD donors out of step 6 STAGE C
# (by_ancestry_qc/cohort_<ANC>_qc: post-QC and post AF-concordance filter, the exact variant set
# step 7 tested) into a staging tree. Runs after step 8; reads in place, modifies nothing upstream.
#
# NO NETWORK: sbatch on biowulf, where compute nodes have no outbound route. Publishing is
# scripts/10_amppd_push.sh, on helix.
#
# One fileset PER STRATUM (QC is per stratum; concatenating would invent a variant set no GWAS ran
# on). NOT re-filtered after subsetting, so the release matches the published sumstats; SUBSET_MAF
# / SUBSET_GENO opt in and are recorded. AMP-PD = source_callset in {wb_dwgs, br_dsnwgs}, read by
# header name from step 6 stage E's retained_samples_manifest.csv — the one definition of callset
# membership. Manifest: size + SHA-256 (coreutils); step 10 adds CRC32C at upload.
#
# KNOBS
#   STAGE_DIR           where the release tree is built (default data/merged/release_amppd)
#   AMPPD_CALLSETS      "wb_dwgs br_dsnwgs"
#   OUT_FMT             pgen (default) | bed | vcf
#   SUBSET_MAF          re-apply --maf inside the subset (default: unset, do not re-filter)
#   SUBSET_GENO         re-apply --geno inside the subset (default: unset)
#   OVERWRITE=1         replace an existing staging tree. Refuses without it.
#
#   ./submit.sh scripts/09_amppd_subset.sh
#
# GUARDRAIL: the deliverable is individual-level AMP-PD genotypes and sample IDs. Human-run, and
# only with confirmation that the recipient is covered by the AMP-PD DUA (README, release
# section). Never read the staging tree's .psam/.fam into a chat. All output below is aggregate.
set -o pipefail

BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

QC_DIR=${MERGED_DIR}/by_ancestry_qc
STAGE_DIR=${STAGE_DIR:-${MERGED_DIR}/release_amppd}
AMPPD_CALLSETS=${AMPPD_CALLSETS:-"wb_dwgs br_dsnwgs"}
OUT_FMT=${OUT_FMT:-pgen}
THREADS=${SLURM_CPUS_PER_TASK:-16}
MANIFEST_TSV=${STAGE_DIR}/MANIFEST.tsv
RELEASE_README=${STAGE_DIR}/README.md
PROV=${STAGE_DIR}/.provenance

banner() { echo "=========================================="; echo "$@"; echo "=========================================="; }

# Stage B's per-stratum tally is an associative array (bash 4+). Without it the tally silently reads
# empty, and an accounting check that cannot run must not look like one that passed.
if [[ -z "${BASH_VERSINFO[0]}" || "${BASH_VERSINFO[0]}" -lt 4 ]]; then
    echo "ERROR: bash 4+ required (this is ${BASH_VERSION:-unknown})." >&2
    echo "  The per-stratum sample tally uses an associative array, and the release accounting" >&2
    echo "  depends on it. Run this on biowulf." >&2
    exit 1
fi

# SHA-256 with coreutils; absence is fatal (a hash-less row looks like one nobody checked).
sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
    else return 127; fi
}


banner "Step 9 — AMP-PD subset of the QC'd association set"
echo "Job ${SLURM_JOB_ID:-interactive} on ${SLURMD_NODENAME:-$(hostname)}   start $(date)"
echo "Source (step 6 stage C): ${QC_DIR}/cohort_<ANC>_qc"
echo "Sample source:           ${RETAINED_MANIFEST}"
echo "AMP-PD callsets:         ${AMPPD_CALLSETS}"
echo "Staging:                 ${STAGE_DIR}"
echo "Format:                  ${OUT_FMT}"
echo "Re-filter inside subset: maf=${SUBSET_MAF:-<none>} geno=${SUBSET_GENO:-<none>}"
echo ""

# ── preconditions: each fails loudly ──
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
# STAGE A — select the AMP-PD IIDs. Columns resolved from the header row, never by position.
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
# STAGE B — one subset per stratum. Keep files join the AMP-PD IIDs to the stratum's OWN .fam
# (headerless, so positional): it supplies the FID, and a sample counts as shipped only if it is
# actually in the fileset.
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

# ── the accounting: fatal, because a release that silently drops samples looks like one with none
# to drop.
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
# STAGE C — manifest + release README, both DERIVED from what stage B wrote; nothing transcribed.
# ─────────────────────────────────────────────────────────────────────────────
echo "--- STAGE C: manifest + release README ---"

# Size and SHA-256 only, from coreutils — no gcloud, no network.
printf 'file\tbytes\tsha256\n' > "${MANIFEST_TSV}"
for F in "${STAGE_DIR}"/amppd_*; do
    B=$(basename "$F")
    case "$B" in *.plink.log) continue ;; esac
    SZ=$(stat -c %s "$F")
    SHA=$(sha256_of "$F")
    [[ -n "$SHA" ]] || { echo "ERROR: sha256sum not found." >&2
                         echo "  Refusing to write a manifest whose hashes are blank." >&2; exit 1; }
    printf '%s\t%s\t%s\n' "$B" "$SZ" "$SHA" >> "${MANIFEST_TSV}"
done

N_OBJ=$(($(wc -l < "${MANIFEST_TSV}" | tr -d ' ') - 1))
TOT_BYTES=$(awk -F'\t' 'NR>1 {s+=$2} END {print s+0}' "${MANIFEST_TSV}")
TOT_GIB=$(awk -v b="$TOT_BYTES" 'BEGIN {printf "%.2f", b/1073741824}')

# ── provenance, machine-readable, consumed by scripts/10_amppd_push.sh ──
{
    echo "built_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "build_job=${SLURM_JOB_ID:-interactive}"
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
    echo "Generated by \`scripts/09_amppd_subset.sh\` on $(date -u +%Y-%m-%d) from the project's"
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
    echo "\`MANIFEST.tsv\` carries per-object size and SHA-256. After downloading:"
    echo ""
    echo "\`\`\`"
    echo "sha256sum -c <(awk -F'\\t' 'NR>1 {print \$3\"  \"\$1}' MANIFEST.tsv)"
    echo "\`\`\`"
    echo ""
    echo "Not MD5: GCS computes no MD5 for composite (parallel-chunked) uploads, so an MD5"
    echo "comparison against these objects silently has nothing to compare."
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
banner "Step 9 complete $(date)  — staging tree ready, nothing uploaded"
echo "NOTHING HAS BEEN UPLOADED — this script has no network code in it at all."
echo "Publishing is scripts/10_amppd_push.sh, and it runs on HELIX, not here:"
echo ""
echo "  ssh helix.nih.gov"
echo "  cd ${BUNDLE}"
echo "  GCS_DEST=gs://<bucket>/<prefix> bash scripts/10_amppd_push.sh"
echo ""
exit 0
