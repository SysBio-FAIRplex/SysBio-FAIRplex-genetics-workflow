#!/bin/bash
#
# STEP 10 (release, not tier 1) — publish step 9's AMP-PD staging tree to a GCS bucket.
#
# Runs on HELIX, never under sbatch: biowulf compute nodes have no outbound network, and an upload
# from one fails looking like an auth problem. Refuses inside a SLURM allocation.
#
#   ssh helix.nih.gov && cd <project root>
#   GCS_DEST=gs://<bucket>/<prefix> bash scripts/10_amppd_push.sh
#
# Builds nothing: reads MANIFEST.tsv, README.md, .provenance and the filesets step 9 wrote, and
# refuses if they are absent. Stage D checks the bucket is private (UBLA on, public access
# prevention enforced, no allUsers/allAuthenticatedUsers); metadata it cannot read is reported
# loudly and the upload proceeds, metadata that reads PUBLIC is refused. Stage E verifies every
# object PRESENT, SIZE equal and CRC32C equal (not MD5: GCS has none for composite uploads).
#
# KNOBS
#   GCS_DEST              gs://bucket/prefix — REQUIRED. No bucket is hardcoded anywhere.
#   SKIP_BUCKET_CHECK=1   skip stage D's privacy checks entirely (they need bucket-level read
#                         permissions that roles/storage.objectAdmin does not grant).
#   STAGE_DIR             the tree step 9 wrote (default data/merged/release_amppd)
#   OVERWRITE=1           allow objects already under GCS_DEST to be replaced. Refuses without it.
#   MOD_GCLOUD            module providing gcloud (default google-cloud-sdk); skipped if gcloud
#                         is already on PATH.
#   ALLOW_COMPUTE_PUSH=1  run inside a SLURM allocation anyway, if you know the node has egress.
#
# GUARDRAIL: uploads individual-level AMP-PD genotypes and sample IDs. Human-run, and only with
# confirmation that the recipient is covered by the AMP-PD DUA (README, release section).
set -o pipefail

BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

STAGE_DIR=${STAGE_DIR:-${MERGED_DIR}/release_amppd}
MANIFEST_TSV=${STAGE_DIR}/MANIFEST.tsv
RELEASE_README=${STAGE_DIR}/README.md
PROV=${STAGE_DIR}/.provenance

banner() { echo "=========================================="; echo "$@"; echo "=========================================="; }


banner "Step 10 — publish the AMP-PD release to GCS"
echo "Host: $(hostname)   start $(date)"

[[ -n "${GCS_DEST}" ]] || {
    echo "ERROR: GCS_DEST is not set, and no bucket is hardcoded in this script." >&2
    echo "  GCS_DEST=gs://<bucket>/<prefix> bash scripts/10_amppd_push.sh" >&2
    exit 2; }
[[ "${GCS_DEST}" == gs://* ]] || { echo "ERROR: GCS_DEST must start with gs:// (got '${GCS_DEST}')" >&2; exit 2; }
GCS_DEST=${GCS_DEST%/}
BUCKET=${GCS_DEST#gs://}; BUCKET=${BUCKET%%/*}

[[ -f "${PROV}" && -f "${MANIFEST_TSV}" ]] || {
    echo "ERROR: no completed build at ${STAGE_DIR} (missing MANIFEST.tsv or .provenance)." >&2
    echo "  Run ./submit.sh scripts/09_amppd_subset.sh on biowulf first." >&2
    exit 1; }

# Compute nodes have no outbound network; refuse rather than fail looking like an auth problem.
if [[ -n "${SLURM_JOB_ID}" && -z "${ALLOW_COMPUTE_PUSH}" ]]; then
    echo "ERROR: this is running inside SLURM job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}." >&2
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

# Stages D and E parse gcloud's --json with python3; check for it now, not three stages later.
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
# STAGE D — destination check, then upload. The reads need storage.buckets.get/.getIamPolicy,
# which roles/storage.objectAdmin lacks, so on a program-managed bucket they often cannot run:
#   metadata UNREADABLE               -> say so loudly, and continue
#   metadata READABLE and says PUBLIC -> refuse
#   neither describable nor listable  -> wrong bucket name; exit
#   SKIP_BUCKET_CHECK=1               -> skip the checks
# ─────────────────────────────────────────────────────────────────────────────
echo "--- STAGE D: destination check ---"

SKIP_D=""
if [[ -n "${SKIP_BUCKET_CHECK}" ]]; then
    echo "  SKIP_BUCKET_CHECK=1 — no check has been made that ${GCS_DEST} is private."
    SKIP_D=1
    DESC=""
else
    DESC=$(gcloud storage buckets describe "gs://${BUCKET}" --format=json 2>/dev/null) || DESC=""
    if [[ -z "$DESC" ]]; then
        # Listing is a different permission: if it works, the bucket exists and only the
        # bucket-level reads are missing.
        if gcloud storage ls "gs://${BUCKET}/" >/dev/null 2>&1; then
            echo "  !! CANNOT READ BUCKET METADATA on gs://${BUCKET}."
            echo "     Objects list fine, so the bucket exists and is reachable; this account"
            echo "     lacks storage.buckets.get / .getIamPolicy (roles/storage.objectAdmin"
            echo "     does not grant them)."
            echo "     NOTHING HERE HAS VERIFIED THAT THIS BUCKET IS PRIVATE — proceeding on"
            echo "     the operator's say-so."
            SKIP_D=1
        else
            echo "ERROR: gs://${BUCKET} can be neither described NOR listed — it does not exist," >&2
            echo "  or this account cannot see it at all. Check the bucket name." >&2
            exit 1
        fi
    fi
fi

if [[ -z "$SKIP_D" ]]; then

# Looked up by name at any depth (gcloud has used snake_case and camelCase). An absent key FAILS.
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

IAM=$(gcloud storage buckets get-iam-policy "gs://${BUCKET}" --format=json 2>/dev/null) || IAM=""
if [[ -n "$IAM" ]]; then
    PUBLIC=$(printf '%s' "$IAM" | grep -cE '"(allUsers|allAuthenticatedUsers)"')
else
    PUBLIC="?"
fi

FAIL=0
printf "  uniform bucket-level access : %s" "$UBLA"
if [[ "$UBLA" == "True" || "$UBLA" == "true" ]]; then echo "  OK"; else echo "  FAIL (must be enabled)"; FAIL=1; fi
printf "  public access prevention    : %s" "${PAP:-<unset>}"
if [[ "$PAP" == "enforced" ]]; then echo "  OK"; else echo "  FAIL (must be 'enforced')"; FAIL=1; fi
printf "  public IAM bindings         : %s" "$PUBLIC"
if [[ "$PUBLIC" == "?" ]]; then echo "  UNREADABLE (not counted as a failure)"
elif [[ "$PUBLIC" -eq 0 ]]; then echo "  OK"
else echo "  FAIL (allUsers/allAuthenticatedUsers present)"; FAIL=1; fi

if [[ "$FAIL" -ne 0 ]]; then
    echo "" >&2
    echo "REFUSING TO UPLOAD. The bucket metadata READ FINE, and it says gs://${BUCKET} is not" >&2
    echo "  private. That is a positive finding rather than a missing one, and this release is" >&2
    echo "  individual-level genotype data under the AMP-PD data use agreement." >&2
    echo "  Fix the bucket configuration, pick a different destination, or SKIP_BUCKET_CHECK=1" >&2
    echo "  if you know something this check does not." >&2
    exit 1
fi
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

# The upload list is MANIFEST.tsv, not a glob (which would sweep up .plink.log files).
UPLOAD_LIST=()
while IFS=$'\t' read -r F _ _; do
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
# STAGE E — verify what landed: PRESENT, SIZE, CRC32C, all mandatory. `cp` exiting 0 cannot
# report a file that was never in its argument list. CRC32C is computed here (step 9 has no
# gcloud); it is the one hash GCS stores.
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "--- STAGE E: verify ---"

# --json, not `ls -L` prose, whose shape changes between gcloud releases.
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
while IFS=$'\t' read -r F SZ SHA; do
    [[ "$F" == "file" ]] && continue
    R=$(awk -F'\t' -v f="$F" '$1==f {print $2 "\t" $3; exit}' "$REMOTE")
    if [[ -z "$R" ]]; then
        echo "  MISSING in bucket: ${F}" >&2; MISMATCH=1; continue
    fi
    RSZ=${R%%$'\t'*}; RCRC=${R##*$'\t'}
    if [[ "$RSZ" != "$SZ" ]]; then
        echo "  SIZE MISMATCH ${F}: staged ${SZ}, bucket ${RSZ}" >&2; MISMATCH=1; continue
    fi
    LCRC=$(gcloud storage hash --skip-md5 "${STAGE_DIR}/${F}" 2>/dev/null \
           | awk -F': *' 'tolower($1) ~ /crc32c/ {print $2; exit}')
    if [[ -z "$LCRC" || "$RCRC" == "-" ]]; then
        # A hash missing on either side leaves CONTENT unverified: a failure, not a warning.
        echo "  CRC32C UNAVAILABLE ${F}: staged '${LCRC:-none}', bucket '${RCRC}'" >&2
        MISMATCH=1; continue
    fi
    if [[ "$LCRC" != "$RCRC" ]]; then
        echo "  CRC32C MISMATCH ${F}: staged ${LCRC}, bucket ${RCRC}" >&2; MISMATCH=1
    fi
done < "${MANIFEST_TSV}"

# README.md and MANIFEST.tsv are not manifest rows; check their presence separately.
for EXTRA in README.md MANIFEST.tsv; do
    awk -F'\t' -v f="$EXTRA" '$1==f {found=1} END {exit !found}' "$REMOTE" \
        || { echo "  MISSING in bucket: ${EXTRA}" >&2; MISMATCH=1; }
done

N_STAGED=$(($(wc -l < "${MANIFEST_TSV}" | tr -d ' ') - 1 + 2))
N_REMOTE=$(wc -l < "$REMOTE" | tr -d ' ')
echo "  objects: ${N_REMOTE} in bucket, ${N_STAGED} expected"
[[ "$N_REMOTE" -eq "$N_STAGED" ]] || { echo "  OBJECT COUNT MISMATCH" >&2; MISMATCH=1; }

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
banner "Step 10 complete $(date)"
echo "  ${GCS_DEST}  — ${N_REMOTE} objects, verified against MANIFEST.tsv"
echo ""
echo "RECORD THIS IN METHODS.md (data availability): destination, object count, byte total, the"
echo "upload window above, and the git commit in ${PROV}. The previous release's exact cp command"
echo "and host were never written down and had to be recovered by audit five weeks later."
echo ""
exit 0
