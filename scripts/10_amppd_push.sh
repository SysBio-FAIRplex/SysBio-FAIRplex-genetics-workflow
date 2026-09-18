#!/bin/bash
#
# STEP 10 (follow-up, not tier 1) — publish the AMP-PD staging tree to a GCS bucket.
#
# RUNS ON HELIX, NOT BIOWULF, AND NOT UNDER SBATCH. Biowulf compute nodes have no general
# outbound network; an upload attempted from one fails in a way that looks like a credentials
# or bucket problem, which is an expensive misdiagnosis. This script refuses to start inside a
# SLURM allocation.
#
#   ssh helix.nih.gov
#   cd <project root>
#   GCS_DEST=gs://<bucket>/<prefix> bash scripts/10_amppd_push.sh
#
# INPUT: the staging tree scripts/09_amppd_subset.sh wrote — MANIFEST.tsv, README.md,
# .provenance and the per-stratum filesets. This script builds nothing and subsets nothing; if
# the tree is absent or incomplete it refuses rather than improvising.
#
# DESTINATION CHECK. This is individual-level genotype data from two controlled-access
# programs; redistribution is governed by the AMP-PD data use agreement. Stage D asks the
# bucket whether it is private — uniform bucket-level access on, public access prevention
# enforced, no allUsers / allAuthenticatedUsers binding.
#
# Those reads need storage.buckets.get and .getIamPolicy, which roles/storage.objectAdmin does
# NOT grant. On a program-managed bucket you can write to but not introspect, the check simply
# cannot run — which is the normal case, not an exception. It is reported loudly and the upload
# proceeds. A bucket whose metadata reads fine and says PUBLIC is still refused: that is a
# positive finding rather than an absent one. See stage D.
#
# INTEGRITY. Step 9's manifest carries size + SHA-256, computed with coreutils and no network.
# This script verifies the upload on three things: every expected object PRESENT, SIZE equal,
# and CRC32C equal between the local file and the bucket. CRC32C rather than MD5 because GCS
# computes no MD5 for composite (parallel-chunked) uploads — 42 of the 46 objects in the
# 2026-08-25 sumstats release came back MD5-less for exactly that reason. `cp` exiting 0 is not
# evidence the release is complete: the failure that misses is a file never in the argument
# list, which no per-object transfer error would report.
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
# GUARDRAIL: uploads sample IDs and genotypes. Run it yourself. The AI must not run it.
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

# Biowulf compute nodes have no general outbound network. A push attempted from one fails in a
# way that looks like a credentials or bucket problem, which is an expensive misdiagnosis.
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
# STAGE D — destination check, then upload.
#
# The three checks are the ones the 2026-09-15 audit of the sumstats bucket used. Reading them
# needs storage.buckets.get and .getIamPolicy, which roles/storage.objectAdmin does NOT grant —
# so on a program-managed bucket you can write to but not introspect, they cannot run at all.
#
# This used to refuse in that case, on the reasoning that an unreadable IAM policy is not
# evidence of a private bucket. True, but it blocked the normal situation while proving nothing,
# and the operator's knowledge of the bucket beats a 403. So:
#
#   metadata UNREADABLE  -> say so loudly, and continue
#   metadata READABLE and says PUBLIC -> refuse; that is a positive finding, not a missing one
#   SKIP_BUCKET_CHECK=1  -> skip stage D's checks entirely
#
# What is NOT relaxed: a bucket that can be neither described nor listed is a wrong name, and
# still exits.
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
        # Separate "cannot introspect" from "does not exist": listing objects is a different
        # permission, so if that works the bucket is real and writable and only the
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

# The upload list comes from MANIFEST.tsv, not from a shell glob. A glob over the staging tree
# would also sweep up the per-stratum .plink.log files, and an unmatched brace pattern (say
# amppd_*.bed on a pgen build) expands to a literal nonexistent path that gcloud then errors on.
# The manifest already is the authoritative list of what this release contains.
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
# STAGE E — verify what landed against what was built.
#
# `cp` exiting 0 is not evidence that the release is complete: the failure this catches is a
# file that was never in the argument list at all, which no per-object transfer error would
# report. Three checks, all mandatory: PRESENT, SIZE, CRC32C.
#
# CRC32C is computed HERE rather than read from the manifest. Step 9 has no gcloud and records
# size + SHA-256 instead — which is the better pair to ship, since a consumer can check SHA-256
# with coreutils after downloading. CRC32C is what GCS itself stores, so it is the only hash
# that can be compared against the bucket without downloading the objects back.
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
while IFS=$'\t' read -r F SZ SHA; do
    [[ "$F" == "file" ]] && continue
    R=$(awk -F'\t' -v f="$F" '$1==f {print $2 "\t" $3; exit}' "$REMOTE")
    if [[ -z "$R" ]]; then
        echo "  MISSING in bucket: ${F}" >&2; MISMATCH=1; continue
    fi
    RSZ=${R%%$'\t'*}; RCRC=${R##*$'\t'}
    if [[ "$RSZ" != "$SZ" ]]; then
        echo "  SIZE MISMATCH ${F}: local ${SZ}, bucket ${RSZ}" >&2; MISMATCH=1; continue
    fi
    LCRC=$(gcloud storage hash --skip-md5 "${STAGE_DIR}/${F}" 2>/dev/null \
           | awk -F': *' 'tolower($1) ~ /crc32c/ {print $2; exit}')
    if [[ -z "$LCRC" || "$RCRC" == "-" ]]; then
        # Not downgraded to a warning: a hash that could not be computed on either side leaves
        # this object's CONTENT unverified, and "size matched" is not the check this release is
        # supposed to get. Absence of a warning is not evidence.
        echo "  CRC32C UNAVAILABLE ${F}: local '${LCRC:-none}', bucket '${RCRC}'" >&2
        MISMATCH=1; continue
    fi
    if [[ "$LCRC" != "$RCRC" ]]; then
        echo "  CRC32C MISMATCH ${F}: local ${LCRC}, bucket ${RCRC}" >&2; MISMATCH=1
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
echo "RECORD THIS IN PROJECT_LOG.md: destination, object count, byte total, the upload window"
echo "above, and the git commit in ${PROV}. The previous release's exact cp command and host"
echo "were never written down and had to be recovered by audit five weeks later."
echo ""
exit 0
