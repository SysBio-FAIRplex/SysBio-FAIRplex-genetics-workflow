#!/bin/bash
#SBATCH --job-name=genotools
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --partition=norm
#
# STEP 1 — per callset: apply the sex file, filter, convert to bed, then GenoTools ancestry + QC.
#
#   ./submit.sh scripts/01_genotools.sh --job-name=genotools_<ds> \
#     --export=PGEN=<raw_pgen_prefix>,SEX_FILE=<sex_update_txt>,OUT_DIR=<out>,DATASET=<ds>
#
# PGEN     raw pgen prefix, before the sex update
# SEX_FILE clinical_core_out/<callset>_update_sex.txt  (config.sh: sex_file)
# OUT_DIR  genotools output directory
# DATASET  short name for logging (wgs_harm, divco_hs, wb_dwgs, br_dsnwgs)
#
# No `set -e` here, so assert that config.sh resolved every path rather than run on empty ones.
BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh" || { echo "ERROR: cannot source ${BUNDLE}/config.sh" >&2; exit 1; }

for _v in REF_PANEL REF_LABELS WGS_ROOT VENV; do
    [[ -n "${!_v:-}" ]] || { echo "ERROR: ${_v} is empty after sourcing config.sh — refusing to run" >&2; exit 1; }
done
[[ -f "${REF_PANEL}.bim" ]] || { echo "ERROR: no reference panel at ${REF_PANEL}.bim" >&2; exit 1; }
[[ -f "${REF_LABELS}" ]]    || { echo "ERROR: no reference labels at ${REF_LABELS}" >&2; exit 1; }

# Derived paths
SEXUPD="${PGEN}_sexupd"
FILTERED="${PGEN}_filtered"
BED="${PGEN}_filtered_bed"
OUT="${OUT_DIR}/FILTERED.${DATASET}"

mkdir -p "${OUT_DIR}"

module load plink/6-alpha
module load python/3.11

# Absolute path: a relative activate would depend on the job's working directory.
source "${VENV}/bin/activate"
cd "${WGS_ROOT}"

echo "=========================================="
echo "GenoTools — ${DATASET}"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       ${SLURMD_NODENAME}"
echo "Input:      ${PGEN}"
echo "Sex file:   ${SEX_FILE}"
echo "Output:     ${OUT}"
echo "Start time: $(date)"
echo "=========================================="

# ── Step 1: Update sex ────────────────────────────────────────────────────────
# --sort-vars is REQUIRED. GenoTools aligns the study matrix to the reference panel BY COLUMN
# POSITION (ancestry.py:247), so a callset in a different variant order is silently mis-projected:
# exit 0, wrong ancestry. Before this step on any new callset:
#   python3 scripts/diag_order.py <panel>.bim <callset>.pvar
echo "Step 1: Updating sex..."
plink2 \
    --pfile "${PGEN}" \
    --update-sex "${SEX_FILE}" \
    --sort-vars \
    --make-pgen \
    --threads 64 \
    --out "${SEXUPD}"

if [[ $? -ne 0 ]]; then
    echo "ERROR: sex update failed" >&2; exit 1
fi
echo "Sex update complete: $(date)"

# ── Step 2: Filter to biallelic PASS SNPs + assign IDs + deduplicate ─────────
# --set-missing-var-ids: assigns CHROM:POS:REF:ALT to any variants with '.' ID (e.g. WB-DWGS)
# --rm-dup exclude-all: removes variants with duplicate IDs (e.g. WGS_Harm scatter interval overlaps)
echo "Step 2: Filtering to biallelic PASS SNPs..."
VAR_ID_TMPL='@:#:$r:$a'
plink2 \
    --pfile "${SEXUPD}" \
    --var-filter \
    --min-alleles 2 \
    --max-alleles 2 \
    --snps-only \
    --set-missing-var-ids "${VAR_ID_TMPL}" \
    --rm-dup exclude-all \
    --make-pgen \
    --threads 64 \
    --out "${FILTERED}"

if [[ $? -ne 0 ]]; then
    echo "ERROR: filtering step failed" >&2; exit 1
fi
echo "Filtering complete: $(date)"

# ── Step 3: Convert filtered pgen → bed ──────────────────────────────────────
echo "Step 3: Converting pgen → bed..."
plink2 \
    --pfile "${FILTERED}" \
    --merge-par \
    --make-bed \
    --threads 64 \
    --out "${BED}"

if [[ $? -ne 0 ]]; then
    echo "ERROR: pgen → bed conversion failed" >&2; exit 1
fi
echo "Conversion complete: $(date)"

# ── Step 4: GenoTools ────────────────────────────────────────────────────────
# Run through genotools_capped.py, NOT the bare `genotools` console script.
# GENOTOOLS_MAX_WORKERS=16 is load-bearing — do not raise it. GenoTools sizes its pool from the
# NODE, and biowulf's per-user `ulimit -u` of 1024 aborts larger pools (192 and 64 workers failed;
# 16 completed, no slower end to end). Thread env vars do not help; the worker count is the lever.
# The grid search fits the 4,008-genome reference panel, so its runtime is the same for any callset.
echo "Step 4: Running GenoTools..."
export GENOTOOLS_MAX_WORKERS="${GENOTOOLS_MAX_WORKERS:-16}"
echo "GENOTOOLS_MAX_WORKERS=${GENOTOOLS_MAX_WORKERS}"
python3 "${WGS_ROOT}/scripts/genotools_capped.py" \
    --bfile "${BED}" \
    --out "${OUT}" \
    --full_output \
    --warn \
    --skip_fails \
    --ancestry \
    --ref_panel "${REF_PANEL}" \
    --ref_labels "${REF_LABELS}" \
    --all_sample

EXIT_CODE=$?
echo "=========================================="
echo "End time: $(date)"
[[ ${EXIT_CODE} -eq 0 ]] && echo "Success." || echo "ERROR: exit code ${EXIT_CODE}" >&2
exit ${EXIT_CODE}