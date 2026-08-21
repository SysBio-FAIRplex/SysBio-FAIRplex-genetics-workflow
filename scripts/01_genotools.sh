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
# KNOWN DEVIATION: the three paths below are ABSOLUTE, so this is the one file in the project that
# is not location-independent. It was reverted to a known-good absolute-path baseline 2026-08-11
# while genotools was failing; all four callsets have since completed, so re-applying the config.sh
# refactor is unblocked and outstanding. See HANDOFF "Known issues".
#
# Fixed paths (shared across all datasets):
REF_PANEL="/data/CARDPB2/sysbio/wgs/data/ref/ref_panel_gp2_prune_rm_underperform_pos_update"
REF_LABELS="/data/CARDPB2/sysbio/wgs/data/ref/ref_panel_ancestry_updated.txt"
WGS_ROOT="/data/CARDPB2/sysbio/wgs"

# Derived paths
SEXUPD="${PGEN}_sexupd"
FILTERED="${PGEN}_filtered"
BED="${PGEN}_filtered_bed"
OUT="${OUT_DIR}/FILTERED.${DATASET}"

mkdir -p "${OUT_DIR}"

module load plink/6-alpha
module load python/3.11

cd "${WGS_ROOT}"
source .venv/bin/activate

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
# --sort-vars is REQUIRED, not cosmetic. GenoTools aligns the study matrix to the reference panel
# BY COLUMN POSITION (ancestry.py:247), and the reorder that would make that safe (:244) sits
# inside `if not self.train`, so it never runs when training from --ref_panel/--ref_labels. A
# callset whose order differs has every column standardized by another variant's mean/SD and
# projected through another variant's loading. Counts match, nothing raises, exit is 0.
#
# BR-DSNWGS was ordered 1,10,11,...,2,20 (per-chromosome files concatenated alphabetically) against
# a numeric panel: 209,068 columns misaligned, all 97 samples labelled CAH. The other three were
# unaffected only because their files happened to be numeric.
# VERIFY ANY NEW CALLSET:  python3 scripts/diag_order.py <panel>.bim <callset>.pvar
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
#
# GENOTOOLS_MAX_WORKERS=16 IS LOAD-BEARING — do not raise it. GenoTools sizes its GridSearchCV
# pool from the NODE (os.cpu_count()), not from this job, and biowulf's RLIMIT_NPROC of 1024 is
# per-user and node-wide. Over-sized pools abort workers on SIGABRT and SLURM records a bare
# `ExitCode 1:0` with no OUT_OF_MEMORY, which is why this read as anything but a resource limit
# for four months. The allocation alone is NOT enough — 64 workers still fails:
#
#     workers | grid-search time | outcome    | aborted workers
#        192  |      36:15       | failed     | 63
#         64  |      42:57       | failed     | 13
#         16  |      62:28       | COMPLETED  | 0
#
# Capping costs nothing: 16 workers ran the same 1080 fits in 62 min against 1:11:59 end to end on
# 128-192. There is no speed/safety tradeoff here.
#
# NO THREAD ENV VARS, deliberately. joblib already caps threads inside each worker; exporting
# OMP/OPENBLAS/MKL/NUMEXPR/NUMBA_NUM_THREADS=1 changes nothing and a run with them removed failed
# identically. The worker COUNT is the lever.
#
# The fit is on the 4,008-genome reference panel, so it does identical work for 97 samples or
# 10,000. Timing: the grid search dies ~36 min after IT starts, not after the job starts — measure
# "Conversion complete" -> "End time", not SLURM Elapsed.
#
# Full chain and ruled-out hypotheses: PROJECT_LOG.md 2026-08-12.
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