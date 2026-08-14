#!/bin/bash
#SBATCH --job-name=genotools
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --partition=norm
#
# ─── REVERTED TO THE PROVEN VERSION — 2026-08-11 ────────────────────────────────
# Byte-identical (executable lines) to scripts_archive's run_genotools.sh, the script that
# actually completed wgs_harm, wb_dwgs and divco_hs on 2026-07-07. Restored as a debugging
# baseline because the config.sh-based rewrite has a 0% success rate: its only two runs were
# br_dsnwgs (26843275, 27057269) and both failed, so "br_dsnwgs is the odd one out" and
# "this script is the odd one out" were never distinguishable.
#
# CONSEQUENCE: the three paths below are absolute again, so this is the one file in the
# project that is not location-independent. That is deliberate and temporary — it removes
# config.sh, $BUNDLE resolution, the venv path and the module names as variables. Re-apply
# the config.sh refactor only once a genotools run has succeeded from this baseline.
#
# Known hazard in the version this replaces: no `set -e`, so a failed `source config.sh`
# left REF_PANEL/REF_LABELS/VENV empty and the script carried on regardless.
#
# Reusable GenoTools runner — submit with dataset-specific variables:
#
#   sbatch \
#     --job-name=genotools_<dataset> \
#     --output=/data/CARDPB2/sysbio/wgs/scripts/logs/genotools_<dataset>.o \
#     --error=/data/CARDPB2/sysbio/wgs/scripts/logs/genotools_<dataset>.e \
#     --export=PGEN=<pgen_prefix>,SEX_FILE=<sex_update_txt>,OUT_DIR=<output_dir>,DATASET=<dataset_name> \
#     /data/CARDPB2/sysbio/wgs/scripts/run_genotools.sh
#
# Required variables:
#   PGEN       — input pgen prefix (raw, before sex update)
#   SEX_FILE   — path to plink-format sex update file in dataset metadata/
#   OUT_DIR    — directory for genotools output
#   DATASET    — short name for logging (e.g. wb_dwgs, wgs_harm, divco_hs)
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
# --sort-vars is REQUIRED, not cosmetic. GenoTools' ancestry step aligns the study genotype
# matrix to the reference panel BY COLUMN POSITION, not by variant name:
#     ancestry.py:247   raw_geno.columns = col_names      # col_names = the PANEL's SNP order
# The reorder that would make that safe (ancestry.py:244) sits inside `if not self.train`, so
# it never runs when training a model from --ref_panel/--ref_labels. If this callset's variant
# order differs from the panel's, every column is standardized by another variant's mean/SD and
# projected through another variant's loading. Counts still match, nothing raises, exit is 0.
#
# BR-DSNWGS was ordered 1,10,11,...,19,2,20,21,22,3,...  — per-chromosome files concatenated in
# ALPHABETICAL order, because "10" sorts before "2". The panel is numeric. Result: all 209,068
# shared columns misaligned, ancestry signal destroyed, and all 97 samples labelled CAH (job
# 27211436). The tell was the projected PCs' sd being FLAT across PC1-PC10 (3.4/1.9/1.8/...)
# where the panel decays 103.9 -> 18.1 — flat spread on every axis is noise, not genomes.
#
# wgs_harm/divco_hs/wb_dwgs were unaffected only because their files happened to be numeric.
# Verify any callset with:  python3 scripts/diag_order.py <panel>.bim <callset>.pvar
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
# GenoTools sizes its GridSearchCV worker pool from the NODE rather than from this job:
#     ancestry.py:516   n_jobs = min(os.cpu_count(), max_workers_by_ram)
# where os.cpu_count() and psutil.virtual_memory() both report node totals and ignore the cgroup.
# On a 192-core node that is 192 loky workers however few CPUs we asked for — job 27057269
# requested 16 and still got 128. biowulf's RLIMIT_NPROC is 1024, PER-USER and node-wide across
# every job you have on that node, and 192 full interpreters plus their threads do not fit.
# pthread_create returns EAGAIN, the C++ layer throws an uncaught std::runtime_error, the worker
# aborts on SIGABRT, joblib raises TerminatedWorkerError, and the parent exits 1 — so SLURM
# records a bare `ExitCode 1:0` with no OUT_OF_MEMORY and no MaxRSS. That is why this read as
# anything but a resource limit for four months. The wrapper reports SLURM_CPUS_PER_TASK and
# SLURM_MEM_PER_NODE instead.
#
# BUT THE ALLOCATION ALONE IS NOT ENOUGH — 64 workers STILL FAILS. An earlier version of this
# comment claimed 64 "leaves ~3x headroom"; the scan below refutes it:
#
#     workers | grid-search time | outcome    | aborted workers
#     --------|------------------|------------|----------------
#        192  |      36:15       | failed     | 63
#         64  |      42:57       | failed     | 13
#         16  |      62:28       | COMPLETED  | 0
#
# So GENOTOOLS_MAX_WORKERS=16 is exported below and is load-bearing. Job 27211436 — the only
# BR-DSNWGS run ever to complete, in 8 attempts — used exactly this, passed by hand via
# --export. Baking it in is what makes that run reproducible from the repo.
#
# The bracket needs no guessing: 64 workers failing means >1024/64 = 16 tasks per worker;
# 16 succeeding means <1024/16 = 64. TASKS_PER_WORKER is therefore in (16, 64], and 32 with the
# half-budget rule reproduces the winning config exactly: (1024 // 2) // 32 = 16.
#
# Capping costs nothing in throughput. 16 workers finished the same 1080 fits in 62 min; the
# July successes on 128-192 workers took 1:11:59 end to end. The oversubscription was buying
# nothing — so there is no speed/safety tradeoff to weigh here.
#
# NO THREAD ENV VARS HERE, deliberately — an earlier version of this script exported
# OMP/OPENBLAS/MKL/NUMEXPR/NUMBA_NUM_THREADS=1 and it was pure superstition. joblib ALREADY caps
# threads inside each worker: its MAX_NUM_THREADS_VARS list includes NUMBA_NUM_THREADS, and the
# per-worker default is max(cpu_count() // n_jobs, 1), which is 1 when n_jobs equals the core
# count. Setting them by hand changes nothing joblib was not already doing. Job 27162204 ran with
# them removed and failed identically. The worker COUNT is the lever, not threads per worker.
#
# Not specific to br_dsnwgs: only 3 of 18 genotools jobs have ever completed, and each of
# wgs_harm / wb_dwgs / divco_hs needed 4-5 attempts. Whether a run survived came down to how much
# of the 1024 budget was free on that node at that moment. The fit is on the 4,008-genome
# reference panel, so it does identical work for 97 samples or 10,000.
#
# Timing: the grid search dies ~36 min after IT starts, not after the job starts — on 27162204
# the plink steps took 36 SECONDS. Measure "Conversion complete" -> "End time", not SLURM Elapsed.
#
# Expect this to run LONGER than the 1:11 successes: a third of the workers doing the same 1080
# fits. --time=48:00:00 has ample room.
#
# See PROJECT_LOG.md 2026-08-12 for the full chain and the hypotheses that were ruled out.
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