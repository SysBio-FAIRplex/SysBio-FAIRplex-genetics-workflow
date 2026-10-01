#!/bin/bash
#SBATCH --job-name=ctrl_ctrl_filter
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=norm
#
# STEP 8 — cross-reference every contrast against the control-vs-control scan.
# ANNOTATES and writes a filtered COPY; never modifies the primary sumstats (METHODS.md §8).
#
#   ./submit.sh scripts/08_ctrl_ctrl_filter.sh
#
# Knobs: CTRL_P (1e-5) · GWSIG (5e-8) · OUT_SUBDIR (gwas)
# Reads sumstats only — no genotypes, no sample IDs. Human-run.
set -o pipefail

BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

GWAS_DIR=${MERGED_DIR}/by_ancestry_qc/${OUT_SUBDIR:-gwas}
CTRL_P=${CTRL_P:-1e-5}
GWSIG=${GWSIG:-5e-8}

module load "${MOD_PYTHON}"
source "${VENV}/bin/activate"

echo "=========================================="
echo "Step 8 — control-vs-control cross-reference"
echo "Job ${SLURM_JOB_ID:-interactive} on ${SLURMD_NODENAME:-$(hostname)}   start $(date)"
echo "GWAS dir: ${GWAS_DIR}   ctrl P < ${CTRL_P}"
echo "=========================================="

[[ -d "$GWAS_DIR" ]] || { echo "ERROR: no GWAS output at $GWAS_DIR — run step 7 first" >&2; exit 1; }

python3 "${BUNDLE}/scripts/08_ctrl_ctrl_filter.py" \
    --gwas-dir "${GWAS_DIR}" \
    --ctrl-p "${CTRL_P}" \
    --gwsig "${GWSIG}"
RC=$?

echo ""
echo "=========================================="
echo "Done $(date)  (exit ${RC})"
echo "=========================================="
exit ${RC}
