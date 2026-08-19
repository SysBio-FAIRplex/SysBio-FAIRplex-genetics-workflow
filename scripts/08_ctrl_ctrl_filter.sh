#!/bin/bash
#SBATCH --job-name=ctrl_ctrl_filter
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=norm
#
# STEP 8 — cross-reference every contrast against the control-vs-control scan.
#
# `control@amppd` vs `control@ampad` has no true disease signal by construction, so anything
# significant there is cohort artifact that survived step 6a. This is the demo's step 13.
#
# It ANNOTATES and writes a filtered COPY; it never modifies the primary sumstats. The two control
# arms are screened for different diseases, so APOE is expected to reach significance there for a
# real reason — filtering destructively would delete it. See the docstring in 08_ctrl_ctrl_filter.py.
#
# Reads summary statistics only, so it is cheap and re-runnable without redoing the GWAS.
#
#   ./submit.sh scripts/08_ctrl_ctrl_filter.sh
#
# Knobs: CTRL_P (1e-5, the demo's threshold) · GWSIG (5e-8) · OUT_SUBDIR (gwas)
#
# GUARDRAIL: reads sumstats only — no genotypes, no sample IDs. Safe for the AI to reason about;
# still run by the user for consistency with the rest of the bundle.
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
echo "Job ${SLURM_JOB_ID:-local} on ${SLURMD_NODENAME:-$(hostname)}   start $(date)"
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
