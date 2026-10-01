#!/bin/bash
#SBATCH --job-name=excludelist
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32g
#SBATCH --partition=norm
#
# STEP 5 (sbatch wrapper) — runs 05_excludelist.py: exclude list + retained manifest from the
# ancestry-split KING pairs and common-set call rates. Writes text files only; prints aggregates.
#
#   ./submit.sh scripts/05_excludelist.sh

set -o pipefail

BUNDLE="${BUNDLE:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "${BUNDLE}/config.sh"

module load "${MOD_PYTHON}"
source "${VENV}/bin/activate"

echo "=========================================="
echo "Step 5 — excludelist — Job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}"
echo "Start: $(date)"
echo "=========================================="

python3 "${BUNDLE}/scripts/05_excludelist.py"

echo "=========================================="
echo "End: $(date)"
