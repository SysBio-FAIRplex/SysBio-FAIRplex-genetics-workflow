#!/bin/bash
#SBATCH --job-name=excludelist
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32g
#SBATCH --partition=norm
#
# STEP 5 (sbatch wrapper) — run 05_excludelist.py to build the excludelist + retained manifest
# from the ancestry-split KING pairs + common-set call rates. It only writes text outputs and
# prints an aggregate summary; NO genotype file is modified here.
#
# It is small enough to run on an interactive node instead — see README "Run order", tier 1.
#
#   ./submit.sh scripts/05_excludelist.sh

set -o pipefail

# Resolve the bundle root and load the one file that knows where anything lives.
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
