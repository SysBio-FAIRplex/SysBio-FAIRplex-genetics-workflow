#!/bin/bash
#
# Build the Python environment. Run once on biowulf (login or interactive node, not sbatch):
#   bash scripts/00_setup_env.sh
#
# Load-bearing: setuptools pinned to 67.8.0 (later versions break GenoTools' pkg_resources, but
# only under sbatch), and `module load python` BEFORE the venv is created.

set -euo pipefail

BUNDLE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${BUNDLE}/config.sh"

echo "=========================================="
echo "Building ${VENV}"
echo "Python module: ${MOD_PYTHON}"
echo "=========================================="

module load "${MOD_PYTHON}"

if [[ -d "${VENV}" ]]; then
    echo "NOTE: ${VENV} already exists — installing into it."
    echo "      Delete it first for a clean rebuild."
else
    mkdir -p "$(dirname "${VENV}")"
    python3 -m venv "${VENV}"
fi

source "${VENV}/bin/activate"

# The setuptools pin must be in place before anything else resolves against it.
pip install --upgrade pip
pip install 'setuptools==67.8.0'
pip install -r "${BUNDLE}/scripts/requirements.txt"

# Only for executing notebooks headless; unpinned, outside the validated set, imported by no step.
pip install nbconvert ipykernel

echo
echo "=========================================="
python3 -c "import genotools, pandas, numpy; print('genotools', genotools.__version__)" 2>/dev/null \
    || echo "WARNING: genotools did not import — check the pip log above"
echo "setuptools: $(python3 -c 'import setuptools; print(setuptools.__version__)')"
echo
echo "Done. Every later step activates this itself; to use it by hand:"
echo "  module load ${MOD_PYTHON} && source ${VENV}/bin/activate"
echo "=========================================="
