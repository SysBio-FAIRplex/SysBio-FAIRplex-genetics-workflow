#!/bin/bash
#
# STEP 0 — build the Python environment. Run once, on a login or interactive node.
#
#   bash scripts/00_setup_env.sh
#
# Not an sbatch job: it is a pip install, and it needs to be re-runnable while you watch it.
#
# Two things here are load-bearing and are the reason this is a script rather than a README line:
#
#   1. setuptools is pinned to 67.8.0. Later versions removed the pkg_resources behaviour GenoTools
#      depends on, and it fails ONLY under sbatch — interactive runs look fine, then the batch job
#      dies with an opaque import error.
#   2. `module load python` must happen BEFORE the venv is created. Building against system Python
#      produces an environment that fails the same way.

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

# Order matters: the setuptools pin must be in place before anything else resolves against it.
pip install --upgrade pip
pip install 'setuptools==67.8.0'
pip install -r "${BUNDLE}/scripts/requirements.txt"

# Needed only to execute clinical_core.ipynb headless (`jupyter nbconvert --execute`). Deliberately
# unpinned and installed separately — they are not part of the version set the pipeline was
# validated against, and nothing in steps 1-7 imports them.
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
