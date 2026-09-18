#!/bin/bash
# 01 — pull the pseudobulk sample lists off biowulf.
#
# Recorded so the transfer is repeatable, not because it needs to be rerun often.
# Structure is preserved, never flattened: two of the six cohorts (CMD, rasle) put the
# cell type in the DIRECTORY and name every file samples.tsv/.csv, while the other four
# put it in the FILENAME. Flattening would collide and would lose the cohort/tissue/
# celltype labels that step 02 attributes each donor to.
#
# Idempotent — rerun to refresh.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REMOTE="${REMOTE:-${USER}@biowulf.nih.gov}"
SRC="${SRC:-/data/CARD/sysbio/data/beta1/pseudobulk/}"
DEST="${HERE}/pseudobulk/"

mkdir -p "${DEST}"

rsync -av --prune-empty-dirs \
  --include='*/' \
  --include='*samples.tsv' \
  --include='*samples.csv' \
  --exclude='*' \
  "${REMOTE}:${SRC}" "${DEST}"

echo "files: $(find "${DEST}" -name '*samples.[tc]sv' | wc -l | tr -d ' ')  (expected 179)"
