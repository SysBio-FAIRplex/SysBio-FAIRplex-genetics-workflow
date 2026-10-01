# ad-pd-gwas — configuration. Sourced, never executed; every script sources it.
#
# No machine-specific paths: everything derives from this file's location. Any value can be
# overridden from the environment, e.g.  GRAIN=/other/grain.csv ./submit.sh scripts/07_gwas.sh

# ── the one root ─────────────────────────────────────────────────────────────
# Code and data share a root on biowulf; data/ may be a symlink to a larger allocation.
# WGS_ROOT is an exported alias the shell steps cd into and the Python helpers accept.
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export WGS_ROOT="${WGS_ROOT:-${PROJECT_ROOT}}"

DATA="${PROJECT_ROOT}/data"
REF_DIR="${DATA}/ref"                  # fetched, large, gitignored — see REF_* below
SHIP_REF="${PROJECT_ROOT}/ref"         # ships with the code, versioned, small
MERGED_DIR="${DATA}/merged"
VENV="${PROJECT_ROOT}/.venv"

# ── the phenotype boundary ───────────────────────────────────────────────────
# clinical_core.py / analysis_grain.py write here and the genotype steps read here. Nothing is
# copied into data/, so a handoff file cannot go stale against the run that produced it.
CLINICAL_OUT="${CLINICAL_OUT:-${PROJECT_ROOT}/clinical_core_out}"

# ── the four source callsets ─────────────────────────────────────────────────
# The AMP-AD callsets keep their pgens in pgen/, WB-DWGS in joint_calls/.
DIR_WGS="${DATA}/amp-ad-genomics/WGS_Harmonization"
DIR_DC="${DATA}/amp-ad-genomics/DivCo_HS"
DIR_WB="${DATA}/amp-pd-genomics/WB-DWGS"
DIR_BR="${DATA}/amp-pd-genomics/BR-DSNWGS"        # postmortem AMP-PD, 97 donors

# RAW_* — step 1's INPUT. Step 1 appends _sexupd/_filtered itself, so pass RAW_*, never PF_*.
RAW_WGS="${DIR_WGS}/pgen/wgs_harm_hg38"
RAW_DC="${DIR_DC}/pgen/divco_hs_hg38"
RAW_WB="${DIR_WB}/joint_calls/all_chrs_merged"
RAW_BR="${DIR_BR}/pgen/br_dsnwgs_hg38"

# PF_* — step 1's filtered output, and step 2's input.
PF_WGS="${RAW_WGS}_filtered"
PF_DC="${RAW_DC}_filtered"
PF_WB="${RAW_WB}_filtered"
PF_BR="${RAW_BR}_filtered"

# Step 2 writes these; step 3 merges them.
NORM_WGS="${DIR_WGS}/pgen/wgs_harm_hg38_norm_bed"
NORM_DC="${DIR_DC}/pgen/divco_hs_hg38_norm_bed"
NORM_WB="${DIR_WB}/joint_calls/all_chrs_merged_norm_bed"
NORM_BR="${DIR_BR}/pgen/br_dsnwgs_hg38_norm_bed"

# Step 1 writes these; steps 4 and 6 read them back.
LBL_WGS="${DIR_WGS}/genotools/FILTERED.wgs_harm_ancestry_umap_linearsvc_predicted_labels.txt"
LBL_DC="${DIR_DC}/genotools/FILTERED.divco_hs_ancestry_umap_linearsvc_predicted_labels.txt"
LBL_WB="${DIR_WB}/genotools/FILTERED.wb_dwgs_ancestry_umap_linearsvc_predicted_labels.txt"
LBL_BR="${DIR_BR}/genotools/FILTERED.br_dsnwgs_ancestry_umap_linearsvc_predicted_labels.txt"

# ── reference data ───────────────────────────────────────────────────────────
# REF_DIR: fetched, large, not redistributable, not in the repo. SHIP_REF: small, versioned.
# How to obtain REF_DIR's contents: README.md, Reference data. The ancestry panel cannot be
# re-downloaded; check it with `md5sum -c ref/ref_panel.md5`.
REF_FASTA="${REF_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa.zst"
REF_PANEL="${REF_DIR}/ref_panel_gp2_prune_rm_underperform_pos_update"
REF_LABELS="${REF_DIR}/ref_panel_ancestry_updated.txt"
HIGHLD_BED="${SHIP_REF}/highld_exclude_hg38.bed"   # PCA input only — ships in ref/
REFFLAT="${SHIP_REF}/refFlat.txt"                  # locus coordinates for gene_annot.py

# ── what the clinical side hands the genotype steps ──────────────────────────
# §10  -> per-callset sex files,  read by step 1.
# §12a -> sample_annot.csv,       read by step 6's AF-concordance stage.
# §12  -> analysis_grain.csv       (§11-13 live in analysis_grain.py)
GRAIN="${GRAIN:-${CLINICAL_OUT}/analysis_grain.csv}"

# The grain's PC-free half (IID, source_callset, pheno, dx_detailed). Stage B needs only
# (callset, dx), so it can run before PCs exist — which is what keeps step 6 a single pass.
ANNOT="${ANNOT:-${CLINICAL_OUT}/sample_annot.csv}"

# Step 1's per-callset sex file:  SEX_FILE="$(sex_file wgs_harm)"
sex_file() { echo "${CLINICAL_OUT}/${1}_update_sex.txt"; }

# §13's phenotype/covariate/contrast files — step 7's inputs, and the sole definition of who is
# a case. Overridable as a set so a sensitivity run (e.g. EXCLUDE_DUAL, rerun to a different
# CLINICAL_OUT) points step 7 at a recorded alternative generation.
PHENO_SRC="${PHENO_SRC:-${CLINICAL_OUT}/pheno}"
COVAR_SRC="${COVAR_SRC:-${CLINICAL_OUT}/covar}"
CONTRASTS_CSV="${CONTRASTS_CSV:-${CLINICAL_OUT}/contrasts.csv}"

# Step 6 stage E writes it; the review tooling and step 9 read it. §12 does not (it globs the
# eigenvecs directly).
RETAINED_MANIFEST="${RETAINED_MANIFEST:-${MERGED_DIR}/by_ancestry_qc/retained_samples_manifest.csv}"

# ── modules ──────────────────────────────────────────────────────────────────
# plink1.9 is needed by step 3 (plink2's non-concatenating --pmerge-list is unimplemented) and
# by step 6 stage B (plink --assoc).
MOD_PLINK2="${MOD_PLINK2:-plink/6-alpha}"
MOD_PLINK1="${MOD_PLINK1:-plink/1.9.0-beta4.4}"
MOD_PYTHON="${MOD_PYTHON:-python/3.11}"
