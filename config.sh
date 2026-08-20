# ad-pd-gwas — configuration.
#
# Contains NO machine-specific paths. Everything below is derived from this file's own
# location, so the project runs wherever it is copied to and there is nothing to edit
# when it moves. Every script sources this; nothing else resolves a location.
#
# Sourced, never executed. Every value can be overridden from the environment:
#   GRAIN=/some/other/grain.csv ./submit.sh scripts/07_gwas.sh

# ── the one root ─────────────────────────────────────────────────────────────
# Code and data share a root: this file sits at the top of the project, and data/ sits
# beside it. On the cluster that is /data/CARDPB2/sysbio/wgs; on a laptop it is wherever
# the repo was cloned. `data/` may be a symlink to a larger allocation — paths below
# resolve through it either way.
#
# WGS_ROOT is kept as an exported alias because the shell steps `cd "${WGS_ROOT}"` and the
# two Python helpers accept it as an override. They self-locate identically by default, so
# nothing breaks if it is unset.
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export WGS_ROOT="${WGS_ROOT:-${PROJECT_ROOT}}"

DATA="${PROJECT_ROOT}/data"
REF_DIR="${DATA}/ref"                  # fetched, large, gitignored — see REF_* below
SHIP_REF="${PROJECT_ROOT}/ref"         # ships with the code, versioned, small
MERGED_DIR="${DATA}/merged"
VENV="${PROJECT_ROOT}/.venv"

# ── the phenotype boundary ───────────────────────────────────────────────────
# clinical_core.py writes everything the genotype steps need into ONE directory, and the
# steps read it from there. Nothing is copied into data/ — a handoff file has exactly one
# location, so it cannot go stale against the run that produced it.
CLINICAL_OUT="${CLINICAL_OUT:-${PROJECT_ROOT}/clinical_core_out}"

# ── the three source callsets ────────────────────────────────────────────────
# Note the asymmetry: the two AMP-AD callsets keep their pgens in pgen/, WB-DWGS in joint_calls/.
DIR_WGS="${DATA}/amp-ad-genomics/WGS_Harmonization"
DIR_DC="${DATA}/amp-ad-genomics/DivCo_HS"
DIR_WB="${DATA}/amp-pd-genomics/WB-DWGS"
DIR_BR="${DATA}/amp-pd-genomics/BR-DSNWGS"        # postmortem AMP-PD, 97 donors

# RAW_* — step 0's output, and step 1's INPUT. Step 1 derives ${PGEN}_sexupd, ${PGEN}_filtered
# and ${PGEN}_filtered_bed from it, so passing a stem that already ends in _filtered would
# produce _filtered_filtered. Pass RAW_*, not PF_*, to step 1.
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
# Two homes, along the same line as everything else: REF_DIR is fetched and large and is
# not in the repo (see README §2 — none of it is redistributable); SHIP_REF is small,
# versioned, and travels with the code. Nothing is copied between them.
REF_FASTA="${REF_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa.zst"
REF_PANEL="${REF_DIR}/ref_panel_gp2_prune_rm_underperform_pos_update"
REF_LABELS="${REF_DIR}/ref_panel_ancestry_updated.txt"
HIGHLD_BED="${SHIP_REF}/highld_exclude_hg38.bed"   # PCA input only — ships in ref/
REFFLAT="${SHIP_REF}/refFlat.txt"                  # locus coordinates for gene_annot.py

# ── what clinical_core.py hands the genotype steps ───────────────────────────
# §10  -> per-callset sex files,  read by step 1.
# §12a -> sample_annot.csv,       read by step 6's AF-concordance stage.
# §12  -> analysis_grain.csv,     read by step 7.  (§11-13 live in analysis_grain.py)
# Overridable so a sensitivity run can point step 7 at an alternative grain.
GRAIN="${GRAIN:-${CLINICAL_OUT}/analysis_grain.csv}"

# sample_annot.csv is the grain's PC-FREE half: IID, source_callset, pheno, dx_detailed.
# It exists because the AF-concordance stage needs (callset, dx) per sample and NOTHING else —
# it never reads a PC — while the grain cannot be built until step 6 has produced PCs. Splitting
# the file splits the dependency, which is what lets step 6 run as ONE pass instead of two.
# §12a writes it from §7's crosswalk and §4's reconciliation, so it is available before step 1.
ANNOT="${ANNOT:-${CLINICAL_OUT}/sample_annot.csv}"

# Step 1 takes SEX_FILE per callset; this is the naming rule it follows.
#   SEX_FILE="$(sex_file wgs_harm)"
sex_file() { echo "${CLINICAL_OUT}/${1}_update_sex.txt"; }

# NOTE: §13 also writes clinical_core_out/{pheno,covar}/, but step 7 does NOT read them —
# it rebuilds both from ${GRAIN} in awk, and defines its own PHENO_DIR/COVAR_DIR under its
# output dir. Two implementations of the same files; the awk one is what actually runs.
# Deliberately not named here so nothing shadows step 7's. See README "known duplication".

# Written by the genotype side (step 6 stage E), read by the review tooling. NOT read by §12,
# which globs the eigenvecs directly — see the docstring in scripts/ancestry_qc_manifest.py.
RETAINED_MANIFEST="${RETAINED_MANIFEST:-${MERGED_DIR}/by_ancestry_qc/retained_samples_manifest.csv}"

# ── modules ──────────────────────────────────────────────────────────────────
# plink2's non-concatenating --pmerge-list is unimplemented in this build, which is why
# step 3 alone drops to plink1.9.
MOD_PLINK2="${MOD_PLINK2:-plink/6-alpha}"
MOD_PLINK1="${MOD_PLINK1:-plink/1.9.0-beta4.4}"
MOD_PYTHON="${MOD_PYTHON:-python/3.11}"
