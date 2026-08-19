# WGS Harmonization Pipeline

## Project Overview

This pipeline processes whole genome sequencing (WGS) data from multiple cohorts into a unified, harmonized GRCh38 format for ancestry prediction, sample and variant QC, and GWAS.

---

## Infrastructure

- **Cluster**: NIH Biowulf (SLURM)
- **Transfer node**: NIH Helix (used for Synapse downloads)
- **Tools**: plink2 v2.00a6LM, plink1.9, bcftools, liftOver, GenoTools v1.3.6, Synapse CLI
- **Python environment**: `.venv` at `/data/CARDPB2/sysbio/wgs/.venv` (Python 3.11, module `python/3.11`)
- **setuptools**: must be pinned to 67.8.0 for `pkg_resources` to work in sbatch

---

## Directory Structure

```
/data/CARDPB2/sysbio/wgs/
  bin/
    liftOver                                         ← UCSC liftOver binary
    hg19ToHg38.over.chain.gz                         ← liftover chain file
  data/
    ref/
      ref_panel_gp2_prune_rm_underperform_pos_update.{bed,bim,fam}
      ref_panel_ancestry_updated.txt
    plots/                                           ← interactive HTML PCA plots per dataset
    merged/
      cohort_merged.{bed,bim,fam}                   ← final merged cohort (3 datasets)
      merge_list.txt
      tmp_merge/
    amp-ad-genomics/
      WGS_Harmonization/
        joint_calls/                                 ← 24 b37 GATK scatter-interval VCFs + .tbi
        metadata/
          wgs_harm_hg38_update_sex.txt               ← plink sex update file
          MayoRNAseq_individual_metadata_harmonized.csv
          MSBB_individual_metadata_harmonized.csv
          ROSMAP_clinical_harmonized.csv
        pgen/
          intermediate/
            wgs_harm_filtered_b37.vcf.gz             ← filtered sorted b37 VCF
            wgs_harm_b37.{pgen,pvar,psam}            ← b37 pgen
          wgs_harm_hg38.{pgen,pvar,psam}             ← GRCh38 pgen
          wgs_harm_hg38_filtered_bed.{bed,bim,fam}   ← filtered bed (genotools input/byproduct)
        genotools/
          FILTERED.wgs_harm.json                     ← ancestry + QC results
          FILTERED.wgs_harm_ancestry_*.{pgen,psam,pvar,samples}
          FILTERED.wgs_harm_{ANCESTRY}.{pgen,psam,pvar}  ← per-ancestry QC'd pgens
      DivCo_HS/
        joint_calls/
          merged.deduped.vcf.gz                      ← full merged GRCh38 VCF (1,019 samples)
          merged.deduped.vcf.gz.tbi
        metadata/
          divco_hs_hg38_update_sex.txt               ← plink sex update file (3-path crosswalk)
          AMP-AD_DiverseCohorts_individual_metadata.csv
          AMP-AD_DiverseCohorts_biospecimen_metadata.csv
          AMP-AD_DiverseCohorts_assay_WGS_metadata.csv
        pgen/
          intermediate/
            divco_hs_filtered_hg38.vcf.gz            ← filtered VCF (no PASS filter — all '.')
          divco_hs_hg38.{pgen,pvar,psam}             ← GRCh38 pgen
          divco_hs_hg38_filtered_bed.{bed,bim,fam}   ← filtered bed (genotools input/byproduct)
        genotools/
          FILTERED.divco_hs.json                     ← ancestry + QC results
          FILTERED.divco_hs_ancestry_*.{pgen,psam,pvar,samples}
          FILTERED.divco_hs_{ANCESTRY}.{pgen,psam,pvar}
    amp-pd-genomics/
      WB-DWGS/
        joint_calls/
          all_chrs_merged.{pgen,pvar,psam}           ← original GRCh38 pgen (unfiltered)
          all_chrs_merged_filtered_bed.{bed,bim,fam} ← filtered bed (genotools input/byproduct)
        metadata/
          amp_pd_wbdwgs_update_sex.txt               ← plink sex update file
        genotools/
          FILTERED.wb_dwgs.json                      ← ancestry + QC results
          FILTERED.wb_dwgs_ancestry_*.{pgen,psam,pvar,samples}
          FILTERED.wb_dwgs_{ANCESTRY}.{pgen,psam,pvar}
      BR-DSNWGS/
        joint_calls/                                 ← NOT YET AVAILABLE
        metadata/
  scripts/
    README.md                                        ← this file
    archive/                                         ← superseded scripts
    logs/                                            ← all sbatch .o/.e files
    swarm/                                           ← swarm files and logs
```

---

## Datasets

| Dataset | Cohort | Samples | Build | Source Format | GenoTools | Merge |
|---|---|---|---|---|---|---|
| WGS_Harmonization | AMP-AD (MayoRNAseq, MSBB, ROSMAP) | 1,894 | b37 → GRCh38 | 24 GATK scatter VCFs | ✅ Complete | ✅ Included |
| WB-DWGS | AMP-PD | 10,418 | GRCh38 | Merged pgen | ✅ Complete | ✅ Included |
| DivCo_HS | AMP-AD diverse cohorts | 1,019 | GRCh38 | Single merged VCF | ✅ Complete | ✅ Included |
| BR-DSNWGS | AMP-PD | TBD | TBD | TBD | ⏳ Pending data | ⏳ Pending |

---

## Scripts

```
scripts/
  # ── Data acquisition ──────────────────────────────────────────────────────
  helix_synapse_download.sh              ← template for Synapse downloads on Helix
  synapse_download_sbatch.sh             ← alternative sbatch-based Synapse download

  # ── WGS_Harmonization (AMP-AD, b37 → GRCh38) ─────────────────────────────
  wgs_harm_s1_filter_concat.sh           ← sbatch: bcftools concat + filter + sort → b37 VCF
  wgs_harm_s2_b37_to_pgen.sh             ← sbatch: plink2 b37 VCF → b37 pgen
  wgs_harm_s3_liftover.sh                ← sbatch: liftOver BED approach → GRCh38 pgen

  # ── DivCo_HS (AMP-AD, already GRCh38) ────────────────────────────────────
  divco_hs_s1_filter.sh                  ← sbatch: bcftools filter + annotate → filtered VCF
  divco_hs_s2_vcf_to_pgen.sh             ← sbatch: plink2 filtered VCF → GRCh38 pgen

  # ── Shared / reusable ─────────────────────────────────────────────────────
  liftover_pgen.sh                       ← reusable liftover worker (called by wgs_harm_s3)
  update_sex.py                          ← Python: join pgen psam with harmonized metadata,
                                            run plink2 --update-sex for each dataset;
                                            sex update files stored in each dataset's metadata/ dir
  run_update_sex.sh                      ← sbatch wrapper for update_sex.py
  run_genotools.sh                       ← reusable sbatch: sex update → filter → bed → GenoTools
                                            submit with --export=PGEN=...,SEX_FILE=...,OUT_DIR=...,DATASET=...
                                            also produces filtered bed as byproduct used in merge

  # ── Merge ─────────────────────────────────────────────────────────────────
  wgs_merge_s1_merge_bed.sh              ← sbatch: plink1.9 union merge of all 3 filtered beds

  # ── Analysis ──────────────────────────────────────────────────────────────
  plot_ancestry_pca.py                   ← plotly 3D PCA plots from genotools JSON → HTML
```

---

## Pipeline Order

### WGS_Harmonization

> **Key context**: The 24 VCFs are GATK scatter-interval files, NOT chromosome-split. Each file contains variants from multiple chromosomes. The pipeline processes them as a genome-wide dataset.

```bash
# 1. Download data (run on Helix in tmux)
synapse get -r syn11707420 --downloadLocation .../WGS_Harmonization/joint_calls/

# 2. Reindex VCFs (stale .tbi files from download)
bash scripts/wgs_harm_reindex_vcfs.sh   # generates + submits as swarm

# 3. Filter + concat all 24 scatter VCFs into single sorted b37 VCF
JID1=$(sbatch scripts/wgs_harm_s1_filter_concat.sh | awk '{print $NF}')

# 4. Convert filtered b37 VCF → b37 pgen
JID2=$(sbatch --dependency=afterok:${JID1} scripts/wgs_harm_s2_b37_to_pgen.sh | awk '{print $NF}')

# 5. Liftover b37 pgen → GRCh38 pgen
JID3=$(sbatch --dependency=afterok:${JID2} scripts/wgs_harm_s3_liftover.sh | awk '{print $NF}')

# 6. Run GenoTools (sex update + filter + bed + ancestry/QC)
sbatch \
  --job-name=genotools_wgs_harm \
  --output=scripts/logs/genotools_wgs_harm.o \
  --error=scripts/logs/genotools_wgs_harm.e \
  --export=PGEN=/data/CARDPB2/sysbio/wgs/data/amp-ad-genomics/WGS_Harmonization/pgen/wgs_harm_hg38,SEX_FILE=/data/CARDPB2/sysbio/wgs/data/amp-ad-genomics/WGS_Harmonization/metadata/wgs_harm_hg38_update_sex.txt,OUT_DIR=/data/CARDPB2/sysbio/wgs/data/amp-ad-genomics/WGS_Harmonization/genotools,DATASET=wgs_harm \
  scripts/run_genotools.sh
# Note: produces wgs_harm_hg38_filtered_bed.{bed,bim,fam} as byproduct used in merge
```

### WB-DWGS

> **Key context**: Already GRCh38 merged pgen. GenoTools handles sex update + filtering + bed conversion in one job.

```bash
# 1. Run GenoTools (sex update + filter + bed + ancestry/QC)
sbatch \
  --job-name=genotools_wb_dwgs \
  --output=scripts/logs/genotools_wb_dwgs.o \
  --error=scripts/logs/genotools_wb_dwgs.e \
  --export=PGEN=/data/CARDPB2/sysbio/wgs/data/amp-pd-genomics/WB-DWGS/joint_calls/all_chrs_merged,SEX_FILE=/data/CARDPB2/sysbio/wgs/data/amp-pd-genomics/WB-DWGS/metadata/amp_pd_wbdwgs_update_sex.txt,OUT_DIR=/data/CARDPB2/sysbio/wgs/data/amp-pd-genomics/WB-DWGS/genotools,DATASET=wb_dwgs \
  scripts/run_genotools.sh
# Note: produces all_chrs_merged_filtered_bed.{bed,bim,fam} as byproduct used in merge
```

### DivCo_HS

> **Key context**: Already GRCh38 and pre-merged. FILTER column is all '.' (no PASS flag applied). Sex crosswalk required 3-path ID resolution — see sex crosswalk note below.

```bash
# 1. Download data (run on Helix in tmux)
synapse get -r syn68259948 --downloadLocation .../DivCo_HS/joint_calls/

# 2. Delete scatter-interval VCFs (only merged.deduped.vcf.gz needed)
rm .../DivCo_HS/joint_calls/chr*.vcf.gz{,.tbi}

# 3. Filter + annotate IDs → filtered GRCh38 VCF (no --apply-filters PASS)
JID1=$(sbatch scripts/divco_hs_s1_filter.sh | awk '{print $NF}')

# 4. Convert filtered VCF → GRCh38 pgen
JID2=$(sbatch --dependency=afterok:${JID1} scripts/divco_hs_s2_vcf_to_pgen.sh | awk '{print $NF}')

# 5. Run GenoTools (sex update + filter + bed + ancestry/QC)
sbatch \
  --job-name=genotools_divco_hs \
  --output=scripts/logs/genotools_divco_hs.o \
  --error=scripts/logs/genotools_divco_hs.e \
  --export=PGEN=/data/CARDPB2/sysbio/wgs/data/amp-ad-genomics/DivCo_HS/pgen/divco_hs_hg38,SEX_FILE=/data/CARDPB2/sysbio/wgs/data/amp-ad-genomics/DivCo_HS/metadata/divco_hs_hg38_update_sex.txt,OUT_DIR=/data/CARDPB2/sysbio/wgs/data/amp-ad-genomics/DivCo_HS/genotools,DATASET=divco_hs \
  scripts/run_genotools.sh
# Note: produces divco_hs_hg38_filtered_bed.{bed,bim,fam} as byproduct used in merge
```

### Merge (all datasets)

> **Key context**: Filtered bed files are a byproduct of each dataset's genotools run. No separate conversion step needed. plink1.9 used for union merge (plink2 non-concatenating merge not implemented in v2.00a6LM).

```bash
# Merge all three filtered bed files (union — missing genotypes for absent variants)
sbatch scripts/wgs_merge_s1_merge_bed.sh
# Output: data/merged/cohort_merged.{bed,bim,fam}
```

### Downstream Plan

Once merge completes and BR-DSNWGS becomes available:

1. **Cross-dataset relatedness** — run genotools `--related` on merged dataset to identify duplicates and cryptic relatives across cohorts
2. **Split by ancestry** — use genotools `ancestry_labels` from individual dataset JSON files
3. **Exclude QC fails** — apply per-ancestry genotools QC outputs (`{ANCESTRY}_pass_fail`) to merged dataset
4. **Add covariates** — PCs from merged dataset, sex, age, APOE, dataset batch
5. **GWAS** — per-ancestry, case/control assignment from harmonized metadata

### BR-DSNWGS

> Not yet available. Pipeline TBD based on format and build when data arrives.

---

## Sex Crosswalk Notes

### DivCo_HS (3-path crosswalk)
psam sample IDs use three different formats requiring separate resolution paths:
- **Direct match** (620 samples): numeric psam IID matches `individualID` in individual metadata directly
- **Suffix stripping** (293 samples): specimen IDs like `<individualID>_DLPFC_WGS` — strip `_DLPFC_WGS` suffix to recover numeric `individualID`
- **Biospecimen lookup** (106 samples): `-D` suffixed specimen IDs like `<specimenID>-D` — match via `specimenID` in biospecimen metadata to get numeric `individualID`

Result: 1019/1019 (100%) matched. Script: `update_sex.py` (see DivCo_HS section).

### WGS_Harmonization
Sex file at `metadata/wgs_harm_hg38_update_sex.txt`. 1886/1894 (99.6%) matched; 8 samples remain unknown sex.

### WB-DWGS
Sex file at `metadata/amp_pd_wbdwgs_update_sex.txt`. 10418/10418 (100%) matched.

---

## Key Gotchas

- **GATK scatter VCFs are not chromosome-split** — WGS_Harmonization VCFs each contain variants from multiple chromosomes. Never assume file number = chromosome number.
- **`#` in swarm files is always a comment** — Biowulf's swarm executor strips `#` even inside quotes. Never use `--set-missing-var-ids @:#` in a swarm file; use `bcftools annotate --set-id` instead. In sbatch scripts, `#` works fine.
- **plink2 v2.00a6LM `--pmerge-list` non-concatenating merge not implemented** — use plink1.9 for multi-cohort merging.
- **PAR1/PAR2 chromosome codes not recognized by plink1.9** — always use `--merge-par` when converting pgen → bed if `--split-par` was used during VCF import.
- **chrX requires `--split-par` + `--update-sex`** — use `hg19` for b37 data, `hg38` for GRCh38. Sex file must be 2-column (IID SEX) if psam has no FID column.
- **liftOver requires `chr`-prefixed chromosome names** — b37 uses bare numbers; add `chr` prefix when building BED from pvar.
- **`set -o pipefail` must come after all `#SBATCH` directives** — SLURM stops parsing `#SBATCH` lines at the first non-comment line. Placing `set -o pipefail` before `#SBATCH` silently drops all resource requests.
- **DivCo_HS FILTER column is all `.`** — no VQSR was applied; omit `--apply-filters PASS` or all variants will be dropped.
- **venv Python must use `module load python/3.11`** — the default system Python (`/usr/local/bin/python`) points to an Anaconda env on compute nodes that breaks `pkg_resources`. Always load `python/3.11` before activating `.venv` in sbatch scripts.
- **setuptools must be pinned to 67.8.0** — newer versions break `pkg_resources` import needed by genotools. Fix: `pip install setuptools==67.8.0`.
- **`subprocess.run()` with `capture_output=True` hides plink2 output** — use without `capture_output` or add `stderr=subprocess.STDOUT` so output appears in sbatch logs.
- **Always use `--dependency=afterok` when chaining sbatch jobs** — avoids downstream jobs running on incomplete input.
- **Synapse downloads may be interrupted** — always rerun `synapse get -r` to resume; it skips completed files and resumes incomplete ones.
- **tabix `.tbi` files may be stale after download** — use `tabix -f -p vcf` to force reindex before running bcftools with `--regions`.