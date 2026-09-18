# WGS Harmonization Pipeline

## Project Overview

This pipeline processes whole genome sequencing (WGS) data from multiple cohorts into a unified,
harmonized GRCh38 format for ancestry prediction, sample and variant QC, and GWAS. The analysis it
feeds is **AD vs PD across AMP-AD and AMP-PD, with AD defined by neuropathology**.

**This file is the reference: what the tools are, where files live, in what order they run, and
which gotchas bite.** It does not carry project state or scientific rationale — those have their
own homes. One job each:

| doc | answers |
|---|---|
| `README.md` (this file) | what the pieces are, where files live, and **in what order things run** |
| `METHODS.md` | **why** each choice was made — draft methods section for the paper |
| `HANDOFF.md` | what is true right now: status and open issues |
| `PROJECT_LOG.md` | what we did and what was ruled out. Append-only; grep it before investigating |
| `RUNLOG.md` | which jobs ran and how they ended. Generated: `bash scripts/runlog.sh --md > RUNLOG.md` |
| `wgs_core.ipynb` | the **executable** run order — it runs on biowulf and submits every step |

**One root.** Code and data share a root — on biowulf `/data/CARDPB2/sysbio/wgs/amp-ad-pd-wgs-gwas`, on a laptop
wherever the repo was cloned. Nothing in the project holds an absolute path: `config.sh` derives
the root from its own location, and `clinical_common.py` does the same. To move it, copy the folder.

**Code reaches the cluster by `git pull`.** There is a remote as of 2026-08-21. Do not rsync code:
`rsync` without `--delete` cannot express a deletion, which caused three stale-artifact bugs.

---

## Infrastructure

- **Cluster**: NIH Biowulf (SLURM)
- **Transfer node**: NIH Helix (used for Synapse downloads)
- **Tools**: plink2 v2.00a6LM, plink1.9, bcftools, liftOver, GenoTools v1.3.6, Synapse CLI
- **Python environment**: `.venv` at `/data/CARDPB2/sysbio/wgs/amp-ad-pd-wgs-gwas/.venv` (Python 3.11, module `python/3.11`)
- **setuptools**: must be pinned to 67.8.0 for `pkg_resources` to work in sbatch
- **`requirements.lock.txt`**: 115 pinned packages, frozen from the working cluster venv
  2026-09-17. **It does not contain `setuptools`** — `pip freeze` omits it — so rebuilding is
  `python3 -m venv .venv && pip install setuptools==67.8.0 && pip install -r requirements.lock.txt`,
  in that order. Before this file existed the environment lived only as installed bytes on one
  disk; a venv cannot be moved, so that made the project root effectively immovable.

---

## Directory Structure

```
/data/CARDPB2/sysbio/wgs/amp-ad-pd-wgs-gwas/
  bin/
    liftOver                                         ← UCSC liftOver binary
    hg19ToHg38.over.chain.gz                         ← liftover chain file
  data/
    ref/
      ref_panel_gp2_prune_rm_underperform_pos_update.{bed,bim,fam}
      ref_panel_ancestry_updated.txt
    plots/                                           ← interactive HTML PCA plots per dataset
    merged/
      cohort_merged.{bed,bim,fam}                   ← merged cohort, ALL FOUR callsets
      merge_list.txt
      tmp_merge/
      relatedness/
        retained_manifest.csv                       ← step 5 output; steps 6 and §11-12 read it
        exclude_reasons.tsv
      exclude_af_concordance.txt                    ← step 6 stage B builds it, stage C applies it
      af_concordance/                               ← stage B scratch (per-cell .afreq, .snplist)
      by_ancestry_qc/                               ← STEP 6 OUTPUT
        unfiltered/                                 ← stage A: the BASELINE, kept permanently
          cohort_{ANC}_qc.{bed,bim,fam}
          cohort_{ANC}_pca.{eigenvec,eigenval}
          retained_samples_manifest.csv             ← "before" for the AF-filter comparison
        cohort_{ANC}_qc.{bed,bim,fam}               ← stage C: the ASSOCIATION set (step 7 reads)
        cohort_{ANC}_pca.{eigenvec,eigenval}        ← stage D: the COVARIATES (§12 reads)
        retained_samples_manifest.csv               ← "after"
        step6_summary.txt
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
        joint_calls/
          AMPPD_postmortem_joint_gt_call_97donors.vcf.gz   ← 97 postmortem donors
        metadata/
        pgen/
          br_dsnwgs_hg38.{pgen,pvar,psam}            ← GRCh38 pgen (see --sort-vars gotcha)
        genotools/
  clinical_core_out/                                 ← THE PHENOTYPE BOUNDARY. One location each;
    {callset}_update_sex.txt                         ←   nothing is copied into data/, so a handoff
    sample_annot.csv                                 ←   file cannot go stale against its producer.
    individual_core.csv                              ← §9 audit table -> analysis_grain.py
    genome_crosswalk.csv                             ← §9 audit table -> analysis_grain.py
    analysis_grain.csv                               ← the grain. Step 7 does NOT read it directly
    pheno/  covar/  contrasts.csv                    ←   — it reads THESE (§13), as of 2026-08-20
    qc_outcomes.csv
  results/                                           ← plots + eta^2 evidence (mostly gitignored)
  scripts/                                           ← see Scripts below
    logs/                                            ← all sbatch .o/.e files
  ref/                                               ← ships with the code: highld BED, refFlat
  README.md   METHODS.md   HANDOFF.md   PROJECT_LOG.md   RUNLOG.md
  config.sh   submit.sh
  clinical_common.py   clinical_core.py   analysis_grain.py
  wgs_core.ipynb                                     ← the orchestrator; runs ON biowulf
```

---

## Datasets

| Dataset | Cohort | Samples | Build | Source Format | GenoTools | Merge |
|---|---|---|---|---|---|---|
| WGS_Harmonization | AMP-AD (MayoRNAseq, MSBB, ROSMAP) | 1,894 | b37 → GRCh38 | 24 GATK scatter VCFs | ✅ Complete | ✅ Included |
| WB-DWGS | AMP-PD | 10,418 | GRCh38 | Merged pgen | ✅ Complete | ✅ Included |
| DivCo_HS | AMP-AD diverse cohorts | 1,019 | GRCh38 | Single merged VCF | ✅ Complete | ✅ Included |
| BR-DSNWGS | AMP-PD (postmortem) | 97 | GRCh38 | Single joint-call VCF | ✅ Complete | ✅ Included |

Merged cohort: **172,497,055 variants × 13,334 samples**; 12,495 retained after step 5.
Per-step state lives in `HANDOFF.md`, not here.

---

## Scripts

Every step is a numbered script in `scripts/`, submitted through `submit.sh`, which resolves the
root and sources `config.sh`. The scripts named in older revisions of this file
(`run_genotools.sh`, `wgs_merge_s1_merge_bed.sh`, `update_sex.py`, the per-callset `*_s1_*.sh`
pairs) were consolidated into these and no longer exist.

```
scripts/
  # ── the pipeline, in order ────────────────────────────────────────────────
  01_genotools.sh              ← per callset: sex update -> filter -> bed -> ancestry + QC
  02_normalize.sh              ← per callset: chr naming + variant IDs, to a common convention
  03_merge.sh                  ← plink1.9 union merge of all four -> cohort_merged
  04_relatedness.sh            ← KING, report-only
  05_excludelist.{py,sh}       ← duplicates, relatives, QC fails -> retained_manifest.csv
  06_ancestry_qc.sh            ← per-ancestry variant QC + PCA. FIVE STAGES, ONE SUBMISSION:
                                   A unfiltered QC+PCA   B build AF exclusion list from A
                                   C apply it            D filtered QC+PCA
                                   E both sample manifests
  07_gwas.sh                   ← plink2 --glm per ancestry per contrast. Reads §13's pheno/covar/
                                 contrasts.csv — it does NOT rebuild arms from the grain.
  08_ctrl_ctrl_filter.{py,sh}  ← control-vs-control artifact scan. ANNOTATES, never subtracts.

  # ── tier 1b, release. After step 8; not part of every run ─────────────────
  09_amppd_subset.sh           ← subset the AMP-PD donors out of step 6 stage C into a staging
                                 tree. sbatch on BIOWULF. No network code in it at all.
  10_amppd_push.sh             ← publish that tree to GCS. Runs on HELIX, not under sbatch.

  # ── called by step 6, not submitted directly ──────────────────────────────
  af_concordance_build.{py,sh} ← stage B. The .sh is a wrapper for re-tuning knobs only.
  ancestry_qc_manifest.py      ← stage E. retained_samples_manifest.csv, once per generation.

  # ── tier 2, preflight (see Run order) ─────────────────────────────────────
  00_setup_env.sh              ← venv + GenoTools install
  diag_order.py                ← variant order vs the reference panel. MANDATORY on any callset
                                 not already cleared, before step 1. Read-only.
  hooks-pre-commit             ← install into .git/hooks/ after a clone

  # ── tier 3, on demand ─────────────────────────────────────────────────────
  gene_annot.py                ← what gene is this variant in; locus coordinates from ref/refFlat.txt
  runlog.sh                    ← regenerates RUNLOG.md from the SLURM accounting log
  nb_guard.py                  ← refuses to commit a notebook carrying stored outputs

  # ── called by step 1, not run directly ────────────────────────────────────
  genotools_capped.py          ← GenoTools entry point that respects the SLURM allocation

review/                        ← figures. Run on the cluster against outputs in place (the
                                 sumstats are ~23 GB); pull the PNGs down, not the inputs.
  plot_af_filter_effect.py     ← the AF-filter before/after proof figure + eta^2 tables
  plot_pcs_by_callset.py       ← PC scatter by callset + eta^2 per PC per stratum
  plot_gwas.py                 ← QQ + Manhattan, with lambda_GC cross-checked against step 7's
  drift_check.py               ← does the working tree differ from a revision in anything that
                                 RUNS? AST for .py, quote-aware comment strip for .sh, so a
                                 rewritten docstring or realigned trailing comment reads as
                                 identical. Read-only. `python3 review/drift_check.py [<rev>]`
  methods_numbers.py           ← re-derives every number in METHODS.md from the artifacts and
                                 diffs it against the doc. Read-only. Run it before any write-up
                                 edit; --strict to fail on anything it could not check.
                                 --discordance measures §6.4's duplicate-pair enrichment (cluster).
```

The clinical side is three files at the project root:

| file | runs | writes |
|---|---|---|
| `clinical_common.py` | imported, never run | — paths, readers, reconciliation rules |
| `clinical_core.py` | **once, before step 1** | sex files, audit tables, `sample_annot.csv` |
| `analysis_grain.py` | **once, after step 6** | `analysis_grain.csv`, pheno/covar, contrasts |

Two runs, not one, and that is irreducible: step 1 applies the sex files, and the grain carries
step 6's PCs because step 7 reads them as covariates. They are two *files* rather than one script
run twice so that nothing re-derives, and so the second half cannot rewrite the sex files step 1
already consumed.

---

## Run order

Everything is grouped by **when it runs and what triggers it**, not by what it is allowed to do.
The distinction matters: `diag_order.py` used to be filed under "read-only diagnostics" alongside
postmortem tools, which made a mandatory gate look optional — and skipping it is exactly how all 97
BR-DSNWGS samples came back mislabelled with a clean exit and 0.97 nominal model accuracy.

`wgs_core.ipynb` is the executable copy of tier 1 and runs on biowulf. This section is the map.

### Tier 1 — the pipeline. Runs every time, in this order.

```
clinical_core.py                §1-10 sex files, §12a sample_annot.csv   (cluster only)
  01 genotools   (x4)           sex update -> filter -> bed -> ancestry + QC
  02 normalize   (x4)           chr naming + variant IDs to a common convention
  03 merge                      plink1.9 union of all four -> cohort_merged
  04 relatedness                KING, report-only
  05 excludelist                duplicates, relatives, QC fails -> retained_manifest.csv
  06 ancestry QC                ONE submission, five stages:
                                  A unfiltered QC+PCA   B build the AF exclusion list from A
                                  C apply it            D filtered QC+PCA
                                  E both sample manifests
analysis_grain.py               §11-13 on step 6's PCs                   (cluster only)
  07 gwas                       plink2 --glm per ancestry per contrast
  review/plot_gwas.py           QQ + Manhattan
  08 ctrl-vs-ctrl               annotates; never subtracts
```

Two hard prerequisites, both of which fail loudly rather than silently:

- **`analysis_grain.py` must run before any step-7 run.** A `contrasts.csv` predating 2026-08-21
  carries no callset-skew columns and step 7 exits 3 naming the missing one. Rerun the grain; do
  not patch the CSV.
- **`clinical_core.py` and `analysis_grain.py` run on the CLUSTER ONLY.** The two machines hold
  different clinical inputs and a laptop run silently produces a smaller, wrong table. See
  `HANDOFF.md` open issue 1.

`review/plot_af_filter_effect.py` runs after step 6, off the two manifests that step wrote.

### Tier 1b — release. Runs after step 8, when data is being handed to someone else.

```
  09 amppd subset               sbatch on BIOWULF  — subset + manifest. Never touches the network
  10 amppd push                 run on HELIX       — upload + verify
```

Not part of every run, and deliberately not renumbered into tier 1: it consumes the pipeline's
output and produces nothing the pipeline reads back. It is listed here rather than under tier 3
because it is a **step in the data process with a fixed position** — after step 8, on artifacts
step 6 wrote — not a tool you reach for at an arbitrary moment.

| | |
|---|---|
| what it ships | `data/merged/by_ancestry_qc/cohort_<ANC>_qc` — step 6 **stage C**, the exact variant set step 7 tested — restricted to `source_callset` in {`wb_dwgs`, `br_dsnwgs`} |
| shape | one fileset **per ancestry stratum**, because the QC is per stratum. Concatenating them would invent a variant set no GWAS in this study ran on |
| who is AMP-PD | read by header name from step 6 stage E's `retained_samples_manifest.csv`. **Not** re-derived from the four genotools label files — `ancestry_qc_manifest.py` is the one implementation of callset membership |
| not re-filtered | dropping the AMP-AD donors leaves some variants monomorphic within the subset; they stay. Re-applying `--maf` would produce a fileset that disagrees with the published sumstats. `SUBSET_MAF`/`SUBSET_GENO` opt in, and the release README says which was used |
| destination | **no bucket is hardcoded.** `GCS_DEST` is required at push time |

**Two scripts on two hosts, and the split is not cosmetic.** Biowulf compute nodes have no
general outbound network, so a push from one fails looking like an auth problem; step 10 refuses
to run inside a SLURM allocation. Keeping them separate also means step 9 has no gcloud
dependency — it records size and SHA-256 with coreutils, which is the better pair to ship anyway
since a consumer can check SHA-256 after download. CRC32C is what GCS stores, so step 10 computes
it at upload time and compares against the bucket.

```bash
# on biowulf
./submit.sh scripts/09_amppd_subset.sh

# then on helix.nih.gov, after step 9 reports "reconciled"
GCS_DEST=gs://<bucket>/<prefix> bash scripts/10_amppd_push.sh
```

**Three things it refuses to do, none of them overridable by a flag:**

- Ship a **partial** release. Every AMP-PD sample in the retained manifest must land in exactly
  one released fileset; the reconciliation is fatal, and it names the strata that came up short.
  A stratum with AMP-PD donors but no step-6 fileset is a hole, and it looks identical to a
  clean run unless something counts.
- Upload to a bucket whose metadata **reads fine and says it is public** — uniform bucket-level
  access off, public access prevention not enforced, or an `allUsers` binding. Those reads need
  `storage.buckets.get`/`.getIamPolicy`, which `roles/storage.objectAdmin` does not grant, so on
  a program-managed bucket the check often cannot run at all; that case is reported loudly and
  the upload proceeds on the operator's say-so. `SKIP_BUCKET_CHECK=1` skips it outright. A bucket
  that can be neither described nor listed is a wrong name and still exits.
- Report success on an **unverified** upload. `cp` exiting 0 does not prove every object arrived
  — the failure it misses is a file that was never in the argument list. Stage E re-lists the
  bucket and compares presence, size and CRC32C against `MANIFEST.tsv`.

Not MD5 anywhere: GCS computes no MD5 for composite (parallel-chunked) uploads, which is why 42
of the 46 objects in the sumstats release came back MD5-less (`PROJECT_LOG.md` 2026-09-15).

### Tier 2 — preflight. Runs when a precondition changed, not every time.

| check | trigger | if it fails |
|---|---|---|
| `bash scripts/00_setup_env.sh` | new machine, or a fresh clone with no `.venv` | nothing else will run |
| `cp scripts/hooks-pre-commit .git/hooks/pre-commit && chmod +x .git/hooks/pre-commit` | fresh clone | notebook outputs can reach git; one already carried a subject-level table |
| `python3 scripts/diag_order.py <panel>.bim <callset>.pvar` | **any callset not already cleared, BEFORE step 1** | non-zero rank drops -> add `--sort-vars` at import, or ancestry output is garbage that exits 0 |

`wgs_harm`, `divco_hs`, `wb_dwgs` and `br_dsnwgs` have all been cleared. A fifth callset has not.

### Tier 3 — on demand. No position in the sequence.

| tool | use |
|---|---|
| `python3 scripts/gene_annot.py --at chr19:44908684` | what gene does this variant sit in — this is what named APOE and HLA-DQB1 in step 8's results |
| `python3 scripts/gene_annot.py CR1 SNCA LRRK2` | gene coordinates, optionally `--flank` |
| `bash scripts/runlog.sh --md > RUNLOG.md` | regenerate the job table from SLURM accounting |
| `python3 scripts/nb_guard.py --staged` | what the pre-commit hook runs; `--strip <nb>` to clean one |

---

## Sex Crosswalk Notes

### DivCo_HS (3-path crosswalk)
psam sample IDs use three different formats requiring separate resolution paths:
- **Direct match** (620 samples): numeric psam IID matches `individualID` in individual metadata directly
- **Suffix stripping** (293 samples): specimen IDs like `<individualID>_DLPFC_WGS` — strip `_DLPFC_WGS` suffix to recover numeric `individualID`
- **Biospecimen lookup** (106 samples): `-D` suffixed specimen IDs (`<specimenID>-D`) — match via `specimenID` in biospecimen metadata to get numeric `individualID`

Result: 1019/1019 (100%) matched. The crosswalk now lives in `clinical_core.py` §7, not the
retired `update_sex.py`.

### WGS_Harmonization
1886/1894 (99.6%) matched; 8 samples remain unknown sex.

### WB-DWGS
10418/10418 (100%) matched.

### BR-DSNWGS
97/97 resolved via the AMP-PD sample inventory (`clinical_core.py` §7). 60M / 37F.

### Where the sex files actually live — NOT in `metadata/`
All four are written by `clinical_core.py` §10 to **`clinical_core_out/{callset}_update_sex.txt`**,
and step 1 reads them from there directly (`config.sh: sex_file`). Nothing is copied into each
callset's `metadata/`. That is the point: one location per handoff file, so it cannot go stale
against the run that produced it. Earlier revisions of this file described the `metadata/` layout,
and those paths are gone.

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

### Learned since (each of these cost a run)

- **GenoTools aligns to the reference panel by COLUMN POSITION, not by variant name.** BR-DSNWGS's
  pgen was ordered `1,10,11,…,19,2,20,…` (per-chromosome files concatenated alphabetically) while
  the panel is numeric, so every shared column was projected as a *different* variant and all 97
  samples came back CAH. Fixed with `--sort-vars` in step 1. **Run
  `python3 scripts/diag_order.py <panel>.bim <callset>.pvar` on any new callset before trusting its
  ancestry output** — the other three passed only because their files happened to be numeric.
- **GenoTools sizes its worker pool from the NODE, not the SLURM allocation.** `os.cpu_count()` on
  a shared node exceeded biowulf's per-user `ulimit -u` of 1024; `pthread_create` returned EAGAIN
  and SLURM logged a bare `ExitCode 1:0`. Needs both `scripts/genotools_capped.py` **and**
  `GENOTOOLS_MAX_WORKERS=16` — the allocation alone is not enough. Pinning `OMP_/NUMBA_NUM_THREADS`
  is *not* the fix; joblib already does that, and the worker count is the lever.
- **Step 4's `COMMON_GENO` must stay below the smallest callset's share of the cohort.** BR is
  97/13,334 = 0.0073, so at `--geno 0.05` every variant BR lacked stayed in the common set and all
  97 BR samples read as 50% missing — which step 5 then consumed as a duplicate tie-break signal.
  `0.005` fixes it. Re-check if a callset under ~0.5% is ever added.
- **A file COUNT is not a file LIST.** The laptop/cluster clinical-input check compared counts,
  passed, and was wrong: `WGS_Harmonization/metadata` had 5 files on both sides but different ones.
  Compare names.
- **`rsync` without `--delete` cannot express a deletion.** Renamed or retired scripts must be
  removed from the cluster by hand; this has caused three separate stale-artifact bugs. Local
  absence does not imply cluster absence — a stale `exclude_af_concordance.txt` existed only on the
  cluster and was silently applied.
- **A joint call emits nothing at sites monomorphic in its own donors.** BR-DSNWGS therefore
  contributes to only about half the association set. Expected, not a defect — but it is a methods
  limitation, and BR sits entirely on the AMP-PD side of the primary contrast.
- **venv Python must use `module load python/3.11`** — the default system Python (`/usr/local/bin/python`) points to an Anaconda env on compute nodes that breaks `pkg_resources`. Always load `python/3.11` before activating `.venv` in sbatch scripts.
- **setuptools must be pinned to 67.8.0** — newer versions break `pkg_resources` import needed by genotools. Fix: `pip install setuptools==67.8.0`.
- **`subprocess.run()` with `capture_output=True` hides plink2 output** — use without `capture_output` or add `stderr=subprocess.STDOUT` so output appears in sbatch logs.
- **Always use `--dependency=afterok` when chaining sbatch jobs** — avoids downstream jobs running on incomplete input.
- **Synapse downloads may be interrupted** — always rerun `synapse get -r` to resume; it skips completed files and resumes incomplete ones.
- **tabix `.tbi` files may be stale after download** — use `tabix -f -p vcf` to force reindex before running bcftools with `--regions`.