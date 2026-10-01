# WGS Harmonization Pipeline

Processes whole-genome sequencing from multiple cohorts into one harmonized GRCh38 set for ancestry
prediction, sample and variant QC, and GWAS. The analysis it feeds is **AD vs PD across AMP-AD and
AMP-PD, with AD defined by neuropathology**.

| doc | answers |
|---|---|
| `README.md` (this file) | what the pieces are, where files live, in what order things run, and the known gaps a rerun inherits |
| `METHODS.md` | **why** each choice was made — draft methods section for the paper |
| `wgs_core.ipynb` | the **executable** run order — it runs on biowulf and submits every step |
| `SUMSTATS_README.md` | the published summary statistics: layout, columns, provenance |
| `CLAUDE.md` | the rules a Claude Code session in this repo follows |

`RUNLOG.md` (which jobs ran and how they ended) is generated on the cluster, not versioned:
`bash scripts/runlog.sh --md > RUNLOG.md`.

---

## Where things run

**Everything runs on NIH biowulf** — the genotype steps, both clinical scripts, the review figures.
The data is controlled-access (AMP-AD via Synapse, AMP-PD via GCP) and both data use agreements
confine it to the cluster.

- A clone off the cluster is for **editing code only**. It holds no data and runs nothing.
- **Nothing sample-level leaves the cluster**, in either direction: `data/`, `clinical_core_out/`,
  `results/` (except the tracked aggregate `results/pca/*.csv`), sumstats, manifests, psams,
  eigenvecs, step logs. Only aggregate text and figures come off it.
- **No participant IDs in code, comments, docs or commit messages.** `.gitignore` covers data
  paths; the git hooks refuse stored notebook outputs (`nb_guard.py`) and ID-shaped text in files
  and commit messages (`id_guard.py`). Purely numeric IDs have no shape to match, so those are on you.
- Bulk downloads go through the **helix** transfer node; compute nodes have no outbound route to
  Synapse or GCS. `wgs_core.ipynb` §0 holds every acquisition command.

**One root.** Code and data share a root on biowulf: `/data/CARDPB2/sysbio/wgs/amp-ad-pd-wgs-gwas`.
Nothing in the project holds an absolute path — `config.sh`, `clinical_common.py` and
`wgs_core.ipynb` derive the root from their own location. To move it, copy the folder.

**Code reaches the cluster through GitHub.** It is pushed from the development clone and pulled on
biowulf. The repo is public, so the cluster checkout holds **no GitHub credentials and stays
read-only**: nothing on the cluster can push, so no data can leave it by git. Aggregate outputs
that belong in the repo (the eta² tables) are copied off the cluster and committed from the
development clone.

```bash
# biowulf login node (or helix), in the project root
git pull --ff-only
```

Run git on the login node or helix, not in a batch job; compute nodes reach the internet only
through a proxy. Code moves only through git: any other copy cannot express a deletion, and a
retired file left behind on the cluster is silently run. `git rev-parse HEAD` on the cluster is the commit that ran — step 9
stamps it into every release.

---

## Infrastructure

- **Cluster**: NIH Biowulf (SLURM). **Transfer node**: NIH Helix.
- **Tools**: plink2 v2.00a6LM, plink1.9, bcftools, liftOver, GenoTools v1.3.6, Synapse CLI
- **Python**: `.venv` in the project root, Python 3.11. Always
  `module load python/3.11 && source .venv/bin/activate` — the system Python on compute nodes is an
  Anaconda py3.9 that breaks `pkg_resources`.
- **`requirements.lock.txt`**: 115 pinned packages, frozen from the working cluster venv. It does
  not contain `setuptools` (`pip freeze` omits it), and setuptools must be **67.8.0** for
  GenoTools' `pkg_resources` import. Rebuild in this order:
  `python3 -m venv .venv && pip install setuptools==67.8.0 && pip install -r requirements.lock.txt`

---

## Directory Structure

```
/data/CARDPB2/sysbio/wgs/amp-ad-pd-wgs-gwas/
  bin/
    liftOver                                         ← UCSC liftOver binary
    hg19ToHg38.over.chain.gz                         ← liftover chain file
  data/
    ref/                                             ← see Reference data
      GRCh38_full_analysis_set_plus_decoy_hla.fa.zst
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
        gwas/                                       ← STEP 7 OUTPUT: sumstats, gwas_summary.csv, plots/
    amp-ad-genomics/
      WGS_Harmonization/
        joint_calls/                                 ← 24 b37 GATK scatter-interval VCFs + .tbi
        metadata/
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
      metadata/                                      ← AMP-PD participants, case/control, demographics
      WB-DWGS/
        joint_calls/
          all_chrs_merged.{pgen,pvar,psam}           ← original GRCh38 pgen (unfiltered)
          all_chrs_merged_filtered_bed.{bed,bim,fam} ← filtered bed (genotools input/byproduct)
        metadata/
        genotools/
          FILTERED.wb_dwgs.json                      ← ancestry + QC results
          FILTERED.wb_dwgs_ancestry_*.{pgen,psam,pvar,samples}
          FILTERED.wb_dwgs_{ANCESTRY}.{pgen,psam,pvar}
      BR-DSNWGS/
        joint_calls/
          AMPPD_postmortem_joint_gt_call_97donors.vcf.gz   ← 97 postmortem donors
        metadata/
        pgen/
          br_dsnwgs_hg38.{pgen,pvar,psam}            ← GRCh38 pgen (see --sort-vars)
        genotools/
  clinical_core_out/                                 ← THE PHENOTYPE BOUNDARY. One location each;
    {callset}_update_sex.txt                         ←   nothing is copied into data/, so a handoff
    sample_annot.csv                                 ←   file cannot go stale against its producer.
    individual_core.csv                              ← §9 audit table -> analysis_grain.py
    genome_crosswalk.csv                             ← §9 audit table -> analysis_grain.py
    analysis_grain.csv                               ← the grain. Step 7 does NOT read it directly
    pheno/  covar/  contrasts.csv                    ←   — it reads THESE (§13)
    qc_outcomes.csv
  results/                                           ← figures + eta^2 tables (only the eta^2 CSVs are versioned)
  scripts/                                           ← see Scripts below
    logs/                                            ← all sbatch .o/.e files
  ref/                                               ← ships with the code: highld BED, refFlat, panel checksums
  README.md   METHODS.md   SUMSTATS_README.md   CLAUDE.md
  config.sh   submit.sh
  clinical_common.py   clinical_core.py   analysis_grain.py
  wgs_core.ipynb                                     ← the orchestrator; runs ON biowulf
```

---

## Datasets

| Dataset | Cohort | Samples | Build | Source format |
|---|---|---|---|---|
| WGS_Harmonization | AMP-AD (MayoRNAseq, MSBB, ROSMAP) | 1,894 | b37 → GRCh38 | 24 GATK scatter VCFs |
| WB-DWGS | AMP-PD | 10,418 | GRCh38 | Merged pgen |
| DivCo_HS | AMP-AD diverse cohorts | 1,019 | GRCh38 | Single merged VCF |
| BR-DSNWGS | AMP-PD (postmortem) | 97 | GRCh38 | Single joint-call VCF |

The current run's merged cohort is 172,497,055 variants × 13,334 samples, 12,495 retained after
step 5. Where every input came from: `wgs_core.ipynb` §0.

---

## Scripts

Every step is a numbered script in `scripts/`, submitted through `submit.sh`, which resolves the
root and sources `config.sh`.

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
  hooks-pre-commit             ← install both hooks into .git/hooks/ after a clone
  hooks-commit-msg

  # ── tier 3, on demand ─────────────────────────────────────────────────────
  gene_annot.py                ← what gene is this variant in; locus coordinates from ref/refFlat.txt
  runlog.sh                    ← regenerates RUNLOG.md from the SLURM accounting log
  nb_guard.py                  ← refuses to commit a notebook carrying stored outputs
  id_guard.py                  ← refuses to commit text or a message shaped like a participant ID

  # ── called by step 1, not run directly ────────────────────────────────────
  genotools_capped.py          ← GenoTools entry point that respects the SLURM allocation

review/                        ← figures and checks, run on biowulf against outputs in place.
                                 Only the eta^2 tables and QQ/Manhattan PNGs leave the cluster;
                                 the PC scatters plot one point per participant and stay.
  plot_af_filter_effect.py     ← the AF-filter before/after proof figure + eta^2 tables
  plot_pcs_by_callset.py       ← PC scatter by callset + eta^2 per PC per stratum
  plot_gwas.py                 ← QQ + Manhattan, with lambda_GC cross-checked against step 7's
  drift_check.py               ← does the working tree differ from a revision in anything that
                                 RUNS? AST for .py, quote-aware comment strip for .sh. Read-only.
                                 `python3 review/drift_check.py [<rev>]`
  methods_numbers.py           ← re-derives every number in METHODS.md from the artifacts and
                                 diffs it against the doc. Read-only. Run it before any write-up
                                 edit; --strict to fail on anything it could not check.
                                 --discordance measures §6.4's duplicate-pair enrichment.
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

Grouped by **when it runs and what triggers it**, not by what it is allowed to do: `diag_order.py`
is a mandatory gate, and filing it with optional diagnostics is how all 97 BR-DSNWGS samples once
came back mislabelled with a clean exit and 0.97 nominal model accuracy.

`wgs_core.ipynb` is the executable copy of tier 1. This section is the map. All of it runs on
biowulf.

### Tier 1 — the pipeline. Runs every time, in this order.

```
clinical_core.py                §1-10 sex files, §12a sample_annot.csv
  01 genotools   (x4)           sex update -> filter -> bed -> ancestry + QC
  02 normalize   (x4)           chr naming + variant IDs to a common convention
  03 merge                      plink1.9 union of all four -> cohort_merged
  04 relatedness                KING, report-only
  05 excludelist                duplicates, relatives, QC fails -> retained_manifest.csv
  06 ancestry QC                ONE submission, five stages:
                                  A unfiltered QC+PCA   B build the AF exclusion list from A
                                  C apply it            D filtered QC+PCA
                                  E both sample manifests
analysis_grain.py               §11-13 on step 6's PCs
  07 gwas                       plink2 --glm per ancestry per contrast
  review/plot_gwas.py           QQ + Manhattan
  08 ctrl-vs-ctrl               annotates; never subtracts
```

- **`analysis_grain.py` must run before any step-7 run.** A `contrasts.csv` without the
  callset-skew columns makes step 7 exit 3 naming the missing one. Rerun the grain; do not patch
  the CSV.
- **Read step 6 stage B's log before step 7.** Stage C applies the exclusion list later in the same
  job, so nothing in stage B can gate. Check the BY CALLSET PAIR table and the per-cell HWE ratio
  table, and resolve anything odd against the per-cell `.afreq`/`.snplist` in `af_concordance/`
  (by header name). Call rate is usually the tell, not frequency.
- `review/plot_af_filter_effect.py` runs after step 6, off the two manifests that step wrote.

### Tier 1b — release. Runs after step 8, when data is being handed to someone else.

```
  09 amppd subset               sbatch on BIOWULF  — subset + manifest. Never touches the network
  10 amppd push                 run on HELIX       — upload + verify
```

**Do not run either without written confirmation that the recipient is covered by the AMP-PD data
use agreement.** This ships individual-level genotypes. The script can prove a bucket is private;
it cannot prove that the people with access to it are approved users.

It consumes the pipeline's output and produces nothing the pipeline reads back, so it is not
numbered into tier 1 — but it is a step with a fixed position (after step 8, on artifacts step 6
wrote), not a tool for an arbitrary moment.

| | |
|---|---|
| what it ships | `data/merged/by_ancestry_qc/cohort_<ANC>_qc` — step 6 **stage C**, the exact variant set step 7 tested — restricted to `source_callset` in {`wb_dwgs`, `br_dsnwgs`} |
| shape | one fileset **per ancestry stratum**, because the QC is per stratum. Concatenating them would invent a variant set no GWAS in this study ran on |
| who is AMP-PD | read by header name from step 6 stage E's `retained_samples_manifest.csv`. **Not** re-derived from the four genotools label files — `ancestry_qc_manifest.py` is the one implementation of callset membership |
| not re-filtered | dropping the AMP-AD donors leaves some variants monomorphic within the subset; they stay. Re-applying `--maf` would produce a fileset that disagrees with the published sumstats. `SUBSET_MAF`/`SUBSET_GENO` opt in, and the release README says which was used |
| destination | **no bucket is hardcoded.** `GCS_DEST` is required at push time |

**Two scripts on two hosts, and the split is not cosmetic.** Biowulf compute nodes have no general
outbound network, so a push from one fails looking like an auth problem; step 10 refuses to run
inside a SLURM allocation. Keeping them separate also means step 9 has no gcloud dependency — it
records size and SHA-256 with coreutils, which a consumer can check after download. CRC32C is what
GCS stores, so step 10 computes it at upload time and compares against the bucket. Not MD5
anywhere: GCS computes no MD5 for composite (parallel-chunked) uploads.

```bash
# on biowulf
./submit.sh scripts/09_amppd_subset.sh

# then on helix.nih.gov, after step 9 reports "reconciled"
GCS_DEST=gs://<bucket>/<prefix> bash scripts/10_amppd_push.sh
```

**Three things it refuses to do, none of them overridable by a flag:**

- Ship a **partial** release. Every AMP-PD sample in the retained manifest must land in exactly
  one released fileset; the reconciliation is fatal, and it names the strata that came up short.
- Upload to a bucket whose metadata **reads fine and says it is public** — uniform bucket-level
  access off, public access prevention not enforced, or an `allUsers` binding. Those reads need
  `storage.buckets.get`/`.getIamPolicy`, which `roles/storage.objectAdmin` does not grant, so on
  a program-managed bucket the check often cannot run at all; that case is reported loudly and
  the upload proceeds on the operator's say-so. `SKIP_BUCKET_CHECK=1` skips it outright. A bucket
  that can be neither described nor listed is a wrong name and still exits.
- Report success on an **unverified** upload. Stage E re-lists the bucket and compares presence,
  size and CRC32C against `MANIFEST.tsv`.

Neither script has been run against real data. Unverified until the first real build: the gcloud
module name on helix (`MOD_GCLOUD` defaults to `google-cloud-sdk`, falling back to `gcloud` on
`PATH`), and whether every AMP-PD sample reconciles — donors in a `PRUNE_FAIL` stratum would trip
it.

### Tier 2 — preflight. Runs when a precondition changed, not every time.

| check | trigger | if it fails |
|---|---|---|
| `bash scripts/00_setup_env.sh` | a fresh cluster checkout with no `.venv` | nothing else will run |
| `cp scripts/hooks-pre-commit .git/hooks/pre-commit && cp scripts/hooks-commit-msg .git/hooks/commit-msg && chmod +x .git/hooks/pre-commit .git/hooks/commit-msg` | any fresh clone | notebook outputs and participant IDs can reach git; both have happened |
| `python3 scripts/diag_order.py <panel>.bim <callset>.pvar` | **any callset not already cleared, BEFORE step 1** | non-zero rank drops -> add `--sort-vars` at import, or ancestry output is garbage that exits 0 |

`wgs_harm`, `divco_hs`, `wb_dwgs` and `br_dsnwgs` have all been cleared. A new callset has not.

### Tier 3 — on demand, on biowulf.

| tool | use |
|---|---|
| `python3 scripts/gene_annot.py --at chr19:44908684` | what gene does this variant sit in — this is what named APOE and HLA-DQB1 in step 8's results |
| `python3 scripts/gene_annot.py CR1 SNCA LRRK2` | gene coordinates, optionally `--flank` |
| `bash scripts/runlog.sh --md > RUNLOG.md` | regenerate the job table from SLURM accounting |
| `python3 scripts/nb_guard.py --staged` | what the pre-commit hook runs; `--strip <nb>` to clean one |
| `python3 scripts/id_guard.py --all` | audit every tracked file for ID-shaped text |

---

## Standing constraints

A run breaks, or silently goes wrong, if any of these changes.

- **`COMMON_GENO=0.005` in step 4 is load-bearing.** It must stay below the smallest callset's
  share of the cohort (BR is 97/13,334 = 0.0073). At `--geno 0.05` every variant BR lacked stayed
  in the common set and all 97 BR samples read as 50% missing, which step 5 then consumed as a
  duplicate tie-break. Re-check if a callset under ~0.5% is ever added.
- **`plink1.9` is a hard dependency of step 6 stage B.** The frequency test is `plink --assoc`
  (the 1-df allelic chi-square); plink2 dropped `--assoc` for `--glm`, which is only asymptotically
  equivalent. If the module goes away the fallback is `--glm` plus a fresh equivalence check.
- **GenoTools aligns to the reference panel by COLUMN POSITION, not by variant name.** A pgen
  ordered `1,10,11,…,19,2,20,…` against a numeric panel projects every shared column as a
  different variant. Hence `--sort-vars` in step 1 and `diag_order.py` before it on any new callset.
- **GenoTools sizes its worker pool from the NODE, not the SLURM allocation**, and exceeds
  biowulf's per-user `ulimit -u` of 1024 (a bare `ExitCode 1:0`). Needs both
  `scripts/genotools_capped.py` and `GENOTOOLS_MAX_WORKERS=16`; 64 still fails.
- **Only step 6a may delete variants; step 8 annotates.** Different entitlements, not different
  degrees of caution — `METHODS.md` §6.2 and §8.
- **A joint call emits nothing at sites monomorphic in its own donors**, so BR-DSNWGS contributes
  to only about half the association set. Expected, and a methods limitation: BR sits entirely on
  the AMP-PD side of the primary contrast.

---

## Reference data

Not in the repo. `wgs_core.ipynb` §0 fetches the public files on helix.

| file | lives in | source | used by |
|---|---|---|---|
| `GRCh38_full_analysis_set_plus_decoy_hla.fa.zst` | `data/ref/` | 1000 Genomes GRCh38 reference, `ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/`, zstd-compressed (plink2 `--fa` reads it directly) | step 2 |
| `ref_panel_gp2_prune_rm_underperform_pos_update.{bed,bim,fam}`, `ref_panel_ancestry_updated.txt` | `data/ref/` | GenoTools ancestry reference panel — **archived upstream**, see below | step 1 |
| `hg19ToHg38.over.chain.gz`, `liftOver` | `bin/` | UCSC, `hgdownload.soe.ucsc.edu` | GRCh37 callsets only |
| `highld_exclude_hg38.bed`, `refFlat.txt` | `ref/` | versioned with the code | step 6, `gene_annot.py` |

**The ancestry panel cannot be re-downloaded.** It is a GenoTools reference panel that was archived
upstream after this study ran; `genotools-download` now serves `1kg_30x_hgdp_ashk_ref_panel`
instead. The copy in `data/ref/` is the copy of record: keep a backup, and verify it with
`md5sum -c ref/ref_panel.md5` from the project root. Step 1 trains GenoTools' ancestry model from
this panel on every run, so:

- to reproduce this study, or add a callset comparable with it, use the archived panel;
- a new study can use the current panel, but its ancestry labels will not match this study's, and
  `diag_order.py` must clear every callset against the new panel before step 1.

---

## Known gaps

Open items a rerun or the write-up inherits.

- **GRCh37 callsets have no script.** `wgs_harm`'s route (per-chromosome filter, merge, `liftOver`
  to hg38) is described in `wgs_core.ipynb` §1 but was never scripted; a new b37 callset needs one.
- **Two of the five DivCo metadata files have no recorded provenance**; `wgs_core.ipynb` §0 has
  syn IDs for three.
- **DivCo's source VCF (`merged.deduped.vcf.gz`) is 0 bytes on the cluster.** The pgen was built
  before it was truncated; re-deriving DivCo needs a re-pull from Synapse.
- **Step 0 (VCF→pgen for BR-DSNWGS) has no script.** `wgs_core.ipynb` §1 is the whole record.
- **No age covariate.** AMP-AD gives age at death, AMP-PD age at baseline; one column would mix two
  quantities. Needs a harmonized age from the phenotype side. `METHODS.md` §10.
- **The callset-skew columns are not in the current results.** `analysis_grain.py` §13 emits
  `max_callset_delta` / `worst_callset` / `callset_one_sided`, but step-7 job 28004190 predates
  them. Rerun the grain, then step 7.
- **`METHODS.md` §6.4's "~7× enriched for duplicate-pair discordance" is unsourced** — it is the
  only non-circular evidence that the flagged variants are technical. Measure it with
  `python3 review/methods_numbers.py --discordance`, or delete the sentence.
- **The step-6 job ID behind the live 4,187-variant list is not recorded.** Recover it with
  `bash scripts/runlog.sh --md`.
- **Fewer than 20 EUR `wgs_harm` controls differ between step 6a's HWE cell and §13.** The `.keep`
  holds 328; 6a ran before §13 became the sole definition of a case. Check whether any of them are
  now cases.
- **Donors sequenced in two callsets may carry different ancestry labels** (always AFR vs EUR).
  Reported by a side analysis against a superseded grain; unverified on the current one. Check
  which label each retained genome carries.
- **Sex updates for the current run, unverified:** `wgs_harm`, `divco_hs` and `wb_dwgs` may have
  been sex-updated from `<callset>/metadata/` rather than `clinical_core_out/`. A rerun reads
  `clinical_core_out/` (`config.sh: sex_file`).
- **`gwas_summary.csv` is not versioned.** It lives on the cluster
  (`data/merged/by_ancestry_qc/gwas/`) and in the sumstats release bucket, and every §7 number traces
  to it. It is aggregate (one row per contrast); decide whether to version it.
  Superseded copies sit under `by_ancestry_qc_pre_pcfix_20260726/` and `plotdata/`.
- **GenoTools fixes not upstreamed** (`dvitale199/GenoTools`): positional column alignment
  (`ancestry.py:247`), worker sizing, and having `get_common_snps` report its overlap.

---

## Sex crosswalk

All four sex files are written by `clinical_core.py` §10 to
**`clinical_core_out/{callset}_update_sex.txt`**, and step 1 reads them there
(`config.sh: sex_file`). Nothing is copied into a callset's `metadata/`.

- **DivCo_HS** — psam IDs come in three formats, resolved separately in `clinical_core.py` §7:
  direct `individualID` match (620), suffix stripping of `<individualID>_DLPFC_WGS` (293), and
  biospecimen lookup of `<specimenID>-D` (106). 1,019/1,019 matched.
- **WGS_Harmonization** — over 99% of 1,894; the few remaining are of unknown sex.
- **WB-DWGS** — 10,418/10,418.
- **BR-DSNWGS** — 97/97, via the AMP-PD sample inventory.

---

## Key gotchas

- **GATK scatter VCFs are not chromosome-split** — each WGS_Harmonization VCF holds variants from
  several chromosomes. Never assume file number = chromosome number.
- **`#` in swarm files is always a comment** — swarm strips it even inside quotes. Never use
  `--set-missing-var-ids @:#` in a swarm file; use `bcftools annotate --set-id`. In sbatch, `#` is fine.
- **plink2 v2.00a6LM's non-concatenating `--pmerge-list` is not implemented** — multi-cohort merges
  use plink1.9.
- **plink1.9 does not recognize PAR1/PAR2** — use `--merge-par` when converting pgen → bed if
  `--split-par` was used at VCF import.
- **chrX needs `--split-par` + `--update-sex`** — `hg19` for b37, `hg38` for GRCh38. The sex file
  must be 2-column (IID SEX) if the psam has no FID column.
- **liftOver needs `chr`-prefixed chromosomes** — b37 uses bare numbers.
- **`set -o pipefail` must come after all `#SBATCH` lines** — SLURM stops parsing directives at the
  first non-comment line, silently dropping every resource request after it.
- **DivCo_HS FILTER is all `.`** — no VQSR; `--apply-filters PASS` drops every variant.
- **`subprocess.run(capture_output=True)` hides plink2's errors** — use `stderr=subprocess.STDOUT`
  so they reach the sbatch log.
- **Chain sbatch jobs with `--dependency=afterok`**, so nothing runs on incomplete input.
- **Synapse downloads get interrupted** — rerun `synapse get -r`; it resumes.
- **tabix `.tbi` files can be stale after download** — `tabix -f -p vcf` before `bcftools --regions`.
- **A file count is not a file list.** Compare names when checking that two directories match.
