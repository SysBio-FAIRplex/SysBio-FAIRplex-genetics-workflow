# ad-pd-gwas

AD-vs-PD GWAS across AMP-AD and AMP-PD, restricted to donors with WGS, with **AD defined by
neuropathology**.

## The goal

1. One harmonized clinical dataset from AMP-AD + AMP-PD. Primary contrast AD vs PD; secondary
   AD and PD vs other dx.
2. WGS QC pipeline on biowulf.
3. GWAS across all contrasts + plots.
4. Follow-up PCA to confirm the callsets align.

## One root

Code and data share a root. On biowulf that root is `/data/CARDPB2/sysbio/wgs`; on a laptop it is
wherever the repo was cloned. **Nothing in this project contains an absolute path** — `config.sh`
derives the root from its own location, and `clinical_core.py`, `05_excludelist.py` and
`ancestry_qc_manifest.py` each do the same. To move the project, copy the folder.

That collapses three things that used to be separate and disagree: `WGS_ROOT`, `PROJECT_ROOT`, and
the old `BUNDLE` (which pointed one directory level away from where the code actually was).

**Every handoff file has exactly one location.** The clinical side writes to `clinical_core_out/`
and the genotype steps read it from there — nothing is copied into `data/`, so nothing can go
stale against the run that produced it.

| | |
|---|---|
| `clinical_core.py` | goal 1. 5 cohorts, 4 callsets. Reads `data/`, writes `clinical_core_out/`. |
| `wgs_core.ipynb` | goal 2 orchestrator — ssh/sbatch driver. Complete through step 2. |
| `config.sh` | every path, derived from its own location. `submit.sh` — sbatch wrapper. |
| `scripts/00`–`07` | env, genotools, normalize, merge, relatedness, excludelist, ancestry QC, GWAS. |
| `scripts/ancestry_qc_manifest.py` | step 6 → clinical boundary: `retained_samples_manifest.csv`. |
| `scripts/gene_annot.py` + `ref/refFlat.txt` | locus coordinates. CLI: `python3 scripts/gene_annot.py CR1 SNCA LRRK2` |
| `scripts/diag_order.py`, `scripts/diag_cah.sh` | read-only. Variant order vs the panel; postmortem of a genotools output dir. |
| `review/` | goals 3 and 4: QQ/Manhattan, ctrl-vs-ctrl mask, eta² of callset on each PC. |
| `data/**/metadata/` | the 11 clinical files `clinical_core.py` opens. Gitignored — controlled access. |
| `PROJECT_LOG.md` | append-only: what we did, why, and what was ruled out. This file is state; that one is history. |

Genotype psams are read **in place** from each callset's `pgen/` or `joint_calls/`. There is no
copy under `metadata/`.

## Run order

Everything runs on the cluster, so there are no round trips. `clinical_core.py` runs twice —
§10 must precede genotools, §11–13 must follow it — and the waiting sections skip cleanly and
say so.

```
git pull                     # cluster takes the code
python3 clinical_core.py     # §1-10 -> clinical_core_out/*_update_sex.txt  (§11-13 print SKIPPED)
./submit.sh scripts/01_genotools.sh ...   # steps 1-6; step 1 reads the sex files in place
python3 clinical_core.py     # §11 QC outcomes, §12 -> analysis_grain.csv now fill in
./submit.sh scripts/07_gwas.sh            # reads $GRAIN from clinical_core_out/
review/plot_gwas.py          # + mask_cohort_artifacts.py
```

## Status

**All four callsets are merged.** `cohort_merged` is now 172,497,055 variants × 13,334 samples
(job 27429821, 2026-08-15). Everything downstream of the merge — relatedness, excludelist,
per-ancestry QC sets and PCs — was built from the old 3-callset cohort and is stale.

| Step | BR-DSNWGS |
|---|---|
| 0 VCF→pgen | done — `br_dsnwgs_hg38.{pgen,psam}` — **but built by notebook cells, not a script** |
| 1 genotools filter | done — `br_dsnwgs_hg38_filtered*` |
| 1 genotools ancestry | done 2026-08-14 — 77 EUR / 14 AJ / 3 AMR / 1 CAH / 1 AFR / 1 AAC |
| 2 normalize | done 2026-08-14 — job 27429119 |
| 3 merge | done 2026-08-15 — 97/97 samples in, variant arithmetic closes exactly |
| 4 relatedness | **next** — blocked on `LBL_BR` (see below) |

**`04_relatedness.sh` does not know BR-DSNWGS exists.** Lines 67 and 100 list `LBL_WGS`,
`LBL_WB` and `LBL_DC` only. Unfixed, the 97 BR samples enter the merged cohort, get no
ancestry label, drop out of every stratum, and are never tested for relatedness — with no
error raised. Fix before running step 4.

`FILTERED.br_dsnwgs_ancestry_umap_linearsvc_predicted_labels.txt` is `LBL_BR` in `config.sh` and
is read by steps 4 and 6. It now holds a real ancestry spread and is usable. (Note the filename
says `linearsvc` but the model is XGBoost — a GenoTools naming artifact, not a description.)

**`analysis_grain.csv`'s BR-DSNWGS rows are wrong and must be regenerated.** The grain carries
19 AFR / 67 EUR over 86 rows; the correct run gives 1 AFR / 77 EUR over 97 samples. Not
reconcilable — regenerate the PCs and labels before step 7.

**Blocker 1 — genotools crashed. FIXED 2026-08-12.** GenoTools sized its GridSearchCV worker
pool from the *node* (`os.cpu_count()`), not the SLURM allocation, and biowulf's per-user
`ulimit -u` of 1024 could not hold that many worker processes — `pthread_create` returned
EAGAIN, workers aborted on SIGABRT, and SLURM logged a bare `ExitCode 1:0`. Fixed by
`scripts/genotools_capped.py`, which reports the allocation instead, **plus
`GENOTOOLS_MAX_WORKERS=16`** — the allocation alone is not enough, 64 workers still aborts 13 of
them. Never specific to BR-DSNWGS: only 3 of 18 genotools jobs had ever completed. **Note:** an
earlier version of this paragraph credited pinning `OMP_/NUMBA_NUM_THREADS` to 1 — that was
wrong, joblib already does it, and the worker *count* was the lever.

**Blocker 2 — all 97 samples predicted CAH. FIXED 2026-08-14.** Not missing data and not a model
problem. BR-DSNWGS's pgen was ordered `1,10,11,…,19,2,20,21,22,3,…` — per-chromosome files
concatenated alphabetically — while the reference panel is numeric, and GenoTools aligns the
study matrix to the panel **by column position, not by name** (`ancestry.py:247`; the reorder at
`:244` that would prevent it sits inside `if not self.train` and never runs when training). Every
shared column was standardized and projected as a *different* variant, destroying the ancestry
signal; the resulting noise clump sits near the centre of the map, which is exactly the CAH rule.
Fixed with `--sort-vars` in `01_genotools.sh` step 1.

**The other three callsets are unaffected — measured, not assumed.** `scripts/diag_order.py`
reports 0 rank drops for `wgs_harm`, `divco_hs` and `wb_dwgs`; they worked only because their
files happened to be numeric. Their labels, the merge, and the grain are sound.

Run `python3 scripts/diag_order.py <panel>.bim <callset>.pvar` on any new callset before trusting
its ancestry output. Full chain and every hypothesis ruled out are in `PROJECT_LOG.md` — read
"Where we are right now" at the top of that file first.

The clinical side is already ahead of it: §7 resolves all 97 BR samples via the AMP-PD sample
inventory, and §10 writes `br_dsnwgs_update_sex.txt` (60M / 37F).

## Known issues

1. **Two implementations of pheno/covar.** `clinical_core.py` §13 writes
   `clinical_core_out/{pheno,covar}/`, but `07_gwas.sh` does not read them — it rebuilds both from
   `$GRAIN` in awk. The awk version is what actually runs. Pick one.

2. **Two implementations of the ctrl-vs-ctrl artifact filter.** `review/mask_cohort_artifacts.py`
   is here; `08_ctrl_ctrl_filter.{py,sh}` is in `scripts_archive/` on the cluster. The archived one
   carries the argument that APOE is *expected* to reach significance in the control-vs-control
   scan (AMP-PD controls are screened for PD, not AD), so subtracting destructively would delete
   the study's strongest true locus. Read it before trusting the `review/` one.

3. **`06_ancestry_qc.sh` reads `$MERGED_DIR/exclude_af_concordance.txt`**, which nothing here
   builds — `06a_af_concordance.{py,sh}` is in `scripts_archive/` on the cluster. The step prints a
   `NONE` note and continues *without the exclusion*. That filter existed because callset
   separation dominated PC1 in AJ (eta² 0.98). Promote those two files if goal 4 resurfaces.

4. **`all_chrs_merged.psam` is not in the local `data/`**, so `wb_dwgs` skips when
   `clinical_core.py` runs on a laptop (3,010 genomes instead of 13,428). Everything else is
   unaffected. To run the full thing locally:
   `rsync biowulf.nih.gov:$ROOT/data/amp-pd-genomics/WB-DWGS/joint_calls/all_chrs_merged.psam data/amp-pd-genomics/WB-DWGS/joint_calls/`

5. **The cluster's `README.md` (30 KB) has not been merged in.** `06_ancestry_qc.sh` references
   its §2 (reference data acquisition) and §3 (per-step commands) by number. Those sections
   should be folded into this file.

6. **DivCo's source VCF is 0 bytes** on the cluster (`merged.deduped.vcf.gz`). The pgen was derived
   before it was truncated, so nothing is blocked, but DivCo cannot be re-derived from source
   without re-pulling from Synapse.

## Provenance

`scripts_archive/ad-pd-gwas-jul28/` on the cluster holds the July version of every script, plus 15
files that were never carried over: the AF-concordance pair, the ctrl-vs-ctrl step, and the
`diag_*` / `probe_*` diagnostics. Of the 12 scripts shared with this repo, 8 are byte-identical;
`config.sh`, `06_ancestry_qc.sh`, `07_gwas.sh` and `05_excludelist.py` differ, with this copy
taken as authoritative.

Not committed anywhere yet — `git init` has been run, but there is no remote.
