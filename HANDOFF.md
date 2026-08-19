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

Everything runs on the cluster, so there are no round trips. **Code reaches the cluster by
`rsync`, not `git pull`** — there is no remote (see Provenance). Because `rsync` without
`--delete` cannot express a deletion, renamed or retired scripts must be removed by hand; that
has caused three separate stale-artifact bugs, so check before assuming.

`clinical_core.py` runs **three** times and step 6 runs **twice**:

```
rsync -av --exclude='logs/' --exclude='__pycache__/' scripts/ helix:$ROOT/scripts/

module load python/3.11 && source .venv/bin/activate   # NOT the system Anaconda py3.9
python3 clinical_core.py                  # §1-10 -> *_update_sex.txt   (§11-13 SKIP)

./submit.sh scripts/01_genotools.sh ...   # once per callset; reads the sex files in place
./submit.sh scripts/02_normalize.sh ...   # once per callset
./submit.sh scripts/03_merge.sh           # all four -> cohort_merged
./submit.sh scripts/04_relatedness.sh     # KING, report-only
./submit.sh scripts/05_excludelist.sh     # -> retained_manifest.csv
./submit.sh scripts/06_ancestry_qc.sh     # PASS 1, unfiltered
python3 scripts/ancestry_qc_manifest.py   # -> retained_samples_manifest.csv (NOT called by 06)
python3 clinical_core.py                  # §11-13 -> analysis_grain.csv

python3 review/plot_pcs_by_callset.py --manifest <manifest>   # eta^2: is the AF filter needed?
./submit.sh scripts/af_concordance_build.sh # -> exclude_af_concordance.txt (needs $GRAIN)
./submit.sh scripts/06_ancestry_qc.sh     # PASS 2, picks the list up automatically
python3 scripts/ancestry_qc_manifest.py   # PCs changed
python3 clinical_core.py                  # grain must be rebuilt on the new PCs

./submit.sh scripts/07_gwas.sh            # reads $GRAIN from clinical_core_out/
python3 review/plot_gwas.py               # + mask_cohort_artifacts.py
```

**The two-pass shape is forced, not stylistic:** `af_concordance_build` needs step 6's
`by_ancestry_qc` output *and* the grain, and step 6 cannot depend on the grain without a circular
ordering (grain ← §12 ← manifest ← step 6). Step 6 refuses an exclusion list older than
`cohort_merged.bed`, so a stale one now fails fast instead of being applied silently.

## Status

**All four callsets are through step 6 pass 1, and the grain is rebuilt.** `cohort_merged` is
172,497,055 variants × 13,334 samples (job 27429821); `analysis_grain.csv` is 12,495 rows × 22 cols,
built 2026-08-18 from the current manifest.

| Step | state |
|---|---|
| 0 VCF→pgen | done — `br_dsnwgs_hg38.{pgen,psam}` — **but built by notebook cells, not a script** |
| 1 genotools | done all four — BR: 77 EUR / 14 AJ / 3 AMR / 1 CAH / 1 AFR / 1 AAC |
| 2 normalize | done all four — BR job 27429119 |
| 3 merge | done — 97/97 BR in, variant arithmetic closes exactly |
| 4 relatedness | done — 13,334/13,334 labelled, 11 strata |
| 5 excludelist | done — **12,495 retained** (839 excluded), 95/97 BR |
| 6 ancestry QC **pass 1** | done unfiltered, job 27602590 — 6 of 11 strata have PCs |
| `ancestry_qc_manifest.py` | done — 12,495-row `retained_samples_manifest.csv` |
| `clinical_core.py` §11–13 | done 2026-08-18 — grain rebuilt, 17 of 44 contrasts viable |
| `af_concordance_build.sh` | done 2026-08-18, job 27697096 — **4,415 variants**, correct grain |
| premise test (eta² before/after) | done 2026-08-19 — **EUR 0.757 → 0.036; AJ 0.984 → 0.962** |
| 6 ancestry QC **pass 2** | next — `mv by_ancestry_qc` FIRST (one-way door), then submit |
| `clinical_core.py` again | then — PCs change, so the grain must be rebuilt |
| 7 GWAS | after that |

**The filter is warranted, and its effect is now MEASURED rather than cited.** Unfiltered eta² on
the four-callset cohort: AJ PC1 **0.984**, EUR PC2 **0.757**. Removing BR changes neither
(0.984 / 0.758), so this is real `wgs_harm`↔`wb_dwgs` structure. Applying the 4,415-variant list:

| stratum | max eta² unfiltered | max eta² filtered | reduction |
|---|---|---|---|
| EUR | 0.757 (PC2) | **0.036** (PC6) | **95.2%** |
| AJ | 0.984 (PC1) | **0.962** (PC1) | 2.2% |

**EUR: solved.** 0.058% of its variants carried essentially all the callset structure in its PCs.

**AJ: NOT solved, and the AF filter cannot solve it.** No AJ cell was ever evaluated — `AJ/control`
is `wgs_harm`=44, below any GWAS-relevant floor — so the list holds no AJ-derived flags, and
EUR-derived flags do not transfer. AJ also shows **no HWE excess** (0.24× chance expectation).
Leading hypothesis: AJ's split is **real sub-continental structure**, in which case PC1 belongs in
the model and the cohort/phenotype collinearity is a methods limitation, not a bug. Untested — see
`PROJECT_LOG.md` 2026-08-19 for the reference-panel projection that would decide it.

**The "0.748" was EUR's, mislabelled as AJ's — this is settled, and the old note here was wrong.**
Today's EUR baseline is 0.757, within 0.009 of it; today's AJ baseline is 0.984, off by 0.236. The
claimed 0.748 → 0.055 collapse matches EUR's behaviour and contradicts AJ's. Nothing regressed and
the four-callset merge degraded nothing.

**AAC 0.877 and AFR 0.708 are NOT callset structure — they are one BR sample each**, at 39.7σ and
19.6σ. Remove that single point and they fall to 0.016 and 0.121. Cause is BR's ~50% missingness in
a stratum where it is n=1, so `--geno 0.05` cannot protect (1/226 = 0.4%). **No action needed:**
every AAC and AFR contrast is `viable_ge100 = 0`, so those PCs never reach the GWAS. The 17 viable
contrasts are AJ (3) and EUR (14). Revisit only if the ≥100 threshold moves or strata get pooled.
A `--mind` guard would be the wrong fix — BR is ~50% missing everywhere, so it would drop all 95.

**Six of eleven strata have PCs:** EUR, AJ, AAC, AMR, AFR, CAH. CAS/EAS/FIN/MDE/SAS `PRUNE_FAIL` —
plink2 will not LD-prune at those sample sizes. MDE fell 50 → 45 after exclusions and lost its PCs
this run; it had them before.

**Step 4's `COMMON_GENO` is load-bearing.** It must stay below the smallest callset's share of
the cohort (BR is 97/13,334 = 0.0073). At the old `--geno 0.05`, every variant BR lacked stayed
in the common set and all 97 BR samples read as 50% missing — which step 5 then consumed as a
duplicate tie-break signal. `0.005` fixes it. Re-check if a callset under ~0.5% is ever added.

**BR-DSNWGS contributes to roughly half the association set**, because a 97-sample joint call
emits nothing at sites monomorphic in its own donors. Expected, not a defect — but it is a
methods limitation, and BR sits entirely on the AMP-PD side of the primary contrast.

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

3. **Step 6 runs TWICE, with `af_concordance_build.sh` between the passes.** Step 6 reads
   `$MERGED_DIR/exclude_af_concordance.txt`; `af_concordance_build.{py,sh}` writes it, and needs
   step 6's `by_ancestry_qc` output plus the grain to do so. Pass 1 therefore runs *without* the
   exclusion, printing a `NONE` note. That is by design, not a gap — step 6 cannot depend on the
   grain without a circular ordering (grain ← §12 ← manifest ← step 6).

   ```
   6  ancestry QC + PCA (pass 1, unfiltered)
      review/plot_pcs_by_callset.py        -> eta² per PC per stratum. IS the filter needed?
      af_concordance_build.sh              -> exclude_af_concordance.txt  (only if it is)
      read the BY CALLSET PAIR table + the SENTINEL_LOCI tripwire in its log
   6  ancestry QC + PCA (pass 2, picks the list up automatically)
      clinical_core.py §12                 -> PCs changed, so the grain must be rebuilt
   ```

   Why the filter exists — **measured on THIS cohort, 2026-08-19.** Applying the 4,415-variant list
   takes EUR's worst callset eta² from **0.757 to 0.036 (95.2%)**, by excluding 0.058% of its
   variants. The effect is concentrated, not diffuse. Method: the three step-6 QC filters are
   per-variant and independent, so `cohort_<ANC>_qc` minus the blacklisted IDs is *exactly* pass 2's
   output — the test needs no step-6 rerun and no overwrite. See `PROJECT_LOG.md` 2026-08-19.

   **It does NOT fix AJ (0.984 → 0.962), and the earlier claim that it would was wrong.** The
   docstring formerly cited `diag_af_crossstratum`: 147 EUR-flagged variants taking *AJ's* eta² from
   0.748 to 0.055, with 147 random variants leaving it at 0.738. That script exists in no commit,
   and the claim is now contradicted by direct measurement — 93.6% of the July list survives in
   today's, so those 147 are almost certainly inside today's 4,415, and removing them plus 4,268
   others moves AJ by 0.022. The 0.748 was EUR's baseline, not AJ's. **Cite the EUR result above
   instead; it is first-party and reproducible.**
   The files were formerly `06a_af_concordance.*`; the number implied they ran before step 6.

   **The cluster had a stale list, and step 6 applied it.** A 4,587-variant
   `exclude_af_concordance.txt` from 2026-07-28 (3-callset cohort) was picked up by the first
   four-callset pass, so that pass was not the unfiltered baseline it was meant to be. Step 6 now
   **refuses** an exclusion list older than `cohort_merged.bed`. Note the earlier wording of this
   item — "nothing here builds it" — was true of the repo and false of the cluster, which is
   precisely why it went unnoticed. Do not assume local absence means cluster absence.

4. **`clinical_core.py` runs on the CLUSTER ONLY. Do not run it on a laptop.** The two machines
   hold different clinical inputs, and that is how the previous `analysis_grain.csv` acquired
   BR rows no cluster run could have produced:

   | directory | laptop | cluster |
   |---|---|---|
   | `amp-pd-genomics/metadata` | 5 files | 5 (pushed up 2026-08-18; was absent) |
   | `WGS_Harmonization/metadata` | 5 files | 5 — **matched on count, not contents** |
   | `DivCo_HS/metadata` | **2 files** | **5 files** |

   The laptop is three DivCo files short and lacks `all_chrs_merged.psam`, so `wb_dwgs` silently
   drops to 3,010 genomes instead of 13,428. **A file count is not a file list** — that check was
   run, passed, and was wrong; `WGS_Harmonization` was missing `MSBB_biospecimen_metadata.csv` and
   `WGS_sample_QC_info.csv`. Do not "fix" DivCo by pushing the laptop's 2 files up: the cluster's
   set is larger and correct.

   Run it as: `module load python/3.11 && source .venv/bin/activate && python3 clinical_core.py`.
   The system Anaconda py3.9 is on `PATH` and is not the pinned environment.

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
