# ad-pd-gwas — state

AD-vs-PD GWAS across AMP-AD and AMP-PD, restricted to donors with WGS, with **AD defined by
neuropathology**.

**This file is state: what has run, what is true now, and what is still open.** It carries no run
order (`README.md`, and `wgs_core.ipynb` is the executable copy), no rationale (`METHODS.md`), and
no history (`PROJECT_LOG.md` — grep its resolved-anomaly index before investigating anything).

## The goal

1. One harmonized clinical dataset from AMP-AD + AMP-PD. Primary contrast AD vs PD; secondary
   AD and PD vs other dx.
2. WGS QC pipeline on biowulf.
3. GWAS across all contrasts + plots.
4. Follow-up PCA to confirm the callsets align.

## One root

Code and data share a root: on biowulf `/data/CARDPB2/sysbio/wgs`, on a laptop wherever the repo
was cloned. **Nothing in this project contains an absolute path** — `config.sh`,
`clinical_common.py`, `05_excludelist.py`, `ancestry_qc_manifest.py` and `wgs_core.ipynb` each
derive the root from their own location. To move the project, copy the folder. (The last two
exceptions closed 2026-08-21: `01_genotools.sh` had carried three hardcoded paths since 2026-08-11
as a debugging baseline, and the notebook hardcoded the cluster root.)

**Every handoff file has exactly one location.** The clinical side writes to `clinical_core_out/`
and the genotype steps read it there — nothing is copied into `data/`, so nothing can go stale
against the run that produced it. Genotype psams are read **in place** from each callset's `pgen/`
or `joint_calls/`.

| | |
|---|---|
| `clinical_core.py` | goal 1, first half. Runs **once, before step 1**. Reads `data/`, writes `clinical_core_out/`. |
| `analysis_grain.py` | goal 1, second half — §11-13. Runs **once, after step 6**, on its PCs. |
| `clinical_common.py` | paths, readers, and the reconciliation rules both halves share. Imported, never run. |
| `wgs_core.ipynb` | the orchestrator. **Runs ON biowulf** — direct sbatch, no ssh except to helix for transfers. |
| `config.sh` | every path, derived from its own location. `submit.sh` — sbatch wrapper. |
| `scripts/01`–`08` | genotools, normalize, merge, relatedness, excludelist, ancestry QC, GWAS, ctrl-vs-ctrl. |
| `scripts/09_amppd_release.sh` | **release, not pipeline.** Subsets the AMP-PD donors out of step 6 stage C and publishes them to GCS. `MODE=build` on biowulf (sbatch, no network), then `MODE=push` on helix. Written 2026-09-17, **never run**. |
| `scripts/af_concordance_build.{py,sh}` | step 6 **stage B**, in-job. The `.sh` is for re-tuning knobs only. |
| `scripts/ancestry_qc_manifest.py` | step 6 stage E: `retained_samples_manifest.csv`, once per generation. |
| `scripts/gene_annot.py` + `ref/refFlat.txt` | the single source of locus coordinates; a read-only CLI with no pipeline callers. `--at chr19:44908684` |
| `scripts/diag_order.py` | read-only **preflight gate** — variant order vs the reference panel, before step 1 on any uncleared callset. |
| `review/` | QQ/Manhattan, eta² of callset on each PC, the AF-filter before/after. Runs on the cluster against outputs in place. |
| `review/methods_numbers.py` | re-derives every numeric claim in `METHODS.md` from the artifacts and diffs it against the doc. Read-only. Seven checks need cluster files and report `????` on a laptop. |
| `data/**/metadata/` | the 11 clinical files `clinical_core.py` opens. Gitignored — controlled access. |

## Run order

**`README.md` "Run order" is the map; `wgs_core.ipynb` is the executable copy.** Not duplicated
here — two copies of a run order is how the wrong one gets followed.

Two prerequisites fail a whole run rather than a step:

- **`analysis_grain.py` before any step-7 rerun.** A `contrasts.csv` predating 2026-08-21 has no
  callset-skew columns and step 7 exits 3 naming the missing one. Rerun the grain, don't patch it.
- **Both clinical scripts are cluster-only** — open issue 1.

Code reaches the cluster by **git bundle** — there are no GitHub credentials on biowulf. See
`CLAUDE.md` rule 5. (The cluster was **not a git checkout at all** until 2026-09-17, despite this
line having claimed `git pull` since 2026-08-21; `PROJECT_LOG.md` 2026-09-17.)

## Status

**COMPLETE through step 8, figures made, as of 2026-08-21.** What remains is the write-up.
A release step (`scripts/09_amppd_release.sh`, 2026-09-17) exists for handing the AMP-PD subset
to a bucket; it has **not been run**, so no genotype data has left the cluster — open issue 10.

**The cluster became a git checkout on 2026-09-17** (it never was one before, despite the docs).
Four files had drifted. Every step that produced the released results — 04, 06, 06a, 07, 08 — was
verified **logic-identical** to the repo by AST comparison, so the 4,187 list, the association set
and job 28004190 all ran the code in git. `METHODS.md` §6 and §8 stand. `PROJECT_LOG.md` 2026-09-17.
`METHODS.md` is drafted and its numbers are now derived rather than transcribed — run
`python3 review/methods_numbers.py`. A slide deck exists (`ad_pd_wgs_gwas_methods.pptx`, untracked)
and is unreviewed. Of the 14 write-up flags it raised, 9 are closed (`PROJECT_LOG.md` 2026-08-24);
the rest are open issues 6, 8 and 9 below.

| | |
|---|---|
| `cohort_merged` | 172,497,055 variants × 13,334 samples (job 27429821) |
| retained | **12,495** of 13,334 (839 excluded: 499 relative, 319 duplicate, 21 sex) |
| exclusion list | **4,187** variants, HWE excess-over-chance gate ON (was 4,415 ungated) |
| `analysis_grain.csv` | 12,495 rows × 22 cols, **17 viable contrasts** — EUR 14, AJ 3 |
| step 7 GWAS | job **28004190** — 44/44 contrasts ran, λ_GC **1.0175–1.0549** |
| step 8 ctrl-vs-ctrl | job **28037485** — annotated; primary sumstats untouched |
| figures | 17 viable contrasts; scipy λ matches step 7's awk λ to 3 dp on every one |

Step 6 ran as a single pass twice: job **27857727** (7m13s, ungated 4,415) then a rerun with the
gate on (2026-08-20 21:37) producing 4,187 and regenerating stages C/D/E. **That rerun's PCs and
manifests are the live ones.** *(Its job ID is not recorded here — recover it from
`bash scripts/runlog.sh --md`.)*

**The two headline results.** The AF-concordance filter takes EUR's callset eta² from 0.757 (PC2)
to 0.036 (PC6) by excluding 4,187 of EUR's 7,538,809 post-QC variants — 0.06%; and step 8's
annotate-never-subtract licence is what kept APOE-ε4 (P=3.55e-15) in the association set. Both are
written up in `METHODS.md` §6 and §8. `python3 review/methods_numbers.py` re-derives both from the
artifacts.

### Per-callset state

All four cleared genotools ancestry, entered the merge, and survive in the retained set.
`br_dsnwgs` is the one worth stating in full, because its rows were wrong until 2026-08-20:

| | |
|---|---|
| step 1, pre-exclusion | 77 EUR / 14 AJ / 3 AMR / 1 CAH / 1 AFR / 1 AAC over 97 |
| in merge | 97/97, variant arithmetic closes exactly |
| retained after step 5 | **95** — 76 EUR / 13 AJ / 3 AMR / 1 AAC / 1 AFR / 1 CAH |

The old grain's 19 AFR / 67 EUR rows are gone. The laptop's copy is still that stale one — see
Provenance.

`FILTERED.br_dsnwgs_ancestry_umap_linearsvc_predicted_labels.txt` is `LBL_BR` in `config.sh`, read
by steps 4 and 6. The filename says `linearsvc` but the model is XGBoost — a GenoTools naming
artifact, not a description.

**Six of eleven strata have PCs:** EUR 10,135 · AJ 1,518 · AAC 226 · AFR 191 · AMR 185 · CAH 106.
CAS/EAS/FIN/MDE/SAS `PRUNE_FAIL` — plink2 will not LD-prune at those sample sizes. MDE fell 50 → 45
after exclusions and lost its PCs this run; it had them before.

**AJ is closed on power and design, not on eta².** `AJ/AD` is 97 (`wgs_harm`=95 + `divco_hs`=2),
under the 100 floor, so AJ cannot field the primary contrast at all. Its three viable contrasts are
all `within_cohort` — ~96% `wb_dwgs` on both arms — and callset↔phenotype collinearity cannot bias a
contrast whose arms share a callset. The unsolved 0.962 has no downstream consumer. The
sub-continental-structure hypothesis is untested and blocks nothing.

## Standing constraints

True now, and a run breaks if any of them changes. Not open work.

- **`COMMON_GENO=0.005` in step 4 is load-bearing.** It must stay below the smallest callset's share
  of the cohort (BR is 97/13,334 = 0.0073). At `--geno 0.05` every variant BR lacked stayed in the
  common set and all 97 BR samples read as 50% missing — which step 5 then consumed as a duplicate
  tie-break. Re-check if a callset under ~0.5% is ever added.
- **`plink1.9` is a hard dependency of step 6 stage B.** The frequency test is `plink --assoc` (the
  1-df allelic chi-square); plink2 dropped `--assoc` for `--glm`, which is only asymptotically
  equivalent. `MOD_PLINK1` loads unconditionally in `06_ancestry_qc.sh` and
  `af_concordance_build.sh`. If that module goes away the fallback is `--glm` plus a fresh
  equivalence check. `07_gwas.sh` deliberately does *not* rely on plink1.9.
- **`GENOTOOLS_MAX_WORKERS=16` in step 1.** 64 still fails; see the script header for the scan.
- **`--sort-vars` in step 1, and `diag_order.py` before it on any new callset.** README tier 2.
- **Only step 6a may delete variants; step 8 annotates.** Different entitlements, not different
  degrees of caution — `METHODS.md` §6.2 and §8.
- **BR-DSNWGS contributes to roughly half the association set**, because a 97-donor joint call emits
  nothing at sites monomorphic in its own donors. Expected; reported as a limitation.

## Open issues

1. **The clinical side runs on the CLUSTER ONLY — both halves.** The two machines hold different
   clinical inputs, and that is how a previous `analysis_grain.csv` acquired BR rows no cluster run
   could have produced:

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

   Run as `module load python/3.11 && source .venv/bin/activate && python3 clinical_core.py`. The
   system Anaconda py3.9 is on `PATH` and is not the pinned environment.

2. **Two DivCo metadata files have no recorded provenance.** The cluster holds 5; `wgs_core.ipynb`
   §0 records syn IDs for 3. Surfaced 2026-08-21 when the acquisition commands moved out of
   `clinical_core.py`'s header into the notebook. A reproducibility gap, not a blocker.

3. **Nothing in this project records how to obtain the reference data.** Not `config.sh`, not
   either README, not `wgs_core.ipynb`. A fresh clone cannot run step 1 or step 6.

   Missing: the GenoTools ancestry panel (`ref_panel_gp2_prune_rm_underperform_pos_update`,
   a pruned GP2 panel), `ref_panel_ancestry_updated.txt`, the GRCh38 fasta, and the
   GRCh37→GRCh38 liftover chain (`hg19ToHg38.over.chain.gz`). Every source lists **where they
   go** and none says **where they come from**.

   **This issue used to read "the cluster's `README.md` (30 KB) has not been merged in," and
   every part of that was wrong** — established 2026-09-17 once the cluster became a checkout
   and the two copies could actually be compared:

   | the claim | the fact |
   |---|---|
   | cluster README is 30 KB | **15,297 bytes.** The 30 KB was the *diff* size (30,753) |
   | cluster's copy is the larger one | repo's is 26 KB — the cluster's is **smaller** |
   | it holds a §2 the repo lacks | **neither copy has numbered sections.** Cluster headings: Project Overview / Infrastructure / Directory Structure / Datasets / Scripts / Pipeline Order / Sex Crosswalk Notes / Key Gotchas |
   | merging §2 closes it | there is no §2 to merge. The content exists nowhere and must be **written**, not reconciled |

   **Five citations pointed at that phantom numbering** — `config.sh` and `06_ancestry_qc.sh`
   at "README §2", `submit.sh` and `05_excludelist.{sh,py}` at "README §3". All five predate
   the oldest README on disk. Repointed 2026-09-17 (comment-only; `review/drift_check.py`
   confirms LOGIC IDENTICAL on all five).

   The cluster README's one unique asset is its `## Pipeline Order` section — per-callset
   download/reindex/filter/liftover command sequences. Most of that now lives in
   `wgs_core.ipynb`; the chain file and the ancestry panel are the parts nothing covers.
   Preserved on the cluster as `README.cluster.30k.bak` (misleading name, kept as written)
   and `~/README.cluster.diff`.

   **To close this:** write the acquisition down — URLs or `gcloud`/`synapse` commands for the
   four artifacts above — and delete the two backups.

4. **DivCo's source VCF is 0 bytes** on the cluster (`merged.deduped.vcf.gz`). The pgen was derived
   before it was truncated, so nothing is blocked, but DivCo cannot be re-derived from source
   without re-pulling from Synapse.

5. **No age covariate, and it cannot simply be added.** AMP-AD supplies age at **death**, AMP-PD age
   at **baseline/analysis**; forcing one column would silently mix two quantities. Age is the
   dominant confounder for both diseases, and the primary contrast is exactly where the two
   variables are least comparable. Fixing it needs the phenotype track to emit one harmonized age,
   which nothing currently does. Every run prints `age covariate: NOT FOUND in grain`.
   `METHODS.md` §10.

6. **The callset-skew columns are emitted but have NOT been run — and the cause is now known.**
   §13 writes `max_callset_delta` / `worst_callset` / `callset_one_sided` into `contrasts.csv` as of
   2026-08-21, and step 7 carries them into `gwas_summary.csv` — but job 28004190 predates them, so
   the summary on disk has none of the three.

   **The cluster's `analysis_grain.py` had no `callset_skew()` function at all** until 2026-09-17:
   written on the laptop 2026-08-21, never rsynced up. That made a **deadlock** nobody had hit,
   because step 7 has not been rerun since — the cluster's `07_gwas.sh` *is* logic-identical to the
   repo's, so it exits 3 on the missing columns, while the cluster's grain could not emit them.
   Fixed by the checkout; regenerate the grain before any step-7 rerun. `PROJECT_LOG.md` 2026-09-17.

   `confound_tag` itself is unchanged and still means *program*, not callset; `METHODS.md` §10
   carries the measurement and the two rejected designs.

7. **Step 0 (VCF→pgen for BR-DSNWGS) has no script.** It ran as notebook cells; `wgs_core.ipynb` §1
   is the whole record.

8. **Seven `METHODS.md` numbers can only be settled on the cluster.** Run
   `python3 review/methods_numbers.py --strict` there; on a laptop those seven report `????` and the
   rest come back clean. **The script itself was never on the cluster** — it arrived 2026-09-17 with
   the git conversion, which is why this has stayed open. It can now actually run. Each one is read-only, and the script names the file it could not open:

   | what it settles | source it needs |
   |---|---|
   | the **94-genome gap** between §1's table (13,428) and §3's merge (13,334), per callset | the `.fam` files named in `merge_list.txt` |
   | §6.4's **exact** excluded share of EUR variants — every doc now says 0.06%, sourced but rounded | `by_ancestry_qc/unfiltered/cohort_EUR_qc.bim` vs the stage-C bim |
   | the live exclusion list is 4,187 | `exclude_af_concordance.txt` |
   | **§6.3's gate table omits EUR `divco_hs`** — evaluated and below the bar, or never evaluated for want of controls? Silence is not an answer | the cluster `analysis_grain.csv`; the laptop's stale 11,918-row copy is refused by name |
   | the **step-6 job ID** behind the live 4,187 list, which every §6 number traces to | `bash scripts/runlog.sh --md` |

   Nothing here is a pipeline defect. They are write-up provenance, and they block the paper, not a
   run.

9. **§6.4's "~7× enriched for duplicate-pair discordance" is unsourced — measure it or delete it.**
   It is the only *non-circular* evidence that the 4,187 flagged variants are technical rather than
   real; every other §6 number is measured on the same frequency signal the filter was built from.
   No numerator, denominator, pair count, script or job ID exists anywhere in the repo.
   `python3 review/methods_numbers.py --discordance` computes it on the cluster (`plink2
   --sample-diff` over the KING duplicate pairs at kinship ≥ 0.354, flagged sites vs the rest).
   **If it cannot be run, the sentence comes out of §6.4** — an unsourced multiplier carrying the
   filter's whole justification is worse than no sentence.

10. **`scripts/09_amppd_release.sh` has never been run, and two of its assumptions are
    unverified on the cluster.** The script is written, syntax-checked, and exercised end-to-end
    against a stubbed `plink2` — its selection, join, accounting and manifest logic all work, and
    all three fatal guards fire. What a laptop cannot settle:

    - **The gcloud module name on biowulf.** `MOD_GCLOUD` defaults to `google-cloud-sdk` and the
      script falls back to whatever `gcloud` is already on `PATH`. If neither resolves, the build
      exits 1 rather than writing a manifest with no CRC32C — `ALLOW_NO_HASH=1` accepts size-only
      verification, which is a weaker check than this release should get. Settle it with
      `module spider google-cloud-sdk` before the first build.
    - **Whether every AMP-PD sample reconciles.** The build refuses to ship unless each one lands
      in exactly one stratum fileset. AMP-PD donors in a stratum step 6 never wrote a fileset for
      (the `PRUNE_FAIL`/sub-2-sample strata) would trip it. Nothing on the laptop can say whether
      any exist, because the laptop's `analysis_grain.csv` is the stale 11,918-row one. The first
      real build answers it in one line.

    No bucket is hardcoded anywhere; `GCS_DEST` is required at push time. Redistribution of
    individual-level AMP-PD genotypes is governed by the DUA — the script will not upload to a
    bucket it cannot prove is private, and that refusal has no override flag.

**Closed since 2026-08-19**, all with their reasoning in `PROJECT_LOG.md`'s index: the pheno/covar
duplication (§13 is the sole definition), the ctrl-vs-ctrl duplication (`review/mask_cohort_artifacts.py`
deleted), step 6's two-pass shape (one pass), the HWE excess-over-chance gate (4,415 → 4,187), and
the sentinel tripwire (deleted, both call sites).

## Provenance

`scripts_archive/ad-pd-gwas-jul28/` on the cluster holds the July version of every script, plus 15
files never carried over. Of the 12 scripts shared with this repo, 8 are byte-identical;
`config.sh`, `06_ancestry_qc.sh`, `07_gwas.sh` and `05_excludelist.py` differ, with this copy
authoritative.

**Remote added 2026-08-21: `SysBio-FAIRplex/amp-ad-pd-wgs-gwas` (private).** Before that the entire
history was 12 commits on one laptop disk. The first push was preceded by a history rewrite — a
notebook had stored a subject-level output table — so no controlled-access data has ever reached the
remote. `scripts/nb_guard.py` plus a `pre-commit` hook now refuse any notebook carrying stored
outputs; install the hook after a clone (README tier 2).

`clinical_core_out/` and `results/` remain gitignored. **The laptop's `analysis_grain.csv` is the
stale 11,918-row one (BR-DSNWGS as 19 AFR) and must never be rsynced upward** over the cluster's
correct 12,495-row grain.

**`results/pca/*_eta2.csv` and `af_filter_effect*.csv` are the one exception in `.gitignore`**, with
the reason written there. They are aggregate (one row per stratum, no IIDs, no genotypes), and the
0.984 / 0.757 baseline previously survived being overwritten only because the numbers had been typed
into `PROJECT_LOG.md` by hand. "Regenerable from code + data" was true in principle and false in
practice: regenerating a *baseline* means re-running the step in a state that no longer exists. The
negation is narrow on purpose — `results/retained_samples_manifest.csv` is one row per sample and
stays out.
