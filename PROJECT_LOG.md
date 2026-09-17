# Project log — ad-pd-gwas

Append-only. Newest entry first. **Nothing here is ever rewritten or deleted**, including
entries that later turn out to be wrong — a hypothesis that was ruled out is the most
expensive thing to rediscover, and the only place it can be recorded is a log.

Three documents, one job each. Keep them apart:

| doc | answers | written by |
|---|---|---|
| `README.md` | what the pieces are and in what order they run | rewritten in place |
| `METHODS.md` | why each choice was made | rewritten in place |
| `HANDOFF.md` | what is true now: status and open issues | rewritten in place |
| `RUNLOG.md` | which jobs ran, and how they ended | generated — `bash scripts/runlog.sh --md > RUNLOG.md` |
| `PROJECT_LOG.md` (this file) | what we did, why, and what we ruled out | appended |

Entry format: date, a one-line title, then whichever of **Did / Found / Changed / Next**
apply. Record dead ends explicitly — they are the point.

**"Where we are right now" stays at the TOP** and is the only part of this file that gets
rewritten. Everything below it is append-only history.

---

## Resolved-anomaly index

**Grep this table before investigating anything.** Every row is a symptom that was already chased
to a verdict; re-deriving one is the most expensive failure mode in this project. Search by locus,
number, or symptom, then read the dated entry.

| symptom | verdict | entry |
|---|---|---|
| CR1 `chr1:207521012:T:C` excluded by 6a | correct — `divco_hs` calls it in **156 of 242 alleles** (64.5%) while the same samples are 242/242 at the LRRK2 site. Site-specific dropout, not a frequency difference. Passes HWE in all three testable callsets | 2026-08-20 (later) |
| LRRK2 `chr12:40227079:C:T` excluded by 6a | **wrong, and un-excluded by the HWE gate.** Never flagged on frequency (max spread 0.034 < 0.05); excluded on one chance-level HWE hit in AJ/`wb_dwgs`, while the same callset with 3,064 EUR controls passes | 2026-08-20 (later), (night, later) |
| all 97 BR-DSNWGS samples labelled CAH, 0.97 model accuracy, exit 0 | variant ORDER mismatch. GenoTools aligns to the panel by column POSITION; BR's pgen was `1,10,11,…,2,20` from alphabetical concatenation. Fixed with `--sort-vars` | 2026-08-14, 2026-08-12 |
| genotools dies with a bare SLURM `ExitCode 1:0`, no OOM, no MaxRSS | worker pool sized from the NODE, not the allocation; `RLIMIT_NPROC` 1024 is per-user node-wide. **64 workers still fails**; 16 completes. Not thread env vars — those were superstition | 2026-08-12 (three entries) |
| eta² of 0.748 attributed to AJ | **EUR's baseline, mislabelled.** Today's EUR is 0.757 (within 0.009); AJ is 0.984 (off by 0.236). The `diag_af_crossstratum` script it was cited from exists in no commit | 2026-08-19 (premise), 2026-08-20 (later) |
| "~748 duplicate genomes" | never a genotype-duplicate count — it is the **enrollment** overlap (846 donors in >1 cohort, 748 DivCo↔ROSMAP). 309 have >1 genome; KING finds 302 clusters; 319 dropped | 2026-08-18, 2026-08-17 |
| every BR-DSNWGS sample reads 50% missing | expected at `--geno 0.05`: BR is 97/13,334 = 0.0073 of the cohort, so every variant it lacks stayed in the common set. `COMMON_GENO=0.005` fixes it. Step 5 had been consuming it as a duplicate tie-break | 2026-08-16, 2026-08-17 |
| step 6 pass 1 was not an unfiltered baseline | a 4,587-variant `exclude_af_concordance.txt` from 2026-07-28 existed **only on the cluster** and was silently applied. Local absence ≠ cluster absence | 2026-08-17 |
| AJ PC1 eta² 0.984, barely moved by the filter (→0.962) | **open, and the AF filter cannot fix it.** No AJ cell was ever evaluated: 6a needs TWO callsets above `MIN_CELL=100` in one cell, and AJ's best is `AJ/control` at `wb_dwgs=638` but `wgs_harm=44`. So the list holds no AJ-derived flags and EUR-derived ones do not transfer. AJ also shows no HWE excess (0.24× chance). Leading hypothesis: real sub-continental structure. Closed on power and design, not on eta² | 2026-08-19 (6a ran clean), (premise) |
| AAC eta² 0.877, AFR 0.708 | **one BR sample each**, at 39.7σ and 19.6σ. Remove that point and they fall to 0.016 and 0.121. No AAC/AFR contrast is viable, so neither PC reaches a GWAS | 2026-08-19, 2026-08-18 |
| APOE reaches significance in the control-vs-control scan | **expected, not an artifact.** The control arms are differentially screened — AMP-AD as cognitively normal, AMP-PD for PD and not AD. Subtracting would delete ε4 at P=3.55e-15 from the study's cleanest AD contrast | 2026-08-21 (step 8) |
| `EUR PD_vs_DLB` is `within_cohort` (Δ=0.0) yet lost 353,068 variants to differential missingness | `confound_tag` measures PROGRAM and pools `wb_dwgs` with `br_dsnwgs`. BR contributes **no DLB at all** and sits at ~50% missingness. `n_diffmiss_excluded` is the signal; one-sidedness alone fires on 5 of 5 and predicts nothing | 2026-08-21 (step 7) |
| 245 of 279 sentinel-window exclusions are MHC | per-variant technical: the two natively-called callsets agree near 0 while the lifted one does not (`chr6:32474706` 0.008/0.005 vs 0.207–0.213). Licenses deletion under 6a. The flag RATE is a methods number, not a decision input | 2026-08-20 (night, MHC) |
| the `--assoc` swap flagged one variant fewer than the Python test | `.assoc` prints CHISQ to four SIGNIFICANT figures, so 25.005 lands as "25" and `25.0 > 25` is false. Threshold on **P**, which is exponential and ~1000× finer at the boundary. Cost exactly `chr6:32555808:T:C` | 2026-08-20 (later still), (evening) |
| a bare `plink2` dies with PermissionError mid-job | on biowulf `plink` and `plink2` are one module family: loading plink/1.9 UNLOADS plink/6-alpha. Call both by absolute path | 2026-08-20 (later still) |
| `07_gwas.sh` put BR-DSNWGS on the AMP-AD side | found 2026-08-11, fixed since — `AMPPD_CALLSETS = {wb_dwgs, br_dsnwgs}` | 2026-08-11 |
| `clinical_core_out/` missing on the cluster | it had never been pushed | 2026-08-11 |
| the sentinel tripwire | **deleted 2026-08-20, both call sites.** It gated nothing (stage C applies the list later in the same job) and its 9-gene scope made both its positive and negative results uninformative. Do not reintroduce | 2026-08-20 (night), (evening) |
| `review/mask_cohort_artifacts.py`, `diag_cah.sh`, `diag_threads.py`, `compare_pcs.py` | **deleted 2026-08-21.** Reasons per file in that entry | 2026-08-21 (git remote) |
| subject-level data in git | one notebook carried 18 stored outputs, one of them an `individual_id` × callset table. Scrubbed from history before the first push; `scripts/nb_guard.py` now blocks it | 2026-08-21 (git remote) |
| three AJ contrasts report `lambda_gc` exactly **0.0000** | **non-convergence, not deflation.** `AD_vs_DLB`, `PD_vs_AD` and `control_amppd_vs_control_ampad` are AJ's cross-callset contrasts; every variant returns `ERRCODE=UNFINISHED` with *P* ≈ 1, so the median chi-square is ~0. AJ PC1 is η² 0.962 on callset and enters as a covariate, so a covariate separates the arms almost perfectly. All three are below the viability floor and none is reported | 2026-08-24 |
| METHODS §10's "roughly six times any cross-program contrast" | **false as written.** EUR `PD_vs_AD` removed 387,207 to differential missingness, more than `PD_vs_DLB`'s 353,068, and was omitted from the comparison set. The mechanism argument is unaffected; only the multiplier failed. Rewritten as "largest within-program removal, same magnitude as the primary cross-program contrast" | 2026-08-24 |
| METHODS §6.4's "0.058% of variants" vs PROJECT_LOG's "0.05%" | both are the same quantity with no stated denominator. 0.058% = 4,385 / 7,538,809, measured on the **ungated 4,415** list; the live gated list is 4,187. Both round to 0.06% of EUR's post-QC variants, which is now what every doc says | 2026-08-24 |

---

## Where we are right now

**The pipeline is COMPLETE through step 8 and the figures are made, as of 2026-08-21.** Steps 0–8
have all run to completion on the four-callset cohort, and `review/plot_gwas.py` has produced QQ +
Manhattan for all 17 viable contrasts (scipy λ matches step 7's awk λ to 3 dp on every one). What
remains is the write-up. `METHODS.md` is drafted and, as of 2026-08-24, its numbers are derived by
`review/methods_numbers.py` rather than transcribed; a slide deck exists and is unreviewed, and the
14 flags it raised are 9 closed, the rest `HANDOFF.md` issues 6, 8 and 9. Read `README.md` for the run order
and `HANDOFF.md` for status and open issues; this block is the short version.

### The chain that ran

| step | job | outcome |
|---|---|---|
| 3 merge | 27429821 | 172,497,055 variants × 13,334 samples |
| 5 excludelist | — | **12,495 retained** (839 excluded: 499 relative, 319 duplicate, 21 sex) |
| 6 ancestry QC | 27857727, then a rerun | one pass, five stages. Ungated list 4,415 → **gated 4,187** |
| grain (§11–13) | — | 12,495 rows × 22 cols, **17 viable contrasts** (EUR 14, AJ 3) |
| 7 GWAS | **28004190** | **44/44 contrasts ran**, λ_GC **1.0175–1.0549** |
| 8 ctrl-vs-ctrl | **28037485** | annotated; primary sumstats untouched |

Six of eleven strata have PCs (EUR, AJ, AAC, AMR, AFR, CAH); CAS/EAS/FIN/MDE/SAS `PRUNE_FAIL` at
their sample sizes. Only EUR and AJ field viable contrasts.

### The two results that matter

**1. The AF-concordance filter works, and it is measured on the live list.** Applying the
4,187-variant exclusion takes EUR's worst callset eta² from **0.757 (PC2) → 0.036 (PC6), 95.3%** by
removing **0.06% of EUR's 7,538,809 post-QC variants** — and the whole profile collapses, not just
the peak: all ten PCs land at or below 0.036 and PC2 itself falls to 0.011. (The "0.058%" quoted
before 2026-08-24 was measured on the *ungated* 4,415 list, where 4,385 of them were present in
EUR: 7,538,809 → 7,534,424. Both round to 0.06%; the exact gated figure is the one
`review/methods_numbers.py` derives on the cluster.) The HWE excess-over-chance gate cost this
nothing — the ungated 4,415 list
gave 0.036 / 95.2%, so the frequency channel carries the whole effect. Tables are versioned at
`results/pca/af_filter_effect_gated4187*.csv`, the first files ever committed under `.gitignore`'s
eta² negations. **AJ is NOT solved (0.984 → 0.962) and this filter cannot solve it** — no AJ cell was
ever powered enough to evaluate, so the list holds no AJ-derived flags. Leading hypothesis is real
sub-continental structure, which would make PC1 a legitimate covariate rather than a bug. Untested,
and it blocks nothing: all three AJ contrasts are `within_cohort`, where callset↔phenotype
collinearity cannot bias the comparison.

**2. Step 8's annotate-never-subtract licence was confirmed, and it saved the study's headline
finding.** All 4 flagged genome-wide hits in `EUR AD_ampad_vs_control_ampad` — the *within-cohort* AD
contrast, the cleanest in the study — are APOE, including **rs429358/ε4 at P=3.55e-15** and
rs7412/ε2. Two of the three genome-wide control-vs-control hits are `HLA-DQB1`, 30 kb from
`HLA-DRB1/DRB5`. Both are the documented screening asymmetry (AMP-AD controls assessed cognitively
normal; AMP-PD controls screened for PD only), not batch. **`.ccfilt.tsv` must not be used for the AD
contrasts — it deletes ε4.** Meanwhile every high-signal contrast came through with zero
control-flagging (`AD_vs_PSP` 2,319/0, `PD_vs_control` 2,523/0, `PD_amppd_vs_control_amppd` 2,508/0),
which is the corroboration that 6a already removed the cohort artifact. Full detail: `HANDOFF.md`
known issue 2 and the 2026-08-21 entry below.

### Next, in order

1. **Run `python3 review/methods_numbers.py` on biowulf.** It settles the seven checks a laptop
   cannot, including the 94-genome gap and §6.4's exact denominator. Read-only.
2. **Measure or delete §6.4's "~7× duplicate-pair discordance."** It is the only non-circular
   evidence that the 4,187 flagged variants are technical, and it has no numerator, denominator or
   job ID anywhere. `HANDOFF.md` issue 9.
3. **Name two variants** for the write-up: `python3 scripts/gene_annot.py --at chr19:44912456 --at
   chr3:106666502`. The first decides whether the APOE signal is described as extending into APOC1;
   the second is the one ctrl-vs-ctrl hit with no AD/PD story and the only plausible true artifact.
4. **Write-up.** The MHC flag-rate command in `HANDOFF.md` is a methods number, deliberately not a
   gate — see the 2026-08-20 decision entry.

### Open, and honestly open

- **No age covariate**, and it cannot be forced: AMP-AD gives age at death, AMP-PD age at baseline.
  The dominant confounder for both diseases. `HANDOFF.md` issue 5.
- **`confound_tag` measures program, not callset.** `within_cohort` is not a clean bill of health —
  `EUR PD_vs_DLB` wears it while losing 353,068 variants to differential missingness, because
  BR-DSNWGS is 71 PD / 0 DLB at ~50% missingness. Annotated as of 2026-08-21; `HANDOFF.md` issue 6.
- **BR-DSNWGS contributes to ~half the association set** — a 97-donor joint call emits nothing at
  sites monomorphic in its own donors. Expected, but a methods limitation, and BR sits entirely on
  the AMP-PD side of the primary contrast.
- **The order check is still not a pre-flight gate in `01_genotools.sh`.** The CAH bug exited 0,
  wrote output, and scored 0.97 model accuracy while corrupting every sample. Run
  `python3 scripts/diag_order.py <panel>.bim <callset>.pvar` by hand on any new callset until it is.
- **Upstream GenoTools fixes unlanded** in `dvitale199/GenoTools`: the positional column alignment
  (`ancestry.py:247`), worker sizing, and making `get_common_snps` report its overlap.
- **Cross-callset sex-file inconsistency, unverified:** the other three callsets were sex-updated
  from `<dataset>/metadata/` rather than the corrected `clinical_core_out/`. Never re-checked.
- **The cluster's 30 KB `README.md` still is not merged in**, and its §2 (reference-data acquisition)
  is the part the repo copy lacks — `06_ancestry_qc.sh` cites that section by number.
- ~~**One implementation duplication left:** `review/mask_cohort_artifacts.py`~~ — CLOSED
  2026-08-21, the file is deleted. `scripts/08_ctrl_ctrl_filter.py` is the only implementation.
  The pheno/covar duplication was resolved 2026-08-20. No known duplications remain.

### Traps that have each cost a run

- **The laptop's `analysis_grain.csv` is the stale 11,918-row one** (BR as 19 AFR). The correct
  12,495-row grain exists **only on the cluster**, and `clinical_core_out/` sits at the project root
  beside the code — so **any laptop→cluster rsync of the root must name files explicitly**, never
  sync the directory. Fourth face of the stale-artifact bug.
- **`rsync` without `--delete` cannot express a deletion.** Retired scripts must be removed by hand;
  this has caused three stale-artifact bugs. Local absence ≠ cluster absence.
- **Both clinical scripts are cluster-only.** The laptop holds different clinical inputs (DivCo is
  3 files short) and silently produces a smaller, wrong table. A file *count* is not a file *list*.
- **`Bash` keeps its working directory between calls.** On 2026-08-20 a stale `cd scripts` made
  `ls`, `git ls-files` and `git check-ignore` all return true answers to the wrong question, and
  produced a confident, wrong claim that `ref/` was untracked. Use absolute paths to check existence.
- **Read columns by header name.** `plink2 --freq` puts `PROVISIONAL_REF?` at column 5, and a
  positional read of an `.afreq` briefly hid the CR1 call rate. This also bit a script reading a
  summary *it had written itself* (2026-08-21), the moment two columns were inserted.
- **Duplicates are settled and the "~748" is sourced:** `clinical_core.py` §8's multi-cohort
  ENROLLMENT count (846 donors in >1 cohort, 748 DivCo↔ROSMAP), never a genotype-duplicate count.
  Chain: 846 enrolled twice → 309 have >1 genome → KING finds 302 clusters over 621 genomes → 319
  drops. Enrollment overlap is a superset of sequencing overlap.
- **Step 4's `COMMON_GENO=0.005` is load-bearing** — it must stay below the smallest callset's share
  of the cohort (BR is 97/13,334 = 0.0073). Re-check if a callset under ~0.5% is ever added.

---

## 2026-09-17 (later) — the cluster was never a git checkout. Four weeks of docs asserting otherwise. Four files drifted; the results are unaffected.

**Did.** Tried to `git pull` on `/data/CARDPB2/sysbio/wgs` to deliver step 9 and got
`fatal: not a git repository`. **The cluster directory has never been a checkout.** Converted it in
place, measured the drift against the repo, and corrected the three documents that asserted the
opposite.

**The documentation bug, which is the real finding.** `CLAUDE.md` rule 5 said the remote "was
created 2026-08-21 **and retired the two-rsync dance**"; `HANDOFF.md` line 60 said "Code reaches the
cluster by `git pull` (remote added 2026-08-21)"; the 2026-08-21 log entry said the same. All three
were true about the *laptop* — the remote exists and holds all 50 tracked files — and false about
the cluster, where the workflow was never adopted. **A rule written in the past tense about an
intention is indistinguishable from a rule written in the past tense about a fact.** Nothing in the
project could have caught this, because every check that would have run `git` on the cluster was
itself premised on the cluster being a repo.

**No GitHub credentials on biowulf**, and GitHub no longer accepts password auth, so fetching from
the remote there fails outright. Transport is now a **bundle** — `git bundle create` on the laptop,
`scp` to helix, `git fetch ~/adpd.bundle main:refs/remotes/origin/main` on the cluster. This is not
the loose-file `scp` rule 5 forbids: a bundle carries the whole commit graph, so the cluster ends up
a genuine checkout at the true commit, and `git rev-parse HEAD` — which `09_amppd_release.sh` stamps
into every release — returns the real hash. Loose-file `scp` into the checkout would have been worse
than useless: `rev-parse` would still succeed and stamp releases with a commit that does **not**
contain the script that built them.

**Two traps in the conversion itself, both worth keeping.**

- **`git show origin/main:.gitignore > .gitignore` creates the file before the command runs.** The
  fetch had silently failed, `git show` errored, and the redirect still left a **zero-byte
  `.gitignore`** — on top of a 4.2 TB data tree. Harmless under `status`, which reports `data/` as a
  single untracked directory, and a loaded gun under any `git add -A`. Use
  `git show … > f.new && mv f.new f`, and gate on `git rev-parse --verify origin/main` first.
- **`git reset origin/main` (mixed, no `--hard`) is the whole conversion.** It sets HEAD and the
  index and touches **no** working-tree file, so the cluster's copies survive and `git status`
  becomes a drift report. `--hard` here would have destroyed the 30 KB `README.md` of issue 3 and
  four files' worth of divergence before anyone could look at them.

**The drift: 16 modified files. Method matters, because two naive checks both gave wrong answers.**

Line counts first suggested the drift was one commit — `ccebea8` "Strip narrative commentary from
pipeline scripts" — because **seven files matched its per-file line counts exactly** (151, 126, 212,
59, 203, 16, 12). Suggestive, not proof, and the second check contradicted it: stripping whole-line
comments reported `07_gwas.sh`, `08_ctrl_ctrl_filter.py` and `af_concordance_build.py` as code
changes. Both were wrong in opposite directions.

What settled it, per file type:

| type | test | why the cruder test failed |
|---|---|---|
| `.py` | **AST comparison with docstrings stripped** | comments are not in an AST at all; `ccebea8` rewrote module docstrings, which a text diff counts as code |
| `.sh` | comment strip **including trailing** comments, quote-aware | `ccebea8` realigned trailing comments on `07_gwas.sh`'s eight knob lines. Values byte-identical (`20`, `0.05`, `1e-4`, `0.02`, `5e-8`, `5`); a whole-line-only stripper called it CODE CHANGED |

By AST, `ccebea8` is **comment-and-docstring-only for all four Python files it touched**, including
`af_concordance_build.py`.

**The verdict, and it is the one that mattered: every step that produced the released results ran
the repo's logic.**

| file | verdict |
|---|---|
| `scripts/06_ancestry_qc.sh` | **LOGIC IDENTICAL** |
| `scripts/af_concordance_build.py` | **LOGIC IDENTICAL** — step 6a, the only licensed deletion |
| `scripts/07_gwas.sh` | **LOGIC IDENTICAL** (trailing comments only) |
| `scripts/08_ctrl_ctrl_filter.py` | **LOGIC IDENTICAL** |
| `scripts/04_relatedness.sh`, `af_concordance_build.sh`, `clinical_core.py`, `submit.sh` | LOGIC IDENTICAL |

So the 4,187 exclusion list, the association set, job 28004190 and the ctrl-vs-ctrl annotation were
all produced by the code in git. **`METHODS.md` §6 and §8 stand.** This was not the expected answer
and it was worth the work to establish rather than assume.

**FOUND — `analysis_grain.py` on the cluster has no `callset_skew()` at all, and this is the cause
of `HANDOFF.md` issue 6.** The function, `ONE_SIDED_MIN`, and the call site are all present in the
repo and absent on the cluster: §13's callset-skew columns were written 2026-08-21 and the file was
never rsynced up.

The consequence is a **deadlock**, undiscovered only because step 7 has not been rerun since. The
cluster's `07_gwas.sh` is logic-identical to HEAD, so it *does* exit 3 when `contrasts.csv` lacks
`max_callset_delta`/`worst_callset`/`callset_one_sided` — while the cluster's `analysis_grain.py`
cannot emit them. Regenerate the grain, step 7 exits 3; regenerate again, identical failure. Issue 6
recorded the symptom ("emitted but have NOT been run") and attributed it to job 28004190 predating
the columns. That was true and incomplete: the code to produce them was never on the machine.

**Also never on the cluster: `review/methods_numbers.py`.** It came down as a *deletion* in the
conversion, meaning it had never been transferred. `HANDOFF.md` issue 8 instructs
`python3 review/methods_numbers.py --strict` **on the cluster** to settle seven `METHODS.md`
numbers. The script was not there. Same for `review/plot_af_filter_effect.py`, `scripts/nb_guard.py`
and `scripts/hooks-pre-commit` — so the pre-commit notebook guard never existed cluster-side either.

**Changed.** Cluster converted to a checkout at `4f73ac3`. Twelve deleted files restored (step 9's
script, `methods_numbers.py`, `plot_af_filter_effect.py`, `nb_guard.py`, `hooks-pre-commit`,
`METHODS.md`, `wgs_core.ipynb`, the `SUMSTATS_*` pair, `demo_sample_check/`). Four drifted files
checked out. `CLAUDE.md` rule 5 rewritten to say what is true and to record the bundle transport.

**Left alone deliberately: `README.md`.** The cluster's 30 KB copy is issue 3 — its §2
(reference-data acquisition) is what the repo version lacks, and `06_ancestry_qc.sh` cites that
section by number. Backed up to `README.cluster.30k.bak`; the merge is still open.

**Next.**
1. Merge the cluster README's §2 into the repo copy — closes issue 3, and it is now a readable
   `git diff` rather than an assertion.
2. `python3 review/methods_numbers.py --strict` on the cluster. It can finally run; issue 8.
3. Regenerate the grain with the repo's `analysis_grain.py` before any step-7 rerun; issue 6.
4. A PAT or SSH key on biowulf retires the bundle step.

**Dead end, recorded so it is not retried:** comparing line counts against a suspected commit. Seven
exact matches out of nine felt conclusive and was not — `analysis_grain.py` matched nothing (137 vs
79) and turned out to be the only file with a substantive difference.

---

## 2026-09-17 — a release step for the AMP-PD subset: `scripts/09_amppd_release.sh`. Written, tested against a stub, NOT run.

**Did.** Wrote step 9 — subset the AMP-PD donors out of the final QC'd association set and publish
them to a GCS bucket. Documented it in `README.md` as **tier 1b, release**, and in `HANDOFF.md` as
open issue 10. Nothing upstream changed; no pipeline artifact is read back by it.

**What it ships, and why that fileset.** `data/merged/by_ancestry_qc/cohort_<ANC>_qc` — step 6
**stage C**, post per-ancestry QC and post the AF-concordance exclusion list. That is the exact
variant set step 7 tested, which is the whole reason it is the thing worth releasing. Deliberately
not `cohort_merged` (pre-QC) and not `unfiltered/` (the baseline twin, which exists to measure the
filter, not to be shipped).

**One fileset per stratum, not one merged fileset.** The QC is per stratum, so each stratum has its
own variant set; concatenating them would produce a fileset no GWAS in this study ran on. This was
the choice put to the user and the one taken.

**Not re-filtered after subsetting, and this is the decision most likely to be second-guessed.**
Dropping the AMP-AD donors leaves some retained variants monomorphic or low-MAF within the AMP-PD
subset. They stay. Re-applying `--maf`/`--geno` would yield a variant set that silently disagrees
with the published sumstats — a released pgen and a released sumstat that do not describe the same
variants is exactly the kind of drift `review/methods_numbers.py` was written to stop. `SUBSET_MAF`
and `SUBSET_GENO` opt in, and the generated release README states which was used either way.

**Who counts as AMP-PD: read, not re-derived.** `source_callset` in {`wb_dwgs`, `br_dsnwgs`}, taken
**by header name** from step 6 stage E's `retained_samples_manifest.csv`. That file is already the
one implementation of per-sample callset membership. Re-deriving it from the four genotools label
files would have been a second one, which is the duplication rule, and it would have been the
easier thing to write.

**Two hosts, because one of them has no network.** `MODE=build` is an sbatch job on biowulf and
never touches the network. `MODE=push` runs on helix and refuses to start inside a SLURM
allocation — a push attempted from a compute node fails looking like an auth or bucket problem,
which is an expensive misdiagnosis to walk into.

**Three refusals, none of them overridable:**

| refusal | the failure it exists for |
|---|---|
| partial release | every AMP-PD sample in the retained manifest must land in exactly one fileset. A stratum with AMP-PD donors but no step-6 fileset is a hole that looks *identical* to a clean run unless something counts. The check names the short strata |
| unverified-private bucket | uniform bucket-level access on, public access prevention enforced, no `allUsers`/`allAuthenticatedUsers`. **If the IAM policy cannot be read, that is also a refusal** — absence of a public binding in a listing that failed is not evidence of absence (rule 3) |
| unverified upload | `cp` exiting 0 does not prove every object arrived; the failure it misses is a file never in the argument list. Stage E re-lists the bucket and compares presence, size and CRC32C against `MANIFEST.tsv` |

The bucket refusal has no flag on purpose. This is individual-level genotype data under the AMP-PD
DUA, and that decision does not belong to a shell variable.

**No bucket is hardcoded.** `GCS_DEST` is required at push time. The sumstats bucket
(`gs://sysbio-gwas/results`, 2026-09-15 entry) holds aggregate statistics; this is a different data
class and the destination is a deliberate choice each time.

**CRC32C, not MD5** — carried straight from the 2026-09-15 finding that 42 of the 46 sumstats
objects came back MD5-less because they were composite uploads. The build records CRC32C and
SHA-256; stage E compares CRC32C. If `gcloud` cannot be found the build **exits 1** rather than
writing a manifest with no hash — `ALLOW_NO_HASH=1` accepts size-only, and says plainly that size
agreement is not content agreement.

**Three parse bugs found and fixed before they could run.** The first two are rule-4 shaped:

- The upload list was a brace-expanded glob (`amppd_*.{pgen,pvar,psam,bed,...}`). With `nullglob`
  off, an unmatched pattern — `amppd_*.bed` on a pgen build — expands to a literal nonexistent path
  that `gcloud` then errors on, and the glob would also have swept up the per-stratum `.plink.log`
  files. Now built from `MANIFEST.tsv`, which is already the authoritative list.
- Stage E parsed `gcloud storage ls -L` prose positionally. That output is human-facing, has
  changed shape between `gsutil` and `gcloud storage`, and renders sizes as `12345 (12.05 KiB)` —
  a positional read of which yields `KiB)`. Both it and the bucket preflight now read
  `--json`/`--format=json` through a small python parser that looks its fields up **by name at any
  depth**, so the camelCase/snake_case/`iamConfiguration`-nested variants gcloud has emitted across
  releases all resolve. An empty listing after a successful `cp` is a loud failure, not
  "everything missing".
- **The third is the one worth remembering.** Both embedded python parsers are single-quoted shell
  arguments, so they cannot contain a single quote — and the f-strings I wrote used `\"` to get
  double quotes inside the expression, which is a `SyntaxError`. It survived my first test only
  because the test harness unescaped the source before running it, i.e. **the test was not running
  the code the script would run.** Caught by the stubbed push instead, where the preflight printed
  a traceback and then reported all three checks as FAIL.

  It failed **closed** — a broken parser yields empty values, empty values fail the check, and the
  upload is refused. That is the right direction and it is not an accident; it is what "a check
  that cannot run must say so loudly" buys you. But it would equally have refused a perfectly
  configured bucket, and the traceback was the only thing distinguishing the two. Both parsers now
  build their output by concatenation, with no f-string and no escaped quote anywhere in the file.

**Tested against a stubbed `plink2` and a stubbed `gcloud`** — three mock strata, a mock manifest,
AMP-PD in two of them. Nine cases, all passing:

| # | case | result |
|---|---|---|
| 1 | staging tree exists, no `OVERWRITE` | refused |
| 2 | AMP-PD sample in a stratum with no step-6 fileset | refused, `MDE manifest 1 shipped 0 gap 1` |
| 3 | `AMPPD_CALLSETS` names a callset that does not exist | refused, **and prints the observed `source_callset` values** — that error is a name typo, not an empty cohort, and the distinction is the whole message |
| 4 | bucket not private (UBLA off, PAP inherited, `allUsers` bound) | refused on all three |
| 5 | private bucket, clean upload | preflight OK, 8 objects |
| 6 | verify with matching CRC32C | success path |
| 7 | one object never arrives | `MISSING in bucket`, object-count mismatch, refused |
| 8 | push over an existing release, no `OVERWRITE` | refused |
| 9 | `MODE=push` inside a SLURM allocation | refused, points at helix |

Selection, the `.fam` join, the accounting arithmetic, the manifest and the generated README all
work. The CRC32C hash parse was exercised against a **real** `gcloud storage hash`, and the four
bucket-JSON shapes and four listing-JSON shapes against the fixed parsers directly.

**Next — two things a laptop cannot settle, both in `HANDOFF.md` issue 10.**

1. The gcloud module name on biowulf. `MOD_GCLOUD` defaults to `google-cloud-sdk`; confirm with
   `module spider google-cloud-sdk` before the first build.
2. Whether every AMP-PD sample reconciles. AMP-PD donors sitting in a stratum step 6 wrote no
   fileset for would trip the accounting guard. Nothing local can say whether any exist — the
   laptop's `analysis_grain.csv` is the stale 11,918-row one. The first real build answers it in
   one line.

**And when it is pushed: record the destination, object count, byte total, upload window and git
commit here.** The previous release's exact `cp` command and host were never written down and had
to be recovered by audit five weeks later (2026-09-15). The push stage prints exactly this list on
success.

---

## 2026-09-15 — the sumstats release destination, recovered and now manifested: `gs://sysbio-gwas/results`

**Did.** Recovered the release destination, which nothing in the repo had recorded, and wrote
`SUMSTATS_MANIFEST.tsv` — 46 objects with size, CRC32C, MD5 where available, component count,
generation and upload timestamp. `SUMSTATS_README.md` §Access now names the bucket.

**Found — the upload is intact and matches the release doc exactly.**

| | |
|---|---|
| destination | `gs://sysbio-gwas/results`, US-CENTRAL1 |
| uploaded | **2026-08-25 02:51:40–02:52:09 UTC** — a 29-second window, so one `cp` invocation |
| contents | 44 `.tsv` + `gwas_summary.csv` + `README.md` = 46 objects, 22.57 GiB |
| doc vs bucket | the 44 contrasts in `SUMSTATS_README.md` and the 44 in the bucket are the **same 44** — no file in one and not the other |
| access | uniform bucket-level access **on**, public access prevention **enforced**, no `allUsers` or `allAuthenticatedUsers` binding. Not public |

**The bucket's `README.md` is byte-identical to the committed `SUMSTATS_README.md`** — 18,401 bytes,
MD5 `5zYm+34zEwvN9kUaaXWGDw==`, matching `git show HEAD:SUMSTATS_README.md` exactly. That is the
provenance link the push never wrote down: the release doc served to the browser is this repo's, at
the state it was committed in. (Today's edit naming the bucket makes the local copy newer; re-upload
it or accept the one-line drift.)

**42 of 46 objects carry no MD5** — they were composite uploads (`Component-Count` 14–15, parallel
chunked `cp`), and GCS does not compute MD5 for composite objects. CRC32C is present on all 46 and
is what any integrity check has to use. Not a defect, but it means `gsutil cp -c` style MD5
comparison will not work against these.

**FOUND — the three non-converged AJ contrasts are visible in the file sizes, and they shipped.**
Median `.tsv` is 531 MB. Two files are wild outliers:

- `gwas_AJ_PD_vs_AD.filtered.tsv` — **10.9 MB, 2.1% of median**
- `gwas_AJ_AD_vs_DLB.filtered.tsv` — **104 MB, 19.6% of median**
- `gwas_AJ_control_amppd_vs_control_ampad.filtered.tsv` — 417 MB, the smallest full-size AJ file

These are exactly the three λ_GC = 0.0000 contrasts from the 2026-08-24 entry. The mechanism is the
release filter: rows are `ADD` only with |BETA| ≤ 5, and a Firth fit that returned
`ERRCODE=UNFINISHED` mostly fails that bound, so the non-convergence shows up as **missing rows**
rather than as a warning. The size ratio is an independent confirmation of the diagnosis, arrived at
from the bucket rather than from the logs.

**Consequence for the browser:** all three are in `gs://sysbio-gwas/results` and will render as
selectable results. They are below the 100-per-arm floor and their λ is meaningless. They should be
suppressed or labelled at the UI layer — the files themselves are correctly filtered, so nothing
downstream flags them.

**Still not recorded:** the exact `gcloud storage cp` command and the host it ran from. The 29-second
window over 22.57 GiB implies a high-bandwidth source, i.e. biowulf or a GCP VM rather than the
laptop, but that is inference and not a record.

---

## 2026-08-24 — the write-up numbers are now derived, not transcribed. `review/methods_numbers.py`. Three METHODS claims were wrong.

**Did.** Worked `METHODS_FLAGS.md` — the 14 items raised while building the slide deck. (That file
was a fourth doc competing with the three-doc split, and two of its items already existed as
`HANDOFF.md` issues. Folded into `HANDOFF.md` issues 6, 8 and 9 and **deleted the same day**; its
ground rules — do not invent a number, do not edit pipeline scripts, record fresh measurements in
this log with their job ID — were followed and are worth keeping in mind.) Wrote
`review/methods_numbers.py`, a read-only checker that re-derives every numeric claim in
`METHODS.md` from the artifacts and diffs it against the value the doc asserts. On a laptop it reports **zero
mismatches and seven checks it could not run**, all of which need cluster artifacts.
Negative-tested by corrupting an asserted value and confirming it reports `FAIL`.

**This is the real fix, and the flags were the symptom.** Every number in `METHODS.md` had been
typed by hand from one of six artifacts. Nothing tied the doc to the files, so drift was invisible
until someone summed a table. There is now one command that catches it:
`python3 review/methods_numbers.py --strict`.

**FOUND — three claims did not survive the derivation.**

1. **§10's "roughly six times any genuinely cross-program contrast" is false.** EUR `PD_vs_AD` —
   `cross_cohort`, and the primary contrast — removed **387,207** variants to differential
   missingness, *more* than `PD_vs_DLB`'s 353,068, and had been silently omitted from the
   comparison set. The mechanism argument is untouched; only the multiplier failed. Rewritten as
   the sharper true statement: `PD_vs_DLB` is the largest removal of any *within-program* contrast
   and sits at the same magnitude as the primary cross-program one **despite scoring Δ = 0.0**.
2. **§6.3's HWE denominator was wrong.** "~1,540 of ~9,800 EUR samples" — 1,540 is exact
   (`wgs_harm` membership in `results/retained_samples_manifest.csv`, counted by splitting the
   pipe-joined `source_callset`), but EUR is **10,135**, as §5 and §6.4 both already said. 9,800
   appears nowhere else in the repo.
3. **§7's λ_GC 1.0175–1.0549 was quoted immediately after "44 contrasts ran"** and reads as
   covering all 44. It covers the 17 viable. Now scoped, with the non-viable spread stated.

**FOUND — why three AJ contrasts return λ_GC exactly 0.0000. Not previously logged.**
`AJ AD_vs_DLB`, `AJ PD_vs_AD` and `AJ control_amppd_vs_control_ampad` are AJ's three *cross-callset*
contrasts. Every variant in all three carries `ERRCODE=UNFINISHED` (100% of the first 300k rows in
two of them; the third is `UNFINISHED` + `INVALID_RESULT`), with *P* ≈ 0.99996 — so the median
chi-square is ~0 and λ collapses to 0. **The logistic model did not converge; the result is not
deflated.** The cause is visible in §6.4: AJ's `AD` arm is 95/97 `wgs_harm` while the PD/DLB/AMP-PD
control arms are ~96% `wb_dwgs`, and AJ PC1 is η² 0.962 on callset and enters as a covariate — so a
covariate separates the arms almost perfectly. The contemporaneous AJ contrasts whose arms *share* a
callset (`PD_vs_control`, etc.) come back `ERRCODE=.` and λ ≈ 1.03. All three λ = 0 contrasts are
below the 100-per-arm floor and none is reported, so nothing downstream is affected. **This is the
association-test face of the AJ collinearity that §6.4 and §10 already describe in PC space.**

**RECONCILED — "0.058%" vs "0.05%", the same quantity in two docs, neither with a denominator.**
0.058% = 4,385 / 7,538,809, and both halves belong to the **ungated 4,415** list: 4,385 is how many
of it were present in EUR's post-QC set (the 2026-08-19 premise test, EUR 7,538,809 → 7,534,424).
The live list is 4,187. Every doc now states the denominator and says **0.06% of EUR's 7,538,809
post-QC variants**; `methods_numbers.py` derives the exact gated figure from
`by_ancestry_qc/unfiltered/cohort_EUR_qc.bim` when run on the cluster.

**Also written into §6.4, all from `af_filter_effect_gated4187{,_per_pc}.csv`, none of it new data.**
The before/after cells are different components (PC2 → PC6) — worst-axis-before vs worst-axis-after,
which is the right comparison but was left for the reader to notice. The table is a **maximum over
PC1–PC10** and is upward-biased as such; now labelled. **The EUR result was understated:** before,
four PCs carried structure (PC2 0.757, PC3 0.155, PC1 0.080, PC7 0.038); after, all ten sit at or
below 0.036 and PC2 itself is 0.011 — the profile collapsed, the peak did not move. And η² is not
interpretable to two decimals at *n* = 106–226: under a filter that only *removes* variants, CAH PC6
rose 0.098 → 0.123 and AFR PC6 0.005 → 0.050. Noise, and now said to be.

**Fixed stale cross-references** left by the 617 → 210 line `HANDOFF.md` cut: `README.md` cited
"known issue 4" for the cluster-only clinical scripts (issue 1), and this file's top block cited
issues 9 and 10 for the age covariate and `confound_tag` (5 and 6).

**Next — the 7 checks that cannot run off a laptop.** Each is one read-only command on biowulf and
`methods_numbers.py` reports the missing path for it. Two are load-bearing for the paper: the
**94-genome gap** between §1's table (13,428) and §3's merge (13,334), which the `.fam` files named
in `merge_list.txt` settle exactly; and **§6.4's "~7× enriched for duplicate-pair discordance"**,
which is the only non-circular evidence that the 4,187 are technical and has no numerator,
denominator or job ID anywhere. That one gets measured or deleted — it is not staying as an
unsourced multiplier. `HANDOFF.md` issue 9.

---

## 2026-08-21 — git remote created. Subject-level data found in notebook outputs and scrubbed from history. Four retired scripts deleted.

**The repo has a remote for the first time:** `SysBio-FAIRplex/amp-ad-pd-wgs-gwas` (private). Until
today the entire history was 12 commits on one laptop disk.

**Found before the first push: `demo_sample_check/scripts/demo_sample_check.ipynb` carried 18 stored
cell outputs, one of them subject-level** — a rendered DataFrame of `individual_id` × callset
membership (real ROSMAP donor IDs, enumerated in the output). The other 17 were aggregate. This is
controlled-access data under the AMP-AD/AMP-PD DUA, and `.gitignore` already excluded the *file*
form of exactly this content (`demo_sample_check/out/` — "carries donor IDs joined to phenotype").
It leaked because `.gitignore` governs paths and this was inside a tracked file.

**Fixed by history rewrite, not by a follow-up commit.** The notebook entered at `1277fc4` and was
never modified after, so `git filter-repo --path … --invert-paths` cost no real source history. The
cleaned notebook was re-added as a fresh commit (`046e723`) and the remote created only afterwards,
so the blob was never published. Verified: `git log --all -- <path>` and
`git rev-list --objects --all | grep` both empty before re-adding.

**Also checked and clean:** `wgs_core.ipynb` has zero stored outputs; both tracked
`results/pca/af_filter_effect_gated4187*.csv` are aggregate (one row per stratum, no IIDs).

**Guard added — `scripts/nb_guard.py` + a `pre-commit` hook** (source at `scripts/hooks-pre-commit`,
installed by hand after a clone). It refuses any staged `.ipynb` carrying stored outputs. It does
**not** try to judge which outputs are aggregate: a judgement is precisely what failed here, and
17-of-18-safe is what makes that judgement look reliable. Verified by negative test — a planted
notebook with one output was refused, exit 1. Rule 3: the check proves it can fire.

**Deleted, with reasons, so the grep finds them here rather than the files:**

| file | why it went |
|---|---|
| `scripts/diag_threads.py` | Blocker 1 (GenoTools worker-pool sizing) closed 2026-08-12; the finding — "the worker *count* is the lever, not the thread env vars" — is in README and at 2026-08-12 below |
| `scripts/diag_cah.sh` | Blocker 2 (all-97-CAH) closed 2026-08-14 by `--sort-vars` in step 1; postmortem of a genotools output dir with no remaining subject |
| `review/compare_pcs.py` | forensic tool for "did the PCs change between two manifests", written during the stale-exclusion-list investigation (2026-08-17). Handled sign flips and PC reordering. No callers; `plot_af_filter_effect.py`'s eta² tables answer the live version of the question |
| `review/mask_cohort_artifacts.py` | the last known duplicate implementation — see the entry above. Its docstring cited `gwas_per_ancestry.sh` (retired) and `docs/METHODS.md` (never existed) |

Deleted rather than moved to an in-repo `archive/`: a second copy of a dead concept is the rule-2
failure mode, and with a remote now in place git history is a real archive. The cluster's
`scripts_archive/ad-pd-gwas-jul28/` is untouched — that is July provenance, not a graveyard for
these.

**Rule 2's other half is NOT done for these four.** `rsync` without `--delete` cannot express a
deletion, so all four still exist on biowulf until removed by hand. That `rm` is pending — see
HANDOFF.

---

## 2026-08-21 — step 8 RAN (job 28037485). APOE and HLA-DQB1 flagged: the annotate-never-subtract licence is CONFIRMED, not just argued.

**Did.** Step 8 completed, exit 0. `.ccannot.tsv` + `.ccfilt.tsv` written per contrast; primary
sumstats untouched. First run with the sentinel block removed — gene naming came from the
`gene_annot.py --at` CLI instead, which is the replacement path working as intended.

**FOUND — the docstring's central claim is now measured.** `08_ctrl_ctrl_filter.py` has argued since
it was written that "APOE is expected to reach significance in the control-vs-control scan for an
entirely real reason" and that subtracting "would delete the strongest true locus in the study."
Both are now confirmed. All **4** flagged genome-wide hits in `EUR AD_ampad_vs_control_ampad` — the
*within-cohort* AD contrast, Δ_amppd = 0.0, the cleanest AD comparison in the study — are APOE:

| variant | identity | P (AD vs ctrl) | CTRL_P |
|---|---|---|---|
| chr19:44908684:T:C | **rs429358, APOE-ε4** | 3.55e-15 | 1.29e-07 |
| chr19:44912456:G:A | ~3 kb past APOE's 3′ end (toward APOC1) | 9.17e-12 | 4.28e-06 |
| chr19:44906745:G:A | inside APOE | 2.27e-11 | 6.20e-06 |
| chr19:44908822:C:T | **rs7412, APOE-ε2** | 6.68e-09 | 2.91e-07 |

APOE's refFlat extent is 19:44,905,796–44,909,395, so three are inside the gene. **`.ccfilt.tsv`
would delete APOE-ε4 at P=3.55e-15 from the cleanest AD contrast.** Mechanism is the documented
screening asymmetry: AMP-AD controls were assessed cognitively normal, depleting ε4 carriers who
would have converted, while AMP-PD controls were screened for PD only and carry ε4 at population
frequency. **Do not use `.ccfilt.tsv` for the AD contrasts.** Report the primary and name APOE as
expected-and-retained.

**FOUND — 2 of the 3 genome-wide ctrl-vs-ctrl hits are `HLA-DQB1`** (chr6:32,661,554 and 32,661,570),
MHC class II, ~30 kb from `HLA-DRB1/DRB5` — the study's other real finding, and an established PD
locus. Same mechanism: AMP-PD controls screened against PD are depleted of PD risk alleles. The
third hit, `chr3:106666502`, is **intergenic, 443 kb from the nearest gene** (`LINC00882`), with no
AD/PD story — that one is the plausible genuine residual batch artifact. So of 3 genome-wide
ctrl-vs-ctrl hits: **2 real biology, 1 likely technical.**

**FOUND — the high-signal contrasts are untouched, which is the corroboration that 6a worked.**
`AD_vs_PSP` 2,319 hits / **0 flagged**; `PD_vs_control` 2,523 / **0**;
`PD_amppd_vs_control_amppd` 2,508 / **0**; `PD_vs_DLB` 318 / 3. Thousands of true-signal hits with
zero control-flagging means step 8 is confirming rather than correcting. Flagging concentrates in
the cross-cohort AD comparisons (`PD_vs_AD` 49/10, `AD_vs_control` 34/6, `AD_vs_DLB` 33/6) and, as
above, lands on real biology even there.

**Caveat to state, not bury: `PD_vs_MCI` is 3 hits, 3 flagged, 0 kept** — its `.ccfilt.tsv` is
empty. With n_ctrl = 159 that contrast was thin regardless, but "every hit flagged" must be reported
as such rather than presented as a null result.

**Open, and worth a methods sentence: `--ctrl-p 1e-5` is what catches APOE.** All four APOE CTRL_P
values are 1.3e-07 – 6.2e-06, i.e. **above** 5e-8 — APOE was not among the 3 genome-wide
ctrl-vs-ctrl hits. The demo's default 1e-5 threshold is aggressive enough to sweep up a locus that
is depleted by study design. At `--ctrl-p 5e-8` the AD contrasts would lose nothing. Not changed:
the flagging is annotation, and a wider net with correct labels is the safer default. But the
threshold, not the mechanism, is what put APOE in the flagged set.

---

## 2026-08-21 — step 7 RAN clean (job 28004190). Found: `within_cohort` hides callset asymmetry.

**Did.** Step 7 completed, 44/44 contrasts `ran=yes`, first run of the rewritten script that reads
§13's pheno/covar. Both new staleness guards passed. λ_GC spans **1.0175–1.0549** across every
contrast — tight control for a design this confounded.

**Found — two results worth carrying.** `control_amppd_vs_control_ampad`, the batch-artifact scan
with no true signal by construction, produced **3** genome-wide hits at λ=1.0319: residual cohort
artifact is close to gone, which is 6a's filter working. And `AD_vs_PSP` produced **2,319** hits
against the **2,318** recorded in `gene_annot.py`'s old MAPT±500kb flank measurement — independent
confirmation the pipeline finds the 17q21.31 block it should. Highest λ is `PD_vs_AD` (the primary,
`cross_cohort`) at 1.0549; lowest is `AD_ampad_vs_control_ampad` at 1.0175.

**FOUND — the confound tag is blind to the callset split inside AMP-PD, and it matters.**
`AMPPD_CALLSETS = {"wb_dwgs", "br_dsnwgs"}`, so `delta_amppd` pools them. `EUR PD_vs_DLB` scored
Δ=0.0 → `within_cohort` — the label `07_gwas.sh` calls "the confound-free trusted backbone" — while
the diffmiss filter removed **353,068** variants from it, ~6× any genuinely cross-program contrast
(~60k). Mechanism, from the grain rather than assumed: BR's 95 retained are 71 PD / 21 control /
3 null, so **zero DLB**, and BR is ~50% missing on the common set. Control:
`PD_amppd_vs_control_amppd` has BR on both arms (71 v 21) → 5,814 removed. ~55 EUR genomes produced
more technical asymmetry than the whole AMP-AD/AMP-PD split.

**Not a correctness problem — a labelling one.** Diffmiss is pre-association and caught it (λ=1.0402,
318 hits). But the summary ranks `PD_vs_DLB` as cleaner than `AD_vs_DLB`, which is backwards.

**Changed.** §13 emits `max_callset_delta`, `worst_callset`, `callset_one_sided`; step 7 carries them
into `gwas_summary.csv`, annotates `within_cohort` rows, and ends by ranking `within_cohort` rows by
`n_diffmiss_excluded`. `review/plot_gwas.py` shows the one-sided callset in the figure subtitle
(via `.get`, so pre-2026-08-21 summaries still plot). `confound_tag` is **unchanged** — it means
"program", and redefining it would reinterpret every row on record. Step 7's summary schema grew by
three columns, so the final viable-rows block was converted from positional `$11/$13/$14` to
header-name indexing; it had silently pointed at the wrong columns the moment they were inserted,
which is rule 4 applying to a file the script wrote itself.

**RULED OUT — a percentage-point threshold on callset share.** The first implementation returned only
`max_callset_delta` and reused the confound tag's 20pp boundary, on the reasoning that it introduced
no new magic number. It is wrong: BR is ~55 of a 2,595-sample arm, so its delta is **2.1pp** and a
20pp rule misses the single case the function was written to catch. Caught by unit-testing the
function against the real arm compositions before it shipped, not by reasoning about it.

**RULED OUT — one-sidedness as a warning.** The replacement flagged any callset with ≥10 samples in
one arm and 0 in the other. It fires on **5 of 5** real contrasts, which is exactly the objection
that retired the sentinel tripwire the day before: a flag that fires on everything carries no
information. The real predictor is one-sidedness of a **sparse** callset — `PD_vs_control` is
one-sided (`wgs_harm` 0 v 407) and removed only 6,331, whereas BR one-sidedness removed 353,068,
because a 97-donor joint call emits nothing at sites monomorphic in its own donors. §13 holds no
missingness data and cannot compute sparsity. **Resolution: `n_diffmiss_excluded` is already the
measurement, so the columns are descriptive and the ranking replaces the threshold.** Two orders of
magnitude separate the BR rows from the rest; no cutoff needed.

**Next.** Step 8 (`08_ctrl_ctrl_filter.sh`) — note it will find few flags, since the ctrl-vs-ctrl
scan yielded 3 hits. Then `review/plot_gwas.py`. §13's new columns mean **`analysis_grain.py` must be
rerun before step 7 is ever rerun** — an existing `contrasts.csv` lacks them and step 7 exits 3
naming the missing column.

---

## 2026-08-20 (night, last) — eta² re-measured on the gated list. The HWE gate cost the filter nothing.

**Did.** Pulled both step-6 manifests to the laptop and ran `review/plot_af_filter_effect.py
--label gated4187`. Wrote `results/pca/af_filter_effect_gated4187.{png,csv}` and
`..._per_pc.csv` — **the first eta² tables ever actually committed**, see below.

**Found — the prediction held, and the headline claim is safe on the live list.**

| stratum | n | before | after | reduction | on the ungated 4,415 |
|---|---|---|---|---|---|
| EUR | 10,135 | 0.757 (PC2) | **0.036** (PC6) | **95.3%** | 0.036 / 95.2% |
| AJ | 1,518 | 0.984 (PC1) | 0.962 (PC1) | 2.2% | identical |
| AAC | 226 | 0.877 | 0.876 | 0.1% | — |
| AFR | 191 | 0.708 | 0.713 | **−0.6%** | — |
| CAH | 106 | 0.309 | 0.199 | 35.6% | — |
| AMR | 185 | 0.203 | 0.170 | 16.0% | — |

Removing the 228 HWE-only variants moved EUR's reduction by 0.001. That is the direct test of the
"HWE-only additions are general variant QC, not cohort-artifact removal" reading, and it passed:
the frequency channel is carrying the effect, and the gate touched only the channel that was not.
**Quote 95.3% against the 4,187 list; the 95.2% figure was the 4,415 list and both are now on
record.** AJ is unchanged to three decimals, as expected — the list holds no AJ-derived flags.

**Ruled out — AFR's negative reduction is not a regression.** AAC and AFR eta² is one BR sample each
at 39.7σ / 19.6σ (2026-08-18 entry), so those rows are dominated by a single point and move by noise
in either direction. CAH and AMR are new numbers, not previously in the two-row table; none of the
four strata fields a `viable_ge100` contrast, so none of these PCs reaches a GWAS.

**Two operational notes worth keeping.** The laptop's system python3 has pandas/numpy but **no
matplotlib** — this script needs `.venv/bin/python`, and the failure is an import traceback at line
40, before any output. And both manifests share the basename `retained_samples_manifest.csv`, so
they must land in separate directories; pulling them into one would compare a file against itself
and report a ~0 reduction that looks exactly like "the filter does nothing".

**Closed the gap between `.gitignore` and reality.** The negations `!results/pca/*_eta2.csv` and
`!results/pca/af_filter_effect*.csv` were added 2026-08-19 with a written rationale — a *baseline*
cannot be regenerated later because that needs a state that no longer exists — but no file had ever
been produced under them, so `results/` was entirely untracked and the 0.757/0.036 numbers survived
only as hand-typed text in two markdown files. That is precisely the fragility the negation was
meant to fix, left in place for a day. These three files are the first to land there.

---

## 2026-08-20 (night, later still) — pheno/covar duplication RESOLVED. Step 7 reads §13's files.

**Changed.** `07_gwas.sh` now reads `clinical_core_out/{pheno,covar}/` and `contrasts.csv` instead of
rebuilding both in awk from `$GRAIN`; the awk blocks (covar builder, pheno builder + arm counting,
FID map, confound-tag arithmetic, age-column sniffing) are deleted. `config.sh` gained
`PHENO_SRC` / `COVAR_SRC` / `CONTRASTS_CSV` and lost the NOTE that documented the duplication.
Known issue 1 is closed; `CLAUDE.md`'s rule-2 corollary now lists one live duplication, not two.

**Why this one and not the other direction.** §13's docstring already asserted "Step 7 reads these
files and runs `plink2 --glm`; it does not parse the grain and does not rebuild an arm in awk, so the
definition exists once and cannot drift between the two." That was false, and a false claim in the
authoritative implementation is worse than no claim. Making it true also puts the definition of "who
is a case" next to the reconciliation logic that produces it, and makes the pheno files inspectable
before a 3-hour job starts rather than materialising inside it.

**Verified before committing to it: the formats were already identical** — covar `#FID IID SEX [AGE]
PC1..PC10` tab-delimited with `NA` for missing, pheno `#FID IID pheno` with case=2/ctrl=1, same
`<case>_vs_<ctrl>` tag convention with `@`→`_`. And §13's `FAM_DIR = PCA_DIR = MERGED/by_ancestry_qc`
is the same directory step 7 takes `--bfile` from, so the FIDs come from the same `.fam`. This was a
wiring change, not a reimplementation.

**Three behaviour changes, all deliberate, all recorded in HANDOFF issue 1.** `MIN_ARM` can now only
*raise* the 20-per-arm floor (§13 wrote no file below it), so `MIN_ARM=0` no longer forces every
contrast. `EXCLUDE_DUAL` is gone from step 7 and lives only in §13 — which is what §13 argued for,
and step 7 now **refuses** the old `--export=EXCLUDE_DUAL=1` rather than silently ignoring it, since
a sensitivity run returning as the primary is the worse failure. The age decision is §13's.

**Found — the old "PREREQ" was a comment where a check belonged.** The header said re-running step 6
invalidates the grain and "this silently runs on stale covariates", which is rule 3 exactly: the
warning was unactionable and nothing enforced it. Replaced with two guards — `covar_<ANC>.txt` must
post-date `cohort_<ANC>_pca.eigenvec` (exit 4), and every pheno IID must appear in the `.fam` being
tested (exit 5). The second is not redundant: a covar can be newer than the eigenvec and still list
samples the current fileset does not contain, which mtime cannot see.

**`contrasts.csv` is read by header name** (rule 4), emitting a fixed-order tab stream for the shell
loop. Tested against a deliberately column-shuffled CSV — parsed correctly — and against a renamed
column, which aborts naming the missing field rather than mislabelling a contrast. Arms containing
`@` survive the round-trip. Not yet run on the cluster.

---

## 2026-08-20 (night, later) — HWE gate RUN. List 4,415 → 4,187. Both predictions held.

**Did.** rsynced the sentinel removal, reran step 6 with `HWE_REQUIRE_EXCESS` at its default (1),
reran `analysis_grain.py`. Full job, not just stage B — stage D's PCs and stage E's manifest both
regenerated (21:37).

**Found — the gate behaves exactly as designed, and the arithmetic closes.** The list is **4,187**,
down 228 from 4,415. `AJ`/`wb_dwgs` (0.24×) and `EUR`/`wb_dwgs` (0.35×) both print `WITHHELD`;
`EUR`/`wgs_harm` (4.5×) voted. Provenance records `require-excess : True`. HANDOFF predicted those
two cells "still contributed up to 229 of the 4,415" — the observed 228 means 228 were unique to
them and exactly 1 was independently flagged by another channel. Both settled predictions held:

| prediction | outcome |
|---|---|
| `chr12:40227079:C:T` (LRRK2) un-excluded | **confirmed** — absent from the list |
| CR1 `chr1:207521012:T:C` still excluded | **confirmed** — present (frequency test + 64.5% `divco_hs` call rate) |
| no sentinel output anywhere | **confirmed** — `step6_summary.txt` has none |

**Found — the grain is unchanged in shape on the new PCs:** 12,495 rows × 22 cols, 17 viable
contrasts, `within_cohort=17 cross_cohort=15 partial=12`, BR's 95 retained as 76 EUR / 13 AJ / 3 AMR
/ 1 AAC / 1 AFR / 1 CAH. So changing the list by 228 variants moved no contrast across the ≥100
viability floor.

**Ruled out, by grepping the log first rather than investigating (rule 1 working as intended).**
`analysis_grain.py` printed `214 dx conflicts` and `10 pheno conflicts`; those exact counts are
already recorded at the 2026-08-17 entry as expected and reconciled by §12. No diagnostic was run.

**Found — the age covariate is absent, by design, and that design was documented only in code.**
`analysis_grain.py` printed `age covariate: NOT FOUND in grain — running WITHOUT age`. Not a
regression: `07_gwas.sh:52-57` and `analysis_grain.py:240-241` both record the reason — AMP-AD
supplies age at death, AMP-PD age at baseline/analysis, and they are not the same variable, so it
cannot be forced. Both detect an age column by header name (`age`, `age_analysis`, `agedeath`,
`age_death`, `age_baseline`, `age_cov`) and warn when absent. **Age is the dominant confounder for
both AD and PD, so this is a live methods limitation** and it was reachable only by reading two
scripts — now HANDOFF known issue 9. Note it is also a third face of known issue 1: the same
six-name detection is implemented twice, independently.

**Still open — the eta² table is now measured on a superseded list.** The
`EUR 0.757 → 0.036 (95.2%)` figure is the 4,415 measurement and the live list is 4,187. The 228
withheld are HWE-only additions (general variant QC, not cohort-artifact removal), so EUR should
barely move — but that is a prediction, and this table is the project's first-party proof the filter
works, cited in three places. `review/plot_af_filter_effect.py` against this run's two manifests is
the next action; HANDOFF carries a warning on the table until then.

---

## 2026-08-20 (night) — the sentinel is DELETED from both call sites, and the 08 sub-question is decided

**Did.** Executed the removal decided earlier today. From `scripts/af_concordance_build.py`: the
`gene_annot` import, `sentinel_detail()`, `keep_n()`, the `---- sentinel loci ----` block including
its refFlat-failure banner, `arm_index` and the loop populating it, and the per-arm `.keep` writes in
`assoc_pair()`. From `scripts/08_ctrl_ctrl_filter.py`: the import and the KNOWN LOCI block. From
`scripts/gene_annot.py`: `SENTINEL_GENES`, `SENTINEL_FLANK`, `SENTINEL_VARIANTS`, `sentinel_hits()`.
All three compile; the CLI is smoke-tested in both modes. **Not yet run on the cluster, and not yet
rsynced.**

**Decided — the 08 sub-question: remove both.** HANDOFF framed this as genuinely open, because 08's
report served the opposite purpose (a known locus in the ctrl-vs-ctrl scan is the screening
asymmetry showing up, i.e. an argument *against* subtracting) and step 8 never deletes. The
volume critique also does not transfer: 08's flagged set is variants at ctrl-vs-ctrl P < 1e-5, tens
to low hundreds, not 4,415, so it would not have produced a 1,700-line census.

What settled it is the **negative**, not the volume. Line 110 of the old file could print `no known
AD/PD locus among the flagged` off an uncited 9-gene list, and in this file that sentence is
actively dangerous — it reads as reassurance that the flags are safe to subtract into `.ccfilt.tsv`.
A flagged variant sitting on a real locus outside those 9 produced byte-identical output. That is
rule 3 in the exact shape 08's own docstring was written to fix (it once carried a hardcoded window
table covering 37% of CR1 and naming no HLA gene, so it could report "no known locus flagged" while
sitting on HLA-DRB1); resolving coordinates from refFlat narrowed the hole without closing it. The
±500 kb mislabelling applies identically — the 21 "LRRK2±500kb" variants were in SLC2A13/C12orf40.

**Nothing was lost, and this is why the removal is not a downgrade.** The substantive warning 08's
block existed to deliver does not depend on a gene list, so it is now printed unconditionally, for
every flagged variant rather than for nine: the control arms differ by disease *screening*, a real
AD locus is EXPECTED here, judge hits individually in `.ccannot.tsv`. The per-locus naming that was
genuinely useful moved to `gene_annot.py --at <chr:pos>`, which is unarbitrary and correctly
labelled.

**Changed, incidentally, both to stop deleted concepts from lingering (rule 2).**
`sentinel_regions` → `gene_regions`: it was never sentinel-specific — it resolves whatever symbols
it is handed — and the name would have outlived the set it was named for. And the `--at` CLI now
parses through `parse_variant_id` rather than its own inline `replace("chr","")`/`partition`, which
(a) stops `parse_variant_id` becoming dead code now that `sentinel_hits` is gone, (b) accepts full
`chr12:40227079:C:T` variant IDs so IDs from an exclusion list paste straight in, and (c) errors
cleanly on unparseable input instead of raising.

**Changed — `CLAUDE.md`, three places that asserted a check which no longer exists.** Rule 6 said
"the sentinel tripwire is the check on 6a's licence"; that was the most load-bearing staleness, since
it told the next session to rely on something absent. It now says there is deliberately no automated
check and describes what replaces it (stage B's BY CALLSET PAIR and HWE ratio tables, the per-cell
`.afreq` intermediates read by header name, call rate as the usual tell). Rule 3's "both sentinel
call sites now print a banner instead" became the stronger general lesson: **a check whose negative
result is not evidence should be deleted, not annotated with a caveat.** Rule 2's example now
records that wiring the replacement is what made the concept testable, and the test killed it.

**Kept, deliberately, as comments where the code was.** Three things that cost a run each and would
otherwise be re-derived: that call rate rather than frequency is what resolves these hits (CR1:
156/242 alleles in `divco_hs`, 242/242 at the LRRK2 site — dropout, a mechanism); that reading
`.afreq` positionally hid that number on the first attempt; and the MAPT±500kb measurement (2,286 of
2,318 EUR AD-vs-PSP hits, the other 32 beyond 46,528,333 in the 17q21.31 inversion tail) which is
why 500 kb was a defensible proxy for step 8's question and a poor one for 6a's.

**Ruled out — a claim made and retracted within this entry: "`ref/` is not in git."** It is. Both
`ref/refFlat.txt` and `ref/highld_exclude_hg38.bed` are tracked, `.gitignore` exists, and its lines
61-64 record the deliberate decision to ship them ("small, public, and the reason a fresh clone can
run gene_annot.py with no setup", with the large non-redistributable references left in `data/ref/`).
So `README.md`'s "ships with the code" is accurate and nothing needs committing.

**The cause is worth more than the claim: `Bash` keeps its working directory between calls.** An
earlier `cd scripts` persisted, so `ls -a`, `git ls-files`, `grep .gitignore` and `git check-ignore`
all ran against `scripts/` — where there is indeed no `ref/` and no `.gitignore` — and every one
returned a true answer to the wrong question. Same failure class as rule 4 (`.afreq` column 5) and
the file-count check: the command succeeded, so nothing looked wrong. **Pass absolute paths, or
prefix `cd <root> &&`, when checking whether a file exists.**

**Next.** 1. Both rsyncs. 2. Step 6 with the HWE gate on (~7 min); read the withheld count and the
new list size. 3. `analysis_grain.py` (~1 min). 4. Step 7.

---

## 2026-08-20 (night) — the MHC flag-rate measurement is OFF the critical path (decision)

**Changed.** Dropped the MHC flag-rate measurement as a prerequisite for step 7. It was item 2 of 5
in the previous entry's "Next, in order" and is now a write-up task. `HANDOFF.md`'s status row and
its "one aggregate question worth keeping" section were reframed to match.

**Why — the decision it was meant to inform was already made, by the entry directly below this
one.** The measurement existed to answer "if the enrichment is large, should the MHC be annotated
rather than subtracted (the step 8 pattern)?" But that same entry's per-arm finding settles it
per-variant: the two natively-called callsets agree near 0 while the lifted one does not —
`chr6:32474706` `divco_hs` 0.008 / `wb_dwgs` 0.005 vs `wgs_harm` 0.207–0.213; `chr6:32588203`
0.000 / 0.000 vs 0.088–0.107. A both-native pair agreeing against the lifted callset is precisely
the licence condition for 6a to delete: disease is held constant inside the cell, so the gap is
technical. An enrichment ratio cannot overturn a per-variant mechanism, so the number could not
have changed the action — which is the definition of not a gate.

**What is NOT being claimed.** The MHC concentration is still real and still worth stating: 245 of
the 279 sentinel hits were MHC, this filter reaches the association set (unlike the high-LD BED,
which is PCA-input-only), and `HLA-DRB1/DRB5` is one of this study's two real findings. So the
magnitude belongs in the methods as a limitation. The command survives in `HANDOFF.md` for that
purpose. What is retired is the idea that it blocks anything.

**Next, in order — revised.** 1. Remove the sentinel from `af_concordance_build.py`, deciding the
`08_ctrl_ctrl_filter.py` sub-question first. 2. Run step 6 with the HWE gate on; read the withheld
count and the new list size. 3. `analysis_grain.py`. 4. Step 7. (The MHC rate moves to write-up.)

---

## 2026-08-20 (evening) — the --assoc swap is VERIFIED. Decision: delete the sentinel entirely.

**Did.** Ran the swapped frequency test against the pre-change list. Fixed a precision bug it
exposed. Ran the sentinel with corrected gene coordinates for the first time, read the result, and
decided to remove the sentinel from the pipeline. **No code was changed for that removal — it is
written up as the next action item in `HANDOFF.md` and is the first thing to do next.**

**THE SWAP IS VERIFIED, on every cell.** `HWE_REQUIRE_EXCESS=0` (gate off, so the only difference
from the committed behaviour is `plink --assoc` replacing the Python z-test) reproduced the
exclusion list **byte-identically: 4,415 variants, `diff` empty**. Because the gate was off, this
tests the swap across all three powered cells, not just the EUR/AD one checked by hand earlier.
`plink --assoc` is now the frequency test, and it demonstrably changed nothing.

**A precision bug, found by that verification and worth remembering.** The first attempt came back
4,414 — one variant short, `chr6:32555808:T:C`, in the MHC. Both implementations had computed
chi-square = **25.005**; `.assoc` prints CHISQ to four SIGNIFICANT figures, so it lands on disk as
`"25"`, and the reader's `chisq > 25` then evaluated `25.0 > 25` → false. Fixed by thresholding on
**P** instead (`P < erfc(zmin/sqrt2)`, the identical rule), because P is exponential and pins
chi-square to ~0.001 near the boundary where `"25"` pins it only to 0.005.

This is the same failure family as the `.afreq` column mislabel earlier the same day: `CLAUDE.md`
rule 4 says read columns by NAME, and the lesson extends to reading them at **sufficient
precision**. A rounded column is not a wrong column, and it is just as capable of flipping a
decision.

**Also fixed: plink and plink2 are ONE MODULE FAMILY here.** `module load plink/1.9.0-beta4.4`
UNLOADS `plink/6-alpha` and leaves a non-executable `plink2` earlier on PATH, so a bare `plink2`
dies with `PermissionError`. The first `--assoc` run failed that way, three stages in, after the
frequency stage had already succeeded. Both binaries are now resolved to absolute paths while each
module is loaded and passed as `--plink1` / `--plink2`; the Python fails fast if either is not
executable. **`06_ancestry_qc.sh` also reloads `MOD_PLINK2` after stage B** — stages C and D call
bare `plink2`, so leaving plink1.9 active there would have failed the exclusion apply and the PCA
minutes after stage B reported success. Verified that plink2 runs correctly by absolute path with
the plink1 module active, and separately that `module` does work in a non-interactive subshell (the
`module load X && exec plink...` shim approach is a viable fallback if a future build needs its own
environment). This bug was LATENT before the swap: the old `MISHAP != 0` conditional load had the
same defect and escaped notice only because `MISHAP` defaults to 0.

**FOUND — the sentinel tripwire does not work, and the decision is to remove it rather than tune
it.** With the hardcoded coordinate table replaced by refFlat extents ±500 kb, it produced **279
hits across 8 loci and ~1,700 lines of per-arm detail**. Reasons it goes:

1. **279 hits is a census, not a tripwire.** Nothing is actionable at that volume, and nothing is
   auto-whitelisted, so it has no mechanical effect either — it is a report that fires *after*
   stage C has already applied the list in the same job.
2. **The 9-gene list is arbitrary and uncited** (established earlier this session — `git log
   --follow` gives one commit, no reference anywhere in the repo). It scrutinises a hand-picked
   handful while ~4,100 other flagged variants get none. Either every deletion needs justification
   or none does.
3. **The ±500 kb flank actively mislabels.** It was sized for step 8's question, where LD blocks
   matter. For 6a it names flank hits after the gene: **none of the 21 "LRRK2±500kb" variants were
   in LRRK2** — `SLC2A13` and `C12orf40`, 240–435 kb away. SNCA's was intergenic, GBA1's in `DAP3`,
   APOE's in `ZNF285`/`ZNF229`. That is worse than the under-reporting it replaced, because it
   manufactures alarm and could hide a genuine in-gene hit inside a list of 21 flank hits.
4. **Its per-variant detail was right for 2 hits and wrong for 279.** My design error.

**Honest accounting of what it bought before being retired.** Two hits, one useful. CR1 turned out
to be correctly excluded (64.5% call rate in `divco_hs` — dropout), so that hit changed nothing.
LRRK2 led to the HWE excess-over-chance defect, which was real — but that defect was already
visible in the ratio table printed on every run (4.5× / 0.35× / 0.24×). The sentinel supplied the
motivation to look, not the evidence. A per-locus **rate** table was considered as a replacement and
also rejected: it would not have prompted the LRRK2 question either, and it is one more thing to
maintain. The aggregate question survives as a one-off measurement instead (below).

**Ruled out: making it configurable** (`SENTINEL_DETAIL=1`). A knob preserves the maintenance cost
and the misleading labels while guaranteeing nobody sets it.

**FOUND — 245 of the 279 hits were MHC**, and this is the one thing worth carrying forward.
`HLA-DRB1±500kb` 213 + `HLA-B±500kb` 32, roughly chr6:30.8–33.1 Mb. Expected in the sense that the
MHC is the hardest region in the genome to align and call — but `06_ancestry_qc.sh`'s own header
argues the MHC is a real AD locus and that masking it from association would delete signals we most
expect to see (which is why the high-LD BED is PCA-input-only), and `HLA-DRB1/DRB5` is one of this
study's two real findings. **This filter reaches the association set.** The flag-rate command with
the correct denominator is in `HANDOFF.md`; run it before step 7 and decide deliberately whether the
MHC should be annotated rather than subtracted there.

**FOUND — the per-arm tables settle the liftover verdict per-variant, which the pairwise rate table
never could.** Repeatedly, the two natively-called callsets agree near 0 while the lifted one does
not: `chr6:32474706` divco_hs 0.008 / wb_dwgs 0.005 vs wgs_harm 0.207–0.213; `chr6:32588203`
0.000 / 0.000 vs 0.088–0.107; `chr12:40509081` absent / 0.000 vs 0.131–0.199. The 2026-08-19 entry
correctly flagged that the BY CALLSET PAIR table could not decide the liftover question because the
both-native pair was never measured. This is per-variant, unambiguous, and points at `wgs_harm`.
Several arms also show `OBS_CT=0` — the variant is absent from `divco_hs` entirely, not discordant.

**State for the next session.** Two commits landed during this session — `16d2c9a` (step 6 as one
pass, the clinical split, the `gene_annot` sentinel wiring, `.gitignore`/README/config) and
`ad73810` (`CLAUDE.md`, `analysis_grain.py`, `clinical_common.py`). **Everything from the `--assoc`
swap onward is still UNCOMMITTED**: `scripts/af_concordance_build.py` (+321 lines — `assoc_pair`,
the HWE gate, `sentinel_detail`, the provenance sidecar, the `--plink1`/`--plink2` args),
`scripts/af_concordance_build.sh` and `scripts/06_ancestry_qc.sh` (module handling + knob
passthrough), plus `review/plot_af_filter_effect.py` which is still untracked. Worth committing the
verified swap before the sentinel removal, so the two changes are separable in history.

The HWE gate is written and defaults ON but has **not been run** — no run has yet produced the
gated list, so 4,415 is still the live number.
`${MERGED_DIR}/exclude_af_concordance.4415.bak` holds the pre-change list;
`exclude_test_assoc.txt` is the verified swap output (also 4,415). The association set and both
manifests are still job 27857727's — untouched, because every test ran stage B alone through
`af_concordance_build.sh`, which is precisely the job that wrapper exists for.

**Next, in order.** 1. Remove the sentinel (HANDOFF has the exact list of what to delete, plus the
open sub-question about `08_ctrl_ctrl_filter.py`, whose sentinel serves the opposite purpose since
step 8 never subtracts). 2. Run the MHC flag-rate measurement. 3. Run step 6 with the HWE gate on
and read the withheld count and the new list size. 4. `analysis_grain.py`. 5. Step 7.

---

## 2026-08-20 (later still) — the frequency test is plink's now, and HWE cells must show excess.

**Did.** Replaced the hand-rolled frequency test with `plink --assoc`, added an excess-over-chance
condition to the HWE channel, and added two diagnostics (per-arm call rate under the sentinel
report, and a provenance sidecar). Not yet run.

**The swap was verified BEFORE it was made, and it changes nothing.** On the largest cell (EUR/AD,
`divco_hs`=121 vs `wgs_harm`=657, 7,538,809 variants):

| | flagged |
|---|---|
| Python rule (`\|dAF\| > 0.05` AND `z > 5`) | 1,698 |
| `plink --assoc`, CHISQ > 25, same `\|dAF\|` floor | **1,698** |
| symmetric difference | **0 in both directions** |

`plink --assoc` at CHISQ > 25 alone gives 2,008 — the extra 310 are variants that clear
significance but fail the effect-size floor, predominantly low-frequency ones. At MAF ~20% in this
cell `z > 5` already requires |dAF| ~ 0.14 so the floor never binds; at MAF ~2% it is reachable at
~0.035, and the floor is what stops the list filling with low-frequency noise. Worth stating in the
methods that the floor is **absolute**, so a 0.035 gap at MAF 2% — a large *relative* discordance —
is kept anyway.

**Why swap at all, given it changes nothing.** Communication. The Python version was a two-sample
test of proportions with pooled variance — textbook, and algebraically the same statistic
(z² = χ²) — but a hand-rolled implementation has to be taken on trust, while `plink --assoc` is a
named 1-df allelic chi-square a collaborator recognises on sight. The verification above is what
made the swap safe to do rather than a change with an unattributable effect.

**Cost, and it contradicts a stated preference:** `--assoc` is plink1.9 only (plink2 dropped it for
`--glm`, which is logistic regression on dosage — asymptotically equivalent, NOT identical, so it
would have moved the flags and cost us the verification). `MOD_PLINK1` is therefore now an
**unconditional** load in both `06_ancestry_qc.sh` and `af_concordance_build.sh`, where it used to
load only for the optional mishap stage. `07_gwas.sh:49` records a deliberate preference not to
depend on plink1.9; step 7 still doesn't, but step 6 stage B now does.

**The HWE gate.** A cell contributes exclusions only if its rejection count exceeds the number
expected by chance at `--hwe`. No multiplier — the bar is "more than chance explains". E/O is the
Benjamini-Hochberg FDR estimate for that cell's rejections, so at or below 1.0× there is nothing to
attribute. `HWE_REQUIRE_EXCESS=0` restores the old unconditional union.

Considered and rejected as over-built for the problem: BH at q=0.05 on `plink2 --hardy` p-values
(needs a different code path, an unverified column assumption, and a regression test, and it moves
the *good* cell too), and moving the HWE channel to a post-hoc annotation on sumstats. The latter
is conceptually right — for the association set, striking a variant from the sumstats is identical
to never testing it, since `--glm` tests variants independently and 5e-8 is a fixed convention — and
it stays available. **Note the one part that is genuinely not post-hoc-able: the PCA input, because
PCs are covariates in every test.** That asymmetry is the reason the frequency channel must stay
pre-applied whatever happens to the HWE one.

**Expected:** list drops by at most 229; LRRK2's `chr12:40227079:C:T` comes off (AJ HWE was its only
channel); CR1 stays (frequency channel, dropout mechanism); EUR eta² stays ~0.036 since
`wgs_harm`'s 1,680 are retained.

**Next.** `HWE_REQUIRE_EXCESS=0` must reproduce 4,415 exactly — which doubles as the regression test
for the swap across ALL three cells, not just the one verified by hand. Then the default run.

---

## 2026-08-20 (later) — both sentinel verdicts settled. The 2026-08-19 HWE claim was wrong.

**Did.** Read the two sentinel variants out of job 27857727's `.afreq` (by header name this time)
and `.snplist` intermediates. This resolves the contradiction flagged in the entry below.

**CR1 `chr1:207521012:T:C` — the dropout reading is CONFIRMED, exactly.** `divco_hs` EUR/AD shows
**OBS_CT 156 against 242 possible (121 samples × 2) = 64.5% call rate**, matching the 2026-08-19
entry's number precisely. Every other arm is complete: `wgs_harm` AD 1314/1314, `wb_dwgs` control
6124/6128, `wgs_harm` control 656/656, both `other` arms 100%.

**And the decisive control: the SAME `divco_hs` samples are 242/242 — fully called — at the LRRK2
variant.** So the dropout is specific to this site, not to those samples. Non-random dropout plus an
ALT frequency inflated to 0.365 against 0.170–0.212 everywhere else is lost reference calls.
**CR1 stays excluded**, on the frequency test (|dAF| 0.196, z ≈ 7) plus this mechanism.

**But the HWE corroboration claimed on 2026-08-19 does not exist.** That entry says CR1 "fails HWE
in **all three** testable callsets"; it **PASSES all three** (`AJ/wb_dwgs`, `EUR/wb_dwgs`,
`EUR/wgs_harm`). It also says LRRK2 failed "both testable callsets"; LRRK2 fails **AJ/`wb_dwgs`
only**. Aggregate HWE counts are byte-identical across jobs 27697096 and 27857727 (97 / 132 /
1,680), so the sets never changed and the earlier reading was simply wrong. **Nothing about the
conclusions changes** — CR1 was always excluded on frequency, not HWE — but "het excess in three
independent callsets is paralog/CNV collapse" must not be repeated. Per the rules in `CLAUDE.md`
this is a new entry; the 2026-08-19 entry stands as written.

**LRRK2 `chr12:40227079:C:T` — confirmed a false positive, and the reason is now airtight.** Not
frequency-flagged (max spread 0.058−0.023 = 0.034, under the 0.05 threshold). Call rate is complete
in every arm. Its only basis is one HWE failure in AJ/`wb_dwgs` — the 0.24×-of-chance row. **The
same callset, tested in EUR with 3,064 controls instead of 638, PASSES.** A real het-excess
mechanism in `wb_dwgs` would appear more strongly in the 4.8×-larger sample, not vanish from it.
This is stronger evidence than the ratio argument alone.

**So the HWE excess-over-chance gate is fully determined**: it leaves CR1 excluded (frequency-
derived, untouched) and un-excludes the LRRK2 variant. No ambiguity left.

**Incidental:** `AJ_*.afreq` does not exist — the glob failed. Expected, and a third confirmation
that no AJ frequency comparison ever ran: no AJ cell had two callsets above `MIN_CELL`.

---

## 2026-08-20 — step 6 ran as ONE pass. Grain rebuilt. The sentinel table was the stale fix.

**Did.** Ran the restructured step 6 (job **27857727**, 7m13s) and rebuilt the grain. Wired both
sentinel call sites onto `gene_annot.py` and deleted the hardcoded tables. Added `CLAUDE.md`.

**The refactor is behavior-preserving, measured two ways.** Stage B rebuilt the exclusion list
through `sample_annot.csv` and got **4,415 variants** — identical to job 27697096's count from the
grain. §12a's self-check reports **12,495 shared IIDs, 0 callset mismatches, 0 dx mismatches, 0
absent**. The annot split changed nothing.

**Stage A reuse worked: 7m13s against ~2.5 h.** All 11 strata printed `reused`, so the migrated
pass-1 filesets were adopted and `cohort_merged` was never re-scanned. Less work than the old pass
2 alone, as predicted. Filtered set: EUR 7,538,809 → 7,534,424 variants; both manifests written at
18:02 with 12,495 rows and `with PCs: 12,361`.

**Grain: 12,495 rows × 22 cols, 17 viable contrasts, and the BR-DSNWGS item is CLOSED.** The 95
retained BR samples are 76 EUR / 13 AJ / 3 AMR / 1 AAC / 1 AFR / 1 CAH — the corrected labels
minus one EUR and one AJ duplicate. The 19 AFR / 67 EUR grain is gone.

**AJ is closed, and not by the eta² argument.** `AJ/AD` is `wgs_harm`=95 + `divco_hs`=2 = 97,
under the 100 floor, so AJ cannot field the primary AD-vs-PD contrast at all. Its three viable
contrasts (`PD_amppd_vs_control_amppd`, `PD_vs_DLB`, `PD_vs_control`) all come back
`within_cohort` — ~96% `wb_dwgs` on both arms. eta² measures callset↔phenotype collinearity, which
cannot bias a contrast whose arms share a callset. **So the unsolved 0.962 has no consumer.** This
retires the AJ question on power and design rather than on the untested sub-continental-structure
hypothesis, which remains untested and no longer blocks anything.

**FOUND — the sentinel table was wrong, and the fix for it had been written and never wired in.**
`gene_annot.py` exists precisely to replace the hardcoded `[(name, chrom, pos, window)]` list, and
its docstring already measured the damage. Re-measured against `ref/refFlat.txt`:

| sentinel | window covered of gene | missed |
|---|---|---|
| CR1 | **37.0%** | 91,765 bp |
| LRRK2 | **52.0%** | 69,284 bp |
| SNCA | 66.5% | 38,304 bp |
| BIN1 | 96.7% | 1,977 bp |

The set also contained **no HLA gene**, while `HLA-DRB1` is one of this study's two real findings.
So `08_ctrl_ctrl_filter.py` could print "no known AD/PD locus among the flagged" while sitting on
it. `SENTINEL_LOCI` and `parse_id` were duplicated verbatim across both scripts — three copies of
one concept, and the consolidated version was inert in the same directory.

**Changed.** Both call sites now `from gene_annot import sentinel_hits` (full refFlat transcript
extents ±500 kb, build-verified via APOE, warns on unresolvable symbols); both hardcoded tables and
both `parse_id` copies deleted; both print an unmissable banner if refFlat cannot be read, because
"none found" must not be printable when the check did not run. Verified: the new windows catch both
of this run's hits as `CR1±500kb` and `LRRK2±500kb`.

**The tripwire is a report, not a gate, and that is new.** In the two-pass shape a human read it
between submissions. Stage C now applies the list four minutes later in the same job, so the
warning arrives after the fact. Left as a report with the timing stated in its own output;
converting it to a hard gate is still open.

**METHODS POINT, unresolved — the HWE stage has no excess-over-chance gate.** `hwe_failed |= fail`
unions every failure in regardless of whether that (stratum × callset) row shows any excess, and
the union is applied to **all** strata. This run: EUR/`wgs_harm` 1,680 vs 377 expected = **4.5×**
(real, and from the smallest sample of the three — consistent with mismapping in the lifted
callset); EUR/`wb_dwgs` 0.35×; AJ/`wb_dwgs` 0.24×. The two below-chance rows still contributed up
to 229 variants. Underdispersion is expected at these control counts, which is a second reason a
below-chance row should not vote. Gating on observed > 2× expected would keep EUR/`wgs_harm` and
drop the other two.

**Also found — a documentation contradiction to settle before that gate is built.** The 2026-08-19
entry says CR1 fails HWE in all three testable callsets; this run's `.snplist` files say it passes
all three, and LRRK2 fails only in AJ. The aggregate counts are byte-identical across the two jobs
(97 / 132 / 1,680), so the underlying sets are the same and one reading is wrong. Which one decides
whether the HWE gate would un-exclude the LRRK2 variant. **Unresolved.**

**Process finding, and the reason `CLAUDE.md` now exists.** Both sentinel hits had already been
resolved in the 2026-08-19 entry, with a *better* CR1 explanation (64.5% call rate in `divco_hs` —
dropout, not a frequency difference) than the one re-derived here from z ≈ 7. The re-derivation
also mislabelled `.afreq` columns positionally, hiding the very call-rate number that settles it.
Two rules written down: grep this log before investigating an anomaly, and never let a replacement
sit unwired.

**Next.** Settle the HWE contradiction (one command against the `.snplist` files), then decide the
HWE excess gate. Both sentinel verdicts stand for now: CR1 excluded on dropout, LRRK2 excluded on a
chance-level hit and would return if the gate lands. Then step 7.

---

## 2026-08-19 — step 6 collapsed to ONE pass. The circular ordering was never real.

**Did.** Restructured step 6 from two submissions with a hand-run script and a mandatory `mv`
between them into a single job with five stages. Split `sample_annot.csv` out of the grain,
switched `af_concordance_build` onto it, added `review/plot_af_filter_effect.py` and a
before/after section to `wgs_core.ipynb`.

**Found — the cycle that forced two passes does not exist, and it is worth being precise about
why, because the wrong version was written into four files.** The claim was
`grain <- §12 <- manifest <- step 6`, so step 6 could not depend on the grain. Two independent
facts dissolve it:

1. **`af_concordance_build` never needed the grain.** It reads `IID -> (source_callset,
   dx_detailed)` and nothing else. Its own `read_grain` also parsed `ancestry` — and the value
   was never used: `cell[(g[2], g[0])]` keys on (dx, callset), and the stratum comes from which
   `cohort_<ANC>_qc` fileset a sample is in. Both fields it does use are pure clinical output,
   available from §7's crosswalk and §4's reconciliation before step 1 runs. They arrived via the
   grain only because the grain is the file that happens to carry dx and the PCs together. A
   file-layout accident was read as a data dependency.

2. **The QC pass never needed re-running.** `--geno`, `--maf` and `--hwe` are per-variant on a
   fixed sample set (no `--mind`), so they commute with `--exclude`. This is the same identity the
   premise test relied on to measure the filter's effect without a rerun — it was used as a
   measurement trick and not recognised as a statement about the pipeline's structure. Pass 2 was
   re-scanning 527 GiB to recompute numbers it already had.

**Also found:** `ancestry_qc_manifest.py`'s docstring says §12 reads its output. §12 does not — it
globs `cohort_*_pca.eigenvec` and merges onto step 5's `retained_manifest.csv` directly. Its only
consumer is `plot_pcs_by_callset.py`. So the manifest is not on the grain's path either, and the
middle link of the claimed cycle was wrong as well.

**Changed.**

| | |
|---|---|
| `clinical_core.py` | new §12a writes `sample_annot.csv`, unconditionally, before any genotype step. Self-checks against an existing grain: shared IIDs must agree on `source_callset` and `dx_detailed`. |
| `af_concordance_build.py` | `--annot` replaces `--grain` (kept as an alias). Reader keys on column NAMES, so it takes either file — which is what makes the two directly comparable. |
| `06_ancestry_qc.sh` | five stages: A unfiltered QC+PCA -> `unfiltered/`, B build the list from A, C apply, D filtered prune+PCA, E both manifests. Stage A is skipped when a valid fileset postdates the merge. |
| `af_concordance_build.sh` | demoted to a knob-tuning wrapper; step 6 calls the `.py` directly. |
| `review/plot_af_filter_effect.py` | new. Dumbbell of max eta^2 per stratum before/after, plus PC scatters for the two largest strata, shared limits per row. |
| `.gitignore` | `results/pca/*_eta2.csv` and `af_filter_effect*.csv` now versioned. |

**Two failure modes the old shape had, now gone by construction rather than by guard.** The
`mv by_ancestry_qc` before pass 2 — the baseline is a permanent named output. And the stale-list
refusal — stage B rebuilds the list in-job from this merge, so an old file cannot be applied. The
mtime guard survives only on `SKIP_AF_BUILD=1`, the one path by which a list this job did not
build can still reach stage C.

**Also fixed a bug introduced while writing this.** The first `.gitignore` draft un-ignored
`results/` wholesale to version the eta^2 tables, which swept in `results/retained_samples_manifest.csv`
— one row per sample. Narrowed to `results/pca/` with explicit negations and verified with
`git check-ignore`.

**Next.** Nothing has been run yet. From the current cluster state the migration is a rename:
job 27602590's output IS stage A's output, so notebook §3 moves it to `by_ancestry_qc/unfiltered/`
and step 6 skips the scan. What runs is stage B, the exclusion, and the filtered prune+PCA —
less work than the old pass 2 alone.

**The check that decides whether this refactor was behaviour-preserving:** stage B must produce
**4,415** variants, matching job 27697096, which built the list from the grain. §12a's self-check
is the finer-grained version of the same question. A mismatch means the two reconciliation paths
disagree and the exclusion list would move for a non-genotype reason — stop there, do not
proceed to step 7.

**Removed `wgs_core_draft.ipynb`.** 14 markdown cells, zero code — a runbook covering steps 0-8,
added in `1277fc4` after the executable `wgs_core.ipynb` already existed. Its status tables and
eta^2 prose were a near-verbatim fourth copy of `HANDOFF.md` and this file's header, and its
step-6 cell was titled "**runs twice**", so the restructuring above made it actively wrong. The
content that was genuinely unique — the step-0 bcftools/plink2 commands for BR-DSNWGS, and the
explicit `--export=` invocations for steps 1 and 2, which `HANDOFF.md` had been eliding as `...` —
is folded into HANDOFF's Run order. Earlier log entries still reference the file by name; they are
history and stay as written.

**Not addressed:** AJ. The filter still does not touch it (0.984 -> 0.962) and this change does
not pretend to. The reference-panel projection that would decide whether AJ's split is real
sub-continental structure is still unrun.

**Split the clinical side in two, same day.** `clinical_core.py` used to hold §1-14 and be run
twice; reaching §11-13 meant re-executing §1-10, which re-read eleven clinical files and
**rewrote the sex-update files step 1 had already consumed**. Nothing checked they still matched
what genotools applied — and given the laptop/cluster clinical-input divergence in known issue 4,
that is a live hazard, not a theoretical one.

Two *executions* are irreducible and this is worth stating plainly, because the shape looks
exactly like the step-6 cycle above and is not the same thing. Step 1 applies the sex files; the
grain carries step 6's PCs because step 7 reads them as covariates. There is no column to split
out — the round trip is real. What was avoidable is running the same script twice:

| file | runs | reads |
|---|---|---|
| `clinical_common.py` | imported, never run | — |
| `clinical_core.py` | once, before step 1 | clinical files + psams |
| `analysis_grain.py` | once, after step 6 | §9's audit tables + manifest + eigenvecs |

`analysis_grain.py` re-derives nothing: §9 already writes `individual_core.csv` and
`genome_crosswalk.csv`, and those carry every column §11-13 use. Verified before moving anything —
§13 turned out to need only the `grain` frame and the `.fam` files, and `AMPPD_CALLSETS` was the
one constant that had to move with it. **Section numbers are unchanged**, so every `§12`-style
cross-reference in `HANDOFF.md`, `config.sh` and `scripts/` still resolves; only the file moved.
§12a stays in `clinical_core.py` despite its number, because it must run before step 1.

**Also updated `README.md`** — it was dated 2026-08-10 and described a three-callset merge, sex
files living in each callset's `metadata/`, and a script list (`run_genotools.sh`,
`wgs_merge_s1_merge_bed.sh`, `update_sex.py`) that no longer exists. Now carries the four-callset
state, the `scripts/00`-`08` list, the real handoff locations, and the gotchas learned since
August. It does **not** carry a run order — that lives in `HANDOFF.md`, and two copies of a run
order is how the wrong one gets followed.

---

## 2026-08-19 — PREMISE TESTED. The filter works in EUR (95%) and does nothing in AJ (2%). The 0.748 was EUR's.

**Did.** Measured eta² before vs after applying the 4,415-variant blacklist, for AJ and EUR.
Non-destructively: `--geno`/`--maf`/`--hwe` are per-variant and mutually independent, so
`cohort_<ANC>_qc` minus the blacklisted IDs **is exactly** what step 6 pass 2 would produce. Three
plink2 calls per stratum to `data/merged/premise_test/`, nothing in `by_ancestry_qc/` touched, no
pipeline rerun, the unfiltered baseline preserved. **This is the measurement no run had ever made.**

Variants actually removed: AJ 7,997,769 → 7,993,751 (**4,018**); EUR 7,538,809 → 7,534,424
(**4,385**). Both strata matched 100% of eigenvec rows to the manifest.

**EUR — the filter works, and the effect is large.**

| PC | unfiltered | filtered |
|---|---|---|
| PC2 | **0.757** | **0.011** |
| PC3 | 0.155 | 0.008 |
| PC1 | 0.080 | 0.032 |
| **max over PCs** | **0.757** | **0.036** |

**95.2% reduction in the worst PC**, and no PC above 0.036 afterward. **0.058% of EUR's variants
carried essentially all of the callset structure in its PCs.** PC5/PC6 rose slightly (0.004→0.034,
0.003→0.036) — residual variance reshuffling once the dominant axis is gone, negligible in
absolute terms and not a finding.

**AJ — the filter does nothing. 0.984 → 0.962, a 2.2% reduction.** PC1 remains the cohort label.
Predicted in advance from the structure: no AJ cell was ever evaluated, so the blacklist contains
no AJ-derived frequency flags, and AJ's HWE stage returned 97 failures against 400 expected by
chance. The 4,018 variants removed from AJ were EUR-derived, and **they do not transfer.**

**RESOLVED — the 0.748 was EUR's, mislabelled as AJ's.** This closes the "UNEXPLAINED" item in
the entry below.

| | baseline | after |
|---|---|---|
| claimed for **AJ** (147 variants, \|dAF\|>0.10) | 0.748 | 0.055 |
| measured **EUR** today (4,385 variants) | **0.757** | 0.011 |
| measured **AJ** today (4,018 variants) | 0.984 | 0.962 |

The claimed baseline matches EUR to within 0.009 and misses AJ by 0.236; the claimed collapse
matches EUR's and misses AJ's non-effect entirely. A second, independent route to the same
conclusion: 93.6% of the July list survives in today's, so those 147 variants are almost certainly
a subset of today's 4,415 — and removing them **plus 4,268 others** moves AJ by 0.022. It is
arithmetically impossible for the 147 alone to have taken AJ from 0.748 to 0.055.

**So nothing regressed and the four-callset merge degraded nothing.** The earlier "AJ is worse than
the 0.748 that motivated the filter" was comparing AJ against an EUR number.

**REFUTED — `af_concordance_build.py:5`'s foundational claim.** "diag_af_crossstratum proved
non-circularly that the callset separation in PC space is variant-intrinsic: dropping the 147
variants EUR flagged at |dAF| > 0.10 took AJ's eta^2 from 0.748 to 0.055." The cross-stratum
transfer is the entire non-circularity argument, and it **fails by direct measurement** with 30×
more variants. This is no longer merely unsourced; it is contradicted. **The filter's real
justification is the one it demonstrates in EUR, which is more than sufficient** — and it is a
first-party, reproducible measurement rather than a citation to a missing script.

**CONSEQUENCE — there are two problems, not one, and only one is solved.**
1. **EUR: solved.** Apply the blacklist. 95% of the callset signal in the PCs is removed by
   excluding 0.058% of variants.
2. **AJ: unsolved, and the AF filter is the wrong instrument.** It cannot measure AJ (the
   `wgs_harm` arm is 44, under any GWAS-relevant floor) and EUR-derived flags do not transfer.

**LEADING HYPOTHESIS for AJ, untested — the split is REAL sub-continental structure, not technical
breakage.** Two independent observations point the same way: AJ shows no HWE excess (0.24× chance
expectation, i.e. below noise), and EUR-derived flags have no effect. AJ is `wb_dwgs` 1,324 vs
`wgs_harm` 174; Ashkenazi is a bottlenecked population and the label is a classifier output, so
the two cohorts' AJ sets may genuinely differ in origin or degree of admixture. If so PC1 is doing
exactly what a PC should do and belongs in the model — the cost is the cohort/phenotype
collinearity, which is a **study-design limitation to state in the methods, not a bug to filter.**

**Next test for AJ:** project its samples onto the GenoTools reference panel and check whether
`wgs_harm`-AJ and `wb_dwgs`-AJ separate on a *reference* axis. Reference axes are built from known
populations and are independent of our callsets — separation there means real ancestry; separation
only in the within-stratum PCA means technical.

**Method note, reusable.** The three step-6 QC filters are per-variant and independent of one
another, so a variant-exclusion experiment needs no step-6 rerun and no overwrite: subtract the IDs
from `cohort_<ANC>_qc`, re-prune, re-PCA to a scratch prefix. Exact, not approximate. This matters
because **re-running step 6 is a one-way door** — it overwrites `cohort_<ANC>_qc.*`, the eigenvecs,
and (via `tee` without `-a`) `step6_summary.txt`, and it poisons any future rebuild of the
blacklist, since `af_concordance_build` derives frequencies from those same files.

---

## 2026-08-19 — 6a ran clean and answered ONE of its three questions. Sentinels cleared.

**Did.** Read job 27697096's log end to end, resolved both sentinel-loci hits against the
per-cell `.afreq`/`.snplist` intermediates, ran the July-list overlap, and tested three
hypotheses about AJ's eta². Two of the three died. No pipeline code changed.

**Provenance of the run is CLEAN, by a six-minute margin.** `sacct`: submitted 16:54:16, started
**16:55:51**, ended 17:02:32, 6m41s, COMPLETED 0:0. The grain was rewritten at **16:49** and the
`.o` prints `grain: 12,495 samples`, so it read the rebuilt grain. This was luck, not design —
the earlier attempt that PROJECT_LOG recorded as RUNNING at 14:58 was resubmitted by hand after
the rebuild. A 16:40 submission would have produced a plausible ~4,400-variant list from the
11,918-row grain, exit 0, with nothing in the log to say so.

**FOUND — AJ was never evaluated. Not one AJ frequency comparison ran.** Every AJ cell reads
`fewer than 2 callsets at >= 100, skipped`; the binding case is `AJ/control` at `wb_dwgs=638,
wgs_harm=44` against `MIN_CELL=100`. The per-cell detail table has three rows, all EUR. AJ's only
contribution is its HWE stage: **97 failures against 400 expected by chance** — below chance, so
no signal.

**So the 4,415-variant list is EUR-derived by construction, and this run says nothing about
whether it moves AJ's 0.984.** That was the whole question. The only evidence it would is
`diag_af_crossstratum`, which does not exist (below).

**FOUND — the BY CALLSET PAIR table is missing the arm that decides the liftover question.**

| pair | cells | shared | flagged | % |
|---|---|---|---|---|
| `wb_dwgs`\|`wgs_harm` | 2 | 15,077,618 | 4,018 | 0.027 |
| `divco_hs`\|`wgs_harm` | 1 | 7,538,809 | 1,698 | 0.023 |
| `divco_hs`\|`wb_dwgs` | **absent** | — | — | — |

No cell has both natively-called callsets above 100 (`EUR/control` is `divco_hs`=46). The script
prints its decision rule — "comparable rate in the both-native pair ⇒ calling/mapping; only
`wgs_harm` pairs ⇒ liftover-specific" — directly above a table that cannot support either branch.
Read casually it looks like the liftover-specific verdict. **The both-native rate is not low, it
is unmeasured.** Same shape as the 2026-08-17 incident: a conclusion resting on the absence of an
expected thing, there a log string, here a table row.

**Both sentinel hits RESOLVED, and neither is the filter eating real biology.**

`CR1 chr1:207521012:T:C` — frequency-flagged in `EUR/AD`. Four independent arms agree at
0.170–0.212; `divco_hs` sits at **0.365** with **OBS_CT 156 of 242 (64.5% call rate)** while every
other arm is complete. Inflated ALT frequency plus non-random dropout is lost reference calls, not
a frequency difference. Confirming: it fails HWE in **all three** testable callsets
(`AJ/wb_dwgs`, `EUR/wb_dwgs`, `EUR/wgs_harm`). With `keep-fewhet` those are het-EXCESS failures,
and het excess in three independent callsets is paralog/CNV collapse — consistent with CR1's known
copy-number variability. Excluding this variant does not remove the CR1 locus; neighbouring
well-behaved variants carry the signal.

`LRRK2 chr12:40227079:C:T` — **not frequency-flagged at all.** Max spread across arms is
0.058−0.023 = 0.034, under the 0.05 threshold. It entered via the HWE stage only, failing in both
testable callsets.

**The cell design held.** The tripwire fired, we looked, and the answer is two intrinsically
unreliable variants. Uniform HWE failure across every testable callset is the *opposite* of a
cohort artifact — a cohort artifact fails in one callset.

**But it exposes a composition point for the methods:** the 1,069 HWE-only additions are general
variant QC, not cohort-artifact removal, and they ride into the same union. The docstring frames
HWE as mechanism-based confirmation of the frequency test, which is true of the 840 overlapping
variants and not of the other 1,069. Both sentinel hits came in at least partly through the HWE
side. Separate them when writing this up — the HWE set would survive any change to cohort
composition, and the frequency set would not.

**HWE fingers `wgs_harm`, and pointedly not AJ.**

| stratum | callset | controls | fail | exp by chance | ratio |
|---|---|---|---|---|---|
| EUR | `wgs_harm` | 328 | 1,680 | 377 | **4.5×** |
| EUR | `wb_dwgs` | 3,064 | 132 | 377 | 0.35× |
| AJ | `wb_dwgs` | 638 | 97 | 400 | 0.24× |

`wgs_harm` clears its expectation 4.5× with the *smallest* sample of the three — real het excess,
consistent with mismapping in the lifted callset, and the strongest single result in the run.
**AJ shows no HWE excess anywhere**, so the mechanism-based test does not corroborate technical
breakage in the stratum whose eta² motivates the entire filter. Caveat: `AJ/wgs_harm` controls are
n=44, under `MIN_HWE_CONTROLS=50`, so this is an absence in the one AJ callset large enough to
test rather than across AJ. That keeps "AJ's PC1 split is real sub-continental structure, not
artifact" alive — and filtering real structure removes real signal.

**KILLED — "`divco_hs` is the broken one."** CR1's `divco_hs` outlier reading suggested the native
callset was at fault, which would have undercut the liftover story. Tested by using `wb_dwgs` as a
third-party outgroup across the 1,692 `EUR/AD` flags with all three AFs present:

| | `divco_hs` | `wgs_harm` |
|---|---|---|
| further from `wb_dwgs` (count) | 887 (52.4%) | 805 (47.6%) |
| mean \|AF − `wb_dwgs`\| | 0.0758 | **0.1146** |
| mean OBS_CT | 206 / 242 (85%) | 1,304 / 1,314 (99%) |

`divco_hs` is further slightly more often, but `wgs_harm` is further by **1.5× the magnitude** and
accounts for ~60% of total deviation from the native outgroup. CR1 was a genuine `divco_hs`
dropout case that does not generalize. **The liftover story survives, with two independent
supports** (this, and HWE above). Caveat on the method: this crosses disease (`divco_hs` AD vs
`wb_dwgs` control), which the cell design forbids for *building* a list — acceptable as a
technical diagnostic over 1,692 variants because real disease effects are sparse, but it must
never become a filter input.

**NEW, unresolved — `divco_hs` averages 85% call rate at flagged variants against `wgs_harm`'s
99%.** `z` uses each arm's per-variant `OBS_CT`, so it corrects *precision*; non-random dropout
biases the AF estimate itself and no z correction reaches bias. Some unknown fraction of the 1,698
`EUR/AD` flags is `divco_hs` missingness rather than frequency disagreement.

**`divco_hs` and `br_dsnwgs` are essentially unexamined by this run.** No `hwe_divco_hs.*` file
exists — `EUR/control` `divco_hs`=46 against `MIN_HWE_CONTROLS=50`. So `divco_hs` got one
frequency comparison and zero HWE tests, `br_dsnwgs` got neither. Of four callsets the run
meaningfully examined two, and both exclusions are threshold misses of 4 samples.

**The list is STABLE against the fourth callset.** Overlap with `july28_3callset`: **4,292 of
4,587** July variants survive, and 4,292 of the new 4,415 were already July's — **97.2% of the new
list is the old list**, only 123 genuinely new. Exactly what BR at 0.7% of the cohort predicts,
and worth stating in the methods as a reproducibility result: the artifact set is a property of
the callsets, not of cohort composition. (Counted with the awk hash from 2026-08-15, not
`sort`+`comm`.)

**And that sharpens the AJ puzzle rather than softening it:** substantially the SAME list, yet AJ
now starts at 0.984 rather than 0.748. A near-identical intervention against a very different
baseline is more consistent with 0.748 having been a different quantity than with the merge having
degraded anything.

**The 0.748 / 0.055 / 147 numbers have NO derivation in this repo.** They came from
`diag_af_crossstratum`, which is in neither the working tree nor **any commit** — `git log --all
-S "0.748"` and `-S "0.055"` both return empty. The numbers survive only as prose in
`af_concordance_build.py:5`, `HANDOFF.md:195`, `wgs_core_draft.ipynb`, and this log. **This is
structurally the same object as the "~748 Rush↔ROSMAP duplicates"** — a docstring number with no
derivation, which this project already disproved once. Treat it as an assertion.

**KILLED — "0.984 is inflated because the manifest now has more callset levels."** The manifest
carries five `source_callset` levels where July had three, and eta² cannot decrease as levels are
added. Recomputed AJ PC1 five ways on `results/retained_samples_manifest.csv`:

| grouping | levels | n | PC1 eta² |
|---|---|---|---|
| as published | 5 | 1,518 | 0.984 |
| drop `br_dsnwgs` | 4 | 1,505 | 0.984 |
| drop BR + fused rows | 3 | 1,504 | 0.984 |
| fused → `divco_hs` | 3 | 1,505 | 0.983 |
| fused → `wgs_harm` | 3 | 1,505 | 0.982 |

Flat. The 0.984 is real PC geometry, and this independently reproduces the "without BR = 0.984"
figure by a different route.

**UNCONFIRMED, recorded because it is cheap to test and would dissolve the mystery:** today's
**EUR** PC2 is 0.757 as published and **0.747** with the fused rows folded into `wgs_harm`.
Against a remembered 0.748. If the July measurement was EUR's and got attributed to AJ in the
docstring, then nothing regressed — EUR is unchanged at ~0.75 and AJ was simply never measured
unfiltered in July. **This is a coincidence of one number and nothing more**; it cannot be
confirmed without `diag_af_crossstratum`. Do not propagate it into `HANDOFF.md` as a finding.

**The validation gap is wider than previously stated.** The 147-variant experiment used
`|dAF| > 0.10`. The production list is 30× larger at 0.05, plus 1,069 HWE-only variants and 23
discordance rows that were never in that experiment — and its negative control (147 *random*
variants leaving eta² at 0.738) only calibrates at that size. Reconstructing the diagnostic on
current data is the one measurement that settles whether the filter earns its keep in AJ: extract
only the strict-threshold subset, recompute AJ's PCs, measure eta², with a size-matched random set
as control.

**Scope limit worth a methods sentence:** 7,538,809 variants could be evaluated by some powered
cell, against 20,642,023 in some stratum's QC set — **13,103,214 (63.5%) were never evaluated and
are kept by default.** The filter is far narrower than "genome-wide".

**KILLED — the `.afreq` parse hazard.** The `.afreq` files looked to have inconsistent column
counts under `grep`, which would be a silent-wrong-answer risk if parsed positionally. All headers
are identical (`#CHROM ID REF ALT PROVISIONAL_REF? ALT_FREQS OBS_CT`) and
`af_concordance_build.py:112` resolves columns by name via `hdr.index()`. Terminal rendering, not
a defect.

**Guard design — two corrections to what was proposed first, both found by looking at the actual
artifacts.**
1. **mtime against the merge is the WRONG check for the grain.** Step 6's guard compares the
   exclusion list to `cohort_merged.bed`, which is right there because the list derives from the
   merge. The grain does not — it derives from the manifest. Grain 2026-08-18 vs merge 2026-08-15
   passes in every failure scenario, including the one that nearly happened. **The load-bearing
   check is `grain rows == manifest rows`**: 12,495 == 12,495 now, and 11,918 != 12,495 at 15:00.
2. **The provenance stamp must be a SIDECAR, not a header.** `wc -l` on the list is 4,415,
   exactly the stated total, so the file is bare variant IDs. Writing `##` lines into it feeds
   them to `plink2 --exclude`. Write `exclude_af_concordance.txt.prov` and have step 6 validate
   that.

Note the `.py` already computes the grain's sample count — it prints `grain: 12,495 samples` to a
log nobody diffs. Moving that value into a sidecar is a two-line change.

**Minor, noted not chased.** The `EUR/AD` flags file yielded 1,706 `chr`-prefixed rows against the
log's reported 1,698 flagged (1,692 with all three AFs + 14 incomplete). An 8-row gap, possibly
duplicate IDs. Also `3,326 + 1,069 + 23 = 4,418` against a stated total of 4,415, consistent with
3 discordance rows already being in the union — this project's arithmetic usually closes exactly,
so both are worth a glance.

---

## 2026-08-18 — the 748/846 mystery is SOLVED, and the grain is rebuilt

**748 and 846 are `clinical_core.py` §8's multi-cohort ENROLLMENT counts.** Today's cluster run
prints them directly:

```
846 donors are enrolled in more than one cohort
    amp_ad_divco + amp_ad_rosmap    748
    amp_ad_divco + amp_ad_mayo       98
```

DivCo's largest contributing group is Rush (875 donors) and ROSMAP is Rush-run, so
"Rush↔ROSMAP" names the 748 exactly. **Neither number was ever a genotype-duplicate count.**

**The full chain, all four numbers now sourced:**

| quantity | count | what it measures |
|---|---|---|
| donors enrolled in >1 cohort | 846 (748 + 98) | clinical enrollment |
| donors with >1 genome | 309 | crosswalk, i.e. actually sequenced twice |
| KING duplicate clusters | 302 | genotype identity |
| genomes in those clusters | 621 | " |
| duplicate drops | 319 | 621 − 302 |

Enrollment overlap is a **superset** of sequencing overlap: being in two cohorts does not mean
being sequenced twice. 319 was never a shortfall against 748 — the docstring compared two
different populations. This supersedes the 2026-08-17 entry's "not currently decidable"; it was
decidable, it just needed a cluster-side `clinical_core.py` run, which needed the metadata fix
below.

**Grain rebuilt: `analysis_grain.csv` 12,495 rows × 22 cols**, from today's manifest. All checks
zero violations. Reconciliations that hold: grain rows = manifest rows = 12,495; ancestry totals
match the manifest stratum-for-stratum; and the status table's 12,582 retained genomes minus the
**87 fused** `divco_hs|wgs_harm` rows (two genomes, one merged sample) = 12,495.

**BR-DSNWGS in the grain:** 71 PD, 21 control, 3 null across its 95 retained — so BR is
predominantly PD cases, which is what makes its ~50% missingness matter for the primary contrast
specifically rather than diffusely.

**Contrast viability:** 17 of 44 contrasts reach ≥100 per arm (within_cohort 17, cross_cohort 15,
partial 12). AJ now clears it for `PD@amppd vs control@amppd`, `PD vs DLB` and `PD vs control`;
previously EUR only.

**FOUND — the two machines were reading DIFFERENT clinical inputs, and that is how the stale grain
got BR rows no cluster run could produce.** Neither had the full set:

| directory | laptop | cluster (before) |
|---|---|---|
| `amp-pd-genomics/metadata` | 5 | **0 — directory absent** |
| `WGS_Harmonization/metadata` | 5 | 5, **but not the same 5** |
| `DivCo_HS/metadata` | **2** | 5 |

The AMP-PD metadata was entirely missing cluster-side, and `WGS_Harmonization` matched on *count*
while lacking `MSBB_biospecimen_metadata.csv` and `WGS_sample_QC_info.csv`. **A file count is not
a file list** — that check was run and passed and was wrong. Fixed by pushing the AMP-PD directory
and those two files up. DivCo was deliberately NOT pushed: the cluster's set is larger and §1–6
reads it correctly, so overwriting with the laptop's 2 could only regress it.

**Consequence to enforce:** `clinical_core.py` must run on the CLUSTER only. The laptop is short
three DivCo files and `all_chrs_merged.psam` (`HANDOFF.md` issue #4), so a locally-built grain is
not comparable to a cluster-built one. Issue #4 currently offers the local workaround; it should
forbid local runs instead.

**Minor, noted not chased:** 4 `wgs_harm` genomes unresolved to a donor; 214 dx conflicts and 10
pheno conflicts among multi-cohort donors (expected — §12 reconciles); 2 cross-source sex conflicts
left for step 1's sex-check; and `age covariate: NOT FOUND`, the known limitation from
`07_gwas.sh`'s header — AMP-AD gives age at death, AMP-PD age at baseline, and they are not
commensurable.

---

## 2026-08-18 — unfiltered eta² baseline: TWO different problems, not one

**Did.** Step 6 pass 1 re-run genuinely unfiltered (job 27602590, `AF-concordance exclusion:
NONE`), then `ancestry_qc_manifest.py` → 12,495-row manifest (12,361 with PCs; 134 in the five
sub-50 strata), then `review/plot_pcs_by_callset.py`.

**Callset composition of the retained manifest:** wb_dwgs 9,919 · wgs_harm 1,701 · divco_hs 693 ·
br_dsnwgs 95 · fused `divco_hs|wgs_harm` 87. Sums to 12,495. (The `!! 87 rows with MULTI source`
warning is that fused group — expected, not a finding.) MDE fell 50 → 45 after exclusions and so
lost its PCs; five strata now have none (CAS/EAS/FIN/MDE/SAS), previously four.

**The eta² table splits cleanly once you remove BR and re-measure:**

| stratum | nBR | PC | eta² all | eta² without BR | BR max \|z\| |
|---|---|---|---|---|---|
| AJ | 13 | PC1 | 0.984 | **0.984** | 0.3 |
| EUR | 76 | PC2 | 0.757 | **0.758** | 1.2 |
| CAH | 1 | PC4 | 0.309 | 0.308 | 0.4 |
| AMR | 3 | PC8 | 0.203 | 0.207 | 1.5 |
| AAC | 1 | PC2 | 0.877 | **0.016** | **39.7** |
| AFR | 1 | PC1 | 0.708 | **0.121** | **19.6** |

**PROBLEM 1 — real callset structure (AJ, EUR).** Removing BR changes nothing, and the scatter
shows two disjoint clouds along AJ's PC1: `wgs_harm` versus `wb_dwgs`, the documented
AMP-AD↔AMP-PD axis. **This is what `af_concordance_build` exists for, and it is warranted.**

**PROBLEM 2 — two single-sample artifacts (AAC, AFR).** One BR sample each, at 39.7σ and 19.6σ.
Drop that single point and the apparent callset signal collapses 0.877 → 0.016 and 0.708 → 0.121.
**The AF filter cannot touch this** — it is not a variant-frequency effect. Mechanism is the
50%-missingness finding with a concrete consequence: where BR is n=1, that sample contributes
1/226 = 0.4% missingness per variant, far under `--geno 0.05`, so every variant BR lacks stays in
and the one sample is missing at half of them with no same-callset neighbours to average against.
It projects to the edge of the map. **This refines the 2026-08-17 assessment** that BR's
missingness was "a limitation, not a bug" — true for EUR and AJ, false in the n=1 strata, where it
is actively producing bad covariates. AAC (226) and AFR (191) are both above the ≥100-case
viability line, so those PCs would reach the GWAS.

**A `--mind` guard in step 6 is the WRONG instrument** — BR is ~50% missing in *every* stratum, so
any threshold that catches these two drops all 95. Excluding the two specific samples is
proportionate: 2 of 12,495.

**UNEXPLAINED — AJ's 0.984 is worse than the 0.748 that originally motivated the filter.** Nothing
we changed touches it (step 6 reads `cohort_merged` directly; `COMMON_GENO` affects step 4 only).
Either the July baseline was measured differently or something about the four-callset merge made it
worse. **Do not assume the 147-variant fix will reproduce its 0.055.**

**FOUND, blocking `af_concordance_build` — `$GRAIN` is stale and nothing checks it.**
`af_concordance_build.sh:76` tests only `[[ -f "$GRAIN" ]]`. `analysis_grain.csv` was built from
the OLD 3-callset manifest, so its cell membership and PCs predate today's run. Since the
(stratum × dx) cells come from the grain, and "disease held constant within a cell" is the sole
reason the filter does not delete APOE, running against a stale grain would produce a
plausible-looking and wrong exclusion list. **Third instance of the same shape:** a step consuming
an artifact by path with no provenance check. `clinical_core.py` §11–13 must run first.

---

## 2026-08-17 — step 6 pass 1 silently applied a STALE AF exclusion list. Guard added.

**What happened.** Step 6 was run over the four-callset merge expecting an unfiltered pass — the
baseline that decides whether the AF-concordance filter is still needed. Its log instead read:

```
AF-concordance exclusion (assoc + PCA): 4587 variants from .../exclude_af_concordance.txt
```

That file is dated **2026-07-28 23:13** and was built by the July `06a` run on the **THREE**-callset
cohort. It was never in the local repo — `HANDOFF.md` issue #3 said "nothing here builds it", which
was true locally and false on the cluster, and that is exactly why nobody looked. Exit 0, ~3 h of
compute, and a log line that reads like a success.

**How it was caught.** `grep -i "NONE at"` over the log returned nothing. The note is printed
unconditionally at `06_ancestry_qc.sh:110`, so its absence meant `AFX_NOTE` had taken the *other*
branch — i.e. the file existed. The absence of an expected string, not the presence of an error.

**Why the obvious recovery does not work.** Rebuilding the list from this run is
self-reinforcing: `af_concordance_build` derives frequencies from `cohort_<ANC>_qc.*`, which no
longer contain the 4,587 excluded variants. They cannot be re-flagged or un-flagged, so the new
list would inherit the old one's decisions and look clean regardless of the truth. Pass 1 has to
be redone genuinely unfiltered. The July list is preserved as
`exclude_af_concordance.txt.july28_3callset` — comparing it to the eventual rebuild measures how
much the fourth callset changed the answer.

**Changed — provenance guard in `06_ancestry_qc.sh`.** The step now REFUSES if the exclusion list
is older than `cohort_merged.bed`:

```bash
if [[ -f "$AF_EXCLUDE" && "$AF_EXCLUDE" -ot "${MERGED}.bed" ]]; then ... exit 1
```

A list predating the merge it is applied to is always wrong. Refusing rather than warning is
deliberate: the previous behaviour was effectively a warning, and it was not read.

**The general pattern, now three for three.** Every step in this pipeline that consumes an
artifact by *path* rather than by provenance has been bitten by a stale one: genotools refusing
to overwrite `relatedness_*`; `02a`/`02b` and `06a_*` lingering on the cluster after renames
because `rsync` without `--delete` cannot express a deletion; and this. The first two announced
themselves. **This one did not, and it silently altered the association set.** The exclusion
list carries no cohort, threshold, or date stamp — worth having
`af_concordance_build` write a header line, so the guard can check content rather than mtime.

**Also relevant, from the same log:** 6 of 11 strata produced PCs (EUR, AJ, AAC, AMR, AFR, CAH);
EAS/MDE/FIN/CAS/SAS hit `PRUNE_FAIL` because plink2 will not LD-prune at those sample sizes.
Consistent with prior runs. Long-range-LD exclusion applied correctly (15 regions).

---

## 2026-08-17 — the "~748 duplicates" was never a real number. RESOLVED, deferring to KING.

**Question.** `05_excludelist.py`'s docstring claimed duplicates are "where the ~748
Rush↔ROSMAP different-IID dups finally get caught," and also referenced "846 duplicates" on the
phenotype side. Step 5 reported 319. Was that a regression?

**No — 748 has no derivation anywhere in the repo.** It and 846 appear in that docstring and
nowhere else: not `clinical_core.py`, not `README.md`, not this log, not any script. (The only
other `748` in the tree is a coincidental eta² value in `06a_af_concordance.py`.)

**And 748 matches no computable quantity.** Duplicate drops relate to clusters as
`drops = members − clusters`, so 319 drops caps membership at 638 (all clusters of size 2). Any
cluster of 3 pushes membership *down*. 748 is unreachable as drops, members, or pairs (~336).

**Measured instead** — two independent counts, which agree:

| | count |
|---|---|
| KING clusters (step 5) | **302** clusters / **621** genomes / 319 drops |
| clinical crosswalk, donors with >1 genome | **309** donors / **620** genomes |

`621 − 302 = 319` exactly, so step 5 is internally consistent, and KING recovers essentially
everything the crosswalk knows about. The docstring was comparing against a number nobody had
checked.

**One real difference, and the decision on it.** KING's clusters are slightly FEWER and LARGER
than the crosswalk predicts — 17 clusters of 3+ against an implied 2 — meaning ~7 clusters unify
samples the crosswalk treats as distinct donors. Expected: `DUP_CUTOFF=0.354` cannot separate a
duplicate from an MZ twin pair, and the script's own comment says so. Candidates are MZ twins,
sample swaps, or one person enrolled under two donor IDs (plausible in a longitudinal cohort
like ROSMAP).

**DECISION: defer to the KING numbers and move on.** Genotypic independence is what the GWAS
model requires, and keeping one member per cluster is correct whether the pair is a twin or a
crosswalk error. At ~7 clusters in 12,495 retained samples it cannot move a result. Recorded
rather than chased.

**Changed.** Replaced the unsourced 748/846 in `05_excludelist.py`'s docstring with the measured
302 / 621 / 319 and 309 / 620, the MZ-twin caveat, and this decision — so the next person
compares against measurements instead of an assertion. The old number is explicitly flagged as
underived so it does not get resurrected.

**Small residual for the clinical side, not blocking.** The two BR-DSNWGS duplicates are
cross-PROGRAM (retained partners `<AMP-AD IID>` AJ and `<AMP-AD IID>` EUR carry AMP-AD-pattern numeric IIDs).
Which genome to keep is settled by the above; which *label* the surviving donor carries is not,
since one copy may be AD-by-neuropathology and the other PD. Two samples, so it cannot move a
GWAS, but it is a genuine AD-vs-PD label conflict and belongs in the §12 reconciliation.

---

## 2026-08-17 — `COMMON_GENO=0.005` confirmed; it changed the metric, not the outcome

**Did.** Re-ran steps 4 and 5 at `COMMON_GENO=0.005`.

**Fix works.** BR's mean F_MISS on the common set went **0.5030 → 0.0009**, against 0.0001 for
the other 13,237. Same order of magnitude now, so the call-rate denominator is finally the fair
one the step 4 header claims, and BR competes on merit in step 5's duplicate tie-break.

**But the prediction in yesterday's entry was wrong.** It said "expect these to move." They did
not:

| | 0.05 | 0.005 |
|---|---|---|
| total excluded | 839 | 839 |
| relative_2nd_deg | 498 | 499 |
| duplicate | 319 | 319 |
| sex | 21 | 21 |
| het | 1 | 0 |
| retained | 12,495 | 12,495 |

The single change is self-consistent: AAC `rel_drop` 2→3 and recorded `het` 1→0 are the same
DivCo AAC sample, now caught as a relative, which takes precedence over its het fail. `AAC pairs`
4→7 and `EUR pairs` 927→928 prove step 4 genuinely re-ran on a different common set — KING on
identical input is deterministic, so the pairs could not otherwise have moved.

**What this means, stated plainly:** close relatives are detected robustly whether KING sees ~5M
or ~10M variants, so the BR handicap never actually distorted the exclusion decisions. That
**retroactively validates the earlier 0.05 run's output** rather than condemning it. The fix was
still worth making — the metric now means what it says, and the next callset added at <1% of the
cohort will not silently inherit the problem — but it bought correctness of reasoning, not a
different answer.

**Worth keeping in mind:** an expected-and-benign property (a small joint call emitting fewer
sites) was being consumed as a quality signal, and the only reason it surfaced was that two
dropped samples' call rates matched to four decimal places. Nothing in the pipeline flagged it.

**Still open, unaffected by this change:** the `duplicate` count of 319 against
`05_excludelist.py`'s "~748 Rush↔ROSMAP" expectation, and the two cross-program BR duplicates
(`<AMP-AD IID>`, `<AMP-AD IID>`) whose phenotype labels need resolving in `genome_crosswalk.csv`.

---

## 2026-08-16 — every BR-DSNWGS sample is 50% missing on the common set (EXPECTED; the step 4 threshold is not)

**Found.** Steps 4 and 5 ran clean, but the two BR samples dropped as duplicates both carried
`callrate=0.497`. Two independently bad samples do not land within 0.0001 of each other, so it
was checked across the callset:

```
other: mean F_MISS 0.0001 over 13237 samples
BR:    mean F_MISS 0.5030 over 97 samples
```

**This is expected behaviour from BR, not a defect.** A 97-sample joint call emits nothing at
sites monomorphic in its 97 donors, and step 0's `--min-alleles 2` drops any site with no ALT.
In a merge, absent reads as *missing*, not hom-ref. The magnitude is about what theory predicts:
segregating sites scale with the harmonic number, and H(96)/H(13236) ≈ 5.14/10.07 ≈ 0.51
against an observed 0.497 present. Right number, right reason. **BR's data is fine.** (Rough
agreement only — the common set is already filtered, so it is not a random draw of sites.)

**What IS wrong is step 4's threshold.** `--geno 0.05` on the union cannot see a callset that is
0.73% of the cohort (97/13,334), so every variant BR lacks stayed in the common set. The script's
own header promises "the same call-rate denominator for every sample, so cross-dataset KING is
fair" — with BR at 50% missing it does not deliver that. Two consequences, both silent:

- **KING ran on a matrix where BR is half missing**, so relatedness involving BR is underpowered.
  BR relatives were more likely *missed* than falsely found — the direction you never notice.
- **BR can never win step 5's duplicate tie-break**, which picks by common-set call rate. Correct
  outcome for the two actual cases (partners at 0.99996 / 0.99995), wrong mechanism.

**Changed.** `04_relatedness.sh` now sets `COMMON_GENO=0.005`. Any threshold below the smallest
callset's cohort share (here 0.0073) forces the common set to variants present in **all four**
callsets, which is what "common" was meant to mean. The reasoning is in the script, including the
instruction to re-check the bound if a callset under ~0.5% of the cohort is ever added.

**Deliberately NOT changed — step 6.** Its per-stratum `--geno 0.05` builds the *association*
set, where requiring BR presence would discard half the variants for all 13,334 samples to
accommodate 97. That trade goes the other way. The consequence is a **limitation to state in the
methods, not a bug**: BR contributes to roughly half the association set, and it sits entirely on
the AMP-PD side of the primary contrast, with step 7's differential-missingness filter as the
only mitigation.

**Also found — the 2 BR duplicates are cross-PROGRAM.** Retained partners are `<AMP-AD IID>` (AJ) and
`<AMP-AD IID>` (EUR) — numeric IIDs, which is the AMP-AD pattern, not AMP-PD. So the same donor appears
in AMP-AD and in AMP-PD's postmortem callset. That is a **phenotype conflict in the primary
contrast**, not merely a genotype duplicate: one copy may carry an AD label and the other PD.
Resolve against `genome_crosswalk.csv` before the grain is rebuilt. BR IIDs are
`PM-MS_<5d>-BLM0-PVC-DWGS` (cluster dup19, n=3) and `PM-MS_<5d>-BLM0-PVC-DWGS` (dup27, n=2);
the `PM-MS_` prefix and the `gs://amp-pd-receipt-2026/mssm/` source both point at Mount Sinai,
which also contributes MSBB on the AMP-AD side.

**Superseded numbers.** The step 5 run below was on the 0.05 common set: 839 excluded (498
relative, 319 duplicate, 21 sex, 1 het), 12,495 retained, 95 of 97 BR retained. Steps 4 and 5
must both be re-run at 0.005; expect these to move.

**Unresolved from that run.** `05_excludelist.py`'s docstring anticipates "~748 Rush↔ROSMAP
different-IID dups"; the run dropped 319 duplicates total. Possibly two different definitions
(donor-level pairs vs genotype drops), possibly a regression — nothing in this log records what
the 3-callset run produced, so it is not currently decidable. Check before trusting the grain.

---

## 2026-08-15 — step 3 merge COMPLETE with all four callsets (job 27429821)

**Did.** `03_merge` COMPLETED 0:0. `cohort_merged` is now four callsets.

**Verified, three ways.** 13,334 samples (13,237 + BR's 97, exact); `97 of 97` BR IIDs present
in the merged fam; and the variant arithmetic closes exactly:

```
171,249,226   previous 3-callset merge
  1,247,829   BR-only IDs measured before the merge
─────────────
172,497,055   observed
```

A clean union — nothing dropped, nothing duplicated. The pre-merge ID check was
`15,537,603 of 16,785,432` BR IDs (92.6%) already present in `cohort_merged`, confirming the
step 2 normalization put BR in the same `chr:pos:REF:ALT` namespace as the other three.

**Method note — `comm` lied, and it took a second pass to notice.** The first overlap check
piped `sort -u` into `comm -12` and returned 804,919 with a `comm: file 2 is not in sorted
order` warning. Real answer: 15,537,603, an undercount by a factor of 19. Cause is locale
collation — under a UTF-8 locale `sort` orders punctuation differently than `comm`'s ordering
check expects, and these IDs are mostly punctuation. **Do not sort-and-`comm` variant IDs.**
Hash the smaller file in awk and stream the larger one: no sort, no locale, one pass, and it
prints numerator and denominator together:

```bash
awk 'NR==FNR{a[$2];next} ($2 in a) && !($2 in seen){seen[$2];n++}
     END{print n, "of", length(a)}' small.bim large.bim
```

**Everything downstream of the merge is now stale** and must be regenerated: relatedness, the
excludelist, the per-ancestry QC sets and PCs, and the BR rows of `analysis_grain.csv`.

**Then fixed `04_relatedness.sh`** — it listed `LBL_WGS`/`LBL_WB`/`LBL_DC` in two places and
would have carried the 97 BR samples through the merge with no ancestry label, dropping them
from every stratum without an error. Both sites now read one `LABEL_FILES=(...)` array, which
is what made the omission survivable in the first place: the input-existence check and the
label union were separate literal lists, so neither could catch the other being short.
`--threads 32` (three sites) is now `${SLURM_CPUS_PER_TASK:-32}` — same fix as steps 2 and 3.

**Not a concern for step 4's dedup:** 13,334 = 13,237 + 97 exactly, so no BR IID collided with
an existing sample. Had one, plink1.9 would have fused them the way it fused the 87 dual-source
DivCo/WGS_Harm genomes, and the total would have come up short.

---

## 2026-08-14 — step 2 run for BR-DSNWGS; step 3 rewritten to include it

**Did.** `norm_br` (job 27429119) COMPLETED 0:0. All four callsets now have normalized beds.
Rewrote `scripts/03_merge.sh`, 109 → 66 lines.

**Found — the merge list did not contain BR-DSNWGS.** `03_merge.sh` hardcoded `NORM_WB` and
`NORM_DC` as secondaries in a heredoc, so running it unchanged would have rebuilt the same
3-callset cohort and silently produced a valid `cohort_merged` with no BR samples in it.
The callsets are now a `SECONDARY=(...)` array at the top; adding one is a one-line edit.

**Changed in `03_merge.sh`:**
- `--threads` and `--memory` are now passed. plink1.9 sizes both from the *node*: its memory
  default is half of system RAM, which on a largemem node is approximately the whole 1500G
  request with no headroom left over. `--memory` is set to 90% of `SLURM_MEM_PER_NODE`. This
  is the same defect class as the step 2 `--threads` omission and the GenoTools worker bug —
  a tool reading the machine instead of the allocation.
- **Both merge passes now stage into `tmp_merge/` and are moved to `cohort_merged` only on
  success.** Previously pass 1 staged but the flip-retry pass wrote straight to the output, so
  a failed *retry* destroyed a good existing `cohort_merged` while a failed first pass did not.
- No `set -e`, deliberately, and now said so in the script: pass 1 is *allowed* to fail, since
  that is what triggers the flip retry. Added a `: "${NORM_WGS:?}"` guard instead, which is the
  part `set -e` was buying — catching a `source config.sh` that silently did nothing.

**Noted, not changed — the flip retry is probably now unreachable.** It fires on a plink
`.missnp`, which requires two filesets to disagree on the alleles at a shared variant ID. Step 2
encodes REF/ALT into every ID (`chr:pos:REF:ALT`), so a shared ID carries identical alleles by
construction. Kept as a cheap net (~10 lines after the rewrite, was ~35), with the reasoning in
the script so the next reader does not have to re-derive it.

**Next — `04_relatedness.sh` has the same omission one step later.** Lines 67 and 100 list
`LBL_WGS`, `LBL_WB`, `LBL_DC` and not `LBL_BR`: line 67 is an existence check, line 100 builds
the union of per-callset ancestry labels. Unfixed, the 97 BR-DSNWGS samples would enter
`cohort_merged`, carry no ancestry label, fall out of every stratum, and never be tested for
relatedness — without an error. Fix before running step 4.

---

## 2026-08-14 — step 2 consolidated to one script; the chr22 gate retired (with its result)

**Did.** Merged `02a_normalize_check.sh` and `02b_normalize.sh` into a single
`scripts/02_normalize.sh`, and deleted `02a`. This entry exists because deleting `02a` would
otherwise delete the only record of what it measured — the number below appears nowhere else
in the repo.

**What `02a` measured, and why it mattered.** Before normalization the three callsets' variant
IDs did not agree, so a union merge would have joined almost nothing. `02a` extracted chr22
from each filtered pgen, ran the *exact* pass `02b` would run genome-wide, and recomputed the
3-way exact-ID overlap:

| | chr22 3-way exact-ID overlap |
|---|---|
| before normalization (chr-stripped comparison) | ~83,905 |
| after `--ref-from-fa force` + `--set-all-var-ids '@:#:$r:$a'` | ~331,015 |
| allele-concordant positions (the ceiling) | ~331,015 |

The pass recovered essentially the whole ceiling, which is what authorized the genome-wide
re-normalize and the re-merge on 13,428 genomes. `--ref-from-fa force` was doing most of the
work: roughly 75% of the disagreement was REF/ALT swaps, not position mismatches.

**Also carried forward from `02a`'s header, because it is a live constraint:** plink2's
`--set-all-var-ids` defaults to `--new-id-max-allele-len 23 error`, which *aborts* on any
allele over 23 bp. WGS callsets are full of long indels. `1000 missing` is what lets the run
finish; variants beyond 1000 bp get a fallback `chr:pos` ID and simply do not align. This flag
must stay identical across every callset or the IDs stop being comparable.

**Why retiring the gate is safe.** A gate is worth its cost when the thing it guards is
expensive. It guarded a genome-wide pass plus a re-merge across 13,428 genomes. BR-DSNWGS is
97 samples and `02_normalize.sh` runs in minutes, so there is nothing left to gate — and `02a`
was hardcoded to `PF_WGS`/`PF_WB`/`PF_DC`, so it could not have taken BR-DSNWGS without being
rewritten anyway.

**Changed in `02_normalize.sh`** (was `02b`, 79 lines → 44):
- **`--threads` is now passed.** This was the only plink2 step in the pipeline that omitted it
  (step 1 passes 64, step 6 passes `$THREADS`). plink2 sizes its thread pool from the *node*,
  so an 8-CPU allocation on a 192-core node was spawning 192 threads. Same shape as the
  GenoTools worker bug, far less severe — threads, not processes, so `RLIMIT_NPROC` was never
  in play — but there was no reason for it.
- **`set -e`.** Closes the hazard named in `01_genotools.sh`'s header: without it a failed
  `source config.sh` leaves `REF_FASTA` empty and the script carries on. Replaces the manual
  `RC=$?` / error-echo block, which said nothing plink2's own stderr and the exit code did not.
- Dropped the `REF=` alias, `du -sh`, and the 9-line banner. No flag changed.

**Next.** Run step 2 on BR-DSNWGS, then check its normalized IDs against `cohort_merged.bim`
before step 3 — if that overlap is near zero, plink1.9 will merge BR-DSNWGS as ~all-new
variants, produce a valid file, and the damage will only surface later as `--geno` failures or
bad PCs. This callset has been the odd one out twice; the check costs seconds.

---

## 2026-08-14 — post-genotools pipeline order: KEEPING the current design (decision)

**Question asked.** Four candidate orderings for everything downstream of genotools:
(1) genotools per dataset with phenotype and `--all_variant`, then merge; (2) merge first,
genotools without `--all_variant` and no phenotype; (3) merge first, preliminary phenotype
(0=control, 1=AD/PD), `--all_variant`; (4) keep as-is.

**Decision: (4).** Reasoning, from the GenoTools source rather than its docs:

**What `--all_variant` actually is.** `main.py:87-93` expands it to `geno=0.05`,
`case_control=1e-4`, `haplotype=1e-4`, `hwe=1e-4`, `filter_controls=True`, **`ld=None`**.
`main.py:131-137` builds the run list only from non-`None` truthy args, so LD pruning does
**not** run and the association set is never rewritten.

**Overlap with step 6 is 2 of ~9 responsibilities**, and the least important 2:

| | step 6 | genotools `--all_variant` |
|---|---|---|
| variant missingness | `--geno 0.05` in stratum | `geno 0.05` in stratum — **exact duplicate** |
| HWE | `1e-6 keep-fewhet`, all samples | `1e-4`, controls only — **conflicting** |
| MAF | `--maf 0.01` | none |
| differential missingness | none | `case_control 1e-4` |
| haplotype mishap | none | `haplotype 1e-4` |
| LD prune | `1000kb 1 0.1`, long-range LD excluded | not run |
| PCA | `--pca 10` per stratum | not run |
| autosomes only | yes | no |

The HWE row is a collision, not a duplicate: `1e-4` removes strictly more than `1e-6`, and
without `keep-fewhet` it also removes heterozygote-deficient variants that step 6 keeps
deliberately. Whichever ran second would leave the harsher filter's result with no record of why.

**Step 6's product is the PCs, not the variant QC.** The filtering exists to make the
per-ancestry eigenvecs clean; step 7 consumes the eigenvecs. `--all_variant` produces none of
that, and the inputs to it are exactly what GenoTools cannot express: a **kb** LD window rather
than `--indep-pairwise 50 5 0.5`, and the long-range-LD/inversion exclusion. Without both, MHC
and 17q21.31 dominate a top PC and inversion haplotype rides into every GWAS as a covariate.

**Why (1) is disqualified on its own.** Rule-based per-callset filtering is fine — step 1's
`--var-filter --min-alleles 2 --snps-only --rm-dup` applies the same rule everywhere. *Data-
dependent* per-callset filtering is not: the merge is a union, so a variant dropped by HWE in
DivCo but kept in WGS_Harm becomes **missing for every DivCo sample** in `cohort_merged`. That
manufactures cohort-correlated missingness in a study where cohort is nearly collinear with
phenotype (AD from AMP-AD, PD from AMP-PD). Step 6 already handles this the right way round:
`--geno 0.05` computed *within* the stratum drops variants not genotyped across the callsets
present. Separately, within-callset ancestry strata are tiny — BR-DSNWGS is 1 AFR, 1 AAC,
1 CAH, 3 AMR — so per-callset HWE and case/control tests there are meaningless.

**Why (2) is worse than as-is.** It moves sample QC (callrate, sex, het) onto merged data. Call
rate distributions differ by callset and depth, so one threshold across the merge preferentially
drops samples from the lower-depth callset. Per-callset is the correct scope, and that is what
the current design does. Without phenotype, `case_control` also cannot run and `filter_controls`
has no control set — so `--all_variant` would be partly inert, which is the coherent reason (3)
pairs a phenotype with it.

**Why (3) is defensible but not taken.** Its real offer is the two filters step 6 lacks, plus
HWE-in-controls for free. The objection is that with cohort ≈ phenotype, `--test-missing` flags
variants differing in call rate *between callsets* and removes them under a case/control label;
some of that is batch cleanup, some is real signal at loci with differing coverage, and the
output cannot distinguish them. The right instrument for batch here is control-vs-control /
AF concordance, which this project already reasoned about — and which is currently missing from
the repo.

**Corrections made during this discussion, recorded because the wrong version was argued first.**
It was claimed that `--all_variant` LD-prunes and therefore destroys the association set, and
that HWE would run on all samples unless `--filter_controls` were passed. **Both are wrong** —
`ld=None` and `filter_controls=True` are set by that same expansion. The wrong version made
options (1) and (3) look disqualified for a reason that does not exist. The decision did not
change, but the reasoning behind it did, and only the corrected reasoning is load-bearing.

**Next, and it does not require restructuring.** `case_control` and `haplotype` are individually
flaggable (`--case_control 1e-4 --haplotype 1e-4`), so they can be run on the merged
per-ancestry sets **as diagnostics** — look at what they would remove before wiring either in as
a filter, precisely because the cohort/phenotype collinearity makes the removals uninterpretable
on their own. Promote `06a_af_concordance` regardless.

---

## 2026-08-14 — BR-DSNWGS genotools SUCCEEDED, real ancestry spread (sort fix confirmed)

**Did.** Re-ran step 1 with `--sort-vars` in step 1 and `GENOTOOLS_MAX_WORKERS=16`, both now
in the repo. The `.o` shows `Using 16 parallel workers for GridSearchCV (RAM: 256.0GB,
CPUs: 16)` — the wrapper reporting the allocation rather than the node, as intended. (Job ID
recoverable from `logs/submissions.tsv`.)

**Found — the labels are real:**

| label | n |
|---|---|
| EUR | 77 |
| AJ | 14 |
| AMR | 3 |
| CAH | 1 |
| AFR | 1 |
| AAC | 1 |

97 total, with a single CAH — the rate you would expect, against 97/97 before.

**The model is the control that proves the point.** Training balanced accuracy 0.9707
(95% CI 0.9490–0.9924), test 0.9838 (0.9751–0.9925) — statistically indistinguishable from the
broken run's 0.970/0.981. Identical model quality, completely different output. The classifier
was never the variable; only the study samples' column alignment was.

**Found — `analysis_grain.csv` is WRONG for BR-DSNWGS, not merely stale.** The grain carries
19 AFR / 67 EUR over 86 rows. The correct run gives 77 EUR / 14 AJ / 3 AMR / 1 AFR / 1 AAC /
1 CAH over 97 samples. Those are not reconcilable — 19 AFR vs 1 AFR is not a sampling
difference. Whatever produced the grain's BR-DSNWGS PCs and labels was not this pipeline in a
correct state, so both must be regenerated before step 7. This upgrades the standing open
question ("either from an earlier round, or values that will not match") to a definite defect.

**Minor, recorded to prevent confusion.** `Best Parameters` includes `xgb__lambda`, so the
classifier is XGBoost — while every output file is still named `..._umap_linearsvc_...`. The
filename is a GenoTools naming artifact and does not describe the model.

---

## 2026-08-14 — CAH root cause: variant ORDER mismatch (DIAGNOSED and FIXED)

**Found.** BR-DSNWGS's pgen lists chromosomes `1,10,11,…,19,2,20,21,22,3,…,9` — per-chromosome
files concatenated in **alphabetical** order, because `"10"` sorts before `"2"`. The reference
panel is numeric. GenoTools aligns the study matrix to the panel **by column position**:

```python
raw_geno.columns = col_names          # ancestry.py:247 — col_names is the PANEL's SNP order
```

The reorder that would make that safe —

```python
if not self.train:
    geno_snps = geno_snps[ref_snps.columns]     # ancestry.py:244
```

— is inside `if not self.train`, so it never runs when training from `--ref_panel`/`--ref_labels`.
All 209,068 shared columns were therefore standardized by another variant's mean/SD and
projected through another variant's loading. Column **counts** match, so nothing raises; the run
exits 0 and writes plausible output.

**The evidence chain, in the order it actually resolved:**

1. **The PC spread said "noise" before any file was examined.** Projected study samples vs the
   reference panel, sd of PC1–PC10:

   | | PC1 | PC2 | PC3 | PC5 | PC8 | PC10 |
   |---|---|---|---|---|---|---|
   | reference | 103.9 | 93.9 | 59.7 | 33.3 | 21.7 | 18.1 |
   | study (97) | 3.4 | 1.9 | 1.8 | 3.9 | 1.8 | 1.8 |

   The panel decays as eigenvalues must. The study is **flat** — equal spread on every axis is
   isotropic noise, not genomes. Real samples inherit the decay they are projected onto.
2. **Not collapsed to the origin** — study means are PC2 `24.47`, PC7 `12.05`, PC10 `-10.62`.
   A *fixed bias vector*: the same permutation hit all 97 samples identically. That distinguishes
   a scramble from missing data, which would have pulled everything to zero.
3. **The chromosome order then explained the offset exactly.** Walking the callset against the
   panel's ranks gives **exactly 2 rank drops**, at `19→2` and at `22→3` — precisely what a
   lexicographic block permutation predicts. (Adjacent-inversion *percentage* is the wrong
   metric here: 2 drops rearrange all 209,068 columns across 22 blocks. `diag_order.py` now
   reports drops and misplaced-column count instead.)

**Why CAH specifically.** `predict_admixed_samples` (`ancestry.py:779-842`) relabels any sample
closer to the overall centroid than to any ancestry centroid. A tight clump of noise near the
middle of the map trips that rule for every sample. CAH was never a claim about ancestry — it is
what the code says when the coordinates are meaningless.

**Ruled out, each with a measurement:**

- **hom-ref encoded as missing** (the previous entry's leading hypothesis). `--geno 0.1` removed
  **225,241 of 16,785,432 variants — 1.34%**. A variant-only VCF would have been gutted. The
  matrix is dense, which also rehabilitates the `call_rate` 0.979 column in `analysis_grain.csv`
  that the previous entry called doubtful. **The `plink2 --missing` command that headed this file
  is answered — do not run it.**
- **Collapsed SNP intersection.** 209,036 of the panel's 209,517 = **99.77%**, read from
  `FILTERED.br_dsnwgs_ancestry_umap_linearsvc_ancestry_model.common_snps` — the number GenoTools
  actually computed, not the one measured on `_filtered_bed`, which GenoTools discards.
- **plink2 drift from the cluster update.** Both the July success and the August failure ran
  **PLINK v2.00a3.3LM (3 Jun 2022)** — GenoTools resolves its *own* plink2 via
  `dependencies.check_plink2()`, so the module never reaches `ancestry.py`. The module *did* move
  (now `v2.00a6LM`, 17 Jun 2024) and it is irrelevant to this. Hypothesis raised and killed in
  one line of a log banner.
- **"The two problems share a cause."** They do not. Blocker 1 was the platform; this bug has been
  latent in GenoTools the whole time and bites only callsets whose order disagrees with the panel.
  It would have failed identically in January. Recorded because the opposite was argued at length
  earlier in this session — the cluster update is a coincidence, and reasoning from "both started
  after the update" was the wrong move. Two independent bugs surfacing together is not evidence
  of a common cause.

**Scope — all four callsets measured, not assumed:**

| callset | chrom ordering | rank drops | verdict |
|---|---|---|---|
| br_dsnwgs | **LEXICOGRAPHIC** | 2 | broken |
| wgs_harm | numeric | 0 | safe |
| divco_hs | numeric | 0 | safe |
| wb_dwgs | numeric | 0 | safe |

The other three worked only because their files happened to be numeric. Their July labels, the
merge, and the grain are unaffected — nothing downstream needs redoing.

**Changed.**
- `01_genotools.sh:89` — `--sort-vars` added to step 1, with the mechanism recorded inline.
  Applies to every callset, so the precondition need not be remembered.
- `01_genotools.sh:191` — `GENOTOOLS_MAX_WORKERS=16` baked in and echoed to the `.o`. The comment
  claiming 64 workers "leaves ~3x headroom" was **wrong** and is corrected: 64 failed with 13
  aborted workers. Job 27211436 is now reproducible from the repo.
- Added `scripts/diag_order.py` (variant order vs the panel) and `scripts/diag_cah.sh`
  (read-only postmortem of a genotools output dir).

**Two artifacts of the tooling, recorded so they are not misread again:**
- The ancestry step's `out_path` carries an `_ancestry` infix (`pipeline.py:86-97`), so the files
  are `FILTERED.<ds>_ancestry_*`. A first pass of `diag_cah.sh` used the pipeline prefix and
  reported every ancestry artifact as MISSING — that was a script bug, not evidence.
- `_common_snps.raw` and `.bim` are deleted by `clean_up()` (`ancestry.py:1087`) at the end of a
  successful run. Their absence is not a failure signal, and the NA fraction of the projected
  matrix cannot be measured post-hoc.
- `concat_logs` **appends**: `_all_logs.log` is cumulative across all 8 BR-DSNWGS attempts, not a
  record of one run. Variant counts are still safe to read (same input every time), but the file
  is not a single job's log.

**Upstream issues to file** (`dvitale199/GenoTools`), in value order:
1. `ancestry.py:247` positional column rename — silent data corruption with a clean exit. Hoist
   the `:244` reorder out of `if not self.train` and add a set-equality check that fails loudly.
2. `get_common_snps` (`utils.py:404`) prints `"Getting Common SNPs"` and never reports its overlap.
3. Worker sizing (`ancestry.py:516`) — `len(os.sched_getaffinity(0))`, an `RLIMIT_NPROC` term, and
   the cgroup memory limit rather than `psutil.virtual_memory().total`.

---

## 2026-08-12 — all 97 samples labelled CAH (LIVE ISSUE, not diagnosed)

**Found.** Job 27211436 exited 0 and produced ancestry labels, but every BR-DSNWGS sample is
**CAH**. That is what a broken projection looks like, and `--skip_fails` plus a silent
`get_common_snps` meant nothing in the pipeline objected.

**The wrapper is excluded as a cause — by evidence, not argument.** Two independent reasons:

1. **The model is excellent.** `Training Balanced Accuracy: 0.970` (95% CI 0.948–0.992),
   `Balanced Accuracy on Test Set: 0.981`. A degenerate classifier cannot score that.
2. **CAH is assigned downstream of the classifier**, by post-hoc geometry:
   ```python
   projected['label'] = np.where(projected['min_distance_ancestry'] == 'ALL', 'CAH', projected['label'])
   ```
   GenoTools computes a PC centroid per reference ancestry plus one overall centroid, and
   relabels any sample sitting closer to the **overall** centroid than to any ancestry
   centroid. Pure Euclidean distance in 50-dim PC space — no UMAP, no LinearSVC, no
   hyperparameters. `n_jobs` cannot move a sample's position in PC space.

So the classifier trained fine on the reference panel; the study samples are being projected
into that space badly, and all 97 land at the global centroid.

**Ruled out — variant ID mismatch.** BR-DSNWGS IDs are `1:66466:A:G`; the reference panel uses
rsIDs (`rs3748597`). ID overlap is **0**. But it is *also* 0 for `wgs_harm`, which worked in
July — because `get_common_snps` matches on composite `chr:pos:a1:a2` keys **in both allele
orientations**, never on IDs. Dead end; the control callset killed it immediately.

**Ruled out — positions or build.** Against the panel's 209,517 positions:

| callset | variants | position overlap |
|---|---|---|
| br_dsnwgs | 16,785,432 | **209,068** (99.8%) |
| wgs_harm | 53,214,908 | 207,706 (99.1%) |

BR-DSNWGS overlaps *better* than the callset that worked. Coordinates and build are fine.

**This is a regression, not a property of the data.** `analysis_grain.csv` already carries
19 AFR / 67 EUR for BR-DSNWGS, with PC spread indistinguishable from every other callset
(sd PC1 0.02074 vs 0.01977–0.02018 elsewhere). Something, at some point, projected these 97
genomes correctly.

**Leading hypothesis, UNTESTED — hom-ref encoded as missing.** If BR-DSNWGS's VCF is
variant-only rather than joint-called with reference blocks, hom-ref genotypes become
`missing` rather than `0/0`. At a panel of ~209k common SNPs each sample then has calls only
where it is non-reference; missing dosages are mean-imputed during projection; every sample
collapses to the global centroid — which *is* the CAH rule. It also explains why WGS_Harm is
unaffected (properly joint-called). BR-DSNWGS's pgen is the one input built by
`wgs_core.ipynb` f-strings rather than by a script, so it is the one conversion that never got
reviewed.

**Caveat against that hypothesis:** `analysis_grain.csv` shows `call_rate` 0.9794 / 0.9744 for
BR-DSNWGS samples. If those came from this callset, missingness is low and the story is wrong.
That column's provenance is doubtful (same grain carries PCs for a callset that had not
cleared step 1), but it is a real reason to measure before building on the theory.

**Next.** Run the `plink2 --missing` command in "Where we are right now". `F_MISS ≈ 0.75`
confirms; `≈ 0.02` refutes and moves the suspect to allele orientation (99.8% of positions
overlap, but `get_common_snps` needs REF/ALT to match in one of two orientations, and a
convention mismatch would empty the usable intersection while leaving position overlap
perfect).

**Upstream fix this exposes.** `get_common_snps()` prints `"Getting Common SNPs"` and nothing
else — no overlap count, no warning on an empty or tiny intersection, no handling of the empty
case. A collapsed intersection therefore produces a clean exit and plausible-looking output.
That is arguably more valuable to fix than the worker-sizing bug: a silent zero here is
indistinguishable from success.

---

## 2026-08-12 — why this bug never surfaced before (asked: "why now?")

**The full heuristic**, from `ancestry.py` on current upstream `main`:

```python
available_ram_gb   = psutil.virtual_memory().total / (1024**3)
gb_per_worker      = 3          # Conservative estimate for ~4k samples
max_workers_by_ram = max(1, int((available_ram_gb * 0.8) / gb_per_worker))
n_jobs             = min(os.cpu_count(), max_workers_by_ram)
```

It reproduces every observed worker count, so the model predicts rather than post-hoc explains:

| node | RAM | cores | `max_workers_by_ram` | predicted | observed |
|---|---|---|---|---|---|
| cn4304 | 503.4 GB | 128 | 134 | 128 | 128 ✓ |
| cn0041 | 755.3 GB | 192 | 201 | 192 | 192 ✓ |
| wrapper | 256 GB | 64 | 68 | 64 | 64 ✓ |

**Why no user has reported it.** It is a *memory* guard with no concept of process or thread
limits, and that is correct on a machine you own. Three conditions have to coincide:

- **A very large core count.** `n_jobs` scales with the machine, so a 16-core VM spawns 16
  workers and never sees this. The bug gets *worse on better hardware*, which inverts the
  usual intuition and is why nobody tests for it.
- **A low per-user `RLIMIT_NPROC`.** 1024 here. Cloud VMs, Terra, containers and workstations
  are typically tens of thousands or unlimited.
- **Shared nodes**, so the budget is also being spent by your other jobs.

Dataset size is irrelevant — the comment says `~4k samples` because the fit is on the
4,008-genome reference panel, identical for 97 samples or 10,000. Years of running this on
many cohorts genuinely could not have surfaced it.

Anyone who *did* hit it would see `SIGABRT`, `TerminatedWorkerError` and SLURM's bare
`ExitCode 1:0`, conclude the cluster was flaky, and retry until it worked — exactly what
happened here for four months at 3 successes in 18 attempts. **Silence is not evidence of
absence.** GP2 processes data "in multiple localities around the world"; some of those are
university clusters with biowulf's shape.

**Best answer to "why now": the code did not change, the platform did.**

---

## 2026-08-12 — BR-DSNWGS genotools SUCCEEDED (job 27211436)

**Winning configuration:** `GENOTOOLS_MAX_WORKERS=16` + `NUMBA_NUM_THREADS=1`, passed via
`--export`, on top of `scripts/genotools_capped.py`. Node cn4318, elapsed **01:03:02**, exit 0,
`Success.` printed. 8th BR-DSNWGS attempt (26843275, 27057269, 27103380, 27119797, 27162204,
27192530, 27202690, 27211436) and the first to complete.

**The scan that located the limit:**

| workers | grid-search time | outcome | aborted workers |
|---|---|---|---|
| 192 | 36:15 | failed | 63 |
| 64 | 42:57 | failed | 13 |
| 16 | 62:28 | **completed** | 0 |

**`TASKS_PER_WORKER` is bracketed by those runs**, no guessing required: 64 workers failing
means >1024/64 = 16 tasks per worker; 16 workers succeeding means <1024/16 = 64. So
16 < TASKS_PER_WORKER <= 64. A value of **32** with the half-budget rule reproduces exactly the
configuration that worked: `(1024 // 2) // 32 = 16`. Use 32 in the upstream patch;
`scripts/diag_threads.py` measures it directly if confirmation is wanted.

**The parallelism was never buying speed.** 16 workers finished the same 1080 fits in 62 min;
the July successes on ~128-192 workers took 1:11:59 end to end. Oversubscription was costing
time as well as risking the process limit — so capping workers has no performance downside to
trade off against, which strengthens the case for the upstream change.

**Correction — repo direction.** `GP2code/GenoTools` is a **fork**; `dvitale199/GenoTools` is the
source (GitHub API: `fork: true`, `parent: dvitale199/GenoTools`). An earlier entry had this
backwards. Upstream `main` (pushed Jul 2026) still contains the unmodified block, so there is no
upgrade path out of this — 1.3.6 and current main behave identically and the fix must be written.

**Upstream fix drafted** (see conversation for the full block): replace `os.cpu_count()` with
`len(os.sched_getaffinity(0))`, add a third limit from `resource.getrlimit(RLIMIT_NPROC)`, and
print all three terms so the binding constraint is visible. The RAM term has the same defect —
`psutil.virtual_memory().total` read 755GB on a node where the job held 256GB — and should read
the cgroup limit (`/sys/fs/cgroup/memory.max`, v1 fallback) the same way.

**Still to verify:** that
`FILTERED.br_dsnwgs_ancestry_umap_linearsvc_predicted_labels.txt` (`LBL_BR`) exists with ~97 rows.
`--skip_fails` means exit 0 does not by itself prove ancestry wrote its output.

---

## 2026-08-12 — the thread env vars were superstition; the worker COUNT is the lever

**Found — `ulimit -u` is 1024**, per-user and node-wide, shared across every job that user has
on the node. Measured with `srun --cpus-per-task=2 --mem=4g -t 5 bash -c 'ulimit -u'`. 192
worker processes is ~19% of the whole budget before counting a single thread, and each worker
is a full interpreter with numba and llvmlite loaded. 63 of 192 workers aborted — what hitting
the wall partway through spawning looks like.

**Found — joblib already caps threads inside each worker, so the five exports did nothing.**
From `joblib/_parallel_backends.py`:

```python
MAX_NUM_THREADS_VARS = [
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS", "NUMEXPR_NUM_THREADS",
]
default_n_threads = max(cpu_count() // n_jobs, 1)
```

`NUMBA_NUM_THREADS` **is** in the list (added by joblib PR #951, merged). With `n_jobs=192` and
`cpu_count()=192` the per-worker default is `max(1, 1)` = 1. So `NUMBA_NUM_THREADS=1` was
already set in every worker on job 27162204 — the run that had the exports commented out. The
joblib *docs page* omits NUMBA from the list; only the source is correct.

**Corrected — the "36,900 threads" figure from the previous entry was wrong.** `_launch_threads`
appeared at the abort site because it was the first thing to request a thread after the limit
was already reached, not because it was requesting 192 of them. Reading "named in the traceback"
as "cause" was the error. Two positions were taken and both were wrong before this one: the
wrapper was called a workaround, then the env vars were called load-bearing. What settled it was
measuring `ulimit -u` rather than reasoning about it — that measurement was available from the
start and would have skipped both detours.

**Changed.**
- Added `scripts/genotools_capped.py` — patches `os.cpu_count` and `psutil.virtual_memory` to
  report `SLURM_CPUS_PER_TASK` / `SLURM_MEM_PER_NODE`, then calls genotools' `handle_main()`.
  Drop-in for the console script; `GENOTOOLS_MAX_WORKERS` overrides. No `site-packages` edits.
- `01_genotools.sh` step 4 now calls the wrapper, and the five thread exports are **removed**
  with the reason recorded so nobody adds them back.

**Version audit (asked: did the venv drift?).** Partly, but it is not this bug. GenoTools 1.3.6
declares `umap_learn==0.5.3` (exact) but only floors elsewhere: `numba>=0.57.1`,
`numpy>=1.23.5`, `scikit_learn>=1.3.0`, `joblib>=1.3.0`. So pip resolved **numpy 2.0.2**, while
umap-learn 0.5.3 predates numpy 2 and is reported incompatible with numpy >= 1.25
(lmcinnes/umap#1007, discussion #1191). Pinning umap exactly while leaving numpy unbounded
permits a combination upstream never tested — worth adding `numpy<2` to `requirements.txt`
independently of this investigation. It cannot explain the intermittency: `ancestry.py` was
installed Jul 2 and the three successes were Jul 7, identical bytes and versions.

**Upstream fix worth making** (we maintain GenoTools): `ancestry.py:516` should use
`len(os.sched_getaffinity(0))` or honour `SLURM_CPUS_PER_TASK` instead of `os.cpu_count()`.
That makes the tool correct on any shared scheduler. `LOKY_MAX_CPU_COUNT` does not help because
`n_jobs` is passed explicitly, and `taskset` does not help because `os.cpu_count()` ignores CPU
affinity on Python 3.11. Once upstream lands, delete the wrapper.

---

## 2026-08-12 — control experiment: caps removed, failure reproduced at the exact source line

**Did.** Ran job 27162204 with the five thread-cap exports commented out — deliberately back
to the pre-fix configuration — to test whether the caps were necessary. Node cn0041,
64 CPUs / 256G allocated, node reporting 192 CPUs / 755.3GB. FAILED, ExitCode 1:0.

**Found — the abort site, named.** The worker tracebacks give the thread-creation call
directly:

```
numba/np/ufunc/parallel.py:522  in _launch_threads
numba/np/ufunc/parallel.py:649  in get_num_threads
umap/umap_.py:2295              in fit
```

UMAP asks numba for its thread count; numba lazily launches a pool sized to `os.cpu_count()`;
it does this **once per worker process**. 192 workers × 192 threads ≈ 36,900 against the
per-user thread limit. `what(): pthread_create has failed: Resource temporarily unavailable`,
**63 workers aborted** (vs 14 on the 128-core node — worse on the bigger machine, as the model
predicts). 53 leaked semaphores.

This is direct confirmation, not inference: the failing call is exactly the one
`NUMBA_NUM_THREADS` controls, and it was unset for this run.

**Corrected — capping `n_jobs` alone would not have been enough.** Earlier in this
investigation the wrapper (pinning `os.cpu_count` to `SLURM_CPUS_PER_TASK`) was called the
principled fix and the env vars a workaround. That was backwards on the arithmetic: 64 workers
× 64 numba threads is still ~4,100 threads. `NUMBA_NUM_THREADS=1` attacks the multiplier and is
the fix; the wrapper is optional belt-and-braces on top.

**Corrected — the "1:04" reading.** `01:04:32` in the `.o` is a wall-clock end time, not a
duration; `sacct` gives Elapsed **00:36:51**. Steps 1–3 took **36 seconds** (97 samples), so
essentially the whole job was grid search. This run did *not* get further than the others — it
sits in the same 31–43 min band as every previous failure. Timing signature to use going
forward: **~36 min measured from "Conversion complete", not from job start.**

**Changed.** Uncommented all five exports in `scripts/01_genotools.sh`, with the traceback and
the arithmetic recorded in the comment so the next person does not repeat the experiment.

**Next.** Sync the script to the cluster (via helix — `scp` to biowulf fails, see 2026-08-11),
confirm `active exports: 5`, resubmit.

---

## 2026-08-11 — `07_gwas.sh` misclassifies BR-DSNWGS as AMP-AD (found, NOT yet fixed)

**Found.** `scripts/07_gwas.sh:247` is `ispd = ($3=="wb_dwgs") ? 1 : 0`, written when there
were three callsets. BR-DSNWGS is AMP-PD, so it is wrong now. Feeding `srcok()` on line 242,
this both **excludes** BR-DSNWGS from every `@amppd` arm and — because the test is `!ispd` —
**admits** it into every `@ampad` arm.

Impact measured against `clinical_core_out/analysis_grain.csv` (86 BR-DSNWGS rows: 63 PD,
20 control, 3 no pheno):

| contrast | current behaviour |
|---|---|
| `PD@amppd:control@amppd` | drops 63 PD cases + 20 controls (AFR 10/8, EUR 53/12) |
| `AD@ampad:control@ampad` | admits 20 AMP-PD brain controls into the AMP-AD control arm |
| bare-dx contrasts | `case_pd+=ispd` undercounts → `case_pct_amppd`, `delta_amppd` and `confound_tag` all understated |

The `@ampad` contamination is the worse half: the script's own header calls those
within-cohort contrasts "confound-free — the trusted backbone".

`clinical_core.py:521` already has this right — `AMPPD_CALLSETS = {"wb_dwgs", "br_dsnwgs"}`,
used by §13 via `.isin()`. So the producer is correct and the awk consumer disagrees with it.
This is the "two implementations of pheno/covar" issue in `HANDOFF.md` §1, except one copy
has now silently gone wrong.

**Not yet urgent** — BR-DSNWGS is not in any merged genotype set, so no existing sumstats are
affected. It must be fixed before the merge, not after.

**Next.** One-line fix (`$3 in amppd` against a set declared once at the top of the script)
as a stopgap, then the consolidation below.

**Consolidation plan (agreed direction, not started).** §13 already writes
`clinical_core_out/{pheno,covar}/` with exactly the layout step 7 rebuilds — 28 pheno files
= 14 contrasts × 2 ancestries, matching step 7's `CONTRASTS` default, same headers. So this
is a deletion (~40 lines of awk), not a reimplementation.

- The only real coupling is FID. §13 *guesses* (`FID==0` for AMP-AD, `FID==IID` for AMP-PD);
  the awk resolves FID from the real `${QC}.fam` because plink2 matches on FID+IID and
  AMP-PD samples carry a non-zero FID.
- IIDs are globally unique (0 duplicates in 11,918 grain rows), so **if** plink2 accepts an
  `#IID`-only header and matches on IID alone, the whole FID problem disappears.
  **VERIFY THIS AGAINST plink 6-alpha FIRST** — it decides the design. Fallback is a ~5-line
  join against `${QC}.fam`.
- The summary can read `n_case` / `n_ctrl` / `case_pct_amppd` / `delta_amppd` /
  `confound_tag` / `viable_ge100` from `contrasts.csv` instead of recomputing them.
- Resulting boundary: **step 7 reads genotypes; it does not derive phenotypes.**
- Must decide where four runtime knobs live: `CONTRASTS` and `EXCLUDE_DUAL` would have to
  move into `clinical_core.py`; `MIN_ARM` and ancestry discovery stay in step 7.
- Free oracle: `contrasts.csv` is produced by the correct code, so step 7's counts must match
  it row for row. Today they do not (AFR `PD@amppd:control@amppd` is 679/809 there, 669/801
  from the awk). After the fix they must agree exactly.

**Open question.** Every BR-DSNWGS row in the grain already carries PC1–PC10 (86/86
populated), but BR-DSNWGS has not cleared step 1, and PCs come from step 6 on the merged
callset. Either those PCs are from an earlier pipeline round or the grain is carrying values
that will not match the run about to be produced. `07_gwas.sh:77` warns a stale grain
"silently runs on stale covariates". Resolve before step 7 runs.

---

## 2026-08-11 — `clinical_core_out/` was never on the cluster

**Found.** First run after the fix (job 27119797) died in ~5 seconds at step 1:
`Error: Failed to open --update-sex file`. The `.o` showed the path arrived intact via
`--export`; the file simply did not exist. `clinical_core.py` had only ever been run locally.

**Note for future debugging:** plink2 does *not* echo the filename in that error message, so
it reads as though the variable were empty. It is not. Check the `.o` header block, which
prints every resolved variable.

**Changed.** Copied `clinical_core_out/` to the cluster via helix (not biowulf — helix is the
transfer node). Verified: 98 / 1020 / 10419 / 1888 rows in the four `*_update_sex.txt`.

**Flagged, not resolved — cross-callset sex-file inconsistency.** The July runs sex-updated
`wgs_harm`, `wb_dwgs` and `divco_hs` from files under `<dataset>/metadata/`, not from
`clinical_core_out/`. One of those inputs (`WB-DWGS/metadata/all_chrs_merged_sexupd.psam`) is
the circular input this session's `clinical_core.py` fix removed — step 1 was reading its own
output. So BR-DSNWGS is about to be sex-updated from the *corrected* source while the other
three carry assignments from the old one. Resolving this likely means re-running step 1 for
all four callsets once BR-DSNWGS works.

**Also noted.** `scp` to biowulf fails with `Received message too long 1165128303`. That
integer is `0x4572726F` = the ASCII bytes `Erro` — biowulf's shell rc prints an error to
stdout on non-interactive sessions, corrupting the transfer protocol. Fix at the source with
`case $- in *i*) ;; *) return;; esac` at the top of `~/.bashrc`; until then, transfer via
helix or pipe through `ssh 'cat > dest'`.

---

## 2026-08-11 — GenoTools step 4 root cause: thread exhaustion (FIXED)

**Did.** Read the logs for the first time in four attempts. Job 27103380's `.e` gave the
whole chain.

**Found — the mechanism, end to end:**

1. `genotools/ancestry.py:516` — `n_jobs = min(os.cpu_count(), max_workers_by_ram)`, where
   `os.cpu_count()` and `psutil.virtual_memory()` (line 513) both report **node** totals and
   ignore the SLURM cgroup. On cn4304 that is 128 CPUs / 503.4 GB.
2. GridSearchCV (line 527, `n_jobs` passed explicitly) spawns that many loky workers.
3. Each worker starts its own numba/OpenMP pool, *also* sized to the node — ~16k threads.
4. `pthread_create` returns `EAGAIN` against the per-user thread limit.
5. The C++ layer throws an uncaught `std::runtime_error`, the worker `abort()`s on `SIGABRT`.
6. joblib raises `TerminatedWorkerError`; the parent exits 1.

Confirming string: `what():  pthread_create has failed: Resource temporarily unavailable`.

**Dead ends, recorded so they are not re-walked:**

- **"It is not an OOM, because SLURM would report OUT_OF_MEMORY."** The premise is wrong —
  SLURM reports that only when the cgroup killer reaps the job *step*. Here a child worker
  died and the parent exited cleanly, so SLURM logged a bare `ExitCode 1:0`. (The conclusion
  happened to be right, for the wrong reason: it was threads, not memory.)
- **`std::bad_alloc` / liblinear memory exhaustion.** Wrong — the exception is
  `std::runtime_error`, which is not what a C++ OOM looks like.
- **"BR-DSNWGS is the odd one out."** No. The fit is on the 4,008-genome **reference panel**,
  not on study samples, so it does identical work for 97 samples or 10,000. Only 3 of 18
  genotools jobs have ever completed; `wgs_harm`, `wb_dwgs` and `divco_hs` each needed 4–5
  attempts. BR-DSNWGS has simply lost the same coin flip four times.
- **"The `config.sh` refactor broke it."** No — reverting to the proven script changed
  nothing, because the script was never the variable.
- **"Bigger nodes fail, smaller nodes succeed."** No. cn0047 (success) and cn0044 (failure)
  are identical hardware, 2×48 cores / 765 GB. The variable is transient contention on a
  *shared* node — `AllocMem` was ~758/765 GB with `FreeMem` as low as 1.9 GB — not geometry.
- **Raising `--mem` or `--cpus-per-task` cannot help.** Job 27057269 requested 16 CPUs /
  128 G and failed identically, because the allocation never enters the calculation.

**Diagnostic that generalises:** every failure sits at 31–43 min, which is where GridSearchCV
begins; the three successes are 1:11:59, 1:11:59 and 3:56:00. **Dying at ~40 min means this
bug. Clearing 40 min means the fix held.**

**Changed.**

- `scripts/01_genotools.sh` — `OMP_NUM_THREADS` / `OPENBLAS` / `MKL` / `NUMEXPR` /
  `NUMBA_NUM_THREADS` all pinned to 1, immediately before step 4 so plink2 steps 1–3 are
  untouched. Collapses ~16k threads to roughly one per worker. Nothing is lost: joblib
  already parallelises at the outer level, so the inner threads were pure oversubscription
  even on the runs that survived.
- `submit.sh` — `--output` / `--error` now carry `%j`. Fixed log names are why three failures
  left no evidence and why this took four runs to diagnose.
- `--skip_fails` **kept**. It continues the pipeline past a failed QC step, which is wanted,
  and is unrelated to the ancestry method.

**Not done — the belt-and-braces option.** `n_jobs` still tracks the node's core count. If a
run dies at ~40 min again, the next move is a small wrapper in `scripts/` that pins
`os.cpu_count` / `psutil.virtual_memory` to `SLURM_CPUS_PER_TASK` / `SLURM_MEM_PER_NODE`
before calling genotools' `handle_main()` — versioned in the repo, no edits to
`site-packages`.

**Process note.** Old `.o`/`.e` for the three July successes are gone (`scripts/logs/` no
longer exists), so their worker counts are unrecoverable. Preserve logs before resubmitting;
`%j` now does this automatically.
