# Project log — ad-pd-gwas

Append-only. Newest entry first. **Nothing here is ever rewritten or deleted**, including
entries that later turn out to be wrong — a hypothesis that was ruled out is the most
expensive thing to rediscover, and the only place it can be recorded is a log.

Three documents, one job each. Keep them apart:

| doc | answers | written by |
|---|---|---|
| `HANDOFF.md` | what is true now | rewritten in place |
| `RUNLOG.md` | which jobs ran, and how they ended | generated — `bash scripts/runlog.sh --md > RUNLOG.md` |
| `PROJECT_LOG.md` (this file) | what we did, why, and what we learned | appended |

Entry format: date, a one-line title, then whichever of **Did / Found / Changed / Next**
apply. Record dead ends explicitly — they are the point.

**"Where we are right now" stays at the TOP** and is the only part of this file that gets
rewritten. Everything below it is append-only history.

---

## Where we are right now

**Both step-1 blockers are fixed, and step 1 is DONE for all four callsets.** The
post-`--sort-vars` run produced a real ancestry spread for BR-DSNWGS — 77 EUR / 14 AJ / 3 AMR /
1 CAH / 1 AFR / 1 AAC across 97 samples.

The CAH failure was a **variant ORDER mismatch**: BR-DSNWGS's pgen was ordered
`1,10,11,…,19,2,20,21,22,3,…` (per-chromosome files concatenated alphabetically) while the
reference panel is numeric, and GenoTools aligns the study matrix to the panel **by column
position, not by name**. Fixed with `--sort-vars` in step 1. See the entries below.

**Scope of the ORDER bug is BR-DSNWGS only — measured, not assumed.** All four callsets were
checked with `scripts/diag_order.py`: `wgs_harm`, `divco_hs` and `wb_dwgs` are numeric with
**0 rank drops**, so their ancestry labels and the merge are unaffected and need no redoing.
**This is not the same question as whether the grain is correct** — `analysis_grain.csv`'s
BR-DSNWGS rows are wrong for a separate reason and are item 1 below.

**Steps 1–5 are DONE for all four callsets.** `cohort_merged` is 172,497,055 variants × 13,334
samples; step 5 retains **12,495** (839 excluded: 499 relative, 319 duplicate, 21 sex).
95 of 97 BR-DSNWGS samples retained.

**The blacklist is BUILT, VALIDATED, and ready to apply.** `af_concordance_build` job 27697096
produced 4,415 variants from the correct grain, and the premise test (2026-08-19 entry below)
measured what it does: **EUR max eta² 0.757 → 0.036, a 95% reduction. Apply it.**

**AJ is a SEPARATE, UNSOLVED problem and the AF filter cannot touch it** — 0.984 → 0.962. Do not
wait on AJ before proceeding; it needs a different instrument, and the leading hypothesis is that
it is real sub-continental structure rather than artifact.

**Immediate next step — step 6 pass 2, but PRESERVE the baseline first. Step 6 is a one-way door.**

```bash
cd /data/CARDPB2/sysbio/wgs && source config.sh
mv ${MERGED_DIR}/by_ancestry_qc ${MERGED_DIR}/by_ancestry_qc_unfiltered   # instant, same fs
./submit.sh scripts/06_ancestry_qc.sh                                     # recreates the dir
```

Without the `mv`, pass 2 overwrites `cohort_<ANC>_qc.*`, every eigenvec, and (via `tee` without
`-a`) `step6_summary.txt` — and it poisons any future rebuild of the blacklist, because
`af_concordance_build` derives its frequencies from those same files.

Then, in order: `ancestry_qc_manifest.py`, `clinical_core.py` again (PCs changed → grain must be
rebuilt), then step 7.

**WITHDRAWN — the `MIN_CELL=40 MIN_HWE_CONTROLS=40` rerun proposed earlier.** `MIN_CELL=100` is
the GWAS-relevance floor (cf. "17 of 44 contrasts reach ≥100 per arm"), not a noise control, and a
cell that cannot field an analysis arm should not vote on a genome-wide exclusion list. The z gate
already handles noise and is sample-size aware.

**Still unmeasured: the both-native liftover verdict.** `divco_hs`↔`wb_dwgs` never ran — they share
only `control` and `other` as diagnoses, and `divco_hs` has 46 EUR controls. It costs **one
`plink2 --freq`** against `cohort_EUR_qc`, not a pipeline rerun, and then the power-matched pairwise
table is arithmetic on the `.afreq` files that already exist. Note the existing pair rates are NOT
power-matched: detection floors are 0.08 / 0.14 / 0.15 across the three cells, so they cannot be
compared as printed.

**The grain in this repo (laptop) is the WRONG one — 11,918 rows, dated 2026-08-10, BR-DSNWGS as
67 EUR / 19 AFR.** The correct 12,495-row grain exists only on the cluster, which is right
(`clinical_core.py` must run there). But `clinical_core_out/` sits in the project root beside the
code, so **any `rsync` of the project laptop→cluster overwrites the good grain with the bad one**
unless `clinical_core_out/` is excluded. Fourth face of the stale-artifact bug.

**The grain in this repo (laptop) is the WRONG one — 11,918 rows, dated 2026-08-10, BR-DSNWGS as
67 EUR / 19 AFR.** The correct 12,495-row grain exists only on the cluster, which is right
(`clinical_core.py` must run there). But `clinical_core_out/` sits in the project root beside the
code, so **any `rsync` of the project laptop→cluster overwrites the good grain with the bad one**
unless `clinical_core_out/` is excluded. Fourth face of the stale-artifact bug.

**Duplicates are settled, and the "~748" is sourced.** It was `clinical_core.py` §8's
multi-cohort ENROLLMENT count (846 donors in >1 cohort, 748 of them DivCo↔ROSMAP) — never a
genotype-duplicate count. Chain: 846 enrolled twice → 309 have >1 genome → KING finds 302
clusters over 621 genomes → 319 drops. Enrollment overlap is a superset of sequencing overlap.

**One small item for the clinical side, not blocking:** the two BR-DSNWGS duplicates are
cross-PROGRAM (partners `<AMP-AD IID>`, `<AMP-AD IID>` carry AMP-AD-pattern IIDs), so the surviving donor's
AD-vs-PD label needs picking in §12. Two samples — it cannot move a GWAS.

**Known limitation to state in the methods (not a bug):** BR-DSNWGS is a 97-sample joint call,
so it emits nothing at sites monomorphic in its own donors and contributes to roughly half the
association set. Step 6's per-stratum `--geno 0.05` cannot see a callset at 0.7% of the cohort,
and tightening it would cost half the variants for all 13,334 samples. BR sits entirely on the
AMP-PD side of the primary contrast; step 7's differential-missingness filter is the mitigation.

**Then, in order:**
1. **Regenerate the BR-DSNWGS rows in `analysis_grain.csv`.** The stored labels (19 AFR /
   67 EUR over 86 rows) are contradicted by the correct run — see the entry below. Those PCs
   and labels are wrong, not merely stale.
2. Add the order check as a pre-flight gate in `01_genotools.sh` — deferred until a run
   succeeded from the reverted baseline, which has now happened. This bug exited 0, wrote its
   output, and scored 0.97 model accuracy while corrupting every sample; a gate is what stops
   the next callset repeating it silently.
3. Promote `06a_af_concordance.{py,sh}` out of `scripts_archive/`. `06_ancestry_qc.sh` reads
   `exclude_af_concordance.txt`, nothing builds it, and the step continues **without** the
   cohort-artifact exclusion. Largest concrete gap in the QC path (`HANDOFF.md` issue #3).
4. Before step 7: fix the `ispd` bug (`07_gwas.sh:247`) — it misclassifies all 86 BR-DSNWGS
   rows and starts to matter the moment this callset reaches the merge.
5. Land the upstream GenoTools fixes in `dvitale199/GenoTools`: the positional column rename
   (`ancestry.py:247`), worker sizing, and making `get_common_snps` report its overlap.

**Also unresolved, from earlier entries:** the cross-callset sex-file inconsistency (the other
three callsets were sex-updated from `<dataset>/metadata/`, not from the corrected
`clinical_core_out/`); and the BR-DSNWGS PCs and ancestry labels already sitting in
`analysis_grain.csv` for a callset that had never cleared step 1.

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
