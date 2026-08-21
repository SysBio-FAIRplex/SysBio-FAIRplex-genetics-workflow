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
derives the root from its own location, and `clinical_common.py`, `05_excludelist.py` and
`ancestry_qc_manifest.py` each do the same. To move the project, copy the folder.

That collapses three things that used to be separate and disagree: `WGS_ROOT`, `PROJECT_ROOT`, and
the old `BUNDLE` (which pointed one directory level away from where the code actually was).

**Every handoff file has exactly one location.** The clinical side writes to `clinical_core_out/`
and the genotype steps read it from there — nothing is copied into `data/`, so nothing can go
stale against the run that produced it.

| | |
|---|---|
| `clinical_core.py` | goal 1, first half. Runs **once, before step 1**. Reads `data/`, writes `clinical_core_out/`. |
| `analysis_grain.py` | goal 1, second half — §11-13. Runs **once, after step 6**, on its PCs. |
| `clinical_common.py` | paths, readers, and the reconciliation rules both halves share. Imported, never run. |
| `wgs_core.ipynb` | goal 2 orchestrator — ssh/sbatch driver. §3 step 6, §4 the before/after proof. |
| `config.sh` | every path, derived from its own location. `submit.sh` — sbatch wrapper. |
| `scripts/00`–`08` | env, genotools, normalize, merge, relatedness, excludelist, ancestry QC, GWAS, ctrl-vs-ctrl. |
| `scripts/af_concordance_build.{py,sh}` | step 6 **stage B**, in-job. The `.sh` is for re-tuning knobs only. |
| `scripts/ancestry_qc_manifest.py` | step 6 stage E: `retained_samples_manifest.csv`, once per generation. |
| `review/plot_af_filter_effect.py` | the collaborator-facing before/after figure + eta² tables. |
| `scripts/gene_annot.py` + `ref/refFlat.txt` | **the single source of locus coordinates**, and now a **read-only CLI with no pipeline callers** — the sentinel API both pipeline scripts imported was deleted 2026-08-20. `python3 scripts/gene_annot.py CR1 SNCA LRRK2` / `--at chr12:40227079` |
| `scripts/diag_order.py` | read-only **preflight gate** — variant order vs the reference panel. Run on any callset not already cleared, BEFORE step 1. |
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

The clinical side runs at two points and step 6 runs **once**:

```bash
# BOTH lines. The clinical side lives at the project ROOT, not under scripts/ — an rsync of
# scripts/ alone silently leaves the cluster on the old single-file clinical_core.py.
# The host is the FQDN: there is no `helix` ssh alias, so a bare `helix:` fails to resolve.
# Transfers go through helix, never biowulf (scp to biowulf fails — PROJECT_LOG 2026-08-11).
ROOT=/data/CARDPB2/sysbio/wgs
rsync -av --exclude='logs/' --exclude='__pycache__/' scripts/ helix.nih.gov:$ROOT/scripts/
rsync -av clinical_common.py clinical_core.py analysis_grain.py config.sh helix.nih.gov:$ROOT/

cd /data/CARDPB2/sysbio/wgs && source config.sh

module load python/3.11 && source .venv/bin/activate   # NOT the system Anaconda py3.9
python3 clinical_core.py                  # §1-10  -> *_update_sex.txt
                                          # §12a   -> sample_annot.csv

# 0. VCF -> pgen. BR-DSNWGS only; the other three shipped as pgen or were built earlier.
#    NO SCRIPT — this ran as notebook cells, and these two commands are the whole record.
#    --sort-vars is not optional: see "Blocker 2" below.
bcftools view --apply-filters 'PASS,.' --min-alleles 2 --max-alleles 2 --type snps \
    $DIR_BR/joint_calls/AMPPD_postmortem_joint_gt_call_97donors.vcf.gz \
  | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' \
  | bgzip -@ 8 > $DIR_BR/pgen/intermediate/br_dsnwgs_filtered.vcf.gz
plink2 --vcf $DIR_BR/pgen/intermediate/br_dsnwgs_filtered.vcf.gz \
    --chr 1-22,X,Y --split-par hg38 --update-sex placeholder_sex.txt \
    --make-pgen --out $RAW_BR

# 1. genotools — once per callset. Reads the sex files in place from clinical_core_out/.
./submit.sh scripts/01_genotools.sh --job-name=genotools_wgs_harm \
  --export=PGEN=$RAW_WGS,SEX_FILE=$(sex_file wgs_harm),OUT_DIR=$DIR_WGS/genotools,DATASET=wgs_harm
./submit.sh scripts/01_genotools.sh --job-name=genotools_divco_hs \
  --export=PGEN=$RAW_DC,SEX_FILE=$(sex_file divco_hs),OUT_DIR=$DIR_DC/genotools,DATASET=divco_hs
./submit.sh scripts/01_genotools.sh --job-name=genotools_wb_dwgs \
  --export=PGEN=$RAW_WB,SEX_FILE=$(sex_file wb_dwgs),OUT_DIR=$DIR_WB/genotools,DATASET=wb_dwgs
./submit.sh scripts/01_genotools.sh --job-name=genotools_br_dsnwgs \
  --export=PGEN=$RAW_BR,SEX_FILE=$(sex_file br_dsnwgs),OUT_DIR=$DIR_BR/genotools,DATASET=br_dsnwgs

# 2. normalize — once per callset.
./submit.sh scripts/02_normalize.sh --job-name=norm_wgs --export=PFILE=$PF_WGS,OUT=$NORM_WGS,TAG=wgs
./submit.sh scripts/02_normalize.sh --job-name=norm_dc  --export=PFILE=$PF_DC,OUT=$NORM_DC,TAG=dc
./submit.sh scripts/02_normalize.sh --job-name=norm_wb  --export=PFILE=$PF_WB,OUT=$NORM_WB,TAG=wb
./submit.sh scripts/02_normalize.sh --job-name=norm_br  --export=PFILE=$PF_BR,OUT=$NORM_BR,TAG=br

./submit.sh scripts/03_merge.sh           # all four -> cohort_merged
./submit.sh scripts/04_relatedness.sh     # KING, report-only
./submit.sh scripts/05_excludelist.sh     # -> retained_manifest.csv
./submit.sh scripts/06_ancestry_qc.sh     # A unfiltered QC+PCA · B build AF list · C apply
                                          # D filtered QC+PCA   · E both manifests
python3 analysis_grain.py                 # §11-13 -> analysis_grain.csv on the new PCs

# the proof the filter works — reads both manifests step 6 wrote in that one job.
# Runs LOCALLY in practice (notebook §4 rsyncs the two manifests down first); it needs no
# genotypes, only the manifests, so either machine works.
python3 review/plot_af_filter_effect.py \
    --before $MERGED_DIR/by_ancestry_qc/unfiltered/retained_samples_manifest.csv \
    --after  $MERGED_DIR/by_ancestry_qc/retained_samples_manifest.csv

# Step 7 reads §13's pheno/covar/contrasts files, NOT $GRAIN (changed 2026-08-20, known issue 1).
# So analysis_grain.py above is a hard prerequisite, and a contrasts.csv predating 2026-08-21 lacks
# the callset-skew columns -> step 7 exits 3 naming the missing one. Rerun the grain, don't patch it.
./submit.sh scripts/07_gwas.sh            # reads $PHENO_SRC/$COVAR_SRC/$CONTRASTS_CSV
python3 review/plot_gwas.py               # QQ + Manhattan, local, on downloaded sumstats
./submit.sh scripts/08_ctrl_ctrl_filter.sh   # annotates; never subtracts (see known issue 2)
```

**Two clinical invocations, and unlike step 6's two passes this one is real.** Step 1 applies the
sex files; the grain carries step 6's PCs because step 7 reads them as covariates. There is no
column to split out. What changed is that they are now two *files* rather than one script run
twice — `analysis_grain.py` reads §9's audit tables instead of re-deriving them, so it no longer
rewrites the sex files step 1 already consumed. Section numbers are unchanged.

**Step 6 used to run twice, and the reason was wrong.** The claimed cycle was
`grain ← §12 ← manifest ← step 6`, so step 6 could not depend on the grain. But
`af_concordance_build` reads exactly `IID → (source_callset, dx_detailed)` — it never reads a PC,
and it never reads ancestry either, because the stratum comes from which fileset a sample is in.
Both fields are pure clinical output, so §12a now writes them as `sample_annot.csv` before step 1
runs. Separately, `--geno/--maf/--hwe` are per-variant on a fixed sample set and therefore
**commute with `--exclude`**, so pass 2 never needed to re-scan `cohort_merged` at all — stage C
applies the list to stage A's output instead. See the header of `scripts/06_ancestry_qc.sh`.

The stale-list guard is now structural rather than a check: stage B rebuilds the list inside the
job from this merge. The mtime refusal survives only on the `SKIP_AF_BUILD=1` path, which is the
one way a list this job did not build can still reach stage C.

## Status

**THE PIPELINE IS COMPLETE THROUGH STEP 8 (2026-08-21).** Steps 0–8 have all run on the four-callset
cohort. What is left is `review/plot_gwas.py` — local, on downloaded sumstats — and the write-up.

| | |
|---|---|
| `cohort_merged` | 172,497,055 variants × 13,334 samples (job 27429821) |
| retained | **12,495** of 13,334 (839 excluded: 499 relative, 319 duplicate, 21 sex) |
| exclusion list | **4,187** variants, HWE excess-over-chance gate ON (was 4,415 ungated) |
| `analysis_grain.csv` | 12,495 rows × 22 cols, **17 viable contrasts** — EUR 14, AJ 3 |
| step 7 GWAS | job **28004190** — 44/44 contrasts ran, λ_GC **1.0175–1.0549** |
| step 8 ctrl-vs-ctrl | job **28037485** — annotated; primary sumstats untouched |

Step 6 ran twice as a single pass: job **27857727** (7m13s, ungated) then a rerun with the gate on
(2026-08-20 21:37), which produced 4,187 and regenerated stages C/D/E. **That rerun's PCs and
manifests are the live ones.** *(Its job ID is not recorded here — recover it from
`bash scripts/runlog.sh --md`.)*

**λ_GC is well controlled everywhere**, 1.0175 (`AD_ampad_vs_control_ampad`) to 1.0549 (`PD_vs_AD`,
the primary and most confounded). No contrast shows inflation needing intervention.

**Read the two headline results in `PROJECT_LOG.md`'s "Where we are right now"**: the filter takes
EUR's callset eta² 0.757 → 0.036, and step 8's licence boundary is what kept APOE-ε4 (P=3.55e-15) in
the association set. Known issue 2 carries the second in full.

---

## The sentinel tripwire is DELETED — done and RUN 2026-08-20/21, both call sites

**Removed from `scripts/af_concordance_build.py`** (the import, `sentinel_detail()`, `keep_n()`, the
`---- sentinel loci ----` block, `arm_index`, and the per-arm `.keep` writes in `assoc_pair()` —
those duplicated the `.pheno` file's membership and existed only to feed `sentinel_detail`), **and
from `scripts/08_ctrl_ctrl_filter.py`.** `SENTINEL_GENES` / `SENTINEL_FLANK` / `SENTINEL_VARIANTS` /
`sentinel_hits()` are gone from `scripts/gene_annot.py` with them; `load_genes` / `gene_region` /
`Annotator` / `parse_variant_id` stay. All three files compile and the CLI is smoke-tested.

**The 08 sub-question was decided: remove both.** 08's report served the opposite purpose (a known
locus in the ctrl-vs-ctrl scan is the screening asymmetry, an argument *against* subtracting), so it
was not obviously the same decision. What settled it is the **negative** it could print — "no known
AD/PD locus among the flagged", off an uncited 9-gene list, reads as reassurance that the flags are
safe to subtract, and a flagged variant on a real locus outside those 9 produced identical output.
That is this project's rule 3 ("absence of a warning is not evidence") in the exact shape 08's own
docstring was written to fix, and resolving coordinates from refFlat had narrowed the hole without
closing it. Nothing was lost: the substantive warning 08 needed does not depend on a gene list and
is now printed unconditionally for every flagged variant.

Two smaller changes fell out, both to keep dead concepts from lingering (rule 2):
`gene_annot.sentinel_regions` → **`gene_regions`** (it was never sentinel-specific, and the name
outlived the set), and the `--at` CLI now parses through `parse_variant_id` instead of its own
inline parser — so it accepts full `chr12:40227079:C:T` IDs, errors cleanly on junk, and
`parse_variant_id` does not become dead code.

**Kept: `scripts/gene_annot.py` as the CLI for "what gene is this variant in."** This is the honest
per-variant replacement, and it is what exposed the mislabelling: `--at chr12:40227079` → `LRRK2`.

**Why the ±500 kb windows made it worse, recorded so it is not repeated.** The flank was sized for
step 8's question ("did my association signal land in a known locus?", where LD blocks matter — it
measured MAPT±500kb catching 2,286 of 2,318 hits). For 6a's question ("am I deleting a variant that
is part of a known locus?") it labels flank hits with the gene's name, and **none of the 21
"LRRK2±500kb" variants were in LRRK2** — they were in `SLC2A13` and `C12orf40`, 240–435 kb away.
Same for SNCA (intergenic), GBA1 (`DAP3`), APOE (`ZNF285`/`ZNF229`).

**All of this has since RUN**, and the removal is confirmed live: `step6_summary.txt` from the gated
rerun contains no sentinel output, and step 8's log (job 28037485) prints the unconditional
screening-asymmetry note in place of the deleted known-loci block. The `gene_annot.py --at` CLI is
what named APOE and HLA-DQB1 in step 8's results — the replacement path doing the job the 9-gene
list was supposed to.

### The MHC aggregate question — NOT a gate on step 7. Decided 2026-08-20.

Dropping the tripwire does not dispose of "is this filter disproportionately hitting known loci?",
and the 279-hit dump did surface something: **245 of the 279 were MHC** (`HLA-DRB1±500kb` 213,
`HLA-B±500kb` 32), roughly chr6:30.8–33.1 Mb.

**But the decision that question was meant to inform is already settled by the per-arm tables**
(`PROJECT_LOG.md` 2026-08-20, "the per-arm tables settle the liftover verdict per-variant"). Those
MHC flags have a per-variant technical mechanism: the two natively-called callsets agree near 0
while the lifted one does not — `chr6:32474706` `divco_hs` 0.008 / `wb_dwgs` 0.005 vs `wgs_harm`
0.207–0.213; `chr6:32588203` 0.000 / 0.000 vs 0.088–0.107. That is exactly the pattern that
licences 6a to delete (disease held constant, so a between-callset gap is technical). "Annotate the
MHC rather than subtract it" is therefore argued *against* on per-variant evidence, not left open.

So the flag rate is a **methods-reporting number**, not a decision input, and it does not gate step
7. Run it when writing up, to state the magnitude — "the filter removes X% of MHC variants against
Y% genome-wide" — not before:

```bash
cd /data/CARDPB2/sysbio/wgs && source config.sh
B=${MERGED_DIR}/by_ancestry_qc/unfiltered/cohort_EUR_qc.bim
L=${MERGED_DIR}/exclude_test_assoc.txt
tested_mhc=$(awk '{split($2,p,":"); if ((p[1]=="chr6"||p[1]=="6") && p[2]>=30800000 && p[2]<=33100000) n++} END {print n+0}' "$B")
tested_all=$(wc -l < "$B")
flag_mhc=$(awk -F: '($1=="chr6"||$1=="6") && $2>=30800000 && $2<=33100000 {n++} END {print n+0}' "$L")
flag_all=$(wc -l < "$L")
python3 -c "
tm,ta,fm,fa=$tested_mhc,$tested_all,$flag_mhc,$flag_all
print(f'MHC    : {fm:,}/{tm:,} = {100*fm/tm:.4f}%')
print(f'genome : {fa:,}/{ta:,} = {100*fa/ta:.4f}%')
print(f'enrichment: {(fm/tm)/(fa/ta):.1f}x')"
```

**Why it is still worth reporting, even though it decides nothing:** `06_ancestry_qc.sh`'s own
header argues the MHC is a real AD locus and that masking it from *association* would delete
signals we most expect to see — which is why the high-LD BED is applied to the PCA input only. This
filter reaches the association set, and `HLA-DRB1/DRB5` is one of this study's two real findings.
So the magnitude of what was removed from that region belongs in the methods as a stated limitation.
What it is *not* is an open question about whether to annotate instead of subtract there — the
per-arm evidence above answers that, and the answer is subtract.

---

**The BR-DSNWGS grain item is CLOSED.** The 95 retained BR samples are 76 EUR / 13 AJ / 3 AMR /
1 AAC / 1 AFR / 1 CAH. The old 19 AFR / 67 EUR rows are gone.

**AJ is closed on power and design, not on eta².** `AJ/AD` is 97 (`wgs_harm`=95 + `divco_hs`=2),
under the 100 floor, so AJ cannot field the primary AD-vs-PD contrast at all. Its three viable
contrasts are all `within_cohort` — ~96% `wb_dwgs` on both arms — and callset↔phenotype
collinearity cannot bias a contrast whose arms share a callset. The unsolved 0.962 has no
downstream consumer. The sub-continental-structure hypothesis is still untested and no longer
blocks anything.

| Step | state |
|---|---|
| 0 VCF→pgen | done — `br_dsnwgs_hg38.{pgen,psam}` — **but built by notebook cells, not a script** |
| 1 genotools | done all four — BR: 77 EUR / 14 AJ / 3 AMR / 1 CAH / 1 AFR / 1 AAC |
| 2 normalize | done all four — BR job 27429119 |
| 3 merge | done — 97/97 BR in, variant arithmetic closes exactly |
| 4 relatedness | done — 13,334/13,334 labelled, 11 strata |
| 5 excludelist | done — **12,495 retained** (839 excluded), 95/97 BR |
| 6 ancestry QC, unfiltered | done, job 27602590 — 6 of 11 strata have PCs. **This is now stage A's output** |
| `ancestry_qc_manifest.py` | done — 12,495-row `retained_samples_manifest.csv` |
| §11–13 (now `analysis_grain.py`) | done 2026-08-18 — grain rebuilt, 17 of 44 contrasts viable |
| `af_concordance_build` | done 2026-08-18, job 27697096 — **4,415 variants**, correct grain |
| premise test (eta² before/after) | done 2026-08-19 — **EUR 0.757 → 0.036; AJ 0.984 → 0.962** |
| step 6 rebuilt as ONE pass | done 2026-08-19, **RUN 2026-08-20** — job 27857727, 7m13s |
| clinical side split in two | done 2026-08-19, **RUN 2026-08-20** — §12a self-check 0 mismatches |
| 6 ancestry QC, single pass | **done** — job 27857727. Stage B rebuilt 4,415 (matches 27697096) |
| `analysis_grain.py` | **done 2026-08-20** — 12,495 rows × 22 cols, 17 viable, BR fixed |
| frequency test → `plink --assoc` | done + **VERIFIED IDENTICAL** 2026-08-20, job `af_swap_test2` |
| HWE excess-over-chance gate | **done + RUN 2026-08-20** — list is **4,187** (was 4,415); LRRK2 un-excluded, CR1 still excluded |
| **remove the sentinel entirely** | **done + RUN 2026-08-20** — both call sites + the `gene_annot` API; no sentinel output in `step6_summary.txt` |
| step 6 rerun on the gated list | **done 2026-08-20** — 4,187-variant list applied; stage D/E regenerated |
| grain on the gated PCs | **done 2026-08-20** — 12,495 × 22, 17 viable, unchanged shape |
| eta² re-measurement on 4,187 | **done 2026-08-20** — EUR 0.757 → 0.036 (**95.3%**), i.e. the gate cost the filter nothing |
| known issue 1 (pheno/covar) | **RESOLVED + RUN 2026-08-21** — step 7 reads §13's files; awk rebuild deleted |
| MHC flag-rate measurement | **not a gate — deferred to write-up 2026-08-20.** Per-arm tables already settled subtract-vs-annotate; the rate is a methods number |
| callset-skew columns (issue 10) | **added 2026-08-21, NOT YET RUN** — §13 emits them; needs `analysis_grain.py` before any step-7 rerun |
| 7 GWAS | **done + RUN 2026-08-21** — job 28004190, 44/44 ran, λ 1.0175–1.0549 |
| 8 ctrl-vs-ctrl | **done + RUN 2026-08-21** — job 28037485. APOE + HLA-DQB1 flagged = screening asymmetry, see known issue 2 |
| `review/plot_gwas.py` | **NEXT ACTION** — runs locally on downloaded sumstats |
| write-up | after that. `PROJECT_LOG.md`'s top block lists the two headline results |

**Running it from here costs less than the old pass 2 did.** Stage A's inputs (`cohort_merged`,
the step-5 manifest, the locked thresholds) have not changed, so job 27602590's output *is* stage
A's output. Notebook §3 moves it to `by_ancestry_qc/unfiltered/`, step 6's reuse rule picks it up,
and the ~13-14 min × 11 strata scan is skipped. What actually runs is stage B, the exclusion, and
the filtered prune+PCA.

**Stage B is the regression test for the split.** It rebuilds the list through `sample_annot.csv`
instead of the grain. Job 27697096 got 4,415 variants from the grain; the same count means the
refactor changed nothing. §12a also prints a direct self-check — every IID shared with an existing
grain must agree on `source_callset` and `dx_detailed`. A disagreement means the exclusion list
would move for a reason that has nothing to do with genotypes; stop there.

**The filter is warranted, and its effect is now MEASURED rather than cited.** Unfiltered eta² on
the four-callset cohort: AJ PC1 **0.984**, EUR PC2 **0.757**. Removing BR changes neither
(0.984 / 0.758), so this is real `wgs_harm`↔`wb_dwgs` structure. Applying the 4,415-variant list:

**Measured on the LIVE 4,187-variant gated list, 2026-08-20** (`results/pca/af_filter_effect_gated4187.csv`):

| stratum | n | max eta² unfiltered | max eta² filtered | reduction |
|---|---|---|---|---|
| EUR | 10,135 | 0.757 (PC2) | **0.036** (PC6) | **95.3%** |
| AJ | 1,518 | 0.984 (PC1) | **0.962** (PC1) | 2.2% |
| AAC | 226 | 0.877 (PC2) | 0.876 (PC2) | 0.1% |
| AFR | 191 | 0.708 (PC1) | 0.713 (PC1) | **−0.6%** |
| CAH | 106 | 0.309 (PC4) | 0.199 (PC4) | 35.6% |
| AMR | 185 | 0.203 (PC8) | 0.170 (PC8) | 16.0% |

**The HWE gate did not cost the filter anything.** On the ungated 4,415 list EUR was 0.036 / 95.2%;
on the gated 4,187 it is 0.036 / 95.3%, and AJ is identical to three decimals. That is what the
"the 228 withheld are general variant QC, not cohort-artifact removal" reading predicted, and it is
now measured rather than argued. **The 95% claim is safe to quote against the live list.**

AFR's reduction is slightly *negative*. Not a regression: AAC and AFR eta² is one BR sample each at
39.7σ / 19.6σ, not callset structure (see below), so these two rows are dominated by a single point
and move by noise. No AAC/AFR/CAH/AMR contrast is `viable_ge100`, so none of those four PCs reaches
a GWAS — the 17 viable contrasts are AJ (3) and EUR (14).

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

**~~`analysis_grain.csv`'s BR-DSNWGS rows are wrong~~ — CLOSED 2026-08-20, regenerated.** The old
grain carried 19 AFR / 67 EUR over 86 rows, which was not reconcilable with the post-`--sort-vars`
run. The cluster's grain is now 12,495 rows × 22 cols on step 6's filtered PCs, and its 95 retained
BR samples are 76 EUR / 13 AJ / 3 AMR / 1 AAC / 1 AFR / 1 CAH — the same item recorded as closed at
the top of this file. (Step 1's pre-exclusion spread is 77 EUR / 1 AFR over 97; the difference is
the 2 samples step 5 drops, not a disagreement.) **The laptop's copy is still the stale 11,918-row
one** and must never be rsynced upward — see Provenance.

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

1. **~~Two implementations of pheno/covar~~ — RESOLVED 2026-08-20. §13 is the only one.**
   `07_gwas.sh` now reads `clinical_core_out/{pheno,covar}/` and `contrasts.csv`; its awk rebuild of
   both files is deleted. §13's docstring already claimed step 7 read them ("the definition exists
   once and cannot drift between the two") — the claim is now true rather than aspirational.
   `config.sh` exports `PHENO_SRC` / `COVAR_SRC` / `CONTRASTS_CSV`.

   Three behaviour changes that follow, none of them accidental:
   - **`MIN_ARM` can only raise the floor.** §13 writes no file below 20 per arm, so `MIN_ARM=0` no
     longer means "force literally all". Lower it in §13 if that is ever wanted.
   - **`EXCLUDE_DUAL` is gone from step 7** and now lives only in §13, which is what §13 argued for:
     a sensitivity run should be a recorded artifact, not a flag remembered at GWAS time. Step 7
     **refuses** the old `--export=EXCLUDE_DUAL=1` rather than ignoring it. Run it by setting
     `EXCLUDE_DUAL=True` in `analysis_grain.py`, rerunning to a separate `CLINICAL_OUT`, then
     pointing `PHENO_SRC`/`COVAR_SRC`/`CONTRASTS_CSV` at that generation.
   - **The age decision moved to §13**, which emits an `AGE` column or does not. Step 7 reports what
     the covar file actually contains instead of sniffing the grain header (issue 9).

   **Two stale-generation guards replaced a header comment.** The old script only warned, in a
   comment, that re-running step 6 invalidates the grain and "this silently runs on stale
   covariates". Now: each `covar_<ANC>.txt` must be newer than its `cohort_<ANC>_pca.eigenvec`
   (exit 4), and every IID in a pheno file must appear in the `.fam` being tested (exit 5) — which
   catches a mismatched generation that mtime alone would miss. `contrasts.csv` is read **by header
   name**, so a reordered column aborts with the missing name instead of mislabelling a contrast.

2. **~~Two implementations of the ctrl-vs-ctrl artifact filter~~ — RESOLVED 2026-08-21.
   `scripts/08_ctrl_ctrl_filter.{py,sh}` is the only one.** `review/mask_cohort_artifacts.py`, the
   older post-hoc version, was deleted: it lacked the licence argument below, and its docstring
   still cited `gwas_per_ancestry.sh` (retired) and a `docs/METHODS.md` that never existed. Step 8
   carries the argument that APOE is *expected* to reach significance in the control-vs-control
   scan — the two control arms are differentially screened, AMP-AD's as cognitively normal and
   AMP-PD's for PD and not AD — so subtracting destructively would delete the study's strongest
   true locus. It therefore annotates (`CTRL_P` in `.ccannot.tsv`) and writes `.ccfilt.tsv`
   separately, never touching the primary.

   This is also the licence boundary against step 6a: **6a may delete variants** because disease is
   held constant inside each of its cells; **step 8 may not**, because its control definitions
   differ across programs. Not duplicated reasoning — different entitlements.

   **CONFIRMED EMPIRICALLY 2026-08-21, job 28037485 — this is no longer an argument, it is a
   measurement.** All 4 flagged genome-wide hits in `EUR AD_ampad_vs_control_ampad` (the
   *within-cohort* AD contrast, Δ=0.0, the cleanest one in the study) are APOE, including both
   causal variants: `chr19:44908684:T:C` = **rs429358/ε4 at P=3.55e-15**, `chr19:44908822:C:T` =
   rs7412/ε2, plus `chr19:44906745:G:A` inside the gene and `chr19:44912456:G:A` ~3 kb past its 3′
   end. **`.ccfilt.tsv` would delete ε4 at 3.55e-15 from the study's cleanest AD contrast — so do
   not use it for the AD contrasts.** Report the primary; name APOE as expected-and-retained.

   Separately, 2 of the 3 genome-wide ctrl-vs-ctrl hits are `HLA-DQB1` (chr6:32,661,554/32,661,570),
   ~30 kb from `HLA-DRB1/DRB5` — the same screening asymmetry at the study's other real finding. The
   third, `chr3:106666502`, is intergenic 443 kb from any gene and is the one plausible genuine
   artifact. **2 of 3 are real biology.**

   **The corroboration that 6a worked:** every high-signal contrast came through with zero
   control-flagging — `AD_vs_PSP` 2,319 hits / 0 flagged, `PD_vs_control` 2,523 / 0,
   `PD_amppd_vs_control_amppd` 2,508 / 0. Step 8 is confirming, not correcting.

   **Caveat to report rather than bury:** `PD_vs_MCI` is 3 hits, 3 flagged, 0 kept — its
   `.ccfilt.tsv` is empty (n_ctrl = 159, thin to begin with).

   **Methods sentence worth writing:** it is the demo's default `--ctrl-p 1e-5` that catches APOE.
   All four APOE CTRL_P values are 1.3e-07 – 6.2e-06, i.e. *above* 5e-8, so APOE is not among the 3
   genome-wide ctrl-vs-ctrl hits. At `--ctrl-p 5e-8` the AD contrasts would lose nothing. Left at
   1e-5 deliberately — flagging is annotation, and a wider net with correct labels is safer — but
   the threshold, not the mechanism, is what put APOE in the flagged set.

3. **~~Step 6 runs TWICE~~ — RESOLVED 2026-08-19. It is one pass.** The two-pass shape was
   justified by a cycle that did not exist, and the write-up is in `scripts/06_ancestry_qc.sh`.
   Short version: the AF build never needed the grain (only `IID → (source_callset, dx_detailed)`,
   now written by §12a as `sample_annot.csv` before step 1 runs), and `--geno/--maf/--hwe` commute
   with `--exclude`, so pass 2's re-scan of `cohort_merged` was recomputing what it already had.

   ```
   6  A  unfiltered per-ancestry QC + PCA   -> by_ancestry_qc/unfiltered/   (the BASELINE, kept)
      B  af_concordance_build              -> exclude_af_concordance.txt   (from A, in-job)
         read the BY CALLSET PAIR table and the per-cell HWE ratio table in its log
         (the sentinel tripwire that used to print here was deleted 2026-08-20)
      C  apply the exclusion               -> cohort_<ANC>_qc              (the ASSOCIATION set)
      D  prune + PCA on the filtered set   -> cohort_<ANC>_pca             (the COVARIATES)
      E  both retained_samples_manifest.csv files
      analysis_grain.py §12                -> PCs changed, so the grain must be rebuilt
   ```

   Two things the old shape needed and this one does not: the `mv by_ancestry_qc` before pass 2
   (the baseline is now a permanent named output, not something you must remember to preserve),
   and the mtime refusal on a stale exclusion list (stage B rebuilds it in-job from this merge; the
   guard survives only on the `SKIP_AF_BUILD=1` path).

   Why the filter exists — **measured on THIS cohort, 2026-08-19.** Applying the 4,415-variant list
   takes EUR's worst callset eta² from **0.757 to 0.036 (95.2%)**, by excluding 0.058% of its
   variants. The effect is concentrated, not diffuse. That measurement is no longer a one-off:
   stage D produces PCs for both generations every run, and `review/plot_af_filter_effect.py`
   renders the before/after — figure and eta² tables — from the pair. See `PROJECT_LOG.md`
   2026-08-19. The eta² CSVs are now versioned (`.gitignore` carries the reason).

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

4. **The clinical side runs on the CLUSTER ONLY — both halves. Do not run either on a laptop.** The two machines
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

   Run them as: `module load python/3.11 && source .venv/bin/activate && python3 clinical_core.py`
   (before step 1) and `python3 analysis_grain.py` (after step 6). The system Anaconda py3.9 is on
   `PATH` and is not the pinned environment.

5. **The cluster's `README.md` (30 KB) still has not been merged in.** The repo's `README.md` was
   brought up to date 2026-08-19 — four callsets, the current script list, the real handoff
   locations, and the gotchas learned since August — but it was rewritten from the repo's own
   state, not merged with the cluster's longer copy. `06_ancestry_qc.sh` references that copy's §2
   (reference data acquisition) by number, and §2 is the part the repo version still lacks:
   nothing here says how to obtain the reference panel or the liftover chain. Diff the two before
   trusting either.

6. **DivCo's source VCF is 0 bytes** on the cluster (`merged.deduped.vcf.gz`). The pgen was derived
   before it was truncated, so nothing is blocked, but DivCo cannot be re-derived from source
   without re-pulling from Synapse.

7. **~~The HWE stage has no excess-over-chance gate~~ — FIXED AND RUN 2026-08-20. CLOSED.** A
   (stratum × callset) cell now contributes exclusions only if its rejection count exceeds the
   number expected by chance at `--hwe`. No multiplier: E/O is the Benjamini-Hochberg FDR estimate
   for that cell's rejections, so at or below 1.0× there is nothing to attribute.
   `HWE_REQUIRE_EXCESS=0` restores the old unconditional union, which is what reproduces the
   4,415-variant list.

   **Result: the list is 4,187 variants** (from 4,415 — 228 removed). Both below-chance cells were
   `WITHHELD` and `EUR`/`wgs_harm` voted, exactly as designed; the provenance sidecar records
   `require-excess : True`. Both settled predictions held: `chr12:40227079:C:T` (LRRK2) is **gone**
   from the list, and CR1 `chr1:207521012:T:C` **remains** excluded on the frequency test plus its
   64.5% `divco_hs` call rate. The 228 is one short of the "up to 229" estimate below, so 228 of
   those variants were unique to the two withheld cells and 1 was independently flagged by another
   channel. `step6_summary.txt` contains no sentinel output, confirming that removal is live too.

   The evidence that motivated the gate, from job 27857727:

   | stratum | callset | controls | fail | exp by chance | ratio |
   |---|---|---|---|---|---|
   | EUR | `wgs_harm` | 328 | 1,680 | 377 | **4.5×** |
   | EUR | `wb_dwgs` | 3,064 | 132 | 377 | 0.35× |
   | AJ | `wb_dwgs` | 638 | 97 | 400 | 0.24× |

   `wgs_harm` clears expectation 4.5× from the *smallest* sample of the three — real het excess,
   consistent with mismapping in the lifted callset. The other two are at or below chance, and
   underdispersion is expected at these control counts (the HWE p-distribution is discrete and
   conservative at n in the hundreds), which is a second reason a below-chance row should not vote.
   Those two rows still contributed up to 229 of the 4,415.

   **What it costs concretely:** `chr12:40227079:C:T`, inside LRRK2, was never flagged by the
   frequency test (max spread 0.034, under the 0.05 threshold) and is excluded from every stratum
   on the strength of one chance-level HWE hit in AJ — a stratum that contributes nothing to the
   primary contrast. Gating on observed > 2× expected keeps EUR/`wgs_harm` and drops the other two.

   **SETTLED 2026-08-20 — the gate is fully determined.** CR1 **passes** HWE in all three testable
   callsets (the 2026-08-19 entry claiming it failed all three was misread; see the later log entry)
   and is excluded on the frequency test plus a **64.5% call rate in `divco_hs`** — 156 of 242
   alleles, where the same samples are 242/242 at the LRRK2 site, so the dropout is site-specific.
   CR1 is unaffected by this gate. LRRK2 fails **AJ/`wb_dwgs` only**, and the same callset tested in
   EUR with 3,064 controls instead of 638 **passes** — a real het-excess mechanism would show up
   more strongly in the larger sample, not vanish. So the gate leaves CR1 excluded and un-excludes
   LRRK2. Build it, rerun step 6 (~7 min) and the grain (~1 min).

   Related composition point for the methods: the 1,069 HWE-only additions are general variant QC,
   not cohort-artifact removal, and they ride into the same union. The docstring's "mechanism-based
   confirmation" framing is true of the 840 overlapping variants and not of the other 1,069.

   **Still open on this, deliberately not built:** moving the HWE channel to a post-hoc annotation
   on sumstats rather than a pre-association deletion. For the association set the two are
   identical — `--glm` tests variants independently and 5e-8 is a fixed convention, so striking a
   variant from the sumstats is the same as never testing it — and `scripts/08_*` already
   annotates-never-subtracts for exactly this reason. **The PCA input is the part that cannot move:
   PCs are covariates in every test, so the frequency channel must stay pre-applied regardless.**

8. **`plink1.9` is now a hard dependency of step 6 stage B.** As of 2026-08-20 the frequency test is
   `plink --assoc` (the 1-df allelic chi-square) rather than a hand-rolled two-sample test of
   proportions — verified identical before the swap: 1,698 vs 1,698 on EUR/AD with a symmetric
   difference of 0 both ways. plink2 dropped `--assoc` in favour of `--glm`, which is logistic
   regression on dosage and only asymptotically equivalent, so `MOD_PLINK1` is now loaded
   **unconditionally** in `06_ancestry_qc.sh` and `af_concordance_build.sh` where it previously
   loaded only for the optional mishap stage. Note `07_gwas.sh:49` records a deliberate preference
   *not* to rely on plink1.9 being present; step 7 still doesn't, but step 6 stage B now does. If
   that module ever goes away, the fallback is `--glm` plus a fresh equivalence check.

9. **No age covariate, and it cannot simply be added.** Every run prints `age covariate: NOT FOUND
   in grain — running WITHOUT age`. This is by design, not a regression: AMP-AD supplies age at
   **death**, AMP-PD age at **baseline/analysis**, and they are not the same variable, so forcing
   one column would silently mix two quantities. Both `analysis_grain.py` §13 and `07_gwas.sh`
   detect an age column by header name (`age`, `age_analysis`, `agedeath`, `age_death`,
   `age_baseline`, `age_cov`) and warn when absent.

   **Age is the dominant confounder for both AD and PD, so state this as a limitation in the
   methods.** AD-by-neuropathology donors are older at collection than the AMP-PD arms, and the
   primary contrast is exactly the one where the two age variables are least comparable — the
   `within_cohort` contrasts are the only ones where an age term would mean the same thing on both
   arms. Fixing it properly needs the phenotype track to emit one harmonized age, which nothing
   currently does. Recorded here 2026-08-20 because until then it lived only in two code comments.

   Also a third instance of known issue 1: that six-name detection is implemented independently in
   both places.

10. **`confound_tag` measures PROGRAM, not callset — and `within_cohort` is not a clean bill of
    health.** `AMPPD_CALLSETS = {"wb_dwgs", "br_dsnwgs"}`, so `delta_amppd` pools the two AMP-PD
    callsets. A contrast can therefore score `delta_amppd = 0.0`, be tagged `within_cohort` — the
    label `07_gwas.sh`'s header calls "the confound-free trusted backbone" — and still have one arm
    carrying a callset the other lacks entirely.

    **Measured on job 28004190 (2026-08-21).** `EUR PD_vs_DLB` is `within_cohort`, Δ=0.0, and step
    7's differential-missingness filter removed **353,068** variants from it — roughly 6× any
    genuinely cross-program contrast (`AD_vs_DLB` 65,245, `AD_vs_control` 63,393,
    `control_amppd_vs_control_ampad` 62,344). Cause: BR-DSNWGS's 95 retained samples are 71 PD /
    21 control / 3 null, so it contributes **no DLB at all**, and BR sits at ~50% missingness on the
    common set. ~55 EUR genomes generated more technical asymmetry than the entire AMP-AD/AMP-PD
    split. `PD_amppd_vs_control_amppd` has BR on both arms (71 v 21) and removed only 5,814,
    which confirms the mechanism rather than assuming it.

    **The results are protected; the label is what misleads.** The diffmiss filter is pre-association
    and caught it — that contrast came back at λ=1.0402 with 318 hits. But anyone reading the summary
    would rank `PD_vs_DLB` as cleaner than `AD_vs_DLB`, which is backwards.

    **What was added, and what was deliberately NOT.** §13 now emits `max_callset_delta`,
    `worst_callset` and `callset_one_sided` into `contrasts.csv`, and step 7 carries all three into
    `gwas_summary.csv` and annotates `within_cohort` rows. `confound_tag` itself is unchanged — it
    means "program", and redefining it would reinterpret every row already on record.

    **Two rejected designs, recorded so they are not retried.** (1) *A percentage-point threshold.*
    Reusing the tag's own 20pp yardstick looked tidy and fails: BR is ~55 of a 2,595-sample arm, so
    its share delta is **2.1pp** and a 20pp rule misses the one case the check exists for. (2)
    *One-sidedness as a warning.* It fires on 5 of 5 real contrasts, which is the objection that
    retired the sentinel tripwire. One-sidedness is not the predictor — one-sidedness of a **sparse**
    callset is. `PD_vs_control` is one-sided (`wgs_harm` 0 v 407) and removed only 6,331. BR is
    sparse because a 97-donor joint call emits nothing at sites monomorphic in its own donors, and
    §13 has no missingness data, so it cannot compute sparsity. **So `n_diffmiss_excluded` is the
    signal and `callset_one_sided` only names which callset explains it.** Step 7 ends by ranking
    `within_cohort` rows by that count — threshold-free, because a two-order-of-magnitude outlier is
    obvious side by side.

## Provenance

`scripts_archive/ad-pd-gwas-jul28/` on the cluster holds the July version of every script, plus 15
files that were never carried over: the AF-concordance pair, the ctrl-vs-ctrl step, and the
`diag_*` / `probe_*` diagnostics. Of the 12 scripts shared with this repo, 8 are byte-identical;
`config.sh`, `06_ancestry_qc.sh`, `07_gwas.sh` and `05_excludelist.py` differ, with this copy
taken as authoritative.

**Committed as of 2026-08-19 (`1277fc4`) — but there is still no remote.** That commit is the first
to carry `PROJECT_LOG.md`, `af_concordance_build.{py,sh}`, `08_ctrl_ctrl_filter.{py,sh}` and the
`diag_*` scripts; none of them had a second copy before it. The work now survives an accidental
`rm` or a bad checkout, and does **not** survive disk or laptop loss. One machine, one disk, with
history.

`clinical_core_out/` and `results/` remain gitignored, deliberately. One consequence still to know:
the laptop's `analysis_grain.csv` is the stale 11,918-row one (BR-DSNWGS as 19 AFR) and must never
be rsynced upward over the cluster's correct 12,495-row grain.

**`results/pca/*_eta2.csv` is no longer unversioned** — as of 2026-08-19 those tables and
`af_filter_effect*.csv` are the one exception in `.gitignore`, with the reason written there. They
are aggregate (one row per stratum, no IIDs, no genotypes), so they are safe to version, and the
0.984 / 0.757 baseline previously survived being overwritten only because the numbers had been
typed into `PROJECT_LOG.md` by hand. "Regenerable from code + data" was true in principle and
false in practice: regenerating a *baseline* means re-running the step in a state that no longer
exists. The negations are narrow on purpose — `results/retained_samples_manifest.csv` is one row
per sample and stays out.
