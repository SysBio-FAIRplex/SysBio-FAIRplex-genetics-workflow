# Working rules for this project

Read `HANDOFF.md` first (what is true now), then the "Where we are right now" block at the top of
`PROJECT_LOG.md` (what we did and what was ruled out). Where either disagrees with a prompt, they
win. `RUNLOG.md` is generated — `bash scripts/runlog.sh --md > RUNLOG.md`.

## 1. Search the log before investigating anything

`PROJECT_LOG.md` is append-only and holds every hypothesis that was ruled out. Before running a
diagnostic on an anomaly — a sentinel hit, an eta² value, a count that looks off — grep it:

```bash
grep -n -i '<locus|number|symptom>' PROJECT_LOG.md
```

**This is not a nicety, it is the most expensive failure mode in the project.** On 2026-08-20 both
sentinel-loci hits were re-derived from scratch and the conclusion was *weaker* than the entry
already sitting at `PROJECT_LOG.md` 2026-08-19: that entry had the `divco_hs` call rate (156 of
242 alleles — dropout, not a frequency difference) and the cross-callset HWE pattern. Re-deriving
also produced a contradiction with it that then had to be resolved separately. A finding that
disagrees with the log is worth reporting; a finding that *rediscovers* the log is waste.

Record every resolved anomaly back into `PROJECT_LOG.md`, including dead ends. That is what makes
the grep work next time.

## 2. A replacement is not done until the thing it replaces is deleted

`scripts/gene_annot.py` was written to replace a hardcoded `SENTINEL_LOCI` coordinate table, and
its docstring measured exactly how wrong that table was (37% of CR1, 52% of LRRK2, no HLA gene at
all). It was then never imported. Both copies of the bad table stayed in place and kept running,
so the tripwire under-reported for weeks with the fix sitting in the same directory. Wired in
2026-08-20 — and then, once wired, its first correct run showed the whole check was unsound and
**the sentinel was deleted outright the same day** (both call sites, plus the `sentinel_*` API in
`gene_annot.py`). Worth keeping in mind when applying this rule: wiring the replacement is what
made the concept testable, and the test killed it. `gene_annot.py` survives as a read-only CLI with
no pipeline callers.

So: **write the replacement, wire every call site, and delete the original in the same change.**
Dead-code fixes read as solved and are not. If a call site cannot be converted yet, say so in
`HANDOFF.md` "Known issues" with the reason — an unwired replacement that nobody flagged is
indistinguishable from a bug.

Corollary — one implementation per concept. Both known duplications are now closed: pheno/covar
2026-08-20 (`analysis_grain.py` §13 is the sole definition of who is a case, and `07_gwas.sh` reads
its files instead of rebuilding them in awk), and the ctrl-vs-ctrl filter 2026-08-21
(`review/mask_cohort_artifacts.py` deleted; `scripts/08_*` is the only implementation). Do not add
a third.

## 3. Absence of a warning is not evidence

Several bugs here exited 0 and looked like clean results: the CAH ancestry bug scored 0.97 model
accuracy while corrupting every sample; a stale exclusion list was applied silently; a "file count
matched" check passed while two files were missing. When a check cannot run, it must say so
loudly — never let "none found" be printable when the thing that would have found it was absent.

The sentinel tripwire used to print a banner for exactly this reason, and on 2026-08-20 it was
deleted instead — which is the stronger form of the same rule. Its refFlat-failure banner covered
the case where the check *could not run*, but nothing covered the case where the check ran fine and
its 9-gene scope was too narrow to mean anything: "no known AD/PD locus among the flagged" was
printable off an uncited list, and read as reassurance. **A check whose negative result is not
evidence should be deleted, not annotated with a caveat** — a warning nobody can act on is the same
failure wearing a different label.

## 4. Verify by name, not by position

Read columns by header name. `plink2 --freq` emits `PROVISIONAL_REF?` as column 5, so a positional
`$5`/`$6` read of an `.afreq` mislabels the frequency as the allele count — that happened on
2026-08-20 in an ad-hoc command and briefly hid the CR1 call rate. `read_afreq` and `read_annot`
both index by header name; ad-hoc `awk` should too.

Related: a file *count* is not a file *list*. That check was run against the clinical metadata
directories, passed, and was wrong.

## 5. Cluster etiquette

- **Hand cluster commands to the user; do not ssh.** They run them and paste back.
- Hosts are FQDNs: `helix.nih.gov` (transfers), `biowulf.nih.gov` (compute). No ssh alias exists.
- **Code reaches the cluster by `git`, never `rsync`.** The remote
  (`SysBio-FAIRplex/amp-ad-pd-wgs-gwas`, private) was created 2026-08-21. `rsync` without
  `--delete` cannot express a deletion, so retired scripts had to be removed by hand and three
  stale-artifact bugs came of it; git expresses deletions. **Local absence still ≠ cluster
  absence** for anything predating the checkout, and for `results/` and `clinical_core_out/`,
  which stay gitignored.

  **This rule said "retired the two-rsync dance" for four weeks while the cluster was not a
  git repository at all** (`git pull` → `fatal: not a git repository`). Converted in place
  2026-09-17; four files had drifted. A rule that describes an intention in the past tense reads
  exactly like a rule that describes a fact — see `PROJECT_LOG.md` 2026-09-17.

  **There are no GitHub credentials on biowulf** and GitHub refuses password auth, so `git fetch`
  against the remote fails there. The transport is a **bundle**: `git bundle create ~/x.bundle main`
  on the laptop, `scp` it to helix, then `git fetch ~/x.bundle main:refs/remotes/origin/main` on
  the cluster. Unlike loose-file `scp` this carries the real commit graph, so `git rev-parse HEAD`
  on the cluster returns the true commit — which `09_amppd_release.sh` stamps into every release.
  A PAT or SSH key on biowulf would retire the bundle step; neither exists yet.
- Never copy `clinical_core_out/` upward by any means — the laptop's `analysis_grain.csv` is stale
  and wrong (11,918 rows against the cluster's correct 12,495).
- Both clinical scripts (`clinical_core.py`, `analysis_grain.py`) are **cluster-only**; the laptop
  holds different clinical inputs and silently produces a smaller, wrong table.
- Python on the cluster is `module load python/3.11 && source .venv/bin/activate`. The system
  Anaconda py3.9 is on `PATH` and is not the pinned environment.
- Ask before adding behavior to `scripts/*.sh` (or the pipeline `.py`). Diagnostics and docs are
  fine unprompted.

## 6. Deleting variants is a licensed act

Only step 6a (`af_concordance_build.py`) may delete, and only because disease is held constant
inside each (stratum × dx) cell, so a between-callset frequency gap has to be technical. Step 8
may not delete — its control arms are differentially screened across programs, so APOE is
*expected* to reach significance there and subtracting would remove the study's strongest true
locus. It annotates instead.

**There is no automated check on 6a's licence, deliberately — as of 2026-08-20 you are it.** A
sentinel-loci tripwire held that role and was deleted: it gated nothing (stage C applies the list
later in the same job), and its 9-gene scope made both its positive and its negative results
uninformative. What replaces it is reading stage B's log — the BY CALLSET PAIR table and the
per-cell HWE ratio table — and resolving anything that looks wrong against the per-cell
`.afreq`/`.snplist` intermediates in `af_concordance/`, indexing by **header name** (rule 4: column
5 of an `.afreq` is `PROVISIONAL_REF?`). Call rate is usually the tell, not frequency: CR1 was
resolved by `divco_hs` calling it in 156 of 242 alleles while the same samples were 242/242 at the
LRRK2 site. To name the gene a variant sits in, use `python3 scripts/gene_annot.py --at <chr:pos>`.

Note the timing: step 6 applies the list at stage C, later in the same job that builds it at stage
B. So nothing about stage B can gate — resolve anomalies before step 7 reads the association set.
