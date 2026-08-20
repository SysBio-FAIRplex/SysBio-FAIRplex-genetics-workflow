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
2026-08-20.

So: **write the replacement, wire every call site, and delete the original in the same change.**
Dead-code fixes read as solved and are not. If a call site cannot be converted yet, say so in
`HANDOFF.md` "Known issues" with the reason — an unwired replacement that nobody flagged is
indistinguishable from a bug.

Corollary — one implementation per concept. Known live duplications, each already in
`HANDOFF.md`: pheno/covar (`analysis_grain.py` §13 vs `07_gwas.sh`'s awk — the awk one runs), and
the ctrl-vs-ctrl filter (`scripts/08_*` is authoritative over `review/mask_cohort_artifacts.py`).
Do not add a third.

## 3. Absence of a warning is not evidence

Several bugs here exited 0 and looked like clean results: the CAH ancestry bug scored 0.97 model
accuracy while corrupting every sample; a stale exclusion list was applied silently; a "file count
matched" check passed while two files were missing. When a check cannot run, it must say so
loudly — never let "none found" be printable when the thing that would have found it was absent.
Both sentinel call sites now print a banner instead.

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
- Code reaches the cluster by `rsync`, not `git pull` — there is no remote. **Two rsyncs**: the
  clinical files live at the project root, not under `scripts/`.
- `rsync` without `--delete` cannot express a deletion. Retired scripts must be removed by hand;
  that has caused three stale-artifact bugs. **Local absence ≠ cluster absence.**
- Never rsync `clinical_core_out/` upward — the laptop's `analysis_grain.csv` is stale and wrong.
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
locus. It annotates instead. The sentinel tripwire is the check on 6a's licence; a hit means stop
and resolve it against the per-cell `.afreq`/`.snplist` intermediates.

Note the timing: step 6 applies the list at stage C, later in the same job that builds it at stage
B. The tripwire therefore reports rather than gates — resolve hits before step 7 reads the
association set.
