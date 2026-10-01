# Working rules for this project

Code is written in a development clone; it is executed **only on NIH biowulf**. The data is controlled-access —
AMP-AD (Synapse) and AMP-PD (GCP) — and both data use agreements confine it to the cluster.
`README.md` is the operator's guide; `METHODS.md` is the rationale.

## 1. The data boundary — above every other rule

- **Only code exists outside biowulf.** Controlled data, and anything derived from it at the
  sample level, exists only on biowulf and never leaves it, in either direction: `data/`,
  `clinical_core_out/`, `results/` (except the tracked `results/pca/*.csv`), sumstats, manifests,
  psams, eigenvecs, and step logs.
- **Nothing that reads project data runs off the cluster.** In the development clone you edit
  code and run git and static checks on code (syntax parsing, `scripts/nb_guard.py`) — nothing
  else. If a data file ever turns up there, do not open it; stop and tell the user.
- **Every run is a command you hand to the user**, who runs it on biowulf and pastes back the
  output. Do not ssh.
- **What comes back must be aggregate.** Commands you hand over print counts, summary tables,
  λ, file paths — never IIDs or donor IDs, per-sample rows, genotypes, or any count or frequency
  from a cell of fewer than 20 samples. A diagnostic that needs subject-level rows writes them to
  a file on the cluster and prints how many. Whatever is pasted into this session leaves the
  cluster.
- **Off the cluster: aggregate text and figures only. Onto the cluster: code, through GitHub.**
- **No participant IDs in code, comments, docs or commit messages** — not even as worked examples.
  Describe ID *formats* with placeholders (`R<7d>`, `PM-<site>_<n>`). `scripts/id_guard.py` enforces
  the known shapes at commit time; purely numeric IDs it cannot see.

## 2. Getting code to the cluster

Hosts are FQDNs: `helix.nih.gov` (transfers), `biowulf.nih.gov` (compute). Code is pushed to
GitHub from the development clone and pulled on the cluster. The repo is public, so the cluster reads it with
**no credentials, and stays read-only** — nothing on biowulf can push, so nothing there can leave
by git:

```bash
# biowulf login node (or helix), in the project root
git pull --ff-only
```

Code moves only through git; nothing else copies it. `09_amppd_subset.sh` stamps `git rev-parse HEAD` into every release,
so the cluster checkout must be the commit that ran.

## 3. Writing code here

- Ask before adding behavior to `scripts/*.sh` or the pipeline `.py`. Docs and read-only
  diagnostics are fine unprompted.
- Python on the cluster is `module load python/3.11 && source .venv/bin/activate`; the system
  Anaconda py3.9 on `PATH` is not the pinned environment.
- **Read columns by header name, never by position.** `plink2 --freq` puts `PROVISIONAL_REF?` at
  column 5.
- **When a check cannot run, it must say so loudly.** Never let "none found" print when the thing
  that would have found it was absent. A check whose negative result is not evidence should be
  deleted, not annotated.
- **A replacement is not done until the original is deleted**, with every call site converted in
  the same change. One implementation per concept.
- **Only step 6a may delete variants; step 8 annotates** (`METHODS.md` §6.2, §8).
