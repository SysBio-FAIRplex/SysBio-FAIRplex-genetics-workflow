#!/usr/bin/env python3
"""STEP 5 — build the cross-dataset exclude list from the ancestry-split KING pairs.

Consumes the report-only relatedness outputs from 04_relatedness.sh and makes the picks ourselves.
Phenotype reconciliation of the 846 duplicates happens separately, on the phenotype side — this
script is genotype-side only (call-rate-driven).

Logic:
  1. DUPLICATES  (KINSHIP >= --dup-cutoff, default 0.354): union-find the dup pairs into clusters
     (a donor genotyped in >2 callsets forms one cluster); within each cluster KEEP the highest
     common-set call rate, DROP the rest.                         -> reason 'duplicate'
     (This is where the ~748 Rush<->ROSMAP different-IID dups finally get caught.)
  2. RELATIVES   (--rel-cutoff < KINSHIP < --dup-cutoff, default 0.0884..0.354) among still-retained
     samples: greedy maximal-unrelated-set — repeatedly drop the highest-degree node, tie-broken by
     LOWEST call rate (so we keep the better-genotyped sample).   -> reason 'relative_2nd_deg'
  3. GENOTOOLS QC-FAILS — recovered from each callset's genotools JSON `<ANC>_pass_fail` manifest
     (see README §3 step 5). Genotools QC is a per-ancestry STEP CHAIN; the
     JSON records, per ancestry, one {status, input, output} entry per step, where input/output are
     the bed-fileset STEMS bracketing that step:
         callrate: input=..._ancestry_<ANC>       -> output=..._<ANC>_callrate
         sex:      input=..._<ANC>_callrate        -> output=..._<ANC>_callrate_sex
         het:      input=..._<ANC>_callrate_sex    -> output=..._<ANC>_callrate_sex_het   (status=False
                   for small ancestries where het is SKIPPED -> that output fileset does not exist)
         related:  ...  -> output=..._<ANC>  (FINAL)   [IGNORED — see below]
     Each step's fail set = (IIDs in `input` fileset) MINUS (IIDs in `output` fileset), read from the
     .fam (genotools ran on --bfile; .psam fallback). het is computed ONLY when status is True.
     related-pruned is IGNORED on purpose: step 4 owns relatedness cross-dataset, so the QC-fail
     union is callrate ∪ sex ∪ het ONLY — no double-pruning.        -> reason 'genotools_qc_fail:<step>'
     WHY the JSON manifest and not `pruned_samples`: probing showed `pruned_samples` is incomplete
     (never lists callrate fails; absent entirely for WB-DWGS). `<ANC>_pass_fail` is uniform across all
     3 callsets and its input/output stems give the authoritative fileset paths (no globbing).
     Pass the 3 callset JSONs via --genotools-json (repeatable; defaults below). QC-fail IIDs are
     matched to the common-set sample list by IID (genotools .fam IIDs == cohort_merged IIDs — same
     filtered source). (--qc-fail also accepts a pre-assembled flat id file, unioned as
     'genotools_qc_fail'.)

GUARDRAIL: run by a human (reads sample-id rows). Writes files (which contain IIDs — the deliverable)
and prints ONLY aggregate summaries to stdout. The AI must not run this.

Outputs (into --out-dir):
  excludelist.txt        FID<TAB>IID                         (plink2 --remove format; the drops)
  exclude_reasons.tsv    FID  IID  reason  detail
  retained_manifest.csv  FID,IID,ancestry,call_rate,dup_cluster_id   (the seed step 6 extends with
                         source_callset + PC1..PCk)
"""
from pathlib import Path
from collections import defaultdict, Counter
import argparse
import os
import glob
import json
import re
import pandas as pd

# Paths follow config.sh, and resolve the same way it does: from this file's own location
# (scripts/ sits directly under the project root). Set WGS_ROOT to relocate, or pass explicit
# flags. No absolute path is hardcoded here — that is what let this file drift from config.sh.
WGS_ROOT = os.environ.get("WGS_ROOT") or str(Path(__file__).resolve().parent.parent)
DATA_ROOT = f"{WGS_ROOT}/data"
D_DEFAULT = f"{DATA_ROOT}/merged/relatedness"
GENOTOOLS_JSON_DEFAULT = [
    f"{DATA_ROOT}/amp-ad-genomics/WGS_Harmonization/genotools/FILTERED.wgs_harm.json",
    f"{DATA_ROOT}/amp-pd-genomics/WB-DWGS/genotools/FILTERED.wb_dwgs.json",
    f"{DATA_ROOT}/amp-ad-genomics/DivCo_HS/genotools/FILTERED.divco_hs.json",
    f"{DATA_ROOT}/amp-pd-genomics/BR-DSNWGS/genotools/FILTERED.br_dsnwgs.json",
]
QC_STEPS = ("callrate", "sex", "het")   # 'related' intentionally excluded — we own relatedness


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rel-dir", default=f"{D_DEFAULT}/by_ancestry",
                    help="dir containing relatedness_<ANC>_pairs.related")
    ap.add_argument("--smiss", default=f"{D_DEFAULT}/cohort_common.smiss")
    ap.add_argument("--labels", default=f"{D_DEFAULT}/common_ancestry_labels.txt",
                    help="FID<TAB>IID<TAB>ancestry for every common-set sample")
    ap.add_argument("--genotools-json", action="append", default=None,
                    help="per-callset genotools FILTERED.<dataset>.json (repeatable). "
                         "Defaults to the 3 callset JSONs. Use --no-genotools to skip.")
    ap.add_argument("--no-genotools", action="store_true",
                    help="skip genotools QC-fail assembly entirely")
    ap.add_argument("--qc-fail", default=None,
                    help="optional PRE-ASSEMBLED flat file of QC-fail ids (FID IID, or IID per line); "
                         "unioned in addition to --genotools-json")
    ap.add_argument("--dup-cutoff", type=float, default=0.354)
    ap.add_argument("--rel-cutoff", type=float, default=0.0884)
    ap.add_argument("--out-dir", default=D_DEFAULT)
    return ap.parse_args()


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def read_fileset_iids(stem):
    """IID set from a plink bed/pgen fileset given its STEM. genotools ran on --bfile, so prefer
    <stem>.fam (IID = whitespace column 2, no header); fall back to <stem>.psam (header-aware).
    Returns (iids, source_ext) or (None, None) if neither file exists."""
    fam = f"{stem}.fam"
    if Path(fam).exists():
        iids = set()
        with open(fam) as fh:
            for line in fh:
                toks = line.split()
                if len(toks) >= 2:
                    iids.add(toks[1])          # FID IID PAT MAT SEX PHENO
        return iids, "fam"
    psam = f"{stem}.psam"
    if Path(psam).exists():
        iids, idx = set(), 0
        with open(psam) as fh:
            for line in fh:
                s = line.rstrip("\n")
                if not s.strip():
                    continue
                toks = s.split()
                if s.startswith("#"):
                    hdr = [t.lstrip("#") for t in toks]
                    idx = hdr.index("IID") if "IID" in hdr else 0
                    continue
                if idx < len(toks):
                    iids.add(toks[idx])
        return iids, "psam"
    return None, None


def assemble_genotools_qc(json_paths, known_iids, iid_index):
    """Assemble callrate ∪ sex ∪ het fails from each callset's genotools `<ANC>_pass_fail` manifest.

    For each ancestry and each step in QC_STEPS, fails = (IIDs in `input` fileset) − (IIDs in
    `output` fileset). het is computed only when its status is True (skipped for small ancestries).
    'related' is never read. Returns (qc_hits, diag):
      qc_hits: {(FID,IID) -> step}  first step wins (callrate < sex < het) for a given sample.
      diag:    per-callset/-ancestry/-step counts + any file-read problems (aggregate; guardrail-safe).
    """
    qc_hits = {}
    raw_by_step = Counter()
    diag = []

    def flag(iid, step):
        for key in iid_index.get(iid, []):
            qc_hits.setdefault(key, step)

    for jpath in json_paths:
        if not jpath:
            continue
        cd = {"json": jpath, "exists": Path(jpath).exists(), "ancestries": {},
              "step_totals": Counter(), "problems": []}
        if not cd["exists"]:
            diag.append(cd)
            continue
        with open(jpath) as fh:
            data = json.load(fh)
        pf_keys = sorted(k for k in data if k.endswith("_pass_fail"))
        for key in pf_keys:
            anc = key[:-len("_pass_fail")]
            pf = data[key] or {}
            adiag = {}
            for step in QC_STEPS:
                cell = pf.get(step)
                if not isinstance(cell, dict):
                    adiag[step] = "no_cell"
                    continue
                if step == "het" and not cell.get("status"):
                    adiag[step] = "skipped"           # het not run for this (small) ancestry
                    continue
                in_iids, _ = read_fileset_iids(cell.get("input", ""))
                out_iids, _ = read_fileset_iids(cell.get("output", ""))
                if in_iids is None or out_iids is None:
                    adiag[step] = "MISSING_FILESET"
                    cd["problems"].append(f"{anc}:{step}: input/output fileset not found")
                    continue
                fails = in_iids - out_iids
                matched = {i for i in fails if i in known_iids}
                for i in matched:
                    flag(i, step)
                raw_by_step[step] += len(matched)
                cd["step_totals"][step] += len(matched)
                adiag[step] = {"fail": len(fails), "matched": len(matched),
                               "unmatched": len(fails) - len(matched)}
            cd["ancestries"][anc] = adiag
        diag.append(cd)
    return qc_hits, {"per_callset": diag, "raw_by_step": raw_by_step}


def main():
    a = parse_args()

    # ── call rate from .smiss (call_rate = 1 - F_MISS), keyed on (FID, IID) ──
    sm = pd.read_csv(a.smiss, sep=r"\s+")
    sm.columns = [c.lstrip("#") for c in sm.columns]
    sm["F_MISS"] = pd.to_numeric(sm["F_MISS"], errors="coerce").fillna(1.0)
    callrate = {(str(f), str(i)): 1.0 - m
                for f, i, m in zip(sm["FID"], sm["IID"], sm["F_MISS"])}

    def cr(k):
        return callrate.get(k, 0.0)   # unknown call rate -> treat as worst

    # ── ancestry labels (every common-set sample) ──
    lab = pd.read_csv(a.labels, sep=r"\s+", header=None,
                      names=["FID", "IID", "ancestry"], dtype=str)
    ancestry = {(f, i): anc for f, i, anc in zip(lab["FID"], lab["IID"], lab["ancestry"])}
    all_samples = list(ancestry.keys())

    # ── load all per-stratum KING pair tables ──
    frames = []
    for f in sorted(glob.glob(f"{a.rel_dir}/relatedness_*_pairs.related")):
        anc = Path(f).name[len("relatedness_"):-len("_pairs.related")]
        df = pd.read_csv(f)
        if df.empty:
            continue
        df.columns = [c.lstrip("#") for c in df.columns]
        df["stratum"] = anc
        frames.append(df)
    cols = ["FID1", "IID1", "FID2", "IID2", "KINSHIP", "stratum"]
    pairs = (pd.concat(frames, ignore_index=True) if frames
             else pd.DataFrame(columns=cols))
    pairs["KINSHIP"] = pd.to_numeric(pairs["KINSHIP"], errors="coerce")
    pairs["k1"] = list(zip(pairs["FID1"].astype(str), pairs["IID1"].astype(str)))
    pairs["k2"] = list(zip(pairs["FID2"].astype(str), pairs["IID2"].astype(str)))

    dropped = {}          # (FID,IID) -> (reason, detail)
    dup_cluster_id = {}   # (FID,IID) -> cluster id (for retained + dropped members)

    # ── 1. duplicates: cluster, keep highest call rate ──
    dup = pairs[pairs["KINSHIP"] >= a.dup_cutoff]
    uf = UnionFind()
    for k1, k2 in zip(dup["k1"], dup["k2"]):
        uf.union(k1, k2)
    clusters = defaultdict(list)
    for k in set(dup["k1"]).union(set(dup["k2"])):
        clusters[uf.find(k)].append(k)
    for idx, (_, members) in enumerate(sorted(clusters.items()), start=1):
        cid = f"dup{idx}"
        keep = max(members, key=lambda k: (cr(k), k))   # highest call rate wins
        for k in members:
            dup_cluster_id[k] = cid
            if k != keep:
                dropped[k] = ("duplicate",
                              f"{cid};n={len(members)};callrate={cr(k):.4f};kept_callrate={cr(keep):.4f}")

    # ── 2. relatives: greedy maximal-unrelated-set among still-retained ──
    rel = pairs[(pairs["KINSHIP"] > a.rel_cutoff) & (pairs["KINSHIP"] < a.dup_cutoff)]
    adj = defaultdict(set)
    for k1, k2 in zip(rel["k1"], rel["k2"]):
        if k1 in dropped or k2 in dropped:
            continue
        adj[k1].add(k2)
        adj[k2].add(k1)
    while True:
        active = [(k, len(v)) for k, v in adj.items() if v]
        if not active:
            break
        maxdeg = max(d for _, d in active)
        victim = min((k for k, d in active if d == maxdeg), key=lambda k: (cr(k), k))
        dropped[victim] = ("relative_2nd_deg", f"deg={len(adj[victim])};callrate={cr(victim):.4f}")
        for nb in list(adj[victim]):
            adj[nb].discard(victim)
        adj[victim].clear()

    # ── 3. genotools per-callset QC-fails (callrate ∪ sex ∪ het; related-prune IGNORED) ──
    iid_index = defaultdict(list)
    for (f, i) in all_samples:
        iid_index[i].append((f, i))
    known_iids = set(iid_index)

    gt_diag = None
    if not a.no_genotools:
        jpaths = a.genotools_json if a.genotools_json is not None else GENOTOOLS_JSON_DEFAULT
        jpaths = [g for g in jpaths if g]
        if jpaths:
            qc_hits, gt_diag = assemble_genotools_qc(jpaths, known_iids, iid_index)
            for k, step in qc_hits.items():
                dropped.setdefault(k, (f"genotools_qc_fail:{step}", ""))

    # optional pre-assembled flat QC-fail file (unioned in addition)
    if a.qc_fail:
        qc = pd.read_csv(a.qc_fail, sep=r"\s+", header=None, dtype=str)
        for _, row in qc.iterrows():
            vals = [v for v in row.tolist() if isinstance(v, str) and v.strip()]
            keys = []
            if len(vals) >= 2:
                keys = [(vals[0], vals[1])]
            elif len(vals) == 1:
                keys = iid_index.get(vals[0], [])
            for k in keys:
                dropped.setdefault(k, ("genotools_qc_fail", ""))

    # ── write outputs ──
    out = Path(a.out_dir)
    with open(out / "excludelist.txt", "w") as fx, open(out / "exclude_reasons.tsv", "w") as fr:
        fr.write("FID\tIID\treason\tdetail\n")
        for (fid, iid), (reason, detail) in sorted(dropped.items()):
            fx.write(f"{fid}\t{iid}\n")
            fr.write(f"{fid}\t{iid}\t{reason}\t{detail}\n")
    with open(out / "retained_manifest.csv", "w") as fm:
        fm.write("FID,IID,ancestry,call_rate,dup_cluster_id\n")
        for k in all_samples:
            if k in dropped:
                continue
            fid, iid = k
            fm.write(f"{fid},{iid},{ancestry.get(k,'')},{cr(k):.5f},{dup_cluster_id.get(k,'')}\n")

    # ── aggregate summary (guardrail-safe: counts only, no ids) ──
    by_reason = Counter(v[0] for v in dropped.values())
    n_total = len(all_samples)
    print("=== exclude-list summary (aggregate) ===")
    print(f"common-set samples:   {n_total}")
    print(f"total excluded:       {len(dropped)}")
    for reason, n in by_reason.most_common():
        print(f"    {reason:22s} {n}")
    print(f"retained:             {n_total - len(dropped)}")
    print("\n-- per stratum: pairs / dup-drops / relative-drops --")
    for anc in sorted(set(pairs["stratum"])) if len(pairs) else []:
        p = pairs[pairs["stratum"] == anc]
        ndup = sum(1 for k, v in dropped.items() if v[0] == "duplicate" and ancestry.get(k) == anc)
        nrel = sum(1 for k, v in dropped.items() if v[0] == "relative_2nd_deg" and ancestry.get(k) == anc)
        print(f"    {anc:4s} pairs={len(p):5d}  dup_drop={ndup:4d}  rel_drop={nrel:4d}")

    # ── genotools QC-fail diagnostics (aggregate; verify before applying) ──
    if gt_diag is None:
        print("\nNOTE: genotools QC assembly SKIPPED (--no-genotools). Excludelist is relatedness-only.")
    else:
        rbs = gt_diag["raw_by_step"]
        print("\n-- genotools QC-fails from <ANC>_pass_fail input-minus-output diffs "
              "(common-set samples flagged, pre-precedence) --")
        print(f"    callrate={rbs['callrate']:5d}  sex={rbs['sex']:5d}  het={rbs['het']:5d}"
              f"   (recorded reasons may differ: dup/relative take precedence)")
        for cd in gt_diag["per_callset"]:
            name = Path(cd["json"]).parent.parent.name
            if not cd["exists"]:
                print(f"    !! {name}: JSON NOT FOUND ({cd['json']}) — QC-fails MISSING for this callset")
                continue
            st = cd["step_totals"]
            print(f"    {name}: {len(cd['ancestries'])} ancestries  "
                  f"callrate={st['callrate']} sex={st['sex']} het={st['het']}")
            for anc, d in sorted(cd["ancestries"].items()):
                parts, warn = [], ""
                for step in QC_STEPS:
                    v = d.get(step)
                    if isinstance(v, dict):
                        parts.append(f"{step}={v['fail']}")
                        if v["unmatched"]:
                            warn += f"  [{step}_unmatched_to_common={v['unmatched']}]"
                    else:
                        parts.append(f"{step}:{v}")   # 'skipped' / 'no_cell' / 'MISSING_FILESET'
                        if v == "MISSING_FILESET":
                            warn += f"  <-- {step} fileset missing"
                print(f"        {anc:4s} {'  '.join(parts)}{warn}")
            for p in cd["problems"]:
                print(f"        !! {p}")
        print("    (het:skipped is expected for small ancestries; MISSING_FILESET is a real problem.)")
        print("    CROSS-CHECK: WGS_Harm sex total should ≈ pruned_samples sex=18 (probe finding).")
    print(f"\nwrote: {out}/excludelist.txt, exclude_reasons.tsv, retained_manifest.csv")


if __name__ == "__main__":
    main()
