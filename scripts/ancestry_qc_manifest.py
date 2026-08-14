#!/usr/bin/env python3
"""STEP 6 (boundary) — assemble retained_samples_manifest.csv.

Joins the step-5 retained manifest + the step-6 per-ancestry PCs + per-callset source into the
single artifact the phenotype side consumes (clinical_core.ipynb §12 reads it to build the grain):
    IID, source_callset, ancestry, call_rate, dup_cluster_id, PC1..PCk

Sources (all cluster-side, so the manifest is emitted entirely on the genetics track):
  --manifest   data/merged/relatedness/retained_manifest.csv   (FID,IID,ancestry,call_rate,dup_cluster_id)
  --pca-dir    data/merged/by_ancestry_qc/                      (cohort_<ANC>_pca.eigenvec; #FID IID PC1..)
  --label      the 3 per-callset genotools label files -> source_callset by IID membership
               (FID<TAB>IID<TAB>ancestry_label; an IID appears in exactly its own callset's file)

PCs are blank for samples in strata that had no PCA (the <50-sample strata plink2 refused to LD-prune —
EAS/MDE/FIN/CAS/SAS); those are below the >=100-cases GWAS-viability cut anyway.

GUARDRAIL: run by a human (reads IID rows). Writes the manifest (IIDs = the deliverable) and prints
ONLY aggregate counts. The AI must not run it.
"""
from pathlib import Path
from collections import Counter, defaultdict
import argparse
import csv
import glob
import os

# Paths follow config.sh, and resolve the same way it does: from this file's own location
# (scripts/ sits directly under the project root). Set WGS_ROOT to relocate, or pass explicit
# flags. No absolute path is hardcoded here — that is what let this file drift from config.sh.
WGS_ROOT = os.environ.get("WGS_ROOT") or str(Path(__file__).resolve().parent.parent)
GT = f"{WGS_ROOT}/data"
LABELS_DEFAULT = [
    f"wgs_harm={GT}/amp-ad-genomics/WGS_Harmonization/genotools/FILTERED.wgs_harm_ancestry_umap_linearsvc_predicted_labels.txt",
    f"wb_dwgs={GT}/amp-pd-genomics/WB-DWGS/genotools/FILTERED.wb_dwgs_ancestry_umap_linearsvc_predicted_labels.txt",
    f"divco_hs={GT}/amp-ad-genomics/DivCo_HS/genotools/FILTERED.divco_hs_ancestry_umap_linearsvc_predicted_labels.txt",
    f"br_dsnwgs={GT}/amp-pd-genomics/BR-DSNWGS/genotools/FILTERED.br_dsnwgs_ancestry_umap_linearsvc_predicted_labels.txt",
]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=f"{WGS_ROOT}/data/merged/relatedness/retained_manifest.csv")
    ap.add_argument("--pca-dir", default=f"{WGS_ROOT}/data/merged/by_ancestry_qc")
    ap.add_argument("--label", action="append", default=None,
                    help="name=path (repeatable); defaults to the 3 callset label files")
    ap.add_argument("--out", default=f"{WGS_ROOT}/data/merged/by_ancestry_qc/retained_samples_manifest.csv")
    ap.add_argument("--npc", type=int, default=10)
    return ap.parse_args()


def read_eigenvec(path):
    """{IID -> [PC1..PCn as strings]} from a plink2 .eigenvec (#FID IID PC1..)."""
    out = {}
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
        iid_i = header.index("IID")
        pc_i = [k for k, c in enumerate(header) if c.upper().startswith("PC")]
        for line in fh:
            t = line.split()
            if len(t) > iid_i:
                out[t[iid_i]] = [t[k] for k in pc_i]
    return out


def read_labels_source(spec):
    """spec 'name=path' -> {IID -> name} for every IID in that callset's label file (IID = col 2)."""
    name, path = spec.split("=", 1)
    ids = set()
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if len(t) >= 2 and not line.startswith("#"):
                ids.add(t[1])
    return name, ids


def main():
    a = parse_args()

    # ── source_callset by IID membership across the 3 per-callset label files ──
    label_specs = a.label if a.label is not None else LABELS_DEFAULT
    source = defaultdict(list)          # IID -> [callset names it appears in]
    for spec in label_specs:
        name, ids = read_labels_source(spec)
        for i in ids:
            source[i].append(name)

    def source_callset(iid):
        s = source.get(iid, [])
        if not s:
            return "unknown"
        return s[0] if len(s) == 1 else "|".join(sorted(set(s)))

    # ── per-ancestry PCs ──
    pcs = {}                            # IID -> [PC strings]
    npc_seen = 0
    for f in sorted(glob.glob(f"{a.pca_dir}/cohort_*_pca.eigenvec")):
        ev = read_eigenvec(f)
        for iid, v in ev.items():
            pcs[iid] = v
            npc_seen = max(npc_seen, len(v))
    npc = min(a.npc, npc_seen) if npc_seen else a.npc

    # ── join over the retained manifest ──
    pc_cols = [f"PC{i}" for i in range(1, npc + 1)]
    n = 0
    n_with_pc = 0
    src_counts = Counter()
    unknown_src = 0
    multi_src = 0
    pc_by_anc = defaultdict(lambda: [0, 0])   # ancestry -> [with_pc, total]

    with open(a.manifest) as fin, open(a.out, "w", newline="") as fout:
        r = csv.DictReader(fin)
        w = csv.writer(fout)
        w.writerow(["IID", "source_callset", "ancestry", "call_rate", "dup_cluster_id"] + pc_cols)
        for row in r:
            iid = row["IID"]
            anc = row.get("ancestry", "")
            sc = source_callset(iid)
            src_counts[sc] += 1
            if sc == "unknown":
                unknown_src += 1
            elif "|" in sc:
                multi_src += 1
            pv = pcs.get(iid)
            has_pc = pv is not None
            pc_out = (pv + [""] * npc)[:npc] if has_pc else [""] * npc
            w.writerow([iid, sc, anc, row.get("call_rate", ""), row.get("dup_cluster_id", "")] + pc_out)
            n += 1
            n_with_pc += has_pc
            pc_by_anc[anc][0] += has_pc
            pc_by_anc[anc][1] += 1

    # ── aggregate summary (guardrail-safe) ──
    print("=== retained_samples_manifest summary ===")
    print(f"rows written:        {n}   -> {a.out}")
    print(f"PC columns:          PC1..PC{npc}")
    print(f"with PCs:            {n_with_pc}    without PCs: {n - n_with_pc} (sub-50-sample strata)")
    print(f"source_callset:      {dict(src_counts)}")
    if unknown_src:
        print(f"  !! {unknown_src} rows with UNKNOWN source (IID not in any label file) — investigate")
    if multi_src:
        print(f"  !! {multi_src} rows with MULTI source (IID in >1 label file) — investigate")
    print("PC coverage per ancestry (with_pc / total):")
    for anc in sorted(pc_by_anc):
        wpc, tot = pc_by_anc[anc]
        print(f"    {anc:4s} {wpc:6d} / {tot:6d}")


if __name__ == "__main__":
    main()
