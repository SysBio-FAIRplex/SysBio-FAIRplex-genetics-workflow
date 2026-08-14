#!/usr/bin/env python3
"""Mask cohort/batch-artifact variants in cross-cohort disease sumstats — POST-HOC (no plink re-run).

A control-vs-control GWAS (AMP-PD controls vs AMP-AD controls, same ancestry — the
`control@amppd:control@ampad` contrast in gwas_per_ancestry.sh) has NO true disease signal by
construction, so every variant it flags is a cohort/batch frequency artifact. This script reads that
scan's sumstats, takes the flagged variants as an artifact set, and applies it as a mask to the
EXISTING cross-cohort disease sumstats (e.g. PD-vs-AD) — annotating a `cc_artifact` column and writing
an artifact-removed copy, and reporting lambda_GC and genome-wide-hit counts before vs after.

Why post-hoc (not `plink2 --exclude` + re-run): the disease GWAS already ran; masking is a cheap
variant-level filter on its output, fully transparent and reversible.

CAVEATS (see docs/METHODS.md): the cc scan is power-limited by the smaller control arm (EUR AMP-AD ctrls ~406),
so it only catches LARGE cross-cohort frequency differences — a coarse filter, not a complete scrub.
It can also drop a real locus that happens to differ between the cohorts' controls (conservative).
Filtering is on p (power-dependent); interpret alongside the within-cohort scans, not as a lone fix.

GUARDRAIL: variant-level association statistics only (no subject rows) -> safe to run locally.

USAGE:
  python3 review/mask_cohort_artifacts.py \
      --cc      results/gwas/gwas_EUR_control_amppd_vs_control_ampad.pheno.glm.logistic.hybrid \
      --disease results/gwas/gwas_EUR_PD_vs_AD.*.hybrid \
                results/gwas/gwas_EUR_AD_vs_control.*.hybrid \
      --threshold 1e-4
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

CHI2_MED = stats.chi2.ppf(0.5, 1)                 # 0.4549
THRESH_GRID = [5e-8, 1e-6, 1e-5, 1e-4, 1e-3]      # sensitivity table


def find_col(hdr, cands):
    return next((c for c in cands if c in hdr), None)


def read_ids_p(path):
    """ADD-test rows of a plink2 sumstats as a df[ID, P] (numeric P, 0<P<=1)."""
    hdr = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    cid = find_col(hdr, ["ID", "SNP"])
    cp = find_col(hdr, ["P", "P_VALUE", "PVAL"])
    ct = find_col(hdr, ["TEST"])
    if not (cid and cp):
        sys.exit(f"ERROR: {path} has no ID/P columns ({hdr})")
    df = pd.read_csv(path, sep="\t", usecols=[c for c in (cid, cp, ct) if c])
    if ct and ct in df:
        df = df[df[ct] == "ADD"]
    df = df.rename(columns={cid: "ID", cp: "P"})
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df[df["P"].notna() & (df["P"] > 0) & (df["P"] <= 1)]
    return df[["ID", "P"]]


def lam_gc(p):
    p = np.asarray(p, dtype=float)
    p = p[(p > 0) & (p <= 1)]
    return float(stats.chi2.isf(np.median(p), 1) / CHI2_MED) if p.size else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cc", required=True, help="control-vs-control sumstats (the artifact source)")
    ap.add_argument("--disease", nargs="+", required=True, help="cross-cohort disease sumstats to mask")
    ap.add_argument("--threshold", type=float, default=1e-4, help="cc P below which a variant is an artifact (default 1e-4)")
    ap.add_argument("--sig", type=float, default=5e-8, help="genome-wide significance for the hit counts (default 5e-8)")
    ap.add_argument("--out-dir", default=None, help="default <dir-of-first-disease-file>/cc_masked")
    a = ap.parse_args()

    out = a.out_dir or os.path.join(os.path.dirname(a.disease[0]) or ".", "cc_masked")
    os.makedirs(out, exist_ok=True)

    # ── the control-vs-control artifact source ──
    cc = read_ids_p(a.cc)
    print(f"control-vs-control scan: {os.path.basename(a.cc)}")
    print(f"  variants tested: {len(cc):,}    lambda_GC: {lam_gc(cc['P'].to_numpy()):.3f}"
          f"   (its own inflation = how much cohort structure lives among controls)")
    print("  artifact variants by threshold:")
    for t in THRESH_GRID:
        n = int((cc["P"] < t).sum())
        mark = "  <- using" if t == a.threshold else ""
        print(f"    P < {t:>8.0e} : {n:>10,}{mark}")
    artifacts = set(cc.loc[cc["P"] < a.threshold, "ID"])
    print(f"  -> masking {len(artifacts):,} variants (cc P < {a.threshold:g})\n")

    # ── apply the mask to each disease scan ──
    hdr = f"{'disease scan':44}{'variants':>11}{'masked':>9}{'gwsig':>7}{'gwsig_masked':>13}{'lam_before':>11}{'lam_after':>10}"
    print(hdr); print("-" * len(hdr))
    for f in a.disease:
        cols = pd.read_csv(f, sep="\t", nrows=0).columns.tolist()
        cid = find_col(cols, ["ID", "SNP"])
        cp = find_col(cols, ["P", "P_VALUE", "PVAL"])
        df = pd.read_csv(f, sep="\t")
        pnum = pd.to_numeric(df[cp], errors="coerce")
        flag = df[cid].isin(artifacts)
        df["cc_artifact"] = flag.astype(int)

        lb, la = lam_gc(pnum.to_numpy()), lam_gc(pnum[~flag].to_numpy())
        gwsig = int((pnum < a.sig).sum())
        gwsig_masked = int(((pnum < a.sig) & flag).sum())

        base = os.path.basename(f)
        stem = base.replace(".hybrid", "").replace(".gz", "")
        df.to_csv(os.path.join(out, stem + ".annotated.tsv"), sep="\t", index=False)
        df[~flag].to_csv(os.path.join(out, stem + ".ccfiltered.tsv"), sep="\t", index=False)
        print(f"{base[:44]:44}{len(df):>11,}{int(flag.sum()):>9,}{gwsig:>7}{gwsig_masked:>13}{lb:>11.3f}{la:>10.3f}")

    print(f"\nwrote *.annotated.tsv (+cc_artifact col) and *.ccfiltered.tsv (artifacts removed) -> {out}")
    print("read: 'masked' = how many of the scan's variants are batch artifacts; 'gwsig_masked' = how")
    print("      many of its genome-wide hits were artifacts; lam_after < lam_before => the mask removed inflation.")


if __name__ == "__main__":
    main()
