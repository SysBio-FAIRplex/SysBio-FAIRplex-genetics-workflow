#!/usr/bin/env python3
"""QQ + Manhattan plots for the per-ancestry per-contrast GWAS (pipeline step 7, review side).

Runs LOCALLY on the sumstats rsync'd back from biowulf. For each contrast that ran, draws a stacked
QQ (top) + Manhattan (bottom) figure, annotated with the arm counts, cohort-confound tag, and genomic
inflation. Genomic inflation is recomputed here from the full P distribution and cross-checked against
the lambda_gc the batch job already wrote to gwas_summary.csv (they should match to ~2 decimals).

GUARDRAIL: sumstats are VARIANT-LEVEL association statistics (no subject rows), so this is safe to run
locally / by the agent. It reads only #CHROM/POS/P (+A1_FREQ) and the summary CSV — never a genotype.

INPUT  step 7's output dir: data/merged/by_ancestry_qc/gwas/, which is where this looks by default
— no copying, because code and data share a root. Running on a laptop instead, put a copy of that
directory at results/gwas/ and it is found there as a fallback.
    (needs gwas_summary.csv + gwas_<ANC>_<CONTRAST>.filtered.tsv, or the raw
     gwas_<ANC>_<CONTRAST>.pheno.glm.logistic.hybrid — either naming is accepted)

USAGE (from anywhere — paths resolve off this file, not the CWD):
    python3 review/plot_gwas.py                          # everything in <bundle>/results/gwas
    python3 review/plot_gwas.py --viable-only            # only >=100/arm contrasts
    python3 review/plot_gwas.py --gwas-dir results/gwas_v2 --out-dir results/figures
    python3 review/plot_gwas.py --sig 5e-8 --min-maf 0.01

Figures default to <gwas-dir>/plots/qqman_<ANC>_<CONTRAST>.png.
"""
import argparse
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

CHI2_MED = stats.chi2.ppf(0.5, 1)          # 0.4549 — median of the null chi-square(1)
# Resolved from this file, not the CWD, so the script works from anywhere. Prefer step 7's real
# output dir; fall back to results/gwas for the case where the dir was copied to a laptop. The
# old default was results/gwas unconditionally, which never resolved — results/ is regenerable
# and is not carried in the repo.
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_CLUSTER_DIR = os.path.join(PROJECT_ROOT, "data", "merged", "by_ancestry_qc", "gwas")
_LOCAL_DIR = os.path.join(PROJECT_ROOT, "results", "gwas")
DEFAULT_DIR = _CLUSTER_DIR if os.path.isdir(_CLUSTER_DIR) else _LOCAL_DIR
# alternating chromosome colors (2-tone, colorblind-safe)
CHR_COLORS = ("#3b6ea5", "#9ec1e3")
SIG = 5e-8
SUGG = 1e-5


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gwas-dir", default=DEFAULT_DIR, help=f"dir with sumstats + gwas_summary.csv (default {DEFAULT_DIR})")
    ap.add_argument("--out-dir", default=None, help="output dir for PNGs (default <gwas-dir>/plots)")
    ap.add_argument("--viable-only", action="store_true", help="only plot contrasts flagged viable_ge100=yes")
    ap.add_argument("--sig", type=float, default=SIG, help=f"genome-wide significance line (default {SIG})")
    ap.add_argument("--min-maf", type=float, default=0.0, help="drop variants below this MAF (uses A1_FREQ; default 0)")
    ap.add_argument("--max-null-points", type=int, default=200_000,
                    help="Manhattan: max points kept among P>0.01 (all P<=0.01 always kept); default 200k")
    return ap.parse_args()


def find_col(hdr, cands):
    return next((c for c in cands if c in hdr), None)


def load_sumstats(path, min_maf):
    """Return a tidy df with CHR(int), POS(int), P(float) [MAF-filtered], or None if unreadable."""
    hdr = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    c_chr = find_col(hdr, ["#CHROM", "CHROM", "#CHR", "CHR"])
    c_pos = find_col(hdr, ["POS", "BP"])
    c_p = find_col(hdr, ["P", "P_VALUE", "PVAL"])
    c_test = find_col(hdr, ["TEST"])
    c_frq = find_col(hdr, ["A1_FREQ", "A1FREQ", "FREQ"])
    if not (c_chr and c_pos and c_p):
        print(f"  !! {os.path.basename(path)}: missing CHROM/POS/P columns {hdr}", file=sys.stderr)
        return None
    use = [c for c in (c_chr, c_pos, c_p, c_test, c_frq) if c]
    df = pd.read_csv(path, sep="\t", usecols=use)
    if c_test and c_test in df:
        df = df[df[c_test] == "ADD"]
    df = df.rename(columns={c_chr: "CHR_RAW", c_pos: "POS", c_p: "P"})
    # P -> numeric, keep 0<P<=1
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df[df["P"].notna() & (df["P"] > 0) & (df["P"] <= 1)]
    if min_maf > 0 and c_frq and c_frq in df:
        f = pd.to_numeric(df[c_frq], errors="coerce")
        maf = np.minimum(f, 1 - f)
        df = df[maf >= min_maf]
    # chrom -> int (autosomes; strip 'chr', map X/Y just in case)
    chrmap = {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}
    cr = df["CHR_RAW"].astype(str).str.replace("chr", "", case=False, regex=False)
    df["CHR"] = cr.replace(chrmap)
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df = df[df["CHR"].notna()]
    df["CHR"] = df["CHR"].astype(int)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df = df[df["POS"].notna()]
    df["POS"] = df["POS"].astype(np.int64)
    return df[["CHR", "POS", "P"]].reset_index(drop=True)


def lambda_gc(p):
    """Genomic inflation factor from a P vector (median chi2 / 0.4549); monotonic => use median P."""
    med = np.median(p)
    return float(stats.chi2.isf(med, 1) / CHI2_MED)


def lambda_1000(lam, n_case, n_ctrl):
    """Rescale lambda_GC to a 1000/1000 study. None when the smaller arm < 200 (rescale unstable) or
    the inputs are degenerate — matches the guard in gwas_per_ancestry.sh."""
    if lam is None or lam <= 0 or n_case <= 0 or n_ctrl <= 0 or min(n_case, n_ctrl) < 200:
        return None
    l = 1 + (lam - 1) * ((1 / n_case) + (1 / n_ctrl)) / (2 / 1000)
    return l if l > 0 else None


def qq_panel(ax, p, lam_data, lam_sum, lam1k_sum):
    n = len(p)
    obs = -np.log10(np.sort(p))
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    # thin the dense null end for a light figure; keep the informative tail (obs>2) fully
    keep = np.ones(n, dtype=bool)
    dense = exp < 2.0
    if dense.sum() > 50_000:
        idx = np.where(dense)[0]
        drop = np.random.default_rng(0).choice(idx, size=idx.size - 50_000, replace=False)
        keep[drop] = False
    ax.scatter(exp[keep], obs[keep], s=6, c="#333333", edgecolors="none", rasterized=True)
    lim = max(exp.max(), obs.max()) * 1.05
    ax.plot([0, lim], [0, lim], color="#c0392b", lw=1)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$")
    ax.set_ylabel(r"Observed $-\log_{10}(p)$")
    txt = fr"$\lambda_{{GC}}$ = {lam_data:.3f} (data)"
    if lam_sum is not None:
        txt += f"\nsummary $\\lambda_{{GC}}$ = {lam_sum:.3f}"
    if lam1k_sum is not None:
        txt += f"\n$\\lambda_{{1000}}$ = {lam1k_sum:.3f}"
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", ha="left",
            fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="#999999", alpha=0.9))
    ax.set_title("QQ", fontsize=10, loc="left")


def manhattan_panel(ax, df, sig, max_null):
    df = df.sort_values(["CHR", "POS"])
    # downsample the null cloud (P>0.01), keep every P<=0.01
    hi = df[df["P"] <= 0.01]
    lo = df[df["P"] > 0.01]
    if len(lo) > max_null:
        lo = lo.sample(n=max_null, random_state=0)
    d = pd.concat([hi, lo]).sort_values(["CHR", "POS"])
    d["mlog"] = -np.log10(d["P"])
    # cumulative x offset per chromosome
    x = np.zeros(len(d), dtype=np.float64)
    ticks, labels = [], []
    off = 0.0
    order = sorted(d["CHR"].unique())
    pos_all = d["POS"].to_numpy()
    chr_all = d["CHR"].to_numpy()
    fill = np.zeros(len(d))
    for i, c in enumerate(order):
        m = chr_all == c
        cpos = pos_all[m]
        span = cpos.max() - cpos.min() if cpos.size else 0
        fill[m] = off + (cpos - (cpos.min() if cpos.size else 0))
        ticks.append(off + span / 2.0)
        labels.append(str(c))
        off += span + span * 0.02 + 1
    d = d.assign(x=fill)
    for i, c in enumerate(order):
        m = d["CHR"] == c
        ax.scatter(d.loc[m, "x"], d.loc[m, "mlog"], s=4,
                   c=CHR_COLORS[i % 2], edgecolors="none", rasterized=True)
    ax.axhline(-np.log10(sig), color="#c0392b", lw=0.9, ls="--")
    ax.axhline(-np.log10(SUGG), color="#e08e0b", lw=0.8, ls=":")
    ax.set_xticks(ticks); ax.set_xticklabels(labels, fontsize=7)
    ax.set_xlim(0, off)
    ymax = max(d["mlog"].max(), -np.log10(sig) + 1)
    ax.set_ylim(0, ymax * 1.05)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(p)$")
    nhit = int((df["P"] <= sig).sum())
    ax.set_title(f"Manhattan   (genome-wide hits P<={sig:g}: {nhit})", fontsize=10, loc="left")


def main():
    a = parse_args()
    gdir = a.gwas_dir
    out_dir = a.out_dir or os.path.join(gdir, "plots")
    os.makedirs(out_dir, exist_ok=True)

    # summary (for counts/confound/lambda annotations); optional but expected
    summary = {}
    spath = os.path.join(gdir, "gwas_summary.csv")
    if os.path.exists(spath):
        sdf = pd.read_csv(spath)
        for _, r in sdf.iterrows():
            summary[(r["ancestry"], r["contrast"])] = r
    else:
        print(f"note: no gwas_summary.csv in {gdir} — plotting without annotations", file=sys.stderr)

    # Two naming forms, both from step 7: the raw plink2 output and the ADD-only |BETA|-filtered
    # copy. The .filtered.tsv is what gets rsync'd back, because the hybrids are several times
    # larger. Both carry the same ADD rows, so lambda from either matches gwas_summary.csv's.
    files = sorted(glob.glob(os.path.join(gdir, "gwas_*.pheno.glm.logistic*"))
                   + glob.glob(os.path.join(gdir, "gwas_*.filtered.tsv")))
    if not files:
        print(f"No sumstats found in {gdir} — expected gwas_<ANC>_<CONTRAST>.filtered.tsv or "
              f"gwas_<ANC>_<CONTRAST>.pheno.glm.logistic*", file=sys.stderr)
        sys.exit(1)

    made = 0
    print(f"{'ANC':4} {'contrast':16} {'nSNP':>10} {'lam_data':>9} {'lam_summary':>12}  status")
    for f in files:
        base = os.path.basename(f)
        m = re.match(r"gwas_([A-Za-z]+)_(.+?)\.(?:pheno\.glm\.logistic|filtered\.tsv)", base)
        if not m:
            continue
        anc, contrast = m.group(1), m.group(2)
        row = summary.get((anc, contrast))
        if a.viable_only and row is not None and str(row.get("viable_ge100")) != "yes":
            continue

        df = load_sumstats(f, a.min_maf)
        if df is None or df.empty:
            print(f"{anc:4} {contrast:16} {'-':>10} {'-':>9} {'-':>12}  EMPTY/unreadable")
            continue
        lam = lambda_gc(df["P"].to_numpy())
        lam_sum = float(row["lambda_gc"]) if row is not None and pd.notna(row.get("lambda_gc")) else None
        nc = int(float(row["n_case"])) if row is not None and pd.notna(row.get("n_case")) else 0
        nk = int(float(row["n_ctrl"])) if row is not None and pd.notna(row.get("n_ctrl")) else 0
        # compute lambda_1000 from the local (scipy) lambda + arm counts, guarded — do NOT trust the
        # summary's raw value (unguarded runs wrote nonsense like -5.09 / 0.26 for small arms).
        lam1k = lambda_1000(lam, nc, nk)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 9),
                                       gridspec_kw={"height_ratios": [1, 1.1]})
        qq_panel(ax1, df["P"].to_numpy(), lam, lam_sum, lam1k)
        manhattan_panel(ax2, df, a.sig, a.max_null_points)

        # figure title with the review-relevant metadata
        if row is not None:
            # callset_one_sided / n_diffmiss_excluded are read with .get so an older
            # gwas_summary.csv (pre-2026-08-21, no callset-skew columns) still plots.
            one_sided = str(row.get("callset_one_sided", "") or "")
            skew = ""
            if one_sided not in ("", "none", "NA", "nan"):
                skew = f"   one-sided callset={one_sided} (diffmiss removed {row.get('n_diffmiss_excluded', '?')})"
            sub = (f"case={row['case_arm']} (n={row['n_case']}, {row['case_pct_amppd']}% AMP-PD)   "
                   f"ctrl={row['ctrl_arm']} (n={row['n_ctrl']}, {row['ctrl_pct_amppd']}% AMP-PD)   "
                   f"confound={row['confound_tag']} (Δ={row['delta_amppd']})   "
                   f"viable≥100={row['viable_ge100']}{skew}")
        else:
            sub = f"{len(df):,} variants"
        fig.suptitle(f"{anc} — {contrast.replace('_', ' ')}", fontsize=13, y=0.995)
        fig.text(0.5, 0.955, sub, ha="center", fontsize=8.5, color="#444444")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        outp = os.path.join(out_dir, f"qqman_{anc}_{contrast}.png")
        fig.savefig(outp, dpi=150)
        plt.close(fig)
        made += 1

        flag = ""
        if lam_sum is not None and abs(lam - lam_sum) > 0.02:
            flag = "  !! lambda mismatch vs summary"
        print(f"{anc:4} {contrast:16} {len(df):>10,} {lam:>9.3f} "
              f"{(lam_sum if lam_sum is not None else float('nan')):>12.3f}  -> {os.path.basename(outp)}{flag}")

    print(f"\n{made} figure(s) -> {out_dir}")


if __name__ == "__main__":
    main()
