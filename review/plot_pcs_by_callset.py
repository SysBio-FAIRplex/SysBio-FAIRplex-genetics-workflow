#!/usr/bin/env python3
"""Do the callsets separate in PC space? Scatter + a quantitative answer.

Reads retained_samples_manifest.csv (IID, source_callset, ancestry, PC1..PC10) and,
per ancestry, plots PC1 vs PC2 coloured by source callset.

The plot answers "is it visible"; the eta-squared table answers "how much". eta^2 is the
share of a PC's variance explained by callset membership — 0 means the PC is blind to
which cohort a sample came from, 1 means the PC *is* the cohort label. That number is the
point: cohort structure in the PCs means the GWAS covariates are partly absorbing batch,
which is protective for confounding but costs power on any cross-cohort contrast.

Runs LOCALLY on a downloaded manifest — no cluster, no genotypes.

Usage:
    python3 plot_pcs_by_callset.py --manifest path/to/retained_samples_manifest.csv
    python3 plot_pcs_by_callset.py --manifest OLD.csv --label pre_pcfix   # compare runs
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Categorical slots 1-3 of the validated palette; all-pairs CVD dE 9.2, normal-vision 16.6.
# A scatter is an ALL-PAIRS form, and only the first three slots clear the all-pairs floors —
# slot 4 puts yellow beside orange and fails (normal-vision 13.7). So values past the third take
# the neutral plus a distinct marker: secondary encoding, not a fourth hue.
#
# TWO neutrals do not fit either — MEASURED 2026-08-18 against the L 0.43-0.77 band on this
# surface: #52514e scores normal-vision dE 14.0 against #7a7973 and #5f5e5a scores 9.4 (both
# below the 15 floor); #adaca6 clears it at 15.6 but drops CVD to 6.1 protan and contrast to
# 2.22; #b8b7b1 falls outside the band at 0.779. So br_dsnwgs and the fused rows SHARE the
# neutral and are separated by marker shape alone.
#
# That is the better treatment here regardless of palette: the marker branch below also raises
# size and alpha, which is what makes a 95-sample callset visible against 10,780 EUR samples.
# A hue would not have.
COLORS = {
    "wb_dwgs":            "#2a78d6",   # AMP-PD blood
    "wgs_harm":           "#eb6834",   # AMP-AD trio
    "divco_hs":           "#1baf7a",   # Diverse Cohorts
    "divco_hs|wgs_harm":  "#7a7973",   # fused at merge (n~87)
    "br_dsnwgs":          "#7a7973",   # AMP-PD postmortem (n=95) — neutral + '^'
}
MARKERS = {"divco_hs|wgs_harm": "x", "br_dsnwgs": "^"}
LABELS = {
    "wb_dwgs": "AMP-PD blood (wb_dwgs)",
    "wgs_harm": "AMP-AD trio (wgs_harm)",
    "divco_hs": "Diverse Cohorts (divco_hs)",
    "divco_hs|wgs_harm": "fused (both)",
    "br_dsnwgs": "AMP-PD postmortem (br_dsnwgs)",
}

INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#e3e2dd"
SURFACE = "#fcfcfb"


def eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    """Share of variance in `values` explained by `groups` (one-way ANOVA effect size)."""
    grand = values.mean()
    ss_total = ((values - grand) ** 2).sum()
    if ss_total == 0:
        return 0.0
    ss_between = sum(
        (groups == g).sum() * (values[groups == g].mean() - grand) ** 2
        for g in np.unique(groups)
    )
    return float(ss_between / ss_total)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    # Resolved from this file, not the CWD — downloaded outputs live under <bundle>/results/
    # so the bundle stays self-contained.
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent.parent / "results" / "pca"))
    ap.add_argument("--label", default="")
    ap.add_argument("--npc", type=int, default=10)
    args = ap.parse_args()

    df = pd.read_csv(args.manifest, dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    pcs = [f"PC{i}" for i in range(1, args.npc + 1) if f"PC{i}" in df.columns]
    for c in pcs:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    have = df.dropna(subset=["PC1", "PC2"])
    print(f"{len(df):,} rows, {len(have):,} with PCs "
          f"({len(df) - len(have):,} without — sub-50 strata)\n")

    ancs = [a for a in have.ancestry.unique() if (have.ancestry == a).sum() >= 2]
    ancs = sorted(ancs, key=lambda a: -(have.ancestry == a).sum())

    print("samples per callset per ancestry:")
    print(pd.crosstab(have.ancestry, have.source_callset, margins=True).to_string())

    # ---- quantitative: how much of each PC is cohort? ----
    rows = []
    for anc in ancs:
        g = have[have.ancestry == anc]
        if g.source_callset.nunique() < 2:
            continue
        rows.append({"ancestry": anc, "n": len(g),
                     "callsets": g.source_callset.nunique(),
                     **{pc: round(eta_squared(g[pc].to_numpy(), g.source_callset.to_numpy()), 3)
                        for pc in pcs}})
    eta = pd.DataFrame(rows).set_index("ancestry")
    print("\neta^2 — share of each PC's variance explained by source callset")
    print("(0 = PC is blind to cohort; 1 = PC is the cohort label)")
    print(eta.to_string())
    if len(eta):
        worst = eta[pcs].max(axis=1)
        print("\nstrongest cohort signal in any PC, per ancestry:")
        for a, v in worst.sort_values(ascending=False).items():
            pc = eta.loc[a, pcs].idxmax()
            flag = "  <-- PC is largely cohort" if v >= 0.25 else ""
            print(f"  {a:5} {pc:>4}  eta^2={v:.3f}{flag}")

    # ---- small multiples ----
    ncol = 3
    nrow = int(np.ceil(len(ancs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 4.0 * nrow),
                             facecolor=SURFACE, squeeze=False)

    for ax, anc in zip(axes.ravel(), ancs):
        g = have[have.ancestry == anc]
        ax.set_facecolor(SURFACE)
        # largest group first so it sits underneath and never hides the small ones
        for src in g.source_callset.value_counts().index:
            s = g[g.source_callset == src]
            ax.scatter(s.PC1, s.PC2,
                       s=26 if src in MARKERS else max(4, 2200 / len(g)),
                       c=COLORS.get(src, "#7a7973"),
                       marker=MARKERS.get(src, "o"),
                       alpha=0.85 if src in MARKERS else min(0.8, 220 / len(g) + 0.18),
                       linewidths=0.9 if src in MARKERS else 0,
                       label=LABELS.get(src, src), rasterized=True)
        e1 = eta.loc[anc, "PC1"] if anc in eta.index else float("nan")
        e2 = eta.loc[anc, "PC2"] if anc in eta.index else float("nan")
        ax.set_title(f"{anc}   n={len(g):,}", color=INK, fontsize=11, loc="left", pad=8)
        ax.set_xlabel(f"PC1   $\\eta^2$={e1:.2f}", color=INK_MUTED, fontsize=9)
        ax.set_ylabel(f"PC2   $\\eta^2$={e2:.2f}", color=INK_MUTED, fontsize=9)
        ax.tick_params(colors=INK_MUTED, labelsize=8)
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)

    for ax in axes.ravel()[len(ancs):]:
        ax.set_visible(False)

    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False,
                     fontsize=9, labelcolor=INK_MUTED, bbox_to_anchor=(0.5, -0.01),
                     markerscale=2.2)
    # Swatches must be solid — the plot alpha is for overplotting, not for identity.
    for lh in leg.legend_handles:
        lh.set_alpha(1.0)
    suffix = f" — {args.label}" if args.label else ""
    fig.suptitle(f"Population structure by source callset{suffix}",
                 color=INK, fontsize=13, x=0.01, ha="left", y=1.0)
    fig.tight_layout(rect=[0, 0.04, 1, 0.97])

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    stem = f"pcs_by_callset{'_' + args.label if args.label else ''}"
    fig.savefig(out / f"{stem}.png", dpi=150, facecolor=SURFACE, bbox_inches="tight")
    eta.to_csv(out / f"{stem}_eta2.csv")
    print(f"\nwrote {out / (stem + '.png')}")
    print(f"      {out / (stem + '_eta2.csv')}")


if __name__ == "__main__":
    main()
