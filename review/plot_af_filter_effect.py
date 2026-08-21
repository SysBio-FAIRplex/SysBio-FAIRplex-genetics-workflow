#!/usr/bin/env python3
"""Does the AF-concordance filter remove callset structure from the PCs? Before vs after.

This is the proof figure for the filter, and it is built to be shown to someone who was not
in the room. It answers three questions in the order a sceptic asks them:

  1. HOW BIG WAS THE PROBLEM?   eta^2 of source callset on each PC, before filtering.
  2. DID THE FILTER FIX IT?     the same number after, on the same samples and the same
                                pruning/PCA settings — the ONLY difference between the two
                                generations is the exclusion list.
  3. WHERE DID IT NOT?          the strata that barely move are drawn at the same scale as
                                the ones that do, because a filter that works everywhere is a
                                claim this cohort does not support.

eta^2 is the share of a PC's variance explained by callset membership: 0 means the PC is blind
to which cohort a sample came from, 1 means the PC *is* the cohort label. Cohort structure in
the PCs means the GWAS covariates are partly absorbing batch — protective against confounding,
but it costs power on exactly the cross-cohort contrast this study is built around, and it
cannot protect the association test itself, where no PC adjustment reaches.

INPUT — the two manifests step 6 writes in one job:
    by_ancestry_qc/unfiltered/retained_samples_manifest.csv   (before)
    by_ancestry_qc/retained_samples_manifest.csv              (after)

Both come out of the SAME step-6 submission, which matters: an earlier version of this
comparison was assembled from two separate runs, and one of them had silently picked up a stale
exclusion list, so the "unfiltered" baseline was not unfiltered. Pairing them by construction is
what makes the delta trustworthy.

Runs LOCALLY on downloaded manifests — no cluster, no genotypes.

Usage:
    python3 plot_af_filter_effect.py --before UNFILTERED.csv --after FILTERED.csv
    python3 plot_af_filter_effect.py --before B.csv --after A.csv --focus EUR AJ
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_pcs_by_callset import COLORS, LABELS, MARKERS, eta_squared  # noqa: E402

# Text and surface tokens, shared with plot_pcs_by_callset.py.
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#e3e2dd"
SURFACE = "#fcfcfb"

# The effect panel is deliberately ACHROMATIC. Hue in this figure means one thing only —
# which callset a sample came from — and spending two more hues on before/after would put four
# colours on screen carrying two different kinds of meaning. Before/after is encoded by fill
# (hollow -> solid) plus position, with both ends directly labelled, so identity never rests on
# colour at all. The validated categorical slots stay reserved for the scatter panels.
BEFORE_EDGE = "#8f8e88"
AFTER_FILL = INK


def eta_table(manifest: Path, npc: int) -> pd.DataFrame:
    """-> DataFrame indexed by ancestry: n, callsets, PC1..PCk eta^2."""
    df = pd.read_csv(manifest, dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    pcs = [f"PC{i}" for i in range(1, npc + 1) if f"PC{i}" in df.columns]
    for c in pcs:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    have = df.dropna(subset=["PC1", "PC2"])

    rows = []
    for anc in have.ancestry.unique():
        g = have[have.ancestry == anc]
        if len(g) < 2 or g.source_callset.nunique() < 2:
            continue
        rows.append({"ancestry": anc, "n": len(g), "callsets": g.source_callset.nunique(),
                     **{pc: eta_squared(g[pc].to_numpy(), g.source_callset.to_numpy())
                        for pc in pcs}})
    if not rows:
        sys.exit(f"{manifest}: no stratum has PCs and >=2 callsets — nothing to compare")
    return pd.DataFrame(rows).set_index("ancestry"), have, pcs


def scatter(ax, g, anc, e1, e2, title):
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
    ax.set_title(title, color=INK, fontsize=10.5, loc="left", pad=6)
    ax.set_xlabel(f"PC1   $\\eta^2$={e1:.3f}", color=INK_MUTED, fontsize=9)
    ax.set_ylabel(f"PC2   $\\eta^2$={e2:.3f}", color=INK_MUTED, fontsize=9)
    ax.tick_params(colors=INK_MUTED, labelsize=8)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", required=True, help="unfiltered retained_samples_manifest.csv")
    ap.add_argument("--after", required=True, help="filtered retained_samples_manifest.csv")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent.parent / "results" / "pca"))
    ap.add_argument("--focus", nargs="+", default=None,
                    help="strata to draw PC scatters for (default: the 2 largest)")
    ap.add_argument("--npc", type=int, default=10)
    ap.add_argument("--label", default="")
    a = ap.parse_args()

    eb, gb, pcs_b = eta_table(Path(a.before), a.npc)
    ea, ga, pcs_a = eta_table(Path(a.after), a.npc)
    pcs = [p for p in pcs_b if p in pcs_a]

    strata = [s for s in eb.index if s in ea.index]
    if not strata:
        sys.exit("no stratum has PCs in BOTH generations — nothing to compare")
    dropped = sorted(set(eb.index) ^ set(ea.index))
    if dropped:
        print(f"note: {', '.join(dropped)} has PCs in only one generation — omitted\n")

    # ── the table. This is the deliverable; the figure renders it. ──
    summary = pd.DataFrame({
        "n": eb.loc[strata, "n"],
        "callsets": eb.loc[strata, "callsets"],
        "max_eta2_before": eb.loc[strata, pcs].max(axis=1),
        "worst_pc_before": eb.loc[strata, pcs].idxmax(axis=1),
        "max_eta2_after": ea.loc[strata, pcs].max(axis=1),
        "worst_pc_after": ea.loc[strata, pcs].idxmax(axis=1),
    })
    summary["reduction"] = 1 - summary.max_eta2_after / summary.max_eta2_before.replace(0, np.nan)
    summary = summary.sort_values("max_eta2_before", ascending=False)

    print("max eta^2 in ANY PC, per stratum — before vs after the AF-concordance filter")
    print("(0 = PCs blind to cohort; 1 = a PC IS the cohort label)\n")
    print(summary.to_string(float_format=lambda v: f"{v:.3f}"))

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    stem = f"af_filter_effect{'_' + a.label if a.label else ''}"
    per_pc = pd.concat({"before": eb.loc[strata, pcs], "after": ea.loc[strata, pcs]},
                       names=["generation"])
    summary.to_csv(out / f"{stem}.csv")
    per_pc.to_csv(out / f"{stem}_per_pc.csv")

    # Default focus is the two LARGEST strata, not the two worst. Those are the ones carrying
    # the GWAS-viable contrasts, so they are what a reader needs to see — and if the largest
    # stratum is also the one the filter fails on, the figure should say so rather than quietly
    # showing two successes.
    focus = a.focus or list(summary.sort_values("n", ascending=False).index[:2])
    focus = [f for f in focus if f in strata]

    # ── figure: headline dumbbell over per-stratum PC scatters ──
    nrow = 1 + len(focus)
    fig = plt.figure(figsize=(9.6, 3.4 + 4.2 * len(focus)), facecolor=SURFACE)
    gs = fig.add_gridspec(nrow, 2, height_ratios=[1.25] + [1.6] * len(focus),
                          hspace=0.42, wspace=0.26,
                          left=0.135, right=0.975, top=0.90, bottom=0.10)

    # Panel A — dumbbell. One row per stratum, before -> after, sorted by the problem's size.
    axd = fig.add_subplot(gs[0, :])
    axd.set_facecolor(SURFACE)
    ys = np.arange(len(summary))[::-1]
    for y, (anc, r) in zip(ys, summary.iterrows()):
        axd.plot([r.max_eta2_before, r.max_eta2_after], [y, y],
                 color=GRID, linewidth=2.4, solid_capstyle="round", zorder=1)
        axd.scatter(r.max_eta2_before, y, s=64, facecolors="none",
                    edgecolors=BEFORE_EDGE, linewidths=1.8, zorder=3)
        axd.scatter(r.max_eta2_after, y, s=52, color=AFTER_FILL, zorder=4)
        # Direct labels at both ends, each placed on the OUTSIDE of the pair. Stacking them
        # above/below instead collides whenever the two values are close — which is exactly
        # the case the figure most needs to render honestly, since a stratum the filter barely
        # moves is the finding, not a rendering edge case.
        b, af = r.max_eta2_before, r.max_eta2_after
        (b_dx, b_ha), (a_dx, a_ha) = (((8, "left"), (-8, "right")) if af <= b
                                      else ((-8, "right"), (8, "left")))
        axd.annotate(f"{b:.3f}", (b, y), xytext=(b_dx, 0), textcoords="offset points",
                     ha=b_ha, va="center", color=INK_MUTED, fontsize=8.5)
        axd.annotate(f"{af:.3f}", (af, y), xytext=(a_dx, 0), textcoords="offset points",
                     ha=a_ha, va="center", color=INK, fontsize=8.5, fontweight="bold")
        if pd.notna(r.reduction):
            axd.annotate(f"−{100 * r.reduction:.0f}%", (1.20, y), ha="right", va="center",
                         color=INK_MUTED, fontsize=8.5)

    axd.set_yticks(ys)
    axd.set_yticklabels([f"{anc}  n={int(r.n):,}" for anc, r in summary.iterrows()],
                        color=INK, fontsize=9.5)
    axd.set_xlim(-0.16, 1.22)
    axd.set_xticks(np.arange(0, 1.01, 0.2))
    axd.set_ylim(-0.7, len(summary) - 0.3)
    axd.set_xlabel("max $\\eta^2$ of source callset across PC1–PC%d"
                   "          (right column: reduction)" % len(pcs),
                   color=INK_MUTED, fontsize=9)
    axd.tick_params(axis="x", colors=INK_MUTED, labelsize=8)
    axd.tick_params(axis="y", length=0)
    axd.grid(True, axis="x", color=GRID, linewidth=0.6)
    axd.set_axisbelow(True)
    for side in ("top", "right", "left"):
        axd.spines[side].set_visible(False)
    axd.spines["bottom"].set_color(GRID)

    handles = [plt.Line2D([], [], marker="o", linestyle="none", markersize=8,
                          markerfacecolor="none", markeredgecolor=BEFORE_EDGE,
                          markeredgewidth=1.8, label="before — unfiltered"),
               plt.Line2D([], [], marker="o", linestyle="none", markersize=7.5,
                          color=AFTER_FILL, label="after — AF-concordance filter applied")]
    axd.legend(handles=handles, loc="lower left", frameon=False, fontsize=8.5,
               labelcolor=INK_MUTED, ncol=2, bbox_to_anchor=(-0.10, 1.02),
               handletextpad=0.5, columnspacing=1.8)

    # Panels B.. — the same strata in PC space, so the number has a picture behind it.
    # The two panels of a row SHARE limits. Letting each autoscale would rescale the after
    # panel to fill its box and hide the very thing being shown: the clusters collapsing into
    # one cloud. (PC sign and rotation are arbitrary between runs, so read the spread, not the
    # orientation — the eta^2 on each axis is the quantity that is actually comparable.)
    for i, anc in enumerate(focus):
        gens = (("before — unfiltered", gb, eb), ("after — filtered", ga, ea))
        subs = [g[g.ancestry == anc] for _, g, _ in gens]
        xs = pd.concat([s.PC1 for s in subs])
        ys_ = pd.concat([s.PC2 for s in subs])
        pad_x, pad_y = 0.05 * (xs.max() - xs.min()), 0.05 * (ys_.max() - ys_.min())
        xlim = (xs.min() - pad_x, xs.max() + pad_x)
        ylim = (ys_.min() - pad_y, ys_.max() + pad_y)

        for j, ((gen, _, et), g) in enumerate(zip(gens, subs)):
            ax = fig.add_subplot(gs[1 + i, j])
            scatter(ax, g, anc, et.loc[anc, "PC1"], et.loc[anc, "PC2"], f"{anc}   {gen}")
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)

    # One legend for callset identity, gathered from the first scatter and shown once — the
    # four scatter panels share an encoding, so repeating it per panel would be noise.
    ax0 = fig.axes[1]
    h, l = ax0.get_legend_handles_labels()
    leg = fig.legend(h, l, loc="upper center", ncol=min(len(l), 3), frameon=False,
                     fontsize=9, labelcolor=INK_MUTED, bbox_to_anchor=(0.5, 0.058),
                     markerscale=2.2)
    for lh in leg.legend_handles:
        lh.set_alpha(1.0)

    suffix = f" — {a.label}" if a.label else ""
    fig.suptitle(f"Callset structure in the ancestry PCs, before and after the "
                 f"AF-concordance filter{suffix}",
                 color=INK, fontsize=13, x=0.012, ha="left", y=0.985)
    # No tight_layout: the gridspec already carries explicit margins, and tight_layout cannot
    # reconcile them with the two figure-level legends (it warns and then guesses).
    fig.savefig(out / f"{stem}.png", dpi=150, facecolor=SURFACE)

    print(f"\nwrote {out / (stem + '.png')}")
    print(f"      {out / (stem + '.csv')}")
    print(f"      {out / (stem + '_per_pc.csv')}")


if __name__ == "__main__":
    main()
