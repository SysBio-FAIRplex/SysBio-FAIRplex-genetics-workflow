#!/usr/bin/env python3
"""Does the AF-concordance filter remove callset structure from the PCs? Before vs after.

The proof figure for the filter (METHODS.md §6.4): eta^2 of source callset on each PC, before
and after, on the same samples and PCA settings. Every stratum is drawn at the same scale, so
the ones the filter barely moves stay visible.

INPUT — the two manifests ONE step-6 submission writes:
    by_ancestry_qc/unfiltered/retained_samples_manifest.csv   (before)
    by_ancestry_qc/retained_samples_manifest.csv              (after)
They must come from the same job: pairing different runs can pick up a different exclusion
list, and then "before" is not unfiltered.

Runs on biowulf against those manifests in place — no genotypes. The PNG's scatter panels plot one
point per participant, so it stays on the cluster; only the two eta^2 CSVs leave it.

Usage:
    python3 review/plot_af_filter_effect.py                        # step 6's pair, by default
    python3 review/plot_af_filter_effect.py --focus EUR AJ
    python3 review/plot_af_filter_effect.py --before B.csv --after A.csv   # any other pair
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

# The effect panel is achromatic: hue means callset only, so before/after is encoded by fill
# (hollow -> solid) and position, with both ends labelled.
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
    # Defaults are the pair ONE step-6 submission writes (see INPUT above), read in place.
    root = Path(__file__).resolve().parent.parent
    qc = root / "data" / "merged" / "by_ancestry_qc"
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", default=str(qc / "unfiltered" / "retained_samples_manifest.csv"),
                    help="unfiltered retained_samples_manifest.csv")
    ap.add_argument("--after", default=str(qc / "retained_samples_manifest.csv"),
                    help="filtered retained_samples_manifest.csv")
    ap.add_argument("--out", default=str(root / "results" / "pca"))
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

    # Default focus: the two LARGEST strata (they carry the viable contrasts), not the two worst,
    # so a failure in the largest stratum is shown rather than hidden.
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
        # Label each end on the outside of the pair, so close values (a stratum the filter barely
        # moves) never collide.
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

    # Panels B..: the same strata in PC space. A row's two panels share limits so the collapse
    # is visible; PC sign and rotation are arbitrary between runs, so compare spread, not
    # orientation.
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
    # No tight_layout: it cannot reconcile the gridspec's margins with the figure-level legends.
    fig.savefig(out / f"{stem}.png", dpi=150, facecolor=SURFACE)

    print(f"\nwrote {out / (stem + '.png')}")
    print(f"      {out / (stem + '.csv')}")
    print(f"      {out / (stem + '_per_pc.csv')}")


if __name__ == "__main__":
    main()
