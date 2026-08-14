#!/usr/bin/env python3
"""Did the PCs actually change between two retained_samples_manifest.csv files?

Eyeballing one IID is unreliable: an eigenvector's sign is arbitrary, so a PC can be
numerically identical in structure yet print as its own negative. Nearby eigenvalues can
also swap PC order between runs. This merges on IID and reports, per PC:

  exact   — how many samples have a byte-identical value (the "nothing changed" signal)
  r       — Pearson correlation, same PC index, old vs new (sign-sensitive)
  |r|     — magnitude, so a pure sign flip reads as 1.00
  best    — the old PC each new PC most resembles, catching a reordering

Usage:
    python3 review/compare_pcs.py OLD.csv NEW.csv
"""
import sys

import numpy as np
import pandas as pd


def load(path, npc=10):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    d.columns = [c.strip() for c in d.columns]
    pcs = [f"PC{i}" for i in range(1, npc + 1) if f"PC{i}" in d.columns]
    for c in pcs:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    keep = ["IID"] + (["ancestry"] if "ancestry" in d.columns else []) + pcs
    return d[keep].dropna(subset=pcs[:1]), pcs


def compare(m, pcs, header):
    """One PC-by-PC table for an already-subset frame."""
    print(header)
    print(f"{'PC':>5} {'exact':>8} {'r':>8} {'|r|':>7}   best match in old")
    identical = True
    changed = 0
    for pc in pcs:
        a = m[f"{pc}_old"].to_numpy()
        b = m[f"{pc}_new"].to_numpy()
        exact = int((a == b).sum())
        identical &= exact == len(m)
        r = np.corrcoef(a, b)[0, 1]
        if abs(r) < 0.99:
            changed += 1
        best_pc, best_r = max(
            ((p, np.corrcoef(m[f"{p}_old"].to_numpy(), b)[0, 1]) for p in pcs),
            key=lambda t: abs(t[1]),
        )
        note = "" if best_pc == pc else f"  <-- reordered (was {best_pc})"
        print(f"{pc:>5} {exact:>8,} {r:>8.3f} {abs(r):>7.3f}   {best_pc} |r|={abs(best_r):.3f}{note}")

    # Subspace overlap: all ten PCs enter --glm together, so what matters for the model
    # is the SPAN of PC1..PC10, not which axis is which. A rotation between adjacent PCs
    # (equal-ish eigenvalues) leaves the span untouched and adjusts for exactly the same
    # structure. Mean squared canonical correlation: 1.0 = same covariate space, 0 = disjoint.
    A = np.linalg.qr(m[[f"{p}_old" for p in pcs]].to_numpy() -
                     m[[f"{p}_old" for p in pcs]].to_numpy().mean(0))[0]
    B = np.linalg.qr(m[[f"{p}_new" for p in pcs]].to_numpy() -
                     m[[f"{p}_new" for p in pcs]].to_numpy().mean(0))[0]
    sv = np.linalg.svd(A.T @ B, compute_uv=False)
    overlap = float((sv ** 2).mean())
    print(f"      subspace overlap (span of PC1..PC{len(pcs)}): {overlap:.3f}"
          f"   [1.0 = identical covariate space]")
    return identical, changed, overlap


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    old, pcs = load(sys.argv[1])
    new, _ = load(sys.argv[2])
    m = old.merge(new, on="IID", suffixes=("_old", "_new"))
    print(f"old {len(old):,} rows | new {len(new):,} rows | matched on IID {len(m):,}\n")
    if not len(m):
        sys.exit("no shared IIDs — wrong files?")

    anc_col = "ancestry_old" if "ancestry_old" in m.columns else (
        "ancestry" if "ancestry" in m.columns else None)

    if anc_col is None:
        identical, changed, ov = compare(m, pcs, "pooled across all samples")
        groups = [("all", identical, changed, len(m), ov)]
    else:
        # PCs are computed PER ANCESTRY, so a pooled correlation mixes six separate
        # PCAs on unrelated scales and is dominated by the largest stratum. Compare
        # within stratum — that is the only comparison the numbers actually support.
        groups = []
        for anc, g in sorted(m.groupby(anc_col), key=lambda t: -len(t[1])):
            if len(g) < 3:
                continue
            print()
            identical, changed, ov = compare(g, pcs, f"=== {anc}  (n={len(g):,}) ===")
            groups.append((anc, identical, changed, len(g), ov))

    print("\n" + "=" * 60)
    if all(g[1] for g in groups):
        print("VERDICT: every PC is byte-identical — the new manifest was NOT pulled,")
        print("         or step 6 wrote the same PCs. Check before launching step 7.")
    else:
        print("VERDICT: PCs changed. Per stratum — axes differing (|r|<0.99) and span overlap:")
        print(f"  {'anc':6} {'axes':>7}  {'span':>6}   n")
        for anc, _, changed, n, ov in groups:
            print(f"  {anc:6} {changed:2}/{len(pcs):<4} {ov:6.3f}   {n:,}")
        print("\n  |r| near 1.00 with r near -1.00 is only a sign flip, not a real change.")
        print("  Whether the NEW PCs are better is a separate question — compare eta^2")
        print("  from plot_pcs_by_callset.py on both manifests.")


if __name__ == "__main__":
    main()
