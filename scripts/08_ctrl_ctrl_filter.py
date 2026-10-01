#!/usr/bin/env python3
"""STEP 8 — annotate every contrast with the control-vs-control result, and emit a filtered copy.

`control@amppd` vs `control@ampad` has no true disease signal by construction, so a genome-wide
hit there is a candidate cohort artifact that survived step 6a.

ANNOTATES, NEVER OVERWRITES. The control arms are screened for different diseases (AMP-PD controls
for PD only), so AD risk alleles are genuinely depleted in the AMP-AD arm and APOE is EXPECTED to
reach significance here; subtracting would delete the study's strongest true locus (METHODS.md §8).
Writes `.ccannot.tsv` (every row + CTRL_P — judge hits there, per variant) and `.ccfilt.tsv`
(a sensitivity copy) beside the untouched primary.

No gene-list report: an uncited locus list makes a reassuring NEGATIVE printable. Name the gene a
hit sits in with  python3 scripts/gene_annot.py --at chr19:44908684

Reads summary statistics only — no genotypes, no sample IDs. Prints aggregate counts.
"""
from pathlib import Path
import argparse
import sys

CTRL_TAG = "control_amppd_vs_control_ampad"


def read_p(path, pmax):
    """-> {id: p} for rows with P < pmax. Streams; the flagged set is small."""
    flagged = {}
    n = 0
    with open(path) as fh:
        hdr = fh.readline().lstrip("#").rstrip("\n").split("\t")
        try:
            i_id, i_p = hdr.index("ID"), hdr.index("P")
        except ValueError:
            sys.exit(f"{path}: no ID/P column in header: {hdr}")
        for line in fh:
            t = line.rstrip("\n").split("\t")
            n += 1
            v = t[i_p]
            if v in ("NA", ""):
                continue
            try:
                pv = float(v)
            except ValueError:
                continue
            if pv < pmax:
                flagged[t[i_id]] = pv
    return flagged, n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gwas-dir", required=True, help="step 7's output dir")
    ap.add_argument("--ctrl-p", type=float, default=1e-5,
                    help="control-vs-control significance for flagging (default 1e-5)")
    ap.add_argument("--gwsig", type=float, default=5e-8)
    a = ap.parse_args()

    gd = Path(a.gwas_dir)
    ctrl_files = sorted(gd.glob(f"gwas_*_{CTRL_TAG}.filtered.tsv"))
    if not ctrl_files:
        sys.exit(f"no gwas_*_{CTRL_TAG}.filtered.tsv under {gd} — did step 7 run that contrast?")

    for cf in ctrl_files:
        anc = cf.name[len("gwas_"):-(len(CTRL_TAG) + len(".filtered.tsv") + 1)]
        flagged, n_tested = read_p(cf, a.ctrl_p)
        print("=" * 78)
        print(f"{anc}: control@amppd vs control@ampad")
        print(f"  variants tested                       : {n_tested:,}")
        print(f"  flagged at P < {a.ctrl_p:g}                 : {len(flagged):,} "
              f"(expected by chance: {n_tested * a.ctrl_p:,.0f})")
        gw = {v: p for v, p in flagged.items() if p < a.gwsig}
        print(f"  of those, genome-wide significant     : {len(gw):,}")
        print("  There is no true signal in this contrast, so an excess over chance is residual")
        print("  cohort artifact that survived step 6a. Near chance = 6a did its job.")

        # Unconditional: this holds for EVERY flagged variant.
        print("  AMP-AD controls are screened cognitively normal, AMP-PD controls are screened for")
        print("  PD only — so AD risk alleles are genuinely depleted in the AMP-AD arm and a real")
        print("  AD locus is EXPECTED here. Judge flagged variants individually in .ccannot.tsv;")
        print("  `gene_annot.py --at <chr:pos>` names the gene a given hit sits in.")

        # ── annotate + filter every other contrast in this ancestry ──
        others = [p for p in sorted(gd.glob(f"gwas_{anc}_*.filtered.tsv"))
                  if CTRL_TAG not in p.name]
        if not others:
            print("  (no other contrasts for this ancestry)")
            continue
        print(f"\n  {'contrast':34} {'hits':>7} {'flagged':>8} {'kept':>7}")
        for op in others:
            tag = op.name[len(f"gwas_{anc}_"):-len(".filtered.tsv")]
            ann = op.with_suffix("").with_suffix(".ccannot.tsv")
            fil = op.with_suffix("").with_suffix(".ccfilt.tsv")
            n_hit = n_drop = n_keep = 0
            with open(op) as fh, open(ann, "w") as fa, open(fil, "w") as ff:
                hdr = fh.readline().rstrip("\n")
                cols = hdr.lstrip("#").split("\t")
                i_id, i_p = cols.index("ID"), cols.index("P")
                fa.write(hdr + "\tCTRL_P\tCTRL_FLAG\n")
                ff.write(hdr + "\n")
                for line in fh:
                    line = line.rstrip("\n")
                    t = line.split("\t")
                    vid = t[i_id]
                    cp = flagged.get(vid)
                    fa.write(f"{line}\t{'NA' if cp is None else f'{cp:.3g}'}\t"
                             f"{'1' if cp is not None else '0'}\n")
                    if cp is None:
                        ff.write(line + "\n")
                    pv = t[i_p]
                    if pv not in ("NA", ""):
                        try:
                            if float(pv) < a.gwsig:
                                n_hit += 1
                                if cp is not None:
                                    n_drop += 1
                                else:
                                    n_keep += 1
                        except ValueError:
                            pass
            print(f"  {tag:34} {n_hit:7,} {n_drop:8,} {n_keep:7,}")
        print()

    print("=" * 78)
    print("Wrote, per contrast, alongside the untouched primary sumstats:")
    print("  .ccannot.tsv  every row + CTRL_P + CTRL_FLAG   <- use this to judge hits individually")
    print("  .ccfilt.tsv   rows NOT flagged in controls     <- the sensitivity-analysis copy")
    print("\nReport the primary result as primary. The filtered copy is a sensitivity analysis,")
    print("not a replacement — the control arms differ by disease screening, not only by batch.")


if __name__ == "__main__":
    main()
