#!/usr/bin/env python3
"""STEP 8 — annotate every contrast with the control-vs-control result, and emit a filtered copy.

WHY. The primary contrast is confounded with callset: AMP-AD supplies the AD cases, AMP-PD the PD
cases. Step 6a removes the variants where that confound is visible at the genotype level, but no
filter is complete. `control@amppd` vs `control@ampad` has no true disease signal by construction —
both arms are non-cases — so anything genome-wide significant there is cohort artifact that survived
6a. This is the demo's step 13 (`AD_versus_PD_GWAS.sh`), which subtracts those variants from the
case-vs-case results.

WHY IT ANNOTATES RATHER THAN OVERWRITES. The two control groups are differentially screened for the
disease we are testing. AMP-AD controls are assessed as cognitively normal; AMP-PD controls are
screened as "No PD Nor Other Neurological Disorder" — i.e. screened for PD, not for AD. So AD risk
alleles are genuinely depleted in the AMP-AD control arm, and APOE is expected to reach significance
in the control-vs-control scan for an entirely real reason. Subtracting it destructively would
delete the strongest true locus in the study, which is very likely why the demo reports APOE from
its unfiltered sumstats and keeps the filtered set as a separate file.

So this writes `.ccfilt.tsv` alongside the primary and never modifies it, adds a CTRL_P column to a
`.ccannot.tsv` so any hit can be judged individually, and reports explicitly whether known AD/PD
loci are among the flagged — a hit there is the screening asymmetry showing up, not an artifact.

GUARDRAIL: reads summary statistics only — no genotypes, no sample IDs. Prints aggregate counts.
"""
from pathlib import Path
import argparse
import sys

SENTINEL_LOCI = [
    ("APOE/rs429358", "19", 44_908_684, 50_000),
    ("APOE/rs7412",   "19", 44_908_822, 50_000),
    ("TREM2",         "6",  41_160_000, 50_000),
    ("BIN1",          "2", 127_100_000, 50_000),
    ("CR1",           "1", 207_500_000, 50_000),
    ("MAPT/17q21.31", "17", 45_900_000, 200_000),
    ("SNCA",          "4",  89_700_000, 100_000),
    ("GBA1",          "1", 155_230_000, 50_000),
    ("LRRK2",         "12", 40_200_000, 100_000),
]

CTRL_TAG = "control_amppd_vs_control_ampad"


def parse_id(vid):
    p = vid.split(":")
    if len(p) < 2:
        return None
    try:
        return (p[0][3:] if p[0].startswith("chr") else p[0]), int(p[1])
    except ValueError:
        return None


def sentinel_hits(ids):
    pos = [(v, parse_id(v)) for v in ids]
    out = []
    for name, c, p, win in SENTINEL_LOCI:
        near = [v for v, cp in pos if cp and cp[0] == c and abs(cp[1] - p) <= win]
        if near:
            out.append((name, near))
    return out


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
                    help="control-vs-control significance for flagging (demo uses 1e-5)")
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

        sh = sentinel_hits(flagged)
        if sh:
            print("\n  KNOWN LOCI among the flagged — read before filtering on them:")
            for name, near in sh:
                print(f"    {name:16} {len(near)}: {', '.join(sorted(near)[:4])}")
            print("    AMP-AD controls are screened cognitively normal, AMP-PD controls are")
            print("    screened for PD only — so AD risk alleles are genuinely depleted in the")
            print("    AMP-AD arm. A hit here is that asymmetry, NOT an artifact. Filtering it")
            print("    out would delete real signal. This is why .ccfilt.tsv is a separate file.")
        else:
            print("\n  no known AD/PD locus among the flagged")

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
    print("  .ccfilt.tsv   rows NOT flagged in controls     <- the demo's step-13 equivalent")
    print("\nReport the primary result as primary. The filtered copy is a sensitivity analysis,")
    print("not a replacement — the control arms differ by disease screening, not only by batch.")


if __name__ == "__main__":
    main()
