#!/usr/bin/env python3
"""Build the per-callset AF-concordance exclusion list for step 6.

WHY, measured on this cohort (2026-08-19, four-callset merge). Applying the list this script builds
takes EUR's worst callset eta^2 from 0.757 to 0.036 -- a 95.2% reduction from excluding 0.058% of
EUR's variants. So the callset separation in PC space IS variant-intrinsic and concentrated, and the
flagged variants are ~7x enriched for duplicate-pair genotype discordance, where the true genotype
is identical by construction. Technical breakage, not population structure and not disease.

SCOPE OF THAT RESULT -- IT DOES NOT GENERALISE ACROSS STRATA. The same list moves AJ only 0.984 ->
0.962 (2.2%). AJ is never evaluated here: its largest same-dx callset pair is AJ/control at
wgs_harm=44, under MIN_CELL, so the list carries no AJ-derived flags and EUR-derived flags do not
transfer. AJ also shows no HWE excess (0.24x its chance expectation). Do not claim this filter
addresses AJ; the leading hypothesis there is real sub-continental structure. See PROJECT_LOG.md.

An earlier version of this docstring cited `diag_af_crossstratum` for 147 EUR-flagged variants
taking AJ's eta^2 from 0.748 to 0.055. That script is in no commit, and the claim is contradicted
by the measurement above: 93.6% of the July list survives in the current one, so those 147 are
almost certainly a subset of today's 4,415, and removing them plus 4,268 others moves AJ by 0.022.
The 0.748 was EUR's baseline mislabelled as AJ's. Recorded so it is not resurrected.

DISEASE HELD CONSTANT WITHIN EVERY COMPARISON, AND WHY THAT IS THE WHOLE DESIGN. wgs_harm/divco_hs
supply the AD cases and wb_dwgs the PD cases, so the primary contrast is confounded with callset and
a raw between-callset frequency difference has two possible causes: technical breakage, or a real
disease effect. The diagnostic could ignore this because it ran on the LD-pruned PCA input, which
happens to exclude APOE. A GENOME-WIDE filter has no such protection — rs429358 differs sharply
between AD-source and PD-source samples for an entirely real reason, and an all-samples filter would
delete the strongest true locus in the study.

So every comparison here is made WITHIN a (stratum x dx) cell. Inside a cell, disease is constant,
so a true disease effect cannot produce a between-callset difference and only technical breakage
can. Two cells carry most of the weight:

  controls  wb_dwgs vs wgs_harm vs divco_hs   the wgs_harm<->wb_dwgs axis that drives PC1
  AD cases  wgs_harm vs divco_hs              lifted vs NATIVE, disease matched by construction

The flags are unioned across cells. Note this is not selection on the outcome: cells are internally
disease-constant, and the list is applied to a case-vs-case contrast.

NOISE CALIBRATION. A fixed |dAF| threshold is sample-size dependent — a 120-sample cell throws far
more noise flags than a 3,000-sample one, which would over-filter small cells AND make the pairwise
rates incomparable for the liftover question below. A variant is flagged only if it clears BOTH
|dAF| > --thresh AND z > --zmin, where z uses each arm's actual per-variant OBS_CT (so differential
missingness is handled too).

SCOPE. Applied to the ASSOCIATION set, not just the PCA input. A variant mismapped badly enough to
bend PC1 produces a spurious association in the test itself, where no PC adjustment reaches it.

WHAT THE PAIRWISE ROLLUP DECIDES. divco_hs vs wb_dwgs is BOTH-NATIVE. If it diverges as badly as the
wgs_harm pairs, the mechanism is calling/mapping, no liftover change could fix it, and step 0 stays
closed. If only the wgs_harm pairs diverge, the breakage is liftover-specific.

SECOND STAGE — PER-CALLSET HWE (the GenoTools `--all_variant` equivalent). The frequency test above
is effect-based: it sees that two callsets disagree, without knowing why. HWE is mechanism-based. A
mismapped variant pools reads from two near-identical genomic locations, which inflates
heterozygosity and breaks Hardy-Weinberg *in the callset carrying the error*. Step 6 computes HWE
POOLED across callsets at 1e-6, so a deviation confined to wgs_harm (~1,540 of ~9,800 EUR samples)
is averaged against a clean majority before the test sees it. Here it runs within each callset, at
GenoTools' 1e-4.

CONTROLS ONLY, for the same reason as above and the reason GenoTools sets filter_controls=True:
case ascertainment genuinely distorts HWE at a real disease locus, so running this over an
AD-case-heavy callset would flag APOE for an entirely real reason. Restricting to controls removes
that channel. keep-fewhet is retained — mismapping produces het EXCESS, which is still caught, while
het deficiency (which pooling and substructure create artificially) is spared.

The two stages are unioned into one exclusion list, and the overlap is reported: HWE independently
re-finding the frequency-flagged variants is a mechanism-based confirmation of an effect-based test.

GUARDRAIL: human-run. Reads the id-bearing grain and genotypes via plink2, writes the exclusion
list, prints ONLY aggregate counts.
"""
from collections import defaultdict
from pathlib import Path
import argparse
import csv
import subprocess
import sys

# Known loci we must NOT silently delete. GRCh38. Tripwire only — nothing is auto-whitelisted.
# A hit means the cell design did not fully separate technical from real, or that locus is
# genuinely broken in one callset. Either way: stop and look.
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


def run(cmd):
    """plink2 with output captured. stderr=STDOUT because capture_output alone hides plink's
    errors entirely (documented gotcha) — a silent failure here would look like a real result."""
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        sys.exit(f"FAILED: {' '.join(cmd)}\n{r.stdout[-4000:]}")
    return r.stdout


def read_grain(path):
    """-> {iid: (callset, ancestry, dx)}. Columns are positional, per step 7's contract:
    1 IID, 3 source_callset, 4 ancestry, 7 dx_detailed. The PCs (13-22) are deliberately NOT read —
    they go stale the moment step 6 re-runs, and nothing here depends on them."""
    out = {}
    with open(path, newline="") as fh:
        rdr = csv.reader(fh)
        next(rdr, None)
        for row in rdr:
            if len(row) >= 7 and row[0].strip():
                out[row[0].strip()] = (row[2].strip(), row[3].strip(), row[6].strip())
    return out


def read_afreq(path):
    """-> {id: (alt_freq, obs_ct)}. OBS_CT is the ALLELE count, so it already absorbs missingness."""
    out = {}
    with open(path) as fh:
        hdr = fh.readline().lstrip("#").rstrip("\n").split("\t")
        i_id, i_f, i_n = hdr.index("ID"), hdr.index("ALT_FREQS"), hdr.index("OBS_CT")
        for line in fh:
            t = line.rstrip("\n").split("\t")
            try:
                out[t[i_id]] = (float(t[i_f]), int(t[i_n]))
            except ValueError:
                pass          # plink writes NA where an arm has no called genotypes
    return out


def parse_id(vid):
    """IDs embed alleles: chr19:44908684:T:C -> ('19', 44908684). None if unparseable."""
    p = vid.split(":")
    if len(p) < 2:
        return None
    try:
        return (p[0][3:] if p[0].startswith("chr") else p[0]), int(p[1])
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc-dir", required=True, help="step 6's by_ancestry_qc")
    ap.add_argument("--grain", required=True)
    ap.add_argument("--out", required=True, help="exclusion list to write")
    ap.add_argument("--work", required=True)
    ap.add_argument("--thresh", type=float, default=0.05,
                    help="effect-size floor. 0.05 is where the diagnostic's arms flattened: "
                         "0.10->eta2 .055, 0.05->.025, 0.02->.021 for 4,684 more variants")
    ap.add_argument("--zmin", type=float, default=5.0,
                    help="significance floor, so the flag rate does not track cell size")
    ap.add_argument("--min-cell", type=int, default=100,
                    help="samples per callset per cell needed to compare at all")
    ap.add_argument("--control-dx", default="control",
                    help="dx_detailed value marking a control. Stage 2 only — the AF stage needs "
                         "no designated control label, since every cell is disease-constant")
    ap.add_argument("--hwe", type=float, default=1e-4,
                    help="per-callset HWE threshold, controls only. 1e-4 is GenoTools' "
                         "--all_variant default. Set 0 to disable the stage entirely")
    ap.add_argument("--min-hwe-controls", type=int, default=50,
                    help="controls per callset per stratum needed to test HWE there")
    ap.add_argument("--hwe-both-tails", action="store_true",
                    help="drop keep-fewhet, i.e. also exclude het-DEFICIENT variants as plink1.9 "
                         "does. Off by default: mismapping causes het excess, while deficiency is "
                         "what substructure produces artificially")
    ap.add_argument("--mishap", type=float, default=0.0,
                    help="per-callset haplotype-missingness threshold (GenoTools' `haplotype` step, "
                         "default 1e-4 there). 0 = OFF, which is the default here because plink1.9 "
                         "--test-mishap is slow and our own measurement says missingness is not the "
                         "driver: restricting to zero-missingness variants moved AJ eta^2 only "
                         "0.982 -> 0.899. Run it once for completeness, not on the critical path")
    ap.add_argument("--min-mishap", type=int, default=100,
                    help="samples per callset per stratum needed to run --test-mishap there")
    ap.add_argument("--discordance", default=None,
                    help="per_variant_discordance.tsv; variants at >= --disc-rate are unioned in")
    ap.add_argument("--disc-rate", type=float, default=0.50)
    ap.add_argument("--ancs", nargs="+", default=None, help="default: every stratum with a .bed")
    ap.add_argument("--dx", nargs="+", default=None, help="default: every dx meeting --min-cell")
    a = ap.parse_args()

    qcd, work = Path(a.qc_dir), Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    grain = read_grain(a.grain)
    print(f"grain: {len(grain):,} samples")
    print(f"flag rule: |dAF| > {a.thresh} AND z > {a.zmin}, within (stratum x dx) cells "
          f"of >= {a.min_cell} per callset\n")

    ancs = a.ancs or sorted(p.name[len("cohort_"):-len("_qc.bed")]
                            for p in qcd.glob("cohort_*_qc.bed"))
    print(f"strata with a QC fileset: {' '.join(ancs)}\n")

    flagged = set()
    evaluated = set()
    universe = set()
    cells = []                                  # (anc, dx, csA, csB, nA, nB, n_shared, n_flag)
    by_pair = defaultdict(lambda: [0, 0, 0])    # pair -> [comparisons, shared, flagged]
    hwe_failed = set()
    hwe_rows = []                               # (anc, callset, n_controls, n_tested, n_fail)
    mishap_failed = set()
    mishap_rows = []                            # (anc, callset, n, n_ref, n_total_with_flanking)

    for anc in ancs:
        qc = qcd / f"cohort_{anc}_qc"
        if not Path(f"{qc}.bed").exists() or not Path(f"{qc}.fam").exists():
            print(f"{anc}: no fileset — skipped")
            continue

        fid = {}
        with open(f"{qc}.fam") as fh:
            for line in fh:
                t = line.split()
                fid[t[1]] = t[0]
        anc_vars = set()
        with open(f"{qc}.bim") as fh:
            for line in fh:
                anc_vars.add(line.split()[1])
        universe |= anc_vars

        # (dx, callset) -> members, restricted to samples actually in this fileset
        cell = defaultdict(list)
        for iid in fid:
            g = grain.get(iid)
            if g:
                cell[(g[2], g[0])].append(iid)

        dxs = a.dx or sorted({d for d, _ in cell})
        freq_cache = {}
        for dx in dxs:
            arms = {cs: ids for (d, cs), ids in cell.items()
                    if d == dx and len(ids) >= a.min_cell}
            allc = {cs: len(ids) for (d, cs), ids in cell.items() if d == dx}
            desc = ", ".join(f"{cs}={n:,}" for cs, n in sorted(allc.items())) or "none"
            if len(arms) < 2:
                print(f"{anc}/{dx:8}: {desc}  -> fewer than 2 callsets at >= {a.min_cell}, skipped")
                continue
            print(f"{anc}/{dx:8}: {desc}")

            for cs, ids in sorted(arms.items()):
                key = (anc, dx, cs)
                if key not in freq_cache:
                    kf = work / f"{anc}_{dx}_{cs}.keep"
                    kf.write_text("".join(f"{fid[i]}\t{i}\n" for i in ids))
                    stem = work / f"{anc}_{dx}_{cs}"
                    run(["plink2", "--bfile", str(qc), "--keep", str(kf),
                         "--freq", "--out", str(stem)])
                    freq_cache[key] = read_afreq(f"{stem}.afreq")

            names = sorted(arms)
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    x, y = names[i], names[j]
                    fx, fy = freq_cache[(anc, dx, x)], freq_cache[(anc, dx, y)]
                    shared = fx.keys() & fy.keys()
                    evaluated |= shared
                    bad = set()
                    for v in shared:
                        f1, n1 = fx[v]
                        f2, n2 = fy[v]
                        d = f1 - f2
                        if abs(d) <= a.thresh or n1 < 2 or n2 < 2:
                            continue
                        p = (f1 * n1 + f2 * n2) / (n1 + n2)      # pooled, allele-count weighted
                        var = p * (1 - p) * (1.0 / n1 + 1.0 / n2)
                        if var > 0 and abs(d) / var ** 0.5 > a.zmin:
                            bad.add(v)
                    flagged |= bad
                    pair = "|".join(sorted((x, y)))
                    by_pair[pair][0] += 1
                    by_pair[pair][1] += len(shared)
                    by_pair[pair][2] += len(bad)
                    cells.append((anc, dx, x, y, len(arms[x]), len(arms[y]), len(shared), len(bad)))
                    pct = 100 * len(bad) / len(shared) if shared else 0.0
                    print(f"    {x} vs {y}: {len(shared):,} shared, {len(bad):,} flagged ({pct:.3f}%)")

        # ── stage 2: per-callset HWE among CONTROLS of this stratum ──
        if a.hwe > 0:
            for cs, ids in sorted(cell.items()):
                dxv, csn = cs
                if dxv != a.control_dx or len(ids) < a.min_hwe_controls:
                    continue
                kf = work / f"hwe_{anc}_{csn}.keep"
                kf.write_text("".join(f"{fid[i]}\t{i}\n" for i in ids))
                stem = work / f"hwe_{anc}_{csn}"
                cmd = ["plink2", "--bfile", str(qc), "--keep", str(kf),
                       "--hwe", repr(a.hwe)]
                if not a.hwe_both_tails:
                    cmd.append("keep-fewhet")
                cmd += ["--write-snplist", "--out", str(stem)]
                run(cmd)
                passing = set()
                with open(f"{stem}.snplist") as fh:
                    for line in fh:
                        passing.add(line.strip())
                fail = anc_vars - passing
                hwe_failed |= fail
                hwe_rows.append((anc, csn, len(ids), len(anc_vars), len(fail)))
                print(f"    HWE {csn} controls (n={len(ids):,}): {len(fail):,} fail "
                      f"({100*len(fail)/len(anc_vars) if anc_vars else 0:.3f}%)")

        # ── stage 3: per-callset haplotype missingness (GenoTools' `haplotype`, plink1.9) ──
        # Runs on ALL samples of the callset, not just controls: --test-mishap keys on missingness
        # predicted by flanking haplotypes, which has no disease channel to worry about.
        if a.mishap > 0:
            by_cs = defaultdict(list)
            for (dxv, csn), ids in cell.items():
                by_cs[csn].extend(ids)
            for csn, ids in sorted(by_cs.items()):
                if len(ids) < a.min_mishap:
                    continue
                kf = work / f"mishap_{anc}_{csn}.keep"
                kf.write_text("".join(f"{fid[i]}\t{i}\n" for i in ids))
                stem = work / f"mishap_{anc}_{csn}"
                run(["plink", "--bfile", str(qc), "--keep", str(kf), "--maf", "0.05",
                     "--test-mishap", "--allow-no-sex", "--out", str(stem)])
                # .missing.hap: SNP HAPLOTYPE F_0 F_1 M_H1 M_H2 CHISQ P FLANKING
                # GenoTools excludes the reference SNP AND everything in FLANKING.
                bad, n_ref = set(), 0
                with open(f"{stem}.missing.hap") as fh:
                    hdr = fh.readline().split()
                    i_snp, i_p = hdr.index("SNP"), hdr.index("P")
                    i_fl = hdr.index("FLANKING") if "FLANKING" in hdr else None
                    for line in fh:
                        t = line.split()
                        if len(t) <= i_p or t[i_p] in ("NA", ""):
                            continue
                        try:
                            if float(t[i_p]) > a.mishap:
                                continue
                        except ValueError:
                            continue
                        n_ref += 1
                        bad.add(t[i_snp])
                        if i_fl is not None and len(t) > i_fl:
                            bad.update(v for v in t[i_fl].split("|") if v)
                bad &= anc_vars                 # FLANKING can name variants outside this QC set
                mishap_failed |= bad
                mishap_rows.append((anc, csn, len(ids), n_ref, len(bad)))
                print(f"    mishap {csn} (n={len(ids):,}): {n_ref:,} reference SNPs at "
                      f"p<={a.mishap:g}, {len(bad):,} variants with flanking")
        print()

    if not cells and not hwe_rows:
        sys.exit("no (stratum x dx) cell had two callsets above --min-cell, and no stratum had "
                 "enough controls in one callset to test HWE — nothing to filter on")

    n_af = len(flagged)
    n_hwe_new = len(hwe_failed - flagged)
    overlap = len(hwe_failed & flagged)
    flagged |= hwe_failed
    n_mishap_new = len(mishap_failed - flagged)
    mishap_overlap = len(mishap_failed & flagged)
    flagged |= mishap_failed

    # ── union in the duplicate-pair discordance tail ──
    n_disc = 0
    if a.discordance and Path(a.discordance).exists():
        with open(a.discordance) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            i_snp, i_rate = hdr.index("SNP"), hdr.index("rate")
            for line in fh:
                t = line.rstrip("\n").split("\t")
                if float(t[i_rate]) >= a.disc_rate:
                    flagged.add(t[i_snp])
                    n_disc += 1
    elif a.discordance:
        print(f"WARNING: {a.discordance} not found — discordance tail NOT unioned in")

    Path(a.out).write_text("".join(f"{v}\n" for v in sorted(flagged)))

    print("=" * 78)
    print("per-cell detail")
    print(f"{'anc':5} {'dx':9} {'callset A':10} {'callset B':10} {'nA':>6} {'nB':>6} "
          f"{'shared':>11} {'flag':>7} {'%':>7}")
    for anc, dx, x, y, na, nb, ns, nf in cells:
        print(f"{anc:5} {dx:9} {x:10} {y:10} {na:6,} {nb:6,} {ns:11,} {nf:7,} "
              f"{100*nf/ns if ns else 0:7.3f}")

    print("\n" + "=" * 78)
    print("BY CALLSET PAIR — this is the liftover verdict")
    print(f"{'pair':26} {'cells':>6} {'shared':>13} {'flagged':>9} {'%':>8}")
    for pair, (nc, ns, nf) in sorted(by_pair.items(), key=lambda kv: -kv[1][2]):
        print(f"{pair:26} {nc:6} {ns:13,} {nf:9,} {100*nf/ns if ns else 0:8.3f}")
    print("  divco_hs|wb_dwgs is BOTH-NATIVE. Comparable rate there => calling/mapping, no")
    print("  liftover fix exists, step 0 stays closed. Only wgs_harm pairs => liftover-specific.")

    if hwe_rows:
        print("\n" + "=" * 78)
        print(f"PER-CALLSET HWE among controls (p < {a.hwe:g}"
              f"{'' if a.hwe_both_tails else ', keep-fewhet'})")
        print(f"{'anc':5} {'callset':10} {'controls':>9} {'tested':>12} {'fail':>8} {'%':>7} "
              f"{'exp by chance':>14}")
        for anc, csn, nc, nt, nf in hwe_rows:
            # keep-fewhet tests one tail, so roughly half the nominal rate turns into exclusions
            exp = nt * a.hwe * (1.0 if a.hwe_both_tails else 0.5)
            print(f"{anc:5} {csn:10} {nc:9,} {nt:12,} {nf:8,} "
                  f"{100*nf/nt if nt else 0:7.3f} {exp:14,.0f}")
        print("  'exp by chance' is the false-positive count at this threshold. Excess over it is")
        print("  the real signal. A callset far above its expectation is the broken one.")
        print(f"\n  HWE failures also flagged by the frequency test : {overlap:,}")
        print(f"  HWE failures the frequency test missed         : {n_hwe_new:,}")
        print("  Overlap is a mechanism-based confirmation of an effect-based test: HWE knows")
        print("  nothing about the other callsets, only that heterozygosity is wrong.")

    if mishap_rows:
        print("\n" + "=" * 78)
        print(f"PER-CALLSET HAPLOTYPE MISSINGNESS (--test-mishap, p <= {a.mishap:g})")
        print(f"{'anc':5} {'callset':10} {'n':>8} {'ref SNPs':>10} {'with flanking':>14}")
        for anc, csn, n, nr, nb in mishap_rows:
            print(f"{anc:5} {csn:10} {n:8,} {nr:10,} {nb:14,}")
        print(f"\n  also caught by the earlier stages : {mishap_overlap:,}")
        print(f"  new                               : {n_mishap_new:,}")

    print("\n" + "=" * 78)
    print(f"flagged by AF concordance (|dAF| > {a.thresh}, z > {a.zmin}) : {n_af:,}")
    print(f"added by per-callset control HWE (p < {a.hwe:g})              : {n_hwe_new:,}")
    if mishap_rows:
        print(f"added by haplotype missingness (p <= {a.mishap:g})            : {n_mishap_new:,}")
    print(f"added by duplicate-pair discordance >= {a.disc_rate}          : {n_disc:,} rows")
    print(f"TOTAL to exclude                                        : {len(flagged):,}")
    print(f"variants some powered cell could evaluate               : {len(evaluated):,}")
    print(f"variants in some stratum's QC set                       : {len(universe):,}")
    un = len(universe - evaluated)
    print(f"  never evaluated, so kept by default                   : {un:,} "
          f"({100*un/len(universe) if universe else 0:.1f}%)")
    print(f"\nwrote {a.out}")

    # ── tripwire: is the filter reaching known biology? ──
    print("\n---- sentinel loci (expect empty; a hit means stop and look) ----")
    pos = [(v, parse_id(v)) for v in flagged]
    hits = 0
    for name, c, p, win in SENTINEL_LOCI:
        near = [v for v, cp in pos if cp and cp[0] == c and abs(cp[1] - p) <= win]
        if near:
            hits += len(near)
            print(f"  {name:16} {len(near)} flagged within {win//1000}kb:")
            for v in sorted(near)[:10]:
                print(f"      {v}")
    if not hits:
        print("  none — no flagged variant falls near a headline AD/PD locus.")
    else:
        print("\n  Inspect these before running step 6. Nothing is auto-whitelisted.")


if __name__ == "__main__":
    main()
