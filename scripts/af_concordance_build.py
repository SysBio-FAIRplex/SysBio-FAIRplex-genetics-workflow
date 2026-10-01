#!/usr/bin/env python3
"""Step 6 stage B: build the per-callset AF-concordance exclusion list.

Rationale, thresholds and measured effect: METHODS.md §6.

Every comparison is made within a (stratum × dx) cell, where disease is constant, so a
between-callset frequency gap can only be technical. Never widen it beyond the cell: across all
samples APOE differs between callsets for a real reason and would be deleted.

Two channels, unioned: frequency (plink --assoc on callset; |dAF| > --thresh AND z > --zmin) and
HWE (within callset, controls only, gated on excess over chance).

Stage C applies the list later in the same job, so nothing here gates; read the log's BY CALLSET
PAIR and HWE tables before step 7 (README, Run order). Reads IID → (source_callset, dx_detailed)
from sample_annot.csv; prints aggregate counts only.
"""
from collections import defaultdict
from pathlib import Path
import argparse
import csv
import math
import os
import shutil
import subprocess
import sys

def run(cmd):
    """Run plink with stderr merged into stdout, so a failure shows its error instead of passing."""
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        sys.exit(f"FAILED: {' '.join(cmd)}\n{r.stdout[-4000:]}")
    return r.stdout


def read_annot(path):
    """-> {iid: (source_callset, dx_detailed)}.

    Keyed on column names, so sample_annot.csv and analysis_grain.csv both work and must agree.
    Everything else (PCs, ancestry) is ignored: the stratum comes from the fileset, not a column."""
    out = {}
    with open(path, newline="") as fh:
        rdr = csv.reader(fh)
        hdr = next(rdr, None)
        if not hdr:
            sys.exit(f"{path}: empty")
        cols = {c.strip(): i for i, c in enumerate(hdr)}
        missing = [c for c in ("IID", "source_callset", "dx_detailed") if c not in cols]
        if missing:
            sys.exit(f"{path}: missing required column(s) {', '.join(missing)}. "
                     f"Expected a sample_annot.csv or analysis_grain.csv; found: {', '.join(hdr)}")
        i_id, i_cs, i_dx = cols["IID"], cols["source_callset"], cols["dx_detailed"]
        for row in rdr:
            if len(row) > max(i_id, i_cs, i_dx) and row[i_id].strip():
                out[row[i_id].strip()] = (row[i_cs].strip(), row[i_dx].strip())
    return out


def assoc_pair(qc, work, fid, anc, dx, x, ids_x, y, ids_y, a):
    """One callset-pair comparison inside a (stratum x dx) cell -> (tested, flagged).

    plink --assoc with callset membership as the phenotype: the standard 1-df allelic chi-square.
    A variant is flagged only if it clears BOTH the significance floor and the absolute |dAF|
    floor; CHISQ alone scales with n, so without the floor the filter would mean something
    different in every stratum. --assoc is plink1.9 only (plink2's --glm is not identical), so
    both shell wrappers load MOD_PLINK1 unconditionally.
    """
    # Arm membership (2 = x, 1 = y); also the record to check when resolving a flagged variant.
    # Which arm is 2 is arbitrary: |F_A - F_U| and CHISQ are symmetric.
    ph = work / f"{anc}_{dx}_{x}__vs__{y}.pheno"
    ph.write_text("".join(f"{fid[i]}\t{i}\t2\n" for i in ids_x)
                  + "".join(f"{fid[i]}\t{i}\t1\n" for i in ids_y))
    stem = work / f"{anc}_{dx}_{x}__vs__{y}"
    run([a.plink1, "--bfile", str(qc), "--keep", str(ph), "--pheno", str(ph),
         "--assoc", "--allow-no-sex", "--out", str(stem)])

    # Threshold on P, not CHISQ: .assoc prints CHISQ to 4 significant figures, so a true 25.005
    # reads as "25" and fails `> 25`. P keeps the resolution at the boundary.
    p_max = math.erfc(a.zmin / 2 ** 0.5)
    tested, bad = set(), set()
    with open(f"{stem}.assoc") as fh:
        hdr = fh.readline().split()
        i_snp, i_fa, i_fu, i_p = (hdr.index("SNP"), hdr.index("F_A"),
                                  hdr.index("F_U"), hdr.index("P"))
        for line in fh:
            t = line.split()
            if len(t) <= i_p:
                continue
            snp = t[i_snp]
            tested.add(snp)
            try:
                pval = float(t[i_p])
                # A1 is plink's minor allele, not necessarily ALT; |F_A - F_U| is invariant to that.
                daf = abs(float(t[i_fa]) - float(t[i_fu]))
            except ValueError:
                continue                          # plink writes NA where an arm is monomorphic
            if pval < p_max and daf > a.thresh:
                bad.add(snp)
    return tested, bad


# Resolving a flagged variant: call rate is usually the tell, not frequency. Run --freq per arm on
# the af_concordance/ intermediates and read OBS_CT, indexing .afreq BY HEADER NAME (column 5 is
# PROVISIONAL_REF?).


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc-dir", required=True,
                    help="step 6 stage A's UNFILTERED per-ancestry filesets "
                         "(by_ancestry_qc/unfiltered). Must be unfiltered: deriving the list "
                         "from a directory that already has a previous list applied would "
                         "pre-filter this filter's own input, and it would look clean either way")
    ap.add_argument("--annot", "--grain", dest="annot", required=True,
                    help="sample_annot.csv (IID, source_callset, dx_detailed). analysis_grain.csv "
                         "also works — same column names — and the two must agree")
    ap.add_argument("--out", required=True, help="exclusion list to write")
    ap.add_argument("--work", required=True)
    # Both binaries by absolute path: plink and plink2 are one module family on biowulf, so loading
    # one unloads the other and a bare name can resolve to a non-executable file.
    ap.add_argument("--plink1", default="plink",
                    help="plink1.9 binary — used for --assoc and --test-mishap. Pass an absolute "
                         "path; the shell wrappers resolve it while that module is loaded")
    ap.add_argument("--plink2", default="plink2",
                    help="plink2 binary — used for --freq, --hwe. Absolute path, as above")
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
    ap.add_argument("--hwe-require-excess", dest="hwe_require_excess",
                    action="store_true", default=True,
                    help="(default) a stratum x callset cell contributes HWE exclusions only if "
                         "its rejection count exceeds the number expected by chance at --hwe. "
                         "E/O is the BH false-discovery estimate for that cell's rejections, so a "
                         "cell at or below 1.0x has no attributable signal to contribute")
    ap.add_argument("--no-hwe-require-excess", dest="hwe_require_excess",
                    action="store_false",
                    help="union every cell's HWE rejections in regardless of excess — the "
                         "pre-2026-08-20 behaviour, kept so the 4,415-variant list is reproducible")
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

    # Fail fast: a bad path would otherwise surface as a PermissionError stages later.
    for label, exe in (("--plink1", a.plink1), ("--plink2", a.plink2)):
        if not (os.path.isabs(exe) or shutil.which(exe)):
            sys.exit(f"{label}={exe} is not executable and not on PATH. plink and plink2 are ONE "
                     f"module family here, so PATH can hold only one — pass absolute paths "
                     f"(the shell wrappers resolve them while each module is loaded).")
        if os.path.isabs(exe) and not os.access(exe, os.X_OK):
            sys.exit(f"{label}={exe} is not executable.")

    qcd, work = Path(a.qc_dir), Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    annot = read_annot(a.annot)
    print(f"annot: {len(annot):,} samples from {a.annot}")
    print(f"flag rule: plink --assoc P < {math.erfc(a.zmin / 2 ** 0.5):.4g} "
          f"(z > {a.zmin:g}, i.e. chi-square > {a.zmin**2:g}) AND |dAF| > {a.thresh}, "
          f"within (stratum x dx) cells of >= {a.min_cell} per callset\n")

    ancs = a.ancs or sorted(p.name[len("cohort_"):-len("_qc.bed")]
                            for p in qcd.glob("cohort_*_qc.bed"))
    print(f"strata with a QC fileset: {' '.join(ancs)}\n")

    flagged = set()
    evaluated = set()
    universe = set()
    cells = []                                  # (anc, dx, csA, csB, nA, nB, n_shared, n_flag)
    by_pair = defaultdict(lambda: [0, 0, 0])    # pair -> [comparisons, shared, flagged]
    hwe_failed = set()
    hwe_withheld = set()                        # failures from cells with no excess over chance
    hwe_rows = []      # (anc, callset, n_controls, n_tested, n_fail, exp, ratio, voted)
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

        # (dx, callset) -> members of this fileset; the stratum is the fileset itself. Unannotated
        # samples are skipped: with no disease to hold constant they cannot join a cell.
        cell = defaultdict(list)
        n_unannotated = 0
        for iid in fid:
            g = annot.get(iid)
            if g:
                cell[(g[1], g[0])].append(iid)
            else:
                n_unannotated += 1
        if n_unannotated:
            print(f"{anc}: {n_unannotated:,} of {len(fid):,} samples not in the annotation — "
                  f"excluded from every cell")

        dxs = a.dx or sorted({d for d, _ in cell})
        for dx in dxs:
            arms = {cs: ids for (d, cs), ids in cell.items()
                    if d == dx and len(ids) >= a.min_cell}
            allc = {cs: len(ids) for (d, cs), ids in cell.items() if d == dx}
            desc = ", ".join(f"{cs}={n:,}" for cs, n in sorted(allc.items())) or "none"
            if len(arms) < 2:
                print(f"{anc}/{dx:8}: {desc}  -> fewer than 2 callsets at >= {a.min_cell}, skipped")
                continue
            print(f"{anc}/{dx:8}: {desc}")

            names = sorted(arms)
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    x, y = names[i], names[j]
                    shared, bad = assoc_pair(qc, work, fid, anc, dx,
                                             x, arms[x], y, arms[y], a)
                    evaluated |= shared
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
                cmd = [a.plink2, "--bfile", str(qc), "--keep", str(kf),
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
                # A fixed threshold rejects ~N*alpha variants when nothing is wrong (halved:
                # keep-fewhet is one tail), so E/O is the BH false-discovery estimate for this
                # cell, and the cell contributes only if it clears its own chance expectation
                # (METHODS.md §6.3). Stage A's --hwe 1e-6 never needs this; this scan is 100x
                # looser, per callset, and its result applies across strata.
                exp = len(anc_vars) * a.hwe * (1.0 if a.hwe_both_tails else 0.5)
                ratio = len(fail) / exp if exp else 0.0
                voted = (not a.hwe_require_excess) or len(fail) > exp
                if voted:
                    hwe_failed |= fail
                else:
                    hwe_withheld |= fail
                hwe_rows.append((anc, csn, len(ids), len(anc_vars), len(fail), exp, ratio, voted))
                print(f"    HWE {csn} controls (n={len(ids):,}): {len(fail):,} fail "
                      f"({100*len(fail)/len(anc_vars) if anc_vars else 0:.3f}%), {ratio:.2f}x chance"
                      + ("" if voted else "  -> NOT unioned in (no excess over chance)"))

        # ── stage 3: per-callset haplotype missingness (GenoTools' `haplotype`, plink1.9) ──
        # All samples, not just controls: missingness predicted by flanking haplotypes has no
        # disease channel.
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
                run([a.plink1, "--bfile", str(qc), "--keep", str(kf), "--maf", "0.05",
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
              f"{'exp by chance':>14} {'ratio':>7}  {'used?':6}")
        # exp/ratio come from the decision point, so the number that gates is the number printed.
        for anc, csn, nc, nt, nf, exp, ratio, voted in hwe_rows:
            print(f"{anc:5} {csn:10} {nc:9,} {nt:12,} {nf:8,} "
                  f"{100*nf/nt if nt else 0:7.3f} {exp:14,.0f} {ratio:7.2f}  "
                  f"{'yes' if voted else 'WITHHELD':6}")
        print("  'exp by chance' is the false-positive count at this threshold; E/O is the BH")
        print("  false-discovery estimate for a cell's rejections. A callset far above its own")
        print("  expectation is the broken one — and one at or below 1.0x has nothing to")
        print("  attribute, so it contributes nothing (--no-hwe-require-excess to override).")
        if hwe_withheld:
            n_w = len(hwe_withheld - flagged)
            print(f"  WITHHELD from cells with no excess: {len(hwe_withheld):,} rejections, "
                  f"{n_w:,} of which are in no other channel and so are NOT excluded.")
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
    print(f"flagged by AF concordance (--assoc P < {math.erfc(a.zmin / 2 ** 0.5):.3g}, "
          f"|dAF| > {a.thresh}) : {n_af:,}")
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
    print("\nThis list is applied by stage C of step 6, LATER IN THIS SAME JOB. Read the BY CALLSET")
    print("PAIR table above before step 7 consumes the association set; resolve anything that looks")
    print("wrong against the per-cell intermediates in af_concordance/ (index .afreq by HEADER NAME")
    print("— column 5 is PROVISIONAL_REF?, not a frequency), and record the verdict in")
    print("METHODS.md §6 so the next run does not re-derive it.")

    # ── provenance, beside the list ──
    # The list must stay a bare ID list for `plink2 --exclude`, so its settings are recorded here.
    prov = Path(f"{a.out}.provenance.txt")
    prov.write_text(
        f"exclusion list : {a.out}\n"
        f"variants       : {len(flagged):,}\n"
        f"built          : job {os.environ.get('SLURM_JOB_ID', 'interactive')} on "
        f"{os.environ.get('SLURMD_NODENAME', 'unknown host')}\n"
        f"qc-dir         : {a.qc_dir}\n"
        f"annot          : {a.annot}\n"
        f"\nfrequency channel : plink --assoc (1-df allelic chi-square), "
        f"P < {math.erfc(a.zmin / 2 ** 0.5):.4g} (z > {a.zmin:g}) AND |dAF| > {a.thresh}\n"
        f"  min-cell        : {a.min_cell} per callset per (stratum x dx) cell\n"
        f"  flagged         : {n_af:,}\n"
        f"HWE channel       : p < {a.hwe:g}"
        f"{'' if a.hwe_both_tails else ', keep-fewhet'}, controls only, "
        f"min {a.min_hwe_controls} controls\n"
        f"  require-excess  : {a.hwe_require_excess}"
        f"{'' if a.hwe_require_excess else '  (PRE-2026-08-20 BEHAVIOUR)'}\n"
        f"  added           : {n_hwe_new:,}   withheld: {len(hwe_withheld):,}\n"
        f"discordance       : rate >= {a.disc_rate} -> {n_disc:,} rows\n"
        f"mishap            : {'off' if a.mishap <= 0 else f'p <= {a.mishap:g}'}\n"
        "\nper-cell HWE (anc callset controls tested fail exp ratio used)\n"
        + "".join(f"  {r[0]:5} {r[1]:10} {r[2]:7,} {r[3]:12,} {r[4]:8,} "
                  f"{r[5]:8,.0f} {r[6]:6.2f}  {'yes' if r[7] else 'WITHHELD'}\n"
                  for r in hwe_rows))
    print(f"\nprovenance -> {prov}")


if __name__ == "__main__":
    main()
