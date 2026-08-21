#!/usr/bin/env python3
"""Build the per-callset AF-concordance exclusion list for step 6 stage B.

Full rationale, thresholds and measured effect: METHODS.md §6.

THE LICENCE TO DELETE, restated here because this is the only script in the project allowed to
remove variants from the association set. Every comparison is made WITHIN a (stratum x dx) cell.
Inside a cell disease is constant, so a true disease effect cannot produce a between-callset
frequency difference and only technical breakage can. A filter comparing callsets over ALL samples
would have no such protection: wgs_harm/divco_hs supply the AD cases and wb_dwgs the PD cases, so
rs429358 differs sharply between them for an entirely real reason, and an all-samples filter would
delete the strongest true locus in the study. Do not widen the comparison beyond the cell.

There is NO automated check on that licence, deliberately. What replaces it is reading this script's
log — the BY CALLSET PAIR table and the per-cell HWE ratio table — and resolving anything that looks
wrong against the per-cell .afreq/.snplist in af_concordance/, indexing BY HEADER NAME (column 5 of
an .afreq is PROVISIONAL_REF?, not the frequency). Call rate is usually the tell, not frequency.
To name the gene a variant sits in: python3 scripts/gene_annot.py --at <chr:pos>

TIMING: step 6 applies this list at stage C, later in the SAME job that builds it here at stage B.
Nothing here can gate. Resolve anomalies before step 7 reads the association set.

TWO CHANNELS, unioned, with the overlap reported:
  frequency  plink --assoc on callset membership; flag needs BOTH |dAF| > --thresh AND z > --zmin
  HWE        within each callset, controls only, gated on exceeding chance expectation per cell

CONTROLS ONLY for the HWE channel, for the same reason as the cell design and the reason GenoTools
sets filter_controls=True: case ascertainment genuinely distorts HWE at a real disease locus, so
running it over an AD-case-heavy callset would flag APOE for an entirely real reason. keep-fewhet is
retained — mismapping produces het EXCESS, which is still caught, while het deficiency (which
pooling and substructure create artificially) is spared.

SCOPE: applied to the ASSOCIATION set, not just the PCA input. A variant mismapped badly enough to
bend PC1 also produces a spurious association in the test itself, where no PC adjustment reaches it.

READS exactly IID -> (source_callset, dx_detailed), from sample_annot.csv. Never a PC, never an
ancestry column — the stratum comes from which cohort_<ANC>_qc fileset a sample appears in.
--grain is accepted as an alias: the reader keys on column NAMES and the grain carries the same two,
so either file gives an identical result. That equivalence is the regression test for the split.

GUARDRAIL: human-run. Reads the id-bearing sample annotation and genotypes via plink2, writes the
exclusion list, prints ONLY aggregate counts.
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

# NO SENTINEL-LOCUS TRIPWIRE HERE, DELIBERATELY — removed 2026-08-20; do not reintroduce it. A
# hand-picked gene list scrutinises a handful of deletions while thousands get none, and lets a
# NEGATIVE be printable that no evidence supports. Per-variant is the honest form of the question:
# the BY CALLSET PAIR table plus gene_annot.py --at. Full argument in PROJECT_LOG.md 2026-08-20.
def run(cmd):
    """plink2 with output captured. stderr=STDOUT because capture_output alone hides plink's
    errors entirely (documented gotcha) — a silent failure here would look like a real result."""
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        sys.exit(f"FAILED: {' '.join(cmd)}\n{r.stdout[-4000:]}")
    return r.stdout


def read_annot(path):
    """-> {iid: (source_callset, dx_detailed)}.

    Keyed on column NAMES, not positions. That is what lets the same reader take either
    sample_annot.csv (the PC-free half, written by §12a before any genotype step) or
    analysis_grain.csv (which carries the same two names plus ancestry, sex and the PCs) — so the
    two can be diffed against each other to prove the split changed nothing. The old positional
    reader could not: it hardcoded the grain's column order, which is also why the dependency on
    the grain looked structural when it was only a file-layout accident.

    Anything else in the file is ignored on purpose. The PCs go stale the moment step 6 re-runs
    and nothing here depends on them; ancestry is never consulted because the stratum is set by
    which fileset a sample is in."""
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

    THE TEST IS PLINK'S, NOT OURS. `plink --assoc` with callset membership as the phenotype is
    the standard 1-df allelic chi-square on the 2x2 allele-count table. This used to be computed
    here in Python — a two-sample test of proportions with pooled variance, `z = |f1-f2| /
    sqrt(p(1-p)(1/n1+1/n2))` — which is the same statistic by algebra (z^2 = chi^2), but a
    reviewer has to take a hand-rolled implementation on trust and can simply recognise this one.

    VERIFIED IDENTICAL before the swap, 2026-08-20, on the largest cell (EUR/AD, divco_hs=121 vs
    wgs_harm=657, 7,538,809 variants): the Python rule flagged 1,698, `--assoc` at CHISQ > 25 with
    the same |dAF| floor flagged 1,698, and the symmetric difference was **0 in both directions**.
    So this swap changed no result — it removed a statistic from the codebase.

    Both criteria are kept, and the second is not optional: CHISQ alone scales with n, so at
    EUR's thousands a negligible difference clears any threshold while at CAH's hundreds a large
    one may not, and the filter would then mean something different in every stratum (the same
    argument 07_gwas.sh:49 makes for its differential-missingness filter). Concretely, at MAF ~20%
    in this cell z > 5 already requires |dAF| ~ 0.14 so the floor never binds, while at MAF ~2% it
    is reachable at ~0.035 and the floor is what stops the list filling with low-frequency noise.
    The floor is ABSOLUTE, deliberately: a 0.035 gap at MAF 2% is a large *relative* discordance
    and is kept anyway, because a filter licensed to delete should err toward keeping.

    --assoc is plink1.9 only (plink2 dropped it for --glm, which is logistic regression on dosage:
    asymptotically equivalent, not identical, so it would have moved the flags and cost us the
    verification above). Both shell wrappers therefore load MOD_PLINK1 unconditionally.
    """
    # The .pheno file IS the arm membership record (2=x, 1=y, one per line) — read it when a
    # flagged variant needs resolving by hand. `--assoc` needs it and nothing else.
    ph = work / f"{anc}_{dx}_{x}__vs__{y}.pheno"
    # plink1.9 case/control coding: 2 = case, 1 = control. Which arm is which is arbitrary —
    # |F_A - F_U| and CHISQ are both symmetric.
    ph.write_text("".join(f"{fid[i]}\t{i}\t2\n" for i in ids_x)
                  + "".join(f"{fid[i]}\t{i}\t1\n" for i in ids_y))
    stem = work / f"{anc}_{dx}_{x}__vs__{y}"
    run([a.plink1, "--bfile", str(qc), "--keep", str(ph), "--pheno", str(ph),
         "--assoc", "--allow-no-sex", "--out", str(stem)])

    # THRESHOLD ON P, NOT ON CHISQ. The three forms are algebraically the same rule, but .assoc
    # prints CHISQ to four SIGNIFICANT figures, so a true 25.005 lands on disk as "25" and
    # `25.0 > 25` is false. P is exponential and carries ~1000x the resolution at the boundary.
    # Switching back to CHISQ silently drops boundary variants — it cost one on the first run.
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
                # A1 is plink's minor allele, not necessarily ALT — but |F_A - F_U| is invariant
                # to which allele is counted, so the effect-size floor transfers unchanged.
                daf = abs(float(t[i_fa]) - float(t[i_fu]))
            except ValueError:
                continue                          # plink writes NA where an arm is monomorphic
            if pval < p_max and daf > a.thresh:
                bad.add(snp)
    return tested, bad


# RESOLVING A FLAGGED VARIANT: CALL RATE is the tell, not frequency. Run --freq per arm over the
# af_concordance/ intermediates, index BY HEADER NAME (plink2 puts PROVISIONAL_REF? at column 5,
# which is what hid this the first time), and look at OBS_CT. CR1 was excluded because divco_hs
# called it in 156 of 242 alleles while the same samples were 242/242 at the LRRK2 site.


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
    # BOTH BINARIES, BY ABSOLUTE PATH — they cannot both be on PATH. On biowulf plink and plink2
    # are one module family, so loading plink/1.9 UNLOADS plink/6-alpha and leaves a
    # non-executable plink2 earlier on PATH; a bare "plink2" then dies with PermissionError.
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

    # Fail here, not 40 minutes in: the PATH clash above surfaces as a PermissionError from deep
    # inside subprocess, three stages late.
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

        # (dx, callset) -> members, restricted to samples in this fileset. The stratum is `anc`,
        # the loop variable over filesets — never a column. An IID absent from the annotation is
        # SKIPPED, not pooled into a blank cell: an unlabelled sample has no disease to hold
        # constant, so it cannot join a comparison whose validity rests on disease being constant.
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
                # A fixed threshold rejects a PREDICTABLE number of variants when nothing is wrong:
                # N*alpha, halved because keep-fewhet takes one tail. At N~8M and alpha=1e-4 that is
                # ~400. So a cell's rejection count is (false positives + true positives) with the
                # first term known in advance, and E/O is exactly the Benjamini-Hochberg FDR
                # estimate for its rejection set.
                #
                # This cell therefore only votes if it clears its OWN expectation. There is no
                # multiplier to justify: the bar is "more rejections than chance explains".
                # Measured 2026-08-20 (job 27857727): EUR/wgs_harm 1,680 vs 377 = 4.5x — real, and
                # from the SMALLEST control sample of the three, so it is not a power artifact;
                # EUR/wb_dwgs 132 vs 377 = 0.35x; AJ/wb_dwgs 97 vs 400 = 0.24x. Under the old
                # unconditional union those two contributed ~229 variants indistinguishable from
                # null rejections — one of them inside LRRK2, deleted from EVERY stratum on the
                # strength of a cell showing no excess at all. Below-chance is not suspicious, it is
                # underpowered: the exact test is discrete and conservative at n in the hundreds.
                #
                # Why routine pruning never needs this: step 6 stage A already ran --hwe 1e-6 per
                # stratum, where the null contributes ~4 in 7.5M. This scan is at 1e-4 (100x the
                # null burden), PER CALLSET, and its result applies ACROSS strata — it gives up both
                # properties that let a conventional prune ignore the null, in exchange for callset
                # attribution. HWE_REQUIRE_EXCESS=0 restores the old unconditional union.
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
        # exp/ratio are carried from the decision point rather than recomputed here: the number
        # that gates and the number that prints must be the same one.
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
    print("PROJECT_LOG.md so the next run does not re-derive it.")

    # ── provenance, beside the list ──
    # The list itself must stay a bare ID list for `plink2 --exclude`, so it cannot carry a
    # header — which is why "what settings produced this 4,415-variant file?" was unanswerable
    # from disk, the same failure class as the 2026-08-17 stale-list incident.
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
