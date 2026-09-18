#!/usr/bin/env python3
"""Re-derive every numeric claim in METHODS.md from artifacts, and diff against the doc.

Read-only. Nothing here writes, deletes, or submits. It exists because the numbers in
METHODS.md were typed by hand from six different artifacts, and three of them had drifted:
a denominator that appeared nowhere else, a multiplier contradicted by the summary CSV, and
a lambda range quoted over the wrong contrast set.

Each claim below carries the value METHODS.md ASSERTS next to the value this script DERIVES.
Three outcomes, and the third is the one that matters:

    OK           derived == asserted
    MISMATCH     the doc is wrong, or the artifact changed under it
    UNAVAILABLE  the source file is not on this machine, so the check DID NOT RUN

UNAVAILABLE is never silent and never counts as a pass (CLAUDE.md rule 3). Six claims need
cluster artifacts; on a laptop they report UNAVAILABLE with the path that was missing.

    python3 review/methods_numbers.py                  # everything reachable from here
    python3 review/methods_numbers.py --strict         # exit 1 on MISMATCH or UNAVAILABLE
    python3 review/methods_numbers.py --only eta2      # one group
    python3 review/methods_numbers.py --discordance    # the one claim that is a computation

Groups: doc (METHODS.md against itself), manifest, eta2, gwas, cluster.
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

OK, MISMATCH, UNAVAILABLE = "OK", "MISMATCH", "UNAVAILABLE"

# Every claim METHODS.md makes that a file on disk can settle. `asserted` is transcribed from
# the doc; change it here and in METHODS.md together, or this script stops being a check.
results = []


def record(section, claim, asserted, derived, tol=0):
    """Compare an asserted value to a derived one. Numbers within `tol`; everything else ==."""
    if derived is None:
        status = UNAVAILABLE
    elif isinstance(asserted, float) or isinstance(derived, float):
        status = OK if abs(float(asserted) - float(derived)) <= tol else MISMATCH
    else:
        status = OK if asserted == derived else MISMATCH
    results.append((section, claim, asserted, derived, status))


def unavailable(section, claim, asserted, why):
    results.append((section, claim, asserted, f"— {why}", UNAVAILABLE))


def read_csv(path):
    return pd.read_csv(path) if path.exists() else None


def wc(path):
    if not path.exists():
        return None
    with open(path, "rb") as fh:
        return sum(1 for _ in fh)


def first_existing(*paths):
    """The first path that exists, else None.

    Two locations per artifact, deliberately. The pipeline writes under `data/merged/...`
    (`config.sh` and `07_gwas.sh` agree on that), and copies get pulled down to `results/`
    for plotting on a laptop. Checking only `results/` made every §7 claim UNAVAILABLE on the
    cluster — where the real file is — while passing on the laptop copy. Canonical first.
    """
    for p in paths:
        if p.exists():
            return p
    return None


def norm_stems(root):
    """The four step-2 output stems, read from config.sh — the one place that resolves them.

    `merge_list.txt` cannot serve this purpose and an earlier version of this script was wrong
    to try. It holds only the three SECONDARY stems (`wgs_harm` enters the merge via --bfile,
    `03_merge.sh:26`), so a sum over it is structurally short by one callset; and it holds
    absolute paths written at run time, which went stale when the project root moved on
    2026-09-17. Sourcing config.sh is immune to both.
    """
    keys = ("NORM_WGS", "NORM_WB", "NORM_DC", "NORM_BR")
    cfg = root / "config.sh"
    if not cfg.exists():
        return None
    cmd = f'source "{cfg}" && ' + " && ".join(f'echo "${k}"' for k in keys)
    try:
        out = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    vals = [ln.strip() for ln in out.stdout.splitlines() if ln.strip()]
    if out.returncode != 0 or len(vals) != len(keys):
        return None
    return dict(zip(keys, vals))


def id_collisions(stems):
    """FID+IID pairs carried by more than one callset. plink1.9 treats these as one individual
    and fuses them into a single merged sample, so they are exactly the 13,428 → 13,334 gap."""
    seen, dup = set(), set()
    for stem in stems.values():
        fam = Path(stem + ".fam")
        if not fam.exists():
            return None
        ids = set()
        with open(fam) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 2:
                    ids.add((parts[0], parts[1]))
        dup |= seen & ids
        seen |= ids
    return len(dup)


# ── manifest: who is in the study ────────────────────────────────────────────────────────

def check_manifest(root):
    # config.sh:117 — the pipeline writes this under by_ancestry_qc/; results/ holds a copy.
    path = first_existing(
        root / "data" / "merged" / "by_ancestry_qc" / "retained_samples_manifest.csv",
        root / "results" / "retained_samples_manifest.csv",
    )
    m = read_csv(path) if path else None
    if m is None:
        for claim, asserted in [("§4 retained samples", 12495),
                                ("§5 EUR n", 10135), ("§5 AJ n", 1518), ("§5 AAC n", 226),
                                ("§5 AFR n", 191), ("§5 AMR n", 185), ("§5 CAH n", 106),
                                ("§6.3 wgs_harm EUR samples", 1540)]:
            unavailable("manifest", claim, asserted,
                        "missing retained_samples_manifest.csv under "
                        "data/merged/by_ancestry_qc/ or results/")
        return

    record("§4", "retained samples", 12495, len(m))
    for anc, asserted in [("EUR", 10135), ("AJ", 1518), ("AAC", 226),
                          ("AFR", 191), ("AMR", 185), ("CAH", 106)]:
        record("§5", f"{anc} n with PCs", asserted, int((m.ancestry == anc).sum()))

    # A genome sequenced in two callsets carries both, pipe-joined, so membership is counted
    # by splitting rather than by grouping on the raw string.
    eur = m[m.ancestry == "EUR"]
    harm = sum(1 for s in eur.source_callset if "wgs_harm" in str(s).split("|"))
    record("§6.3", "wgs_harm share of EUR (HWE denominator)", 1540, harm)

    record("§4", "excluded = 13,334 − retained", 839, 13334 - len(m))
    record("§4", "excluded = relatives + duplicates + sex", 839, 499 + 319 + 21)


# ── eta2: the headline result ────────────────────────────────────────────────────────────

def check_eta2(root):
    per_stratum = read_csv(root / "results" / "pca" / "af_filter_effect_gated4187.csv")
    per_pc = read_csv(root / "results" / "pca" / "af_filter_effect_gated4187_per_pc.csv")

    table = {  # §6.4: stratum -> (before, worst PC before, after, worst PC after)
        "EUR": (0.757, "PC2", 0.036, "PC6"),
        "AJ": (0.984, "PC1", 0.962, "PC1"),
        "AAC": (0.877, "PC2", 0.876, "PC2"),
        "AFR": (0.708, "PC1", 0.713, "PC1"),
        "CAH": (0.309, "PC4", 0.199, "PC4"),
        "AMR": (0.203, "PC8", 0.170, "PC8"),
    }
    if per_stratum is None:
        for anc, (b, _, a, _) in table.items():
            unavailable("§6.4", f"{anc} max eta2 before/after", f"{b}/{a}",
                        "missing af_filter_effect_gated4187.csv")
    else:
        idx = per_stratum.set_index("ancestry")
        for anc, (b, bpc, a, apc) in table.items():
            row = idx.loc[anc]
            record("§6.4", f"{anc} max eta2 before", b, round(float(row.max_eta2_before), 3), tol=0.0005)
            record("§6.4", f"{anc} worst PC before", bpc, row.worst_pc_before)
            record("§6.4", f"{anc} max eta2 after", a, round(float(row.max_eta2_after), 3), tol=0.0005)
            record("§6.4", f"{anc} worst PC after", apc, row.worst_pc_after)

    if per_pc is None:
        unavailable("§6.4", "EUR profile collapse (all 10 PCs <= 0.036)", True,
                    "missing af_filter_effect_gated4187_per_pc.csv")
        return

    pcs = [f"PC{i}" for i in range(1, 11)]
    idx = per_pc.set_index(["generation", "ancestry"])
    before, after = idx.loc[("before", "EUR")], idx.loc[("after", "EUR")]

    # The claim is that the whole profile collapsed, not that the peak moved: every component
    # after the filter sits at or below what used to be the SIXTH-largest value.
    record("§6.4", "EUR PCs above 0.036 after filtering", 0,
           int(sum(1 for pc in pcs if float(after[pc]) > 0.0365)))
    record("§6.4", "EUR PC2 after filtering", 0.011, round(float(after.PC2), 3), tol=0.0005)
    for pc, asserted in [("PC1", 0.080), ("PC2", 0.757), ("PC3", 0.155), ("PC7", 0.038)]:
        record("§6.4", f"EUR {pc} before filtering", asserted, round(float(before[pc]), 3), tol=0.0005)

    # Small-strata instability: a filter that only REMOVES variants moved these upward.
    for anc, pc, b, a in [("CAH", "PC6", 0.098, 0.123), ("AFR", "PC6", 0.005, 0.050)]:
        record("§6.4", f"{anc} {pc} before (noise floor)", b,
               round(float(idx.loc[("before", anc)][pc]), 3), tol=0.0005)
        record("§6.4", f"{anc} {pc} after (noise floor)", a,
               round(float(idx.loc[("after", anc)][pc]), 3), tol=0.0005)


# ── gwas: contrasts, inflation, differential missingness ─────────────────────────────────

def check_gwas(root):
    # NOT by_ancestry_qc_pre_pcfix_*/ — that is a superseded generation still on the cluster.
    path = first_existing(
        root / "data" / "merged" / "by_ancestry_qc" / "gwas" / "gwas_summary.csv",
        root / "results" / "gwas" / "gwas_summary.csv",
    )
    g = read_csv(path) if path else None
    if g is None:
        unavailable("§7", "all contrast-level claims", "44 ran / 17 viable",
                    "missing gwas_summary.csv under data/merged/by_ancestry_qc/gwas/ "
                    "or results/gwas/")
        return
    print(f"  §7 source: {path.relative_to(root)}")

    viable = g[g.viable_ge100 == "yes"]
    record("§7", "contrasts run", 44, len(g))
    record("§7", "viable contrasts", 17, len(viable))
    record("§7", "viable EUR", 14, int((viable.ancestry == "EUR").sum()))
    record("§7", "viable AJ", 3, int((viable.ancestry == "AJ").sum()))

    record("§7", "lambda_GC min over VIABLE", 1.0175, float(viable.lambda_gc.min()), tol=1e-9)
    record("§7", "lambda_GC max over VIABLE", 1.0549, float(viable.lambda_gc.max()), tol=1e-9)
    record("§7", "contrast holding the viable max", "PD_vs_AD",
           viable.loc[viable.lambda_gc.idxmax(), "contrast"])

    nonviable = g[g.viable_ge100 == "no"]
    record("§7", "lambda_GC max over NON-viable", 1.3247, float(nonviable.lambda_gc.max()), tol=1e-9)
    record("§7", "contrasts returning lambda = 0", 3, int((g.lambda_gc == 0).sum()))
    record("§7", "all lambda = 0 contrasts are AJ", True,
           bool((g[g.lambda_gc == 0].ancestry == "AJ").all()))

    def one(anc, contrast, col):
        row = g[(g.ancestry == anc) & (g.contrast == contrast)]
        return None if row.empty else row.iloc[0][col]

    # §10, differential missingness. PD_vs_DLB is the largest WITHIN-program removal; the
    # doc used to call it "six times any cross-program contrast", which PD_vs_AD contradicts.
    for contrast, asserted in [("PD_vs_AD", 387207), ("PD_vs_DLB", 353068), ("AD_vs_DLB", 65245),
                               ("AD_vs_control", 63393), ("control_amppd_vs_control_ampad", 62344),
                               ("PD_amppd_vs_control_amppd", 5814)]:
        record("§10", f"EUR {contrast} diffmiss removals", asserted,
               int(one("EUR", contrast, "n_diffmiss_excluded")))

    within = g[(g.ancestry == "EUR") & (g.confound_tag == "within_cohort")]
    record("§10", "PD_vs_DLB is the largest within-program removal", "PD_vs_DLB",
           within.loc[within.n_diffmiss_excluded.idxmax(), "contrast"])
    record("§10", "EUR PD_vs_DLB lambda", 1.0402, float(one("EUR", "PD_vs_DLB", "lambda_gc")), tol=1e-9)
    record("§10", "EUR PD_vs_AD arm size (the ~55 br_dsnwgs share)", 2595,
           int(one("EUR", "PD_vs_AD", "n_case")))

    # §8, corroboration: high-signal contrasts came through with zero control-flagged hits.
    for contrast, asserted in [("AD_vs_PSP", 2319), ("PD_vs_control", 2523),
                               ("PD_amppd_vs_control_amppd", 2508)]:
        record("§8", f"EUR {contrast} genome-wide hits", asserted, int(one("EUR", contrast, "n_gwsig")))
    record("§8", "EUR AD_ampad_vs_control_ampad genome-wide hits", 23,
           int(one("EUR", "AD_ampad_vs_control_ampad", "n_gwsig")))

    if "max_callset_delta" in g.columns and g.max_callset_delta.notna().any():
        record("§10", "callset-skew columns populated", True, True)
    else:
        unavailable("§10", "callset-skew columns populated", True,
                    "job 28004190 predates them; needs a step-7 rerun (HANDOFF issue 6)")


# ── cluster: artifacts that live only on biowulf ─────────────────────────────────────────

def check_cluster(root):
    merged = root / "data" / "merged"

    # §1 vs §3: the four normalized filesets sum to 13,428 and the merge reports 13,334. The
    # 94-genome difference is NOT QC — nothing is dropped upstream of the merge. It is 94 donors
    # carried in both wgs_harm and divco_hs under one sample ID, which plink fuses into a single
    # merged sample. 87 survive exclusion as `divco_hs|wgs_harm` rows in the grain
    # (PROJECT_LOG 2026-08-19, 2026-09-17). Three separate claims, so three separate records:
    # the inputs sum, the merge output, and the gap equals the collision count.
    stems = norm_stems(root)
    if stems is None:
        unavailable("§1/§3", "per-callset counts sum to 13,428", 13428,
                    f"could not source {root / 'config.sh'} for the NORM_* stems")
    else:
        total, breakdown = 0, []
        for key, stem in stems.items():
            n = wc(Path(stem + ".fam"))
            breakdown.append((key, Path(stem).name, n))
            if n is not None:
                total += n
        print("  per-callset .fam counts entering the merge:")
        for key, base, n in breakdown:
            print(f"    {key:9s} {base:50s} {n if n is not None else 'MISSING'}")
        if any(n is None for _, _, n in breakdown):
            unavailable("§1/§3", "per-callset counts sum to 13,428", 13428,
                        "a normalized .fam is missing — see the counts above")
        else:
            record("§1/§3", "per-callset counts sum to 13,428", 13428, total)
            n_merged = wc(merged / "cohort_merged.fam")
            if n_merged is None:
                unavailable("§3", "merged samples", 13334,
                            f"missing {merged / 'cohort_merged.fam'}")
            else:
                record("§3", "merged samples", 13334, n_merged)
                record("§3", "genomes fused at merge (13,428 − 13,334)", 94, total - n_merged)
            dups = id_collisions(stems)
            record("§3", "sample IDs shared across callsets", 94, dups)

    # §6.4's denominator. The unfiltered stage-A bim is EUR's post-QC variant count; the
    # difference against the stage-C bim is how many of the 4,187 were actually in EUR.
    unfiltered_bim = merged / "by_ancestry_qc" / "unfiltered" / "cohort_EUR_qc.bim"
    filtered_bim = merged / "by_ancestry_qc" / "cohort_EUR_qc.bim"
    n_unfiltered, n_filtered = wc(unfiltered_bim), wc(filtered_bim)
    if n_unfiltered is None or n_filtered is None:
        unavailable("§6.4", "EUR post-QC variant count (the 0.06% denominator)", 7538809,
                    "missing cohort_EUR_qc.bim under by_ancestry_qc[/unfiltered]")
    else:
        record("§6.4", "EUR post-QC variants (denominator)", 7538809, n_unfiltered)
        removed = n_unfiltered - n_filtered
        pct = 100 * removed / n_unfiltered
        print(f"  EUR: {n_unfiltered:,} → {n_filtered:,}, {removed:,} removed = {pct:.4f}%")
        record("§6.4", "excluded share of EUR variants (%, 2 dp)", 0.06, round(pct, 2), tol=0.005)

    excl = merged / "exclude_af_concordance.txt"
    record("§6.4", "exclusion list size", 4187, wc(excl)) if excl.exists() else unavailable(
        "§6.4", "exclusion list size", 4187, f"missing {excl}")

    # §6.3's gate table omits EUR divco_hs. The HWE channel needs MIN_CELL controls per
    # (stratum x callset); the grain is the only place the control count per cell exists.
    grain = root / "clinical_core_out" / "analysis_grain.csv"
    gr = read_csv(grain)
    if gr is None:
        unavailable("§6.3", "EUR divco_hs control count (missing gate-table row)", "not stated",
                    f"missing {grain}")
    elif len(gr) != 12495:
        unavailable("§6.3", "EUR divco_hs control count (missing gate-table row)", "not stated",
                    f"{grain} has {len(gr)} rows, not 12,495 — this is the stale laptop grain")
    else:
        pheno_col = "pheno" if "pheno" in gr.columns else None
        cs_col = "source_callset" if "source_callset" in gr.columns else None
        if not (pheno_col and cs_col):
            unavailable("§6.3", "EUR divco_hs control count", "not stated",
                        f"grain lacks pheno/source_callset (has {list(gr.columns)})")
        else:
            eur_ctrl = gr[(gr.ancestry == "EUR") & (gr[pheno_col] == "control")]
            print("  EUR control counts per callset (§6.3 gate table, MIN_CELL = 100):")
            for cs in ("wgs_harm", "wb_dwgs", "divco_hs", "br_dsnwgs"):
                n = sum(1 for s in eur_ctrl[cs_col] if cs in str(s).split("|"))
                print(f"    {cs:12s} {n:6d}{'   below MIN_CELL' if n < 100 else ''}")
    # The gate table's Controls column is NOT a grain count, and sourcing it from the grain
    # reported a false MISMATCH on 2026-09-17 (asserted 328, derived 373). Two different
    # populations: 6a ran inside step 6, before analysis_grain.py §13 became the sole
    # definition of who is a case, and its cells are the .keep files it handed plink. The
    # .keep is the artifact of record (CLAUDE.md rule 6) — 328, and the cell's own plink log
    # says "--keep: 328 samples remaining". The grain counts above are printed as context and
    # deliberately not compared: membership-splitting on source_callset gives 373 for EUR
    # wgs_harm, of which 33 are fused divco_hs|wgs_harm rows and 340 sole-source. Neither is
    # 328, and neither is meant to be.
    afc = merged / "af_concordance"
    for callset, asserted in (("wgs_harm", 328), ("wb_dwgs", 3064)):
        keep = afc / f"EUR_control_{callset}.keep"
        n = wc(keep)
        if n is None:
            unavailable("§6.3", f"EUR {callset} controls (HWE cell)", asserted,
                        f"missing {keep.relative_to(root)}")
        else:
            record("§6.3", f"EUR {callset} controls (HWE cell)", asserted, n)

    # §6.4's ~7x duplicate-discordance enrichment. It is a computation, not a lookup —
    # `--discordance` runs it. See measure_discordance().
    unavailable("§6.4", "duplicate-pair discordance enrichment", "~7x",
                "a computation, not an artifact — --discordance on biowulf (HANDOFF issue 9)")

    # §6.4 traceability: which step-6 job produced the live 4,187 list.
    unavailable("§6.4", "job ID behind the live 4,187 list", "not recorded",
                "recover with: bash scripts/runlog.sh --md")


# ── doc: claims METHODS.md makes about itself ────────────────────────────────────────────

def check_doc(root):
    path = root / "METHODS.md"
    if not path.exists():
        unavailable("doc", "METHODS.md self-consistency", "n/a", f"missing {path}")
        return
    text = path.read_text()

    # §1's table of supplied genomes must sum to the number §3 says it merged, minus the
    # per-callset QC loss §3 states. The table summing to something the text contradicts is
    # exactly the discrepancy this file was written to stop.
    rows = [ln for ln in text.split("\n") if ln.startswith("| `") and "AMP-" in ln]
    supplied = 0
    for ln in rows:
        for cell in ln.split("|"):
            cell = cell.strip().replace(",", "")
            if cell.isdigit():
                supplied += int(cell)
                break
    record("§1", "callset table sums to genomes as supplied", 13428, supplied)
    # Not a QC loss — nothing is dropped upstream of the merge. These are cross-callset ID
    # collisions that plink fuses into one sample each; check_cluster() derives the same 94
    # independently from the .fam files.
    record("§1/§3", "supplied − merged = genomes fused at merge", 94, supplied - 13334)

    limitations = text.split("## 10. Limitations")[-1]
    record("§10", "bolded limitations", 7,
           sum(1 for ln in limitations.split("\n") if ln.startswith("**")))


GROUPS = {"doc": check_doc, "manifest": check_manifest, "eta2": check_eta2,
          "gwas": check_gwas, "cluster": check_cluster}


# ── --discordance: the one claim in METHODS.md that is a computation ─────────────────────

def measure_discordance(root, plink2, tmp, dup_cutoff=0.354):
    """Genotype discordance between duplicate pairs, at flagged variants vs everywhere else.

    §6.4 claims the 4,187 flagged variants are "~7x enriched for genotype discordance between
    duplicate sample pairs". That is the only NON-CIRCULAR evidence the list is technical
    rather than real — everything else in §6 is measured on the same frequency signal the
    filter was built from. It has no numerator, denominator or job ID anywhere in the repo;
    HANDOFF.md issue 9 says measure it here or delete the sentence.

    A duplicate pair is the same person sequenced twice, so the true genotype is identical by
    construction and every difference is technical. If the flagged variants are technical
    breakage, they discord at a higher rate than the genome. If they are real biology, they
    discord at the background rate and §6's licence to DELETE is not supported.

    Runs on biowulf only: needs plink2 and the EUR genotypes. Cluster-only, read-only.
    """
    import subprocess

    merged = root / "data" / "merged"
    bfile = merged / "by_ancestry_qc" / "unfiltered" / "cohort_EUR_qc"
    excl = merged / "exclude_af_concordance.txt"
    pairdir = merged / "relatedness" / "by_ancestry"

    for p in (Path(str(bfile) + ".bed"), excl, pairdir):
        if not p.exists():
            print(f"CANNOT RUN: missing {p}")
            print("This is a cluster measurement. Nothing was estimated.")
            return 1

    frames = []
    for f in sorted(pairdir.glob("relatedness_*_pairs.related")):
        df = pd.read_csv(f)
        if df.empty:
            continue
        df.columns = [c.lstrip("#") for c in df.columns]
        frames.append(df)
    if not frames:
        print(f"CANNOT RUN: no relatedness_*_pairs.related under {pairdir}")
        return 1

    pairs = pd.concat(frames, ignore_index=True)
    pairs["KINSHIP"] = pd.to_numeric(pairs["KINSHIP"], errors="coerce")
    dup = pairs[pairs.KINSHIP >= dup_cutoff]

    # Restrict to pairs whose BOTH members are in this fileset, or plink2 errors on the pair.
    fam = pd.read_csv(str(bfile) + ".fam", sep=r"\s+", header=None, dtype=str)
    present = set(fam[1])
    dup = dup[dup.IID1.astype(str).isin(present) & dup.IID2.astype(str).isin(present)]
    if dup.empty:
        print("CANNOT RUN: no duplicate pair has both members in the EUR fileset.")
        return 1

    tmp.mkdir(parents=True, exist_ok=True)
    pairfile = tmp / "dup_pairs.txt"
    pairfile.write_text("".join(f"{a}\t{b}\n" for a, b in zip(dup.IID1, dup.IID2)))
    print(f"{len(dup)} duplicate pairs at kinship >= {dup_cutoff}, both members in EUR.")

    def run(tag, variant_arg):
        out = tmp / tag
        cmd = [plink2, "--bfile", str(bfile), "--sample-diff", "counts-only",
               "id-pairs", f"file={pairfile}", *variant_arg, "--out", str(out)]
        print("  " + " ".join(cmd))
        subprocess.run(cmd, check=True)
        d = pd.read_csv(str(out) + ".sdiff.summary", sep=r"\s+")
        d.columns = [c.lstrip("#") for c in d.columns]
        diff, obs = int(d.DIFF_CT.sum()), int(d.OBS_CT.sum())
        return diff, obs, diff / obs if obs else float("nan")

    flag_diff, flag_obs, flag_rate = run("flagged", ["--extract", str(excl)])
    back_diff, back_obs, back_rate = run("background", ["--exclude", str(excl)])

    print(f"\n  flagged     {flag_diff:>12,} / {flag_obs:>14,} = {flag_rate:.6f}")
    print(f"  background  {back_diff:>12,} / {back_obs:>14,} = {back_rate:.6f}")
    print(f"  enrichment  {flag_rate / back_rate:.2f}x   over {len(dup)} duplicate pairs")
    print("\nRecord both rates, the pair count and the job ID in PROJECT_LOG.md, then cite "
          "them from METHODS.md §6.4 in place of the unsourced multiplier.")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--only", choices=sorted(GROUPS), nargs="+")
    ap.add_argument("--strict", action="store_true",
                    help="exit 1 on any MISMATCH or UNAVAILABLE")
    ap.add_argument("--json", action="store_true", help="emit results as JSON")
    ap.add_argument("--discordance", action="store_true",
                    help="measure §6.4's duplicate-pair discordance enrichment (cluster only)")
    ap.add_argument("--plink2", default="plink2",
                    help="plink2 path; use an absolute one if plink/1.9 is loaded")
    ap.add_argument("--tmp", type=Path, default=Path("sdiff_tmp"))
    args = ap.parse_args()

    if args.discordance:
        return measure_discordance(args.root, args.plink2, args.tmp)

    for name in (args.only or list(GROUPS)):
        print(f"\n── {name} " + "─" * (72 - len(name)))
        GROUPS[name](args.root)

    if args.json:
        print(json.dumps([{"section": s, "claim": c, "asserted": str(a),
                           "derived": str(d), "status": st}
                          for s, c, a, d, st in results], indent=2))
        return 0

    width = max(len(c) for _, c, _, _, _ in results) + 2
    print(f"\n{'':10s}{'claim':{width}s}{'asserted':>14s}   derived")
    print("─" * (26 + width + 30))
    for section, claim, asserted, derived, status in results:
        mark = {OK: "  ok ", MISMATCH: " FAIL", UNAVAILABLE: " ????"}[status]
        print(f"{mark:6s}{section:8s}{claim:{width}s}{str(asserted):>14s}   {derived}")

    bad = [r for r in results if r[4] == MISMATCH]
    unk = [r for r in results if r[4] == UNAVAILABLE]
    print(f"\n{len(results) - len(bad) - len(unk)} ok · {len(bad)} MISMATCH · {len(unk)} could not run")
    if unk:
        print("A check that could not run is NOT a pass. Re-run on biowulf to settle these.")
    return 1 if (args.strict and (bad or unk)) else (1 if bad else 0)


if __name__ == "__main__":
    sys.exit(main())
