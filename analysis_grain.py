#!/usr/bin/env python
# coding: utf-8

# # Analysis grain — the second half of the clinical side
#
# Runs ONCE, **after genetics step 6**. `clinical_core.py` runs once *before* step 1; this is the
# other end of that round trip, and the round trip is irreducible: step 1 applies the sex files,
# and the grain carries step 6's PCs because step 7 reads them as covariates.
#
# | § | Does | Needs |
# |---|---|---|
# | 11 | **in <-** genotools/relatedness QC outcomes, labelled by reason | steps 1-5 |
# | 12 | **in <-** ancestry + PCs, reconciled to one label per genome | step 6 |
# | 13 | **out ->** per-ancestry covariate and per-contrast phenotype files | step 6 |
#
# It reads §9's audit tables and re-derives nothing. `core` and `genomes` below are exactly the
# frames `clinical_core.py` held in memory — which is why this is a separate file rather than
# sections at the end of that one, where reaching them meant rewriting the sex-update files step 1
# had already consumed.
#
# **Guardrail.** Every cell prints aggregates. No cell prints a subject-level row.


import pandas as pd

from clinical_common import (
    AMPPD_CALLSETS, EXCLUDE_REASONS, OUT, PCA_DIR, RETAINED_MANIFEST,
    load_audit_tables, rd, reconcile, rel,
)

pd.set_option("display.width", 130)

core, genomes = load_audit_tables()
print(f"individual_core.csv   {len(core):,} rows")
print(f"genome_crosswalk.csv  {len(genomes):,} rows")



# ## 11. In ← QC outcomes
#
# Attaches steps 1-5's keep/drop reasons to the clinical table, so a loss reads per cohort and per
# phenotype arm rather than per callset. Consumes `05_excludelist.py`'s output; recomputes nothing.
#
# | File | Columns |
# |---|---|
# | `exclude_reasons.tsv` | `FID  IID  reason  detail` |
# | `retained_manifest.csv` | `FID,IID,ancestry,call_rate,dup_cluster_id` |

if not EXCLUDE_REASONS.exists():
    print(f"SKIPPED — no exclude reasons at {rel(EXCLUDE_REASONS)}")
    print("Run genetics steps 1-5, then rsync data/merged/relatedness/ down.")
else:
    reasons = rd(EXCLUDE_REASONS, sep="\t")
    retained = rd(RETAINED_MANIFEST)

    outcome = genomes[["IID", "source_callset", "individual_id", "source_dataset"]].copy()
    outcome = outcome.merge(reasons[["IID", "reason", "detail"]], on="IID", how="left")
    outcome = outcome.merge(retained[["IID", "ancestry", "call_rate"]], on="IID", how="left")
    outcome["reason"] = outcome["reason"].fillna("")
    outcome["status"] = (outcome.IID.isin(set(retained.IID))
                         .map({True: "retained", False: "dropped"}))

    labs = core.set_index(["individual_id", "source_dataset"])[["pheno", "dx_detailed"]]
    outcome = outcome.join(labs, on=["individual_id", "source_dataset"])
    outcome[["pheno", "dx_detailed"]] = outcome[["pheno", "dx_detailed"]].fillna("")

    outcome.to_csv(OUT / "qc_outcomes.csv", index=False)

    print("status x callset:")
    print(pd.crosstab(outcome.source_callset, outcome.status, margins=True).to_string())
    print("\ndrop reason x callset:")
    print(pd.crosstab(outcome.loc[outcome.status == "dropped", "reason"],
                      outcome.loc[outcome.status == "dropped", "source_callset"],
                      margins=True).to_string())
    print("\nretained genomes by phenotype arm:")
    print(pd.crosstab(outcome.loc[outcome.status == "retained", "pheno"].replace("", "(null)"),
                      outcome.loc[outcome.status == "retained", "source_callset"],
                      margins=True).to_string())
    print(f"\nqc_outcomes.csv  {len(outcome):,} rows -> {rel(OUT)}")



# ## 12. In ← ancestry and PCs → the analysis grain
#
# One row per retained genome: the label the GWAS tests and the covariates it adjusts for. Joins the
# retained manifest, step 6's per-ancestry PCs, and §7/§6's audit tables. Nothing is re-derived.
#
# A donor can appear in two source studies — Diverse Cohorts and the ROSMAP/Mayo/MSBB trio —
# phenotyped by programs that do not always agree, so each genome needs one reconciled label:
#
# | Field | Rule |
# |---|---|
# | `pheno` | agree → that value; exactly one labelled → the labelled one; both labelled and differing → **AD-dominant** (either says AD → AD), else the non-control label |
# | `dx_detailed` | pinned to `AD` when reconciled `pheno` is AD; otherwise prefer the trio's clinical instrument, falling back to Diverse Cohorts when the trio is null |
# | `sex` | each genome takes **its own callset's** donor's sex — never reconciled |
#
# AD dominates because it is a neuropathological call and outranks a clinical one. `sex` is NOT
# reconciled: it is a property of the sequenced sample, and disagreement is evidence of a swap.
# Conflicts are flagged, never dropped — `pheno_conflict`/`dx_conflict`/`sex_conflict` travel with
# the row so a sensitivity analysis can exclude them without rebuilding.

eigenvecs = sorted(PCA_DIR.glob("cohort_*_pca.eigenvec")) if PCA_DIR.exists() else []

if not RETAINED_MANIFEST.exists() or not eigenvecs:
    print(f"SKIPPED — need {rel(RETAINED_MANIFEST)} and {rel(PCA_DIR)}/cohort_*_pca.eigenvec")
    print("Run genetics steps 1-6, then rsync data/merged/ down.")
else:
    pcs = pd.concat([rd(f, sep=r"\s+") for f in eigenvecs], ignore_index=True)
    pcs = pcs.rename(columns={"#FID": "FID"})
    pc_cols = [c for c in pcs.columns if c.upper().startswith("PC")][:10]

    manifest = rd(RETAINED_MANIFEST).merge(pcs[["IID"] + pc_cols], on="IID", how="left")

    # donor -> their core rows (a donor in two cohorts has two)
    by_donor = {}
    for r in core.to_dict("records"):
        by_donor.setdefault(r["individual_id"], []).append(r)

    # genome -> every (donor, cohort, callset) it resolved to.  A list, not a single row:
    # the same IID exists in both divco_hs and wgs_harm for the fused dual-source samples.
    by_genome = {}
    for r in genomes.itertuples():
        by_genome.setdefault(r.IID, []).append(r)

    grain_rows = []
    for m in manifest.to_dict("records"):
        hits = by_genome.get(m["IID"], [])
        donor = next((h.individual_id for h in hits if h.individual_id), "")
        rows = by_donor.get(donor, []) if donor else []

        pheno, dx, pconf, dconf = reconcile(rows) if rows else ("", "", False, False)

        # sex comes from the cohort(s) THIS genome resolved to, not the reconciliation
        own_cohorts = {h.source_dataset for h in hits}
        own_sex = {r["sex"] for r in rows if r["source_dataset"] in own_cohorts
                   and r["sex"] in ("1", "2")}
        sex = sorted(own_sex)[0] if own_sex else ""
        sconf = len(own_sex) > 1

        grain_rows.append({
            "IID": m["IID"], "individual_id": donor,
            "source_callset": "|".join(sorted({h.source_callset for h in hits}))
                              or m.get("source_callset", ""),
            "ancestry": m["ancestry"], "sex": sex, "pheno": pheno, "dx_detailed": dx,
            "pheno_conflict": int(pconf), "dx_conflict": int(dconf), "sex_conflict": int(sconf),
            "call_rate": m.get("call_rate", ""), "dup_cluster_id": m.get("dup_cluster_id", ""),
            **{c: m.get(c, "") for c in pc_cols}})

    GRAIN_COLUMNS = ["IID", "individual_id", "source_callset", "ancestry", "sex", "pheno",
                     "dx_detailed", "pheno_conflict", "dx_conflict", "sex_conflict",
                     "call_rate", "dup_cluster_id"] + pc_cols
    grain = pd.DataFrame(grain_rows)[GRAIN_COLUMNS]
    grain.to_csv(OUT / "analysis_grain.csv", index=False)

    print(f"retained genomes: {len(grain):,}")
    print(f"  no donor resolved:  {int((grain.individual_id == '').sum()):,}")
    print(f"  pheno conflicts:    {int(grain.pheno_conflict.sum()):,}")
    print(f"  dx conflicts:       {int(grain.dx_conflict.sum()):,}")
    print(f"  sex conflicts:      {int(grain.sex_conflict.sum()):,}")
    print(f"  missing PCs:        {int((grain[pc_cols[0]] == '').sum()):,} (sub-50-sample strata)")
    print("\ndx_detailed x ancestry:")
    print(pd.crosstab(grain.dx_detailed.replace("", "(null)"), grain.ancestry,
                      margins=True).to_string())
    print(f"\nanalysis_grain.csv  {len(grain):,} rows x {len(GRAIN_COLUMNS)} cols -> {rel(OUT)}")


# ## 13. Out → phenotype, covariate, and contrast-manifest files
#
# Everything `07_gwas.sh` needs. **This section is the sole definition of who is a case.**
# Step 7 reads these files and runs `plink2 --glm`; it does not parse the grain and does not
# rebuild an arm in awk, so the definition exists once and cannot drift between the two.
#
# | File | Format |
# |---|---|
# | `covar_<ANC>.txt` | `#FID IID SEX [AGE] PC1..PC10`, missing as `NA` |
# | `pheno_<ANC>_<CASE>_vs_<CTRL>.txt` | `#FID IID pheno`, case=2 control=1 |
# | `contrasts.csv` | one row per written contrast: arm sizes, cohort composition, confound tag, callset skew |
#
# **FID comes from the genotype fileset**, not the grain: plink2 matches on FID+IID and
# defaults a missing FID to `0`, so the AMP-PD callsets — which carry a non-zero FID — need
# their real one. `cohort_<ANC>_qc.fam` therefore has to exist alongside the PCs.
#
# A contrast arm may be restricted by source: `control@amppd` means AMP-PD callsets only
# (`wb_dwgs` or `br_dsnwgs`), `@ampad` the rest. Bare arms take any source. Both arms must
# reach `MIN_ARM` for a file to be written.
#
# **`delta_amppd`** is each contrast's gap in AMP-PD share — how much of a result could be cohort
# rather than disease. Computed here off `AMPPD_CALLSETS`, so adding a callset updates it in one
# place. It measures PROGRAM, not callset; see `max_callset_delta` below and METHODS.md §10.
#
# **`EXCLUDE_DUAL`** drops samples resolving to both `divco_hs` and `wgs_harm`, for a WGS_Harm-only
# sensitivity run. Here rather than a flag at GWAS time, so the variant is a recorded artifact.

CONTRASTS = [
    ("PD", "AD"), ("PD", "DLB"), ("PD", "MCI"), ("PD", "PSP"), ("PD", "control"), ("PD", "other"),
    ("AD", "MCI"), ("AD", "DLB"), ("AD", "PSP"), ("AD", "control"), ("AD", "other"),
    ("control@amppd", "control@ampad"), ("PD@amppd", "control@amppd"), ("AD@ampad", "control@ampad"),
]
MIN_ARM = 20
FAM_DIR = PCA_DIR          # cohort_<ANC>_qc.fam lives beside the eigenvecs

DUAL_CALLSET = "divco_hs|wgs_harm"   # §12 joins multi-callset genomes with "|"
EXCLUDE_DUAL = False                 # True -> WGS_Harm-only sensitivity variant

# Age is the dominant confounder for both AD and PD, and the grain does not carry it yet:
# AMP-AD records age-at-death, AMP-PD age-at-baseline, and they are not the same variable.
# Emitted into covar the moment a harmonized column appears, so step 7 needs no change.
AGE_COLUMNS = ("age", "age_analysis", "agedeath", "age_death", "age_baseline", "age_cov")


def arm_mask(g, arm):
    """'AD@ampad' -> rows with dx_detailed == AD from a non-AMP-PD callset."""
    dx, _, src = arm.partition("@")
    m = g.dx_detailed == dx
    if src == "amppd":
        return m & g.source_callset.isin(AMPPD_CALLSETS)
    if src == "ampad":
        return m & ~g.source_callset.isin(AMPPD_CALLSETS)
    return m


def amppd_pct(g, mask):
    """Share of an arm that came from an AMP-PD callset. Drives the confound tag."""
    n = int(mask.sum())
    if not n:
        return float("nan")
    return 100.0 * int((mask & g.source_callset.isin(AMPPD_CALLSETS)).sum()) / n


ONE_SIDED_MIN = 10   # samples in one arm, zero in the other, before it is worth naming


def callset_skew(g, cm, km):
    """-> (max |case% - ctrl%| over single callsets, which callset, one-sided callsets).

    WHY THIS EXISTS SEPARATELY FROM delta_amppd. `delta_amppd` pools wb_dwgs with br_dsnwgs
    because both are AMP-PD, which is the right grain for the disease/program confound — but it
    is blind to an asymmetry *inside* a program, and on 2026-08-21 that turned out to matter a
    lot. EUR PD-vs-DLB scored delta_amppd = 0.0 and was tagged `within_cohort`, the label
    07_gwas.sh's own header calls "the confound-free trusted backbone" — while carrying
    BR-DSNWGS on the PD arm and none on the DLB arm (BR's 95 retained samples are 71 PD /
    21 control / 3 null, so it contributes no DLB at all). BR sits at ~50% missingness on the
    common set, and step 7's differential-missingness filter excluded **353,068** variants from
    that contrast — ~6x any genuinely cross-program contrast (~60k). Seventy-one samples
    produced more technical asymmetry than the entire AMP-AD/AMP-PD split.

    WHY ONE-SIDEDNESS AND NOT A PERCENTAGE GAP. The first version of this returned only the
    max percentage-point delta and reused the confound tag's 20pp threshold, which was tidy and
    WRONG: BR is ~55 of the 2,595-sample PD arm, so its delta is **2.1pp** and a 20pp rule misses
    the one case the function was written to catch. The quantity that matters is not how far the
    shares differ but whether a callset is ABSENT from one arm — every site that callset uniquely
    fails to call is then differentially missing by construction, however small its share. So the
    percentage figures are kept as description and the flag is categorical.

    ONE_SIDED_MIN = 10 is a judgment call, not a derived number: enough samples for plink's 2x2
    missingness test to resolve at these arm sizes, low enough to catch a callset the size of BR.

    Reported per callset rather than as a BR-specific flag on purpose: BR is the instance, an
    unbalanced callset is the class, and a fifth callset should be caught without an edit here.
    Deliberately NOT folded into confound_tag — that tag means "program", and changing what it
    measures would reinterpret every row already on record. A contrast can legitimately be
    `within_cohort` AND callset-skewed; step 7 flags that pair rather than reclassifying it.
    """
    n_c, n_k = int(cm.sum()), int(km.sum())
    if not n_c or not n_k:
        return float("nan"), "NA", "NA"
    worst, worst_cs, one_sided = -1.0, "NA", []
    for cs in sorted(set(g.source_callset.dropna())):
        in_cs = g.source_callset == cs
        a, b = int((cm & in_cs).sum()), int((km & in_cs).sum())
        d = abs(100.0 * a / n_c - 100.0 * b / n_k)
        if d > worst:                      # sorted() iteration -> ties resolve alphabetically
            worst, worst_cs = d, cs
        if (a >= ONE_SIDED_MIN and b == 0) or (b >= ONE_SIDED_MIN and a == 0):
            one_sided.append(f"{cs}({a}v{b})")
    return worst, worst_cs, (";".join(one_sided) if one_sided else "none")


def confound_tag(delta):
    """<=20pp apart -> the arms share a cohort base; >=70pp -> they are effectively disjoint."""
    if pd.isna(delta):
        return "NA"
    return "within_cohort" if delta <= 20 else ("cross_cohort" if delta >= 70 else "partial")


if "grain" not in globals():
    print("SKIPPED — §12 has not produced a grain yet.")
else:
    PHENO_DIR = OUT / "pheno"
    COVAR_DIR = OUT / "covar"
    PHENO_DIR.mkdir(parents=True, exist_ok=True)
    COVAR_DIR.mkdir(parents=True, exist_ok=True)

    age_col = next((c for c in grain.columns if c.strip().lower() in AGE_COLUMNS), None)
    print(f"age covariate: {age_col or 'NOT FOUND in grain — running WITHOUT age'}")
    if EXCLUDE_DUAL:
        print(f"EXCLUDE_DUAL on — dropping source_callset == {DUAL_CALLSET!r}")

    written = []
    for anc, g in grain.groupby("ancestry"):
        fam = FAM_DIR / f"cohort_{anc}_qc.fam"
        if not fam.exists():
            print(f"{anc:5} SKIPPED — no {fam.name}")
            continue
        if EXCLUDE_DUAL:
            g = g[g.source_callset != DUAL_CALLSET]
        # .fam is headerless: FID IID PAT MAT SEX PHENO
        f = pd.read_csv(fam, sep=r"\s+", dtype=str, header=None, usecols=[0, 1])
        fid = dict(zip(f[1], f[0]))
        g = g.assign(FID=g.IID.map(lambda i: fid.get(i, "0")))

        cov_cols = ["FID", "IID", "sex"] + ([age_col] if age_col else []) + pc_cols
        cov_names = ["#FID", "IID", "SEX"] + (["AGE"] if age_col else []) + pc_cols
        covar = g[cov_cols].replace("", "NA")
        covar.columns = cov_names
        covar.to_csv(COVAR_DIR / f"covar_{anc}.txt", sep="\t", index=False)

        for case, ctrl in CONTRASTS:
            cm, km = arm_mask(g, case), arm_mask(g, ctrl)
            if cm.sum() < MIN_ARM or km.sum() < MIN_ARM:
                continue
            tag = f"{case}_vs_{ctrl}".replace("@", "_")
            out = pd.concat([
                g.loc[cm, ["FID", "IID"]].assign(pheno="2"),
                g.loc[km, ["FID", "IID"]].assign(pheno="1")])
            out.columns = ["#FID", "IID", "pheno"]
            out.to_csv(PHENO_DIR / f"pheno_{anc}_{tag}.txt", sep="\t", index=False)

            n_case, n_ctrl = int(cm.sum()), int(km.sum())
            case_pct, ctrl_pct = amppd_pct(g, cm), amppd_pct(g, km)
            delta = abs(case_pct - ctrl_pct)
            cs_delta, cs_worst, cs_one_sided = callset_skew(g, cm, km)
            written.append({
                "ancestry": anc, "contrast": tag, "case_arm": case, "ctrl_arm": ctrl,
                "n_case": n_case, "n_ctrl": n_ctrl,
                "case_pct_amppd": round(case_pct, 1), "ctrl_pct_amppd": round(ctrl_pct, 1),
                "delta_amppd": round(delta, 1), "confound_tag": confound_tag(delta),
                # Per-callset asymmetry, which delta_amppd cannot see. The first two describe;
                # callset_one_sided is the one that predicts differential missingness, because a
                # callset absent from an arm makes every site it uniquely fails to call
                # differentially missing regardless of how small its share is.
                "max_callset_delta": round(cs_delta, 1), "worst_callset": cs_worst,
                "callset_one_sided": cs_one_sided,
                "viable_ge100": int(n_case >= 100 and n_ctrl >= 100)})

    if written:
        w = pd.DataFrame(written).sort_values(["ancestry", "contrast"], ignore_index=True)
        # The manifest step 7 loops over. Carrying the arm sizes and cohort composition here
        # means step 7 never has to re-derive an arm to report on it.
        w.to_csv(OUT / "contrasts.csv", index=False)
        print(w.to_string(index=False))
        print(f"\n{len(w)} phenotype files + {w.ancestry.nunique()} covariate files "
              f"+ contrasts.csv -> {rel(OUT)}")
        print(f"viable at >=100 per arm: {int(w.viable_ge100.sum())}")
        print("confound tag: " + "  ".join(f"{k}={v}" for k, v in
                                           w.confound_tag.value_counts().items()))
    else:
        print("no contrast reached MIN_ARM in any ancestry")
        pd.DataFrame(columns=[
            "ancestry", "contrast", "case_arm", "ctrl_arm", "n_case", "n_ctrl",
            "case_pct_amppd", "ctrl_pct_amppd", "delta_amppd", "confound_tag",
            "viable_ge100"]).to_csv(OUT / "contrasts.csv", index=False)

