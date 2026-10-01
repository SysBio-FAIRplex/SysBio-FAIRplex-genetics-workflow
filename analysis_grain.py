#!/usr/bin/env python
# coding: utf-8

# # Analysis grain — the second half of the clinical side
#
# Runs ONCE, **after genetics step 6**; `clinical_core.py` is the first half.
#
# | § | Does | Needs |
# |---|---|---|
# | 11 | **in <-** genotools/relatedness QC outcomes, labelled by reason | steps 1-5 |
# | 12 | **in <-** ancestry + PCs, reconciled to one label per genome | step 6 |
# | 13 | **out ->** per-ancestry covariate and per-contrast phenotype files | step 6 |
#
# Reads §9's audit tables and re-derives nothing.
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
# Steps 1-5's keep/drop reasons per cohort and phenotype arm, from `05_excludelist.py`'s
# `exclude_reasons.tsv` (FID IID reason detail) and `retained_manifest.csv`.

if not EXCLUDE_REASONS.exists():
    print(f"SKIPPED — no exclude reasons at {rel(EXCLUDE_REASONS)}")
    print("Run genetics steps 1-5 first; this reads data/merged/relatedness/ in place.")
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
# One row per retained genome: the label the GWAS tests and its covariates. A donor in two source
# studies gets one reconciled label (clinical_common.reconcile):
#
# | Field | Rule |
# |---|---|
# | `pheno` | agree → that value; exactly one labelled → the labelled one; both labelled and differing → **AD-dominant** (either says AD → AD), else the non-control label |
# | `dx_detailed` | pinned to `AD` when reconciled `pheno` is AD; otherwise prefer the trio's clinical instrument, falling back to Diverse Cohorts when the trio is null |
# | `sex` | each genome takes **its own callset's** donor's sex — never reconciled |
#
# Conflicts are flagged, never dropped, so a sensitivity analysis can exclude them.

eigenvecs = sorted(PCA_DIR.glob("cohort_*_pca.eigenvec")) if PCA_DIR.exists() else []

if not RETAINED_MANIFEST.exists() or not eigenvecs:
    print(f"SKIPPED — need {rel(RETAINED_MANIFEST)} and {rel(PCA_DIR)}/cohort_*_pca.eigenvec")
    print("Run genetics steps 1-6 first; this reads data/merged/ in place.")
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
# Everything `07_gwas.sh` needs. **This section is the sole definition of who is a case**; step 7
# reads these files and never rebuilds an arm.
#
# | File | Format |
# |---|---|
# | `covar_<ANC>.txt` | `#FID IID SEX [AGE] PC1..PC10`, missing as `NA` |
# | `pheno_<ANC>_<CASE>_vs_<CTRL>.txt` | `#FID IID pheno`, case=2 control=1 |
# | `contrasts.csv` | one row per written contrast: arm sizes, cohort composition, confound tag, callset skew |
#
# **FID comes from `cohort_<ANC>_qc.fam`**, not the grain: plink2 matches on FID+IID, and the
# AMP-PD callsets carry a non-zero FID.
#
# `X@amppd` restricts an arm to AMP-PD callsets, `X@ampad` to the rest; both arms must reach
# `MIN_ARM`. `delta_amppd` is the gap in AMP-PD share — PROGRAM, not callset (METHODS.md §10).
# `EXCLUDE_DUAL` writes a WGS_Harm-only sensitivity variant as a recorded artifact.

CONTRASTS = [
    ("PD", "AD"), ("PD", "DLB"), ("PD", "MCI"), ("PD", "PSP"), ("PD", "control"), ("PD", "other"),
    ("AD", "MCI"), ("AD", "DLB"), ("AD", "PSP"), ("AD", "control"), ("AD", "other"),
    ("control@amppd", "control@ampad"), ("PD@amppd", "control@amppd"), ("AD@ampad", "control@ampad"),
]
MIN_ARM = 20
FAM_DIR = PCA_DIR          # cohort_<ANC>_qc.fam lives beside the eigenvecs

DUAL_CALLSET = "divco_hs|wgs_harm"   # §12 joins multi-callset genomes with "|"
EXCLUDE_DUAL = False                 # True -> WGS_Harm-only sensitivity variant

# No harmonized age yet (AMP-AD: at death; AMP-PD: at baseline — METHODS.md §10). Emitted into
# covar as soon as one of these columns appears; step 7 needs no change.
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

    delta_amppd pools AMP-PD callsets, so it cannot see a callset present in one arm and absent
    from the other — and absence, not the size of the share, is what makes every site that
    callset fails to call differentially missing (METHODS.md §10). Hence a categorical one-sided
    flag rather than a percentage threshold. Kept separate from confound_tag, which means program.

    ONE_SIDED_MIN is a judgment call: enough samples for plink's missingness test to resolve,
    few enough to catch a callset the size of BR-DSNWGS.
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
                # Per-callset asymmetry; callset_one_sided is the one that predicts missingness.
                "max_callset_delta": round(cs_delta, 1), "worst_callset": cs_worst,
                "callset_one_sided": cs_one_sided,
                "viable_ge100": int(n_case >= 100 and n_ctrl >= 100)})

    if written:
        w = pd.DataFrame(written).sort_values(["ancestry", "contrast"], ignore_index=True)
        # The manifest step 7 loops over, with arm sizes and composition so it re-derives nothing.
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

