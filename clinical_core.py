#!/usr/bin/env python
# coding: utf-8

# # Harmonized clinical core
#
# Turns per-cohort clinical files into the phenotype inputs the GWAS needs. **Runs ONCE, before
# genetics step 1**, and needs no genotype-pipeline output. Criteria: METHODS.md §2.
#
# | § | Does | Needs |
# |---|---|---|
# | 1–6 | read clinical files, derive `pheno` / `dx_detailed` per donor | clinical files |
# | 7–9 | resolve genotype samples to donors, check, write audit tables | psam files |
# | 10 | **out →** per-callset sex-update files | — |
# | 12a | **out →** `sample_annot.csv`, the grain's PC-free half | — |
# | 14 | appendix: how definition-dependent is the AD arm (read-only) | — |
#
# §11–13 live in `analysis_grain.py`, which runs after step 6; section numbers are shared.
# Acquisition commands for every input: wgs_core.ipynb §0.
#
# **Guardrail.** Every cell prints aggregates. No cell prints a subject-level row.
# Genotype files are read for sample IDs only and are never modified.

import pandas as pd

# Shared with analysis_grain.py so both halves apply identical reconciliation (§12a checks it).
from clinical_common import (
    AD_META, AMPPD_META, DIR_BR, DIR_DC, DIR_WB, DIR_WGS, DIVCO_META,
    DX_VALUES, OUT, PHENO_VALUES, SEX_VALUES,
    rd, reconcile, rel,
)

pd.set_option("display.width", 130)


# ## 1. Clinical sources
#
# One file per cohort at individual grain; AMP-PD takes three, joined on `participant_id`.

clinical_files = {
    "amp_ad_rosmap": AD_META / "ROSMAP_clinical_harmonized.csv",
    "amp_ad_mayo":   AD_META / "MayoRNAseq_individual_metadata_harmonized.csv",
    "amp_ad_msbb":   AD_META / "MSBB_individual_metadata_harmonized.csv",
    "amp_ad_divco":  DIVCO_META / "AMP-AD_DiverseCohorts_individual_metadata.csv",
}

# AMP-PD spans three files. Covers every AMP-PD study, antemortem and postmortem alike.
amp_pd = (rd(AMPPD_META / "amp_pd_participants.csv")
          .merge(rd(AMPPD_META / "amp_pd_case_control.csv"), on="participant_id", how="left")
          .merge(rd(AMPPD_META / "Demographics.csv"), on="participant_id", how="left"))

# the donor-id column differs: AMP-PD calls it participant_id, everyone else individualID
ID_COL = {s: "individualID" for s in clinical_files}
ID_COL["amp_pd"] = "participant_id"

raw = {src: rd(path) for src, path in clinical_files.items()}
raw["amp_pd"] = amp_pd

for src, df in raw.items():
    ids = df[ID_COL[src]]
    print(f"{src:16} {len(df):7,} rows x {df.shape[1]:3} cols   "
          f"unique ids: {ids.nunique():,}")


# ## 2. The raw phenotype instruments
#
# The columns every derived label traces back to. The cohorts do not share an instrument
# (METHODS.md §2), and AMP-PD has no neuropathology at all.

instruments = {
    "amp_ad_rosmap": ["Braak", "amyCerad", "dcfdx_lv"],
    "amp_ad_mayo":   ["Braak", "amyThal", "diagnosis"],
    "amp_ad_msbb":   ["Braak", "amyCerad"],
    "amp_ad_divco":  ["dataContributionGroup", "ADoutcome", "mayoDx"],
    "amp_pd":        ["case_control_other_latest", "diagnosis_latest"],
}

for src, cols in instruments.items():
    print("=" * 70)
    print(src)
    for c in cols:
        vc = raw[src][c].replace("", "(blank)").value_counts()
        head = vc.head(12)
        print(f"\n  {c}  ({len(vc)} distinct)")
        for val, n in head.items():
            print(f"    {str(val)[:52]:54} {n:6,}")
        if len(vc) > len(head):
            print(f"    ... {len(vc) - len(head)} more")
    print()


# ## 3. Criteria — `pheno`
#
# `{AD, PD, control, other}` or null. **AD is neuropathological**; clinical codes never define it.
# Missing inputs yield null, not `other`, and nulls are excluded from every arm. Rationale:
# METHODS.md §2.
#
# | Cohort | Rule |
# |---|---|
# | ROSMAP, MSBB | AD = Braak≥IV & CERAD moderate/frequent; control = Braak≤III & CERAD sparse/none; else other |
# | Mayo | AD = Braak≥IV & Thal≥2; control = Braak≤III & Thal<2; else other — then controls whose `diagnosis` names a disease become other |
# | Diverse Cohorts | pre-adjudicated `ADoutcome`, or `mayoDx` for the Mayo contribution group |
# | AMP-PD | `case_control_other_latest`: Case→PD, Control→control, Other excluded |

MISSING = {"missing or unknown", "NA", ""}

BRAAK_HIGH = {"Stage IV", "Stage V", "Stage VI"}                    # >= IV
BRAAK_LOW  = {"None", "Stage I", "Stage II", "Stage III"}           # <= III

CERAD_HIGH = {"Frequent/Definite/C3", "Moderate/Probable/C2"}
CERAD_LOW  = {"Sparse/Possible/C1", "None/No AD/C0"}

THAL_HIGH  = {"Phase 2", "Phase 3", "Phase 4", "Phase 5"}           # >= 2
THAL_LOW   = {"None", "Phase 1"}                                    # < 2


def pheno_braak_cerad(r):
    """ROSMAP, MSBB."""
    braak, cerad = r["Braak"], r["amyCerad"]
    if braak in MISSING or cerad in MISSING:
        return ""
    if braak in BRAAK_HIGH and cerad in CERAD_HIGH:
        return "AD"
    if braak in BRAAK_LOW and cerad in CERAD_LOW:
        return "control"
    return "other"


def pheno_mayo(r):
    """Mayo: Braak + Thal, then demote impure controls."""
    braak, thal, dx = r["Braak"], r["amyThal"], r["diagnosis"]
    if braak in MISSING or thal in MISSING:
        return ""
    if braak in BRAAK_HIGH and thal in THAL_HIGH:
        return "AD"
    if braak in BRAAK_LOW and thal in THAL_LOW:
        return "control" if dx in {"control", "NA", ""} else "other"
    return "other"


def pheno_divco(r):
    """Diverse Cohorts: pre-adjudicated outcome; Mayo group uses its own column."""
    field = "mayoDx" if r["dataContributionGroup"] == "Mayo" else "ADoutcome"
    return {"AD": "AD", "Control": "control", "Other": "other"}.get(r[field], "")


def pheno_amppd(r):
    """AMP-PD: curated case/control/other flag. "Other" (LBD, MSA, PSP, essential tremor)
    is excluded from pheno."""
    return {"Case": "PD", "Control": "control"}.get(r["case_control_other_latest"], "")


print("pheno rules defined")


# ## 4. Criteria — `dx_detailed`
#
# `{PD, AD, MCI, DLB, PSP, control, other}` or null, for secondary analyses. **If `pheno` is AD,
# `dx_detailed` is AD** in every cohort. Otherwise (METHODS.md §2):
#
# | Cohort | Instrument |
# |---|---|
# | AMP-PD | `diagnosis_latest` → PD / DLB / PSP / control; everything else other |
# | ROSMAP | `dcfdx_lv` → 1 control, 2–3 MCI, 4–6 other (4/5 are clinical AD-dementia, not AD) |
# | Mayo | `diagnosis` == progressive supranuclear palsy → PSP; else keep `pheno` |
# | MSBB, Diverse Cohorts | keep `pheno` |

LATEST_DX_MAP = {
    "Parkinson's Disease": "PD",
    "Idiopathic PD": "PD",
    "LBD": "DLB",
    "Dementia With Lewy Bodies": "DLB",
    "Progressive Supranuclear Palsy": "PSP",
    "No PD Nor Other Neurological Disorder": "control",
}

# ROSMAP clinical consensus codes: 1=NCI, 2/3=MCI, 4/5=AD-dementia, 6=other
DCFDX_MAP = {"1": "control", "2": "MCI", "3": "MCI",
             "4": "other", "5": "other", "6": "other"}


def dx_amppd(r, pheno):
    return LATEST_DX_MAP.get(r["diagnosis_latest"], "other")


def dx_rosmap(r, pheno):
    return DCFDX_MAP.get(r["dcfdx_lv"], "")


def dx_mayo(r, pheno):
    if r["diagnosis"] == "progressive supranuclear palsy":
        return "PSP"
    return pheno


def dx_same_as_pheno(r, pheno):
    return pheno


RULES = {
    "amp_ad_rosmap": (pheno_braak_cerad, dx_rosmap),
    "amp_ad_mayo":   (pheno_mayo,        dx_mayo),
    "amp_ad_msbb":   (pheno_braak_cerad, dx_same_as_pheno),
    "amp_ad_divco":  (pheno_divco,       dx_same_as_pheno),
    "amp_pd":        (pheno_amppd,       dx_amppd),
}


def derive(df, src):
    pheno_fn, dx_fn = RULES[src]
    phenos, dxs = [], []
    for r in df.to_dict("records"):
        p = pheno_fn(r)
        phenos.append(p)
        dxs.append("AD" if p == "AD" else dx_fn(r, p))   # step 1: AD pinned
    return phenos, dxs


print("dx_detailed rules defined")


# ## 5. Derive the labels
#
# Distributions are printed for checking against §2's raw counts.

labels = {}
for src, df in raw.items():
    p, d = derive(df, src)
    labels[src] = pd.DataFrame({"pheno": p, "dx_detailed": d})
    print("=" * 70)
    print(src)
    for col in ("pheno", "dx_detailed"):
        vc = labels[src][col].replace("", "(null)").value_counts()
        print(f"  {col:12} " + "  ".join(f"{k}={v:,}" for k, v in vc.items()))


# ## 6. One donor table
#
# `core` — one row per `(individual_id, source_dataset)`, with `sex` in PLINK2 coding (1/2/0) and
# `projid` for ROSMAP only. The key is the **pair**: Diverse Cohorts re-enrolls donors from the
# other AMP-AD studies and each cohort is sequenced separately, so deduplicating on
# `individual_id` would discard genomes. Sex words differ in case across sources, hence `.lower()`.

CORE_COLUMNS = ["individual_id", "source_dataset", "sex", "projid", "pheno", "dx_detailed"]

SEX = {"male": "1", "female": "2", "missing or unknown": "0"}

frames = []
for src, df in raw.items():
    sex = df["sex"].str.lower().map(SEX).fillna("")
    unmapped = df.loc[sex == "", "sex"].value_counts().to_dict()

    frames.append(pd.DataFrame({
        "individual_id": df[ID_COL[src]],
        "source_dataset": src,
        "sex": sex,
        "projid": df["projid"] if "projid" in df.columns else "",
        "pheno": labels[src]["pheno"].values,
        "dx_detailed": labels[src]["dx_detailed"].values,
    })[CORE_COLUMNS])

    print(f"{src:16} rows={len(df):7,}  unmapped sex: {unmapped or '-'}")

core = pd.concat(frames, ignore_index=True)
print(f"\ncore: {len(core):,} rows x {core.shape[1]} cols")
print(core.source_dataset.value_counts().to_string())

for col in ("sex", "pheno", "dx_detailed"):
    print(f"=== {col} x source_dataset ===")
    print(pd.crosstab(core[col].replace("", "(null)"), core.source_dataset, margins=True).to_string())
    print()


# ## 7. Genotype samples
#
# `genomes` — one row per `.psam` sample across the four callsets, resolved to a donor. A callset
# is not a cohort (`wgs_harm` joint-calls ROSMAP + Mayo + MSBB), so each row carries both
# `source_callset` and `source_dataset`. Resolution rules per callset: METHODS.md §2.
#
# All three `wgs_harm` rules run on every sample, so a collision between studies surfaces as
# `n_rules_hit > 1` instead of being settled by rule order. A missing `.psam` is skipped with a note.

projid2donor = dict(zip(raw["amp_ad_rosmap"]["projid"], raw["amp_ad_rosmap"]["individualID"]))
mayo_donors = set(raw["amp_ad_mayo"]["individualID"])
divco_donors = set(raw["amp_ad_divco"]["individualID"])

rosmap_qc = rd(AD_META / "WGS_sample_QC_info.csv")   # CR-only line terminators; pandas handles it
qc_by_sample = {r["WGS_id"]: r for r in rosmap_qc.to_dict("records")}

msbb_bio = rd(AD_META / "MSBB_biospecimen_metadata.csv")
msbb_by_specimen = {r["specimenID"]: r for r in
                    msbb_bio[msbb_bio["assay"] == "wholeGenomeSeq"]
                    .drop_duplicates("specimenID").to_dict("records")}

divco_bio = rd(DIVCO_META / "AMP-AD_DiverseCohorts_biospecimen_metadata.csv")
divco_by_specimen = {r["specimenID"]: r
                     for r in divco_bio.drop_duplicates("specimenID").to_dict("records")}

br_inventory = rd(AMPPD_META / "wgs_BR-DSNWGS_sample_inventory.csv")
br_sample2donor = dict(zip(br_inventory["sample_id"], br_inventory["participant_id"]))


# Each resolve function: sample id -> list of (source_dataset, rule, individual_id, specimenID, tissue).
# Empty list = unresolved. More than one = an id claimed by two studies.

def resolve_wgs_harm(iid):
    hits = []
    if iid in qc_by_sample:
        q = qc_by_sample[iid]
        hits.append(("amp_ad_rosmap", "rosmap:wgs_qc_manifest",
                     projid2donor.get(q["projid"], ""), "", q["Source.Tissue.Type"]))
    if iid in mayo_donors:
        hits.append(("amp_ad_mayo", "mayo:identity", iid, "", ""))
    if iid in msbb_by_specimen:
        b = msbb_by_specimen[iid]
        hits.append(("amp_ad_msbb", "msbb:biospecimen_specimenID",
                     b["individualID"], iid, b["tissue"] or b["organ"]))
    return hits


def resolve_divco(iid):
    if iid in divco_donors:
        return [("amp_ad_divco", "divco:individualID_direct", iid, "", "")]
    parts = iid.split("_")
    if len(parts) > 1 and parts[0] in divco_donors:
        return [("amp_ad_divco", "divco:specimenID_suffix_strip", parts[0], iid,
                 parts[1] if len(parts) >= 3 else "")]
    if iid in divco_by_specimen:
        b = divco_by_specimen[iid]
        return [("amp_ad_divco", "divco:biospecimen_lookup", b["individualID"], iid,
                 b["tissue"] or b["organ"])]
    return []


def resolve_wb_dwgs(iid):
    return [("amp_pd", "wb_dwgs:psam_identity", iid, "", "whole_blood")]


def resolve_br_dsnwgs(iid):
    donor = br_sample2donor.get(iid, "")
    return [("amp_pd", "br_dsnwgs:sample_inventory", donor, iid, "brain")] if donor else []


# Raw psams (config.sh RAW_*), read in place — never step 1's `_sexupd` output, which this feeds.
CALLSETS = [
    {"name": "wgs_harm",  "psam": DIR_WGS / "pgen/wgs_harm_hg38.psam",
     "iid": "#IID", "fid": None,   "resolve": resolve_wgs_harm},
    {"name": "divco_hs",  "psam": DIR_DC / "pgen/divco_hs_hg38.psam",
     "iid": "#IID", "fid": None,   "resolve": resolve_divco},
    {"name": "wb_dwgs",   "psam": DIR_WB / "joint_calls/all_chrs_merged.psam",
     "iid": "IID",  "fid": "#FID", "resolve": resolve_wb_dwgs},
    {"name": "br_dsnwgs", "psam": DIR_BR / "pgen/br_dsnwgs_hg38.psam",
     "iid": "#IID", "fid": None,   "resolve": resolve_br_dsnwgs},
]

GENOME_COLUMNS = ["IID", "source_callset", "individual_id", "source_dataset",
                  "rule", "specimenID", "tissue", "n_rules_hit"]

psams, rows = {}, []
for cs in CALLSETS:
    if not cs["psam"].exists():
        print(f"{cs['name']:12} SKIPPED — no psam at {rel(cs['psam'])}")
        continue
    psam = rd(cs["psam"], sep="\t")
    psams[cs["name"]] = psam
    for iid in psam[cs["iid"]]:
        hits = cs["resolve"](iid)
        src, rule, donor, spec, tissue = hits[0] if hits else ("", "unresolved", "", "", "")
        rows.append({"IID": iid, "source_callset": cs["name"], "individual_id": donor,
                     "source_dataset": src, "rule": rule, "specimenID": spec,
                     "tissue": tissue, "n_rules_hit": len(hits)})

# Explicit columns so zero psams still yields the right schema for §8's checks.
genomes = pd.DataFrame(rows, columns=GENOME_COLUMNS)

print(f"\ngenome samples: {len(genomes):,}")
print(pd.crosstab(genomes.source_callset,
                  genomes.source_dataset.replace("", "(unresolved)"), margins=True).to_string())
print("\nby rule:")
print(genomes.rule.value_counts().to_string())


# ## 8. Checks
#
# Every row is a violation count; anything nonzero stops the write in §9. Vocabularies come from
# clinical_common.

core_keys = set(zip(core.individual_id, core.source_dataset))
resolved = genomes[genomes.individual_id != ""]

checks = pd.DataFrame([
    ("core: (individual_id, source_dataset) unique",
     core.duplicated(["individual_id", "source_dataset"]).sum()),
    ("core: individual_id never blank", (core.individual_id == "").sum()),
    ("core: sex within vocabulary", (~core.sex.isin(SEX_VALUES + [""])).sum()),
    ("core: pheno within vocabulary", (~core.pheno.isin(PHENO_VALUES + [""])).sum()),
    ("core: dx_detailed within vocabulary", (~core.dx_detailed.isin(DX_VALUES + [""])).sum()),
    ("core: AD identical across pheno and dx_detailed",
     ((core.pheno == "AD") != (core.dx_detailed == "AD")).sum()),
    ("genomes: (IID, source_callset) unique",
     genomes.duplicated(["IID", "source_callset"]).sum()),
    ("genomes: no sample claimed by >1 study rule", (genomes.n_rules_hit > 1).sum()),
    ("genomes: every resolved genome lands on a core donor",
     sum(k not in core_keys for k in zip(resolved.individual_id, resolved.source_dataset))),
], columns=["check", "violations"])

print(checks.to_string(index=False))

unresolved = genomes[genomes.individual_id == ""]
if len(unresolved):
    print(f"\nunresolved samples: {len(unresolved):,}")
    print(unresolved.groupby("source_callset").size().to_string())

xsrc = core.groupby("individual_id")["source_dataset"].nunique()
print(f"\n{int((xsrc > 1).sum()):,} donors are enrolled in more than one cohort "
      f"(expected — §12 reconciles them):")
print(core[core.individual_id.isin(xsrc[xsrc > 1].index)]
      .groupby("individual_id")["source_dataset"]
      .apply(lambda s: " + ".join(sorted(s))).value_counts().to_string())

coverage = pd.DataFrame({
    "donors": core.groupby("source_dataset").size(),
    "genomes": resolved.groupby("source_dataset").size(),
    "donors_with_genome": resolved.groupby("source_dataset")["individual_id"].nunique(),
}).fillna(0).astype(int)
coverage["pct"] = (coverage.donors_with_genome / coverage.donors).map("{:.1%}".format)
print("\ncoverage:")
print(coverage.to_string())


# ## 9. Write the audit tables
#
# The record of how each sample resolved to a donor; analysis_grain.py reads both.

if checks.violations.sum():
    raise AssertionError("invariants violated:\n"
                         + checks.query("violations > 0").to_string(index=False))

OUT.mkdir(parents=True, exist_ok=True)
core.to_csv(OUT / "individual_core.csv", index=False)
genomes.drop(columns=["n_rules_hit"]).to_csv(OUT / "genome_crosswalk.csv", index=False)

print(f"individual_core.csv   {len(core):,} rows")
print(f"genome_crosswalk.csv  {len(genomes):,} rows")
print(f"\nwritten to {rel(OUT)}")


# ## 10. Out → sex-update files
#
# One per callset, applied by `01_genotools.sh` (`plink2 --update-sex`) before any filtering.
# **Sex is per-cohort, never reconciled**: disagreements are emitted as-is so step 1's sex check
# can adjudicate them. Unknown sex is left unwritten. Columns follow each `.psam`.

sex_by_key = {(r.individual_id, r.source_dataset): r.sex for r in core.itertuples()}

n_conflict = int((core[core.sex.isin(["1", "2"])].groupby("individual_id").sex.nunique() > 1).sum())
print(f"cross-source sex conflicts: {n_conflict} donor(s) — emitted per-cohort, "
      f"step 1's sex-check adjudicates\n" + "=" * 72)

for cs in CALLSETS:
    if cs["name"] not in psams:
        continue
    psam = psams[cs["name"]]
    g = genomes[genomes.source_callset == cs["name"]].set_index("IID")
    sex = psam[cs["iid"]].map(
        lambda i: sex_by_key.get((g.individual_id.get(i, ""), g.source_dataset.get(i, "")), ""))
    keep = sex.isin(["1", "2"])

    out = pd.DataFrame({cs["iid"]: psam[cs["iid"]], "SEX": sex})
    cols = [cs["iid"], "SEX"]
    if cs["fid"]:
        out.insert(0, cs["fid"], psam[cs["fid"]])
        cols = [cs["fid"]] + cols

    out.loc[keep, cols].to_csv(OUT / f"{cs['name']}_update_sex.txt", sep="\t", index=False)

    print(f"[{cs['name']:10}] samples {len(psam):6,}  sex written {int(keep.sum()):6,} "
          f"(male {int((sex == '1').sum()):,}, female {int((sex == '2').sum()):,})  "
          f"left as-is {int((~keep).sum()):,}")

    # No check against the psam's SEX: raw callsets ship it as NA. Step 1's sex check covers it.

print(f"\nwritten to {rel(OUT)} — step 1 reads them from here directly (config.sh: sex_file)")


# ## 12a. Out → `sample_annot.csv` — the grain's PC-free half
#
# One row per genome: `IID, source_callset, pheno, dx_detailed`, reconciled exactly as the grain
# is. Step 6's AF-concordance stage needs only these, so writing them before any genotype step
# lets step 6 run as one pass. Runs **unconditionally**.
#
# Where a grain exists, every shared IID must agree on both columns; a disagreement would change
# step 6's exclusion list for a non-genotype reason.

by_genome_annot = {}
for r in genomes.itertuples():
    by_genome_annot.setdefault(r.IID, []).append(r)

by_donor_annot = {}
for r in core.to_dict("records"):
    by_donor_annot.setdefault(r["individual_id"], []).append(r)

annot_rows = []
for iid, hits in by_genome_annot.items():
    donor = next((h.individual_id for h in hits if h.individual_id), "")
    rows = by_donor_annot.get(donor, []) if donor else []
    pheno, dx, _, _ = reconcile(rows) if rows else ("", "", False, False)
    annot_rows.append({
        "IID": iid,
        "source_callset": "|".join(sorted({h.source_callset for h in hits})),
        "pheno": pheno,
        "dx_detailed": dx})

annot = pd.DataFrame(annot_rows)[["IID", "source_callset", "pheno", "dx_detailed"]]
annot = annot.sort_values("IID").reset_index(drop=True)
annot.to_csv(OUT / "sample_annot.csv", index=False)

print(f"sample_annot.csv  {len(annot):,} genomes -> {rel(OUT)}")
print("\ndx_detailed x source_callset:")
print(pd.crosstab(annot.dx_detailed.replace("", "(null)"),
                  annot.source_callset, margins=True).to_string())

# Cells the AF filter could use (>=100 per callset per dx); strata are set in step 6, so an upper bound.
print("\ncells with >=100 genomes per callset (upper bound — stratum splits these further):")
_pairs = (annot[annot.dx_detailed != ""].groupby(["dx_detailed", "source_callset"]).size())
for dx_v, grp in _pairs.groupby(level=0):
    powered = [f"{cs}={n:,}" for (_, cs), n in grp.items() if n >= 100]
    if len(powered) >= 2:
        print(f"    {dx_v:10} {', '.join(powered)}")

# ── regression check against the grain, where one exists ──
_grain_path = OUT / "analysis_grain.csv"
if _grain_path.exists():
    _g = pd.read_csv(_grain_path, dtype=str, keep_default_na=False)
    _cmp = _g[["IID", "source_callset", "dx_detailed"]].merge(
        annot[["IID", "source_callset", "dx_detailed"]], on="IID",
        how="left", suffixes=("_grain", "_annot"), indicator=True)
    _missing = int((_cmp._merge == "left_only").sum())
    _both = _cmp[_cmp._merge == "both"]
    _d_cs = int((_both.source_callset_grain != _both.source_callset_annot).sum())
    _d_dx = int((_both.dx_detailed_grain != _both.dx_detailed_annot).sum())
    print(f"\nvs analysis_grain.csv: {len(_both):,} shared IIDs, "
          f"{_d_cs} callset mismatches, {_d_dx} dx mismatches, {_missing} grain IIDs absent here")
    if _d_cs or _d_dx or _missing:
        print("  !! The two reconciliation paths DISAGREE. Step 6's exclusion list would change")
        print("     for a non-genotype reason — investigate before submitting step 6.")
    else:
        print("  identical on both columns — the AF build is unaffected by the split")


# ## 14. Appendix — how much of the AD arm is definition-dependent?
#
# **Read-only, not part of the build.** How different would the AD arm be if each cohort's
# non-derived diagnosis column were taken at face value? Only ROSMAP's `dcfdx_lv` is clinical;
# Mayo `diagnosis` and DivCo `reag` are other pathology calls; MSBB `CDR` is severity, not
# etiology, so its AD count is an upper bound. METHODS.md §10.

AD_SOURCES = ["amp_ad_rosmap", "amp_ad_mayo", "amp_ad_msbb", "amp_ad_divco"]
ALT_COL = {"amp_ad_rosmap": "dcfdx_lv", "amp_ad_mayo": "diagnosis",
           "amp_ad_msbb": "CDR", "amp_ad_divco": "reag"}

for src in AD_SOURCES:
    col = ALT_COL[src]
    print("=" * 78)
    print(f"{src}   derived pheno (rows)  x  {col} (cols)")
    print(pd.crosstab(labels[src]["pheno"].replace("", "(null)"),
                      raw[src][col].replace("", "(blank)").values, margins=True).to_string())
    print()

# One explicit alternative-AD call per cohort; each is a judgment, MSBB's the weakest.

def alt_ad_rosmap(v):      # clinical AD-dementia (dcfdx_lv 4/5)
    return "AD" if v in {"4", "5"} else ("" if v in MISSING else "not_AD")

def alt_ad_mayo(v):        # pre-made neuropathological consensus
    return "AD" if v == "Alzheimer Disease" else ("" if v in MISSING else "not_AD")

def alt_ad_msbb(v):        # CDR >= 1 = dementia of ANY cause -> upper bound, not a diagnosis
    if v in MISSING:
        return ""
    try:
        return "AD" if float(v) >= 1 else "not_AD"
    except ValueError:
        return ""

def alt_ad_divco(v):       # NIA-Reagan high likelihood
    return "AD" if v == "High Likelihood" else ("" if v in MISSING else "not_AD")

ALT_FN = {"amp_ad_rosmap": alt_ad_rosmap, "amp_ad_mayo": alt_ad_mayo,
          "amp_ad_msbb": alt_ad_msbb, "amp_ad_divco": alt_ad_divco}

has_genome = set(zip(resolved.individual_id, resolved.source_dataset))

rows = []
for src in AD_SOURCES:
    ids = raw[src][ID_COL[src]]
    derived = labels[src]["pheno"].values
    alt = raw[src][ALT_COL[src]].map(ALT_FN[src]).values
    geno = [(i, src) in has_genome for i in ids]

    d_ad, a_ad = derived == "AD", alt == "AD"
    rows.append({
        "cohort": src.replace("amp_ad_", ""), "instrument": ALT_COL[src],
        "derived_AD": int(d_ad.sum()), "alt_AD": int(a_ad.sum()),
        "both": int((d_ad & a_ad).sum()),
        "derived_only": int((d_ad & ~a_ad).sum()), "alt_only": int((a_ad & ~d_ad).sum()),
        "alt_unusable": int((alt == "").sum()), "derived_null": int((derived == "").sum()),
        "derived_AD_geno": int((d_ad & geno).sum()), "alt_AD_geno": int((a_ad & geno).sum()),
    })

summary = pd.DataFrame(rows).set_index("cohort")
print("AD cases under each definition (donors):")
print(summary[["instrument", "derived_AD", "alt_AD", "both", "derived_only", "alt_only"]].to_string())
print("\nunusable under each definition (donors the rule cannot classify):")
print(summary[["derived_null", "alt_unusable"]].to_string())
print("\nAD cases WITH a genome — the number that sets GWAS power:")
print(summary[["derived_AD_geno", "alt_AD_geno"]].to_string())

tot_d, tot_a = summary.derived_AD_geno.sum(), summary.alt_AD_geno.sum()
print(f"\ntotal AD genomes: derived={tot_d:,}  alternative={tot_a:,}  "
      f"difference={tot_a - tot_d:+,}")


# ### 14a. Are we losing samples to the neuropathological definition?
#
# Null donors overstate the loss; null genomes are what cost power. Both are counted.

ad = core[core.source_dataset.isin(AD_SOURCES)].copy()
ad["has_genome"] = [(i, s) in has_genome for i, s in zip(ad.individual_id, ad.source_dataset)]

for label, g in (("all donors", ad), ("donors with >=1 genome", ad[ad.has_genome])):
    print(f"--- pheno x cohort, {label} ---")
    print(pd.crosstab(g.pheno.replace("", "(null)"), g.source_dataset, margins=True).to_string())
    print()

n_donor = int((ad.pheno == "").sum())
n_geno = int((ad.loc[ad.has_genome, "pheno"] == "").sum())
tot_geno = int(ad.has_genome.sum())
print("unclassifiable by the neuropath rule:")
print(f"  {n_donor:,} / {len(ad):,} donors  ({n_donor/len(ad):.1%})")
print(f"  {n_geno:,} / {tot_geno:,} genomes ({n_geno/tot_geno:.1%})   <- the number that matters")


# ### 14b. Could the alternative instrument recover the ones we lose?
#
# For genome-carrying donors with a null `pheno`: how many the alternative column classifies,
# and how many as AD.

recovery = []
for src in AD_SOURCES:
    col, df = ALT_COL[src], raw[src]
    null_geno = [(p == "") and ((i, src) in has_genome)
                 for p, i in zip(labels[src]["pheno"], df[ID_COL[src]])]
    n = sum(null_geno)
    vals = df.loc[null_geno, col]
    usable = int((~vals.isin(MISSING)).sum()) if n else 0
    as_ad = int(vals.map(ALT_FN[src]).eq("AD").sum()) if n else 0
    print(f"{src.replace('amp_ad_', ''):8} {n:4} unclassifiable genomes | "
          f"{col} classifies {usable}, of which AD = {as_ad}")
    recovery.append({"cohort": src.replace("amp_ad_", ""), "null_genomes": n,
                     "recoverable": usable, "as_AD": as_ad})

rec = pd.DataFrame(recovery).set_index("cohort")
print("\n" + rec.to_string())

ad_now = int(((ad.pheno == "AD") & ad.has_genome).sum())
gain = int(rec.as_AD.sum())
print(f"\nAD genomes today: {ad_now:,}")
print(f"recoverable as AD: {gain:,}  ->  {ad_now + gain:,} "
      f"({gain/ad_now:+.1%} change in AD arm size)")
print(f"total genomes recoverable into SOME arm: {int(rec.recoverable.sum()):,}")
