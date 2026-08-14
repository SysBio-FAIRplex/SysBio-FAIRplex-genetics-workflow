#!/usr/bin/env python
# coding: utf-8

# # Harmonized clinical core
#
# Turns per-cohort clinical files into the phenotype inputs the GWAS needs.
#
# The pipeline alternates between this notebook and the cluster. Each cluster round trip
# is one section here; sections waiting on cluster output skip cleanly and say so.
#
# | § | Does | Needs |
# |---|---|---|
# | 1–6 | read clinical files, derive `pheno` / `dx_detailed` per donor | local |
# | 7–9 | resolve genotype samples to donors, check, write audit tables | psam files |
# | 10 | **out →** per-callset sex-update files | — |
# | 11 | **in ←** genotools/relatedness QC outcomes, labelled by reason | steps 1–5 |
# | 12 | **in ←** ancestry + PCs, reconciled to one label per genome | step 6 |
# | 13 | **out →** per-ancestry covariate and per-contrast phenotype files | step 6 |
# | 14 | appendix: how definition-dependent is the AD arm (read-only) | — |
#
# Five cohorts — ROSMAP, Mayo, and MSBB (one AMP-AD harmonized set), AMP-AD Diverse
# Cohorts, and AMP-PD — across four genotype callsets.
#
# **Phenotype criteria** follow the AMP-AD diagnosis-criteria specification (Synapse
# syn51757663; harmonized dictionary syn73713784). The governing principle is that
# **AD is defined by neuropathology**.
#
# **Guardrail.** Every cell prints aggregates. No cell prints a subject-level row.
# Genotype files are read for sample IDs only and are never modified.

# In[1]:


# synapse get syn51757644 --downloadLocation data/amp-ad-genomics/DivCo_HS/metadata/
# synapse get syn51757645 --downloadLocation data/amp-ad-genomics/DivCo_HS/metadata/
# synapse get syn51757646 --downloadLocation data/amp-ad-genomics/DivCo_HS/metadata/

# synapse get syn73713768 --downloadLocation data/amp-ad-genomics/WGS_Harmonization/metadata/
# synapse get syn12178037 --downloadLocation data/amp-ad-genomics/WGS_Harmonization/metadata/
# synapse get syn73713767 --downloadLocation data/amp-ad-genomics/WGS_Harmonization/metadata/
# synapse get syn21893059 --downloadLocation data/amp-ad-genomics/WGS_Harmonization/metadata/
# synapse get syn73713766 --downloadLocation data/amp-ad-genomics/WGS_Harmonization/metadata/


# In[2]:


# gcloud storage cp gs://amp-pd-data/releases/2023_v4release_1027/clinical/Demographics.csv data/amp-pd-genomics/metadata/ --billing-project 8641313829
# gcloud storage cp gs://amp-pd-data/releases/2023_v4release_1027/amp_pd_participants.csv data/amp-pd-genomics/metadata/ --billing-project 8641313829
# gcloud storage cp gs://amp-pd-data/releases/2023_v4release_1027/amp_pd_case_control.csv data/amp-pd-genomics/metadata/ --billing-project 8641313829
# gcloud storage cp gs://amp-pd-data/releases/2023_v4release_1027/wgs_BR-DSNWGS_sample_inventory.csv data/amp-pd-genomics/metadata/ --billing-project 8641313829


# In[3]:


from pathlib import Path

import pandas as pd

# Every path below is anchored to this file, so the script runs identically from any
# working directory and nothing outside this folder is ever read or written. No absolute
# paths, no environment variables, no cluster paths: to move the analysis, copy the folder.
PROJECT_ROOT = Path(__file__).resolve().parent

DATA = PROJECT_ROOT / "data"
OUT = PROJECT_ROOT / "clinical_core_out"

# Per-callset roots. Same names, same layout as config.sh's DIR_* — one vocabulary across
# the clinical and genotype sides. Note the asymmetry config.sh also documents: the two
# AMP-AD callsets keep their pgens in pgen/, WB-DWGS in joint_calls/.
DIR_WGS = DATA / "amp-ad-genomics/WGS_Harmonization"
DIR_DC = DATA / "amp-ad-genomics/DivCo_HS"
DIR_WB = DATA / "amp-pd-genomics/WB-DWGS"
DIR_BR = DATA / "amp-pd-genomics/BR-DSNWGS"

# Clinical metadata — downloaded from Synapse / GCP (provenance at the top of this file).
AD_META = DIR_WGS / "metadata"
DIVCO_META = DIR_DC / "metadata"
AMPPD_META = DATA / "amp-pd-genomics/metadata"
WB_META = DIR_WB / "metadata"
BR_META = DIR_BR / "metadata"

# Written by the cluster, rsynced back down into this folder between round trips
# (see §11, §12, §13). Absent on a fresh copy — those sections skip and say so.
MERGED = DATA / "merged"
EXCLUDE_REASONS = MERGED / "relatedness/exclude_reasons.tsv"
RETAINED_MANIFEST = MERGED / "relatedness/retained_manifest.csv"
PCA_DIR = MERGED / "by_ancestry_qc"


def rel(path):
    """Display form: a path shown relative to PROJECT_ROOT, so a run's log reads the
    same on every machine."""
    return Path(path).relative_to(PROJECT_ROOT)


def rd(path, sep=","):
    """Read everything as stripped strings. Missing stays "" (never NaN) so nulls
    show up in value_counts and crosstabs instead of being silently dropped."""
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    return df.apply(lambda s: s.str.strip())


pd.set_option("display.width", 130)


# ## 1. Clinical sources
#
# Seven files across five cohorts, each at individual grain. ROSMAP, Mayo, MSBB, and
# Diverse Cohorts are one file each; AMP-PD takes three, joined on `participant_id`.
# ROSMAP/Mayo/MSBB share the AMP-AD harmonized schema; the other two use their own.

# In[4]:


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
          f"unique ids: {ids.nunique():,}   e.g. {ids.iloc[0]}")


# ## 2. The raw phenotype instruments
#
# The columns the derivation reads. Every label produced later traces back to one of these.
#
# The cohorts do not share an instrument. ROSMAP and MSBB carry Braak staging and CERAD
# plaque scores; Mayo carries Braak and Thal phase plus a neuropathological `diagnosis`;
# Diverse Cohorts ships a pre-adjudicated outcome; AMP-PD is clinical throughout and has
# no neuropathology at all.

# In[5]:


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
# `pheno` is the primary case/control label: `{AD, PD, control, other}`, or null.
#
# **AD is neuropathological.** Clinical diagnosis codes never define AD.
#
# **Missing inputs yield null, not `other`.** `other` means classified and neither case
# nor control; null means not classifiable.
#
# | Cohort | Rule |
# |---|---|
# | ROSMAP, MSBB | Braak + CERAD. AD = Braak≥IV & CERAD moderate/frequent; control = Braak≤III & CERAD sparse/none; else other |
# | Mayo | Braak + Thal. AD = Braak≥IV & Thal≥2; control = Braak≤III & Thal<2; else other — then the control-purity screen below |
# | Diverse Cohorts | pre-adjudicated `ADoutcome`, or `mayoDx` for the Mayo contribution group |
# | AMP-PD | curated `case_control_other_latest`: Case→PD, Control→control, Other excluded |
#
# **Mayo's control-purity screen.** Braak and Thal are both AD-specific axes, so a PSP
# brain — a 4R-tauopathy with little amyloid — scores as a control on both and would
# contaminate the control arm. Any neuropathological control whose `diagnosis` names a
# disease is demoted to `other`. The screen is one-directional: it never rescues a null
# and never overrides AD. It applies to Mayo alone, because Mayo is the only cohort with
# an independent neuropathological diagnosis alongside the staging axes.
#
# **Why Diverse Cohorts keeps its pre-made call** rather than recomputing from Braak+Thal:
# the harmonized Thal column is a lossy subset of the plaque data the original
# adjudication rested on, so recomputing would null out donors who have Braak but no Thal.
#
# **Never Braak alone** for any cohort — it misclassifies amyloid-negative age-related
# tauopathy as AD.

# In[6]:


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
# A second label for secondary analyses: `{PD, AD, MCI, DLB, PSP, control, other}`, or null.
#
# **Step 1, every cohort: if `pheno` is AD, `dx_detailed` is AD.** This pins the AD set
# to be identical across both fields — a donor cannot be AD in one analysis and MCI in
# another.
#
# **Step 2, otherwise, use the cohort's specific instrument:**
#
# | Cohort | Instrument |
# |---|---|
# | AMP-PD | `diagnosis_latest` → PD / DLB / PSP / control; everything else other |
# | ROSMAP | `dcfdx_lv` → 1 control, 2–3 MCI, 4–6 other |
# | Mayo | `diagnosis` == progressive supranuclear palsy → PSP; else keep `pheno` |
# | MSBB, Diverse Cohorts | keep `pheno` — no independent specific diagnosis available |
#
# ROSMAP codes 4 and 5 are clinical AD-dementia. They map to `other`, not AD, because AD
# here is neuropathological. `dcfdx_lv` is used over `cogdx` because it is far more complete.
#
# The two fields disagree for non-AD labels by design — they use different instruments.
# An analysis that needs the strict set should use `pheno`.

# In[7]:


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
# Applied per cohort, with the resulting distributions printed so each one can be checked
# against the raw instrument counts in §2.

# In[8]:


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
# `core` — one row per `(individual_id, source_dataset)`.
#
# | Column | Notes |
# |---|---|
# | `individual_id` | native donor id, unnamespaced so it joins back to source files |
# | `source_dataset` | cohort |
# | `sex` | PLINK2 coding — `1` male, `2` female, `0` unknown |
# | `projid` | ROSMAP's donor identifier; blank elsewhere |
# | `pheno` | primary label |
# | `dx_detailed` | secondary label |
#
# The key is the **pair**, not `individual_id` alone: Diverse Cohorts re-enrolls donors
# from the other AMP-AD studies, and each cohort is sequenced separately, so both rows
# correspond to real data. Deduplicating on `individual_id` would discard genomes.
#
# All five sources record sex as a word in a column named `sex` — lowercase in the AMP-AD
# files, title-case in AMP-PD — so the map reads a lowercased value.

# In[9]:


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


# In[10]:


for col in ("sex", "pheno", "dx_detailed"):
    print(f"=== {col} x source_dataset ===")
    print(pd.crosstab(core[col].replace("", "(null)"), core.source_dataset, margins=True).to_string())
    print()


# ## 7. Genotype samples
#
# `genomes` — one row per `.psam` sample across all four callsets, resolved to a donor.
#
# A callset is not a cohort. `wgs_harm` is a joint call over ROSMAP + Mayo + MSBB, and the
# same ROSMAP donor can also appear in `divco_hs`. So every row carries both
# `source_callset` (which file it was sequenced and joint-called in) and `source_dataset`
# (which study phenotyped the donor).
#
# | Callset | Resolution |
# |---|---|
# | `wgs_harm` | ROSMAP: `#IID` == `WGS_id` → `projid` → `individualID`. Mayo: identity. MSBB: `#IID` == `specimenID` → `individualID` |
# | `divco_hs` | direct id, else strip tissue suffix (`<individualID>_DLPFC_WGS`→`<individualID>`), else biospecimen lookup |
# | `wb_dwgs` | `IID` == donor id (identity) |
# | `br_dsnwgs` | `sample_id` → `participant_id` via the AMP-PD sample inventory |
#
# All three `wgs_harm` rules are tested against **every** sample rather than
# short-circuiting, so a collision between studies surfaces as `n_rules_hit > 1` instead of
# being settled silently by rule order.
#
# Adding a callset is one entry in `CALLSETS` plus one resolve function. A callset whose
# `.psam` is not down yet is skipped with a note.

# In[11]:


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


# Each psam is the RAW genotype fileset's own — the same file config.sh names as RAW_*, read
# in place. Nothing is copied into metadata/. That matters for wb_dwgs in particular: the
# previous path read `all_chrs_merged_sexupd.psam`, which is step 1's OUTPUT, so this section
# was reading a file produced by the pipeline it feeds. Raw psams change only when a callset
# is rebuilt, so there is no round trip to keep in sync.
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

# AMP-PD callsets — drives the @amppd / @ampad contrast tags in §13.
AMPPD_CALLSETS = {"wb_dwgs", "br_dsnwgs"}

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

# Columns given explicitly so that zero available psams yields an empty frame with the right
# schema rather than a column-less one, which would turn §8's checks into an AttributeError.
genomes = pd.DataFrame(rows, columns=GENOME_COLUMNS)

print(f"\ngenome samples: {len(genomes):,}")
print(pd.crosstab(genomes.source_callset,
                  genomes.source_dataset.replace("", "(unresolved)"), margins=True).to_string())
print("\nby rule:")
print(genomes.rule.value_counts().to_string())


# ## 8. Checks
#
# Every row is a violation count. Zero is a pass; anything else stops the write in §9.

# In[12]:


SEX_VALUES = ["0", "1", "2"]
PHENO_VALUES = ["AD", "PD", "control", "other"]
DX_VALUES = ["PD", "AD", "MCI", "DLB", "PSP", "control", "other"]

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
# Neither file is read by the pipeline — they are the durable record of how a sample
# resolved to a donor, and the basis for anything produced later.

# In[13]:


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
# One file per callset, applied by `01_genotools.sh` via `plink2 --update-sex` before it
# filters anything.
#
# **Sex is per-cohort, never reconciled.** A genome takes the sex its *own* source study
# recorded. Cross-source disagreements are emitted as-is — the genotype sex-check in step 1
# adjudicates them, and it can only do that if it sees the claimed value. Only unknown sex
# (`0` or blank) is left unwritten, so those samples keep whatever their `.psam` says.
#
# No sample is removed here. Column layout follows each `.psam`: callsets without an FID
# get `#IID SEX`, `wb_dwgs` gets `#FID IID SEX`.

# In[14]:


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

    # No comparison against the psam's own SEX column: every raw callset ships SEX as NA, so
    # there is nothing to disagree with. (The check that used to live here read wb_dwgs's
    # _sexupd psam — step 1's output — and so only ever confirmed that plink had applied the
    # file we are writing right now.) Genotype-inferred sex is checked against these values by
    # genotools in step 1, which is where the finding belongs.

print(f"\nwritten to {rel(OUT)} — step 1 reads them from here directly (config.sh: sex_file)")


# ## 11. In ← QC outcomes
#
# After genetics steps 1–5, every sample has been kept or dropped for a stated reason:
# genotools per-ancestry QC (call rate, sex-check, heterozygosity), duplicate resolution,
# or 2nd-degree relatedness pruning. This attaches those outcomes to the clinical table so
# a loss can be read per cohort and per phenotype arm rather than per callset alone.
#
# Reads two files produced by `05_excludelist.py`:
#
# | File | Columns |
# |---|---|
# | `exclude_reasons.tsv` | `FID  IID  reason  detail` |
# | `retained_manifest.csv` | `FID,IID,ancestry,call_rate,dup_cluster_id` |
#
# The picks themselves stay on the genotype side — they are call-rate-driven and operate on
# KING kinship output. This section consumes them; it does not recompute them.

# In[15]:


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
# One row per retained genome, carrying the label the GWAS tests and the covariates it
# adjusts for. Joins the retained manifest, the per-ancestry PCs from step 6, the genome
# crosswalk from §7, and the donor table from §6.
#
# **Why reconciliation is needed.** A donor can appear in two source studies — once in
# Diverse Cohorts and once in the ROSMAP/Mayo/MSBB trio — phenotyped by two different
# programs that do not always agree. The grain is per genome, so each genome needs one label.
#
# | Field | Rule |
# |---|---|
# | `pheno` | agree → that value; exactly one labelled → the labelled one; both labelled and differing → **AD-dominant** (either says AD → AD), else the non-control label |
# | `dx_detailed` | pinned to `AD` when reconciled `pheno` is AD; otherwise prefer the trio's clinical instrument, falling back to Diverse Cohorts when the trio is null |
# | `sex` | each genome takes **its own callset's** donor's sex — never reconciled |
#
# AD dominates because the AD label is neuropathological (§3) and a pathology call outranks
# a clinical one. `sex` is deliberately not reconciled: it is a property of the sequenced
# sample, and disagreement is evidence of a sample swap worth keeping as a flag.
#
# Conflicts are flagged, never dropped — `pheno_conflict`, `dx_conflict`, `sex_conflict`
# travel with the row so a sensitivity analysis can exclude them without rebuilding.

# In[16]:


DIVCO = "amp_ad_divco"
TRIO = {"amp_ad_msbb", "amp_ad_mayo", "amp_ad_rosmap"}


def reconcile_pheno(d, t):
    """Diverse Cohorts label vs trio label (each a value or "") -> (reconciled, conflict?)."""
    if d == t or not d or not t:      # agree, or only one of them is labelled
        return (d or t), False
    if "AD" in (d, t):                # AD-dominant: neuropathology outranks clinical
        return "AD", True
    non_control = {d, t} - {"control"}
    if "control" in (d, t) and len(non_control) == 1:
        return non_control.pop(), True
    return t, True                    # not observed in the data; deterministic + flagged


def reconcile_dx(dx_d, dx_t, pheno):
    """Pin to AD when pheno is AD, else prefer the trio, falling back to Diverse Cohorts."""
    conflict = bool(dx_d and dx_t and dx_d != dx_t)
    if pheno == "AD":
        return "AD", conflict
    return (dx_t or dx_d), conflict


def reconcile(rows):
    """Core rows for one donor -> (pheno, dx_detailed, pheno_conflict, dx_conflict)."""
    if len(rows) == 1:
        return rows[0]["pheno"], rows[0]["dx_detailed"], False, False
    dvc = next((r for r in rows if r["source_dataset"] == DIVCO), None)
    tri = next((r for r in rows if r["source_dataset"] in TRIO), None)
    if not (dvc and tri):             # >1 row but not the divco+trio pattern — not observed
        first = next((r for r in rows if r["pheno"]), rows[0])
        return first["pheno"], first["dx_detailed"], False, False
    pheno, pconf = reconcile_pheno(dvc["pheno"], tri["pheno"])
    dx, dconf = reconcile_dx(dvc["dx_detailed"], tri["dx_detailed"], pheno)
    return pheno, dx, pconf, dconf


# In[17]:


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
# | `contrasts.csv` | one row per written contrast: arm sizes, cohort composition, confound tag |
#
# **FID comes from the genotype fileset**, not the grain: plink2 matches on FID+IID and
# defaults a missing FID to `0`, so the AMP-PD callsets — which carry a non-zero FID — need
# their real one. `cohort_<ANC>_qc.fam` therefore has to exist alongside the PCs.
#
# A contrast arm may be restricted by source: `control@amppd` means AMP-PD callsets only
# (`wb_dwgs` or `br_dsnwgs`), `@ampad` the rest. Bare arms take any source. Both arms must
# reach `MIN_ARM` for a file to be written.
#
# **The confound tag** measures how far apart the two arms are in AMP-PD share. The primary
# contrast is confounded with callset by construction — AMP-AD supplies the AD cases, AMP-PD
# the PD cases — so `delta_amppd` is the number that says how much of a result could be
# cohort rather than disease. It is computed here, off `AMPPD_CALLSETS`, so adding a callset
# updates it in one place.
#
# **`EXCLUDE_DUAL`** drops the samples that resolved to both `divco_hs` and `wgs_harm`, for a
# WGS_Harm-only sensitivity run. It lives here rather than as a flag at GWAS time so that the
# sensitivity variant is a recorded artifact you can point at, not a switch someone remembers.

# In[18]:


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
            written.append({
                "ancestry": anc, "contrast": tag, "case_arm": case, "ctrl_arm": ctrl,
                "n_case": n_case, "n_ctrl": n_ctrl,
                "case_pct_amppd": round(case_pct, 1), "ctrl_pct_amppd": round(ctrl_pct, 1),
                "delta_amppd": round(delta, 1), "confound_tag": confound_tag(delta),
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


# ## 14. Appendix — how much of the AD arm is definition-dependent?
#
# **Read-only. Nothing below changes any table, and this section is not part of the build.**
#
# `pheno` defines AD neuropathologically. Every AD cohort also carries a *non-derived*
# diagnosis column we do not use. This asks one question: if we took those columns at face
# value instead, how different would the AD arm be?
#
# The four cohorts are not comparable on this axis:
#
# | Cohort | Alternative column | Complete | What it actually measures |
# |---|---|---|---|
# | ROSMAP | `dcfdx_lv` | 100% | **clinical** cognitive consensus — the only true clinical contrast |
# | Mayo | `diagnosis` | 98.5% | **neuropathological** consensus — pre-made vs our recomputation |
# | MSBB | `CDR` | 100% | dementia **severity**, not etiology (`diagnosis` is entirely NA) |
# | DivCo | `reag` | 79% | NIA-Reagan likelihood — also neuropathological |
#
# So only ROSMAP answers "clinical vs neuropathology." Mayo answers "someone else's
# pathology call vs ours." MSBB has no diagnosis at all — CDR says how impaired a donor was,
# not what caused it, so treat any AD count derived from it as an upper bound. DivCo's
# `ADoutcome` is already the call `pheno` uses, making `reag` a third pathology instrument.
#
# Why bother, given the definition is fixed: the completeness gap is real. Mayo's `amyThal`
# is ~44% complete against 98.5% for `diagnosis`, which is why our rule nulls 347 of 620
# Mayo donors. That is a power argument for a sensitivity arm, not a reason to redefine the
# primary phenotype.

# In[19]:


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


# In[20]:


# One explicit "AD by the alternative instrument" call per cohort. Each mapping is a
# judgement — stated here rather than buried, because the MSBB one in particular is weak.

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
# The rule returns null whenever its staging inputs are missing, and those donors fall out
# of both arms. Counting null *donors* overstates it — a donor with no genome contributes
# nothing to a GWAS either way. The cell below counts both, and the gap is the point.

# In[21]:


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
# Only the unclassifiable genomes are worth chasing. For each cohort this takes the
# genome-carrying donors with a null `pheno` and asks what its alternative column says
# about them: how many become classifiable, and into which arm.

# In[22]:


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
