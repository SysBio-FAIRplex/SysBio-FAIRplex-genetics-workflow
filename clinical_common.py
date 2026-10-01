#!/usr/bin/env python
# coding: utf-8
"""Paths, readers and donor-level reconciliation shared by clinical_core.py and analysis_grain.py.

The clinical side runs twice, and that is irreducible: clinical_core.py (§1-10, §12a) before
step 1, which applies the sex files; analysis_grain.py (§11-13) after step 6, because the grain
carries step 6's PCs. The second half reads §9's audit tables and re-derives nothing.

GUARDRAIL: imported, never run. Defines no output and prints nothing.
"""
from pathlib import Path

import pandas as pd

# Every path is anchored to this file; nothing outside the project root is read or written.
PROJECT_ROOT = Path(__file__).resolve().parent

DATA = PROJECT_ROOT / "data"
OUT = PROJECT_ROOT / "clinical_core_out"

# Per-callset roots, named as config.sh's DIR_*. AMP-AD keeps pgens in pgen/, WB-DWGS in joint_calls/.
DIR_WGS = DATA / "amp-ad-genomics/WGS_Harmonization"
DIR_DC = DATA / "amp-ad-genomics/DivCo_HS"
DIR_WB = DATA / "amp-pd-genomics/WB-DWGS"
DIR_BR = DATA / "amp-pd-genomics/BR-DSNWGS"

# Clinical metadata — pulled from Synapse / GCP on helix into data/ (provenance: wgs_core.ipynb §0).
AD_META = DIR_WGS / "metadata"
DIVCO_META = DIR_DC / "metadata"
AMPPD_META = DATA / "amp-pd-genomics/metadata"
WB_META = DIR_WB / "metadata"
BR_META = DIR_BR / "metadata"

# Written by genetics steps 1-6; analysis_grain.py's sections skip, and say so, until they exist.
MERGED = DATA / "merged"
EXCLUDE_REASONS = MERGED / "relatedness/exclude_reasons.tsv"
RETAINED_MANIFEST = MERGED / "relatedness/retained_manifest.csv"
PCA_DIR = MERGED / "by_ancestry_qc"

# §9's audit tables: the handoff from clinical_core.py to analysis_grain.py.
INDIVIDUAL_CORE = OUT / "individual_core.csv"
GENOME_CROSSWALK = OUT / "genome_crosswalk.csv"

# Vocabularies, shared so the two scripts cannot drift on what a valid label is.
SEX_VALUES = ["0", "1", "2"]
PHENO_VALUES = ["AD", "PD", "control", "other"]
DX_VALUES = ["PD", "AD", "MCI", "DLB", "PSP", "control", "other"]

# AMP-PD callsets. Drives §13's @amppd/@ampad arms and the confound percentages beside them.
AMPPD_CALLSETS = {"wb_dwgs", "br_dsnwgs"}


def rel(path):
    """A path relative to PROJECT_ROOT, for logs."""
    return Path(path).relative_to(PROJECT_ROOT)


def rd(path, sep=","):
    """Read everything as stripped strings. Missing stays "" (never NaN) so nulls
    show up in value_counts and crosstabs instead of being silently dropped."""
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    return df.apply(lambda s: s.str.strip())


# ── donor-level reconciliation ───────────────────────────────────────────────
# A donor in both Diverse Cohorts and the ROSMAP/Mayo/MSBB trio can carry two labels; each genome
# needs one. AD dominates: neuropathology outranks a clinical call. Sex is never reconciled — a
# disagreement is evidence of a sample swap. Shared so §12 and §12a apply identical rules.
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


def load_audit_tables():
    """§9's two audit tables -> (core, genomes), as clinical_core.py wrote them."""
    missing = [p for p in (INDIVIDUAL_CORE, GENOME_CROSSWALK) if not p.exists()]
    if missing:
        raise SystemExit(
            "missing §9 audit table(s): " + ", ".join(str(rel(p)) for p in missing)
            + "\nRun clinical_core.py first — it writes both, and needs no genotype-pipeline output.")
    return rd(INDIVIDUAL_CORE), rd(GENOME_CROSSWALK)
