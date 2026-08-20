#!/usr/bin/env python
# coding: utf-8
"""Shared between `clinical_core.py` and `analysis_grain.py`.

The clinical side runs at TWO points in the pipeline, and that is irreducible:

    clinical_core.py     §1-10, §12a    BEFORE step 1   — step 1 applies the sex files
    <steps 1-6 on the cluster>
    analysis_grain.py    §11-13         AFTER step 6    — the grain carries step 6's PCs

Unlike step 6's two passes — which came from one CSV carrying two unrelated things, and
collapsed once it was split — this round trip is real. The grain must carry PCs because step 7
reads them as covariates, and PCs require the genotype pipeline, which requires the sex files.
There is no column to split out.

What WAS avoidable is running the same script twice. `clinical_core.py` used to re-execute
§1-10 on its second invocation purely to rebuild `core` and `genomes` in memory — which also
silently rewrote the sex-update files that step 1 had already consumed. If any clinical input
had changed between the two runs, the files on disk would have stopped matching what genotools
actually applied, with nothing to flag it. §9's audit tables already hold exactly what the later
sections need, so `analysis_grain.py` reads those instead and nothing is re-derived.

This module holds only what both genuinely share: where things live, how files are read, and
the donor-level reconciliation rules.

GUARDRAIL: imported, never run. Defines no output and prints nothing.
"""
from pathlib import Path

import pandas as pd

# Every path below is anchored to this file, so both scripts run identically from any working
# directory and nothing outside this folder is ever read or written. No absolute paths, no
# environment variables, no cluster paths: to move the analysis, copy the folder.
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

# Clinical metadata — downloaded from Synapse / GCP (provenance at the top of clinical_core.py).
AD_META = DIR_WGS / "metadata"
DIVCO_META = DIR_DC / "metadata"
AMPPD_META = DATA / "amp-pd-genomics/metadata"
WB_META = DIR_WB / "metadata"
BR_META = DIR_BR / "metadata"

# Written by the cluster, rsynced back down into this folder between round trips.
# Absent on a fresh copy — analysis_grain.py's sections skip and say so.
MERGED = DATA / "merged"
EXCLUDE_REASONS = MERGED / "relatedness/exclude_reasons.tsv"
RETAINED_MANIFEST = MERGED / "relatedness/retained_manifest.csv"
PCA_DIR = MERGED / "by_ancestry_qc"

# §9's audit tables. Written by clinical_core.py, read by analysis_grain.py — this pair IS the
# handoff between the two halves, and it is why the second half re-derives nothing.
INDIVIDUAL_CORE = OUT / "individual_core.csv"
GENOME_CROSSWALK = OUT / "genome_crosswalk.csv"

# Vocabularies, shared so the two scripts cannot drift on what a valid label is.
SEX_VALUES = ["0", "1", "2"]
PHENO_VALUES = ["AD", "PD", "control", "other"]
DX_VALUES = ["PD", "AD", "MCI", "DLB", "PSP", "control", "other"]

# Which callsets are AMP-PD. Drives the @amppd / @ampad contrast tags in §13, and the
# cohort-confound percentages beside them. Shared because it is a fact about the callsets,
# not about either script: adding a callset should change both halves at once.
AMPPD_CALLSETS = {"wb_dwgs", "br_dsnwgs"}


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


# ── donor-level reconciliation ───────────────────────────────────────────────
#
# A donor can appear in two source studies — once in Diverse Cohorts and once in the
# ROSMAP/Mayo/MSBB trio — phenotyped by two different programs that do not always agree. The
# grain is per genome, so each genome needs one label.
#
# AD dominates because the AD label is neuropathological and a pathology call outranks a
# clinical one. `sex` is deliberately NOT reconciled here: it is a property of the sequenced
# sample, and disagreement is evidence of a sample swap worth keeping as a flag.
#
# These live in the shared module because BOTH halves apply them — §12a for sample_annot.csv
# and §12 for the grain — and the whole point of the split is that the two agree exactly.
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
    """§9's two audit tables -> (core, genomes). The handoff into analysis_grain.py.

    Returns them exactly as clinical_core.py held them in memory, minus `n_rules_hit`, which
    §9 drops on write and nothing downstream reads. Everything §11-13 need is here: core
    supplies (individual_id, source_dataset, sex, projid, pheno, dx_detailed) and genomes
    supplies (IID, source_callset, individual_id, source_dataset, rule, specimenID, tissue).
    """
    missing = [p for p in (INDIVIDUAL_CORE, GENOME_CROSSWALK) if not p.exists()]
    if missing:
        raise SystemExit(
            "missing §9 audit table(s): " + ", ".join(str(rel(p)) for p in missing)
            + "\nRun clinical_core.py first — it writes both, and needs no cluster output.")
    return rd(INDIVIDUAL_CORE), rd(GENOME_CROSSWALK)
