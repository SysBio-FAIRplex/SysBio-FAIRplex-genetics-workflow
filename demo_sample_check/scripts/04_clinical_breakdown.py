#!/usr/bin/env python3
"""04 — sex / phenotype / ancestry breakdown of the overlapping donors.

    python3 scripts/04_clinical_breakdown.py

Joins out/overlap.tsv to clinical_core_out/analysis_grain.csv on individual_id.

Two things the join has to get right:

1. The grain is POST-QC. A donor can match a callset in step 03 and still be absent here,
   because ancestry QC, relatedness or call-rate dropped it. That gap is the interesting
   number — it is the difference between "we have a genotype file with this donor in it"
   and "this donor is in the GWAS" — so it is reported rather than silently lost to an
   inner join.

2. 178 individual_ids carry more than one genotype IID (the same person sequenced in two
   callsets, e.g. divco_hs and wgs_harm). Rows are collapsed to one per donor so the
   breakdown counts people, not samples. Disagreement between a donor's rows on sex,
   pheno or ancestry is reported instead of being resolved by whichever sorted first.

Writes out/clinical_breakdown.tsv — one row per overlapping donor, with the pseudobulk
cohorts it appears in and its clinical fields.
"""
import csv
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parent
OUT = HERE / "out"

GRAIN = Path(os.environ.get("GRAIN", REPO / "clinical_core_out/analysis_grain.csv"))

SEX = {"1": "male", "2": "female", "0": "unknown", "": "unknown"}
FIELDS = ("sex", "pheno", "ancestry", "dx_detailed")


def table(title: str, rows: dict[str, Counter], order: list[str]) -> None:
    """cohort-by-value counts, printed as a matrix."""
    vals = sorted({v for c in rows.values() for v in c}, key=lambda x: (x == "", x))
    if not vals:
        return
    w = max([len(v) for v in vals] + [6]) + 2
    print(f"\n{title}")
    print(f"{'cohort':<16}{'n':>6}  " + "".join(f"{v or '(blank)':>{w}}" for v in vals))
    print("-" * (24 + w * len(vals)))
    for c in order:
        if c not in rows:
            continue
        tot = sum(rows[c].values())
        print(f"{c:<16}{tot:>6}  " + "".join(f"{rows[c].get(v, 0):>{w}}" for v in vals))


def main() -> int:
    ov_path = OUT / "overlap.tsv"
    if not ov_path.exists():
        sys.exit(f"{ov_path} missing — run 03_overlap.py first")
    if not GRAIN.exists():
        sys.exit(f"{GRAIN} missing — it is written by clinical_core.py §12")

    overlap = list(csv.DictReader(ov_path.open(), delimiter="\t"))
    # donor -> pseudobulk cohorts, and donor -> callsets it was found in
    donor_cohorts: dict[str, str] = {}
    donor_callsets: dict[str, set[str]] = defaultdict(set)
    for r in overlap:
        donor_cohorts[r["donor"]] = r["cohort"]
        donor_callsets[r["donor"]].add(r["callset"])
    donors = set(donor_cohorts)

    grain_rows: dict[str, list[dict]] = defaultdict(list)
    for r in csv.DictReader(GRAIN.open()):
        if r["individual_id"] in donors:
            grain_rows[r["individual_id"]].append(r)

    conflicts: list[tuple[str, str, list[str]]] = []
    merged: dict[str, dict] = {}
    for d, rows in grain_rows.items():
        rec = dict(rows[0])
        for f in FIELDS:
            seen = {r[f] for r in rows}
            if len(seen) > 1:
                conflicts.append((d, f, sorted(seen)))
        rec["_callsets"] = ",".join(sorted({r["source_callset"] for r in rows}))
        rec["_n_iids"] = len(rows)
        merged[d] = rec

    in_grain = set(merged)
    dropped = donors - in_grain

    print(f"overlapping donors (step 03)      : {len(donors)}")
    print(f"  present in analysis_grain        : {len(in_grain)}")
    print(f"  matched a callset but not in grain: {len(dropped)}   <- removed by QC")
    if conflicts:
        print(f"\n  ! {len(conflicts)} field conflicts across multi-IID donors:")
        for d, f, vals in conflicts[:10]:
            print(f"      {d}  {f}: {vals}")

    cohorts = sorted({donor_cohorts[d] for d in donors})
    order = cohorts + ["ALL"]

    def counts(field, transform=lambda x: x):
        out: dict[str, Counter] = defaultdict(Counter)
        for d, rec in merged.items():
            v = transform(rec[field])
            out[donor_cohorts[d]][v] += 1
            out["ALL"][v] += 1
        return out

    print("\n" + "=" * 60)
    print("Donors in the analysis grain (i.e. in the GWAS), by pseudobulk cohort")
    print("=" * 60)
    table("sex", counts("sex", lambda v: SEX.get(v, v)), order)
    table("pheno", counts("pheno"), order)
    table("ancestry", counts("ancestry"), order)
    table("dx_detailed", counts("dx_detailed"), order)

    drop_by_cohort = Counter(donor_cohorts[d] for d in dropped)
    if dropped:
        print("\nQC-dropped overlapping donors, by cohort")
        for c in cohorts:
            if drop_by_cohort[c]:
                print(f"  {c:<16}{drop_by_cohort[c]:>5}")

    out_path = OUT / "clinical_breakdown.tsv"
    cols = ["donor", "pseudobulk_cohort", "callsets_matched", "in_grain",
            "grain_callset", "n_genotype_iids", *FIELDS]
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for d in sorted(donors):
            rec = merged.get(d)
            w.writerow([
                d, donor_cohorts[d], ",".join(sorted(donor_callsets[d])),
                int(d in in_grain),
                rec["_callsets"] if rec else "",
                rec["_n_iids"] if rec else 0,
                *[(rec[f] if rec else "") for f in FIELDS],
            ])
    print(f"\nwrote {out_path}  ({len(donors)} donors)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
