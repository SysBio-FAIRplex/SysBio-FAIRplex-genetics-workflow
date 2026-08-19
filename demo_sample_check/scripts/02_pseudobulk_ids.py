#!/usr/bin/env python3
"""02 — collapse the 179 pseudobulk sample lists to one donor-level table.

Two things make this less trivial than a `cat`:

1. The six cohorts name their donor column differently, so DONOR_COLS is tried in order
   and the first column present wins. It is an ordered priority list rather than a
   per-cohort mapping so a new cohort needs no code change.

2. The CMD files are CELL-level, not donor-level: `ID` is a 10x barcode and one file can
   carry 144k rows for 11 donors. Taking `ID` would produce nonsense. Every cohort is
   therefore deduplicated to distinct donors per file.

`projid` is carried alongside the donor where the cohort has it (the ROSMAP-space files),
because it is a second, independent key into the genotypes: a handful of AMP-AD
WGS_Harmonization samples are named MAP<projid> and match on nothing else.

Writes two files:
  pseudobulk_donors.tsv         long — one row per (cohort, tissue, celltype, donor)
  pseudobulk_donors_unique.txt  the deduplicated donor list, one ID per line

The long table is what makes step 03 able to say WHERE an overlapping donor came from,
which is the whole point of not flattening the directory tree in step 01.
"""
import csv
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
PSEUDOBULK = HERE / "pseudobulk"
OUT = HERE / "out"

# First match wins. `sample_id` is last because in the 3 CMD files that lack `donor_id`
# it is the donor, but elsewhere a donor can have several biosamples.
DONOR_COLS = ("donor_id", "participant_id", "individualID", "donor", "sample_id")


def celltype_of(path: Path) -> str:
    """CMD and rasle put the cell type in the directory; the rest in the filename."""
    stem = path.stem
    return path.parent.name if stem == "samples" else stem.removesuffix("_samples")


def norm_projid(raw: str | None) -> str:
    """projid is written as a float in rosmap ('3283241.0') and an int in mitrosmap."""
    raw = (raw or "").strip()
    if not raw:
        return ""
    try:
        return str(int(float(raw)))
    except ValueError:
        return raw


def main() -> int:
    files = sorted(PSEUDOBULK.rglob("*samples.[tc]sv"))
    if not files:
        sys.exit(f"no sample files under {PSEUDOBULK} — run 01_pull_pseudobulk.sh first")

    OUT.mkdir(exist_ok=True)
    rows, skipped = [], []

    for f in files:
        rel = f.relative_to(PSEUDOBULK)
        cohort = rel.parts[0].removesuffix("_pseudobulk_out")
        tissue = rel.parts[1]
        celltype = celltype_of(f)
        delim = "\t" if f.suffix == ".tsv" else ","

        with f.open(newline="") as fh:
            reader = csv.DictReader(fh, delimiter=delim)
            fields = reader.fieldnames or []
            col = next((c for c in DONOR_COLS if c in fields), None)
            if col is None:
                skipped.append((str(rel), fields))
                continue
            has_projid = "projid" in fields
            donors = set()
            for r in reader:
                d = (r.get(col) or "").strip()
                if not d:
                    continue
                donors.add((d, norm_projid(r.get("projid")) if has_projid else ""))

        for d, projid in sorted(donors):
            rows.append((cohort, tissue, celltype, d, projid, col))

    if skipped:
        print("WARNING — no donor column found, files skipped:", file=sys.stderr)
        for rel, cols in skipped:
            print(f"  {rel}  cols={cols}", file=sys.stderr)

    long_path = OUT / "pseudobulk_donors.tsv"
    with long_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cohort", "tissue", "celltype", "donor", "projid", "id_column"])
        w.writerows(rows)

    unique = sorted({r[3] for r in rows})
    (OUT / "pseudobulk_donors_unique.txt").write_text("\n".join(unique) + "\n")

    per_cohort = Counter(r[0] for r in rows)
    uniq_per_cohort = {c: len({r[3] for r in rows if r[0] == c}) for c in per_cohort}

    print(f"{len(files)} files -> {len(rows)} (cohort,tissue,celltype,donor) rows")
    print(f"{len(unique)} unique donors overall\n")
    print(f"{'cohort':<16} {'rows':>7} {'unique donors':>14}")
    for c in sorted(per_cohort):
        print(f"{c:<16} {per_cohort[c]:>7} {uniq_per_cohort[c]:>14}")
    print(f"\nwrote {long_path}")
    print(f"wrote {OUT / 'pseudobulk_donors_unique.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
