#!/usr/bin/env python3
"""03 — overlap between the pseudobulk donors and one or more genotype sample lists.

    python3 scripts/03_overlap.py                       # every callset in CALLSETS
    python3 scripts/03_overlap.py name=/path/to.psam    # or an explicit list

Takes .psam / .fam / plain one-ID-per-line. Every callset is reported side by side, so
re-running this against the finished GWAS (plan step 6) means adding a line to CALLSETS,
not editing any logic.

## Why matching is not a set intersection

A genotype IID is not a person. The four callsets were assembled by different groups and
name samples after the specimen, the assay, or a sequencing manifest entry, so the same
donor is `MAP<8d>` in AMP-AD WGS_Harmonization, `R<7d>` in every ROSMAP-derived
pseudobulk, and `<individualID>_DLPFC_WGS` in DivCo_HS. The pseudobulk files, by contrast, are
already donor-level.

The bridge is `clinical_core_out/genome_crosswalk.csv`, which clinical_core.py builds and
the rest of this project already relies on: one row per genotype IID, giving the
`individual_id` it belongs to and the `rule` that established it. This script does NOT
re-derive that mapping — an earlier version of it did, with hand-written regexes, and
recovered 317 donors where the crosswalk finds 716. Most of the difference is ROSMAP,
where `MAP<n>` -> `R<m>` comes from a sequencing manifest and is not obtainable by string
surgery on the ID at all.

So: genotype IID --(crosswalk)--> individual_id <--(equality)-- pseudobulk donor.

`projid` is retained as a secondary key for the ROSMAP-space cohorts. It is expected to
add nothing now that the crosswalk is in play; it is reported separately so that if it
ever does fire, that is visible rather than silently folded into the total.
"""
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parent
OUT = HERE / "out"

CROSSWALK = Path(
    __import__("os").environ.get("CROSSWALK", REPO / "clinical_core_out/genome_crosswalk.csv")
)

# name -> path. Add a line here to check a new callset; nothing else changes.
CALLSETS = {
    # The demo's two arms. Both are EUR-only, so neither can reach the AFR donors that the
    # current callsets do — these columns are a ceiling within EUR, not overall.
    "demo_AMP_PD_EUR": HERE / "ad-pd-gwas-demo/FILTERED.AMP_PD_EUR.psam",
    "demo_AMP_AD_EUR": HERE / "ad-pd-gwas-demo/FILTERED.AMP_AD_EUR.psam",
    "AMP-PD_WB-DWGS": REPO / "data/amp-pd-genomics/WB-DWGS/joint_calls/all_chrs_merged.psam",
    "AMP-PD_BR-DSNWGS": REPO / "data/amp-pd-genomics/BR-DSNWGS/pgen/br_dsnwgs_hg38.psam",
    "AMP-AD_WGS_Harm": REPO / "data/amp-ad-genomics/WGS_Harmonization/pgen/wgs_harm_hg38.psam",
    "AMP-AD_DivCo_HS": REPO / "data/amp-ad-genomics/DivCo_HS/pgen/divco_hs_hg38.psam",
}


def read_ids(path: Path) -> list[str]:
    """IID column of a .psam/.fam, else the first field of a plain list."""
    lines = path.read_text().splitlines()
    if not lines:
        return []
    head = lines[0].lstrip("#").split("\t")
    if "IID" in head:  # .psam
        i, body = head.index("IID"), lines[1:]
    elif path.suffix == ".fam":
        i, body = 1, lines
    else:
        i, body = 0, lines
    out = []
    for ln in body:
        if not ln.strip():
            continue
        parts = ln.split("\t") if "\t" in ln else ln.split()
        if i < len(parts):
            out.append(parts[i].strip())
    return out


def main() -> int:
    long_path = OUT / "pseudobulk_donors.tsv"
    if not long_path.exists():
        sys.exit(f"{long_path} missing — run 02_pseudobulk_ids.py first")
    pb = list(csv.DictReader(long_path.open(), delimiter="\t"))

    if not CROSSWALK.exists():
        sys.exit(
            f"{CROSSWALK} missing — it is built by clinical_core.py and is what makes the\n"
            "genotype IIDs comparable to donor IDs. Set CROSSWALK= to point elsewhere."
        )
    cw = list(csv.DictReader(CROSSWALK.open()))
    iid_to_ind = {r["IID"]: (r["individual_id"], r["rule"]) for r in cw}

    args = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
    callsets = {k: Path(v) for k, v in args.items()} if args else CALLSETS

    # donor -> candidate keys; a donor may sit in several cohorts (13 do), so cohort
    # membership is a set per cohort and the per-cohort rows below can sum to more than
    # the unique total. That is correct, not double counting.
    donor_keys: dict[str, set[str]] = defaultdict(set)
    cohort_donors: dict[str, set[str]] = defaultdict(set)
    projid_of: dict[str, str] = {}
    for r in pb:
        d = r["donor"]
        donor_keys[d].add(d)
        if r.get("projid"):
            donor_keys[d].add(r["projid"])
            projid_of[d] = r["projid"]
        cohort_donors[r["cohort"]].add(d)
    cohorts = sorted(cohort_donors)

    hits: dict[str, dict[str, tuple[str, str]]] = {}  # callset -> donor -> (geno IID, rule)
    sizes: dict[str, int] = {}
    unmapped: dict[str, int] = {}
    for name, path in callsets.items():
        if not path.exists():
            print(f"  ! skipping {name}: {path} not found", file=sys.stderr)
            continue
        gids = read_ids(path)
        sizes[name] = len(gids)
        index: dict[str, tuple[str, str]] = {}
        miss = 0
        for g in gids:
            if g in iid_to_ind:
                ind, rule = iid_to_ind[g]
                index.setdefault(ind, (g, rule))
            else:
                miss += 1
                index.setdefault(g, (g, "not-in-crosswalk"))
        unmapped[name] = miss
        hits[name] = {
            d: index[k] for d, keys in donor_keys.items() for k in keys if k in index
        }

    names = list(hits)
    w = max((len(n) for n in names), default=10) + 2
    hdr = f"{'cohort':<16}{'donors':>7}  " + "".join(f"{n:>{w}}" for n in names)
    print(hdr)
    print("-" * len(hdr))
    for c in cohorts:
        ds = cohort_donors[c]
        print(f"{c:<16}{len(ds):>7}  " + "".join(f"{len(ds & set(hits[n])):>{w}}" for n in names))
    print("-" * len(hdr))
    alld = set(donor_keys)
    print(f"{'ALL unique':<16}{len(alld):>7}  " + "".join(f"{len(hits[n]):>{w}}" for n in names))
    print(f"{'callset size':<16}{'':>7}  " + "".join(f"{sizes[n]:>{w}}" for n in names))
    print(f"{'not in xwalk':<16}{'':>7}  " + "".join(f"{unmapped[n]:>{w}}" for n in names))

    matched = set().union(*(set(h) for h in hits.values())) if hits else set()
    print(f"\nmatched by at least one callset: {len(matched)} / {len(alld)} donors")
    for name in names:
        rules = Counter(r for _, r in hits[name].values())
        if rules:
            print(f"  {name:<20} " + ", ".join(f"{k}={v}" for k, v in rules.most_common()))

    # did projid ever succeed where the donor ID itself did not?
    via_projid = sorted(
        d
        for name in names
        for d in hits[name]
        if d in projid_of and hits[name][d][0] not in (d,) and projid_of[d] == hits[name][d][0]
    )
    print(f"matched only via projid: {len(via_projid)}")

    rows_out = [
        (donor_cohorts, d, name, hits[name][d][0], hits[name][d][1])
        for name in names
        for d in sorted(hits[name])
        for donor_cohorts in [",".join(sorted(c for c in cohorts if d in cohort_donors[c]))]
    ]
    ov = OUT / "overlap.tsv"
    with ov.open("w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(["cohort", "donor", "callset", "genotype_iid", "crosswalk_rule"])
        wr.writerows(rows_out)
    unmatched = sorted(alld - matched)
    (OUT / "unmatched_donors.txt").write_text("\n".join(unmatched) + "\n")
    print(f"\nwrote {ov}  ({len(rows_out)} rows)")
    print(f"wrote {OUT / 'unmatched_donors.txt'}  ({len(unmatched)} donors)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
