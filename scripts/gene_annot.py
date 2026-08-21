#!/usr/bin/env python3
"""Gene coordinates from a UCSC refFlat table — THE single source of locus coordinates.

WHY THIS EXISTS. Locus coordinates used to be hardcoded, as a `SENTINEL_LOCI` list of
`(name, chrom, pos, window)` tuples duplicated verbatim across the scripts that needed them,
including af_concordance_build.py. Independent copies of hand-typed coordinates is a bug waiting to
happen, and it had already happened. Measured against `ref/refFlat.txt` (GRCh38), the hardcoded
windows missed real gene territory:

    CR1     91,765 bp of the gene outside the window   (63% of a 145 kb gene)
    LRRK2   69,284 bp
    SNCA    38,304 bp
    BIN1     1,977 bp
    GBA1    the symbol does not exist in refFlat at all — it is `GBA` there

A sentinel that silently covers 37% of its gene, or that never matches because the symbol is wrong,
reports "no known locus flagged" and looks like a clean result. That is the failure mode this
replaces.

BUILD. refFlat carries no build marker, so `load_genes` verifies GRCh38 by checking APOE against its
known GRCh38 extent (chr19:44.90-44.91 Mb; hg19 would put it at ~45.41 Mb) and raises otherwise. A
silently-hg19 annotation would mislabel every variant in the study by ~500 kb at this locus.

CONVENTIONS. refFlat `txStart` is 0-based half-open, `txEnd` is 1-based inclusive, so a 1-based
variant position `p` lies in a transcript when `txStart < p <= txEnd`. Gene extents are collapsed
across all transcripts of a symbol (min start, max end). Non-primary contigs (`chr*_alt`,
`chr*_random`, `chrUn_*`) are dropped — 10,456 of 88,819 rows.

USAGE
    from gene_annot import load_genes, gene_region, Annotator
    genes = load_genes()                                  # {symbol: (chrom, start1, end)}
    reg   = gene_region("CR1", flank=50_000, genes=genes)  # ('1', 207446157, 207691765)
    ann   = Annotator(genes)
    ann.genes_at("19", 44_888_997)                         # ['NECTIN2']
    ann.nearest("19", 44_890_000)                          # ('NECTIN2', 777) — signed distance 0 if inside

    # CLI
    python3 scripts/gene_annot.py CR1 SNCA LRRK2 GBA1
    python3 scripts/gene_annot.py --at chr19:44888997 --at 6:32000000
    python3 scripts/gene_annot.py --at chr12:40227079:C:T     # full variant IDs work too

THIS CLI IS THE ANSWER TO "what gene is this variant in". It replaced a 9-gene sentinel list that
two pipeline scripts used to check their flagged variants against (deleted 2026-08-20 — see the
comment where SENTINEL_GENES used to live). Per-variant and unarbitrary beats a hand-picked set:
asked this way, the 21 variants once reported as "LRRK2±500kb" turn out to be in SLC2A13 and
C12orf40, 240-435 kb from LRRK2.

GUARDRAIL: reads a public gene annotation only. No genotypes, no sample data.
"""
from bisect import bisect_right
import os
import sys

# refFlat uses older HGNC symbols for some genes. Extend as needed; the CLI warns on any name that
# fails to resolve rather than returning nothing silently.
ALIASES = {
    "GBA1": "GBA",        # HGNC renamed GBA -> GBA1 in 2022; refFlat predates it
    "PVRL2": "NECTIN2",   # older symbol for NECTIN2
    "PVR": "NECTIN2",     # occasional mislabel in older locus tables
}

# GRCh38 sanity anchor. (chrom, lo, hi) that APOE's collapsed extent must fall inside.
_BUILD_ANCHOR = ("19", 44_900_000, 44_915_000)

_PRIMARY = {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}


def refflat_path(explicit=None):
    """Resolve the refFlat table: explicit arg, then $REFFLAT, $REF_DIR/refFlat.txt, <bundle>/ref/."""
    cands = []
    if explicit:
        cands.append(explicit)
    if os.environ.get("REFFLAT"):
        cands.append(os.environ["REFFLAT"])
    if os.environ.get("REF_DIR"):
        cands.append(os.path.join(os.environ["REF_DIR"], "refFlat.txt"))
    bundle = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cands.append(os.path.join(bundle, "ref", "refFlat.txt"))
    for c in cands:
        if c and os.path.exists(c):
            return c
    raise FileNotFoundError(
        "refFlat.txt not found. Looked in: " + ", ".join(str(c) for c in cands) +
        "\nSet $REFFLAT or place it at <bundle>/ref/refFlat.txt (and copy it to $REF_DIR on the cluster)."
    )


def load_genes(path=None, verify_build=True):
    """-> {symbol: (chrom_no_chr_prefix, start_1based, end)}, collapsed across transcripts.

    Primary contigs only. A symbol occurring on more than one primary contig (92 of them, mostly
    paralogous families) keeps its longest span; `load_genes.multi` lists them after the call so a
    caller can check whether a locus it cares about is ambiguous.
    """
    path = refflat_path(path)
    spans = {}      # symbol -> {chrom: [start, end]}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            chrom = f[2]
            if "_" in chrom:
                continue                      # alt / random / unplaced
            c = chrom[3:] if chrom.startswith("chr") else chrom
            if c not in _PRIMARY:
                continue
            try:
                s, e = int(f[4]) + 1, int(f[5])   # txStart is 0-based half-open
            except ValueError:
                continue
            g = f[0]
            per = spans.setdefault(g, {})
            if c in per:
                per[c][0] = min(per[c][0], s)
                per[c][1] = max(per[c][1], e)
            else:
                per[c] = [s, e]

    genes, multi = {}, []
    for g, per in spans.items():
        if len(per) > 1:
            multi.append(g)
        c = max(per, key=lambda k: per[k][1] - per[k][0])   # longest span wins
        genes[g] = (c, per[c][0], per[c][1])
    load_genes.multi = sorted(multi)
    load_genes.path = path

    if verify_build:
        if "APOE" not in genes:
            raise ValueError(f"{path}: no APOE row — cannot verify this is a GRCh38 refFlat table")
        c, s, e = genes["APOE"]
        ac, lo, hi = _BUILD_ANCHOR
        if not (c == ac and lo <= s <= hi and lo <= e <= hi):
            raise ValueError(
                f"{path}: APOE at {c}:{s}-{e}, expected GRCh38 {ac}:{lo}-{hi}. "
                "This looks like the wrong build (hg19 puts APOE near 45,409,039). "
                "Pass verify_build=False only if you are certain."
            )
    return genes


def gene_region(name, flank=0, genes=None):
    """-> (chrom, start, end) for a symbol, widened by `flank` bp. Applies ALIASES. None if absent."""
    if genes is None:
        genes = load_genes()
    key = name if name in genes else ALIASES.get(name, name)
    if key not in genes:
        return None
    c, s, e = genes[key]
    return (c, max(1, s - flank), e + flank)


def gene_regions(names, flank=50_000, genes=None, warn=True):
    """-> [(display_name, chrom, start, end)] for a list of symbols, skipping unresolved ones.

    Unresolved names are reported to stderr rather than dropped in silence — that silence is
    precisely how the old GBA1 entry went unnoticed.

    Named `sentinel_regions` until 2026-08-20. Nothing about it was ever sentinel-specific — it
    resolves whatever symbols it is handed — and the name outlived the sentinel set it was named
    for, which is how a deleted concept quietly stays alive in a codebase.
    """
    if genes is None:
        genes = load_genes()
    out, missing = [], []
    for n in names:
        r = gene_region(n, flank, genes)
        if r is None:
            missing.append(n)
            continue
        out.append((n, r[0], r[1], r[2]))
    if missing and warn:
        print(f"WARNING: no refFlat entry for: {', '.join(missing)} — add an alias in "
              f"gene_annot.ALIASES or fix the symbol", file=sys.stderr)
    return out


# ── NO SENTINEL SET LIVES HERE ANY MORE ───────────────────────────────────────────────────────────
# `SENTINEL_GENES` (a 9-gene AD/PD list), `SENTINEL_FLANK` (500 kb), `SENTINEL_VARIANTS` (the two
# APOE causal positions) and `sentinel_hits()` were deleted 2026-08-20, together with both callers —
# af_concordance_build.py's tripwire and 08_ctrl_ctrl_filter.py's known-loci report. Reasons are in
# those two files and in PROJECT_LOG.md; the short version is that an uncited 9-gene list makes a
# NEGATIVE result printable ("no known locus flagged") that no evidence supports, and the 500 kb
# flanks attach a gene's name to variants in entirely different genes 240-435 kb away.
#
# Do not reintroduce it. The honest form of the same question is per-variant and already available:
# `Annotator.label()` / the `--at` CLI name the gene a position actually falls in, for every variant
# rather than for nine.
#
# Two measurements from the deleted code, kept because they cost a run each and would otherwise have
# to be re-derived: MAPT±500kb caught 2,286 of EUR AD-vs-PSP's 2,318 genome-wide hits, the other 32
# (1.4%) sitting beyond 46,528,333 in the tail of the 17q21.31 inversion LD block — so 500 kb was a
# defensible LD-block proxy for "did my signal land in a known locus", and a poor one for "is this
# variant part of a known locus". Sizing windows from r² decay against an external LD panel is the
# principled fix and was never done. And a uniform flank was chosen over per-locus tuning on purpose:
# tuning from our own hit spans sizes the check from the answer it is meant to protect.


def parse_variant_id(vid):
    """'chr19:44908684:T:C' -> ('19', 44908684). None if unparseable."""
    p = str(vid).split(":")
    if len(p) < 2:
        return None
    c = p[0][3:] if p[0].startswith("chr") else p[0]
    try:
        return (c, int(p[1]))
    except ValueError:
        return None


class Annotator:
    """Position -> gene lookup. Per-chromosome intervals sorted by start, with a running max-end so
    overlapping genes are handled correctly (a plain bisect on starts would miss a long gene that
    began well before the query)."""

    def __init__(self, genes=None):
        if genes is None:
            genes = load_genes()
        self.by_chrom = {}
        for g, (c, s, e) in genes.items():
            self.by_chrom.setdefault(c, []).append((s, e, g))
        for c, iv in self.by_chrom.items():
            iv.sort()
            starts = [x[0] for x in iv]
            maxend, run = [], 0
            for _, e, _g in iv:
                run = max(run, e)
                maxend.append(run)
            self.by_chrom[c] = (iv, starts, maxend)

    @staticmethod
    def _norm(chrom):
        c = str(chrom)
        return c[3:] if c.lower().startswith("chr") else c

    def genes_at(self, chrom, pos):
        """-> sorted list of symbols whose collapsed extent contains `pos` (1-based). [] if none."""
        c = self._norm(chrom)
        rec = self.by_chrom.get(c)
        if not rec:
            return []
        iv, starts, maxend = rec
        i = bisect_right(starts, pos)          # all intervals with start <= pos
        hits = []
        for j in range(i - 1, -1, -1):
            if maxend[j] < pos:                # no interval at or before j can still reach pos
                break
            s, e, g = iv[j]
            if s <= pos <= e:
                hits.append(g)
        return sorted(set(hits))

    def nearest(self, chrom, pos):
        """-> (symbol, distance) for the closest gene; distance 0 when inside. (None, None) if the
        chromosome is absent. Ties broken alphabetically for determinism."""
        inside = self.genes_at(chrom, pos)
        if inside:
            return (inside[0], 0)
        c = self._norm(chrom)
        rec = self.by_chrom.get(c)
        if not rec:
            return (None, None)
        iv, starts, _ = rec
        best, bestd = None, None
        i = bisect_right(starts, pos)
        for j in range(max(0, i - 50), min(len(iv), i + 50)):   # local scan is enough once sorted
            s, e, g = iv[j]
            d = 0 if s <= pos <= e else (s - pos if pos < s else pos - e)
            if bestd is None or d < bestd or (d == bestd and g < best):
                best, bestd = g, d
        return (best, bestd)

    def label(self, chrom, pos, max_dist=100_000):
        """-> a display string: 'GENE' inside, 'GENE(~1.2kb)' nearby, 'intergenic' beyond max_dist.

        The distance is unsigned on purpose. Strand-aware upstream/downstream would need the gene's
        strand and would invite reading a regulatory relationship into what is only proximity."""
        g, d = self.nearest(chrom, pos)
        if g is None or d is None or d > max_dist:
            return "intergenic"
        return g if d == 0 else f"{g}(~{d/1000:.1f}kb)"


def _main(argv):
    args = argv[1:]
    if not args:
        print(__doc__)
        return 0
    ats, names, flank = [], [], 0
    i = 0
    while i < len(args):
        if args[i] == "--at":
            ats.append(args[i + 1]); i += 2
        elif args[i] == "--flank":
            flank = int(args[i + 1]); i += 2
        else:
            names.append(args[i]); i += 1

    genes = load_genes()
    print(f"refFlat: {load_genes.path}", file=sys.stderr)
    print(f"  {len(genes):,} symbols on primary contigs; build verified GRCh38 via APOE "
          f"({genes['APOE'][0]}:{genes['APOE'][1]:,}-{genes['APOE'][2]:,})", file=sys.stderr)
    if load_genes.multi:
        print(f"  {len(load_genes.multi)} symbols span >1 contig (longest span kept)", file=sys.stderr)

    if names:
        for n, c, s, e in gene_regions(names, flank, genes):
            via = "" if n in genes else f"  [via alias -> {ALIASES.get(n)}]"
            print(f"{n:12} {c}:{s:,}-{e:,}   {e - s + 1:,} bp{via}")
    if ats:
        ann = Annotator(genes)
        for a in ats:
            # parse_variant_id, not a second inline parser: it accepts a bare 'chr:pos' and a full
            # 'chr19:44908684:T:C' variant ID alike, so the IDs in an exclusion list or a .assoc can
            # be pasted straight in. It was the only remaining caller-facing use after the sentinel
            # deletion, and a duplicate parse here is what would have made it dead code.
            cp = parse_variant_id(a)
            if cp is None:
                print(f"{a}: cannot parse as chr:pos or a chr:pos:ref:alt variant ID",
                      file=sys.stderr)
                continue
            c, p = cp
            g = ann.genes_at(c, p)
            n, d = ann.nearest(c, p)
            print(f"chr{c}:{p:,}  in={','.join(g) if g else '-'}  nearest={n} ({d:,} bp)  "
                  f"label={ann.label(c, p)}")
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv))
