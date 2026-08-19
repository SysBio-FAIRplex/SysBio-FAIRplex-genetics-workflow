#!/usr/bin/env python3
"""
DIAGNOSTIC — is a callset's variant ORDER the same as the reference panel's?  READ-ONLY.

    python3 scripts/diag_order.py <ref_panel.bim> <callset.pvar|.bim> [more callsets...]

WHY THIS EXISTS
---------------
genotools/ancestry.py:246-247 builds the study genotype matrix and then does:

    raw_geno.columns = col_names        # col_names is the REFERENCE panel's SNP order

That is a POSITIONAL rename. Column i of the study matrix is asserted to be the same variant
as column i of the reference matrix, and nothing checks it. The reorder that would enforce it,

    geno_snps = geno_snps[ref_snps.columns]      # ancestry.py:244

sits inside `if not self.train`, so it does not run when a model is being trained from a
--ref_panel / --ref_labels pair. If the two orders differ, every column is silently mislabeled:
each variant is standardized by another variant's mean/SD and projected through another
variant's loading. Counts still match, no exception is raised, and the run exits 0.

The projected PCs then show a fixed offset (same permutation for every sample) plus isotropic
noise with NO eigenvalue decay across PCs — which is what BR-DSNWGS's projection looks like:
study sd 1.8-3.9 flat across PC1-PC10 against a reference decaying 103.9 -> 18.1.

WHAT IT REPORTS
---------------
  1. chromosome order of appearance   — lexicographic (1,10,11,...,2,20) vs numeric (1,2,3,...)
                                        is the classic per-chromosome-concatenation artifact
  2. within-chromosome position sort   — whether positions ascend inside each chromosome
  3. relative order of the SHARED variants — the one that actually decides it: walk the
     callset in file order, keep variants the panel also has, and check their panel ranks
     ascend. Inversions here are the bug; zero inversions means the positional rename is safe.

Run it on the callset that WORKS as a control. wgs_harm cleared this step in July under the
same code, so its order must agree with the panel; if BR-DSNWGS's does not, that is the fault.
"""

import sys
import os
import re


def norm_chrom(c):
    """chr1 / 1 / CHR1 -> '1'.  Keeps X/Y/MT as-is, uppercased."""
    c = str(c).strip()
    c = re.sub(r'^chr', '', c, flags=re.IGNORECASE)
    return c.upper()


def read_variants(path):
    """Yield (chrom, pos) in FILE ORDER. Handles .bim (CHR ID CM BP ...) and .pvar (#CHROM POS ID ...)."""
    is_pvar = path.endswith('.pvar')
    with open(path) as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                # .pvar header line — confirm column layout is the expected CHROM POS
                continue
            f = line.split(None, 5)
            if len(f) < 4:
                continue
            if is_pvar:
                yield norm_chrom(f[0]), f[1]
            else:
                yield norm_chrom(f[0]), f[3]


def numeric_key(c):
    return (0, int(c)) if c.isdigit() else (1, 0)


def describe(path, ref_rank=None):
    """Single streaming pass. Callset .pvar files run to 16.8M variants, so nothing is
    materialized except the panel's own rank dict (~209k) and the shared-variant ranks."""
    label = os.path.basename(path)
    print(f'\n--- {label} ---')
    if not os.path.isfile(path):
        print('    MISSING')
        return {} if ref_rank is None else None

    building_ref = ref_rank is None
    rank = {} if building_ref else None
    ranks = [] if not building_ref else None

    total = 0
    order, seen = [], set()
    unsorted_chroms = []
    prev_c, prev_p = None, None

    for c, p_raw in read_variants(path):
        p = int(p_raw)

        if c not in seen:
            seen.add(c)
            order.append(c)

        if c != prev_c:
            prev_c, prev_p = c, p
        else:
            if p < prev_p and c not in unsorted_chroms:
                unsorted_chroms.append(c)
            prev_p = p

        if building_ref:
            rank.setdefault((c, p_raw), total)
        elif (c, p_raw) in ref_rank:
            ranks.append(ref_rank[(c, p_raw)])

        total += 1

    print(f'    variants                : {total:,}')
    shown = ' '.join(order[:26]) + (' ...' if len(order) > 26 else '')
    print(f'    chromosome order        : {shown}')

    autosomes = [c for c in order if c.isdigit()]
    lexicographic = autosomes == sorted(autosomes, key=str)
    numeric = autosomes == sorted(autosomes, key=numeric_key)
    if numeric and not lexicographic:
        print('    chrom ordering          : NUMERIC (1,2,3,...)')
    elif lexicographic and not numeric:
        print('    chrom ordering          : LEXICOGRAPHIC (1,10,11,...) <-- per-chrom concat artifact')
    elif numeric and lexicographic:
        print('    chrom ordering          : both (too few autosomes to distinguish)')
    else:
        print('    chrom ordering          : NEITHER — interleaved or unsorted')

    print(f'    positions ascend in-chrom: {"yes" if not unsorted_chroms else "NO -> " + ",".join(unsorted_chroms[:10])}')

    if building_ref:
        return rank

    # THE decisive check: do the shared variants appear in the same relative order?
    print(f'    shared with panel       : {len(ranks):,} of {len(ref_rank):,} panel variants')
    if not ranks:
        print('    relative order          : n/a — no overlap')
        return

    # Count rank DROPS, not a percentage. A chromosome-block permutation rearranges every
    # column while producing only a handful of adjacent drops — lexicographic order
    # (1,10,...,19,2,20,21,22,3,...) drops exactly twice, at 19->2 and at 22->3. Reporting
    # "2 of 209,067 pairs = 0.0%" reads as harmless and is the opposite of the truth, so
    # report the drops and the misplaced-column count instead.
    drops = [i for i, (a, b) in enumerate(zip(ranks, ranks[1:])) if b < a]
    print(f'    rank drops              : {len(drops):,}')
    if not drops:
        print('    VERDICT                 : SAME ORDER — the positional rename at')
        print('                              ancestry.py:247 is safe for this callset.')
        return

    misplaced = sum(1 for i, r in enumerate(ranks) if r != i)
    print(f'    columns out of position : {misplaced:,} of {len(ranks):,} '
          f'({100.0 * misplaced / len(ranks):.1f}%)')
    print('    VERDICT                 : ORDER DIFFERS. ancestry.py:247 renames columns')
    print('                              POSITIONALLY, so each variant is standardized by')
    print('                              another variant\'s mean/SD and projected through')
    print('                              another variant\'s loading. This is the bug.')


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    ref_path, callsets = sys.argv[1], sys.argv[2:]

    print('=' * 62)
    print('REFERENCE PANEL (defines the column order genotools projects onto)')
    print('=' * 62)
    ref_rank = describe(ref_path)

    print()
    print('=' * 62)
    print('CALLSETS')
    print('=' * 62)
    for path in callsets:
        describe(path, ref_rank)

    print()
    print('=' * 62)
    print('DONE — read-only, nothing was modified')
    print('=' * 62)


if __name__ == '__main__':
    main()
