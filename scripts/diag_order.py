#!/usr/bin/env python3
"""
DIAGNOSTIC — is a callset's variant ORDER the same as the reference panel's?  READ-ONLY.

    python3 scripts/diag_order.py <ref_panel.bim> <callset.pvar|.bim> [more callsets...]

Mandatory before step 1 on any callset not already cleared. GenoTools (ancestry.py:247) renames
the study matrix's columns to the panel's SNP order BY POSITION; the reorder that would make that
safe (ancestry.py:244) is skipped when training from --ref_panel/--ref_labels. If the orders
differ, every variant is standardized and projected as a different one — counts match, nothing
raises, exit is 0, and the ancestry labels are wrong.

Reports, per callset:
  1. chromosome order of appearance — lexicographic (1,10,11,...,2,20) is the per-chromosome
     concatenation artifact
  2. whether positions ascend within each chromosome
  3. the decisive check: walking the callset in file order, do the variants it shares with the
     panel keep ascending panel rank? Zero rank drops means the positional rename is safe.
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
                continue                      # .pvar header line
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

    # Count rank DROPS, not a percentage: a chromosome-block permutation misplaces most columns
    # with only a few drops (lexicographic order drops twice), so a percentage reads as harmless.
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
