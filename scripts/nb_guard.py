#!/usr/bin/env python3
"""Refuse to commit a notebook that carries stored cell outputs.

    python3 scripts/nb_guard.py [--staged | <notebook> ...]

WHY. On 2026-08-21 `demo_sample_check.ipynb` was found holding a rendered DataFrame of
individual_id x callset membership -- subject-level controlled-access data, in git, in a
repo about to get its first remote. The history was rewritten before the push. Nothing
had flagged it, because .gitignore governs paths and this leaked inside a tracked file.

Notebook SOURCE is reproducible and belongs in git. Notebook OUTPUTS are a rendering of
whatever the data happened to be, so they inherit the DUA that covers the data. This
enforces that split by refusing outputs entirely rather than trying to judge which ones
are aggregate -- a judgement is exactly what failed here.

Exit 0 clean, 1 on any notebook with outputs, 2 on a usage error.
"""
import json
import subprocess
import sys


def outputs_in(path):
    """-> number of stored outputs. A notebook git cannot parse is a failure, not a pass."""
    with open(path) as fh:
        nb = json.load(fh)
    return sum(len(c.get("outputs", [])) for c in nb.get("cells", []))


def staged_notebooks():
    out = subprocess.run(["git", "diff", "--cached", "--name-only", "--diff-filter=ACM"],
                         capture_output=True, text=True, check=True).stdout
    return [p for p in out.split("\n") if p.endswith(".ipynb")]


def main(argv):
    args = argv[1:]
    if not args:
        print(__doc__)
        return 2
    paths = staged_notebooks() if args == ["--staged"] else args

    dirty = []
    for p in paths:
        try:
            n = outputs_in(p)
        except (OSError, ValueError) as e:
            print(f"nb_guard: cannot read {p}: {e}", file=sys.stderr)
            return 1
        if n:
            dirty.append((p, n))

    if dirty:
        print("nb_guard: REFUSING -- these notebooks carry stored outputs:", file=sys.stderr)
        for p, n in dirty:
            print(f"    {p}  ({n} outputs)", file=sys.stderr)
        print("\nStrip them, then re-stage:\n"
              "    python3 scripts/nb_guard.py --strip " + " ".join(p for p, _ in dirty),
              file=sys.stderr)
        return 1

    print(f"nb_guard: {len(paths)} notebook(s) clean")
    return 0


def strip(paths):
    for p in paths:
        with open(p) as fh:
            nb = json.load(fh)
        n = 0
        for c in nb.get("cells", []):
            if c.get("cell_type") == "code":
                n += len(c.get("outputs", []))
                c["outputs"] = []
                c["execution_count"] = None
        with open(p, "w") as fh:
            json.dump(nb, fh, indent=1, ensure_ascii=False)
            fh.write("\n")
        print(f"stripped {n} outputs from {p}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--strip":
        sys.exit(strip(sys.argv[2:]))
    sys.exit(main(sys.argv))
