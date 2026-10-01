#!/usr/bin/env python3
"""Refuse to commit text that carries a participant or specimen ID.

    python3 scripts/id_guard.py --staged            # pre-commit: the staged content of every file
    python3 scripts/id_guard.py --message <file>    # commit-msg: the commit message
    python3 scripts/id_guard.py --all               # every tracked file, as it is on disk
    python3 scripts/id_guard.py <file> ...

.gitignore governs paths and nb_guard.py governs notebook outputs; neither sees an ID typed into
a doc, a comment or a commit message. This matches the ID shapes the cohorts use. Describe a
format with a placeholder (`R<7d>`, `PM-<site>_<n>`), never a real value.

A pattern check is a floor, not a proof: purely numeric IDs (Mayo, most DivCo individualIDs) have
no shape to match. Hits are reported by file, line and pattern — never the value itself.

Exit 0 clean, 1 on any hit or unreadable input, 2 on a usage error.
"""
import re
import subprocess
import sys

PATTERNS = {
    "ROSMAP individualID": r"\bR\d{7}\b",
    "ROSMAP WGS sample": r"\bMAP\d{7,}\b",
    "AMP-PD participant": r"\b(?:BF|HB|LB|LC|PD|PM|PP|SU|SY)-[A-Z]{0,4}[_-]?\d{3,}",
    "MSBB specimen": r"\bAMPAD_[A-Z]+_\d+",
    "DivCo specimen": r"\b\d{4,6}_DLPFC",
}
RX = [(name, re.compile(p)) for name, p in PATTERNS.items()]


def git(*args):
    return subprocess.run(["git", *args], capture_output=True, check=True).stdout


def scan(label, data):
    """-> [(label, line number, pattern name)]. Binary content is skipped."""
    if b"\0" in data[:8192]:
        return []
    hits = []
    for n, line in enumerate(data.decode("utf-8", errors="replace").splitlines(), 1):
        hits += [(label, n, name) for name, rx in RX if rx.search(line)]
    return hits


def main(argv):
    args = argv[1:]
    if not args:
        print(__doc__)
        return 2
    try:
        if args == ["--staged"]:
            paths = git("diff", "--cached", "--name-only", "--diff-filter=ACMR", "-z").decode().split("\0")
            inputs = [(p, git("show", f":{p}")) for p in paths if p]
        elif args[0] == "--message" and len(args) == 2:
            inputs = [("commit message", open(args[1], "rb").read())]
        elif args == ["--all"]:
            paths = git("ls-files", "-z").decode().split("\0")
            inputs = [(p, open(p, "rb").read()) for p in paths if p]
        else:
            inputs = [(p, open(p, "rb").read()) for p in args]
    except (OSError, subprocess.CalledProcessError) as e:
        print(f"id_guard: cannot read input: {e}", file=sys.stderr)
        return 1

    hits = [h for label, data in inputs for h in scan(label, data)]
    if hits:
        print("id_guard: REFUSING -- text shaped like a participant/specimen ID:", file=sys.stderr)
        for label, n, name in hits:
            print(f"    {label}:{n}  ({name})", file=sys.stderr)
        print("\nReplace each with a format placeholder, e.g. R<7d>, then re-stage.", file=sys.stderr)
        return 1

    print(f"id_guard: {len(inputs)} input(s) clean")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
