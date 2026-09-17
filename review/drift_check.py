#!/usr/bin/env python3
"""Does the working tree differ from HEAD in anything that RUNS?

Read-only. Written 2026-09-17, when the cluster turned out never to have been a git checkout
and sixteen files showed as modified. Two naive tests both gave wrong answers there:

  - matching a file's diff line-count against a suspected commit. Seven of nine matched exactly
    and it was still wrong — analysis_grain.py matched nothing and was the one file with a real
    difference.
  - stripping whole-line comments. That reports a rewritten Python DOCSTRING as a code change,
    and a REALIGNED TRAILING comment in shell as a code change. Both happened; both were noise.

So the test is per language:
    .py   AST comparison with docstrings stripped — comments are not in an AST at all
    .sh   comments removed INCLUDING trailing ones, quote-aware, then compared

Exit 0 if every difference is commentary, 1 if any file differs in code (and it names them).
A file that cannot be parsed is reported as a difference, never as a pass: an unparseable file
is an unanswered question, not a clean one.

    python3 review/drift_check.py              # vs HEAD
    python3 review/drift_check.py <rev>        # vs any revision
"""
import argparse
import ast
import re
import subprocess
import sys

TEXT_EXT = (".sh", ".bash")


def git(*args):
    r = subprocess.run(["git", *args], capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None


def strip_docstrings(tree):
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = node.body
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                node.body = body[1:]
    return tree


def shell_code(text):
    """Executable content of a shell script: no blank lines, no comments, trailing ones included.

    Quote-aware, because a '#' inside a string is not a comment — plink's --out stems and the
    printf format strings in this project both contain them.
    """
    out = []
    for line in text.splitlines():
        if re.match(r"^\s*(#|$)", line):
            continue
        quote, cut = None, len(line)
        for i, ch in enumerate(line):
            if quote:
                if ch == quote:
                    quote = None
            elif ch in "\"'":
                quote = ch
            elif ch == "#":
                cut = i
                break
        stripped = line[:cut].rstrip()
        if stripped:
            out.append(stripped)
    return "\n".join(out)


def compare(path, rev):
    """-> (same, note). `same` is None when the comparison could not be made."""
    old = git("show", f"{rev}:{path}")
    if old is None:
        return None, f"not in {rev}"
    try:
        new = open(path).read()
    except OSError as e:
        return None, str(e)

    if path.endswith(".py"):
        try:
            a = ast.dump(strip_docstrings(ast.parse(old)))
            b = ast.dump(strip_docstrings(ast.parse(new)))
        except SyntaxError as e:
            return None, f"parse error: {e}"
        return a == b, "AST"
    if path.endswith(TEXT_EXT):
        return shell_code(old) == shell_code(new), "shell"
    return old == new, "bytes"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rev", nargs="?", default="HEAD")
    ap.add_argument("--all", action="store_true",
                    help="every tracked file, not only the ones git reports as modified")
    a = ap.parse_args()

    listing = git("ls-files") if a.all else git("diff", "--name-only", a.rev)
    if listing is None:
        print(f"ERROR: cannot list files against {a.rev} — is this a git repository?",
              file=sys.stderr)
        return 2
    files = [f for f in listing.split() if not f.endswith(".md")]
    if not files:
        print(f"No modified non-doc files against {a.rev}.")
        return 0

    differ, unknown = [], []
    for f in sorted(files):
        same, note = compare(f, a.rev)
        if same is None:
            unknown.append((f, note))
            print(f"  CANNOT COMPARE    {f}  ({note})")
        elif same:
            print(f"  LOGIC IDENTICAL   {f}")
        else:
            differ.append(f)
            print(f"  LOGIC CHANGED     {f}")

    print()
    if differ or unknown:
        # Unknowns count as failures on purpose. A file this could not read or parse is an
        # unanswered question; reporting it alongside the passes would make "none found"
        # printable off a check that did not run.
        for f in differ:
            print(f"    git diff {a.rev} -- {f}")
        if unknown:
            print(f"  {len(unknown)} file(s) could NOT be compared — not the same as clean:")
            for f, note in unknown:
                print(f"    {f}  ({note})")
        print(f"\n{len(differ)} changed, {len(unknown)} uncomparable. Read each before overwriting.")
        return 1
    print(f"All {len(files)} modified non-doc file(s) differ only in commentary.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
