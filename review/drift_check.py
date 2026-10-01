#!/usr/bin/env python3
"""Does the working tree differ from a revision in anything that RUNS?

Read-only. Each modified file (.md skipped) is compared against <rev>, per language:
    .py   AST comparison with docstrings stripped (comments are not in an AST)
    .sh   comments removed, trailing ones included, quote-aware, then compared
    other byte comparison
A whole-line comment strip is not enough: it reads a rewritten docstring or a realigned trailing
comment as a code change. Note that #SBATCH directives are comments to this test, so a change to
one is NOT reported.

Exit 0 if every difference is commentary, 1 if any file differs in code (named). A file that
cannot be read or parsed counts as a difference, never as a pass.

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
        # Unknowns count as failures: an unread or unparsed file is not a clean one.
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
