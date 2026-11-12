"""Normalise the line endings of the released tables to LF.

WHY: `pandas.DataFrame.to_csv` writes `os.linesep` as the row terminator
when it is given a path, so on Windows every table a script writes gets
CRLF endings. The repository's `.gitattributes` normalises `*.tsv` to LF,
so the committed blob and the working-tree file differ, and the
SHA-256 of the working tree does not match the SHA-256 a reader computes
from a checkout. That is what made 4 of the 6 fingerprints in
`docs/pre_registration_manifest.json` fail the reproducibility audit.

This script rewrites the affected files with LF endings, in place. It is
idempotent and it never touches the content of a file -- only the row
terminator.

Usage:
    python scripts/normalise_line_endings.py            # report
    python scripts/normalise_line_endings.py --write    # fix
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGETS = ["results/tables", "results/runs_v40", "results/runs_v40_clingen",
           "results/runs_v40_cx", "results/runs_pubcount_adj_v40",
           "review"]


def tracked_files() -> list[Path]:
    out = subprocess.check_output(["git", "ls-files"], cwd=str(ROOT), text=True)
    files = []
    for line in out.splitlines():
        p = ROOT / line
        if p.is_file() and p.suffix in {".tsv", ".csv", ".md", ".json", ".py"}:
            files.append(p)
    return files


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true",
                    help="rewrite the files with LF endings")
    args = ap.parse_args()

    changed = []
    for p in tracked_files():
        b = p.read_bytes()
        if b"\r\n" in b:
            changed.append(p)
            if args.write:
                p.write_bytes(b.replace(b"\r\n", b"\n"))

    if not changed:
        print("all tracked files already use LF endings")
        return 0
    print(f"{len(changed)} file(s) with CRLF endings")
    for p in changed[:40]:
        print(f"  {p.relative_to(ROOT)}")
    if len(changed) > 40:
        print(f"  ... and {len(changed) - 40} more")
    if args.write:
        print("rewritten with LF endings")
    else:
        print("run again with --write to fix them")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
