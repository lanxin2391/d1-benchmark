"""Refresh the audit counts quoted in the documents.

The claim count, the ground-truth size and the test count move every
time a control is added, and three documents quote them. This script
recomputes them and rewrites the quotes, so no document can drift away
from the audit.
"""
from __future__ import annotations

import pathlib
import re
import subprocess
import sys

sys.stdout.reconfigure(encoding="utf-8")
ROOT = pathlib.Path(".")

# ---- 1. the current counts -------------------------------------------------
r = subprocess.run([sys.executable, "scripts/verify_claims.py"],
                   capture_output=True, text=True, cwd=str(ROOT),
                   encoding="utf-8", errors="replace")
m = re.search(r"RESULT: (\d+) PASS / (\d+) FAIL", r.stdout or "")
if not m:
    print("verify_claims did not report a result; nothing to do")
    raise SystemExit(1)
n_claims, n_fail = int(m.group(1)), int(m.group(2))
if n_fail:
    print(f"verify_claims has {n_fail} FAILs; fix them before refreshing")
    raise SystemExit(1)

gt = ROOT / "review/ground_truth.tsv"
n_truth = sum(1 for _ in gt.read_text(encoding="utf-8").splitlines()) - 1

r = subprocess.run([sys.executable, "-m", "pytest", "-q"],
                   capture_output=True, text=True, cwd=str(ROOT),
                   encoding="utf-8", errors="replace")
m = re.search(r"(\d+) passed", r.stdout or "")
n_tests = int(m.group(1)) if m else 0
print(f"the current counts: {n_claims} claim checks, {n_truth} ground-truth "
      f"quantities, {n_tests} tests")

# ---- 2. the documents that quote them --------------------------------------
EDITS = [
    ("docs/research_report_draft.md",
     [(r"`python scripts/verify_claims\.py` \((\d+) PASS / 0 FAIL\)",
       f"`python scripts/verify_claims.py` ({n_claims} PASS / 0 FAIL)"),
      (r"manuscript against that table \((\d+) PASS / 0 FAIL\)",
       f"manuscript against that table ({n_claims} PASS / 0 FAIL)"),
      (r"a ground-truth table of (\d+) recomputed quantities",
       f"a ground-truth table of {n_truth} recomputed quantities"),
      (r"`review/ground_truth\.tsv` \((\d+) recomputed quantities\)",
       f"`review/ground_truth.tsv` ({n_truth} recomputed quantities)")]),
    ("docs/manuscript_a_sections.md",
     [(r"two audit scripts \((\d+) and 71\nchecks\) and a ground-truth "
       r"table of (\d+) recomputed quantities",
       f"two audit scripts ({n_claims} and 71\nchecks) and a ground-truth "
       f"table of {n_truth} recomputed quantities")]),
    ("docs/report_1118_speaker_script_20min.md",
     [(r"\u4e24\u4e2a\u5ba1\u8ba1\u811a\u672c\uff08(\d+) \u9879\u4e0e 71 \u9879\uff09",
       f"\u4e24\u4e2a\u5ba1\u8ba1\u811a\u672c\uff08{n_claims} \u9879\u4e0e 71 \u9879\uff09"),
      (r"(\d+) \u4e2a\u91cd\u7b97\u91cf", f"{n_truth} \u4e2a\u91cd\u7b97\u91cf")]),
    ("scripts/verify_report.py",
     [(r'\("audit", "(\d+) PASS / 0 FAIL", "the claim-verification result"\)',
       f'("audit", "{n_claims} PASS / 0 FAIL", "the claim-verification result")'),
      (r'\("groundtruth", "(\d+) recomputed quantities", "the ground-truth table"\)',
       f'("groundtruth", "{n_truth} recomputed quantities", "the ground-truth table")'),
      (r'\("audit", "(\d+) passed", "the test suite"\)',
       f'("audit", "{n_tests} passed", "the test suite")')]),
]
for path, pairs in EDITS:
    p = ROOT / path
    t = p.read_text(encoding="utf-8")
    for rx, new in pairs:
        t, k = re.subn(rx, new.replace("\\", "\\\\"), t, count=1)
        if k == 0:
            print(f"  {path}: no match for {rx[:48]!r}")
        else:
            print(f"  {path}: {new[:60]}")
    p.write_text(t, encoding="utf-8")
print("the documents quote the current audit")
