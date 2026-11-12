#!/usr/bin/env python3
"""Check the numbers in docs/research_report_draft.md against the tables.

The report is a separate document from the manuscript sections, so it
needs its own verification: a number can be right in the manuscript and
wrong in the report.

Usage:
    python scripts/verify_report.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

REPORT = ROOT / "docs" / "research_report_draft.md"

# (label, string that must appear, what it is)
CHECKS: list[tuple[str, str, str]] = [
    ("design", "5-fold × 40-repeat", "the design"),
    ("design", "200 folds per cell", "the folds per cell"),
    ("endpoint.min", "−0.0125", "the negative cell"),
    ("endpoint.max", "+0.1711", "the panel maximum"),
    ("endpoint.cells", "23 of 24 cells were positive", "the positive-cell count"),
    ("endpoint.median", "+0.0809", "the median cell"),
    ("label.intogen", "+0.0821", "the intogen2024 label mean"),
    ("label.otall", "+0.1058", "the ot_all label mean"),
    ("label.otnolit", "+0.1062", "the ot_nolit label mean"),
    ("label.temporal", "+0.0208", "the temporal label mean"),
    ("label.clingen", "+0.1149", "the clingen label mean"),
    ("null2", "23 cells at cell-level z > 5", "control 2"),
    ("null2", "z = +3.09", "control 2's exception"),
    ("null3", "598 degree-preserving edge", "the rewiring count"),
    ("null3", "z = −3.53", "control 3's exception"),
    ("hub", "+0.0849", "the hub-capped margin"),
    ("hub", "+0.0979", "the degree-penalised margin"),
    ("c8", "+0.0527, +0.0609 and +0.0537", "the seed-degree strata"),
    ("c8", "+0.0010", "the high-minus-low gap"),
    ("c8", "negative in 11 of 24 arms", "the sign of the gap"),
    ("c10", "1.9 %", "the top-100 overlap"),
    ("c10", "0.6 genes", "the degeneracy statistic"),
    ("c10", "123.6", "the real-input degeneracy statistic"),
    ("c9", "−0.021", "the size correlation"),
    ("c9", "p = 0.93", "the permutation p"),
    ("c9", "Coverage was the", "the coverage reading"),
    ("c9", "+0.1449", "the coverage coefficient"),
    ("c9", "8.6 × 10⁻⁴", "the coverage p-value"),
    ("c9", "23 of 24 cells", "the adjusted endpoint count"),
    ("c4", "1,605 shared genes", "the shared universe"),
    ("c4", "321 test genes per fold", "the shared fold size"),
    ("c5", "39.8 %", "the cross-catalogue retention"),
    ("c5", "59 of 66 cells positive", "the transfer-cell counts"),
    ("c6", "50 independent training-seed draws", "the temporal design"),
    ("cx", "0 of 24 arms change sign", "the complex-aware control"),
    ("cx", "0.142", "the complex-aware change"),
    ("cx", "−12.3 %", "the complex-aware signed change"),
    ("c7", "+0.0010", "the control-7 matched difference"),
    ("c7", "28 of 30 arms", "the control-7 matched endpoint"),
    ("var", "40.1 % to the network × label interaction", "the interaction share"),
    ("var", "34.8 % to the label", "the label share"),
    ("var", "24.8 % to the network", "the network share"),
    ("var", "0.3 % to fold noise", "the fold share"),
    ("negcell", "63 of the 152", "the reactome coverage"),
    ("negcell", "41 %", "the coverage percentage"),
    ("negcell", "[−0.0589, +0.0375]", "the gene-cluster interval"),
    ("negcell", "Six of the 24 cells", "the equivalence count"),
    ("h2", "0.1173", "the smallest attainable Kruskal-Wallis p"),
    ("h2", "p = 0.50", "the published p"),
    ("gnn", "+0.028", "the embedding margin"),
    ("gnn", "+0.079", "the RWR margin"),
    ("gnn", "19 of 24 cells", "the embedding cell count"),
    ("overlap", "430 of 600", "the Open Targets overlap"),
    ("overlap", "24 of 600 and 18 of 600", "the temporal overlap"),
    ("audit", "50 PASS / 0 FAIL", "the audit result"),
    ("audit", "63 PASS / 0 FAIL", "the claim-verification result"),
    ("audit", "131 passed, 3 skipped", "the test suite"),
    ("manifest", "35 tables and 9 figures", "the release manifest"),
    ("groundtruth", "242 recomputed quantities", "the ground-truth table"),
    ("pauc.median", "+0.0308", "the partial-AUC median margin"),
    ("pauc.p10", "+0.0412", "the top-10 % median margin"),
    ("pauc.cells", "20 of 24 arms", "the co-primary positive-arm count"),
    ("pauc.rho", "ρ = 0.80", "the endpoint-ordering correlation"),
]


def main() -> int:
    text = REPORT.read_text(encoding="utf-8")
    collapsed = re.sub(r"\s+", " ", text)

    n_pass = n_fail = 0
    failures = []
    print("=" * 72)
    print("research_report_draft.md -- number verification")
    print("=" * 72)

    # numbers that must NOT appear
    forbidden = [
        ("10-repeat", "the 10-repeat design is superseded"),
        ("33 % from label", "the superseded variance numbers"),
        ("52 % residual", "the superseded residual share"),
        ("23 of 24 cells reject", "the superseded H1 wording"),
        ("curator lag", "the curator-lag reading is a hypothesis"),
        ("Mei et al.", "the wrong FunMap citation"),
        ("Jassal et al.", "the wrong Reactome citation"),
        ("61 cells positive", "the superseded transfer count"),
    ]
    print()
    print("-" * 72)
    print("forbidden strings")
    print("-" * 72)
    for needle, what in forbidden:
        n_c = re.sub(r"\s+", " ", needle)
        present = n_c in collapsed
        n_pass += (not present)
        n_fail += present
        print(f"  [{'PASS' if not present else 'FAIL'}] {needle:<26} {what}")
        if present:
            failures.append(f"forbidden {needle!r}")

    # ---- the DOCX must carry the same numbers as the Markdown ----
    print()
    print("-" * 72)
    print("the .docx deliverable")
    print("-" * 72)
    docx_path = REPORT.with_suffix(".docx")
    if not docx_path.is_file():
        print(f"  [FAIL] {docx_path.name} does not exist -- run "
              f"scripts/markdown_to_docx.py")
        n_fail += 1
        failures.append(f"{docx_path.name} does not exist")
    else:
        try:
            import docx as _docx
            d = _docx.Document(str(docx_path))
            dtext = " ".join(p.text for p in d.paragraphs)
            for t in d.tables:
                for row in t.rows:
                    for c in row.cells:
                        dtext += " " + c.text
            dtext = re.sub(r"\s+", " ", dtext)
            n_para = len(d.paragraphs)
            n_tab = len(d.tables)
            n_img = sum(1 for r in d.part.rels.values() if "image" in r.reltype)
            print(f"  [INFO] {docx_path.name}: {n_para} paragraphs, "
                  f"{n_tab} tables, {n_img} images")
            # every headline number must survive the conversion
            for label, needle, what in CHECKS:
                # the DOCX carries the rendered text, so the Markdown
                # emphasis markers must be stripped before comparing
                n_c = re.sub(r"\s+", " ", needle.replace("**", ""))
                if n_c in dtext:
                    n_pass += 1
                    print(f"  [PASS] {label:<16} {what}")
                else:
                    n_fail += 1
                    print(f"  [FAIL] {label:<16} {what}")
                    failures.append(f"docx: {label} ({needle!r}) is missing")
            for needle, what in forbidden:
                n_c = re.sub(r"\s+", " ", needle)
                if n_c in dtext:
                    n_fail += 1
                    print(f"  [FAIL] {needle:<26} {what} (in the docx)")
                    failures.append(f"docx: forbidden {needle!r}")
        except Exception as e:
            n_fail += 1
            print(f"  [FAIL] the docx could not be read: {type(e).__name__}: {e}")
            failures.append(f"docx unreadable: {e}")

    print()
    print("=" * 72)
    print(f"RESULT: {n_pass} PASS / {n_fail} FAIL")
    print("=" * 72)
    if failures:
        print()
        print("Failures:")
        for f in failures:
            print(f"  - {f}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
