#!/usr/bin/env python3
"""Control 7 comparison table: the original label vs the publication-adjusted label.

FIX (2026-11-11). The previous version of this script

* hard-coded an absolute path (``sys.path.insert(0, r'D:\\Bioinformatics\\d1-benchmark')``)
  that does not exist on any other machine, and
* read the **10-repeat** run files (``results/runs/end_to_end_*.tsv`` and
  ``results/runs/pubcount_adj/``), so the table it wrote was on the
  superseded design.

The released design is the 40-repeat one:

* original label      -> ``results/runs_v40/a4_{net}__{label}.tsv``
* adjusted label      -> ``results/runs_pubcount_adj_v40/a4_{net}__{label}_pubcount_adjusted.tsv``

and the comparison table is now written from those, so it covers all
**30** arms (6 networks x 5 label sets) on the design of record.

Usage
-----
    python scripts/compare_pubcount_adj.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.evaluate import margins  # noqa: E402

LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new",
          "clingen_label"]
NETS = ["funmap", "string_full700", "string_phys700", "intact", "reactome",
        "rna_coexp"]

ORIG_DIR = ROOT / "results" / "runs_v40"
ORIG_CLINGEN_DIR = ROOT / "results" / "runs_v40_clingen"
ADJ_DIR = ROOT / "results" / "runs_pubcount_adj_v40"
OUT = ROOT / "results" / "tables" / "pubcount_adj_comparison.tsv"


def _orig_path(label: str, net: str) -> Path:
    """The original-label run file for a (label, network) arm."""
    if label == "clingen_label":
        return ORIG_CLINGEN_DIR / f"a4_{net}__{label}.tsv"
    return ORIG_DIR / f"a4_{net}__{label}.tsv"


def _delta(path: Path) -> tuple[float, float, int]:
    """(mean delta_auroc, sd, n folds) at alpha = 0.5."""
    m = margins(pd.read_csv(path, sep="\t"))
    m = m[m.alpha == 0.5]
    return (float(m.delta_auroc.mean()), float(m.delta_auroc.std(ddof=1)),
            int(len(m)))


def main() -> int:
    rows = []
    for lbl in LABELS:
        for net in NETS:
            orig = _orig_path(lbl, net)
            adj = ADJ_DIR / f"a4_{net}__{lbl}_pubcount_adjusted.tsv"
            if not orig.exists() or not adj.exists():
                print(f"  missing for {lbl} x {net}: "
                      f"orig={orig.exists()}, adj={adj.exists()}")
                continue
            o_d, o_sd, o_n = _delta(orig)
            a_d, a_sd, a_n = _delta(adj)
            rows.append({
                "label": lbl, "network": net,
                "orig_dauroc": o_d, "orig_sd": o_sd,
                "adj_dauroc": a_d, "adj_sd": a_sd,
                "diff": a_d - o_d,
                "orig_n_folds": o_n, "adj_n_folds": a_n,
            })

    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False, float_format="%.6f", lineterminator="\n")

    print(df.to_string(index=False))
    print()
    print(f"arms: {len(df)}   mean diff: {df['diff'].mean():+.4f}   "
          f"adjusted < original in {int((df['diff'] < 0).sum())}/{len(df)}")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
