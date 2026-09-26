"""Run one network x label set x universe arm on real data (used from S1 onward).

Example (the S1 smoke run):
    python scripts/run_arm.py ^
        --network data/processed/networks/funmap.tsv ^
        --splits  data/processed/splits/intogen2024__native_funmap.tsv ^
        --network-id funmap --label-id intogen2024 --universe-id native_funmap ^
        --out results/runs/smoke_funmap_intogen2024.tsv
(In PowerShell use a backtick ` instead of ^ at line ends, or write it on one line.)
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from d1.engine.evaluate import run_arm, summarize  # noqa: E402
from d1.engine.rwr import ALPHA  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--network", required=True, help="C2 network file")
    ap.add_argument("--splits", required=True, help="C4 split file")
    ap.add_argument("--network-id", required=True)
    ap.add_argument("--label-id", required=True)
    ap.add_argument("--universe-id", required=True)
    ap.add_argument("--alpha", type=float, nargs="+", default=[ALPHA])
    ap.add_argument("--out", required=True, help="C5 result table to write")
    a = ap.parse_args()

    t = time.time()
    res = run_arm(a.network, a.splits, a.network_id, a.label_id, a.universe_id,
                  alphas=tuple(a.alpha))
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    res.to_csv(a.out, sep="\t", index=False)
    print(summarize(res).round(4).to_string(index=False))
    print(f"\n{len(res)} rows written to {a.out} in {time.time() - t:.1f} s")


if __name__ == "__main__":
    main()
