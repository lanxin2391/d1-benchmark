#!/usr/bin/env python3
"""Gene-cluster bootstrap CIs for every (network, label) cell.

Work-plan v3.1 task B1. Pre-registration section 9.3 asks for a cluster
bootstrap resampling GENES, not gene-arm rows.

Needs the per-gene score dumps written by task A3
(results/scores_v40/*.npz). For each cell and each alpha this computes the
gene-cluster CI and writes it to
results/tables/stats_bootstrap_ci_gene.tsv.

The three resampling units are reported side by side in
results/tables/stats_ci_all_units.tsv:
    fold        percentile bootstrap over the per-fold deltas
    repeat      percentile bootstrap over the per-repeat mean deltas
                (the primary; calibrated per work-plan A1/A2)
    gene        cluster bootstrap over genes (what the pre-registration
                names; the conservative bound)

Usage
-----
    python scripts/run_gene_cluster_ci.py \\
        --scores-dir results/scores_v40 \\
        --out results/tables/stats_bootstrap_ci_gene.tsv
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from d1.engine.intervals import gene_cluster_bootstrap  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scores-dir", default=str(REPO / "results" / "scores_v40"))
    ap.add_argument("--out", default=str(REPO / "results" / "tables" /
                                         "stats_bootstrap_ci_gene.tsv"))
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261029)
    ap.add_argument("--alphas", default="0.5",
                    help="comma-separated alphas to compute "
                         "(default 0.5, the primary endpoint)")
    args = ap.parse_args(argv)

    scores_dir = Path(args.scores_dir)
    want = [float(a) for a in args.alphas.split(",")]

    rows = []
    t0 = time.time()
    for f in sorted(scores_dir.glob("*.npz")):
        net, lab = f.stem.split("__")
        d = np.load(f, allow_pickle=True)
        alphas = d["alphas"]
        for a in want:
            ai = int(np.argmin(np.abs(alphas - a)))
            res = gene_cluster_bootstrap(
                d["rwr"][ai].T, d["degree"], d["y"].astype(int), d["in_test"],
                n_boot=args.n_boot, seed=args.seed)
            rows.append({
                "network": net, "label": lab, "alpha": round(float(alphas[ai]), 2),
                "n_folds": res["n_folds_used"], "n_boot": res["n_boot"],
                "observed_mean": res["observed_mean"],
                "observed_sd": res["observed_sd"],
                "boot_mean": res["boot_mean"], "boot_sd": res["boot_sd"],
                "ci_low": res["ci_low"], "ci_high": res["ci_high"],
                "excl0": bool(not (res["ci_low"] <= 0 <= res["ci_high"])),
                "width": res["ci_high"] - res["ci_low"],
            })
        print(f"  [{net} x {lab}] done  ({time.time()-t0:.0f}s elapsed)", flush=True)

    out = pd.DataFrame(rows).sort_values(["network", "label"])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    n_ok = int(out.excl0.sum())
    print()
    print(f"gene-cluster CI: {n_ok}/{len(out)} arms whose CI excludes 0")
    print(f"mean CI width: {out.width.mean():.4f}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
