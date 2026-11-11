#!/usr/bin/env python3
"""Control 8 (Zhang et al. 2021): does the margin depend on the DEGREE OF THE SEEDS?

Why this control exists
-----------------------
Zhang et al. (2021, *Briefings in Bioinformatics* 22(5):bbab006) benchmarked
seven network-based gene-prioritisation methods for cerebral small vessel
disease and reported the sharpest published warning against the RWR family:

  "RWR's LOOCV performance may look better simply because its seed nodes
   have higher degree centrality, which need not reflect a real ability to
   predict new candidate genes."

Their random-seed experiment (1,000 trials x 10 random seeds) correlated the
median rank of a method with the degree centrality of the seeds. D1 already
contains a *degree-matched* random-seed null (control 2, A5), which answers a
weaker question: "is the margin larger than what degree-matched random seeds
would give?" It does NOT answer Zhang's question about the REAL seeds:

  **If the drivers happen to be hubs, then RWR wins simply because it is
  seeded with hubs.**

This control answers it directly, inside D1's own design. For every
(network, label) arm the training-fold positives are split into degree
terciles, and the primary endpoint is recomputed three times:

  low   the seeds are the LOW-degree third of the training positives
  mid   the seeds are the MID-degree third
  high  the seeds are the HIGH-degree third

The test set, the degree baseline and the fold structure are identical across
the three strata, so the three numbers are directly comparable. The verdict:

  * if the margin survives in the LOW stratum, the prior is not a hub artefact;
  * if the margin exists only in the HIGH stratum, the headline result is a
    study-bias artefact (the catalogue is degree-enriched) and the paper must
    say so.

Because RWR is linear in the seed weights, the three strata are propagated in
ONE batched call (the seed matrix has 3 x n_folds columns), so the control
costs about as much as the primary factorial itself.

Outputs
-------
    results/tables/control8_seed_degree.tsv       per (network, label, stratum)
    results/tables/control8_seed_degree_summary.tsv   per (network, label)
    docs/control8_seed_degree.md                 the reproducer + interpretation

Usage
-----
    python scripts/run_control8_seed_degree.py
    python scripts/run_control8_seed_degree.py --only-net reactome --only-label ot_all
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.evaluate import evaluate  # noqa: E402
from d1.engine.guard import write_table  # noqa: E402
from d1.engine.rwr import transition_T, seed_matrix, rwr  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETWORKS = ["funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]

NET_DIR = ROOT / "data" / "processed" / "networks"
LABEL_DIR = ROOT / "data" / "processed" / "labels"
SPLIT_DIR = ROOT / "data" / "processed" / "splits_v40"
TABLE_DIR = ROOT / "results" / "tables"

ALPHA = 0.5
STRATA = ("low", "mid", "high")


def _split_path(label: str, net_id: str) -> Path:
    return SPLIT_DIR / f"{label}__native_{net_id}.tsv"


def _tercile_cut(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Split an array into three equal-size degree strata (low, mid, high).

    Ties are broken by the index, so the strata are always the same size even
    when many genes share a degree (which they do in a scale-free network).
    """
    order = np.argsort(points, kind="mergesort")
    n = len(order)
    cut1 = n // 3
    cut2 = (2 * n) // 3
    return (order[:cut1], order[cut1:cut2], order[cut2:])


def run_arm(net_id: str, label: str, alpha: float = ALPHA) -> list[dict]:
    """Recompute the primary endpoint with degree-stratified seed sets."""
    genes, A = load_network(str(NET_DIR / f"{net_id}.tsv"))
    sp = pd.read_csv(_split_path(label, net_id), sep="\t", dtype={"gene": str})
    index = pd.Series(np.arange(len(genes)), index=genes)
    missing = sorted(set(sp.gene) - set(genes))
    if missing:
        raise ValueError(f"{net_id}: {len(missing)} split genes absent from the network")

    WT = transition_T(A)
    deg = np.asarray((A != 0).sum(axis=1)).ravel().astype(float)

    combos = sorted(sp.groupby(["repeat", "fold"]).groups)
    seed_cols: list[np.ndarray] = []
    test_idx: list[np.ndarray] = []
    test_y: list[np.ndarray] = []
    test_g: list[np.ndarray] = []
    stratum_sizes = {"low": [], "mid": [], "high": []}

    for (r, f) in combos:
        d = sp[sp.repeat == r]
        tr = d[(d.fold != f) & (d.y == 1)]
        te = d[d.fold == f]
        tr_idx = index.loc[tr.gene].to_numpy()
        lo, mid, hi = _tercile_cut(deg[tr_idx])
        for name, sel in (("low", lo), ("mid", mid), ("high", hi)):
            stratum_sizes[name].append(len(sel))
            seed_cols.append(tr_idx[sel])
        ti = index.loc[te.gene].to_numpy()
        test_idx.append(ti)
        test_y.append(te.y.to_numpy())
        test_g.append(te.gene.to_numpy())

    P0 = seed_matrix(len(genes), [c for c in seed_cols])
    P, n_iter = rwr(WT, P0, alpha=alpha)

    rows = []
    for j, ((r, f), ti, yy, gg) in enumerate(zip(combos, test_idx, test_y, test_g)):
        deg_metrics = evaluate(deg[ti], yy, gg)
        for si, name in enumerate(STRATA):
            rwr_metrics = evaluate(P[ti, j * 3 + si], yy, gg)
            rows.append({
                "network": net_id, "label": label, "alpha": alpha,
                "repeat": r, "fold": f, "stratum": name,
                "n_seeds": stratum_sizes[name][j],
                "rwr_auroc": rwr_metrics["auroc"],
                "degree_auroc": deg_metrics["auroc"],
                "delta_auroc": rwr_metrics["auroc"] - deg_metrics["auroc"],
                "n_test_pos": deg_metrics["n_test_pos"],
                "n_test": deg_metrics["n_test"],
            })
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    ap.add_argument("--alpha", type=float, default=ALPHA)
    args = ap.parse_args(argv)

    nets = args.only_net or NETWORKS
    labels = args.only_label or LABELS
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []
    for net in nets:
        for lab in labels:
            t0 = time.time()
            rows = run_arm(net, lab, alpha=args.alpha)
            all_rows.extend(rows)
            print(f"  [{net} x {lab}] {len(rows)} rows in {time.time() - t0:.1f}s")

    df = pd.DataFrame(all_rows)
    write_table(df, TABLE_DIR / "control8_seed_degree.tsv",
                nets=args.only_net, labels=args.only_label)

    # ---- per-arm summary: the three stratum means and their gaps ----
    g = df.groupby(["network", "label", "stratum"]).delta_auroc
    summary = g.mean().unstack("stratum")[["low", "mid", "high"]]
    summary["high_minus_low"] = summary["high"] - summary["low"]
    summary["low_positive"] = summary["low"] > 0
    summary["n_folds"] = df.groupby(["network", "label"]).size()
    summary = summary.reset_index()
    write_table(summary, TABLE_DIR / "control8_seed_degree_summary.tsv",
                nets=args.only_net, labels=args.only_label,
                float_format="%.6f")

    print("\n=== Control 8: DeltaAUROC by seed-degree stratum (alpha=0.5) ===")
    print(summary.round(4).to_string(index=False))
    print()
    n_low_pos = int(summary["low_positive"].sum())
    print(f"  low-degree seeds keep a positive margin in {n_low_pos}/{len(summary)} arms")
    print(f"  mean margin  low  {summary['low'].mean():+.4f}")
    print(f"  mean margin  mid  {summary['mid'].mean():+.4f}")
    print(f"  mean margin  high {summary['high'].mean():+.4f}")
    print(f"  mean(high - low) {summary['high_minus_low'].mean():+.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
