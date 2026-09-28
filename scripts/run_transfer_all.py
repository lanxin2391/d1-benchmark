"""W3 Control 5 — runner that executes all 12 transfer arm pairs across
6 networks (= up to 72 arms; missing cells are skipped if the seed
pool or target pool is empty).

Output:
  results/runs/transfer/{label_a}_to_{label_b}__{network}.tsv
  results/tables/results_stats_transfer.tsv (one row per arm:
    network, label_a, label_b, n_target, n_seeds, mean_delta_auroc, sd,
    t, p, full_set_delta, full_set_auroc_rwr, full_set_auroc_deg)

Usage::

    python scripts/run_transfer_all.py
    python scripts/run_transfer_all.py --n-bootstrap 50
    python scripts/run_transfer_all.py --networks funmap string_full700
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

LABELS_DIR = REPO_ROOT / "data" / "processed" / "labels"
NETWORKS_DIR = REPO_ROOT / "data" / "processed" / "networks"
TRANSFER_DIR = REPO_ROOT / "data" / "processed" / "transfer"
RUNS_DIR = REPO_ROOT / "results" / "runs" / "transfer"
TABLES_DIR = REPO_ROOT / "results" / "tables"
RUNS_DIR.mkdir(parents=True, exist_ok=True)
TABLES_DIR.mkdir(parents=True, exist_ok=True)

DEFAULT_LABEL_POOL = ("intogen2024", "intogen_temporal_new", "ot_all",
                      "ot_nolit")
DEFAULT_NETWORKS = ("funmap", "string_full700", "string_phys700",
                    "intact", "reactome", "rna_coexp")


def run_one(label_a: str, label_b: str, network: str,
            n_bootstrap: int, alpha: float = 0.5, seed: int = 20260926,
            ) -> dict | None:
    """Run one transfer arm. Returns a stats row or None if skipped."""
    sp = TRANSFER_DIR / f"{label_a}__to__{label_b}__{network}.tsv"
    if not sp.exists():
        return None
    out = RUNS_DIR / f"{label_a}_to_{label_b}__{network}.tsv"
    np_path = NETWORKS_DIR / f"{network}.tsv"

    # inline-import to keep this script loadable for `--list`
    from d1.engine.rwr import rwr, seed_matrix, transition_T
    from d1.networks.io import load_network
    from sklearn.metrics import roc_auc_score

    transfer_df = pd.read_csv(sp, sep="\t", dtype=str)
    seeds_df = pd.read_csv(LABELS_DIR / f"{label_a}.tsv",
                           sep="\t", dtype=str)[["gene"]].dropna()

    genes, A = load_network(np_path)
    index = pd.Series(np.arange(len(genes)), index=genes)
    target_genes = sorted(set(transfer_df[transfer_df.role == "target"].gene))
    background_genes = sorted(set(transfer_df[transfer_df.role == "background"].gene))
    network_lcc = set(genes)
    seed_pool = sorted((set(seeds_df.gene) & network_lcc) - set(target_genes))

    if not target_genes:
        return None
    if not seed_pool:
        return None

    target_idx = index.loc[target_genes].to_numpy()
    # Build the *real* background pool from the network LCC minus both
    # labels minus the seed pool. The transfer file's 1000 sampled bg
    # is for storage; for evaluation we want every non-(target ∪ seeds
    # ∪ label_a ∪ label_b) gene in the LCC.
    label_a_genes = set(pd.read_csv(LABELS_DIR / f"{label_a}.tsv",
                                    sep="\t", dtype=str).gene)
    label_b_genes = set(pd.read_csv(LABELS_DIR / f"{label_b}.tsv",
                                    sep="\t", dtype=str).gene)
    real_bg_genes = sorted(network_lcc - (label_a_genes | label_b_genes
                                          | set(target_genes)
                                          | set(seed_pool)))
    bg_idx_pool = index.loc[real_bg_genes].to_numpy()
    seed_idx = index.loc[seed_pool].to_numpy()

    WT = transition_T(A)
    P0 = seed_matrix(len(genes), [seed_idx])
    P, n_iter = rwr(WT, P0, alpha=alpha)
    rwr_scores = np.asarray(P[:, 0]).ravel()
    deg = np.asarray(A.sum(axis=1)).ravel()

    target_mask = np.zeros(len(genes), dtype=bool)
    target_mask[target_idx] = True
    bg_pool_size = len(bg_idx_pool)
    bg_size = min(1000, bg_pool_size) if bg_pool_size > 0 else 0

    if bg_size == 0:
        return None

    rng = np.random.default_rng(seed)
    rows = []
    base = {"network": network, "label_a": label_a, "label_b": label_b,
            "alpha": alpha, "n_iter": n_iter,
            "n_target": len(target_idx), "n_seeds": len(seed_idx)}
    for f in range(n_bootstrap):
        sel = rng.choice(bg_idx_pool, size=bg_size, replace=False)
        idx = np.concatenate([target_idx, sel])
        y = np.zeros(len(genes), dtype=int)
        y[target_mask] = 1
        a_rwr = float(roc_auc_score(y[idx], rwr_scores[idx]))
        a_deg = float(roc_auc_score(y[idx], deg[idx]))
        rows.append({**base, "fold": f, "method": "rwr", "auroc": a_rwr,
                     "n_test": int(len(idx))})
        rows.append({**base, "fold": f, "method": "degree", "auroc": a_deg,
                     "n_test": int(len(idx))})

    # full-set
    y_full = np.zeros(len(genes), dtype=int)
    y_full[target_mask] = 1
    a_rwr_full = float(roc_auc_score(y_full, rwr_scores))
    a_deg_full = float(roc_auc_score(y_full, deg))

    out_df = pd.DataFrame(rows)
    out_df.to_csv(out, sep="\t", index=False)

    boot = out_df
    rwr_aurocs = boot[boot.method == "rwr"].auroc.values
    deg_aurocs = boot[boot.method == "degree"].auroc.values
    delta = rwr_aurocs - deg_aurocs
    t_stat, p_val = stats.ttest_1samp(delta, popmean=0, alternative="greater")

    return {
        **base,
        "n_bootstrap": n_bootstrap,
        "n_background": bg_size,
        "mean_delta_auroc": float(delta.mean()),
        "sd_delta_auroc": float(delta.std(ddof=1)),
        "t_statistic": float(t_stat),
        "p_value_one_sided": float(p_val),
        "full_set_auroc_rwr": a_rwr_full,
        "full_set_auroc_deg": a_deg_full,
        "full_set_delta_auroc": a_rwr_full - a_deg_full,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--networks", nargs="+", default=list(DEFAULT_NETWORKS))
    ap.add_argument("--labels", nargs="+", default=list(DEFAULT_LABEL_POOL))
    ap.add_argument("--n-bootstrap", type=int, default=100)
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=20260926)
    ap.add_argument("--out-table", type=Path,
                    default=TABLES_DIR / "results_stats_transfer.tsv")
    args = ap.parse_args(argv)

    pairs = []
    for la in args.labels:
        for lb in args.labels:
            if la != lb:
                pairs.append((la, lb))
    print(f"Running {len(pairs)} ordered label pairs × {len(args.networks)} "
          f"networks = {len(pairs) * len(args.networks)} arms\n")

    rows = []
    t0 = time.time()
    for i, (la, lb) in enumerate(pairs, start=1):
        for net in args.networks:
            try:
                r = run_one(la, lb, net, n_bootstrap=args.n_bootstrap,
                            alpha=args.alpha, seed=args.seed)
            except Exception as e:
                print(f"  [skip] {la}->{lb} on {net}: {e}")
                continue
            if r is None:
                continue
            rows.append(r)
            print(f"  [{i:3d}/{len(pairs) * len(args.networks)}] "
                  f"{la:>22} -> {lb:<22} on {net:<14}: "
                  f"ΔAUROC = {r['mean_delta_auroc']:+.4f} ± "
                  f"{r['sd_delta_auroc']:.4f}  "
                  f"p = {r['p_value_one_sided']:.2g}  "
                  f"(n_target={r['n_target']}, n_seeds={r['n_seeds']})")

    print(f"\nTotal runtime: {time.time() - t0:.1f}s")
    if not rows:
        print("No arms completed")
        return 1
    table = pd.DataFrame(rows)
    table = table.sort_values("mean_delta_auroc", ascending=False)
    table.to_csv(args.out_table, sep="\t", index=False, float_format="%.6f")
    print(f"Wrote {args.out_table} ({len(table)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())