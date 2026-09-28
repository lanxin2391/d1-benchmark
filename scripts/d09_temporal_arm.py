"""D-09 temporal-drift arm.

The protocol's §4.3 says "seeds = training-fold positives". For D-09 we
want to test the harder question: does the network help us **find
newly-discovered drivers** that were not in our training set?

Setup:
* universe      = (intogen2024 ∩ LCC) ∪ (intogen_temporal_new ∩ LCC)
                  = intogen2024 ∩ LCC (since temporal_new ⊂ intogen2024)
* seeds (train) = intogen2024 ∩ LCC minus intogen_temporal_new
                  = the 481 - 357 = 357 known drivers in the LCC, on which
                    we learn a network prior
* test target   = intogen_temporal_new ∩ LCC (104 genes we are trying to
                  recover). Evaluated against a held-out set of LCC
                  negatives that are NOT in either training or target
                  (background of LCC minus any intogen2024 driver).

Because the test target is small (~104) and folds are stratified over the
*full* intogen2024 driver set, we instead do **bootstrap evaluation**: we
sample k (= 100) random draws from the universe minus training drivers,
rank them, and compute AUROC of (target seeds vs sampled background).
This is what the protocol's "control 6" reduces to when the test target
is disjoint from training.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from d1.engine.baselines import degree_scores
from d1.engine.rwr import ALPHA, rwr, seed_matrix, transition_T
from d1.networks.io import load_network

LABELS_DIR = REPO_ROOT / "data" / "processed" / "labels"
RESULTS_DIR = REPO_ROOT / "results" / "runs"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ARGS in cli
def _load_labels(label_id: str) -> set[str]:
    p = LABELS_DIR / f"{label_id}.tsv"
    return set(pd.read_csv(p, sep="\t").gene.tolist())


def _to_idx(genes_sorted: np.ndarray, genes_set: set[str]) -> np.ndarray:
    """Map gene names to integer indices in genes_sorted."""
    pos = {g: i for i, g in enumerate(genes_sorted)}
    return np.asarray([pos[g] for g in genes_set if g in pos], dtype=int)


def temporal_arm(network_id: str,
                 bootstrap_k: int = 100,
                 bootstrap_n: int = 200,
                 alpha: float = ALPHA,
                 seed: int = 20260926) -> pd.DataFrame:
    """One bootstrap-evaluated D-09 arm; returns the long-format C5 result."""
    net_path = REPO_ROOT / "data" / "processed" / "networks" / f"{network_id}.tsv"
    genes, A = load_network(net_path)
    n = len(genes)

    intogen2024 = _load_labels("intogen2024")
    temporal = _load_labels("intogen_temporal_new")
    assert temporal <= intogen2024, "temporal_new must be a subset of intogen2024"

    # restrict to network LCC (genes are the LCC nodes, since A is the LCC adj)
    lcc_genes = set(genes.tolist())
    train_pool = (intogen2024 - temporal) & lcc_genes
    test_target = temporal & lcc_genes
    # background: anything in LCC that is NOT in intogen2024 (so neither train nor test)
    background = lcc_genes - intogen2024

    if not test_target:
        raise ValueError("no temporal_new positives in network LCC")
    if not train_pool:
        raise ValueError("no training drivers in network LCC")

    train_idx = _to_idx(genes, train_pool)
    test_idx = _to_idx(genes, test_target)
    bg_idx_all = _to_idx(genes, background)

    # 1) Compute RWR prior (vector form, single seed vector)
    WT = transition_T(A)
    P0 = seed_matrix(n, [train_idx])      # (n, 1)
    P, n_iter = rwr(WT, P0, alpha=alpha)
    scores = P[:, 0]                       # (n,) per-gene RWR score
    deg = degree_scores(A)
    rows = []
    rng = np.random.default_rng(seed)

    for k in range(bootstrap_k):
        # draw bootstrap_n background genes (with replacement, like sklearn)
        sampled = rng.choice(bg_idx_all, size=bootstrap_n, replace=True)
        idx = np.concatenate([test_idx, sampled])
        y = np.concatenate([np.ones(len(test_idx), dtype=int),
                            np.zeros(len(sampled), dtype=int)])
        # unique gene names for precision_at_k ties
        gene_names = genes[idx]
        auroc = float(roc_auc_score(y, scores[idx]))
        auprc = float(average_precision_score(y, scores[idx]))
        rows.append(dict(network=network_id, label_set="intogen_temporal_new",
                         universe="native_temporal_eval", alpha=alpha,
                         repeat=k, fold=0, method="rwr", n_iter=n_iter,
                         auroc=auroc, auprc=auprc, n_test_pos=int(y.sum()),
                         n_test=int(len(y)),
                         n_train_pool=len(train_pool),
                         n_target=len(test_target)))
        # baseline (degree)
        auroc_deg = float(roc_auc_score(y, deg[idx]))
        auprc_deg = float(average_precision_score(y, deg[idx]))
        rows.append(dict(network=network_id, label_set="intogen_temporal_new",
                         universe="native_temporal_eval", alpha=alpha,
                         repeat=k, fold=0, method="degree", n_iter=np.nan,
                         auroc=auroc_deg, auprc=auprc_deg, n_test_pos=int(y.sum()),
                         n_test=int(len(y)),
                         n_train_pool=len(train_pool),
                         n_target=len(test_target)))

    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--networks", nargs="+",
                    default=["funmap", "string_full700"])
    ap.add_argument("--bootstrap-k", type=int, default=100)
    ap.add_argument("--bootstrap-n", type=int, default=200)
    ap.add_argument("--alpha", type=float, default=ALPHA)
    ap.add_argument("--out", default=str(RESULTS_DIR / "d09_temporal.tsv"))
    args = ap.parse_args()

    all_rows = []
    for net_id in args.networks:
        t = time.time()
        print(f"=== {net_id} ===")
        df = temporal_arm(net_id,
                          bootstrap_k=args.bootstrap_k,
                          bootstrap_n=args.bootstrap_n,
                          alpha=args.alpha)
        all_rows.append(df)
        rwr = df[df.method == "rwr"]
        deg = df[df.method == "degree"]
        delta = rwr.auroc.values - deg.auroc.values
        print(f"  n_train_pool = {rwr.n_train_pool.iloc[0]}, "
              f"n_target = {rwr.n_target.iloc[0]}")
        print(f"  RWR AUROC  = {rwr.auroc.mean():.4f} ± {rwr.auroc.std():.4f}")
        print(f"  degree     = {deg.auroc.mean():.4f} ± {deg.auroc.std():.4f}")
        print(f"  delta      = {delta.mean():.4f} ± {delta.std():.4f}  "
              f"(k={len(delta)} bootstraps)")
        print(f"  time       = {time.time()-t:.1f}s")
    out = pd.concat(all_rows, ignore_index=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
