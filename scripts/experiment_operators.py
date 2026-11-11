"""Experiment: which shrinkage operator actually reduces the per-gene rank SD?

The pre-registered operator (shrink toward the median of the *score*
neighbourhood) turned out to be order-preserving and therefore a no-op on the
ranking. This script compares candidate operators on one arm so that the
amended operator is chosen on evidence rather than on taste.

Operators compared (all on the same folds, same confidence weights):

* ``knn_score``   -- pre-registered: k nearest by score value.
* ``knn_rank``    -- fixed-width window in rank space.
* ``graph_1hop``  -- the gene's **network** neighbours (1-hop) inside the fold.
* ``graph_2hop``  -- 2-hop network neighbourhood.
* ``graph_denoise`` -- the gene's network neighbours, the neighbourhood median
  computed on a *denoised* score vector (median filter over the score
  neighbourhood) before averaging.
* ``global``      -- control C2 (shrink toward the fold's global median).
* ``random``      -- control C1 (matched-strength random neighbourhood).

Usage:
    python scripts/experiment_operators.py --only-net funmap --only-label intogen2024
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from d1.engine.stability import (  # noqa: E402
    crowding_counts, knn_positions, percentile_ranks, weighted_median,
)
from d1.networks.io import load_network  # noqa: E402

SCORE_DIR = REPO_ROOT / "results" / "scores_v40"
SPLIT_DIR = REPO_ROOT / "data" / "processed" / "splits_v40"
NET_DIR = REPO_ROOT / "data" / "processed" / "networks"


def _adj_list(A):
    A = csr_matrix(A)
    return [A.indices[A.indptr[i]:A.indptr[i + 1]] for i in range(A.shape[0])]


def _graph_neighbour_median(scores, adj, rows_idx, conf, order_hint=None):
    """Confidence-weighted median of the 1-hop network neighbours inside the fold."""
    n = len(rows_idx)
    vals = np.full((n, 0), np.nan)
    # build the neighbour lists restricted to the fold
    pos_of = {g: i for i, g in enumerate(rows_idx)}
    neigh = []
    for i, g in enumerate(rows_idx):
        nb = [pos_of[j] for j in adj[g] if j in pos_of]
        neigh.append(np.asarray(nb, dtype=int))
    max_nb = max((len(x) for x in neigh), default=0)
    if max_nb == 0:
        return scores.copy()
    vals = np.full((n, max_nb), np.nan)
    wts = np.full((n, max_nb), np.nan)
    for i, nb in enumerate(neigh):
        if nb.size:
            vals[i, :nb.size] = scores[nb]
            wts[i, :nb.size] = conf[nb]
    # weighted median with missing entries: replace NaN weights by 0
    valid = ~np.isnan(wts)
    w = np.where(valid, wts, 0.0)
    v = np.where(valid, vals, 0.0)
    return weighted_median(v, w)


def _rank_window_median(scores, genes, w):
    order = np.lexsort((np.asarray(genes).astype(str), scores))
    sv = scores[order]
    n = len(sv)
    idx = np.arange(n)[:, None] + np.arange(-w, w + 1)[None, :]
    idx = np.clip(idx, 0, n - 1)
    med = np.median(sv[idx], axis=1)
    out = np.empty(n)
    out[order] = med
    return out


def run(net_id, label_id, alpha=0.5, lams=(0.25, 0.5, 0.75, 1.0), k=20, w=20):
    d = np.load(SCORE_DIR / f"{net_id}__{label_id}__alpha{alpha}.npz", allow_pickle=True)
    P, genes = d["P"], np.asarray(d["genes"], dtype=str)
    y_all, degree = d["y"], d["degree"]
    fold_repeat = d["fold_repeat"]
    n_genes, n_folds = P.shape

    sp = pd.read_csv(SPLIT_DIR / f"{label_id}__native_{net_id}.tsv", sep="\t",
                     dtype={"gene": str})
    gpos = pd.Series(np.arange(n_genes), index=genes)
    sp = sp.assign(idx=gpos.loc[sp.gene].to_numpy()).sort_values(["repeat", "fold"])
    cols = list(sp.groupby(["repeat", "fold"]).groups.keys())
    test_idx = [sp[(sp.repeat == r) & (sp.fold == f)].idx.to_numpy() for (r, f) in cols]

    net_genes, A = load_network(str(NET_DIR / f"{net_id}.tsv"))
    adj = _adj_list(A)

    # per-repeat score matrix -> confidence (leave-one-repeat-out)
    S = np.full((n_folds, n_genes), np.nan)
    for j, te in enumerate(test_idx):
        S[j, te] = P[te, j]
    from d1.engine.stability import leave_one_repeat_out_confidence, score_matrix_by_repeat
    conf = leave_one_repeat_out_confidence(score_matrix_by_repeat(S, fold_repeat))

    acc = {}
    def new_acc():
        return {"sum": np.zeros(n_genes), "sq": np.zeros(n_genes), "cnt": np.zeros(n_genes)}

    def add(a, te, p):
        a["sum"][te] += p
        a["sq"][te] += p * p
        a["cnt"][te] += 1

    def sd(a):
        c = a["cnt"]
        with np.errstate(invalid="ignore", divide="ignore"):
            var = (a["sq"] - a["sum"] ** 2 / c) / (c - 1)
        return np.sqrt(np.maximum(var, 0.0))

    methods = ["baseline", "knn_score", "knn_rank", "graph_1hop", "global", "random"]
    accs = {(m, lam): new_acc() for m in methods for lam in lams}
    accs[("baseline", 0.0)] = new_acc()
    d_auc = {k2: [] for k2 in accs}

    rng = np.random.default_rng(7)
    for j, te in enumerate(test_idx):
        s = P[te, j]
        g = genes[te]
        rep = fold_repeat[j]
        c = conf[rep][te]
        order = np.lexsort((g, s))
        sv = s[order]
        pos = knn_positions(sv, k)
        vals = np.concatenate([sv[:, None], sv[pos]], axis=-1)
        wts = np.concatenate([c[order][:, None], c[order][pos]], axis=-1)
        med_sorted = weighted_median(vals, wts)
        med_knn = np.empty(len(te))
        med_knn[order] = med_sorted
        med_rank = _rank_window_median(s, g, w)
        med_graph = _graph_neighbour_median(s, adj, te, c)
        gmed = float(np.median(s))
        rnd = rng.integers(0, len(te), size=(len(te), k))
        rv = s[rnd]
        med_rand = np.median(rv, axis=1)

        add(accs[("baseline", 0.0)], te, percentile_ranks(s, g))
        d_auc[("baseline", 0.0)].append(roc_auc_score(y_all[te], s))
        for lam in lams:
            for name, med in (("knn_score", med_knn), ("knn_rank", med_rank),
                              ("graph_1hop", med_graph), ("global", np.full(len(te), gmed)),
                              ("random", med_rand)):
                sp_ = s + lam * (med - s)
                add(accs[(name, lam)], te, percentile_ranks(sp_, g))
                d_auc[(name, lam)].append(roc_auc_score(y_all[te], sp_))

    rows = []
    for key, a in accs.items():
        rows.append({"method": key[0], "lam": key[1],
                     "median_rank_sd": float(np.nanmedian(sd(a))),
                     "auroc": float(np.mean(d_auc[key]))})
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    args = ap.parse_args(argv)
    nets = args.only_net or ["funmap", "reactome"]
    labels = args.only_label or ["intogen2024"]
    out = []
    for label in labels:
        for net in nets:
            df = run(net, label, args.alpha)
            df["network"], df["label_set"] = net, label
            out.append(df)
            print(f"\n=== {net} x {label} ===")
            print(df.pivot_table(index="method", columns="lam", values="median_rank_sd")
                  .round(5).to_string())
            print(df.pivot_table(index="method", columns="lam", values="auroc")
                  .round(4).to_string())
    pd.concat(out, ignore_index=True).to_csv(
        REPO_ROOT / "scratch" / "experiment_operators.tsv", sep="\t", index=False, lineterminator="\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
