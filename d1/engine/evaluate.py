"""Metrics and the evaluation harness (A1).

One "arm" = one network x one label set x one universe. For every repeat and fold of the split
file (contract C4):
    seeds = training-fold positives (y == 1 and fold != f)
    test  = all genes of fold f (positives and negatives); seeds are never in it (D-06)
    RWR   = propagate from the seeds on the whole network, score the test genes
    degree baseline = degree of the test genes (control 1)
All folds and repeats are propagated in ONE batched RWR call per alpha, which is why the whole
factorial runs in minutes.

Universe (decision D-12): split genes must be a subset of the network's genes. RWR always runs
on the full network (LCC); seeds and test genes come only from the split file. A `native`
universe = all network genes; the `shared` universe = genes present in all six networks.
"""
import os

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from d1.engine.baselines import degree_scores
from d1.engine.rwr import ALPHA, rwr, seed_matrix, transition_T
from d1.networks.io import load_network

KS = (50, 100, 500)
SPLIT_COLUMNS = ["gene", "repeat", "fold", "y"]
METRICS = ("auroc", "auprc") + tuple(f"p_at_{k}" for k in KS)
ARM_KEYS = ["network", "label_set", "universe"]


# ---------------------------------------------------------------- metrics
def precision_at_k(scores, y, genes, k):
    """Fraction of positives among the top-k genes. Ties: score desc, then symbol asc (D-07)."""
    y = np.asarray(y)
    if len(y) < k:
        return np.nan
    order = np.lexsort((np.asarray(genes).astype(str), -np.asarray(scores, dtype=float)))
    return float(y[order[:k]].mean())


def evaluate(scores, y, genes, ks=KS):
    """AUROC, AUPRC and precision-at-k for one ranked test set."""
    scores = np.asarray(scores, dtype=float)
    y = np.asarray(y).astype(int)
    out = {"n_test_pos": int(y.sum()), "n_test": int(len(y))}
    if 0 < y.sum() < len(y):
        out["auroc"] = float(roc_auc_score(y, scores))
        out["auprc"] = float(average_precision_score(y, scores))
    else:                                   # only one class: metrics undefined
        out["auroc"] = out["auprc"] = np.nan
    for k in ks:
        out[f"p_at_{k}"] = precision_at_k(scores, y, genes, k)
    return out


# ---------------------------------------------------------------- split files
def check_splits(df):
    """Raise ValueError if a split table breaks contract C4."""
    missing = set(SPLIT_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"split table is missing columns {missing} (contract C4)")
    if not set(pd.unique(df.y)) <= {0, 1}:
        raise ValueError("y must be 0 or 1")
    gene_sets = set()
    for r, d in df.groupby("repeat"):
        if d.gene.duplicated().any():
            raise ValueError(f"repeat {r}: a gene appears in more than one fold")
        gene_sets.add(frozenset(d.gene))
    if len(gene_sets) != 1:
        raise ValueError("repeats do not cover the same set of genes")
    if (df.groupby("gene").y.nunique() > 1).any():
        raise ValueError("a gene has different y in different repeats")
    if (df.groupby(["repeat", "fold"]).y.sum() == 0).any():
        raise ValueError("a fold has no positives")
    return True


def read_splits(path):
    df = pd.read_csv(path, sep="\t", dtype={"gene": str})
    check_splits(df)
    return df


# ---------------------------------------------------------------- harness
def run_arm(network, splits, network_id, label_id, universe_id, alphas=(ALPHA,), ks=KS):
    """Evaluate RWR and the degree baseline for one arm.

    network : path to a C2 file, or a (genes, A) tuple from load_network
    splits  : path to a C4 file, or a DataFrame with the C4 columns
    Returns a DataFrame in the C5 result format.
    """
    genes, A = load_network(network) if isinstance(network, (str, os.PathLike)) else network
    sp_df = read_splits(splits) if isinstance(splits, (str, os.PathLike)) else splits
    check_splits(sp_df)

    index = pd.Series(np.arange(len(genes)), index=genes)
    missing = sorted(set(sp_df.gene) - set(genes))
    if missing:
        raise ValueError(f"{len(missing)} split genes are not in network {network_id}, "
                         f"e.g. {missing[:5]}; build splits on a universe inside the network")
    sp_df = sp_df.assign(idx=index.loc[sp_df.gene].to_numpy())

    WT = transition_T(A)
    deg = degree_scores(A)

    combos, seeds, tests = [], [], []
    for (r, f) in sorted(sp_df.groupby(["repeat", "fold"]).groups):
        d = sp_df[sp_df.repeat == r]
        seeds.append(d.loc[(d.fold != f) & (d.y == 1), "idx"].to_numpy())
        tests.append(d[d.fold == f])
        combos.append((r, f))

    base = {"network": network_id, "label_set": label_id, "universe": universe_id}
    rows = []
    for (r, f), t in zip(combos, tests):
        ti = t.idx.to_numpy()
        rows.append({**base, "alpha": np.nan, "repeat": r, "fold": f, "method": "degree",
                     "n_iter": np.nan, **evaluate(deg[ti], t.y.to_numpy(), t.gene.to_numpy(), ks)})

    P0 = seed_matrix(len(genes), seeds)
    for a in alphas:
        P, n_iter = rwr(WT, P0, alpha=a)
        for j, ((r, f), t) in enumerate(zip(combos, tests)):
            ti = t.idx.to_numpy()
            rows.append({**base, "alpha": a, "repeat": r, "fold": f, "method": "rwr",
                         "n_iter": n_iter,
                         **evaluate(P[ti, j], t.y.to_numpy(), t.gene.to_numpy(), ks)})
    cols = ARM_KEYS + ["alpha", "repeat", "fold", "method", "n_iter", *METRICS,
                       "n_test_pos", "n_test"]
    return pd.DataFrame(rows)[cols]


def margins(results, metrics=METRICS):
    """Per-fold margins RWR - degree (the primary endpoint is delta_auroc)."""
    keys = ARM_KEYS + ["repeat", "fold"]
    metrics = list(metrics)
    deg = results[results.method == "degree"].set_index(keys)[metrics]
    out = []
    for a, d in results[results.method == "rwr"].groupby("alpha"):
        d = d.set_index(keys)[metrics]
        dg = deg.loc[d.index]
        m = (d - dg).add_prefix("delta_").join(d.add_prefix("rwr_")).join(dg.add_prefix("degree_"))
        m.insert(0, "alpha", a)
        out.append(m.reset_index())
    return pd.concat(out, ignore_index=True)


def summarize(results):
    """Mean and SD over folds x repeats of AUROC for RWR, degree and their margin."""
    m = margins(results)
    g = m.groupby(ARM_KEYS + ["alpha"])
    return pd.DataFrame({
        "rwr_auroc": g.rwr_auroc.mean(),
        "degree_auroc": g.degree_auroc.mean(),
        "delta_auroc": g.delta_auroc.mean(),
        "delta_auroc_sd": g.delta_auroc.std(),
        "delta_auprc": g.delta_auprc.mean(),
        "n_folds": g.size(),
    }).reset_index()
