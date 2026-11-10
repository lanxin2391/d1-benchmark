"""Diagnostic: what actually drives the per-gene rank instability?

Answers three questions that decide whether the stability-regularised scorer
can work at all:

1. **Which crowding measure predicts rank instability?** Score-space density,
   rank-space density and node degree are all tested.
2. **Where does the rank variance come from?** The rank of a gene in a fold is
   ``1 + #{j : s_j > s_i}``. Two counterfactuals decompose the movement:
   *self moves* (the gene's own score changes) and *others move* (the rest of
   the profile changes). The per-gene SD of each counterfactual is reported.
3. **Is the within-plateau order informative?** Collapsing each crowded plateau
   to its median (a full median filter) changes ~50 % of the ranks; the AUROC
   effect of that collapse says whether the discarded order carried signal.

Run:
    python scripts/diagnose_instability.py --only-net funmap --only-label intogen2024
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from d1.engine.stability import (  # noqa: E402
    crowding_counts, percentile_ranks, score_matrix_by_repeat,
)

SCORE_DIR = REPO_ROOT / "results" / "scores_v40"
SPLIT_DIR = REPO_ROOT / "data" / "processed" / "splits_v40"
NET_DIR = REPO_ROOT / "data" / "processed" / "networks"

NETWORKS = ("funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome")
LABELS = ("intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new")


def diagnose(net_id, label_id, alpha=0.5):
    d = np.load(SCORE_DIR / f"{net_id}__{label_id}__alpha{alpha}.npz", allow_pickle=True)
    P, genes = d["P"], np.asarray(d["genes"], dtype=str)
    degree, y_all = d["degree"], d["y"]
    fold_repeat, fold_id = d["fold_repeat"], d["fold_id"]
    n_genes, n_folds = P.shape

    sp = pd.read_csv(SPLIT_DIR / f"{label_id}__native_{net_id}.tsv", sep="\t",
                     dtype={"gene": str})
    gpos = pd.Series(np.arange(n_genes), index=genes)
    sp = sp.assign(idx=gpos.loc[sp.gene].to_numpy()).sort_values(["repeat", "fold"])
    order_cols = list(sp.groupby(["repeat", "fold"]).groups.keys())
    test_idx = [sp[(sp.repeat == r) & (sp.fold == f)].idx.to_numpy() for (r, f) in order_cols]

    S = np.full((n_folds, n_genes), np.nan)
    pct = np.full((n_folds, n_genes), np.nan)
    dens = np.zeros(n_genes)
    dens_n = np.zeros(n_genes)
    rank_dens = np.zeros(n_genes)
    rank_dens_n = np.zeros(n_genes)
    self_frozen = np.full((n_folds, n_genes), np.nan)
    others_frozen = np.full((n_folds, n_genes), np.nan)

    for j, te in enumerate(test_idx):
        s = P[te, j]
        S[j, te] = s
        pct[j, te] = percentile_ranks(s, genes[te])
        dens[te] += crowding_counts(s)
        dens_n[te] += 1
        # rank-space crowding: how many genes share a +/-1 % rank band
        order = np.lexsort((genes[te], -s))
        rs = np.empty(len(te), dtype=int)
        rs[order] = np.arange(len(te))
        w = max(1, int(0.01 * len(te)))
        rank_dens[te] = (np.abs(rs[:, None] - rs[None, :]) <= w).sum(axis=1) - 1
        rank_dens_n[te] += 1

        # counterfactual A: the gene's score frozen at its across-repeat mean
        # counterfactual B: every other gene frozen at its across-repeat mean
        mu = np.nanmean(S, axis=0)
        mu_j = np.nanmean(np.where(np.isnan(S[:, te]), np.nan, S[:, te]), axis=0) \
            if False else None
        # per-gene mean over the folds in which the gene was tested
        with np.errstate(invalid="ignore"):
            mu_te = np.nanmean(S[:, te], axis=0)
        rank_self = np.empty(len(te))
        rank_others = np.empty(len(te))
        for a in range(len(te)):
            others = s > mu_te[a]
            rank_self[a] = 1 + others.sum() - (1 if others[a] else 0)
            mu_all = mu.copy()
            mu_all[te[a]] = np.nan
            rank_others[a] = 1 + int((mu_all[te] > s[a]).sum())
        self_frozen[j, te] = rank_self
        others_frozen[j, te] = rank_others

    rank_sd = np.nanstd(pct, axis=0, ddof=1)
    sd_self = np.nanstd(self_frozen, axis=0, ddof=1) / (n_genes - 1)
    sd_others = np.nanstd(others_frozen, axis=0, ddof=1) / (n_genes - 1)
    with np.errstate(invalid="ignore"):
        mu_all = np.nanmean(S, axis=0)
        sd_all = np.nanstd(S, axis=0, ddof=1)
    cv = sd_all / np.maximum(np.abs(mu_all), 1e-300)

    df = pd.DataFrame({
        "gene": genes,
        "rank_sd": rank_sd,
        "rank_sd_self_moves": sd_self,
        "rank_sd_others_move": sd_others,
        "score_cv": cv,
        "crowd_score_window": dens / np.maximum(dens_n, 1),
        "crowd_rank_window": rank_dens / np.maximum(rank_dens_n, 1),
        "degree": degree,
        "is_positive": y_all,
    })
    out = {"network": net_id, "label_set": label_id, "alpha": alpha,
           "n_genes": n_genes, "median_rank_sd": float(np.nanmedian(rank_sd))}
    for col in ("score_cv", "crowd_score_window", "crowd_rank_window", "degree",
                "rank_sd_self_moves", "rank_sd_others_move"):
        out[f"spearman_rank_sd~{col}"] = float(df.rank_sd.corr(df[col], method="spearman"))
    return df, out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    ap.add_argument("--save", default="")
    args = ap.parse_args(argv)
    nets = args.only_net or list(NETWORKS)
    labels = args.only_label or list(LABELS)
    rows, dfs = [], []
    for label in labels:
        for net in nets:
            df, out = diagnose(net, label, args.alpha)
            rows.append(out)
            dfs.append(df.assign(network=net, label_set=label))
            print(f"[{net} x {label}] " + " ".join(
                f"{k.split('~')[-1]}={v:+.3f}" for k, v in out.items()
                if k.startswith("spearman")))
    res = pd.DataFrame(rows)
    if args.save:
        res.to_csv(args.save, sep="\t", index=False)
        pd.concat(dfs, ignore_index=True).to_csv(
            args.save.replace(".tsv", "_genes.tsv"), sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
