"""Which per-gene property predicts the ranking instability of a network prior?

The pre-registered mechanism was "crowding causes rank instability" (S2). This
script tests it against the alternatives, per arm, and writes the attribution
table used by the stability figure and the manuscript section.

Predictors (all per gene, all computed on the *unsmoothed* prior):

* ``crowding``        -- number of genes within +/- Scott's bandwidth of the
  gene's score (score-space crowding);
* ``crowding_rank``   -- number of genes within +/- 1 % of ranks;
* ``score_cv``        -- coefficient of variation of the gene's score across
  the 40 repeats (the gene's own score noise);
* ``degree``          -- node degree in the network.

Target: the per-gene SD of the within-fold percentile rank across repeats.

Usage:
    python scripts/run_stability_attribution.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from d1.engine.stability import percentile_ranks  # noqa: E402

SCORE_DIR = REPO_ROOT / "results" / "scores_v40"
SPLIT_DIR = REPO_ROOT / "data" / "processed" / "splits_v40"
TABLE_DIR = REPO_ROOT / "results" / "tables"

NETWORKS = ("funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome")
LABELS = ("intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new")
PREDICTORS = ("crowding", "crowding_rank", "score_cv", "degree")


def rank_space_crowding(pct, genes):
    """Number of genes within +/- 1 % of ranks, for one fold."""
    order = np.lexsort((np.asarray(genes).astype(str), -pct))
    rs = np.empty(len(pct), dtype=int)
    rs[order] = np.arange(len(pct))
    w = max(1, int(0.01 * len(pct)))
    return (np.abs(rs[:, None] - rs[None, :]) <= w).sum(axis=1) - 1


def arm_attribution(net_id, label_id, alpha=0.5):
    d = np.load(SCORE_DIR / f"{net_id}__{label_id}__alpha{alpha}.npz", allow_pickle=True)
    P, genes = d["P"], np.asarray(d["genes"], dtype=str)
    degree = d["degree"]
    fold_repeat = d["fold_repeat"]
    n_genes, n_folds = P.shape

    sp = pd.read_csv(SPLIT_DIR / f"{label_id}__native_{net_id}.tsv", sep="\t",
                     dtype={"gene": str})
    gpos = pd.Series(np.arange(n_genes), index=genes)
    sp = sp.assign(idx=gpos.loc[sp.gene].to_numpy()).sort_values(["repeat", "fold"])
    cols = list(sp.groupby(["repeat", "fold"]).groups.keys())
    test_idx = [sp[(sp.repeat == r) & (sp.fold == f)].idx.to_numpy() for (r, f) in cols]

    pct = np.full((n_folds, n_genes), np.nan)
    crowd_rank = np.zeros(n_genes)
    for j, te in enumerate(test_idx):
        p = percentile_ranks(P[te, j], genes[te])
        pct[j, te] = p
        crowd_rank[te] += rank_space_crowding(p, genes[te])
    crowd_rank /= n_folds

    rank_sd = np.nanstd(pct, axis=0, ddof=1)
    gene = pd.read_csv(TABLE_DIR / "stability_per_gene.tsv", sep="\t")
    gene = gene[(gene.network == net_id) & (gene.label_set == label_id)
                & (gene.alpha == alpha)].reset_index(drop=True)
    assert len(gene) == n_genes, (len(gene), n_genes)
    df = pd.DataFrame({
        "rank_sd": rank_sd,
        "crowding": gene.crowding.to_numpy(),
        "crowding_rank": crowd_rank,
        "score_cv": gene.score_cv.to_numpy(),
        "degree": degree,
    })
    out = {"network": net_id, "label_set": label_id, "alpha": alpha,
           "n_genes": n_genes}
    for p in PREDICTORS:
        out[f"rho_{p}"] = float(df.rank_sd.corr(df[p], method="spearman"))
        out[f"pearson_{p}"] = float(np.corrcoef(df.rank_sd, df[p])[0, 1])
    out["rho_rank_sd_vs_score_cv"] = out["rho_score_cv"]
    return df, out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    args = ap.parse_args(argv)
    nets = args.only_net or list(NETWORKS)
    labels = args.only_label or list(LABELS)
    rows = []
    for label in labels:
        for net in nets:
            _, out = arm_attribution(net, label, args.alpha)
            rows.append(out)
            print(f"[{net} x {label}] " + "  ".join(
                f"{k.replace('rho_',''):>13}={v:+.3f}"
                for k, v in out.items() if k.startswith("rho_")))
    res = pd.DataFrame(rows)
    path = TABLE_DIR / "stability_crowding_attribution.tsv"
    res.to_csv(path, sep="\t", index=False)
    print(f"\nWrote {path}")
    print(res[[c for c in res.columns if c.startswith("rho_")]].describe().round(3).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
