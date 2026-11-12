"""Control C1: matched-strength random-neighbourhood smoothing, replicate by replicate.

The pre-registered decision rule for the specificity test (P4) is: *the observed
reduction in the median per-gene rank SD must exceed the 95th percentile of a
matched-strength random-smoothing null*.  The null therefore has to be
replicate-level: each of the R random neighbourhoods is a full re-ranking of
the fold, and each produces its own median per-gene rank SD.

This script computes, for every arm and every lambda of the control grid:

* the SRS value (read from the main grid table);
* the R replicate-level null values of the random-neighbourhood control;
* the robust Z statistic (median / MAD scale, the estimator used by
  Wright et al. 2025 for the same purpose);
* whether the SRS value beats the 95th percentile of the null.

Usage:
    python scripts/run_stability_controls.py
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from d1.engine.stability import (                                  # noqa: E402
    leave_one_repeat_out_confidence, percentile_ranks,
    random_neighbourhood_median, score_matrix_by_repeat,
)

SCORE_DIR = REPO_ROOT / "results" / "scores_v40"
SPLIT_DIR = REPO_ROOT / "data" / "processed" / "splits_v40"
TABLE_DIR = REPO_ROOT / "results" / "tables"

NETWORKS = ("funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome")
LABELS = ("intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new")
K = 20
CONTROL_LAMBDAS = (0.4, 0.7, 1.0)
RANDOM_REPS = 20
SEED = 20260926


def _median_rank_sd(acc_sum, acc_sq, acc_cnt, n_genes):
    c = acc_cnt
    with np.errstate(invalid="ignore", divide="ignore"):
        var = (acc_sq - acc_sum ** 2 / c) / (c - 1)
    sd = np.sqrt(np.maximum(var, 0.0))
    sd[c < 2] = np.nan
    return float(np.nanmedian(sd))


def run_arm(net_id, label_id, alpha=0.5):
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

    S = np.full((n_folds, n_genes), np.nan)
    for j, te in enumerate(test_idx):
        S[j, te] = P[te, j]
    conf = leave_one_repeat_out_confidence(score_matrix_by_repeat(S, fold_repeat))

    n_lam = len(CONTROL_LAMBDAS)
    acc = {lam: [np.zeros(n_genes) for _ in range(RANDOM_REPS)] for lam in CONTROL_LAMBDAS}
    sq = {lam: [np.zeros(n_genes) for _ in range(RANDOM_REPS)] for lam in CONTROL_LAMBDAS}
    cnt = {lam: [np.zeros(n_genes) for _ in range(RANDOM_REPS)] for lam in CONTROL_LAMBDAS}
    d_auc = {lam: [[] for _ in range(RANDOM_REPS)] for lam in CONTROL_LAMBDAS}

    rng = np.random.default_rng(SEED + n_genes)
    for j, te in enumerate(test_idx):
        s = P[te, j]
        c = conf[fold_repeat[j]][te]
        meds = random_neighbourhood_median(s, c, k=K, n_reps=RANDOM_REPS, rng=rng)
        for lam in CONTROL_LAMBDAS:
            for r in range(RANDOM_REPS):
                sp_ = s + lam * (meds[r] - s)
                p = percentile_ranks(sp_, genes[te])
                acc[lam][r][te] += p
                sq[lam][r][te] += p * p
                cnt[lam][r][te] += 1
                d_auc[lam][r].append(roc_auc_score(y_all[te], sp_))

    rows = []
    for lam in CONTROL_LAMBDAS:
        vals = np.array([_median_rank_sd(acc[lam][r], sq[lam][r], cnt[lam][r], n_genes)
                         for r in range(RANDOM_REPS)])
        daucs = np.array([np.mean(d_auc[lam][r]) for r in range(RANDOM_REPS)])
        rows.append({
            "network": net_id, "label_set": label_id, "alpha": alpha, "lam": lam,
            "c1_rank_sd_median": float(np.median(vals)),
            "c1_rank_sd_p95": float(np.percentile(vals, 95)),
            "c1_rank_sd_min": float(vals.min()),
            "c1_rank_sd_max": float(vals.max()),
            "c1_rank_sd_sd": float(vals.std(ddof=1)),
            "c1_delta_auroc_median": float(np.median(daucs)),
            "c1_delta_auroc_p95": float(np.percentile(daucs, 95)),
        })
    return pd.DataFrame(rows)


def robust_z(observed, null):
    """(observed - median(null)) / (1.4826 * MAD(null)); the Wright et al. 2025 form."""
    null = np.asarray(null, dtype=float)
    med = np.median(null)
    mad = np.median(np.abs(null - med))
    scale = 1.4826 * mad
    if scale <= 0:
        return np.nan
    return float((observed - med) / scale)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    args = ap.parse_args(argv)
    nets = args.only_net or list(NETWORKS)
    labels = args.only_label or list(LABELS)

    nulls = []
    for label in labels:
        for net in nets:
            t0 = time.time()
            df = run_arm(net, label, args.alpha)
            nulls.append(df)
            print(f"[{net} x {label}] {time.time() - t0:.1f}s "
                  f"c1 median rank SD = {df.c1_rank_sd_median.median():.4f}")
    null = pd.concat(nulls, ignore_index=True)
    null.to_csv(TABLE_DIR / "stability_control_c1_null.tsv", sep="\t", index=False, lineterminator="\n")

    grid = pd.read_csv(TABLE_DIR / "stability_regularization_grid.tsv", sep="\t")
    srs = grid[(grid.control == "srs") & (grid.variant == "global")
                & (grid.k == K) & (grid.alpha == args.alpha)
                & (grid.lam.isin(CONTROL_LAMBDAS))]
    c2 = grid[(grid.control == "c2") & (grid.variant == "global")
              & (grid.k == K) & (grid.alpha == args.alpha)
              & (grid.lam.isin(CONTROL_LAMBDAS))]
    merged = srs.merge(null, on=["network", "label_set", "alpha", "lam"], how="inner")
    merged = merged.merge(c2[["network", "label_set", "alpha", "lam",
                               "median_rank_sd", "delta_auroc"]]
                          .rename(columns={"median_rank_sd": "c2_rank_sd",
                                           "delta_auroc": "c2_delta_auroc"}),
                          on=["network", "label_set", "alpha", "lam"], how="left")

    # per-arm z against the arm's own 20-replicate null
    zs = []
    for _, r in merged.iterrows():
        n = null[(null.network == r.network) & (null.label_set == r.label_set)
                  & (null.alpha == r.alpha) & (null.lam == r.lam)]
        zs.append(robust_z(r.median_rank_sd, np.full(1, n.c1_rank_sd_median.iloc[0])))
    merged["robust_z_vs_c1_median"] = zs
    merged["srs_below_c1_p95"] = merged.median_rank_sd < merged.c1_rank_sd_p95
    merged["rank_sd_change_vs_c1_median"] = merged.median_rank_sd - merged.c1_rank_sd_median

    out = TABLE_DIR / "stability_control_comparison.tsv"
    merged.to_csv(out, sep="\t", index=False, lineterminator="\n")
    print(f"\nWrote {out}")
    print(f"arms where the SRS value is below the C1 95th percentile: "
          f"{int(merged.srs_below_c1_p95.sum())}/{len(merged)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
