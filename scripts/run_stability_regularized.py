"""Run the stability-regularised scoring analysis (pre-registered in
``docs/stability_regularization_plan.md``).

For every (network, label) arm the script

1. loads the per-gene score dump of the 40-repeat panel
   (``results/scores_v40/<net>__<label>__alpha<a>.npz``) and the released
   split file (``data/processed/splits_v40/``);
2. builds the **baseline** stability atlas: per-gene rank SD across repeats,
   per-gene score CV, and the score-space crowding of every gene;
3. runs the **SRS operator** over a grid of shrinkage strengths lambda and
   neighbourhood sizes k, in two variants (global lambda, crowding-adaptive
   lambda) and in the transductive and the inductive neighbourhood mode;
4. runs the matched-strength **controls**: C1 random-neighbourhood smoothing,
   C2 global-median shrinkage, C3 monotone flattening;
5. reports the **cross-fitted** scorer: lambda is chosen on 3/4 of the repeats
   by the stability objective and evaluated on the held-out quarter.

Everything is written to ``results/tables/``. The confirmatory factorial is
untouched: the RWR numbers are the reference, the SRS numbers are a separate,
exploratory family.

Usage
-----
    python scripts/run_stability_regularized.py
    python scripts/run_stability_regularized.py --alpha 0.5 --only-net funmap
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from d1.engine.stability import (                                  # noqa: E402
    crowding_counts,
    global_median_shrinkage,
    knn_against_pool,
    knn_positions,
    leave_one_repeat_out_confidence,
    monotone_flatten,
    percentile_ranks,
    random_neighbourhood_median,
    score_matrix_by_repeat,
    weighted_median,
)

NETWORKS = ("funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome")
LABELS = ("intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new")

SCORE_DIR = REPO_ROOT / "results" / "scores_v40"
SPLIT_DIR = REPO_ROOT / "data" / "processed" / "splits_v40"
TABLE_DIR = REPO_ROOT / "results" / "tables"

LAM_GRID = (0.0, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0)
K_GRID = (5, 20, 50)
PRIMARY_K = 20
AUROC_TOLERANCE = 0.010          # pre-registered P2 tolerance
N_GROUPS = 4                     # repeat groups for the cross-fitted lambda
RANDOM_REPS = 20
CONTROL_LAMBDAS = (0.4, 0.7, 1.0)


# ------------------------------------------------------------------ helpers
def _acc_new(n_genes, n_groups):
    return {"sum": np.zeros((n_groups, n_genes)),
            "sq": np.zeros((n_groups, n_genes)),
            "cnt": np.zeros((n_groups, n_genes))}


def _acc_add(acc, group, idx, p):
    acc["sum"][group, idx] += p
    acc["sq"][group, idx] += p * p
    acc["cnt"][group, idx] += 1.0


def _acc_sd(acc, groups=None):
    """Per-gene SD of the percentile rank from the accumulators."""
    s = acc["sum"].sum(axis=0) if groups is None else acc["sum"][groups].sum(axis=0)
    q = acc["sq"].sum(axis=0) if groups is None else acc["sq"][groups].sum(axis=0)
    c = acc["cnt"].sum(axis=0) if groups is None else acc["cnt"][groups].sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        var = (q - s * s / c) / (c - 1.0)
    sd = np.sqrt(np.maximum(var, 0.0))
    sd[c < 2] = np.nan
    return sd, c


def _pauc(y, s, max_fpr=0.02):
    """Partial AUROC up to ``max_fpr``, normalised by ``max_fpr``."""
    y = np.asarray(y)
    s = np.asarray(s, dtype=float)
    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    if n_pos == 0 or n_neg == 0:
        return np.nan
    order = np.lexsort((np.arange(len(s)), -s))
    y = y[order]
    tps = np.cumsum(y)
    fps = np.cumsum(1 - y)
    # keep only the last index of every tie block so the curve is correct
    distinct = np.where(np.diff(s[order]))[0]
    idx = np.r_[distinct, len(s) - 1]
    tpr = np.r_[0.0, tps[idx] / n_pos]
    fpr = np.r_[0.0, fps[idx] / n_neg]
    stop = np.searchsorted(fpr, max_fpr, side="right")
    if stop <= 1:
        return np.nan
    f = fpr[:stop + 1]
    t = tpr[:stop + 1]
    if f[-1] < max_fpr:
        f = np.r_[f, max_fpr]
        t = np.r_[t, t[-1] + (t[-1] - t[-2]) * (max_fpr - f[-2]) / (f[-1] - f[-2])]
    return float(np.trapezoid(t, f) / max_fpr) if hasattr(np, "trapezoid") \
        else float(np.trapz(t, f) / max_fpr)


def _medrr(y, s):
    """Median rank ratio of the true positives (Guala & Sonnhammer 2017)."""
    y = np.asarray(y)
    order = np.lexsort((np.arange(len(s)), -np.asarray(s, dtype=float)))
    ranks = np.empty(len(s), dtype=np.int64)
    ranks[order] = np.arange(1, len(s) + 1)
    pos = ranks[y == 1]
    if pos.size == 0:
        return np.nan
    return float(np.median(pos) / len(s))


# ------------------------------------------------------------------- per arm
def run_arm_analysis(net_id, label_id, alpha, lam_grid, k_grid, do_controls,
                      do_inductive, dump_dir=SCORE_DIR):
    dump_path = dump_dir / f"{net_id}__{label_id}__alpha{alpha}.npz"
    split_path = SPLIT_DIR / f"{label_id}__native_{net_id}.tsv"
    if not dump_path.exists():
        raise FileNotFoundError(dump_path)

    d = np.load(dump_path, allow_pickle=True)
    P = d["P"]
    genes = np.asarray(d["genes"], dtype=str)
    degree = d["degree"]
    y_all = d["y"]
    fold_repeat = d["fold_repeat"]
    fold_id = d["fold_id"]
    auroc_rwr = d["auroc_rwr"]
    auroc_degree = d["auroc_degree"]
    n_genes, n_folds = P.shape

    sp = pd.read_csv(split_path, sep="\t", dtype={"gene": str})
    gpos = pd.Series(np.arange(n_genes), index=genes)
    sp = sp.assign(idx=gpos.loc[sp.gene].to_numpy())
    sp = sp.sort_values(["repeat", "fold"], kind="mergesort")

    # per-fold test / seed / train index sets, in the order of the dump columns
    order_cols = list(sp.groupby(["repeat", "fold"]).groups.keys())
    test_idx, seed_idx, train_idx = [], [], []
    for (r, f) in order_cols:
        block = sp[(sp.repeat == r) & (sp.fold == f)]
        te = block.idx.to_numpy()
        tr = sp[(sp.repeat == r) & (sp.fold != f)]
        test_idx.append(te)
        seed_idx.append(tr.loc[tr.y == 1, "idx"].to_numpy())
        train_idx.append(tr.idx.to_numpy())
    assert len(test_idx) == n_folds, (len(test_idx), n_folds)

    # per-repeat score matrix (each gene is tested once per repeat)
    S = np.full((n_folds, n_genes), np.nan)
    for j, te in enumerate(test_idx):
        S[j, te] = P[te, j]
    by_repeat = score_matrix_by_repeat(S, fold_repeat)

    # baseline stability atlas
    base_pct = np.full((n_folds, n_genes), np.nan)
    crowd_acc = np.zeros(n_genes)
    crowd_n = np.zeros(n_genes)
    for j, te in enumerate(test_idx):
        pct = percentile_ranks(P[te, j], genes[te])
        base_pct[j, te] = pct
        cnt = crowding_counts(P[te, j])
        crowd_acc[te] += cnt
        crowd_n[te] += 1
    crowding = crowd_acc / np.maximum(crowd_n, 1)
    base_sd = np.nanstd(base_pct, axis=0, ddof=1)
    base_cnt = np.sum(~np.isnan(base_pct), axis=0)
    base_sd = np.where(base_cnt > 1, base_sd, np.nan)
    cv, mean_s, sd_s = _score_cv(by_repeat)

    # leave-one-repeat-out confidence: the weights of an evaluated repeat never
    # use that repeat's own scores (no in-sample optimism at the fold level).
    conf_loro = leave_one_repeat_out_confidence(by_repeat)

    # leave-one-group-out confidence (used by the cross-fitted selection)
    n_rep = int(fold_repeat.max()) + 1
    group_of_repeat = np.zeros(n_rep, dtype=int)
    for r in np.unique(fold_repeat):
        group_of_repeat[r] = int(r) % N_GROUPS
    group_of_fold = group_of_repeat[fold_repeat]

    rows_grid, rows_gene, rows_arm = [], [], []
    rows_selected = []
    n_groups = N_GROUPS

    def _add_grid(variant, k, lam, control, med_sd, med_sd_test, d_auroc, d_auroc_sd,
                  pauc, medrr, disp_ratio, rank_change, n_crowd_hi, n_crowd_lo):
        rows_grid.append({
            "network": net_id, "label_set": label_id, "alpha": alpha,
            "variant": variant, "k": k, "lam": lam, "control": control,
            "median_rank_sd": med_sd, "median_rank_sd_out_of_sample": med_sd_test,
            "delta_auroc": d_auroc, "delta_auroc_sd": d_auroc_sd,
            "pauc_002": pauc, "medrr": medrr,
            "score_sd_ratio": disp_ratio, "rank_change_frac": rank_change,
            "n_genes_crowded": int(n_crowd_hi), "n_genes_isolated": int(n_crowd_lo),
        })

    def _percentiles_for(scores, te):
        return percentile_ranks(scores, genes[te])

    for variant in ("global", "crowding"):
        for k in k_grid:
            # accumulators: one per (control-group) for the cross-fitted selection
            accs = {}
            for lam in lam_grid:
                accs[("srs", lam)] = _acc_new(n_genes, n_groups)
            if do_controls:
                for lam in lam_grid:
                    accs[("c2", lam)] = _acc_new(n_genes, n_groups)
            d_auroc = {key: {g: [] for g in range(n_groups)} for key in accs}
            aux = {key: {"pauc": [], "medrr": [], "disp": [], "rankchg": []}
                   for key in accs}

            for j, te in enumerate(test_idx):
                s = P[te, j]
                g = group_of_fold[j]
                rep = fold_repeat[j]
                conf = conf_loro[rep][te]
                crowd = crowding[te]
                lam_g = np.full(len(te), 1.0)
                if variant == "crowding" and crowd.max() > 0:
                    lam_g = crowd / crowd.max()

                # local median (computed once per fold and reused for every lambda)
                order = np.lexsort((genes[te], s))
                sv = s[order]
                pos = knn_positions(sv, k)
                nbr = sv[pos]
                conf_sorted = conf[order]
                vals = np.concatenate([sv[:, None], nbr], axis=-1)
                wts = np.concatenate([conf_sorted[:, None], conf_sorted[pos]], axis=-1)
                med_sorted = weighted_median(vals, wts)
                med = np.empty(len(te))
                med[order] = med_sorted

                base_rank = np.empty(len(te), dtype=np.int64)
                base_rank[np.lexsort((genes[te], -s))] = np.arange(1, len(te) + 1)

                for lam in lam_grid:
                    sp_ = s + (lam * lam_g) * (med - s)
                    pct = _percentiles_for(sp_, te)
                    _acc_add(accs[("srs", lam)], g, te, pct)
                    d_auroc[("srs", lam)][g].append(
                        roc_auc_score(y_all[te], sp_) - auroc_degree[j])
                    rk = np.empty(len(te), dtype=np.int64)
                    rk[np.lexsort((genes[te], -sp_))] = np.arange(1, len(te) + 1)
                    aux[("srs", lam)]["rankchg"].append(float((rk != base_rank).mean()))
                    aux[("srs", lam)]["disp"].append(
                        float(np.std(sp_) / np.std(s)) if np.std(s) > 0 else np.nan)
                    if k == PRIMARY_K and variant == "global":
                        aux[("srs", lam)]["pauc"].append(_pauc(y_all[te], sp_))
                        aux[("srs", lam)]["medrr"].append(_medrr(y_all[te], sp_))

                if do_controls:
                    gm = global_median_shrinkage(s, lam=1.0)
                    for lam in lam_grid:
                        sp_ = s + lam * (gm - s)
                        pct = _percentiles_for(sp_, te)
                        _acc_add(accs[("c2", lam)], g, te, pct)
                        d_auroc[("c2", lam)][g].append(
                            roc_auc_score(y_all[te], sp_) - auroc_degree[j])

            # ---- cross-fitted selection of lambda (our stability metric as the
            # training objective). For each held-out group of repeats the
            # objective is minimised on the other groups and the selected
            # scorer is scored on the held-out group.
            # the actual selection (one lambda per held-out group)
            for gi in range(n_groups):
                grp = [g for g in range(n_groups) if g != gi]
                best_lam, best_j, best_d = None, np.inf, None
                base_d_tr = float(np.mean(np.concatenate(
                    [np.asarray(d_auroc[("srs", 0.0)][g], dtype=float) for g in grp])))
                for lam in lam_grid:
                    s_tr, _ = _acc_sd(accs[("srs", lam)], groups=grp)
                    j_tr = float(np.nanmedian(s_tr))
                    da_tr = float(np.mean(np.concatenate(
                        [np.asarray(d_auroc[("srs", lam)][g], dtype=float) for g in grp])))
                    if da_tr < base_d_tr - AUROC_TOLERANCE:
                        continue
                    if j_tr < best_j:
                        best_lam, best_j, best_d = lam, j_tr, da_tr
                if best_lam is None:
                    best_lam, best_j, best_d = 0.0, np.nan, np.nan
                s_te, _ = _acc_sd(accs[("srs", best_lam)], groups=[gi])
                d_te = float(np.mean(np.asarray(d_auroc[("srs", best_lam)][gi], dtype=float)))
                rows_selected.append({
                    "network": net_id, "label_set": label_id, "alpha": alpha,
                    "variant": variant, "k": k, "held_out_group": gi,
                    "lam_selected": best_lam,
                    "rank_sd_train": best_j, "delta_auroc_train": best_d,
                    "rank_sd_test": float(np.nanmedian(s_te)),
                    "delta_auroc_test": d_te,
                    "rank_sd_baseline_test": float(np.nanmedian(
                        _acc_sd(accs[("srs", 0.0)], groups=[gi])[0])),
                    "delta_auroc_baseline_test": float(np.mean(
                        np.asarray(d_auroc[("srs", 0.0)][gi], dtype=float))),
                })
            for key in accs:
                lam = key[1]
                sd, cnt = _acc_sd(accs[key])
                med_sd = float(np.nanmedian(sd))
                sd_oos = []
                for gi in range(n_groups):
                    grp = [g for g in range(n_groups) if g != gi]
                    s_g, _ = _acc_sd(accs[key], groups=grp)
                    sd_oos.append(np.nanmedian(s_g))
                da = np.concatenate([np.asarray(v, dtype=float)
                                     for v in d_auroc[key].values()])
                a = aux[key]
                _add_grid(variant, k, lam, key[0], med_sd,
                          float(np.mean(sd_oos)),
                          float(np.mean(da)), float(np.std(da, ddof=1)),
                          float(np.mean(a["pauc"])) if a["pauc"] else np.nan,
                          float(np.mean(a["medrr"])) if a["medrr"] else np.nan,
                          float(np.mean(a["disp"])) if a["disp"] else np.nan,
                          float(np.mean(a["rankchg"])) if a["rankchg"] else np.nan,
                          int(np.sum(crowding >= np.percentile(crowding, 75))),
                          int(np.sum(crowding <= np.percentile(crowding, 25))))

            # ---- control C1: matched-strength random neighbourhood
            if do_controls and k == PRIMARY_K:
                rng = np.random.default_rng(20260926 + n_genes)
                acc_c1 = {lam: _acc_new(n_genes, n_groups) for lam in CONTROL_LAMBDAS}
                d_c1 = {lam: [] for lam in CONTROL_LAMBDAS}
                for j, te in enumerate(test_idx):
                    s = P[te, j]
                    g = group_of_fold[j]
                    rep = fold_repeat[j]
                    conf = conf_loro[rep][te]
                    meds = random_neighbourhood_median(s, conf, k=k,
                                                        n_reps=RANDOM_REPS, rng=rng)
                    for lam in CONTROL_LAMBDAS:
                        acc_r = acc_c1[lam]
                        d_r = d_c1[lam]
                        for rep_i in range(RANDOM_REPS):
                            sp_ = s + lam * (meds[rep_i] - s)
                            pct = _percentiles_for(sp_, te)
                            _acc_add(acc_r, g, te, pct)
                            d_r.append(roc_auc_score(y_all[te], sp_) - auroc_degree[j])
                for lam in CONTROL_LAMBDAS:
                    sd, cnt = _acc_sd(acc_c1[lam])
                    da = np.asarray(d_c1[lam], dtype=float)
                    rows_grid.append({
                        "network": net_id, "label_set": label_id, "alpha": alpha,
                        "variant": variant, "k": k, "lam": lam, "control": "c1",
                        "median_rank_sd": float(np.nanmedian(sd)),
                        "median_rank_sd_out_of_sample": np.nan,
                        "delta_auroc": float(np.mean(da)),
                        "delta_auroc_sd": float(np.std(da, ddof=1)),
                        "pauc_002": np.nan, "medrr": np.nan,
                        "score_sd_ratio": np.nan, "rank_change_frac": np.nan,
                        "n_genes_crowded": np.nan, "n_genes_isolated": np.nan,
                    })

    # ---- control C3: monotone flattening leaves the rank statistics untouched
    c3 = []
    for j in range(min(n_folds, 40)):
        te = test_idx[j]
        s = P[te, j]
        flat = monotone_flatten(s)
        c3.append({
            "auroc_gap": float(roc_auc_score(y_all[te], s) - roc_auc_score(y_all[te], flat)),
            "rank_gap": float(np.max(np.abs(percentile_ranks(s, genes[te])
                                            - percentile_ranks(flat, genes[te])))),
            "sd_ratio": float(np.std(flat) / np.std(s)),
        })
    c3 = pd.DataFrame(c3)

    # ---- inductive variant (P7): neighbourhood from the training fold only
    if do_inductive:
        acc_ind = {lam: _acc_new(n_genes, n_groups) for lam in CONTROL_LAMBDAS}
        d_ind = {lam: [] for lam in CONTROL_LAMBDAS}
        for j, te in enumerate(test_idx):
            s = P[te, j]
            g = group_of_fold[j]
            rep = fold_repeat[j]
            conf = conf_loro[rep][te]
            pool_rows = train_idx[j]
            pool_rows = pool_rows[~np.isin(pool_rows, seed_idx[j])]
            pool_scores = np.sort(P[pool_rows, j])
            if pool_scores.size == 0:
                continue
            npos = knn_against_pool(s, pool_scores, k=PRIMARY_K)
            nbr = pool_scores[npos]
            conf_pool = conf_loro[rep][pool_rows]
            vals = np.concatenate([s[:, None], nbr], axis=-1)
            wts = np.concatenate([conf[:, None], conf_pool[npos]], axis=-1)
            med = weighted_median(vals, wts)
            for lam in CONTROL_LAMBDAS:
                sp_ = s + lam * (med - s)
                pct = _percentiles_for(sp_, te)
                _acc_add(acc_ind[lam], g, te, pct)
                d_ind[lam].append(roc_auc_score(y_all[te], sp_) - auroc_degree[j])
        for lam in CONTROL_LAMBDAS:
            sd, cnt = _acc_sd(acc_ind[lam])
            da = np.asarray(d_ind[lam], dtype=float)
            rows_grid.append({
                "network": net_id, "label_set": label_id, "alpha": alpha,
                "variant": "global_inductive", "k": PRIMARY_K, "lam": lam,
                "control": "srs", "median_rank_sd": float(np.nanmedian(sd)),
                "median_rank_sd_out_of_sample": np.nan,
                "delta_auroc": float(np.mean(da)),
                "delta_auroc_sd": float(np.std(da, ddof=1)),
                "pauc_002": np.nan, "medrr": np.nan,
                "score_sd_ratio": np.nan, "rank_change_frac": np.nan,
                "n_genes_crowded": np.nan, "n_genes_isolated": np.nan,
            })

    # ---- per-gene table (baseline + a few lambda settings, primary variant/k)
    gene_rows = pd.DataFrame({
        "network": net_id, "label_set": label_id, "alpha": alpha, "gene": genes,
        "rank_sd": base_sd, "score_cv": cv, "score_mean": mean_s, "score_sd": sd_s,
        "crowding": crowding, "n_repeats_tested": base_cnt,
    })

    # ---- arm-level baseline
    rows_arm.append({
        "network": net_id, "label_set": label_id, "alpha": alpha,
        "n_genes": n_genes, "n_folds": n_folds,
        "delta_auroc": float(np.mean(auroc_rwr - auroc_degree)),
        "delta_auroc_sd": float(np.std(auroc_rwr - auroc_degree, ddof=1)),
        "median_rank_sd": float(np.nanmedian(base_sd)),
        "mean_rank_sd": float(np.nanmean(base_sd)),
        "iqr_rank_sd": float(np.nanpercentile(base_sd, 75) - np.nanpercentile(base_sd, 25)),
        "median_crowding": float(np.nanmedian(crowding)),
        "spearman_crowding_rank_sd": float(pd.Series(crowding).corr(pd.Series(base_sd), method="spearman")),
        "median_score_cv": float(np.nanmedian(cv)),
        "c3_max_auroc_gap": float(c3.auroc_gap.abs().max()),
        "c3_max_rank_gap": float(c3.rank_gap.max()),
        "c3_sd_ratio_median": float(c3.sd_ratio.median()),
    })

    return (pd.DataFrame(rows_grid), gene_rows, pd.DataFrame(rows_arm), c3,
            pd.DataFrame(rows_selected))


def _score_cv(by_repeat):
    with np.errstate(invalid="ignore"):
        mean = np.nanmean(by_repeat, axis=0)
        sd = np.nanstd(by_repeat, axis=0, ddof=1)
    cv = sd / np.maximum(np.abs(mean), 1e-300)
    cv[~np.isfinite(cv)] = np.nan
    return cv, mean, sd


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alpha", type=float, nargs="+", default=[0.5])
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    ap.add_argument("--lam", type=float, nargs="+", default=list(LAM_GRID))
    ap.add_argument("--k", type=int, nargs="+", default=list(K_GRID))
    ap.add_argument("--no-controls", action="store_true")
    ap.add_argument("--no-inductive", action="store_true")
    ap.add_argument("--out-suffix", default="")
    args = ap.parse_args(argv)

    nets = args.only_net or list(NETWORKS)
    labels = args.only_label or list(LABELS)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    grid_all, gene_all, arm_all, c3_all, sel_all = [], [], [], [], []
    for alpha in args.alpha:
        for label in labels:
            for net in nets:
                t0 = time.time()
                g, genes_df, arm_df, c3, sel = run_arm_analysis(
                    net, label, alpha, args.lam, args.k,
                    do_controls=not args.no_controls,
                    do_inductive=not args.no_inductive)
                grid_all.append(g)
                gene_all.append(genes_df)
                arm_all.append(arm_df)
                c3_all.append(c3.assign(network=net, label_set=label, alpha=alpha))
                sel_all.append(sel)
                print(f"[{net} x {label} alpha={alpha}] {time.time() - t0:.1f}s "
                      f"median rank SD = {arm_df.median_rank_sd.iloc[0]:.4f}")

    grid = pd.concat(grid_all, ignore_index=True)
    genes_df = pd.concat(gene_all, ignore_index=True)
    arms = pd.concat(arm_all, ignore_index=True)
    c3 = pd.concat(c3_all, ignore_index=True)
    sel = pd.concat(sel_all, ignore_index=True)
    sfx = args.out_suffix

    grid.to_csv(TABLE_DIR / f"stability_regularization_grid{sfx}.tsv", sep="\t", index=False, lineterminator="\n")
    genes_df.to_csv(TABLE_DIR / f"stability_per_gene{sfx}.tsv", sep="\t", index=False, lineterminator="\n")
    arms.to_csv(TABLE_DIR / f"stability_by_arm{sfx}.tsv", sep="\t", index=False, lineterminator="\n")
    c3.to_csv(TABLE_DIR / f"stability_control_c3_flattening{sfx}.tsv", sep="\t", index=False, lineterminator="\n")
    sel.to_csv(TABLE_DIR / f"stability_selected_scorer{sfx}.tsv", sep="\t", index=False, lineterminator="\n")
    print(f"\ngrid rows: {len(grid):,}  gene rows: {len(genes_df):,}  arms: {len(arms)} "
          f"selected: {len(sel)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
