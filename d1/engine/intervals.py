"""Confidence-interval constructions for the D1 benchmark.

Four fold/repeat-level constructions plus the gene-cluster bootstrap that
pre-registration section 9.3 asks for.

    fold        percentile bootstrap over the per-fold delta values
    repeat_pct  percentile bootstrap over the per-repeat mean deltas
    repeat_t    t interval on the per-repeat mean deltas (n-1 df)
    repeat_bca  bias-corrected and accelerated bootstrap on the repeat means
    gene        cluster bootstrap resampling GENES (the unit the
                pre-registration names), not gene-arm rows

The gene bootstrap needs the per-gene scores of every fold, which the C5
result table does not carry -- see `d1.engine.evaluate.dump_per_gene_scores`
and work-plan v3.1 task A3.

Why the unit matters (work-plan v3.1 task A1)
----------------------------------------------
The 50 (or 200) folds are 10 (or 40) repeats x 5 folds, and inside one
repeat the five folds come from a single StratifiedKFold partition, so a
gene in fold 1's test set is not in fold 2's. The five fold values of a
repeat are negatively dependent (measured intra-repeat ICC = -0.235).
Treating the folds as independent inflates every standard error by about
two-fold, which a coverage simulation confirms: the fold-level interval
covers 100 % of the time, the repeat-level constructions 93.6-95.0 %.
"""
from __future__ import annotations

import numpy as np
from scipy import stats

__all__ = [
    "ci_percentile",
    "ci_t",
    "ci_bca",
    "gene_cluster_bootstrap",
]


# --------------------------------------------------------------- fold/repeat

def ci_percentile(values: np.ndarray, rng: np.random.Generator,
                   n_boot: int) -> tuple[float, float]:
    """Percentile bootstrap CI of the mean, resampling the rows of `values`."""
    values = np.asarray(values, dtype=float)
    n = len(values)
    idx = rng.integers(0, n, size=(n_boot, n))
    means = values[idx].mean(axis=1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return float(lo), float(hi)


def ci_t(values: np.ndarray) -> tuple[float, float]:
    """Two-sided t interval for the mean (n-1 df).

    With only 10-40 clusters this is the interval that is actually
    calibrated; the percentile bootstrap under-covers at n=10.
    """
    values = np.asarray(values, dtype=float)
    n = len(values)
    if n < 2:
        return (float("nan"), float("nan"))
    se = values.std(ddof=1) / np.sqrt(n)
    tcrit = stats.t.ppf(0.975, n - 1)
    m = values.mean()
    return float(m - tcrit * se), float(m + tcrit * se)


def ci_bca(values: np.ndarray, rng: np.random.Generator,
           n_boot: int) -> tuple[float, float]:
    """Bias-corrected and accelerated (BCa) bootstrap CI of the mean."""
    values = np.asarray(values, dtype=float)
    n = len(values)
    if n < 3:
        return ci_t(values)
    obs = values.mean()
    idx = rng.integers(0, n, size=(n_boot, n))
    bm = values[idx].mean(axis=1)
    prop = float(np.mean(bm < obs))
    prop = min(max(prop, 1.0 / (2 * n_boot)), 1.0 - 1.0 / (2 * n_boot))
    z0 = stats.norm.ppf(prop)

    jack = np.array([(values.sum() - v) / (n - 1) for v in values])
    jb = jack.mean()
    num = np.sum((jb - jack) ** 3)
    den = 6.0 * (np.sum((jb - jack) ** 2) ** 1.5)
    a = num / den if den > 0 else 0.0

    zl, zh = stats.norm.ppf(0.025), stats.norm.ppf(0.975)
    al = float(np.clip(stats.norm.cdf(z0 + (z0 + zl) / (1 - a * (z0 + zl))), 0, 1))
    ah = float(np.clip(stats.norm.cdf(z0 + (z0 + zh) / (1 - a * (z0 + zh))), 0, 1))
    lo, hi = np.percentile(bm, [100 * al, 100 * ah])
    return float(lo), float(hi)


# ------------------------------------------------------------ gene cluster

def _weighted_auc_prep(scores: np.ndarray, labels: np.ndarray):
    """Pre-sort one fold's test set for fast weighted AUROC.

    Returns (order, group_start, group_id) where `group_start` marks the
    first position of every distinct score and `group_id[p]` is the group
    index of position p.
    """
    order = np.argsort(scores, kind="mergesort")
    s = scores[order]
    # boundaries of the tied groups
    new = np.empty(len(s), dtype=bool)
    new[0] = True
    new[1:] = s[1:] != s[:-1]
    group_start = np.flatnonzero(new)
    group_id = np.cumsum(new) - 1
    return order, group_start, group_id


def _weighted_auc_batch(s_sorted: np.ndarray, y_sorted: np.ndarray,
                        w: np.ndarray, group_start: np.ndarray,
                        group_id: np.ndarray) -> np.ndarray:
    """Weighted AUROC for a batch of weight vectors.

    s_sorted : ascending scores
    y_sorted : 0/1 labels in the same order
    w        : (B, n) non-negative weights (bootstrap multiplicities)
    Returns (B,) AUROC values.
    """
    wneg = w * (1.0 - y_sorted)[None, :]
    wpos = w * y_sorted[None, :]

    # negatives per tied group, then the exclusive prefix sum
    gneg = np.add.reduceat(wneg, group_start, axis=1)          # (B, G)
    cum_before = np.cumsum(gneg, axis=1) - gneg                 # exclusive
    L = cum_before[:, group_id]                                # (B, n)
    E = gneg[:, group_id]                                      # (B, n)

    num = (wpos * (L + 0.5 * E)).sum(axis=1)
    den = wpos.sum(axis=1) * wneg.sum(axis=1)
    out = np.full(len(num), np.nan)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out


def gene_cluster_bootstrap(scores: np.ndarray, degree: np.ndarray,
                           y: np.ndarray, in_test: np.ndarray, *,
                           n_boot: int = 1000, seed: int = 20261029,
                           chunk: int = 100) -> dict:
    """Cluster bootstrap over GENES, per pre-registration section 9.3.

    scores  : (F, N) RWR score of every gene in every fold
    degree  : (N,)   degree score of every gene
    y       : (N,)   label of every gene (0/1)
    in_test : (F, N) bool, True where the gene is in that fold's test set

    Procedure (one bootstrap replicate):
      1. draw N gene indices with replacement from the universe -- this
         draw is SHARED by all folds, so a gene moves as a block;
      2. for each fold, the resampled test set is the drawn genes that
         fall in that fold's test set, keeping their multiplicities;
      3. recompute AUROC for RWR and for the degree baseline on that
         resampled set and take the difference;
      4. average over folds.

    The resampled statistic is the mean over folds of
    (AUROC_RWR - AUROC_degree), i.e. the same quantity the fold-level
    and repeat-level bootstraps estimate.

    Returns a dict with the 95 % percentile CI of the mean delta, the
    bootstrap distribution's mean and SD, and the per-fold AUROC of the
    observed (unresampled) data for reference.
    """
    F, N = scores.shape
    rng = np.random.default_rng(seed)

    # ---- prepare each fold ----
    prep = []
    for j in range(F):
        idx = np.flatnonzero(in_test[j])
        if idx.size < 2:
            prep.append(None)
            continue
        s = scores[j, idx]
        yy = y[idx]
        d = degree[idx]
        o1, g1, gi1 = _weighted_auc_prep(s, yy)
        o2, g2, gi2 = _weighted_auc_prep(d, yy)
        prep.append((idx, o1, g1, gi1, o2, g2, gi2))

    # ---- observed (unresampled) per-fold deltas ----
    obs_delta = []
    for j, pr in enumerate(prep):
        if pr is None:
            continue
        idx, o1, g1, gi1, o2, g2, gi2 = pr
        a1 = _weighted_auc_batch(scores[j, idx][o1], y[idx][o1],
                                 np.ones((1, idx.size)), g1, gi1)[0]
        a2 = _weighted_auc_batch(degree[idx][o2], y[idx][o2],
                                 np.ones((1, idx.size)), g2, gi2)[0]
        obs_delta.append(a1 - a2)
    obs_delta = np.asarray(obs_delta, dtype=float)

    # ---- bootstrap ----
    dist = np.empty(n_boot, dtype=float)
    done = 0
    while done < n_boot:
        b = min(chunk, n_boot - done)
        draws = rng.integers(0, N, size=(b, N))
        # multiplicity of each gene: count over the draw axis
        mult = np.zeros((b, N), dtype=np.float64)
        rows = np.repeat(np.arange(b), N)
        np.add.at(mult, (rows, draws.ravel()), 1.0)

        # A fold whose resampled test set happens to contain no positive
        # (or no negative) has an undefined AUROC; it must be skipped for
        # that replicate only, not poison the whole replicate.
        acc = np.zeros(b, dtype=float)
        cnt = np.zeros(b, dtype=float)
        for j, pr in enumerate(prep):
            if pr is None:
                continue
            idx, o1, g1, gi1, o2, g2, gi2 = pr
            w = mult[:, idx]
            a1 = _weighted_auc_batch(scores[j, idx][o1], y[idx][o1], w[:, o1], g1, gi1)
            a2 = _weighted_auc_batch(degree[idx][o2], y[idx][o2], w[:, o2], g2, gi2)
            dlt = a1 - a2
            good = np.isfinite(dlt)
            acc[good] += dlt[good]
            cnt[good] += 1.0
        with np.errstate(invalid="ignore", divide="ignore"):
            vals = np.where(cnt > 0, acc / np.maximum(cnt, 1.0), np.nan)
        dist[done:done + b] = vals
        done += b

    dist = dist[np.isfinite(dist)]
    if dist.size == 0:
        return {"ci_low": float("nan"), "ci_high": float("nan"),
                "boot_mean": float("nan"), "boot_sd": float("nan"),
                "observed_mean": float(np.nanmean(obs_delta)),
                "observed_sd": float(np.nanstd(obs_delta, ddof=1)),
                "n_folds_used": 0, "n_boot": 0}
    lo, hi = np.percentile(dist, [2.5, 97.5])
    return {
        "ci_low": float(lo),
        "ci_high": float(hi),
        "boot_mean": float(dist.mean()),
        "boot_sd": float(dist.std(ddof=1)),
        "observed_mean": float(np.nanmean(obs_delta)),
        "observed_sd": float(np.nanstd(obs_delta, ddof=1)),
        "n_folds_used": int(len(prep)),
        "n_boot": int(dist.size),
    }
