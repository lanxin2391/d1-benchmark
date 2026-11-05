"""Stability-regularised scoring (SRS): operators and metrics.

This module implements the analysis pre-registered in
``docs/stability_regularization_plan.md`` (committed *before* any result existed).
It is **exploratory** (work-plan v3.1 E4/R19): it changes the learner's output
and therefore never enters the H1 decision. The confirmatory factorial keeps the
plain RWR scores.

What is here
------------
* :func:`percentile_ranks` -- deterministic within-fold ranks (D-07 tie-break).
* :func:`crowding_counts` -- how many genes share a gene's score neighbourhood.
* :func:`knn_positions` -- the k nearest neighbours **in score space**.
* :func:`weighted_median` -- confidence-weighted local median.
* :func:`clms_smooth` -- the SRS operator (confidence-weighted local-median
  shrinkage), plus its two matched-strength controls.
* :func:`rank_sd_matrix` / :func:`score_cv_matrix` -- the stability metrics.

Conventions inherited from the harness
--------------------------------------
* D-07 tie-break: score descending, then gene symbol ascending.
* The degree baseline is never smoothed.
* No function in this module ever reads a label vector. The operator is a
  function of the score vector alone.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "percentile_ranks",
    "scott_bandwidth",
    "crowding_counts",
    "knn_positions",
    "knn_against_pool",
    "weighted_median",
    "clms_smooth",
    "random_neighbourhood_median",
    "global_median_shrinkage",
    "rank_sd_matrix",
    "score_cv_matrix",
    "score_matrix_by_repeat",
    "leave_one_repeat_out_confidence",
    "monotone_flatten",
]


# ----------------------------------------------------------------- ranking
def percentile_ranks(scores, genes):
    """Within-fold percentile ranks, 0 = best, 1 = worst.

    D-07: score descending, gene symbol ascending. Returns ranks as
    ``(rank - 1) / (n - 1)`` so that folds of different sizes are comparable.
    """
    scores = np.asarray(scores, dtype=float)
    n = scores.size
    if n < 2:
        return np.zeros(n)
    order = np.lexsort((np.asarray(genes).astype(str), -scores))
    rank = np.empty(n, dtype=np.int64)
    rank[order] = np.arange(1, n + 1)
    return (rank - 1) / (n - 1)


def rank_sd_matrix(percentiles):
    """Per-gene SD of the percentile rank across the folds a gene was tested in.

    ``percentiles`` is a (n_folds, n_genes) matrix; entries for folds in which
    the gene was not tested must be NaN. Returns ``(sd, count)``.
    """
    p = np.asarray(percentiles, dtype=float)
    with np.errstate(invalid="ignore"):
        cnt = np.sum(~np.isnan(p), axis=0)
        sd = np.nanstd(p, axis=0, ddof=1)
    sd[cnt < 2] = np.nan
    return sd, cnt


def score_cv_matrix(scores_by_repeat):
    """Per-gene coefficient of variation of the score across repeats.

    ``scores_by_repeat`` is (n_repeats, n_genes); NaN where the gene was never
    tested in that repeat (cannot happen for the primary panel, but the code
    tolerates it).
    """
    s = np.asarray(scores_by_repeat, dtype=float)
    with np.errstate(invalid="ignore"):
        mean = np.nanmean(s, axis=0)
        sd = np.nanstd(s, axis=0, ddof=1)
    cv = sd / np.maximum(np.abs(mean), 1e-300)
    cv[~np.isfinite(cv)] = np.nan
    return cv, mean, sd


# ----------------------------------------------------------------- crowding
def scott_bandwidth(scores):
    """Scott's rule bandwidth for a 1-D sample (Gaussian reference rule)."""
    x = np.asarray(scores, dtype=float)
    n = x.size
    if n < 2:
        return 0.0
    sd = np.std(x, ddof=1)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    sigma = min(sd, iqr / 1.349) if iqr > 0 else sd
    if sigma <= 0:
        return 0.0
    return 1.06 * sigma * n ** (-1.0 / 5.0)


def crowding_counts(scores, width=None):
    """Number of genes within ``+/- width`` of each gene's score.

    ``width`` defaults to Scott's bandwidth of the fold's score distribution, so
    the window adapts to the scale of the prior. This is the crowding signal:
    a gene with many near-tied neighbours has an unstable rank, because a tiny
    change in its score reorders it against them.

    Returns the counts in the input order.
    """
    x = np.asarray(scores, dtype=float)
    n = x.size
    if width is None:
        width = scott_bandwidth(x)
    if width <= 0:
        return np.ones(n, dtype=np.int64)
    srt = np.sort(x)
    lo = np.searchsorted(srt, x - width, side="left")
    hi = np.searchsorted(srt, x + width, side="right")
    return (hi - lo).astype(np.int64)


# ------------------------------------------------------------- neighbourhoods
def knn_positions(sorted_scores, k):
    """Positions (in the ascending-sorted array) of the ``k`` nearest neighbours.

    "Nearest" is by *score value*: |s_j - s_i| smallest, self excluded. The
    candidate window is the 2k positions around i, which provably contains the
    k nearest on each side, so taking the k smallest distances in that window
    gives the global k nearest. Ties in distance are broken by the position in
    the ascending-sorted array (i.e. the smaller score wins, then the symbol,
    because the caller sorts by (score, symbol)).

    Returns an (n, k) integer array of positions.
    """
    s = np.asarray(sorted_scores, dtype=float)
    n = s.size
    k = int(min(k, max(n - 1, 1)))
    pos = np.arange(n)
    left = -np.arange(k, 0, -1)
    right = np.arange(1, k + 1)
    cand = pos[:, None] + np.concatenate([left, right])[None, :]
    valid = (cand >= 0) & (cand < n)
    clipped = np.clip(cand, 0, n - 1)
    d = np.abs(s[clipped] - s[:, None])
    d = np.where(valid, d, np.inf)
    idx = np.argsort(d, axis=1, kind="stable")[:, :k]
    return np.take_along_axis(clipped, idx, axis=1)


def knn_against_pool(query_scores, pool_scores, k):
    """Indices of the ``k`` nearest pool genes (by score value) for each query.

    ``pool_scores`` is the *sorted* array of the neighbourhood pool;
    ``query_scores`` are the genes to be shrunk. Used by the inductive variant,
    whose neighbourhood is the training-fold negatives rather than the other
    test genes. Returns an (n_query, k) array of indices into the pool.
    """
    q = np.asarray(query_scores, dtype=float)
    pool = np.asarray(pool_scores, dtype=float)
    n, m = q.size, pool.size
    k = int(min(k, max(m, 1)))
    if m == 0:
        return np.zeros((n, k), dtype=np.int64)
    pos = np.searchsorted(pool, q)
    left = -np.arange(k, 0, -1)
    right = np.arange(0, k)
    cand = pos[:, None] + np.concatenate([left, right])[None, :]
    valid = (cand >= 0) & (cand < m)
    clipped = np.clip(cand, 0, m - 1)
    d = np.abs(pool[clipped] - q[:, None])
    d = np.where(valid, d, np.inf)
    idx = np.argsort(d, axis=1, kind="stable")[:, :k]
    return np.take_along_axis(clipped, idx, axis=1)


def weighted_median(values, weights):
    """Weighted median along the last axis.

    ``values`` and ``weights`` broadcast against each other on the leading axes.
    The result minimises ``sum_j w_j |v_j - x|``.
    """
    v = np.asarray(values, dtype=float)
    w = np.broadcast_to(np.asarray(weights, dtype=float), v.shape)
    order = np.argsort(v, axis=-1, kind="stable")
    vs = np.take_along_axis(v, order, axis=-1)
    ws = np.take_along_axis(w, order, axis=-1)
    cw = np.cumsum(ws, axis=-1)
    half = 0.5 * cw[..., -1:]
    idx = np.count_nonzero(cw < half, axis=-1)
    idx = np.clip(idx, 0, v.shape[-1] - 1)
    return np.take_along_axis(vs, idx[..., None], axis=-1)[..., 0]


# ----------------------------------------------------------------- operator
def _local_median(scores, genes, confidence, k, neighbour_pos=None):
    """Confidence-weighted local median of every gene, in the input order.

    ``neighbour_pos`` (optional) gives, for every gene, the positions of its
    neighbours **in the ascending-sorted score array**; when it is None the
    score-space k nearest neighbours are computed.
    """
    scores = np.asarray(scores, dtype=float)
    n = scores.size
    order = np.lexsort((np.asarray(genes).astype(str), scores))
    sv = scores[order]
    pos = knn_positions(sv, k) if neighbour_pos is None else neighbour_pos
    nbr = sv[pos]
    conf_sorted = np.asarray(confidence, dtype=float)[order]
    vals = np.concatenate([sv[:, None], nbr], axis=-1)
    wts = np.concatenate([conf_sorted[:, None], conf_sorted[pos]], axis=-1)
    med = weighted_median(vals, wts)
    out = np.empty(n, dtype=float)
    out[order] = med
    return out


def clms_smooth(scores, genes, confidence, k, lam, crowding=None,
                variant="global", neighbour_pos=None):
    """Confidence-weighted local-median shrinkage (the SRS operator).

    ``s'_g = s_g + lambda_g * (m_g - s_g)`` where ``m_g`` is the
    confidence-weighted median of the score neighbourhood of ``g``.

    ``variant="global"``   -> ``lambda_g = lam`` for every gene.
    ``variant="crowding"`` -> ``lambda_g = lam * crowd_g / max(crowd)``:
        crowded genes are shrunk fully, isolated genes barely at all.

    Returns the smoothed scores in the input order.
    """
    scores = np.asarray(scores, dtype=float)
    m = _local_median(scores, genes, confidence, k, neighbour_pos=neighbour_pos)
    if variant == "crowding":
        c = np.asarray(crowding, dtype=float)
        if c.max() <= 0:
            lam_g = np.full(c.shape, lam)
        else:
            lam_g = lam * c / c.max()
    elif variant == "global":
        lam_g = np.full(scores.shape, float(lam))
    else:
        raise ValueError(f"unknown variant {variant!r}")
    return scores + lam_g * (m - scores)


def random_neighbourhood_median(scores, confidence, k, n_reps, rng):
    """Control C1: local median over a *random* neighbourhood of the same size.

    Returns an array of shape (n_reps, n) with the local median for each
    replicate. The neighbourhood is drawn uniformly from the other genes of the
    fold (with replacement), so the operator and its strength are identical to
    :func:`clms_smooth` but the crowding information is destroyed.
    """
    s = np.asarray(scores, dtype=float)
    n = s.size
    k = int(min(k, max(n - 1, 1)))
    conf = np.asarray(confidence, dtype=float)
    # draw uniformly from all genes except self: (r + self + 1) mod n
    draws = (rng.integers(0, n - 1, size=(n_reps, n, k))
             + np.arange(n)[None, :, None] + 1) % n
    vals = np.concatenate([np.broadcast_to(s[None, :, None], (n_reps, n, 1)),
                            s[draws]], axis=-1)
    wts = np.concatenate([np.broadcast_to(conf[None, :, None], (n_reps, n, 1)),
                          conf[draws]], axis=-1)
    return weighted_median(vals, wts)


def global_median_shrinkage(scores, lam):
    """Control C2: shrink every score toward the fold's global median."""
    s = np.asarray(scores, dtype=float)
    med = float(np.median(s))
    return s + lam * (med - s)


def monotone_flatten(scores, strength=0.5):
    """Control C3: a strictly monotone flattening of the score vector.

    ``s' = lo + (hi - lo) * ((s - lo) / (hi - lo)) ** gamma`` with
    ``gamma = strength`` in (0, 1): monotone increasing, and it compresses the
    dynamic range of the scores. Used to demonstrate that rank-based statistics
    (rank SD, AUROC) are invariant to monotone score transforms, so any change
    produced by the SRS operator is a *ranking* change and not a scale change.
    """
    s = np.asarray(scores, dtype=float)
    lo, hi = s.min(), s.max()
    if hi <= lo:
        return s.copy()
    gamma = float(strength)
    if not 0.0 < gamma <= 1.0:
        raise ValueError("strength must be in (0, 1]")
    z = (s - lo) / (hi - lo)
    return lo + (hi - lo) * z ** gamma


# ------------------------------------------------------------- confidence
def score_matrix_by_repeat(scores_by_fold, fold_repeat):
    """(n_repeats, n_genes) matrix of per-repeat scores.

    ``scores_by_fold`` is the (n_folds, n_genes) matrix of scores, with NaN
    where the gene was not in that fold's test set. Every gene is tested
    exactly once per repeat, so collapsing the folds of a repeat (by taking
    the single non-NaN value per gene) gives the per-repeat score of every
    gene.
    """
    s = np.asarray(scores_by_fold, dtype=float)
    fold_repeat = np.asarray(fold_repeat)
    n_rep = int(fold_repeat.max()) + 1
    out = np.empty((n_rep, s.shape[1]))
    for r in range(n_rep):
        block = s[fold_repeat == r]
        with np.errstate(invalid="ignore"):
            out[r] = np.nanmean(block, axis=0)
    return out


def leave_one_repeat_out_confidence(scores_by_repeat, eps=1e-300):
    """Confidence weights estimated leave-one-repeat-out.

    For every repeat r and gene g::

        CV_g^(-r) = sd_{r' != r}(s_g) / max(mean_{r' != r}(s_g), eps)
        c_g^(-r)  = 1 / (1 + CV_g^(-r)^2)

    A gene whose score is stable relative to its level gets weight 1; a gene
    whose score fluctuates gets a low weight. Because the estimate for repeat r
    never uses repeat r, no evaluated repeat contributes to its own weights.

    Returns an (n_repeats, n_genes) array of confidences.
    """
    s = np.asarray(scores_by_repeat, dtype=float)
    n_rep, n_genes = s.shape
    if n_rep < 3:
        return np.ones((n_rep, n_genes))
    total = np.nansum(s, axis=0)
    total_sq = np.nansum(s ** 2, axis=0)
    out = np.empty((n_rep, n_genes))
    for r in range(n_rep):
        rest = n_rep - 1
        mean_o = (total - s[r]) / rest
        var_o = (total_sq - s[r] ** 2 - rest * mean_o ** 2) / (rest - 1)
        var_o = np.maximum(var_o, 0.0)
        cv = np.sqrt(var_o) / np.maximum(np.abs(mean_o), eps)
        out[r] = 1.0 / (1.0 + cv ** 2)
    return out
