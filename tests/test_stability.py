"""Tests for the stability-regularised scoring operators (d1/engine/stability.py).

The properties checked here are the ones the pre-registered analysis relies on:
the neighbourhood is the true score-space k nearest, the weighted median is the
weighted-median minimiser, the operator never sees a label, and rank-based
statistics are invariant to monotone score transforms (the argument behind
control C3).
"""
from __future__ import annotations

import numpy as np
import pytest
from sklearn.metrics import roc_auc_score

from d1.engine.stability import (
    clms_smooth,
    crowding_counts,
    global_median_shrinkage,
    knn_positions,
    leave_one_repeat_out_confidence,
    monotone_flatten,
    percentile_ranks,
    random_neighbourhood_median,
    scott_bandwidth,
    score_matrix_by_repeat,
    weighted_median,
)


def test_percentile_ranks_follows_D07_tie_break():
    # equal scores -> gene symbol ascending decides (D-07)
    genes = np.array(["B", "A", "C"])
    scores = np.array([1.0, 1.0, 2.0])
    p = percentile_ranks(scores, genes)
    # C (score 2) is rank 1 -> percentile 0; A (score 1, symbol first) rank 2; B rank 3
    assert p[genes == "C"] == 0.0
    assert p[genes == "A"] == pytest.approx(1 / 2)
    assert p[genes == "B"] == pytest.approx(1.0)


def test_percentile_ranks_are_a_permutation():
    rng = np.random.default_rng(0)
    genes = np.array([f"G{i:03d}" for i in range(50)])
    scores = rng.random(50)
    p = percentile_ranks(scores, genes)
    assert np.all(np.diff(p) != 0)
    assert p.min() == 0.0 and p.max() == 1.0


def test_knn_positions_returns_the_true_k_nearest_in_score_space():
    rng = np.random.default_rng(1)
    x = np.sort(rng.random(200) * 10)
    for k in (1, 5, 20):
        pos = knn_positions(x, k)
        assert pos.shape == (200, k)
        for i in (0, 1, 37, 100, 199):
            d = np.abs(x - x[i])
            d[i] = np.inf
            expect = np.sort(np.argsort(d, kind="stable")[:k])
            got = np.sort(pos[i])
            # the selected set must have the same distances as the true k nearest
            assert np.allclose(np.sort(d[got]), np.sort(d[expect]))


def test_knn_positions_never_returns_self():
    rng = np.random.default_rng(2)
    x = np.sort(rng.random(60))
    pos = knn_positions(x, 7)
    for i in range(60):
        assert i not in set(pos[i].tolist())


def test_weighted_median_equals_plain_median_for_uniform_weights():
    rng = np.random.default_rng(3)
    v = rng.random((11, 13))
    m = weighted_median(v, np.ones_like(v))
    assert np.allclose(m, np.median(v, axis=-1))


def test_weighted_median_minimises_the_weighted_absolute_deviation():
    rng = np.random.default_rng(4)
    v = rng.random(9)
    w = rng.random(9) + 0.05
    m = weighted_median(v[None, :], w[None, :])[0]
    grid = np.linspace(v.min() - 1, v.max() + 1, 20001)
    obj = np.abs(grid[:, None] - v[None, :]) @ w
    best = grid[np.argmin(obj)]
    assert abs(m - best) < 1e-3


def test_weighted_median_handles_batches():
    rng = np.random.default_rng(5)
    v = rng.random((4, 7, 5))
    w = rng.random((4, 7, 5)) + 0.1
    m = weighted_median(v, w)
    assert m.shape == (4, 7)
    for a in range(4):
        for b in range(7):
            assert m[a, b] == pytest.approx(weighted_median(v[a, b][None, :],
                                                             w[a, b][None, :])[0])


def test_lambda_zero_is_the_identity():
    rng = np.random.default_rng(6)
    genes = np.array([f"G{i}" for i in range(40)])
    s = rng.random(40)
    out = clms_smooth(s, genes, np.ones(40), k=5, lam=0.0)
    assert np.allclose(out, s)


def test_operator_never_reads_labels():
    """Permuting the labels cannot change the smoothed scores: the operator's
    signature does not accept a label vector at all."""
    rng = np.random.default_rng(7)
    genes = np.array([f"G{i}" for i in range(30)])
    s = rng.random(30)
    conf = rng.random(30)
    out = clms_smooth(s, genes, conf, k=4, lam=0.5)
    # same input, genes permuted together with the scores -> the ranking-based
    # quantity that matters is unchanged
    perm = rng.permutation(30)
    out2 = clms_smooth(s[perm], genes[perm], conf[perm], k=4, lam=0.5)
    assert np.allclose(out[perm], out2)


def test_monotone_flatten_leaves_ranks_and_auroc_unchanged():
    """Control C3: a monotone flattening changes neither rank SD nor AUROC."""
    rng = np.random.default_rng(8)
    genes = np.array([f"G{i}" for i in range(200)])
    s = rng.random(200)
    y = (rng.random(200) < 0.3).astype(int)
    flat = monotone_flatten(s)
    assert np.allclose(percentile_ranks(s, genes), percentile_ranks(flat, genes))
    assert roc_auc_score(y, s) == pytest.approx(roc_auc_score(y, flat))
    assert np.std(s) > np.std(flat)          # the scores really were compressed


def test_global_median_shrinkage_moves_everything_towards_the_median():
    s = np.array([0.0, 1.0, 2.0, 10.0])
    out = global_median_shrinkage(s, 0.5)
    med = np.median(s)
    assert np.allclose(out, s + 0.5 * (med - s))
    assert out.max() < s.max()


def test_crowding_counts_counts_neighbours_in_the_window():
    s = np.array([0.0, 0.1, 0.2, 5.0])
    counts = crowding_counts(s, width=0.15)
    assert counts[0] == 2       # 0.0 and 0.1 within 0.15
    assert counts[1] == 3       # 0.0, 0.1, 0.2
    assert counts[3] == 1       # 5.0 alone


def test_scott_bandwidth_is_positive_on_spread_data():
    rng = np.random.default_rng(9)
    assert scott_bandwidth(rng.random(500) * 100) > 0


def test_leave_one_repeat_out_confidence_excludes_the_held_out_repeat():
    # repeat 1 is an outlier. The confidence of repeat 1 is computed from the
    # other repeats only, so the outlier does not lower its own weight
    # (leave-one-repeat-out); the other repeats do see it and lose confidence.
    s = np.array([[1.0, 1.0], [100.0, 1.0], [1.0, 1.0]])
    c = leave_one_repeat_out_confidence(s)
    assert c.shape == s.shape
    assert c[1, 0] == pytest.approx(1.0)     # its own outlier is not in the estimate
    assert c[0, 0] < 0.5                     # repeat 0 sees the outlier
    assert c[2, 0] < 0.5
    assert np.all(c[:, 1] > 0.9)             # the stable gene keeps a high weight


def test_random_neighbourhood_median_shape_and_weighting():
    rng = np.random.default_rng(10)
    s = rng.random(50)
    conf = np.ones(50)
    m = random_neighbourhood_median(s, conf, k=5, n_reps=3, rng=rng)
    assert m.shape == (3, 50)
    # with uniform weights the local median is the plain median of 6 values
    # -> it must lie inside the range of the scores
    assert np.all(m >= s.min()) and np.all(m <= s.max())
    assert not np.allclose(m[0], m[1])   # different replicates


def test_score_matrix_by_repeat_places_each_fold_in_its_repeat():
    # fold 0 and 1 belong to repeat 0, folds 2 and 3 to repeat 1; each gene is
    # tested in exactly one fold of each repeat.
    scores = np.array([[1.0, np.nan], [np.nan, 2.0], [3.0, np.nan], [np.nan, 4.0]])
    fold_repeat = np.array([0, 0, 1, 1])
    m = score_matrix_by_repeat(scores, fold_repeat)
    assert m.shape == (2, 2)
    assert m[0, 0] == 1.0 and m[0, 1] == 2.0
    assert m[1, 0] == 3.0 and m[1, 1] == 4.0
