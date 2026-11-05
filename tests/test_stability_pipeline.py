"""Integration-level tests for the stability analysis (work-plan v3.1, S-phase).

These check the properties the analysis relies on, not the numbers it produced:

* the 40-repeat split panel really is a strict superset of the released
  10-repeat panel (repeat index r is the only thing that changes);
* the SRS operator is a no-op on the ranking when the neighbourhood is the
  k nearest genes in score space (the finding the analysis reports);
* the matched-strength random control makes the ranking *less* stable;
* the monotone-flattening control leaves the rank statistics unchanged.

They are deliberately small and synthetic so that they run in seconds.
"""
from __future__ import annotations

import numpy as np
import pytest

from d1.engine.stability import (
    clms_smooth,
    monotone_flatten,
    percentile_ranks,
    random_neighbourhood_median,
)
from d1.splits import make_splits


def test_more_repeats_extend_the_panel_without_changing_it():
    """make_splits(n_rep=40) must agree with make_splits(n_rep=10) on repeats 0-9."""
    rng = np.random.default_rng(0)
    universe = np.array([f"G{i:03d}" for i in range(200)])
    positives = set(rng.choice(universe, 40, replace=False).tolist())
    a = make_splits(universe, positives, n_rep=10, k=5)
    b = make_splits(universe, positives, n_rep=40, k=5)
    assert a.equals(b[b.repeat < 10].reset_index(drop=True))
    assert b.repeat.nunique() == 40


def test_srs_operator_preserves_the_ranking():
    """The reported finding: shrinking toward the score-neighbourhood median is
    order-preserving up to near-ties, so the ranking -- and therefore the rank
    SD -- barely moves. The random control does not have this property."""
    from scipy.stats import spearmanr
    rng = np.random.default_rng(1)
    n = 400
    genes = np.array([f"G{i:04d}" for i in range(n)])
    scores = rng.random(n) ** 3          # skewed, like an RWR profile
    base = percentile_ranks(scores, genes)
    for k in (5, 20, 50):
        for lam in (0.25, 0.5, 1.0):
            out = clms_smooth(scores, genes, np.ones(n), k=k, lam=lam)
            rho = spearmanr(base, percentile_ranks(out, genes)).statistic
            assert rho > 0.99
    # the matched-strength random control destroys the ordering
    conf = np.ones(n)
    med = random_neighbourhood_median(scores, conf, k=20, n_reps=5,
                                      rng=np.random.default_rng(2))
    rand = scores + 0.7 * (med[0] - scores)
    rho_rand = spearmanr(base, percentile_ranks(rand, genes)).statistic
    assert rho_rand < spearmanr(base, percentile_ranks(
        clms_smooth(scores, genes, conf, k=20, lam=0.7), genes)).statistic


def test_srs_operator_compresses_the_scores():
    """... but it does compress them: the scores change, only the ranks do not."""
    rng = np.random.default_rng(2)
    n = 400
    genes = np.array([f"G{i:04d}" for i in range(n)])
    scores = rng.random(n) ** 3
    out = clms_smooth(scores, genes, np.ones(n), k=20, lam=1.0)
    assert np.std(out) < np.std(scores)
    assert not np.allclose(out, scores)


def test_random_control_is_worse_than_the_observed_operator():
    """The matched-strength random control must degrade the ranking."""
    rng = np.random.default_rng(3)
    n = 300
    genes = np.array([f"G{i:04d}" for i in range(n)])
    scores = rng.random(n) ** 3
    conf = np.ones(n)
    obs = clms_smooth(scores, genes, conf, k=20, lam=0.7)
    med = random_neighbourhood_median(scores, conf, k=20, n_reps=5,
                                      rng=np.random.default_rng(4))
    rand = scores + 0.7 * (med - scores)
    # the random operator changes many more ranks than the score-neighbourhood one
    def changed(a, b):
        return int((percentile_ranks(a, genes) != percentile_ranks(b, genes)).sum())
    changed_random = max(changed(rand[i], scores) for i in range(rand.shape[0]))
    assert changed_random > changed(obs, scores)


def test_monotone_flattening_leaves_rank_statistics_unchanged():
    rng = np.random.default_rng(5)
    n = 250
    genes = np.array([f"G{i:04d}" for i in range(n)])
    scores = rng.random(n)
    flat = monotone_flatten(scores)
    assert np.array_equal(percentile_ranks(scores, genes),
                          percentile_ranks(flat, genes))
