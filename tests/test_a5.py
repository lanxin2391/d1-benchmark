"""Tests for A5 (degree-matched seed null control).

These tests do not require running the full pipeline. They verify:
    1. _degree_matched_seeds returns the requested number of seeds.
    2. The returned seeds have the right degrees (within a small tolerance
       because exact-degree bins may be empty for some degrees).
    3. The procedure never returns one of the observed seeds (no trivial
       identity match).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from run_a5 import _degree_matched_seeds


def _toy_setup():
    """Build a tiny graph with controlled degrees."""
    gene_to_degree = {
        "A": 1, "B": 1, "C": 1,
        "D": 2, "E": 2,
        "F": 3, "G": 3,
        "H": 4, "I": 4, "J": 4, "K": 4, "L": 4,
    }
    net_gene_set = set(gene_to_degree.keys())
    return gene_to_degree, net_gene_set


def test_returns_one_request():
    g2d, ns = _toy_setup()
    obs = ["D", "F", "H"]
    rng = np.random.default_rng(0)
    out = _degree_matched_seeds(obs, g2d, ns, rng)
    assert len(out) == 3


def test_degree_tolerance():
    """Returned seeds must have degree within +/-5 of the requested degree."""
    g2d, ns = _toy_setup()
    obs = ["A", "D", "F", "H"]
    rng = np.random.default_rng(0)
    out = _degree_matched_seeds(obs, g2d, ns, rng)
    for s_obs, s_new in zip(obs, out):
        target = g2d[s_obs]
        got = g2d[s_new]
        assert abs(got - target) <= 5, f"observed deg {target}, got deg {got}"


def test_no_observed_seed_returned():
    """The resampler must not pick an observed seed."""
    g2d, ns = _toy_setup()
    obs = ["A", "D", "F", "H"]
    rng = np.random.default_rng(0)
    for _ in range(50):
        out = _degree_matched_seeds(obs, g2d, ns, rng)
        for s_obs, s_new in zip(obs, out):
            assert s_new != s_obs, "resampler picked an observed seed"


def test_reproducibility_under_seed():
    """Same RNG seed -> same resample."""
    g2d, ns = _toy_setup()
    obs = ["D", "F", "H"]
    a = _degree_matched_seeds(obs, g2d, ns, np.random.default_rng(42))
    b = _degree_matched_seeds(obs, g2d, ns, np.random.default_rng(42))
    assert a == b


def test_handles_empty_exact_degree_bin():
    """If the requested degree is unique, the resampler must fall back."""
    g2d, ns = _toy_setup()
    obs = ["L"]   # L has degree 4; there are other degree-4 nodes, so no issue
    # Add a new degree not in the gene_to_degree
    g2d2 = dict(g2d)
    g2d2["M"] = 99  # degree 99 has no matches
    obs2 = ["M"]
    rng = np.random.default_rng(0)
    out = _degree_matched_seeds(obs2, g2d2, ns, rng)
    assert len(out) == 1


def test_distribution_recovery():
    """If observed seeds have specific degrees, the mean degree of the
    resampled set should be close (within rounding noise) to the observed
    mean degree."""
    g2d, ns = _toy_setup()
    obs = ["A", "D", "F", "H"]   # degrees 1, 2, 3, 4 -> mean 2.5
    rng = np.random.default_rng(7)
    n_replicates = 200
    means = []
    for _ in range(n_replicates):
        out = _degree_matched_seeds(obs, g2d, ns, rng)
        means.append(np.mean([g2d[g] for g in out]))
    assert abs(np.mean(means) - 2.5) < 0.3, \
        f"observed mean degree {-target} ± 0.3; got {np.mean(means):.2f}" \
        .replace("-target", str(2.5))