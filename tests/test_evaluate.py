"""Tests for metrics, split checks and the harness (A1), using toy networks with known answers."""
import numpy as np
import pandas as pd
import pytest
import networkx as nx

from d1.engine.baselines import degree_scores
from d1.engine.evaluate import (check_splits, evaluate, margins, precision_at_k, run_arm,
                                summarize)
from d1.networks.io import finalize_edges, load_network
from d1.toy import make_toy_splits, toy_network


# ---------------- metrics
def test_precision_at_k_tie_break_by_symbol():
    scores = np.array([1.0, 1.0, 1.0, 0.0])
    y = np.array([0, 1, 0, 1])
    genes = np.array(["C", "A", "B", "D"])          # ties -> A, B, C order
    assert precision_at_k(scores, y, genes, 1) == 1.0
    assert precision_at_k(scores, y, genes, 2) == 0.5
    assert np.isnan(precision_at_k(scores, y, genes, 10))


def test_evaluate_perfect_and_inverted():
    y = np.array([1, 1, 0, 0, 0])
    genes = np.array(list("abcde"))
    good = evaluate(np.array([5, 4, 3, 2, 1]), y, genes, ks=(2,))
    bad = evaluate(np.array([1, 2, 3, 4, 5]), y, genes, ks=(2,))
    assert good["auroc"] == 1.0 and good["p_at_2"] == 1.0
    assert bad["auroc"] == 0.0 and bad["p_at_2"] == 0.0


def test_evaluate_single_class_gives_nan():
    out = evaluate(np.array([1.0, 2.0]), np.array([0, 0]), np.array(["a", "b"]))
    assert np.isnan(out["auroc"])


def test_degree_scores():
    G = nx.star_graph(4)
    A = nx.to_scipy_sparse_array(G, nodelist=range(5))
    assert list(degree_scores(A)) == [4, 1, 1, 1, 1]


# ---------------- split checks
def test_check_splits_accepts_toy_and_rejects_broken():
    s = make_toy_splits([f"g{i}" for i in range(100)], {f"g{i}" for i in range(20)})
    assert check_splits(s)
    assert (s.groupby("repeat").size() == 100).all()
    broken = pd.concat([s, s.iloc[[0]]])
    with pytest.raises(ValueError, match="more than one fold"):
        check_splits(broken)


# ---------------- harness on toy scenarios
def _arm(tmp_path, scenario, n=1500, n_pos=100, n_rep=3):
    edges, pos = toy_network(scenario, n=n, n_pos=n_pos, seed=1)
    finalize_edges(edges, f"toy_{scenario}", out_dir=tmp_path)
    genes, A = load_network(tmp_path / f"toy_{scenario}.tsv")
    splits = make_toy_splits(genes, pos & set(genes), n_rep=n_rep)
    res = run_arm((genes, A), splits, f"toy_{scenario}", "toy", "native", alphas=(0.5,))
    return res, summarize(res).iloc[0]


def test_result_table_shape(tmp_path):
    res, _ = _arm(tmp_path, "random", n_rep=2)
    assert len(res) == 2 * 5 * 2                        # 2 repeats x 5 folds x 2 methods
    assert res[res.method == "degree"].alpha.isna().all()
    assert (res[res.method == "rwr"].n_iter < 200).all()
    m = margins(res)
    assert len(m) == 10 and np.allclose(m.delta_auroc, m.rwr_auroc - m.degree_auroc)


def test_random_labels_no_signal(tmp_path):
    _, s = _arm(tmp_path, "random")
    assert abs(s.rwr_auroc - 0.5) < 0.1
    assert abs(s.degree_auroc - 0.5) < 0.1


def test_hub_labels_degree_wins(tmp_path):
    _, s = _arm(tmp_path, "hubs")
    assert s.degree_auroc > 0.85
    assert s.delta_auroc < 0.05                         # RWR adds little beyond degree


def test_planted_module_rwr_beats_degree(tmp_path):
    _, s = _arm(tmp_path, "module")
    assert abs(s.degree_auroc - 0.5) < 0.1              # degrees were preserved
    assert s.delta_auroc > 0.15


def test_split_genes_outside_network_rejected(tmp_path):
    edges, pos = toy_network("random", n=300, n_pos=30, seed=2)
    finalize_edges(edges, "t", out_dir=tmp_path)
    genes, A = load_network(tmp_path / "t.tsv")
    splits = make_toy_splits(list(genes) + ["NOT_IN_NETWORK"], pos, n_rep=1)
    with pytest.raises(ValueError, match="not in network"):
        run_arm((genes, A), splits, "t", "toy", "native")


def test_shared_universe_subset_works(tmp_path):
    """D-12: splits on a subset of genes; propagation still uses the whole network."""
    edges, pos = toy_network("module", n=1500, n_pos=100, seed=3)
    finalize_edges(edges, "t", out_dir=tmp_path)
    genes, A = load_network(tmp_path / "t.tsv")
    subset = set(genes[::2]) | pos                      # half the genes + all positives
    splits = make_toy_splits(subset, pos, n_rep=2)
    res = run_arm((genes, A), splits, "t", "toy", "shared")
    assert res.n_test.groupby(res.repeat).sum().iloc[0] == 2 * len(subset)  # 2 methods
