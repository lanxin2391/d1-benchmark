"""Unit tests for the fixed learner (A0). All must pass before the engine is locked."""
import numpy as np
import pytest
import scipy.sparse as sp
import networkx as nx

from d1.engine.rwr import transition_T, seed_matrix, rwr


def adj(G):
    return sp.csr_matrix(nx.to_scipy_sparse_array(G, nodelist=sorted(G.nodes()), dtype=float))


def test_transition_is_row_stochastic():
    WT = transition_T(adj(nx.karate_club_graph()))
    col_sums = np.asarray(WT.sum(axis=0)).ravel()      # columns of W^T = rows of W
    assert np.allclose(col_sums, 1.0)


def test_mass_is_conserved():
    G = nx.karate_club_graph()
    WT = transition_T(adj(G))
    P0 = seed_matrix(G.number_of_nodes(), [[0, 1, 2], [33], [5, 10]])
    P, _ = rwr(WT, P0, alpha=0.5)
    assert np.allclose(P.sum(axis=0), 1.0, atol=1e-9)
    assert (P >= 0).all()


def test_alpha_one_returns_seeds():
    G = nx.karate_club_graph()
    WT = transition_T(adj(G))
    P0 = seed_matrix(G.number_of_nodes(), [[3, 4]])
    P, n_iter = rwr(WT, P0, alpha=1.0)
    assert np.allclose(P, P0)
    assert n_iter == 1


@pytest.mark.parametrize("alpha", [0.3, 0.5, 0.7])
def test_matches_closed_form(alpha):
    """Fixed point of P = a*P0 + (1-a) W^T P is P = a (I - (1-a) W^T)^-1 P0."""
    G = nx.connected_watts_strogatz_graph(30, 4, 0.3, seed=1)
    WT = transition_T(adj(G))
    P0 = seed_matrix(30, [[0, 7], [12], [3, 19, 25]])
    P, _ = rwr(WT, P0, alpha=alpha)
    exact = alpha * np.linalg.solve(np.eye(30) - (1 - alpha) * WT.toarray(), P0)
    assert np.allclose(P, exact, atol=1e-6)


def test_star_centre_ranks_first():
    G = nx.star_graph(10)                                  # node 0 is the centre
    WT = transition_T(adj(G))
    P, _ = rwr(WT, seed_matrix(11, [list(range(1, 11))]))
    assert P[:, 0].argmax() == 0


def test_signal_stays_in_seeded_community():
    G = nx.barbell_graph(6, 0)                             # two cliques 0-5 and 6-11, one bridge
    WT = transition_T(adj(G))
    P, _ = rwr(WT, seed_matrix(12, [[0, 1]]))
    assert P[2:6, 0].min() > P[6:12, 0].max()


def test_single_vector_input():
    G = nx.karate_club_graph()
    WT = transition_T(adj(G))
    p0 = seed_matrix(34, [[0]])[:, 0]
    p, _ = rwr(WT, p0)
    assert p.shape == (34,)


def test_rejects_degree_zero_nodes():
    A = sp.csr_matrix(np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=float))
    with pytest.raises(ValueError, match="degree 0"):
        transition_T(A)


def test_rejects_asymmetric_and_self_loops():
    with pytest.raises(ValueError, match="symmetric"):
        transition_T(sp.csr_matrix(np.array([[0, 1], [0, 0]], dtype=float)))
    with pytest.raises(ValueError, match="self-loops"):
        transition_T(sp.csr_matrix(np.array([[1, 1], [1, 0]], dtype=float)))


def test_empty_seed_set_rejected():
    with pytest.raises(ValueError, match="empty"):
        seed_matrix(5, [[]])


def test_warns_when_not_converged():
    G = nx.path_graph(200)
    WT = transition_T(adj(G))
    with pytest.warns(RuntimeWarning, match="did not converge"):
        rwr(WT, seed_matrix(200, [[0]]), alpha=0.01, max_iter=3)
