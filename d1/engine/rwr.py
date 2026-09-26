"""Random walk with restart (RWR): the fixed learner of the D1 benchmark.

LOCKED (protocol section 4.3). The learner is deliberately simple and is never tuned:
the study measures what the network prior contributes, and any learner change would confound
that. Any change to this file needs a dated entry in docs/decisions.md.

    adjacency   A: symmetric, unweighted (decision D-01), no self-loops, LCC only
    transition  W = D^-1 A                   (row-stochastic: every row sums to 1)
    iteration   P <- alpha * P0 + (1 - alpha) * W^T @ P
    restart     alpha = 0.5 (primary); 0.3 and 0.7 are sensitivity arms
    stop        worst-column L1 change < 1e-8, or 200 iterations
    seeds       training-fold positives, uniform weights (each column of P0 sums to 1)

Reading the formula: at every step the walker either jumps back to a random seed gene
(probability alpha) or moves to a random neighbour of its current gene (probability 1 - alpha).
P[:, j] is the long-run probability of finding the walker at each gene for seed set j.
Genes close to many seeds (in the network) get high scores.
"""
import warnings

import numpy as np
import scipy.sparse as sp

ALPHA = 0.5
ALPHA_SENSITIVITY = (0.3, 0.7)
TOL = 1e-8
MAX_ITER = 200


def transition_T(A):
    """Return W^T (CSR) for a symmetric sparse adjacency matrix A, where W = D^-1 A.

    Raises ValueError if A is not square, not symmetric, has self-loops, or has
    genes with degree 0 (restrict the network to its LCC first).
    """
    A = sp.csr_matrix(A, dtype=np.float64)
    n, m = A.shape
    if n != m:
        raise ValueError(f"adjacency must be square, got {A.shape}")
    if (abs(A - A.T) > 1e-12).nnz:
        raise ValueError("adjacency must be symmetric (undirected network)")
    if A.diagonal().any():
        raise ValueError("adjacency has self-loops; remove them first")
    deg = np.asarray(A.sum(axis=1)).ravel()
    if (deg == 0).any():
        raise ValueError(f"{int((deg == 0).sum())} nodes have degree 0; "
                         "restrict the network to its largest connected component")
    W = sp.diags(1.0 / deg) @ A
    return sp.csr_matrix(W.T)


def seed_matrix(n, seed_index_lists):
    """Build P0 (n x k): column j is uniform over the seed indices in seed_index_lists[j]."""
    P0 = np.zeros((n, len(seed_index_lists)))
    for j, idx in enumerate(seed_index_lists):
        idx = np.unique(np.asarray(idx, dtype=int))
        if idx.size == 0:
            raise ValueError(f"seed set {j} is empty")
        P0[idx, j] = 1.0 / idx.size
    return P0


def rwr(WT, P0, alpha=ALPHA, tol=TOL, max_iter=MAX_ITER):
    """Run RWR for one or many seed vectors at once.

    WT : sparse (n x n), from transition_T
    P0 : array (n,) or (n x k); each column sums to 1
    Returns (P, n_iter). P has the same shape as P0. n_iter is the number of iterations
    until the slowest column converged.
    """
    if not 0.0 < alpha <= 1.0:
        raise ValueError("alpha must be in (0, 1]")
    P0 = np.asarray(P0, dtype=np.float64)
    squeeze = P0.ndim == 1
    if squeeze:
        P0 = P0[:, None]
    P = P0.copy()
    delta = np.inf
    for it in range(1, max_iter + 1):
        P_new = alpha * P0 + (1.0 - alpha) * (WT @ P)
        delta = np.abs(P_new - P).sum(axis=0).max()
        P = P_new
        if delta < tol:
            break
    else:
        warnings.warn(f"RWR did not converge in {max_iter} iterations "
                      f"(last change {delta:.2e} > tol {tol:.0e})", RuntimeWarning)
    return (P[:, 0] if squeeze else P), it
