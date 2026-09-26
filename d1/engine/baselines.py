"""Baselines the network prior is measured against.

Control 1 (protocol 4.4): degree-only ranking. The primary endpoint is always measured against
this baseline, never against chance: a network prior is only interesting if it beats
"rank genes by how many neighbours they have".
"""
import numpy as np
import scipy.sparse as sp


def degree_scores(A):
    """Unweighted degree of every gene (number of distinct neighbours) in adjacency A."""
    A = sp.csr_matrix(A)
    return np.asarray((A != 0).sum(axis=1)).ravel().astype(float)
