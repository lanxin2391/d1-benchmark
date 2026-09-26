"""Synthetic networks with known answers, for testing the harness before real data arrives.

Three scenarios, each built on a scale-free (Barabasi-Albert) background network:

  random : positives are random genes.
           Expected: RWR ~ 0.5 and degree ~ 0.5 AUROC; margin ~ 0.
  hubs   : positives are drawn from the 10% highest-degree genes.
           Expected: degree AUROC is very high and RWR adds little or nothing
           (the "prior is only a degree proxy" outcome of hypothesis H1).
  module : positives are random genes that are then wired into a module by
           degree-PRESERVING edge swaps, so their degrees are unchanged.
           Expected: degree ~ 0.5, RWR clearly higher; a real margin over degree.

NOTE: `make_toy_splits` is only for toy data. The official split files are produced by
Person B's d1/splits.py (contract C4).
"""
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.model_selection import StratifiedKFold


def _names(n):
    return np.array([f"G{i:05d}" for i in range(n)])


def _plant_module(G, pos, n_swaps, rng):
    """Degree-preserving swaps: (u,u'),(v,v') -> (u,v),(u',v') for positives u,v and
    non-positive neighbours u', v'. Every node keeps its degree."""
    pos_set = set(pos)
    done, tries = 0, 0
    while done < n_swaps and tries < 50 * n_swaps:
        tries += 1
        u, v = rng.choice(pos, 2, replace=False)
        if G.has_edge(u, v):
            continue
        nu = [x for x in G.neighbors(u) if x not in pos_set]
        nv = [x for x in G.neighbors(v) if x not in pos_set]
        if not nu or not nv:
            continue
        u2, v2 = rng.choice(nu), rng.choice(nv)
        if u2 == v2 or G.has_edge(u2, v2):
            continue
        G.remove_edge(u, u2)
        G.remove_edge(v, v2)
        G.add_edge(u, v)
        G.add_edge(u2, v2)
        done += 1
    return done


def toy_network(scenario, n=3000, m=3, n_pos=150, seed=0):
    """Return (edges DataFrame with gene_a, gene_b, weight; set of positive gene names)."""
    rng = np.random.default_rng(seed)
    G = nx.barabasi_albert_graph(n, m, seed=seed)
    nodes = np.arange(n)
    if scenario == "random":
        pos = rng.choice(nodes, n_pos, replace=False)
    elif scenario == "hubs":
        deg = np.array([G.degree(i) for i in nodes])
        top = nodes[np.argsort(-deg)[: n // 10]]
        pos = rng.choice(top, n_pos, replace=False)
    elif scenario == "module":
        pos = rng.choice(nodes, n_pos, replace=False)
        _plant_module(G, list(pos), n_swaps=n_pos // 2, rng=rng)
    else:
        raise ValueError("scenario must be 'random', 'hubs' or 'module'")
    names = _names(n)
    e = np.array(G.edges())
    edges = pd.DataFrame({"gene_a": names[e[:, 0]], "gene_b": names[e[:, 1]], "weight": 1.0})
    return edges, set(names[pos])


def make_toy_splits(universe, positives, n_rep=10, k=5, base_seed=20260926):
    """Stratified 5-fold x 10-repeat split table in the C4 format (toy use only)."""
    genes = np.array(sorted(universe))
    y = np.isin(genes, list(positives)).astype(int)
    rows = []
    for r in range(n_rep):
        skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=base_seed + r)
        for f, (_, test) in enumerate(skf.split(genes, y)):
            rows.append(pd.DataFrame({"gene": genes[test], "repeat": r, "fold": f, "y": y[test]}))
    return pd.concat(rows, ignore_index=True)
