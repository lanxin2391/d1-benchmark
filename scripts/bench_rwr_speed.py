"""Time the RWR engine on a FunMap-sized network (W0, task A0).

The protocol (section 5) measured on the real FunMap network (10,525 genes, 196,800 edges):
    6.2 ms per seed vector (batches of 50), 23 iterations to tol 1e-8, sparse matrix ~5 MB,
    one degree-preserving rewiring (10 x |E| swaps) ~0.7 min.
This script builds a scale-free random network of the same size and measures the same things,
so you know your PC can run the whole factorial. Paste the output into docs/engine_benchmark.md.

Usage:  python scripts/bench_rwr_speed.py            (RWR timing only, ~10 s)
        python scripts/bench_rwr_speed.py --rewire   (also time one rewiring, ~1-2 min)
"""
import argparse
import os
import platform
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402

from d1.engine.rwr import rwr, seed_matrix, transition_T  # noqa: E402

N_GENES, N_EDGES = 10_525, 196_800


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rewire", action="store_true", help="also time one 10x|E| rewiring")
    ap.add_argument("--batch", type=int, default=50)
    args = ap.parse_args()

    print(f"Python {platform.python_version()} on {platform.platform()}, "
          f"{os.cpu_count()} logical CPUs")
    t = time.time()
    G = nx.barabasi_albert_graph(N_GENES, 19, seed=0)          # ~200k edges, scale-free
    rng = np.random.default_rng(0)
    extra = G.number_of_edges() - N_EDGES
    if extra > 0:                                              # trim to exactly 196,800 edges
        edges = list(G.edges())
        for i in rng.choice(len(edges), extra, replace=False):
            G.remove_edge(*edges[i])
    G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    A = sp.csr_matrix(nx.to_scipy_sparse_array(G, dtype=float))
    print(f"network: {A.shape[0]} genes, {int(A.nnz / 2)} edges "
          f"(built in {time.time() - t:.1f} s)")

    WT = transition_T(A)
    mb = (WT.data.nbytes + WT.indices.nbytes + WT.indptr.nbytes) / 1e6
    print(f"sparse transition matrix: {mb:.1f} MB (dense would be "
          f"{A.shape[0] ** 2 * 8 / 1e9:.2f} GB)")

    n = A.shape[0]
    seeds = [rng.choice(n, 400, replace=False) for _ in range(args.batch)]  # ~ training positives
    P0 = seed_matrix(n, seeds)
    rwr(WT, P0)                                                # warm-up
    reps = 5
    t = time.time()
    for _ in range(reps):
        _, n_iter = rwr(WT, P0, alpha=0.5)
    dt = (time.time() - t) / reps
    print(f"RWR alpha=0.5: {n_iter} iterations, {dt * 1000:.0f} ms per batch of {args.batch}, "
          f"{dt * 1000 / args.batch:.2f} ms per seed vector  (protocol: 6.2 ms, 23 iterations)")

    real_arms = 6 * 4 * 50 * dt / args.batch
    print(f"projected real arms (6 networks x 4 labels x 50 folds): {real_arms:.1f} s")

    if args.rewire:
        H = G.copy()
        m = H.number_of_edges()
        t = time.time()
        nx.double_edge_swap(H, nswap=10 * m, max_tries=100 * m, seed=1)
        dt = time.time() - t
        same = all(H.degree(v) == G.degree(v) for v in G)
        print(f"one rewiring (10 x |E| = {10 * m} swaps): {dt / 60:.2f} min, degrees preserved: "
              f"{same}  (protocol: 0.7 min)")
        print(f"projected 100 rewirings x 6 networks on 1 core: {dt * 600 / 3600:.1f} h")


if __name__ == "__main__":
    main()
