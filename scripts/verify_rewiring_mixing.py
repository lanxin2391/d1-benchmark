"""How well mixed is the 2 x |E| rewiring? (work plan P1-1)

The frozen pre-registration (docs/pre_registration.md, line 113) asks
for **10 x |E| double-edge swaps per rewiring**; the implementation does
**2 x |E|** (the 10 x version was projected at >6 h per network) and the
run metadata used to RECORD "10 * |E|" -- i.e. the metadata claimed a
number that was never executed.

This script measures whether the 2 x |E| version is good enough, so the
deviation can be documented with evidence instead of an assertion:

  * the fraction of edges SHARED with the original network after the
    rewiring. A well-mixed rewiring with the same degree sequence shares
    essentially no edge with the original; a badly mixed one keeps many.
  * the fraction of nodes whose NEIGHBOURHOOD changed completely (no
    shared neighbour) -- the local-mixing version of the same question.

The reference point is the null expectation: two graphs with the same
degree sequence drawn independently share about |E| * (d-1)^2 / (2m-1)
of their edges, which is what a fully mixed rewiring gives.

Usage:
    E:\\anaconda\\python.exe scripts/verify_rewiring_mixing.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import networkx as nx
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_rewiring import _read_edge_list, _rewire_one  # noqa: E402

NET_DIR = ROOT / "data" / "processed" / "networks"
TABLES = ROOT / "results" / "tables"
N_SAMPLES = {"intact": 1, "string_full700": 2, "funmap": 2, "rna_coexp": 2,
             "string_phys700": 2, "reactome": 2}
SWAPS_PER_EDGE = 2          # the implementation's setting
PREREG_SWAPS_PER_EDGE = 10  # what the frozen pre-registration asks for


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def mixing_stats(G: nx.Graph, Gr: nx.Graph) -> dict:
    e1 = {frozenset((u, v)) for u, v in G.edges()}
    e2 = {frozenset((u, v)) for u, v in Gr.edges()}
    shared = len(e1 & e2)
    n1, n2 = set(G.nodes()), set(Gr.nodes())
    # nodes with no shared neighbour at all
    nb1 = {u: set(G[u]) for u in n1}
    nb2 = {u: set(Gr[u]) for u in n2}
    common = n1 & n2
    n_changed = sum(1 for u in common
                    if not (nb1[u] & nb2[u]))
    return {
        "n_nodes": len(n1), "n_edges": len(e1),
        "n_shared_edges": shared,
        "frac_shared_edges": shared / max(len(e1), 1),
        "frac_nodes_no_shared_neighbour": n_changed / max(len(common), 1),
    }


def main() -> int:
    global SWAPS_PER_EDGE
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--swaps-per-edge", type=int, default=SWAPS_PER_EDGE)
    ap.add_argument("--out", default="rewiring_mixing_check.tsv")
    args = ap.parse_args()
    SWAPS_PER_EDGE = args.swaps_per_edge

    rows = []
    for net_id, n_samp in N_SAMPLES.items():
        G, meta = _read_edge_list(NET_DIR / f"{net_id}.tsv")
        m = G.number_of_edges()
        # the expectation for two independent graphs with the same degree
        # sequence (the fully mixed reference)
        degs = [d for _, d in G.degree()]
        exp_shared = sum(d * (d - 1) for d in degs) / max(2 * m - 1, 1)
        log(f"{net_id}: {len(G)} nodes, {m} edges; "
            f"expected shared edges under full mixing ~ {exp_shared:.1f}")
        for s in range(n_samp):
            t0 = time.time()
            Gr = _rewire_one(G, s, nswaps_per_edge=SWAPS_PER_EDGE)
            st = mixing_stats(G, Gr)
            st.update({"network": net_id, "sample": s,
                       "swaps_per_edge": SWAPS_PER_EDGE,
                       "n_swaps": SWAPS_PER_EDGE * m,
                       "expected_shared_edges_full_mixing": exp_shared,
                       "seconds": time.time() - t0})
            rows.append(st)
            log(f"  sample {s}: shared edges {st['n_shared_edges']} "
                f"({st['frac_shared_edges']:.4%}, expected "
                f"{exp_shared:.1f}) | nodes with no shared neighbour "
                f"{st['frac_nodes_no_shared_neighbour']:.1%} "
                f"({time.time()-t0:.0f}s)")
    df = pd.DataFrame(rows)
    out = TABLES / args.out
    df.to_csv(out, sep="\t", index=False, lineterminator="\n")
    log(f"wrote {out.name} ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
