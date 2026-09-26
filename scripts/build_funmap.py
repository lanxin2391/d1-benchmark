"""Build the FunMap protein network (A2 task).

Reads:
    data/raw/networks/funmap/funmap.tsv
        headerless, two columns of HGNC symbols, ~196,800 edges, no weights
Writes:
    data/processed/networks/funmap.tsv        (contract C2)
    data/processed/networks/funmap.meta.json  (contract C2b)

The FunMap raw file has no header. pandas would otherwise treat the first
data row as column names, so we pass header=None.

Decision references:
    D-01 unweighted -> weight column set to 1.0 (kept for future use)
    D-02 HGNC symbols are already in the file (no ID mapping needed)
"""
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from d1.networks.io import finalize_edges, check_network_file  # noqa: E402

RAW_PATH = os.path.join(ROOT, "data", "raw", "networks", "funmap", "funmap.tsv")
OUT_DIR = os.path.join(ROOT, "data", "processed", "networks")
NET_ID = "funmap"


def main():
    # header=None: FunMap has no column header row
    edges = pd.read_csv(RAW_PATH, sep="\t", header=None, names=["gene_a", "gene_b"],
                        dtype={"gene_a": str, "gene_b": str})
    print(f"raw funmap edges: {len(edges)}")

    df, meta = finalize_edges(
        edges, net_id=NET_ID,
        meta={"source_file": "funmap.tsv", "threshold": None,
              "description": "FunMap protein-derived network (linkedomics.org)"},
        out_dir=OUT_DIR,
    )
    print(f"after finalize_edges: {len(df)} edges, {meta['n_nodes_lcc']} nodes in LCC")
    print(f"  removed self-loops: {meta['n_self_loops_removed']}")
    print(f"  removed duplicate edges: {meta['n_edges_dedup'] - len(df)} (kept max weight)")

    # Final sanity check
    check_network_file(os.path.join(OUT_DIR, f"{NET_ID}.tsv"))
    print(f"OK: {NET_ID}.tsv passes contract C2")


if __name__ == "__main__":
    main()