"""Reading and writing network files (contract C2 / C2b).

Every network in the panel is written by `finalize_edges`, so all six networks are cleaned in
exactly the same way; every analysis reads networks with `load_network`.
"""
import datetime as _dt
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

NETWORK_COLUMNS = ["gene_a", "gene_b", "weight"]


def _adjacency(a, b, genes):
    """Symmetric 0/1 CSR adjacency for edge lists a, b over the sorted gene array `genes`."""
    index = pd.Series(np.arange(len(genes)), index=genes)
    i = index.loc[a].to_numpy()
    j = index.loc[b].to_numpy()
    n = len(genes)
    A = sp.coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n)).tocsr()
    A = A + A.T
    A.data[:] = 1.0
    return A


def finalize_edges(df, net_id, meta=None, out_dir="data/processed/networks", write=True):
    """Clean an edge list and (optionally) write {net_id}.tsv + {net_id}.meta.json.

    df : DataFrame with columns gene_a, gene_b and optionally weight. Gene columns must
         already be HGNC symbols; NaN means "could not be mapped" and the edge is dropped.
    Steps: drop unmapped -> drop self-loops -> order pair alphabetically -> merge duplicates
    (keep max weight) -> keep largest connected component -> sort.
    Returns (edges DataFrame, meta dict).
    """
    meta = dict(meta or {})
    df = df.copy()
    if "weight" not in df.columns:
        df["weight"] = 1.0
    df = df[NETWORK_COLUMNS]
    n_input = len(df)

    df = df.dropna(subset=["gene_a", "gene_b"])
    df["gene_a"] = df["gene_a"].astype(str)
    df["gene_b"] = df["gene_b"].astype(str)
    n_mapped = len(df)

    df = df[df.gene_a != df.gene_b]
    n_self = n_mapped - len(df)

    lo = np.where(df.gene_a < df.gene_b, df.gene_a, df.gene_b)
    hi = np.where(df.gene_a < df.gene_b, df.gene_b, df.gene_a)
    df = pd.DataFrame({"gene_a": lo, "gene_b": hi, "weight": df.weight.to_numpy()})
    df = df.groupby(["gene_a", "gene_b"], as_index=False, sort=False)["weight"].max()
    n_dedup = len(df)

    genes = np.unique(np.concatenate([df.gene_a.to_numpy(), df.gene_b.to_numpy()]))
    if len(genes) == 0:
        raise ValueError(f"{net_id}: no edges left after cleaning")
    A = _adjacency(df.gene_a.to_numpy(), df.gene_b.to_numpy(), genes)
    n_comp, labels = connected_components(A, directed=False)
    largest = np.bincount(labels).argmax()
    keep = set(genes[labels == largest])
    df = df[df.gene_a.isin(keep)]          # both ends are in the same component
    df = df.sort_values(["gene_a", "gene_b"]).reset_index(drop=True)

    meta.update({
        "net_id": net_id,
        "n_input_edges": int(n_input),
        "n_mapped_edges": int(n_mapped),
        "n_self_loops_removed": int(n_self),
        "n_edges_dedup": int(n_dedup),
        "n_nodes_before_lcc": int(len(genes)),
        "n_components": int(n_comp),
        "n_nodes_lcc": int(len(keep)),
        "n_edges_lcc": int(len(df)),
        "date": _dt.date.today().isoformat(),
        "script": os.path.basename(sys.argv[0]) if sys.argv and sys.argv[0] else "interactive",
    })
    if write:
        os.makedirs(out_dir, exist_ok=True)
        df.to_csv(os.path.join(out_dir, f"{net_id}.tsv"), sep="\t", index=False)
        with open(os.path.join(out_dir, f"{net_id}.meta.json"), "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)
    return df, meta


def load_network(path):
    """Read a C2 network file. Returns (genes: sorted np.ndarray of str, A: symmetric 0/1 CSR).

    Weights are ignored on purpose (decision D-01: primary analysis is unweighted).
    """
    df = pd.read_csv(path, sep="\t", dtype={"gene_a": str, "gene_b": str})
    missing = set(NETWORK_COLUMNS[:2]) - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {missing} (contract C2)")
    genes = np.unique(np.concatenate([df.gene_a.to_numpy(), df.gene_b.to_numpy()]))
    return genes, _adjacency(df.gene_a.to_numpy(), df.gene_b.to_numpy(), genes)


def check_network_file(path):
    """Assert that a file follows contract C2. Returns a small summary dict."""
    df = pd.read_csv(path, sep="\t", dtype={"gene_a": str, "gene_b": str})
    problems = []
    if list(df.columns[:3]) != NETWORK_COLUMNS:
        problems.append(f"header should be {NETWORK_COLUMNS}, got {list(df.columns)}")
    if df[["gene_a", "gene_b"]].isna().any().any():
        problems.append("missing gene names")
    if (df.gene_a == df.gene_b).any():
        problems.append("self-loops present")
    if not (df.gene_a < df.gene_b).all():
        problems.append("some rows have gene_a >= gene_b")
    if df.duplicated(["gene_a", "gene_b"]).any():
        problems.append("duplicate edges")
    genes, A = load_network(path)
    n_comp, _ = connected_components(A, directed=False)
    if n_comp != 1:
        problems.append(f"network has {n_comp} components (should be LCC only)")
    if problems:
        raise ValueError(f"{path} breaks contract C2: " + "; ".join(problems))
    return {"n_nodes": int(len(genes)), "n_edges": int(len(df))}
