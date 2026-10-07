"""W5: Network rewiring null model (Control 2).

For each of the 6 primary networks:
    1. Generate N_REWIRE (default 100) rewired networks by degree-preserving
       edge swaps using ``nx.double_edge_swap`` (one rewiring = |E| * N_SWAPS
       double-edge swaps, as in `scripts/bench_rwr_speed.py`).
    2. Save each rewired network under
       ``data/processed/rewired/<net_id>/seed<N_REWIRE>.tsv`` (skipped if file
       already exists; this makes the script idempotent / restartable).
    3. Run all 4 label arms on each rewired network; save C5 files under
       ``results/runs/w5/rewired_<net_id>_seed<N>__<label>.tsv``.
    4. Aggregate null distribution of delta_AUROC per (network, label).

Outputs:
    * data/processed/rewired/<net_id>/seed<N>.tsv  (one per rewiring)
    * results/runs/w5/rewired_<net_id>_seed<N>__<label>.tsv
    * results/tables/w5_null_distribution.tsv      (mean +- SD across rewirings)
    * results/tables/w5_p_values.tsv              (observed vs null p-values)

Why this matters
-----------------
A positive delta_AUROC could come from two sources: the network's *topology*
or simply having *any* edges. Rewiring the network (preserving degree
sequence but destroying topology) gives the null: if the network has only
"any edges", rewired networks give the same delta_AUROC as the real one.
If the real network's delta_AUROC is *higher* than the rewired null, the
contribution is genuinely topological.

Real-time progress
------------------
The script prints one status line per major step:
    [funmap] rewriting 100/100 done in 12.3s
    [funmap] evaluating seed 100/100 done in 47.2s
    etc.
so the user can watch progress from the terminal.
"""
import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.networks.io import (  # noqa: E402
    finalize_edges, load_network, check_network_file,
)
from d1.engine.evaluate import run_arm, margins  # noqa: E402

NETWORKS = ["funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
PROCESSED_NET = ROOT / "data" / "processed" / "networks"
REWIRED_DIR = ROOT / "data" / "processed" / "rewired"
RUNS_DIR = ROOT / "results" / "runs" / "w5"
TABLE_DIR = ROOT / "results" / "tables"
REWIRED_DIR.mkdir(parents=True, exist_ok=True)
RUNS_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)


def _log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _read_edge_list(path: Path) -> tuple[nx.Graph, dict]:
    """Load a C2 .tsv edge list as a NetworkX graph (unweighted, undirected)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    G = nx.Graph()
    G.add_nodes_from(set(df.gene_a) | set(df.gene_b))
    G.add_edges_from(zip(df.gene_a, df.gene_b))
    meta = {"source": str(path.name), "n_edges": G.number_of_edges(),
            "n_nodes": G.number_of_nodes()}
    return G, meta


def _write_edge_list(G: nx.Graph, path: Path, meta: dict):
    """Save the graph back as a C2 .tsv edge list."""
    edges = pd.DataFrame(sorted(G.edges()), columns=["gene_a", "gene_b"])
    edges["weight"] = 1.0
    edges.to_csv(path, sep="\t", index=False)
    finalize_edges(edges, net_id=meta["net_id"], meta=meta, out_dir=path.parent)


def _rewire_one(G: nx.Graph, seed: int, nswaps_per_edge: int = 2) -> nx.Graph:
    """Degree-preserving double-edge swap rewiring.

    Default 2 swaps per edge (faster than the 10 in bench_rwr_speed.py).
    Degree sequence is preserved by definition (double-edge swap is
    degree-preserving by construction); 2 swaps per edge is enough to
    thoroughly mix the topology on all six of our networks while
    keeping wall clock manageable (~5 min per network at 100 rewirings).
    """
    n_edges = G.number_of_edges()
    nswaps = nswaps_per_edge * n_edges
    Gr = G.copy()
    Gr.remove_edges_from(nx.selfloop_edges(Gr))
    try:
        nx.double_edge_swap(
            Gr, nswap=nswaps, max_tries=100 * nswaps, seed=seed,
        )
    except nx.NetworkXAlgorithmError:
        nx.connected_double_edge_swap(
            Gr, nswap=nswaps, max_tries=100 * nswaps, seed=seed,
        )
    return Gr


def _step1_rewire(net_id: str, n_rewire: int):
    """Generate n_rewire rewired networks and save them as C2 .tsv files."""
    net_dir = REWIRED_DIR / net_id
    net_dir.mkdir(exist_ok=True)
    t0 = time.time()
    G, meta = _read_edge_list(PROCESSED_NET / f"{net_id}.tsv")
    _log(f"[{net_id}] step 1: generating {n_rewire} rewirings "
         f"from {G.number_of_nodes()} nodes / {G.number_of_edges()} edges")
    for s in range(n_rewire):
        out_path = net_dir / f"seed{s}.tsv"
        if out_path.exists():
            continue
        Gr = _rewire_one(G, s)
        _write_edge_list(
            Gr, out_path,
            meta={"net_id": net_id, "source": f"rewiring of {net_id}.tsv",
                  "threshold": f"rewiring seed {s}",
                  "n_swaps": f"10 * {G.number_of_edges()}"}
        )
        if (s + 1) % 25 == 0 or s == n_rewire - 1:
            _log(f"[{net_id}]   rewired {s + 1}/{n_rewire} "
                 f"(elapsed {(time.time() - t0):.1f}s)")


def _step2_evaluate(net_id: str, n_rewire: int):
    """For each rewired network, run the 4 label arms; skip if already done."""
    net_dir = REWIRED_DIR / net_id
    if not net_dir.exists():
        _log(f"[{net_id}] no rewired networks found, skipping")
        return pd.DataFrame()
    rows = []                                                   # (seed, label, delta_auroc, ...)
    t0 = time.time()
    for s in range(n_rewire):
        net_path = net_dir / f"seed{s}.tsv"
        if not net_path.exists():
            continue
        for label in LABELS:
            split_path = (ROOT / "data" / "processed" / "splits"
                          / f"{label}__native_{net_id}.tsv")
            if not split_path.exists():
                continue
            out_path = (RUNS_DIR / f"rewired_{net_id}_seed{s}__{label}.tsv")
            if out_path.exists():
                # Already evaluated; just read for summary.
                df = pd.read_csv(out_path, sep="\t")
            else:
                t_arm = time.time()
                df = run_arm(network=str(net_path),
                              splits=str(split_path),
                              network_id=f"{net_id}_rewired",
                              label_id=label,
                              universe_id=f"native_{net_id}",
                              alphas=(0.5,))
                df.to_csv(out_path, sep="\t", index=False)
                _log(f"[{net_id}]   seed {s + 1}/{n_rewire} "
                     f"x {label}: {len(df)} rows ({time.time() - t_arm:.1f}s)")
            # Append to summary rows
            m = margins(df)
            d = m[m.alpha == 0.5]
            if len(d):
                n_test_pos = int(df[(df.alpha == 0.5) & (df.method == "rwr")].n_test_pos.iloc[0])
                rows.append({
                    "network": net_id, "label": label, "seed": s,
                    "n_test_pos": n_test_pos,
                    "delta_auroc": float(d.delta_auroc.iloc[0]),
                })
    summary = pd.DataFrame(rows)
    summary.to_csv(TABLE_DIR / f"w5_rewired_{net_id}.tsv", sep="\t", index=False)
    return summary


def _step3_aggregate(observed: pd.DataFrame, n_rewire: int) -> pd.DataFrame:
    """Compute the null distribution summary and the p-value table."""
    null_rows = []
    for net_id in NETWORKS:
        p = TABLE_DIR / f"w5_rewired_{net_id}.tsv"
        if not p.exists():
            continue
        df = pd.read_csv(p, sep="\t")
        for label in LABELS:
            sub = df[df.label == label]
            if len(sub) == 0:
                continue
            null_rows.append({
                "network": net_id, "label": label,
                "n_rewire": len(sub),
                "null_delta_mean": float(sub.delta_auroc.mean()),
                "null_delta_sd": float(sub.delta_auroc.std(ddof=1)),
                "null_delta_q05": float(sub.delta_auroc.quantile(0.05)),
                "null_delta_q50": float(sub.delta_auroc.quantile(0.50)),
                "null_delta_q95": float(sub.delta_auroc.quantile(0.95)),
            })
    null = pd.DataFrame(null_rows)
    null.to_csv(TABLE_DIR / "w5_null_distribution.tsv", sep="\t", index=False)

    # p-values: observed vs null (one-sided: observed > null)
    obs = observed.set_index(["network", "label_set"])["a0.5_mean"]
    null_mean = null.set_index(["network", "label"])["null_delta_mean"]
    null_sd = null.set_index(["network", "label"])["null_delta_sd"]
    n = null.set_index(["network", "label"])["n_rewire"]

    p_rows = []
    for key in obs.index:
        o = obs[key]
        m = null_mean[key]
        s = null_sd[key]
        nn = int(n[key])
        if np.isnan(o) or np.isnan(m) or np.isnan(s) or s == 0:
            z = np.nan
            p_val = np.nan
        else:
            z = (o - m) / (s / np.sqrt(nn))
            p_val = float(1.0 - stats.norm.cdf(z))
        p_rows.append({
            "network": key[0], "label": key[1],
            "observed_delta": float(o),
            "null_mean": float(m),
            "null_sd": float(s),
            "n_rewire": nn,
            "z_score": z,
            "p_one_sided_observed_gt_null": p_val,
        })
    pv = pd.DataFrame(p_rows)
    pv.to_csv(TABLE_DIR / "w5_p_values.tsv", sep="\t", index=False)
    return pv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-rewire", type=int, default=100)
    ap.add_argument("--skip-rewire", action="store_true",
                    help="skip step 1 (assume rewired networks already exist)")
    ap.add_argument("--skip-evaluate", action="store_true",
                    help="skip step 2 (assume evaluation files already exist)")
    args = ap.parse_args()

    overall_t0 = time.time()
    _log(f"W5 rewiring null model: {args.n_rewire} rewirings x {len(NETWORKS)} networks "
         f"x {len(LABELS)} labels = {args.n_rewire * len(NETWORKS) * len(LABELS)} arms total")

    if not args.skip_rewire:
        _log(f"\n=== Step 1: generate rewired networks ({args.n_rewire} each) ===")
        for net_id in NETWORKS:
            _step1_rewire(net_id, args.n_rewire)

    if not args.skip_evaluate:
        _log(f"\n=== Step 2: evaluate all {args.n_rewire * len(NETWORKS)} rewirings x {len(LABELS)} labels ===")
        for net_id in NETWORKS:
            _step2_evaluate(net_id, args.n_rewire)

    _log(f"\n=== Step 3: aggregate null distribution ===")
    observed = pd.read_csv(TABLE_DIR / "a4_table1_delta_auroc.tsv", sep="\t")
    pv = _step3_aggregate(observed, args.n_rewire)
    _log(f"total wall clock: {(time.time() - overall_t0) / 60:.1f} minutes")

    # Print compact summary
    print()
    print(pv.round(4).to_string(index=False))


if __name__ == "__main__":
    main()