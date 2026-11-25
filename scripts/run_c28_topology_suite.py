"""Control 28 -- the topology prior suite (is the prior's edge over degree
a generic "topology beats degree" effect?).

The headline is a network prior against a degree baseline, and Control 19
already showed that closeness beats degree in every arm while degree is
never the strongest classical baseline. This control asks the sharper
question: is the prior's increment over degree something any
topology-aware score gets, or is it specific to a diffusion prior?

Three new scorers, all evaluated on the same folds, the same test genes
and the same degree baseline as the primary endpoint:

    betweenness centrality  -- a classical "how central is this gene in
                               the network's paths" score, which is not a
                               degree and not a diffusion; computed once
                               per network (approximate for the large
                               ones, with the sample count recorded),
    the 1-hop seed count    -- how many training seeds are one edge away,
    the 2-hop seed count    -- how many training seeds are within two
                               edges (the sparse square of the adjacency
                               times the seed indicator).

Outputs:
  results/tables/c28_topology_suite.tsv
  results/tables/c28_topology_suite_summary.tsv
  docs/c28_topology_suite.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.baselines import degree_scores  # noqa: E402
from d1.engine.evaluate import evaluate, read_splits  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETS = ROOT / "data/processed/networks"
SPLITS = ROOT / "data/processed/splits_v40"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/c28_topology_suite.md"
NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
BETWEENNESS_K = 400       # the sample count for the approximation


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def betweenness(genes: list[str], A) -> tuple[np.ndarray, str]:
    """The betweenness centrality, approximated for the large networks."""
    import networkx as nx
    G = nx.from_scipy_sparse_array(A)
    n = G.number_of_nodes()
    # the graph from the sparse array carries the integer node labels,
    # so the centrality is read positionally (the gene list is in the
    # same order as the matrix rows)
    if n <= 3000:
        bc = nx.betweenness_centrality(G, normalized=True)
        return (np.array([bc[i] for i in range(n)], dtype=float), "exact")
    bc = nx.betweenness_centrality(G, k=BETWEENNESS_K, normalized=True,
                                   seed=20261125)
    return (np.array([bc[i] for i in range(n)], dtype=float),
            f"approximate (k={BETWEENNESS_K})")


def main() -> int:
    rows = []
    for net in NETWORKS:
        net_file = NETS / f"{net}.tsv"
        if not net_file.exists():
            continue
        genes, A = load_network(str(net_file))
        deg = degree_scores(A)
        t0 = time.time()
        bc, bc_kind = betweenness(genes, A)
        log(f"  {net}: the betweenness ({bc_kind}) in {time.time() - t0:.0f}s")
        A2 = (A @ A).tocsr()
        A2.data = np.minimum(A2.data, 1.0)      # count the paths once
        index = pd.Series(np.arange(len(genes)), index=genes)
        for label in LABELS:
            sp_path = SPLITS / f"{label}__native_{net}.tsv"
            if not sp_path.exists():
                continue
            sp = read_splits(sp_path)
            sp = sp.assign(idx=index.loc[sp.gene].to_numpy())
            for (r, f) in sorted(sp.groupby(["repeat", "fold"]).groups):
                d = sp[sp.repeat == r]
                seeds = d.loc[(d.fold != f) & (d.y == 1), "idx"].to_numpy()
                test = d[d.fold == f]
                if len(seeds) == 0 or len(test) == 0:
                    continue
                ti = test.idx.to_numpy()
                y = test.y.to_numpy()
                gn = test.gene.to_numpy()
                s = np.zeros(len(genes))
                s[seeds] = 1.0
                hop1 = np.asarray(A @ s).ravel()
                hop2 = np.asarray(A2 @ s).ravel()
                ev_d = evaluate(deg[ti], y, gn)
                for tag, sc in (("betweenness", bc),
                                ("the 1-hop seed count", hop1),
                                ("the 2-hop seed count", hop2)):
                    ev = evaluate(sc[ti], y, gn)
                    rows.append({
                        "network": net, "label": label, "repeat": int(r),
                        "fold": int(f), "scorer": tag,
                        "auroc": ev["auroc"],
                        "degree_auroc": ev_d["auroc"],
                        "delta_auroc": ev["auroc"] - ev_d["auroc"],
                    })
            log(f"    {net} x {label}: done")
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "c28_topology_suite.tsv", sep="\t", index=False)

    arm = df.groupby(["network", "label", "scorer"]).delta_auroc.agg(
        ["mean", "std", "size"]).reset_index()
    # the released RWR margin and the closeness margin for the same arms
    rwr = pd.read_csv(TABLES / "a4_table1_delta_auroc_v40.tsv", sep="\t")
    rwr = rwr[["network", "label_set", "a0.5_mean"]].rename(
        columns={"label_set": "label", "a0.5_mean": "rwr_delta"})
    arm = arm.merge(rwr, on=["network", "label"], how="left")
    fam = TABLES / "baseline_family_arms.tsv"
    if fam.exists():
        bf = pd.read_csv(fam, sep="\t")
        cols = [c for c in bf.columns
                if c in ("network", "label", "label_set", "closeness_delta",
                         "delta_auroc", "method", "scorer")]
        arm = arm.merge(bf[cols], on=["network", "label"], how="left")
    arm.to_csv(TABLES / "c28_topology_suite_summary.tsv", sep="\t",
               index=False)
    write_doc(arm)
    log(f"wrote the Control 28 tables ({len(df)} fold rows)")
    return 0


def write_doc(arm: pd.DataFrame) -> None:
    piv = arm.pivot_table(index=["network", "label"], columns="scorer",
                          values="mean").reset_index()
    lines = [
        "# Control 28 — the topology prior suite",
        "",
        "**Date:** 2026-11-25  ",
        "**Script:** `scripts/run_c28_topology_suite.py`  ",
        "**Tables:** `c28_topology_suite.tsv`, "
        "`c28_topology_suite_summary.tsv`",
        "",
        "## The question",
        "",
        "The headline is a network prior against a degree baseline, and "
        "Control 19 already showed that closeness beats degree in every "
        "arm while degree is never the strongest classical baseline. This "
        "control asks the sharper question: **is the prior's increment over "
        "degree something any topology-aware score gets, or is it specific "
        "to a diffusion prior?**",
        "",
        "Three new scorers, on the same folds, the same test genes and the "
        "same degree baseline as the primary endpoint: the betweenness "
        "centrality (a classical 'how central is this gene in the "
        "network's paths' score, approximated for the large networks), the "
        "number of training seeds one edge away, and the number within two "
        "edges.",
        "",
        "## The mean delta-AUROC against degree, per arm",
        "",
        "| network | label | the RWR prior | betweenness | the 1-hop seed "
        "count | the 2-hop seed count |",
        "|---|---|---|---|---|---|",
    ]
    for _, r in piv.iterrows():
        def g(c):
            v = r.get(c)
            return f"{v:+.4f}" if pd.notna(v) else "-"
        lines.append(f"| {r['network']} | {r['label']} | "
                     f"{g('rwr_delta')} | {g('betweenness')} | "
                     f"{g('the 1-hop seed count')} | "
                     f"{g('the 2-hop seed count')} |")
    lines += ["", "## Reading", ""]
    for scorer in ("betweenness", "the 1-hop seed count",
                   "the 2-hop seed count"):
        sub = piv[piv.columns[0].notna()] if False else piv
        v = piv[scorer].dropna()
        lines.append(f"* **{scorer}**: the median margin against degree "
                     f"{v.median():+.4f}, the range {v.min():+.4f} to "
                     f"{v.max():+.4f}, positive in "
                     f"{int((v > 0).sum())} of {len(v)} arms")
    lines += [
        "",
        "## Caveats",
        "",
        "1. The betweenness is approximated with 100 sampled sources for "
        "the networks with more than 3,000 nodes (exact for the rest), so "
        "its values are noisy at the level of the sampling error; the "
        "comparison is at the arm level.",
        "2. The seed-count scores are fold-dependent by construction (the "
        "seeds change), the betweenness is not (it is a property of the "
        "network alone).",
        "3. **Exploratory** (amendment 4): no pre-registered endpoint is "
        "touched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
