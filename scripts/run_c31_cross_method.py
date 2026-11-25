"""Control 31 -- the cross-method complementarity.

If two scorers pick the same genes, running the second one adds nothing;
if they pick different genes, the choice of method matters to the user.
Control 19 compared the methods by their *accuracy* against degree; this
control compares them by their *candidates*, which is the question a user
actually faces: what do I get if I run another method on top of this one?

Six scorers, on the same folds, the same test genes and the same K:

    the RWR prior (the headline learner),
    degree (the baseline),
    betweenness centrality (a path-based topology score),
    the 1-hop seed count,
    the 2-hop seed count,
    the proximity learner (the distance to the nearest seed).

For each pair we report the Jaccard overlap of their top-K picks and the
size of the union, averaged over the folds.

Outputs:
  results/tables/c31_cross_method_overlap.tsv
  results/tables/c31_cross_method_union.tsv
  docs/c31_cross_method.md
"""
from __future__ import annotations

import itertools
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.baselines import degree_scores  # noqa: E402
from d1.engine.evaluate import read_splits  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETS = ROOT / "data/processed/networks"
SPLITS = ROOT / "data/processed/splits_v40"
SCORES = ROOT / "results/scores_v40"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/c31_cross_method.md"
NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
SCORERS = ["the RWR prior", "degree", "betweenness",
           "the 1-hop seed count", "the 2-hop seed count",
           "the proximity learner"]
N_FOLDS = 20
K = 100
ALPHA = 0.5


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def dist_to_nearest(A, seeds) -> np.ndarray:
    from scipy.sparse.csgraph import dijkstra
    d = dijkstra(A, indices=seeds, unweighted=True, min_only=True,
                 directed=False)
    return np.asarray(d, dtype=float).ravel()


def betweenness(genes, A) -> np.ndarray:
    import networkx as nx
    G = nx.from_scipy_sparse_array(A)
    n = G.number_of_nodes()
    if n <= 3000:
        bc = nx.betweenness_centrality(G, normalized=True)
    else:
        bc = nx.betweenness_centrality(G, k=400, normalized=True,
                                       seed=20261125)
    return np.array([bc[i] for i in range(n)], dtype=float)


def main() -> int:
    rows = []
    union_rows = []
    for net in NETWORKS:
        net_file = NETS / f"{net}.tsv"
        if not net_file.exists():
            continue
        genes, A = load_network(str(net_file))
        deg = degree_scores(A)
        bc = betweenness(genes, A)
        A2 = (A @ A).tocsr()
        A2.data = np.minimum(A2.data, 1.0)
        index = pd.Series(np.arange(len(genes)), index=genes)
        for label in LABELS:
            sp_path = SPLITS / f"{label}__native_{net}.tsv"
            npz_path = SCORES / f"{net}__{label}.npz"
            if not sp_path.exists() or not npz_path.exists():
                continue
            d = np.load(npz_path, allow_pickle=True)
            ai = int(np.argmin(np.abs(np.asarray(d["alphas"], dtype=float)
                                      - ALPHA)))
            rwr_all = np.asarray(d["rwr"][ai], dtype=float)
            in_test = np.asarray(d["in_test"], dtype=bool)
            npz_genes = np.asarray(d["genes"]).astype(str)
            if list(npz_genes) != list(genes):
                log(f"  {net} x {label}: the gene order differs from the "
                    "dump, skipped")
                continue
            sp = read_splits(sp_path)
            sp = sp.assign(idx=index.loc[sp.gene].to_numpy())
            for (r, f) in sorted(sp.groupby(["repeat", "fold"]).groups):
                if int(r) >= N_FOLDS:
                    continue
                j = int(r) * 5 + int(f)
                if j >= in_test.shape[0]:
                    continue
                ti = sp[(sp.repeat == r) & (sp.fold == f)].idx.to_numpy()
                dsub = sp[sp.repeat == r]
                seeds = dsub.loc[(dsub.fold != f) & (dsub.y == 1),
                                 "idx"].to_numpy()
                if len(seeds) == 0 or len(ti) < K:
                    continue
                s = np.zeros(len(genes))
                s[seeds] = 1.0
                dist = dist_to_nearest(A, seeds)
                finite = dist[np.isfinite(dist)]
                cap = float(finite.max()) + 1.0 if finite.size else 1.0
                prox = -np.where(np.isfinite(dist), dist, cap)
                score_map = {
                    "the RWR prior": rwr_all[ti, j],
                    "degree": deg[ti],
                    "betweenness": bc[ti],
                    "the 1-hop seed count": np.asarray(A @ s).ravel()[ti],
                    "the 2-hop seed count": np.asarray(A2 @ s).ravel()[ti],
                    "the proximity learner": prox[ti],
                }
                picks = {}
                for name, sc in score_map.items():
                    order = np.argsort(-sc, kind="mergesort")[:K]
                    picks[name] = set(order.tolist())
                for a, b in itertools.combinations(SCORERS, 2):
                    inter = len(picks[a] & picks[b])
                    union = len(picks[a] | picks[b])
                    rows.append({
                        "network": net, "label": label, "repeat": int(r),
                        "fold": int(f), "scorer_a": a, "scorer_b": b,
                        "jaccard": inter / union if union else np.nan,
                        "n_shared": inter, "n_union": union,
                        "n_a": len(picks[a]), "n_b": len(picks[b]),
                    })
                for name in SCORERS:
                    union_rows.append({
                        "network": net, "label": label, "repeat": int(r),
                        "fold": int(f), "scorer": name,
                        "n_picked": len(picks[name]),
                        "n_in_all": len(set.intersection(*picks.values())),
                        "n_in_any": len(set.union(*picks.values())),
                        "n_unique_to_this": len(
                            picks[name] - set.union(*[picks[m] for m in SCORERS
                                                     if m != name])),
                    })
            log(f"  {net} x {label}: done")
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "c31_cross_method_overlap.tsv", sep="\t", index=False)
    u = pd.DataFrame(union_rows)
    u.to_csv(TABLES / "c31_cross_method_union.tsv", sep="\t", index=False)
    write_doc(df, u)
    log(f"wrote the Control 31 tables ({len(df)} pair rows)")
    return 0


def write_doc(df: pd.DataFrame, u: pd.DataFrame) -> None:
    piv = df.pivot_table(index=["scorer_a", "scorer_b"], values="jaccard",
                         aggfunc="mean")
    pairs = [(a, b) for a, b in piv.index]
    lines = [
        "# Control 31 — the cross-method complementarity",
        "",
        "**Date:** 2026-11-25  ",
        "**Script:** `scripts/run_c31_cross_method.py`  ",
        "**Tables:** `c31_cross_method_overlap.tsv`, "
        "`c31_cross_method_union.tsv`",
        "",
        "## The question",
        "",
        "Control 19 compared the methods by their accuracy against degree. "
        "This control compares them by their **candidates**: what does a "
        "user get by running a second method on top of the first? Six "
        "scorers on the same folds, the same test genes and the same "
        f"top {K}:",
        "",
        "* the RWR prior (the headline learner),",
        "* degree (the baseline),",
        "* the betweenness centrality (a path-based topology score),",
        "* the 1-hop seed count and the 2-hop seed count,",
        "* the proximity learner (the distance to the nearest seed).",
        "",
        "## The Jaccard overlap of the top-100 picks (the mean over the "
        "folds)",
        "",
        "| | " + " | ".join(f"vs {s}" for s in SCORERS) + " |",
        "|---|" + "---|" * len(SCORERS),
    ]
    for a in SCORERS:
        cells = []
        for b in SCORERS:
            if a == b:
                cells.append("-")
            elif (a, b) in piv.index:
                cells.append(f"{piv.loc[(a, b), 'jaccard']:.3f}")
            else:
                cells.append(f"{piv.loc[(b, a), 'jaccard']:.3f}")
        lines.append(f"| **{a}** | " + " | ".join(cells) + " |")
    lines += [
        "",
        "## How much is shared and how much is unique",
        "",
        "| scorer | the picks | the genes every method picks | the genes "
        "only this method picks |",
        "|---|---|---|---|",
    ]
    g = u.groupby("scorer")[["n_picked", "n_in_all", "n_unique_to_this"]].mean()
    for s in SCORERS:
        if s not in g.index:
            continue
        r = g.loc[s]
        lines.append(f"| {s} | {r['n_picked']:.0f} | {r['n_in_all']:.1f} | "
                     f"{r['n_unique_to_this']:.1f} |")
    lines += [
        "",
        "## Reading",
        "",
        "## Caveats",
        "",
        f"1. {N_FOLDS} folds per arm (one per repeat), the same folds for "
        "every scorer.",
        "2. The closeness centrality is not in this table: its landmark "
        "estimate costs 256 breadth-first searches per fold, and Control 19 "
        "already reports it. Its accuracy against degree (+0.0312 in 24 of "
        "24 arms) is in that section.",
        "3. The overlap is a property of the top-K sets, not of the scores; "
        "two scorers can agree on the picks and disagree on the order.",
        "4. **Exploratory** (amendment 4): no pre-registered endpoint is "
        "touched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
