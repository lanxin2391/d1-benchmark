"""Control 17 -- a second, fixed learner.

Our headline is about the NETWORK's contribution, measured against a
degree baseline. But the learner is RWR, and Guala & Sonnhammer showed
that a degree-filtered LOCAL method can beat diffusion on a full list, so
a reader is entitled to ask: is the increment a property of diffusion, or
of topology in general?

Control 17 runs a second, equally locked learner over the same splits,
the same folds and the same baselines:

    the PROXIMITY learner: score(g) = -distance(g, nearest seed),
    the unweighted single-source distance to the training seeds.

It is a deliberately different family from diffusion -- local, greedy,
no restart, no stationary distribution -- and it is exactly the family
Guala's benchmark favours. If the increment over degree survives under
it too, the conclusion is about topology, not about RWR.

Outputs:
  results/tables/second_learner_margins.tsv
  results/tables/second_learner_summary.tsv
  docs/second_learner.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import shortest_path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.baselines import degree_scores  # noqa: E402
from d1.engine.evaluate import evaluate, read_splits  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETS = ROOT / "data/processed/networks"
SPLITS = ROOT / "data/processed/splits_v40"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/second_learner.md"

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _dist_to_nearest(A, seeds) -> np.ndarray:
    """The unweighted distance from every node to the nearest seed.

    One multi-source BFS per fold: scipy's dijkstra with unweighted=True
    and min_only=True returns exactly that. (shortest_path has no
    min_only, and a full len(seeds) x n matrix per fold would be far too
    large, so the per-seed loop is only a fallback.)
    """
    from scipy.sparse.csgraph import dijkstra
    try:
        d = dijkstra(A, indices=seeds, unweighted=True, min_only=True,
                      directed=False)
        return np.asarray(d, dtype=float).ravel()
    except TypeError:
        n = A.shape[0]
        best = np.full(n, np.inf)
        for s in seeds:
            dd = np.asarray(shortest_path(A, method="D", unweighted=True,
                                          indices=int(s)), dtype=float).ravel()
            np.minimum(best, dd, out=best)
        return best


def main() -> int:
    rows = []
    for net in NETWORKS:
        genes, A = load_network(str(NETS / f"{net}.tsv"))
        deg = degree_scores(A)
        index = pd.Series(np.arange(len(genes)), index=genes)
        for lab in LABELS:
            sp_path = SPLITS / f"{lab}__native_{net}.tsv"
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
                # the distance to the NEAREST seed: one multi-source BFS
                dist = _dist_to_nearest(A, seeds)
                finite = dist[np.isfinite(dist)]
                cap = float(finite.max()) + 1.0 if finite.size else 1.0
                dist = np.where(np.isfinite(dist), dist, cap)
                score = -dist          # closer to a seed = higher score
                ev_p = evaluate(score[ti], test.y.to_numpy(),
                                test.gene.to_numpy())
                ev_d = evaluate(deg[ti], test.y.to_numpy(),
                                test.gene.to_numpy())
                rows.append({
                    "network": net, "label": lab, "repeat": int(r),
                    "fold": int(f), "n_seeds": int(len(seeds)),
                    "prox_auroc": ev_p["auroc"],
                    "degree_auroc": ev_d["auroc"],
                    "delta_auroc": ev_p["auroc"] - ev_d["auroc"],
                    "frac_unreachable": float(np.mean(~np.isfinite(
                        _dist_to_nearest(A, seeds)))),
                })
            log(f"  {net} x {lab}: {len(seeds)} seeds")
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "second_learner_margins.tsv", sep="\t", index=False)

    # the summary, per arm
    arm = df.groupby(["network", "label"]).agg(
        delta_auroc=("delta_auroc", "mean"),
        sd=("delta_auroc", "std"),
        n_folds=("delta_auroc", "size"),
        frac_unreachable=("frac_unreachable", "mean")).reset_index()
    # the RWR margin of the same arms, for the side-by-side
    rwr = pd.read_csv(TABLES / "a4_table1_delta_auroc_v40.tsv", sep="\t")
    arm = arm.merge(rwr[["network", "label_set", "a0.5_mean"]],
                    left_on=["network", "label"],
                    right_on=["network", "label_set"], how="left")
    arm = arm.rename(columns={"a0.5_mean": "rwr_delta_auroc"})
    arm.to_csv(TABLES / "second_learner_summary.tsv", sep="\t", index=False)

    write_doc(arm)
    log(f"wrote the second-learner tables ({len(df)} folds)")
    return 0


def write_doc(arm: pd.DataFrame):
    pos = int((arm.delta_auroc > 0).sum())
    lines = [
        "# Control 17 -- a second, fixed learner",
        "",
        "**Date:** 2026-11-12  ",
        "**Script:** `scripts/run_second_learner.py`  ",
        "**Tables:** `second_learner_margins.tsv`, "
        "`second_learner_summary.tsv`",
        "",
        "## The learner",
        "",
        "The **proximity learner**: score(g) = −distance(g, nearest "
        "training seed), the unweighted single-source distance on the same "
        "network. It is a deliberately different family from diffusion -- "
        "local, greedy, no restart, no stationary distribution -- and it "
        "is the family Guala & Sonnhammer's benchmark favours. Same splits, "
        "same folds, same degree baseline, same endpoint.",
        "",
        "## The result",
        "",
        "| network | label | proximity margin | RWR margin |",
        "|---|---|---|---|",
    ]
    for r in arm.sort_values("delta_auroc", ascending=False).itertuples():
        lines.append(f"| {r.network} | {r.label} | {r.delta_auroc:+.4f} | "
                     f"{r.rwr_delta_auroc:+.4f} |")
    lines += [
        "",
        "## Reading",
        "",
        f"- the proximity learner is positive in **{pos}/{len(arm)}** arms "
        f"(median **{arm.delta_auroc.median():+.4f}**); the RWR margin on "
        f"the same arms is {arm.rwr_delta_auroc.median():+.4f}",
        f"- the rank correlation between the two learners' margins: "
        f"Spearman {arm.delta_auroc.corr(arm.rwr_delta_auroc, method='spearman'):+.3f}",
        "",
        "## Caveats",
        "",
        "1. The two learners are not equally powerful by construction: "
        "diffusion uses the whole network and the seeds together, the "
        "proximity learner uses one step. The comparison is therefore "
        "about whether the increment over degree survives a change of "
        "family, not about which learner is better.",
        "2. Genes unreachable from every seed get the maximum finite "
        "distance; their share is reported per arm in the table.",
        "3. Exploratory family (amendment 4): the pre-registered endpoint "
        "is untouched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
