#!/usr/bin/env python3
"""Control 10 (Hill et al. 2019): does the learner actually USE its input?

Why this control exists
-----------------------
Hill et al. (2019, *PLoS Computational Biology* 15(12):e1007403) benchmarked 17
network-contextualisation algorithms and noticed that the algorithms differ
sharply in how they behave when the input is MEANINGLESS. They fed 10,000
random gene lists of length 200 and counted, for every node, how often it was
ranked into the top 200. The result was algorithm-specific and sometimes
extreme: ToppNet-HITS put 239 nodes into the top 200 of more than half of the
random runs, SigNet 200, Network Propagation 0, Random Walk 4.

The interpretation Hill et al. draw is the one a reviewer will draw about D1:

    an algorithm that returns almost the same nodes whatever the input is
    is not using the input at all -- it is ranking a fixed structural set, in
    practice the hubs.

D1 already contains a *degree-matched random-seed null* (control 2, A5) that
answers a related question about the MARGIN. It does not measure whether the
learner's RANKING responds to its input, which is the quantity Hill et al.
measured and the quantity a reader can check by eye.

What this script does
---------------------
For every (network, label) arm and every one of the 200 (repeat, fold)
evaluation folds, the locked learner is run twice on the SAME fold universe:

  real input    the seeds are the training-fold positives (the real signal);
  random input  the seeds are a random set of the SAME SIZE, drawn from the
                network's genes EXCLUDING the fold's own test genes -- the
                same exclusion that applies to the real seeds, so the random
                seeds are never the genes being scored.

The decisive quantities
-----------------------
1. INPUT SENSITIVITY (the test the control exists for). For each fold, the
   Spearman correlation and the top-K Jaccard overlap between the ranking
   under the real seeds and the ranking under the random seeds. If the learner
   ignored its input, these would be ~1.0; if it used it, they will be low.

2. HILL'S DEGENERACY STATISTIC. Within the RANDOM-input regime, the number of
   genes that land in the top K in more than half of the independent random
   runs. A method that ignores its input scores high here (Hill et al.:
   ToppNet-HITS 239, SigNet 200, Random Walk 4, Network Propagation 0).

3. HUB CONCENTRATION. The mean degree percentile of the top K and the
   overlap of the top K with the network's own top K by degree, in both
   regimes. This says whether the real input pulls the ranking toward the
   hubs relative to a meaningless input.

4. RANK STABILITY. The mean pairwise Spearman between two runs of the SAME
   regime: high within the real regime (the signal is reproducible), low
   within the random regime (the output is not a fixed structural set).

Outputs
-------
    results/tables/control10_random_input.tsv           per (network, label, regime)
    results/tables/control10_random_input_summary.tsv   per (network, label)
    docs/control10_random_input.md                     the reproducer

Usage
-----
    python scripts/run_control10_random_input.py
    python scripts/run_control10_random_input.py --n-random 5 --only-net reactome
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.rwr import transition_T, seed_matrix, rwr  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETWORKS = ["funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]

NET_DIR = ROOT / "data" / "processed" / "networks"
SPLIT_DIR = ROOT / "data" / "processed" / "splits_v40"
TABLE_DIR = ROOT / "results" / "tables"

ALPHA = 0.5
TOP_K = 100


def _rank_top_k(scores: np.ndarray, k: int) -> np.ndarray:
    """Positions of the k highest-scoring genes (ties broken by position)."""
    order = np.argsort(-scores, kind="mergesort")
    return order[:k]


def _persistent(top_sets: list[np.ndarray]) -> int:
    """Genes that land in the top-K in MORE THAN HALF of the runs.

    The top-K sets are given as GLOBAL network indices, so they are
    comparable across folds and repeats. This is Hill et al.'s degeneracy
    statistic.
    """
    if not top_sets:
        return 0
    counts = np.zeros(max(int(s.max()) for s in top_sets) + 1, dtype=int)
    for s in top_sets:
        counts[s] += 1
    return int((counts > len(top_sets) / 2).sum())


def run_arm(net_id: str, label: str, n_random: int, alpha: float,
            rng: np.random.Generator) -> list[dict]:
    genes, A = load_network(str(NET_DIR / f"{net_id}.tsv"))
    sp = pd.read_csv(SPLIT_DIR / f"{label}__native_{net_id}.tsv", sep="\t",
                     dtype={"gene": str})
    index = pd.Series(np.arange(len(genes)), index=genes)
    WT = transition_T(A)
    deg = np.asarray((A != 0).sum(axis=1)).ravel().astype(float)

    # degree percentile of every gene (0 = lowest, 1 = highest)
    order = np.argsort(deg, kind="mergesort")
    pct = np.empty(len(deg), dtype=float)
    pct[order] = np.arange(len(deg)) / max(len(deg) - 1, 1)
    hub_topk = set(_rank_top_k(deg, TOP_K).tolist())

    # ONE fixed gene ordering for the whole arm, so that the rankings of
    # different folds and different runs are directly comparable.
    univ = np.sort(sp.gene.unique())
    univ_idx = index.loc[univ].to_numpy()

    combos = sorted(sp.groupby(["repeat", "fold"]).groups)
    real_seeds, fold_univ_idx = [], []
    for (r, f) in combos:
        d = sp[sp.repeat == r]
        tr = d[(d.fold != f) & (d.y == 1)]
        real_seeds.append(index.loc[tr.gene].to_numpy())
        te = d[d.fold == f]
        fold_univ_idx.append(index.loc[te.gene].to_numpy())
    n_folds = len(combos)
    all_gene_idx = np.arange(len(genes))

    def _random_seed_matrix() -> list[np.ndarray]:
        """A random seed matrix of the same shape as the real one."""
        out = []
        for j in range(n_folds):
            k = len(real_seeds[j])
            pool = np.setdiff1d(all_gene_idx, fold_univ_idx[j])
            if k == 0 or len(pool) == 0:
                out.append(np.array([], dtype=int))
                continue
            out.append(rng.choice(pool, size=min(k, len(pool)), replace=False))
        return out

    def _propagate(seed_lists: list[np.ndarray]):
        P0 = seed_matrix(len(genes), seed_lists)
        P, _ = rwr(WT, P0, alpha=alpha)
        scores = [P[univ_idx, j] for j in range(n_folds)]
        tops = [univ_idx[_rank_top_k(sc, TOP_K)] for sc in scores]
        return scores, tops

    def _regime(name: str, seed_lists: list[np.ndarray]) -> dict:
        scores, tops = _propagate(seed_lists)
        deg_of_top = float(np.mean([pct[s].mean() for s in tops]))
        hub_ovl = float(np.mean(
            [len(set(s.tolist()) & hub_topk) / TOP_K for s in tops]))
        corrs = []
        for a in range(0, len(scores), 2):
            b = a + 1
            if b < len(scores):
                c = stats.spearmanr(scores[a], scores[b]).statistic
                if np.isfinite(c):
                    corrs.append(c)
        return {"network": net_id, "label": label, "regime": name,
                "n_runs": len(tops),
                "mean_seed_set_size": float(np.mean([len(s) for s in seed_lists])),
                "n_persistent_topk_genes": _persistent(tops),
                "mean_degree_percentile_of_topk": deg_of_top,
                "mean_hub_overlap_topk": hub_ovl,
                "mean_pairwise_rank_spearman": float(np.mean(corrs)) if corrs else np.nan,
                "alpha": alpha}

    # ---------------- the real input ----------------
    real_scores, real_tops = _propagate(list(real_seeds))
    rows = [_regime("real_input", list(real_seeds))]

    # ---------------- the random input (n_random independent replicates) ----
    rand_scores_all, rand_tops_all = [], []
    for rep in range(n_random):
        seeds = _random_seed_matrix()
        sc, tp = _propagate(seeds)
        rand_scores_all.append(sc)
        rand_tops_all.append(tp)
        rows.append(_regime(f"random_input_rep{rep}", seeds))

    # ---------------- the decisive comparison: real vs random, fold by fold --
    # For every random replicate and every fold, compare the ranking under
    # the real seeds with the ranking under the random seeds.
    sens_rho, sens_jac = [], []
    for rep in range(n_random):
        for j in range(n_folds):
            c = stats.spearmanr(real_scores[j], rand_scores_all[rep][j]).statistic
            if np.isfinite(c):
                sens_rho.append(c)
            a = set(real_tops[j].tolist())
            b = set(rand_tops_all[rep][j].tolist())
            sens_jac.append(len(a & b) / len(a | b))
    rows.append({
        "network": net_id, "label": label, "regime": "real_vs_random",
        "n_runs": len(sens_rho),
        "mean_seed_set_size": float(np.mean([len(s) for s in real_seeds])),
        "n_persistent_topk_genes": -1,
        "mean_degree_percentile_of_topk": np.nan,
        "mean_hub_overlap_topk": float(np.mean(sens_jac)),
        "mean_pairwise_rank_spearman": float(np.mean(sens_rho)),
        "alpha": alpha,
    })
    return rows


def main(argv: list[str] | None = None) -> int:
    global TOP_K
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    ap.add_argument("--n-random", type=int, default=5)
    ap.add_argument("--alpha", type=float, default=ALPHA)
    ap.add_argument("--top-k", type=int, default=TOP_K)
    ap.add_argument("--seed", type=int, default=20261111)
    args = ap.parse_args(argv)

    TOP_K = args.top_k
    nets = args.only_net or NETWORKS
    labels = args.only_label or LABELS
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    all_rows = []
    for net in nets:
        for lab in labels:
            t0 = time.time()
            all_rows.extend(run_arm(net, lab, args.n_random, args.alpha, rng))
            print(f"  [{net} x {lab}] {time.time() - t0:.1f}s")

    df = pd.DataFrame(all_rows)
    df.to_csv(TABLE_DIR / "control10_random_input.tsv", sep="\t",
              index=False, float_format="%.6f", lineterminator="\n")

    real = df[df.regime == "real_input"].set_index(["network", "label"])
    rnd = df[df.regime.str.startswith("random_input_rep")].groupby(
        ["network", "label"]).mean(numeric_only=True)
    sens = df[df.regime == "real_vs_random"].set_index(["network", "label"])

    summary = pd.DataFrame(index=real.index)
    summary["n_persistent_topk_genes__real"] = real["n_persistent_topk_genes"]
    summary["n_persistent_topk_genes__random"] = rnd["n_persistent_topk_genes"]
    summary["mean_degree_percentile_of_topk__real"] = \
        real["mean_degree_percentile_of_topk"]
    summary["mean_degree_percentile_of_topk__random"] = \
        rnd["mean_degree_percentile_of_topk"]
    summary["mean_hub_overlap_topk__real"] = real["mean_hub_overlap_topk"]
    summary["mean_hub_overlap_topk__random"] = rnd["mean_hub_overlap_topk"]
    summary["mean_pairwise_rank_spearman__within_real"] = \
        real["mean_pairwise_rank_spearman"]
    summary["mean_pairwise_rank_spearman__within_random"] = \
        rnd["mean_pairwise_rank_spearman"]
    summary["input_sensitivity_spearman_real_vs_random"] = \
        sens["mean_pairwise_rank_spearman"]
    summary["input_sensitivity_topk_jaccard_real_vs_random"] = \
        sens["mean_hub_overlap_topk"]
    summary = summary.reset_index()
    summary.to_csv(TABLE_DIR / "control10_random_input_summary.tsv", sep="\t",
                   index=False, float_format="%.6f", lineterminator="\n")

    print("\n=== Control 10: does the learner use its input? (Hill et al. 2019) ===")
    cols = ["network", "label",
            "n_persistent_topk_genes__real", "n_persistent_topk_genes__random",
            "input_sensitivity_spearman_real_vs_random",
            "input_sensitivity_topk_jaccard_real_vs_random",
            "mean_degree_percentile_of_topk__real",
            "mean_degree_percentile_of_topk__random"]
    print(summary[cols].round(4).to_string(index=False))
    print()
    print(f"  Hill degeneracy statistic (genes in the top-{TOP_K} of >50 % of runs):")
    print(f"    real input    {summary['n_persistent_topk_genes__real'].mean():.1f}")
    print(f"    random input  {summary['n_persistent_topk_genes__random'].mean():.1f}"
          "   (Hill et al.: ToppNet-HITS 239, SigNet 200, Random Walk 4, "
          "Network Propagation 0)")
    print(f"  INPUT SENSITIVITY (Spearman between the real-input and the "
          f"random-input ranking of the same fold):")
    print(f"    {summary['input_sensitivity_spearman_real_vs_random'].mean():.4f}"
          "   (1.0 would mean the input is ignored)")
    print(f"  INPUT SENSITIVITY (top-{TOP_K} Jaccard):")
    print(f"    {summary['input_sensitivity_topk_jaccard_real_vs_random'].mean():.4f}")
    print(f"  mean degree percentile of the top-{TOP_K}: "
          f"real {summary['mean_degree_percentile_of_topk__real'].mean():.4f}  "
          f"random {summary['mean_degree_percentile_of_topk__random'].mean():.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
