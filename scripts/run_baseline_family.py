"""Control 19 -- the baseline family sweep (#1).

Our headline says: the NETWORK adds information beyond DEGREE, measured as
delta-AUROC = AUROC(RWR) - AUROC(degree). Two objections follow, and this
script answers both of them on the same splits, the same folds and the same
test sets as the primary analysis.

  (A) "degree is a strawman baseline."
      Answer: compare degree against the other classical centrality baselines
      on the SAME test folds -- closeness (harmonic, landmark-estimated),
      eigenvector centrality and PageRank. If one of them is systematically
      stronger than degree, our margin is measured against a weak baseline.

  (B) "RWR is a lucky special case."
      Answer: run a second, locked diffusion setting inside the family --
      the literature's default PageRank damping 0.85, which in our parameter-
      isation is RWR with restart 0.15 (outside our pre-registered
      sensitivity range of 0.3-0.7) -- and compare its margin with the
      primary RWR margin. (The proximity learner of Control 17 covers a
      different family entirely; this one covers the same family with a
      different damping.)

Every scorer is evaluated on the identical test genes, so all margins are
paired within a fold. The RWR alpha = 0.5 margin recomputed here is also a
reproducibility check on Table 1: it must match
results/tables/a4_table1_delta_auroc_v40.tsv.

Outputs:
  results/tables/baseline_family_folds.tsv
  results/tables/baseline_family_arms.tsv
  results/tables/baseline_family_summary.tsv
  docs/baseline_family.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sparse
from scipy.sparse.csgraph import dijkstra

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.baselines import degree_scores  # noqa: E402
from d1.engine.evaluate import evaluate, read_splits  # noqa: E402
from d1.engine.rwr import rwr, seed_matrix, transition_T  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETS = ROOT / "data/processed/networks"
SPLITS = ROOT / "data/processed/splits_v40"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/baseline_family.md"

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new",
          "ncg_known"]
PRIMARY_LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]

ALPHA_PRIMARY = 0.5
ALPHA_PPR = 0.15          # PageRank damping 0.85 in our parameterisation
N_LANDMARKS = 256
RNG_SEED = 20260927
PAGERANK_STATIC_ALPHA = 0.85

VALIDATION = []


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ------------------------------------------------------------------ scores
def harmonic_closeness(A, n_landmarks=N_LANDMARKS, seed=RNG_SEED):
    """Closeness centrality estimated from random landmarks.

    The exact quantity needs all-pairs distances (O(V*E)); the standard
    landmark estimator uses the distances from K random sources:
        harmonic(v) ~ mean_l 1 / d(v, l)
    which is what we use. Unreachable pairs contribute 0 (the graph is
    connected by construction -- every network is LCC only).
    """
    n = A.shape[0]
    rng = np.random.default_rng(seed)
    landmarks = rng.choice(n, size=min(n_landmarks, n), replace=False)
    acc = np.zeros(n)
    for l in landmarks:
        d = np.asarray(dijkstra(A, indices=int(l), unweighted=True,
                                 min_only=True, directed=False), dtype=float)
        acc += np.where(np.isfinite(d), 1.0 / np.maximum(d, 1.0), 0.0)
    return acc / len(landmarks)


def exact_harmonic(A):
    """Exact harmonic centrality from all-pairs distances (small graphs only)."""
    from scipy.sparse.csgraph import shortest_path
    d = np.asarray(shortest_path(A, method="D", unweighted=True,
                                  directed=False), dtype=float)
    inv = np.where(np.isfinite(d), 1.0 / np.maximum(d, 1.0), 0.0)
    return inv.sum(axis=1) / A.shape[0]


def eigenvector_centrality(A, iters=200, tol=1e-10):
    """Power-iteration eigenvector centrality (the standard definition)."""
    n = A.shape[0]
    v = np.ones(n) / np.sqrt(n)
    for _ in range(iters):
        w = A @ v
        nw = np.linalg.norm(w)
        if nw < 1e-300:
            break
        v_new = w / nw
        if np.abs(v_new - v).max() < tol:
            v = v_new
            break
        v = v_new
    return v


def pagerank_static(A, alpha=PAGERANK_STATIC_ALPHA, iters=300, tol=1e-12):
    """PageRank with uniform teleportation (no personalization)."""
    n = A.shape[0]
    deg = np.asarray(A.sum(axis=1)).ravel()
    M = sparse.diags(1.0 / deg) @ A          # column-stochastic for symmetric A
    p = np.full(n, 1.0 / n)
    for _ in range(iters):
        p_new = (1.0 - alpha) / n + alpha * (M @ p)
        if np.abs(p_new - p).max() < tol:
            p = p_new
            break
        p = p_new
    return p


# ------------------------------------------------------------------- main
def main() -> int:
    rows = []
    for net in NETWORKS:
        genes, A = load_network(str(NETS / f"{net}.tsv"))
        n = len(genes)
        A = sparse.csr_matrix(A)
        index = pd.Series(np.arange(n), index=genes)

        deg = degree_scores(A)
        log(f"{net}: {n} genes, {A.nnz // 2} edges -- static baselines")
        close = harmonic_closeness(A)
        if n <= 5000:
            ex = exact_harmonic(A)
            rho = float(pd.Series(close).corr(pd.Series(ex), method="spearman"))
            VALIDATION.append({"network": net, "n_genes": n,
                               "landmark_vs_exact_spearman": rho})
            log(f"  {net}: landmark closeness vs exact harmonic, "
                f"Spearman {rho:+.4f}")
        eigen = eigenvector_centrality(A)
        prank = pagerank_static(A)
        WT = transition_T(A)

        for lab in LABELS:
            sp_path = SPLITS / f"{lab}__native_{net}.tsv"
            if not sp_path.exists():
                continue
            sp = read_splits(sp_path)
            sp = sp.assign(idx=index.loc[sp.gene].to_numpy())

            # batch the folds: one column of P0 per (repeat, fold)
            keys = sorted(sp.groupby(["repeat", "fold"]).groups.keys())
            seed_lists, tests = [], []
            for (r, f) in keys:
                d = sp[sp.repeat == r]
                seeds = d.loc[(d.fold != f) & (d.y == 1), "idx"].to_numpy()
                test = d[d.fold == f]
                if len(seeds) == 0 or len(test) == 0:
                    seed_lists.append(np.array([], dtype=int))
                    tests.append(test)
                    continue
                seed_lists.append(seeds)
                tests.append(test)
            P0 = seed_matrix(n, [s for s in seed_lists if len(s) > 0])
            keep = [i for i, s in enumerate(seed_lists) if len(s) > 0]
            P_rwr, _ = rwr(WT, P0, alpha=ALPHA_PRIMARY)
            P_ppr, _ = rwr(WT, P0, alpha=ALPHA_PPR)

            for j, i in enumerate(keep):
                test = tests[i]
                ti = test.idx.to_numpy()
                y = test.y.to_numpy()
                g = test.gene.to_numpy()
                ev = {}
                for name, sc in (("degree", deg), ("closeness", close),
                                 ("eigen", eigen), ("pagerank", prank)):
                    ev[name] = evaluate(sc[ti], y, g)["auroc"]
                ev["rwr05"] = evaluate(P_rwr[:, j][ti], y, g)["auroc"]
                ev["rwr15"] = evaluate(P_ppr[:, j][ti], y, g)["auroc"]
                rows.append({
                    "network": net, "label": lab,
                    "repeat": int(keys[i][0]), "fold": int(keys[i][1]),
                    "n_seeds": int(len(seed_lists[i])),
                    **{f"auroc_{k}": v for k, v in ev.items()},
                    "margin_rwr05_vs_degree": ev["rwr05"] - ev["degree"],
                    "margin_rwr15_vs_degree": ev["rwr15"] - ev["degree"],
                    "margin_rwr05_vs_closeness": ev["rwr05"] - ev["closeness"],
                    "margin_rwr05_vs_eigen": ev["rwr05"] - ev["eigen"],
                    "margin_rwr05_vs_pagerank": ev["rwr05"] - ev["pagerank"],
                    "margin_closeness_vs_degree": ev["closeness"] - ev["degree"],
                    "margin_eigen_vs_degree": ev["eigen"] - ev["degree"],
                    "margin_pagerank_vs_degree": ev["pagerank"] - ev["degree"],
                })
            log(f"  {net} x {lab}: {len(keep)} folds")

    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "baseline_family_folds.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------- arms
    mcols = [c for c in df.columns if c.startswith("margin_")]
    acols = [c for c in df.columns if c.startswith("auroc_")]
    arm = df.groupby(["network", "label"]).agg(
        **{c: (c, "mean") for c in acols + mcols},
        n_folds=("fold", "size")).reset_index()
    arm.to_csv(TABLES / "baseline_family_arms.tsv", sep="\t", index=False)

    # which baseline is strongest per arm, and the margin against it
    base_cols = ["auroc_degree", "auroc_closeness", "auroc_eigen",
                 "auroc_pagerank"]
    arm["strongest_baseline"] = arm[base_cols].idxmax(axis=1).str.replace(
        "auroc_", "", regex=False)
    arm["margin_rwr05_vs_best_baseline"] = arm.apply(
        lambda r: r["auroc_rwr05"] - r[base_cols].max(), axis=1)
    arm.to_csv(TABLES / "baseline_family_arms.tsv", sep="\t", index=False)

    # -------------------------------------------------------------- summary
    summary = []
    mcols = mcols + ["margin_rwr05_vs_best_baseline"]
    for m in mcols:
        prim = arm[arm.label.isin(PRIMARY_LABELS)]
        summary.append({
            "comparison": m.replace("margin_", ""),
            "arms_primary": int(len(prim)),
            "arms_positive": int((prim[m] > 0).sum()),
            "median_primary": float(prim[m].median()),
            "mean_primary": float(prim[m].mean()),
            "min_primary": float(prim[m].min()),
            "max_primary": float(prim[m].max()),
            "arms_positive_all30": int((arm[m] > 0).sum()),
            "median_all30": float(arm[m].median()),
        })
    sdf = pd.DataFrame(summary)
    sdf.to_csv(TABLES / "baseline_family_summary.tsv", sep="\t", index=False)

    # reproducibility check against Table 1
    a4 = pd.read_csv(TABLES / "a4_table1_delta_auroc_v40.tsv", sep="\t")
    chk = arm[["network", "label", "margin_rwr05_vs_degree"]].merge(
        a4[["network", "label_set", "a0.5_mean"]],
        left_on=["network", "label"], right_on=["network", "label_set"],
        how="left")
    diff = (chk["margin_rwr05_vs_degree"] - chk["a0.5_mean"]).abs()
    log(f"Table 1 check: max |recomputed - published| = {diff.max():.2e} "
        f"over {chk['a0.5_mean'].notna().sum()} arms")

    # which baseline is strongest per arm
    write_doc(df, arm, sdf, chk, diff.max())
    log(f"wrote the baseline-family tables ({len(df)} folds, {len(arm)} arms)")
    return 0


def write_doc(df, arm, sdf, chk, max_diff):
    val = VALIDATION[0] if VALIDATION else None
    prim = arm[arm.label.isin(PRIMARY_LABELS)]
    lines = [
        "# Control 19 -- the baseline family sweep (#1)",
        "",
        "**Date:** 2026-11-12  ",
        "**Script:** `scripts/run_baseline_family.py`  ",
        "**Tables:** `baseline_family_folds.tsv`, `baseline_family_arms.tsv`, "
        "`baseline_family_summary.tsv`",
        "",
        "## Question",
        "",
        "Our endpoint is `delta-AUROC = AUROC(RWR) - AUROC(degree)`. Two "
        "objections: (A) *degree is a strawman baseline*; (B) *RWR is a lucky "
        "special case*. Both are answered here on the same splits, the same "
        "folds and the same test genes.",
        "",
        "## Scorers",
        "",
        "| scorer | kind | definition |",
        "|---|---|---|",
        "| degree | baseline (control 1) | number of distinct neighbours |",
        "| closeness | baseline | harmonic centrality, estimated from 256 "
        "random landmarks (the exact quantity needs O(V*E) all-pairs "
        "distances) |",
        "| eigen | baseline | eigenvector centrality by power iteration |",
        "| pagerank | baseline | PageRank, uniform teleportation, damping 0.85 |",
        "| rwr05 | learner (primary) | RWR, restart 0.5 (locked) |",
        "| rwr15 | learner (second damping) | RWR with restart 0.15 = the "
        "literature's default PageRank damping 0.85, outside our "
        "pre-registered sensitivity range of 0.3-0.7 |",
        "",
        "## The baselines against degree (primary 24 arms)",
        "",
        "| baseline | arms where it beats degree | median margin vs degree |",
        "|---|---|---|",
    ]
    for m in ["closeness_vs_degree", "eigen_vs_degree", "pagerank_vs_degree"]:
        r = sdf[sdf.comparison == m].iloc[0]
        lines.append(f"| {m.replace('_vs_degree','')} | "
                     f"{r.arms_positive}/{r.arms_primary} | "
                     f"{r.median_primary:+.4f} |")
    lines += [
        "",
        "## The margins that matter",
        "",
        "| comparison | positive arms | median | min | max |",
        "|---|---|---|---|---|",
    ]
    for m in ["rwr05_vs_degree", "rwr15_vs_degree", "rwr05_vs_closeness",
              "rwr05_vs_eigen", "rwr05_vs_pagerank",
              "rwr05_vs_best_baseline"]:
        r = sdf[sdf.comparison == m].iloc[0]
        lines.append(f"| {m} | {r.arms_positive}/{r.arms_primary} | "
                     f"{r.median_primary:+.4f} | {r.min_primary:+.4f} | "
                     f"{r.max_primary:+.4f} |")
    sp = prim["margin_rwr05_vs_degree"].corr(
        prim["margin_rwr05_vs_closeness"], method="spearman")
    lines += [
        "",
        "## Reading",
        "",
        f"- the strongest baseline per arm: "
        + ", ".join(f"{k} {v}/{len(prim)}"
                    for k, v in prim.strongest_baseline.value_counts().items()),
        f"- Spearman correlation between the RWR-vs-degree margin and the "
        f"RWR-vs-closeness margin across the 24 primary arms: {sp:+.3f}",
        f"- the RWR alpha = 0.5 margin recomputed here differs from Table 1 "
        f"(`a4_table1_delta_auroc_v40.tsv`) by at most {max_diff:.1e}",
        "",
        "## Caveats",
        "",
        "1. The landmark estimate of closeness is an approximation; on "
        "reactome (the smallest panel) the exact harmonic centrality agrees "
        "with it at Spearman "
        + (f"{val['landmark_vs_exact_spearman']:+.4f}" if val else "n/a")
        + ".",
        "2. `rwr15` is the same family as `rwr05` with a different damping, "
        "not an independent learner; the independent-family check is "
        "Control 17 (the proximity learner).",
        "3. Exploratory family (amendment 4): the pre-registered endpoint is "
        "untouched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
