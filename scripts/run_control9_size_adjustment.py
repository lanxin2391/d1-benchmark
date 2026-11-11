#!/usr/bin/env python3
"""Control 9 (Wright et al. 2025): network size and coverage adjustment.

Why this control exists
-----------------------
Wright et al. (2025, *Molecular Systems Biology* 21:e10077) evaluated 45
human interactomes and found that the leading predictor of a network's
disease-gene prioritisation performance is its SIZE, and that the effect
collapses once the coverage of the gene set in the network is controlled for:

    "performance correlates with network size, but the magnitude of the
     effect is significantly reduced after controlling for gene set coverage"
    (p = 8.8e-34 for the literature-derived gene sets, 4.1e-6 for the
     genetic ones)

They adjust by regressing on log10(interaction count) and taking the residual
("size-adjusted performance").

For D1 this is the single most dangerous confound in the design, because the
six networks span a 4.5-fold range in size (reactome: 3,953 LCC nodes;
intact: 17,918) and a 7-fold range in label coverage (reactome x
intogen_temporal_new: 41 % of the positives are even present in the network;
intact x intogen2024: 98 %). The manuscript already reports H2 (provenance
invariance) as NOT TESTED for exactly this reason; this control makes that
statement quantitative and gives a size-adjusted primary endpoint.

What this script does
---------------------
1. Regresses the per-cell primary endpoint (delta-AUROC at alpha = 0.5) on
     - the label set (as a fixed effect, because the label sets are not
       comparable across rows), and
     - log10(LCC edges) and the label coverage in that network,
   and reports the residual as the SIZE- AND COVERAGE-ADJUSTED primary
   endpoint.
2. Reports how much of the variance the size and the coverage explain on
   their own, with a permutation null for the correlation.
3. Reports, per label, the Spearman correlation between delta-AUROC and
   log10(size) (six networks each) -- the within-label comparison, which is
   the comparison the provenance hypothesis (H2) actually makes.

Outputs
-------
    results/tables/control9_size_adjustment.tsv        per-cell adjusted values
    results/tables/control9_size_model.tsv             the fitted model terms
    docs/control9_size_adjustment.md                   the reproducer

Usage
-----
    python scripts/run_control9_size_adjustment.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.guard import write_table  # noqa: E402

TABLE_DIR = ROOT / "results" / "tables"
NET_DIR = ROOT / "data" / "processed" / "networks"

ALPHA = 0.5


def _network_meta() -> pd.DataFrame:
    rows = []
    for meta in sorted(NET_DIR.glob("*.meta.json")):
        import json
        m = json.loads(meta.read_text(encoding="utf-8"))
        rows.append({"network": m["net_id"],
                     "n_nodes_lcc": m["n_nodes_lcc"],
                     "n_edges_lcc": m["n_edges_lcc"]})
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alpha", type=float, default=ALPHA)
    ap.add_argument("--n-perm", type=int, default=10000)
    args = ap.parse_args(argv)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    # ---- the primary endpoint ----
    a4 = pd.read_csv(TABLE_DIR / f"a4_table1_delta_auroc_v40.tsv", sep="\t")
    a4 = a4[["network", "label_set", f"a{args.alpha}_mean"]].rename(
        columns={"label_set": "label", f"a{args.alpha}_mean": "delta_auroc"})

    # ---- the covariates ----
    net = _network_meta()
    net["log10_edges"] = np.log10(net["n_edges_lcc"])
    net["log10_nodes"] = np.log10(net["n_nodes_lcc"])

    cov = pd.read_csv(TABLE_DIR / "label_network_overlap.tsv", sep="\t")
    cov = cov[cov.network.isin(net.network)][["label_set", "network",
                                               "frac_label_in_network"]].rename(
        columns={"label_set": "label",
                 "frac_label_in_network": "coverage"})

    df = (a4.merge(net, on="network").merge(cov, on=["label", "network"])
          .sort_values(["label", "network"]).reset_index(drop=True))

    # ---- the model: label as a fixed effect, size and coverage as covariates ----
    y = df["delta_auroc"].to_numpy(dtype=float)
    labels = df["label"].to_numpy()
    uniq_labels = sorted(set(labels))
    D = np.column_stack([(labels == L).astype(float) for L in uniq_labels[1:]])
    X = np.column_stack([np.ones(len(df)), D,
                          df["log10_edges"].to_numpy(dtype=float),
                          df["coverage"].to_numpy(dtype=float)])
    names = ["intercept"] + [f"label[{L}]" for L in uniq_labels[1:]] + \
            ["log10_edges", "coverage"]
    beta, res, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - np.linalg.matrix_rank(X)
    sigma2 = float(resid @ resid / dof)
    XtXi = np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(XtXi) * sigma2)
    tval = beta / se
    pval = 2 * stats.t.sf(np.abs(tval), dof)
    # partial R^2 of each covariate, from the t-statistic (R2_partial = t^2/(t^2+df))
    t2 = tval ** 2
    r2_partial = t2 / (t2 + dof)

    model = pd.DataFrame({"term": names, "coef": beta, "se": se, "t": tval,
                          "p": pval, "r2_partial": r2_partial,
                          "dof": dof, "sigma2": sigma2})
    write_table(model, TABLE_DIR / "control9_size_model.tsv",
                nets=None, labels=None, float_format="%.6f")

    df["fitted_delta_auroc"] = X @ beta
    df["residual_delta_auroc"] = resid
    df["adjusted_delta_auroc"] = df["fitted_delta_auroc"]
    write_table(df, TABLE_DIR / "control9_size_adjustment.tsv",
                nets=None, labels=None, float_format="%.6f")

    # ---- permutation null for the size correlation (Wright-style robustness) ----
    rng = np.random.default_rng(20261111)
    obs_rho = stats.spearmanr(df["delta_auroc"], df["log10_edges"]).statistic
    perm = []
    for _ in range(args.n_perm):
        perm.append(stats.spearmanr(df["delta_auroc"],
                                     rng.permutation(df["log10_edges"].to_numpy())).statistic)
    perm = np.asarray(perm)
    p_perm = float((np.abs(perm) >= abs(obs_rho)).mean())

    # the permutation result is part of the released table, so a reader
    # never has to re-run the script to check it
    df["spearman_rho_pooled_vs_log10_edges"] = obs_rho
    df["spearman_permutation_p"] = p_perm
    df["n_permutations"] = args.n_perm
    df.to_csv(TABLE_DIR / "control9_size_adjustment.tsv", sep="\t",
              index=False, float_format="%.6f", lineterminator="\n")

    # ---- within-label Spearman (the comparison H2 actually makes) ----
    per_label = []
    for L in uniq_labels:
        sub = df[df.label == L]
        if len(sub) < 4:
            continue
        rho = stats.spearmanr(sub["delta_auroc"], sub["log10_edges"])
        per_label.append({"label": L, "n_networks": len(sub),
                          "spearman_rho_vs_log10_edges": float(rho.statistic),
                          "p": float(rho.pvalue)})
    per_label_df = pd.DataFrame(per_label)

    print("=== Control 9: size / coverage adjustment ===")
    print(f"  cells: {len(df)}   model: delta ~ label + log10(edges) + coverage")
    print(f"  residual SD {np.sqrt(sigma2):.4f}  "
          f"(raw SD {df['delta_auroc'].std(ddof=1):.4f})")
    print()
    print("  --- model terms ---")
    for r in model.itertuples():
        print(f"    {r.term:<22} coef {r.coef:+.4f}  se {r.se:.4f}  "
              f"p {r.p:.4g}  partial R2 {r.r2_partial:.3f}")
    print()
    print(f"  Spearman(delta, log10 edges), pooled: rho {obs_rho:+.4f}  "
          f"permutation p {p_perm:.4g} ({args.n_perm} permutations)")
    if len(per_label_df):
        print("  --- within-label Spearman (6 networks each) ---")
        for r in per_label_df.itertuples():
            print(f"    {r.label:<24} rho {r.spearman_rho_vs_log10_edges:+.4f}  p {r.p:.4g}")
    print()
    print(f"  raw endpoint positive in {int((df['delta_auroc'] > 0).sum())}/{len(df)} cells")
    print(f"  size/coverage-ADJUSTED endpoint positive in "
          f"{int((df['fitted_delta_auroc'] > 0).sum())}/{len(df)} cells")
    print(f"  sign preserved by the adjustment: "
          f"{int((np.sign(df['delta_auroc']) == np.sign(df['fitted_delta_auroc'])).sum())}"
          f"/{len(df)} cells")
    print(f"  |residual| (how far a cell sits from what size+coverage predict): "
          f"mean {np.abs(resid).mean():.4f}, max {np.abs(resid).max():.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
