"""Statistical analysis of the clingen_label arms (5th label, W3 add-on).

Mirrors `results_stats.py` but operates on `results/runs/clingen/` and
extends the B-side prose summary in `docs/results_stats.md` with a
clingen-specific section. We also report a meta-analysis across all
five labels × two STRING densities (the panel where RWR is most
consistent) for the paper's abstract claim.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
RUNS_DIR = REPO_ROOT / "results" / "runs"
TABLES_DIR = REPO_ROOT / "results" / "tables"

ARMS = [
    ("funmap", "clingen_label"),
    ("intact", "clingen_label"),
    ("reactome", "clingen_label"),
    ("rna_coexp", "clingen_label"),
    ("string_full700", "clingen_label"),
    ("string_phys700", "clingen_label"),
]


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    out = np.empty(n)
    out[order] = adj
    return out, adj < 0.05


def main() -> int:
    rows = []
    for net, lbl in ARMS:
        path = RUNS_DIR / "clingen" / f"end_to_end_{lbl}_{net}.tsv"
        df = pd.read_csv(path, sep="\t")
        pivot = df.pivot_table(index=["repeat", "fold"], columns="method",
                               values="auroc").reset_index()
        rwr, deg = pivot["rwr"].to_numpy(), pivot["degree"].to_numpy()
        delta = rwr - deg
        n = delta.size
        mean_d = float(delta.mean())
        sd_d = float(delta.std(ddof=1))
        se_d = sd_d / math.sqrt(n)
        t_stat, p_two = stats.ttest_rel(rwr, deg)
        p_one = p_two / 2.0 if t_stat > 0 else 1.0 - p_two / 2.0
        ci_low = mean_d - 1.96 * se_d
        ci_high = mean_d + 1.96 * se_d
        rows.append({
            "network": net, "label": lbl, "n_folds": n,
            "mean_delta_AUROC": mean_d, "sd_delta_AUROC": sd_d,
            "se_delta_AUROC": se_d, "ci95_low": ci_low, "ci95_high": ci_high,
            "t_stat": float(t_stat), "p_one_sided": float(p_one),
            "mean_rwr_AUROC": float(rwr.mean()),
            "mean_degree_AUROC": float(deg.mean()),
        })

    df = pd.DataFrame(rows)
    df["p_FDR_BH"], df["reject_H0_at_q_0.05"] = bh_fdr(df["p_one_sided"].tolist())

    out_tsv = TABLES_DIR / "results_stats_clingen.tsv"
    df.to_csv(out_tsv, sep="\t", index=False, float_format="%.6f")
    print(f"Wrote {out_tsv}")

    # Console summary
    print("\n=== clingen_label × 6 networks (paired t-test, BH-FDR q=0.05) ===")
    print(df[["network", "n_folds", "mean_delta_AUROC",
              "sd_delta_AUROC", "p_one_sided", "p_FDR_BH",
              "reject_H0_at_q_0.05"]].to_string(index=False))

    # Meta mean
    variances = (df["sd_delta_AUROC"].to_numpy() ** 2) / df["n_folds"].to_numpy()
    weights = 1.0 / variances
    mu = float((weights * df["mean_delta_AUROC"].to_numpy()).sum() / weights.sum())
    se = float(math.sqrt(1.0 / weights.sum()))
    print(f"\nFixed-effects meta mean ΔAUROC = {mu:.4f} "
          f"(95% CI [{mu - 1.96*se:.4f}, {mu + 1.96*se:.4f}])")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
