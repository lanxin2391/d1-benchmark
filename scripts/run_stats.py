"""Statistical significance tests for the D1 benchmark.

For each (network, label) arm at alpha=0.5 (primary endpoint):

    1. One-sided t-test: is delta_AUROC > 0 across the 50 folds?
    2. Multiple-comparison correction: Benjamini-Hochberg FDR across
       the 24 arms of Table 1.

Also produce the same for alpha in {0.3, 0.5, 0.7} as a sensitivity
analysis.

Output:
    results/tables/results_stats.tsv
        rows = (network, label, alpha)
        cols = mean, std, n, t, p_one_sided, p_two_sided, q_FDR
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

NET_ORDER = ["funmap", "rna_coexp", "string_phys700",
             "intact", "string_full700", "reactome"]
LABEL_ORDER = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]


def _bh_q(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted p-values."""
    p = np.asarray(pvals)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    # Enforce monotonicity (from the bottom)
    q_min = np.minimum.accumulate(q[::-1])[::-1]
    q[order] = np.minimum(q, q_min[order])
    return np.clip(q, 0, 1)


def main():
    combined = pd.read_csv(ROOT / "results" / "runs" / "a4_combined.tsv",
                            sep="\t")
    # Per-fold delta_AUROC.
    from d1.engine.evaluate import margins
    m = margins(combined)

    rows = []
    for alpha in (0.3, 0.5, 0.7):
        for net in NET_ORDER:
            for lab in LABEL_ORDER:
                sub = m[(m.alpha == alpha) & (m.network == net)
                         & (m.label_set == lab)]
                if len(sub) == 0:
                    continue
                deltas = sub.delta_auroc.to_numpy()
                # one-sided t-test: is the mean > 0?
                t_one, p_one = stats.ttest_1samp(deltas, popmean=0,
                                                   alternative="greater")
                # two-sided for completeness
                t_two, p_two = stats.ttest_1samp(deltas, popmean=0)
                # Cohen's d (effect size)
                d_cohen = float(deltas.mean() / (deltas.std(ddof=1) or 1.0))
                rows.append({
                    "network": net,
                    "label_set": lab,
                    "alpha": alpha,
                    "n_folds": len(deltas),
                    "mean": float(deltas.mean()),
                    "sd": float(deltas.std(ddof=1)),
                    "t_one_sided": float(t_one),
                    "p_one_sided": float(p_one),
                    "p_two_sided": float(p_two),
                    "cohens_d": d_cohen,
                })

    df = pd.DataFrame(rows)
    # FDR across all rows (24 arms x 3 alpha = 72 tests).
    df["q_FDR"] = _bh_q(df["p_one_sided"].to_numpy())

    out = ROOT / "results" / "tables" / "results_stats.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nWrote {out}")

    # Print summary at alpha=0.5 (primary endpoint)
    primary = df[df.alpha == 0.5].sort_values(["network", "label_set"])
    print("\nPrimary analysis (alpha=0.5):")
    print(primary[["network", "label_set", "n_folds", "mean", "sd",
                     "p_one_sided", "q_FDR"]].to_string(index=False))

    # Count significant arms
    sig_q05 = (primary.q_FDR < 0.05).sum()
    sig_q10 = (primary.q_FDR < 0.10).sum()
    pos = (primary["mean"] > 0).sum()
    print(f"\nSummary at alpha=0.5:")
    print(f"  arms with delta_AUROC > 0: {pos}/{len(primary)}")
    print(f"  arms with q_FDR < 0.05 (significant): {sig_q05}/{len(primary)}")
    print(f"  arms with q_FDR < 0.10 (marginally significant): {sig_q10}/{len(primary)}")


if __name__ == "__main__":
    main()