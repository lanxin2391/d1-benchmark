"""Corrected statistical layer for the D1 benchmark (protocol section 4.6).

Three defects in the original `scripts/run_stats_analysis.py` are fixed
here (see pre-registration amendment 2):

1. ALPHA FILTER. The original `_bootstrap_ci` and `_tost` selected rows
   with `a4[(a4.network == net) & (a4.label_set == label)]` and no
   alpha filter. `pivot_table` then averaged the RWR arm over
   alpha in {0.3, 0.5, 0.7} while the degree arm (whose alpha is NaN)
   stayed single-valued, so every released CI, TOST verdict and A5
   z-score was computed on the mean over the three alphas.

2. TOST VERDICT ORDER. The original code tested
   `if -margin <= mean <= margin` BEFORE the equivalence branch, so any
   cell whose point estimate happened to sit inside the margin was
   labelled BOUND and never reached the equivalence test. 24 of 24
   released rows are BOUND for that reason.

3. NO CI IN THE TOST OUTPUT. The original table carried only two
   one-sided p-values; the equivalence interval was never written, so the
   released table could not be checked against the pre-registration.

The resampling unit is now explicit. Three units are reported side by
side (decision of 2026-10-29, see work plan v3.1 Q1):

    fold    the per-fold delta values
    repeat  the per-repeat mean deltas (PRIMARY; the coverage
            simulation of work plan A1/A2 shows this is the calibrated
            construction, coverage 0.936)
    gene    the gene-cluster bootstrap of pre-registration section 9.3
            (coverage not established; the conservative bound)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

__all__ = [
    "load_primary",
    "unit_values",
    "tost",
]


def load_primary(combined: str | Path, alpha: float = 0.5) -> pd.DataFrame:
    """Load the per-fold deltas at ONE alpha.

    The degree rows carry alpha = NaN (degree does not depend on alpha),
    so they must be kept unconditionally; only the RWR rows are filtered.
    """
    df = pd.read_csv(combined, sep="\t")
    df = df[(df["method"] == "degree") | (df["alpha"] == alpha)]
    sub = df[df["method"].isin(["rwr", "degree"])]
    piv = (sub.pivot_table(index=["network", "label_set", "repeat", "fold"],
                            columns="method", values="auroc")
              .reset_index().dropna(subset=["rwr", "degree"]))
    piv["delta"] = piv["rwr"] - piv["degree"]
    return piv


def unit_values(piv: pd.DataFrame, unit: str) -> np.ndarray:
    """The observations a CI or a test should be computed on.

    fold   : one value per (repeat, fold) pair
    repeat : one value per repeat (the mean of its five folds)
    """
    if unit == "fold":
        return piv["delta"].to_numpy()
    if unit == "repeat":
        return (piv.groupby(["network", "label_set", "repeat"])["delta"]
                   .mean().to_numpy())
    raise ValueError(f"unknown unit {unit!r}")


def tost(values: np.ndarray, margin: float = 0.05, alpha: float = 0.05) -> dict:
    """Two one-sided tests for equivalence to zero within +-margin.

    H0 : |mu| >= margin      (the effect is at least the margin)
    H1 : |mu| <  margin      (the effect is indistinguishable from 0)

    H0 is rejected when BOTH one-sided p-values are below `alpha`; that is
    equivalent to the (1 - 2*alpha) confidence interval lying entirely
    inside [-margin, +margin]. The (1 - alpha) interval is reported as
    well, because the pre-registration's wording refers to a 95 % interval.

    A cell that is EQUIVALENT to zero is a NULL: there is no evidence of
    an effect, in either direction.
    """
    values = np.asarray(values, dtype=float)
    n = len(values)
    mean = float(values.mean())
    sd = float(values.std(ddof=1))
    out = {"mean": mean, "sd": sd, "n": n, "margin": margin, "alpha": alpha}
    if n < 2 or sd == 0:
        out.update({"t_low": np.nan, "p_low": np.nan, "t_high": np.nan,
                    "p_high": np.nan, "ci90_low": np.nan, "ci90_high": np.nan,
                    "ci95_low": np.nan, "ci95_high": np.nan,
                    "verdict": "zero variance" if sd == 0 else "insufficient data"})
        return out

    se = sd / np.sqrt(n)
    df = n - 1
    # lower test: H0: mean <= -margin ; reject when the mean is above it
    t_low = (mean + margin) / se
    p_low = 1.0 - stats.t.cdf(t_low, df)
    # upper test: H0: mean >= +margin ; reject when the mean is below it
    t_high = (mean - margin) / se
    p_high = stats.t.cdf(t_high, df)

    t90 = stats.t.ppf(1 - alpha, df)
    t95 = stats.t.ppf(1 - alpha / 2, df)
    ci90 = (mean - t90 * se, mean + t90 * se)
    ci95 = (mean - t95 * se, mean + t95 * se)

    equivalent = bool(p_low < alpha and p_high < alpha)
    if equivalent:
        verdict = ("EQUIVALENT to zero "
                   f"(mean={mean:.4f}, margin=+-{margin}, "
                   f"90% CI [{ci90[0]:+.4f}, {ci90[1]:+.4f}] inside)")
    else:
        verdict = ("NOT equivalent to zero "
                   f"(mean={mean:.4f}, margin=+-{margin}, "
                   f"90% CI [{ci90[0]:+.4f}, {ci90[1]:+.4f}])")

    out.update({"t_low": float(t_low), "p_low": float(p_low),
                "t_high": float(t_high), "p_high": float(p_high),
                "ci90_low": float(ci90[0]), "ci90_high": float(ci90[1]),
                "ci95_low": float(ci95[0]), "ci95_high": float(ci95[1]),
                "equivalent": equivalent,
                "verdict": verdict})
    return out
