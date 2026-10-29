#!/usr/bin/env python3
"""Per-cell confidence intervals under four constructions, and the H1 verdict.

Work-plan v3.1 task A4. Replaces the ad-hoc comparison that A1 used.

Four interval constructions are computed for every (network, label) cell:

  fold       percentile bootstrap over the 50 (or 200) fold values
             -- treats the folds as independent. A1 measured this at
             100 % coverage: far too conservative, unusable.
  repeat_pct percentile bootstrap over the per-repeat means
             -- this is the construction the pre-registration names
             ("cluster bootstrap"). A1 measured 0.90 coverage at 10
             repeats and 0.936 at 40 repeats, i.e. it becomes calibrated
             once the repeat count is raised.
  repeat_t   t interval on the per-repeat means (n-1 df)
  repeat_bca bias-corrected and accelerated bootstrap on the repeat means

The pre-registered H1 rule is: the 95 % CI excludes 0 in EVERY arm.

Usage
-----
    python scripts/run_ci_variants.py
    python scripts/run_ci_variants.py --combined results/runs_v40/a4_combined.tsv \\
            --out results/tables/stats_ci_variants_v40.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
DEFAULT_COMBINED = REPO / "results" / "runs" / "a4_combined.tsv"
DEFAULT_OUT = REPO / "results" / "tables" / "stats_ci_variants.tsv"
ALPHA = 0.5


def ci_percentile(v: np.ndarray, rng: np.random.Generator,
                   n_boot: int) -> tuple[float, float]:
    n = len(v)
    idx = rng.integers(0, n, size=(n_boot, n))
    return tuple(np.percentile(v[idx].mean(axis=1), [2.5, 97.5]))  # type: ignore


def ci_t(v: np.ndarray) -> tuple[float, float]:
    n = len(v)
    if n < 2:
        return (float("nan"), float("nan"))
    se = v.std(ddof=1) / np.sqrt(n)
    t = stats.t.ppf(0.975, n - 1)
    m = v.mean()
    return float(m - t * se), float(m + t * se)


def ci_bca(v: np.ndarray, rng: np.random.Generator,
           n_boot: int) -> tuple[float, float]:
    n = len(v)
    if n < 3:
        return ci_t(v)
    obs = v.mean()
    idx = rng.integers(0, n, size=(n_boot, n))
    bm = v[idx].mean(axis=1)
    prop = float(np.mean(bm < obs))
    prop = min(max(prop, 1.0 / (2 * n_boot)), 1.0 - 1.0 / (2 * n_boot))
    z0 = stats.norm.ppf(prop)
    jack = np.array([(v.sum() - x) / (n - 1) for x in v])
    jb = jack.mean()
    num = np.sum((jb - jack) ** 3)
    den = 6.0 * (np.sum((jb - jack) ** 2) ** 1.5)
    a = num / den if den > 0 else 0.0
    zl, zh = stats.norm.ppf(0.025), stats.norm.ppf(0.975)
    al = float(np.clip(stats.norm.cdf(z0 + (z0 + zl) / (1 - a * (z0 + zl))), 0, 1))
    ah = float(np.clip(stats.norm.cdf(z0 + (z0 + zh) / (1 - a * (z0 + zh))), 0, 1))
    return tuple(np.percentile(bm, [100 * al, 100 * ah]))  # type: ignore


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--combined", default=str(DEFAULT_COMBINED))
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261029)
    args = ap.parse_args(argv)

    df = pd.read_csv(args.combined, sep="\t")
    # degree rows carry alpha=NaN (degree does not depend on alpha); keep all
    # of them, and keep only the rwr rows at alpha=0.5.
    df = df[(df["method"] == "degree") | (df["alpha"] == ALPHA)]
    sub = df[df["method"].isin(["rwr", "degree"])]
    piv = (sub.pivot_table(index=["network", "label_set", "repeat", "fold"],
                            columns="method", values="auroc")
              .reset_index().dropna(subset=["rwr", "degree"]))
    piv["delta"] = piv["rwr"] - piv["degree"]

    rng = np.random.default_rng(args.seed)
    rows = []
    for (nw, lb), g in piv.groupby(["network", "label_set"]):
        g = g.sort_values(["repeat", "fold"])
        x = g["delta"].to_numpy()
        rep = g.groupby("repeat")["delta"].mean().to_numpy()
        r = {"network": nw, "label": lb,
             "mean": float(x.mean()),
             "n_folds": int(x.size),
             "n_repeats": int(rep.size),
             "mean_repeat": float(rep.mean()),
             "sd_repeat": float(rep.std(ddof=1))}
        for key, fn in (("fold", lambda: ci_percentile(x, rng, args.n_boot)),
                        ("repeat_pct", lambda: ci_percentile(rep, rng, args.n_boot)),
                        ("repeat_t", lambda: ci_t(rep)),
                        ("repeat_bca", lambda: ci_bca(rep, rng, args.n_boot))):
            lo, hi = fn()
            r[f"ci_{key}_lo"] = lo
            r[f"ci_{key}_hi"] = hi
            r[f"excl0_{key}"] = bool(not (lo <= 0 <= hi))
            r[f"width_{key}"] = hi - lo
        rows.append(r)

    out = pd.DataFrame(rows).sort_values(["network", "label"])
    out.to_csv(args.out, sep="\t", index=False)

    print(f"n cells: {len(out)}  (folds/cell: {out.n_folds.unique()}, "
          f"repeats/cell: {out.n_repeats.unique()})")
    print()
    print("=== H1 verdict: arms whose 95% CI excludes 0 (rule: ALL cells) ===")
    for key, lab in (("fold", "fold-level percentile (current published)"),
                     ("repeat_pct", "repeat-level percentile (pre-reg wording)"),
                     ("repeat_t", "repeat-level t"),
                     ("repeat_bca", "repeat-level BCa")):
        n = int(out[f"excl0_{key}"].sum())
        print(f"  {lab:<44} {n}/{len(out)}"
              + ("   <-- H1 holds" if n == len(out) else ""))
    print()
    print(f"  mean CI width: fold {out.width_fold.mean():.4f} | "
          f"repeat_pct {out.width_repeat_pct.mean():.4f} | "
          f"repeat_t {out.width_repeat_t.mean():.4f} | "
          f"repeat_bca {out.width_repeat_bca.mean():.4f}")

    fails = out[~out.excl0_repeat_pct]
    if len(fails):
        print("\n  arms whose repeat-percentile CI contains 0:")
        print(fails[["network", "label", "mean", "ci_repeat_pct_lo",
                     "ci_repeat_pct_hi"]].to_string(index=False))
    else:
        print("\n  No arm's repeat-percentile CI contains 0 -> the "
              "pre-registered rule is met.")

    neg = out[out["mean"] < 0]
    if len(neg):
        print("\n  negative-mean arms (for the record):")
        print(neg[["network", "label", "mean", "ci_repeat_pct_lo",
                   "ci_repeat_pct_hi"]].to_string(index=False))

    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
