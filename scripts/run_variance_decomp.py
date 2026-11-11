#!/usr/bin/env python3
"""Correct variance decomposition for the D1 benchmark primary factorial.

REPLACES the broken `_variance_decomp()` in scripts/run_stats_analysis.py.

Why the old one was wrong
--------------------------
The old code did:

    fit = smf.mixedlm("observed_delta ~ 1", sub, groups=sub["fold_idx"])
    var_resid = fit.scale                 # (1) within-FOLD-INDEX residual
    var_fold  = fit.cov_re.iloc[0, 0]     # (2) variance of 5 fold intercepts
    network_var = groupby("network").mean().var()   # (3) variance of NETWORK MEANS
    label_var   = groupby("label").mean().var()     # (4) variance of LABEL MEANS
    network_label_var = 0                 # (5) interaction zeroed out
    total = var_resid + var_fold + network_var + label_var

Five defects:
  (1) `groups=fold_idx` groups the 1200 observations into only 5 buckets
      (fold 0..4). The within-bucket residual therefore contains EVERY
      between-cell difference plus the within-cell noise -- it is not a
      noise bucket.
  (2)+(3)+(4) var_resid already contains network_var and label_var, so
      summing them double-counts. This is what produced the bogus
      "residual = 52%".
  (3)+(4) the variance of network means is NOT a network variance
      component: each network mean averages over 4 labels x 50 folds, so
      it also carries the label main effect and the interaction.
  (5) the protocol explicitly asks for a network x label interaction term
      (protocol section 4.6). It was set to 0 and never added to `total`.

This script does the textbook balanced random-effects ANOVA instead:

    sigma2_within = MS_within
    sigma2_int    = (MS_int    - MS_within) / n
    sigma2_label  = (MS_label  - MS_within) / (n * a)
    sigma2_network= (MS_network- MS_within) / (n * b)

    with n = 50 folds/cell, a = 6 networks, b = 4 labels.

Two fractions are reported because the protocol's wording matters:

  * "between-arm spread"  = Var(the 24 cell means)
                              = sigma2_net + sigma2_lab + sigma2_int
                                + sigma2_within / n
    This is what protocol section 4.6 asks about ("how much of the
    between-arm spread is attributable to the network, to the labels,
    and to noise").

  * "total variance"       = sigma2_net + sigma2_lab + sigma2_int
                              + sigma2_within
    This includes every within-cell fold fluctuation and is dominated by
    it (~90%), which is why the two framings must not be mixed.

Outputs
-------
  results/tables/stats_variance_decomp.tsv     (rewritten, 4 rows + header)
  results/tables/stats_variance_components.tsv (raw variance components)
  docs/variance_decomposition.md                (written)

Usage
-----
    python scripts/run_variance_decomp.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
A4 = REPO / "results" / "runs" / "a4_combined.tsv"
OUT_TABLE = REPO / "results" / "tables" / "stats_variance_decomp.tsv"
OUT_COMP = REPO / "results" / "tables" / "stats_variance_components.tsv"
OUT_DOC = REPO / "docs" / "variance_decomposition.md"

ALPHA = 0.5


def set_source(a4: str | None) -> None:
    """Point the script at a specific a4_combined.tsv.

    The primary factorial was migrated to the 40-repeat design (work plan
    A2, amendment 2 A2.3d), so the variance decomposition has to be
    computed on the SAME data as the confidence intervals and the H1
    verdict -- otherwise the manuscript quotes between-arm spreads from
    one design and intervals from another.
    """
    global A4
    if a4:
        A4 = Path(a4)
        if not A4.is_absolute():
            A4 = REPO / A4


def load_cell_means() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (cell means 6x4, per-fold deltas long)."""
    df = pd.read_csv(A4, sep="\t")
    # degree rows carry alpha=NaN (degree does not depend on alpha)
    df = df[(df["method"] == "degree") | (df["alpha"] == ALPHA)]
    sub = df[df["method"].isin(["rwr", "degree"])]
    piv = sub.pivot_table(index=["network", "label_set", "repeat", "fold"],
                          columns="method", values="auroc").reset_index()
    piv = piv.dropna(subset=["rwr", "degree"])
    piv["delta"] = piv["rwr"] - piv["degree"]
    cell = piv.groupby(["network", "label_set"])["delta"]
    M = cell.mean().unstack()
    return M, piv


def anova(M: pd.DataFrame, piv: pd.DataFrame) -> dict:
    grand = M.values.mean()
    a, b = M.shape
    n = int(piv.groupby(["network", "label_set"]).size().min())
    r = M.mean(axis=1).values
    c = M.mean(axis=0).values
    resid = M.values - r[:, None] - c[None, :] + grand

    ss_net = n * ((r - grand) ** 2).sum()
    ss_lab = n * ((c - grand) ** 2).sum()
    ss_int = n * (resid ** 2).sum()

    # within-cell: sum over cells of sum over folds of (x - cell mean)^2
    g = piv.groupby(["network", "label_set"])["delta"]
    ss_within = ((g.transform("mean") - g.transform("mean")) ** 2).sum()  # placeholder
    ss_within = 0.0
    for (nw, lb), gg in piv.groupby(["network", "label_set"]):
        ss_within += ((gg["delta"] - gg["delta"].mean()) ** 2).sum()

    df_net, df_lab, df_int = a - 1, b - 1, (a - 1) * (b - 1)
    df_within = a * b * (n - 1)

    ms_net = ss_net / df_net
    ms_lab = ss_lab / df_lab
    ms_int = ss_int / df_int
    ms_within = ss_within / df_within

    s2_within = ms_within
    s2_int = max((ms_int - ms_within) / n, 0.0)
    s2_lab = max((ms_lab - ms_within) / (n * a), 0.0)
    s2_net = max((ms_net - ms_within) / (n * b), 0.0)

    # between-arm spread = variance of the 24 cell means
    s2_between = s2_net + s2_lab + s2_int + s2_within / n
    s2_total = s2_net + s2_lab + s2_int + s2_within

    out = {
        "n": n, "a": a, "b": b,
        "ss_net": ss_net, "ss_lab": ss_lab, "ss_int": ss_int,
        "ss_within": ss_within,
        "df_net": df_net, "df_lab": df_lab, "df_int": df_int,
        "df_within": df_within,
        "ms_net": ms_net, "ms_lab": ms_lab, "ms_int": ms_int,
        "ms_within": ms_within,
        "sigma2_network": s2_net, "sigma2_label": s2_lab,
        "sigma2_interaction": s2_int, "sigma2_within": s2_within,
        "sigma2_between_arm_spread": s2_between,
        "sigma2_total": s2_total,
        "F_network": ms_net / ms_within, "p_network": 1 - stats.f.cdf(ms_net / ms_within, df_net, df_within),
        "F_label": ms_lab / ms_within, "p_label": 1 - stats.f.cdf(ms_lab / ms_within, df_lab, df_within),
        "F_interaction": ms_int / ms_within, "p_interaction": 1 - stats.f.cdf(ms_int / ms_within, df_int, df_within),
    }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a4", default=None,
                    help="path to a4_combined.tsv (default "
                         "results/runs/a4_combined.tsv; use "
                         "results/runs_v40/a4_combined.tsv for the 40-repeat "
                         "primary design)")
    args = ap.parse_args()
    set_source(args.a4)

    M, piv = load_cell_means()
    v = anova(M, piv)

    print("=== cell means (Delta AUROC), alpha=0.5 ===")
    print(M.round(4).to_string())
    print()
    print(f"grid: {v['a']} networks x {v['b']} labels x {v['n']} folds  "
          f"(balanced: {int(piv.groupby(['network','label_set']).size().min())}"
          f"..{int(piv.groupby(['network','label_set']).size().max())})")

    comps = [
        ("network (between)", v["sigma2_network"], v["sigma2_network"]),
        ("label (between)", v["sigma2_label"], v["sigma2_label"]),
        ("network x label (interaction)", v["sigma2_interaction"], v["sigma2_interaction"]),
        ("fold (within cell)", v["sigma2_within"] / v["n"], v["sigma2_within"]),
    ]
    tot_b = v["sigma2_between_arm_spread"]
    tot_t = v["sigma2_total"]

    print("\n=== variance components ===")
    print(f"{'source':<34}{'sigma^2':>14}{'% between-arm':>16}{'% total':>12}")
    for name, s2_b, s2_t in comps:
        print(f"{name:<34}{s2_b:>14.6g}{100*s2_b/tot_b:>15.2f}%{100*s2_t/tot_t:>11.2f}%")
    print(f"{'-- between-arm spread':<34}{tot_b:>14.6g}")
    print(f"{'-- total variance':<34}{tot_t:>14.6g}")
    print("   NOTE: in the 'between-arm' column the noise row is sigma2_within/n")
    print("         (the noise contribution TO A CELL MEAN), not sigma2_within.")

    print("\n=== F tests (vs within-cell noise) ===")
    for k in ("network", "label", "interaction"):
        F = v[f"F_{k}"]
        p = v[f"p_{k}"]
        star = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
        print(f"  {k:<14} F = {F:>8.2f}   p = {p:.3g}  {star}")

    # ---- write tables ----
    rows = []
    for name, s2_b, s2_t in comps:
        rows.append({
            "variance_source": name,
            "sigma2": s2_b,
            "fraction_of_between_arm_spread": s2_b / tot_b,
            "fraction_of_total_variance": s2_t / tot_t,
        })
    df_out = pd.DataFrame(rows)
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(OUT_TABLE, sep="\t", index=False, lineterminator="\n")
    print(f"\nwrote {OUT_TABLE}")

    pd.DataFrame([{
        "grid": f"{v['a']}x{v['b']}x{v['n']}",
        "sigma2_network": v["sigma2_network"],
        "sigma2_label": v["sigma2_label"],
        "sigma2_interaction": v["sigma2_interaction"],
        "sigma2_within_fold": v["sigma2_within"],
        "sigma2_between_arm_spread": tot_b,
        "sigma2_total": tot_t,
        "F_network": v["F_network"], "p_network": v["p_network"],
        "F_label": v["F_label"], "p_label": v["p_label"],
        "F_interaction": v["F_interaction"], "p_interaction": v["p_interaction"],
        "df_within": v["df_within"],
    }]).to_csv(OUT_COMP, sep="\t", index=False, lineterminator="\n")
    print(f"wrote {OUT_COMP}")

    # ---- write the doc ----
    def pct(x):
        return f"{100*x:.2f}%"

    doc = f"""# Variance decomposition (corrected)

**Date:** regenerated by `scripts/run_variance_decomp.py`
**Data:** `{A4.relative_to(REPO).as_posix()}`, alpha = 0.5, {v['a']} networks x {v['b']} labels x {v['n']} folds (balanced)
**Replaces:** the broken `_variance_decomp()` in `scripts/run_stats_analysis.py`

## Why this file exists

The previous `results/tables/stats_variance_decomp.tsv` reported a "residual"
fraction of 51.8% and was wrong. The old code fitted
`mixedlm("observed_delta ~ 1", groups=fold_idx)` with only five groups
(fold 0..4), so `fit.scale` contained every between-cell difference plus
the within-cell noise; it then added `network_var` and `label_var`
(which were variances of group MEANS, themselves already inside that
residual) into the total, double-counting them; and it set the
protocol-required network x label interaction term to exactly 0.

## Result

Balanced random-effects ANOVA on the {v['a']}x{v['b']}x{v['n']} cell means.

Two fractions are reported because they answer **two different questions**.

| source | sigma^2 | % of between-arm spread | % of total variance |
|---|---|---|---|
| network (between) | {v['sigma2_network']:.6g} | {pct(v['sigma2_network']/tot_b)} | {pct(v['sigma2_network']/tot_t)} |
| label (between) | {v['sigma2_label']:.6g} | {pct(v['sigma2_label']/tot_b)} | {pct(v['sigma2_label']/tot_t)} |
| network x label (interaction) | {v['sigma2_interaction']:.6g} | {pct(v['sigma2_interaction']/tot_b)} | {pct(v['sigma2_interaction']/tot_t)} |
| fold (within cell) | {v['sigma2_within']:.6g} | {pct(v['sigma2_within']/v['n']/tot_b)} | {pct(v['sigma2_within']/tot_t)} |
| **between-arm spread** | **{tot_b:.6g}** | 100% | {pct(tot_b/tot_t)} |
| **total variance** | **{tot_t:.6g}** | | 100% |

(In the "between-arm spread" column the noise row is sigma2_within / n --
the noise contribution TO A CELL MEAN -- because averaging n folds shrinks
it by n.)

F tests against the within-cell noise:
- network: F = {v['F_network']:.2f}, p = {v['p_network']:.3g}
- label: F = {v['F_label']:.2f}, p = {v['p_label']:.3g}
- network x label: F = {v['F_interaction']:.2f}, p = {v['p_interaction']:.3g}

## Reading -- and why the old "52% noise" was wrong

There are **two separate problems** and the old table conflated them.

**(1) What makes the 24 cell means differ from each other?**
(between-arm spread, the thing protocol section 4.6 asks about)

- network x label interaction: {pct(v['sigma2_interaction']/tot_b)}
- label main effect: {pct(v['sigma2_label']/tot_b)}
- network main effect: {pct(v['sigma2_network']/tot_b)}
- fold-to-fold sampling noise: **{pct(v['sigma2_within']/v['n']/tot_b)}**

So the between-arm spread is essentially **label and label x network
structure**. Sampling noise contributes almost nothing to it.

**(2) How precise is each cell's own estimate?**
(total variance of the raw per-fold measurements)

- fold-to-fold noise: {pct(v['sigma2_within']/tot_t)}
- network x label interaction: {pct(v['sigma2_interaction']/tot_t)}
- label: {pct(v['sigma2_label']/tot_t)}
- network: {pct(v['sigma2_network']/tot_t)}

This is what sets the width of every confidence interval.

**The old table said "noise ~52%". The corrected numbers are
{pct(v['sigma2_within']/v['n']/tot_b)} (of between-arm spread) or
{pct(v['sigma2_within']/tot_t)} (of total variance), depending on which
question is being asked.** Most of what was called "residual" was the
label main effect and the network x label interaction being
double-counted into a residual bucket.

## Consequence for the noise-reduction plan

- To reduce **(1)** the between-arm spread -> add more and more *independent*
  label sets. Adding ranking models cannot touch it.
- To reduce **(2)** the CI width -> add folds / repeats (the noise is
  sampling error of the estimator, not a defect of the method). Adding
  ranking models cannot touch it either.
- An ensemble of ranking methods reduces neither: it changes the
  estimand rather than the precision, and it is not the
  pre-registered quantity.
"""
    OUT_DOC.write_text(doc, encoding="utf-8")
    print(f"wrote {OUT_DOC}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
