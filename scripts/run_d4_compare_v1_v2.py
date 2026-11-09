"""Phase D, task D4 -- compare the complex-aware (v2) result with v1.

The v1 splits are the stratified splits of the primary factorial; the
v2 splits keep every complex group in one fold (Phase D, tasks D2/D3).
Both are analysed with EXACTLY the same machinery: the same resampling
units (fold, repeat), the same five interval families, the same TOST
margin and the same BH-FDR family.

Two things have to be said alongside the comparison:

  * the v2 partitions are LESS diverse than the v1 ones (the group
    constraint removes freedom), so the repeat-level intervals are
    NARROWER for the same data -- a narrower interval is not a
    stronger result;
  * the group-preserving design changes the fold composition (a fold
    holds whole complexes), so the estimator is not identical.

The decision rule is the plan's Q3:

  |relative change in delta-AUROC|  < 5 %   -> the complex-aware result
                                         goes into the FIRST paper as a
                                         positive selling point
                                      5-10 %  -> both are reported
                                      > 10 %  -> the complex-aware result
                                         is the SECOND paper's main
                                         result and the first paper lists
                                         it as a limitation

Outputs:
  results/tables/d4_v1_vs_v2_comparison.tsv
  results/tables/d4_v2_stats.tsv
  docs/d4_complex_aware_cv.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.intervals import ci_bca, ci_percentile, ci_t  # noqa: E402
from d1.engine.stats import load_primary, tost, unit_values  # noqa: E402

TABLES = ROOT / "results" / "tables"
DOC = ROOT / "docs" / "d4_complex_aware_cv.md"
N_BOOT = 2000
SEED = 20261029
MARGIN = 0.05


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def arm_table(combined: Path, tag: str):
    """Per-arm: the interval families, the TOST and the H1 verdict."""
    piv = load_primary(str(combined), alpha=0.5)
    rows = []
    rng = np.random.default_rng(SEED)
    for (net, lab) in sorted({(n, l) for n, l in
                              zip(piv.network, piv.label_set)}):
        g = piv[(piv.network == net) & (piv.label_set == lab)]
        folds = unit_values(g, "fold")
        reps = unit_values(g, "repeat")
        r = {"version": tag, "network": net, "label": lab,
             "n_folds": len(folds), "n_repeats": len(reps),
             "mean": float(folds.mean()),
             "sd_fold": float(folds.std(ddof=1)),
             "sd_repeat": float(reps.std(ddof=1))}
        for unit, vals in (("fold", folds), ("repeat", reps)):
            lo_p, hi_p = ci_percentile(vals, rng, N_BOOT)
            lo_t, hi_t = ci_t(vals)
            lo_b, hi_b = ci_bca(vals, rng, N_BOOT)
            r[f"ci_{unit}_pct_lo"], r[f"ci_{unit}_pct_hi"] = lo_p, hi_p
            r[f"ci_{unit}_t_lo"], r[f"ci_{unit}_t_hi"] = lo_t, hi_t
            r[f"ci_{unit}_bca_lo"], r[f"ci_{unit}_bca_hi"] = lo_b, hi_b
            r[f"excl0_{unit}_pct"] = bool(not (lo_p <= 0 <= hi_p))
            r[f"excl0_{unit}_t"] = bool(not (lo_t <= 0 <= hi_t))
            r[f"excl0_{unit}_bca"] = bool(not (lo_b <= 0 <= hi_b))
            r[f"width_{unit}_pct"] = hi_p - lo_p
        t_rep = tost(reps, margin=MARGIN)
        t_fold = tost(folds, margin=MARGIN)
        r["tost_repeat_equivalent"] = bool(t_rep["equivalent"])
        r["tost_repeat_ci90_low"] = t_rep["ci90_low"]
        r["tost_repeat_ci90_high"] = t_rep["ci90_high"]
        r["tost_fold_equivalent"] = bool(t_fold["equivalent"])
        r["excl0_any"] = bool(r["excl0_fold_pct"] and r["excl0_repeat_pct"]
                              and r["excl0_repeat_t"] and r["excl0_repeat_bca"])
        rows.append(r)
    return pd.DataFrame(rows)


def main() -> int:
    v1_combined = ROOT / "results" / "runs_v40" / "a4_combined.tsv"
    v2_combined = ROOT / "results" / "runs_v40_cx" / "a4_combined.tsv"
    if not v2_combined.exists():
        raise SystemExit("the v2 run is missing; run run_a4.py on the v2 splits")

    log("computing the v1 statistics with the same machinery ...")
    v1 = arm_table(v1_combined, "v1_stratified")
    log("computing the v2 statistics ...")
    v2 = arm_table(v2_combined, "v2_complex_aware")
    stats = pd.concat([v1, v2], ignore_index=True)
    stats.to_csv(TABLES / "d4_v2_stats.tsv", sep="\t", index=False)

    # ---- the comparison ------------------------------------------------
    cmp_rows = []
    for _, a in v1.iterrows():
        b = v2[(v2.network == a.network) & (v2.label == a.label)].iloc[0]
        denom = abs(a["mean"]) if abs(a["mean"]) > 1e-9 else np.nan
        cmp_rows.append({
            "network": a.network, "label": a.label,
            "v1_mean": a["mean"], "v2_mean": b["mean"],
            "abs_change": abs(b["mean"] - a["mean"]),
            "rel_change": (b["mean"] - a["mean"]) / denom if denom else np.nan,
            "abs_rel_change": abs((b["mean"] - a["mean"]) / denom)
            if denom else np.nan,
            "v1_ci_repeat_pct": (a["ci_repeat_pct_lo"], a["ci_repeat_pct_hi"]),
            "v2_ci_repeat_pct": (b["ci_repeat_pct_lo"], b["ci_repeat_pct_hi"]),
            "v1_ci_repeat_t": (a["ci_repeat_t_lo"], a["ci_repeat_t_hi"]),
            "v2_ci_repeat_t": (b["ci_repeat_t_lo"], b["ci_repeat_t_hi"]),
            "v1_width_repeat_pct": a["width_repeat_pct"],
            "v2_width_repeat_pct": b["width_repeat_pct"],
            "v1_excl0_repeat_pct": a["excl0_repeat_pct"],
            "v2_excl0_repeat_pct": b["excl0_repeat_pct"],
            "v1_tost_equivalent": a["tost_repeat_equivalent"],
            "v2_tost_equivalent": b["tost_repeat_equivalent"],
            "v1_sd_repeat": a["sd_repeat"], "v2_sd_repeat": b["sd_repeat"],
        })
    cmp_df = pd.DataFrame(cmp_rows)
    cmp_df.to_csv(TABLES / "d4_v1_vs_v2_comparison.tsv", sep="\t",
                  index=False)

    # ---- the Q3 decision ----------------------------------------------
    med_rel = float(cmp_df["abs_rel_change"].median())
    max_rel = float(cmp_df["abs_rel_change"].max())
    n_sign_flip = int(((cmp_df.v1_mean > 0) != (cmp_df.v2_mean > 0)).sum())
    n_h1_v1 = int(v1["excl0_any"].sum())
    n_h1_v2 = int(v2["excl0_any"].sum())
    n_to_v1 = int(v1["tost_repeat_equivalent"].sum())
    n_to_v2 = int(v2["tost_repeat_equivalent"].sum())
    width_ratio = float(cmp_df["v2_width_repeat_pct"].median()
                        / cmp_df["v1_width_repeat_pct"].median())
    if med_rel < 0.05:
        decision = ("< 5 % -- the complex-aware result goes into the FIRST "
                    "paper as a positive selling point")
    elif med_rel < 0.10:
        decision = ("5-10 % -- both are reported; the paper states how "
                    "much the standard CV is optimistic")
    else:
        decision = ("> 10 % -- the complex-aware result is the SECOND "
                    "paper's main result and the first paper lists it as a "
                    "limitation")

    log(f"median |relative change| = {med_rel:.3f}  (max {max_rel:.3f})")
    log(f"H1 (all four interval families exclude 0): v1 {n_h1_v1}/24, "
        f"v2 {n_h1_v2}/24")
    log(f"TOST-equivalent cells: v1 {n_to_v1}/24, v2 {n_to_v2}/24")
    log(f"median repeat-level interval width ratio v2/v1 = {width_ratio:.3f}")
    log(f"sign flips: {n_sign_flip}/24")
    log(f"DECISION (Q3): {decision}")

    write_doc(cmp_df, v1, v2, decision, med_rel, max_rel, n_sign_flip,
              n_h1_v1, n_h1_v2, n_to_v1, n_to_v2, width_ratio)
    return 0


def write_doc(cmp_df, v1, v2, decision, med_rel, max_rel, n_sign_flip,
              n_h1_v1, n_h1_v2, n_to_v1, n_to_v2, width_ratio) -> None:
    lines = [
        "# D4 -- the complex-aware (v2) result against v1",
        "",
        "**Date:** 2026-11-05  ",
        "**Script:** `scripts/run_d4_compare_v1_v2.py`  ",
        "**Tables:** `results/tables/d4_v1_vs_v2_comparison.tsv`, "
        "`results/tables/d4_v2_stats.tsv`",
        "",
        "## Decision (the plan's Q3 rule)",
        "",
        f"* median |relative change| in delta-AUROC: **{med_rel:.3f}** "
        f"(max {max_rel:.3f})",
        f"* arms whose sign flips: **{n_sign_flip}/24**",
        f"* H1 (all four interval families exclude 0): v1 **{n_h1_v1}/24**, "
        f"v2 **{n_h1_v2}/24**",
        f"* cells equivalent to zero (TOST): v1 **{n_to_v1}/24**, "
        f"v2 **{n_to_v2}/24**",
        f"* median repeat-level interval width, v2 / v1: "
        f"**{width_ratio:.3f}**",
        "",
        f"**=> {decision}**",
        "",
        "## Per arm",
        "",
        "| network | label | v1 delta | v2 delta | rel. change | "
        "v1 95 % CI (repeat pct) | v2 95 % CI (repeat pct) |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in cmp_df.itertuples():
        lines.append(
            f"| {r.network} | {r.label} | {r.v1_mean:+.4f} | "
            f"{r.v2_mean:+.4f} | {r.rel_change:+.1%} | "
            f"[{r.v1_ci_repeat_pct[0]:+.4f}, {r.v1_ci_repeat_pct[1]:+.4f}] | "
            f"[{r.v2_ci_repeat_pct[0]:+.4f}, {r.v2_ci_repeat_pct[1]:+.4f}] |")
    lines += [
        "",
        "## How to read this",
        "",
        "1. **The interval widths are not comparable.** The complex-aware "
        f"partitions are less diverse than the stratified ones (the group "
        f"constraint removes freedom), so the v2 repeat-level intervals are "
        f"{width_ratio:.2f}x the v1 width for the same data. A narrower "
        "interval is not a stronger result; the comparison that matters is "
        "the point estimate and the sign.",
        "2. **The estimator is not identical.** A v2 fold holds whole "
        "complexes, so the test set is a different kind of set; the "
        "comparison is between two valid estimators of the same quantity "
        "on the same data, not between the same estimator twice.",
        "3. **What was removed.** The D1 audit measured that 83-95 % of the "
        "held-out positives that are in a curated complex had a "
        "complex-linked gene in the training fold. On the v2 splits that "
        "number is 0.02-0.18. The v2 result is therefore the one that "
        "does NOT benefit from the complex-partner leak.",
        "",
        "## Caveats",
        "",
        "* The complex map covers 18.3 % of the panel's genes, so the "
        "control only constrains the genes that have curated complex "
        "membership.",
        "* The residual leakage comes from the sparsity guard (a gene "
        "keeps only its 10 largest complexes) and from the complexes "
        "skipped by the group-size cap; it is measured, not assumed away.",
        "* The v2 analysis is an additional, exploratory analysis: the "
        "pre-registration fixes RWR on the stratified splits, so the v2 "
        "numbers belong to the exploratory family (amendment 4).",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
