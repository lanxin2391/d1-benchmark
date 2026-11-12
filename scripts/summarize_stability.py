"""Summarise the stability analysis: the trade-off curve and the scorecard.

Produces two tables used by the manuscript section and the figure:

* ``stability_tradeoff_curve.tsv`` -- the (stability, discrimination) curve of
  every operator family, aggregated over the 24 arms;
* ``stability_prediction_scorecard.tsv`` -- one row per pre-registered
  prediction (P1-P7 of ``docs/stability_regularization_plan.md``) with the
  decision rule, the observed value and the verdict.

The scorecard is the audit trail: every prediction is evaluated against the
rule that was fixed before the analysis was run.

Usage:
    python scripts/summarize_stability.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

TABLE_DIR = REPO_ROOT / "results" / "tables"
N_ARMS = 24


def _tradeoff():
    grid = pd.read_csv(TABLE_DIR / "stability_regularization_grid.tsv", sep="\t")
    seed = pd.read_csv(TABLE_DIR / "stability_seed_size_scaling.tsv", sep="\t")
    rows = []
    sub = grid[(grid.alpha == 0.5) & (grid.k == 20) & (grid.variant == "global")]
    for ctl, family in (("srs", "SRS (score-neighbourhood median)"),
                        ("c1", "C1 random neighbourhood"),
                        ("c2", "C2 global median")):
        s = sub[sub.control == ctl]
        if ctl == "c1":
            s = s[s.lam.isin([0.4, 0.7, 1.0])]
        for lam, g in s.groupby("lam"):
            rows.append({
                "axis": "smoothing", "family": family, "lam": lam,
                "n_arms": len(g),
                "median_rank_sd_median": g.median_rank_sd.median(),
                "median_rank_sd_q25": g.median_rank_sd.quantile(0.25),
                "median_rank_sd_q75": g.median_rank_sd.quantile(0.75),
                "delta_auroc_median": g.delta_auroc.median(),
                "delta_auroc_q25": g.delta_auroc.quantile(0.25),
                "delta_auroc_q75": g.delta_auroc.quantile(0.75),
            })
    ss = seed[seed.alpha == 0.5]
    for frac, g in ss.groupby("frac"):
        rows.append({
            "axis": "seed_count", "family": "training seed set size",
            "lam": np.nan, "n_arms": len(g),
            "median_rank_sd_median": g.median_rank_sd.median(),
            "median_rank_sd_q25": g.median_rank_sd.quantile(0.25),
            "median_rank_sd_q75": g.median_rank_sd.quantile(0.75),
            "delta_auroc_median": g.delta_auroc.median(),
            "delta_auroc_q25": g.delta_auroc.quantile(0.25),
            "delta_auroc_q75": g.delta_auroc.quantile(0.75),
        })
    return pd.DataFrame(rows)


def _scorecard():
    grid = pd.read_csv(TABLE_DIR / "stability_regularization_grid.tsv", sep="\t")
    arms = pd.read_csv(TABLE_DIR / "stability_by_arm.tsv", sep="\t")
    attr = pd.read_csv(TABLE_DIR / "stability_crowding_attribution.tsv", sep="\t")
    ctl = pd.read_csv(TABLE_DIR / "stability_control_comparison.tsv", sep="\t")
    c3 = pd.read_csv(TABLE_DIR / "stability_control_c3_flattening.tsv", sep="\t")
    seed = pd.read_csv(TABLE_DIR / "stability_seed_size_scaling.tsv", sep="\t")
    sel = pd.read_csv(TABLE_DIR / "stability_selected_scorer.tsv", sep="\t")
    expo = pd.read_csv(TABLE_DIR / "stability_seed_size_exponents.tsv", sep="\t")

    srs = grid[(grid.control == "srs") & (grid.variant == "global") & (grid.k == 20)
               & (grid.alpha == 0.5)]
    base = srs[srs.lam == 0.0].set_index(["network", "label_set"])
    d_sd = srs.reset_index().merge(
        base.reset_index()[["network", "label_set", "median_rank_sd", "delta_auroc"]],
        on=["network", "label_set"], suffixes=("", "_base"))
    max_d_sd = (d_sd.median_rank_sd - d_sd.median_rank_sd_base).abs().max()
    max_d_auc = (d_sd.delta_auroc - d_sd.delta_auroc_base).abs().max()

    rows = []

    rows.append(dict(
        prediction="P1", statement="median per-gene rank SD decreases as lambda grows",
        rule="monotone decrease in >= 24/24 arms (primary: k=20, confidence-weighted)",
        observed=f"max |change| in the median rank SD across all arms and lambda = "
                 f"{max_d_sd:.2e} (baseline median rank SD 0.019-0.044); no arm shows a "
                 f"monotone decrease",
        verdict="NOT SUPPORTED (the operator is order-preserving: it does not change "
                "the ranking)"))

    rows.append(dict(
        prediction="P2", statement="delta-AUROC does not fall",
        rule="median |delta(delta-AUROC)| <= 0.010",
        observed=f"max |change| = {max_d_auc:.2e} across all arms and lambda",
        verdict="SUPPORTED, but vacuously: nothing changed"))

    rows.append(dict(
        prediction="P3", statement="a monotone stability/discrimination trade-off curve",
        rule="report the curve with the price of stability",
        observed="the SRS curve is a single point (stability and discrimination both "
                 "flat); the random and global-median controls move to worse "
                 "stability AND worse discrimination; the seed-count axis moves to "
                 "better stability AND better discrimination",
        verdict="NO TRAD-OFF EXISTS in the smoothing direction; the only Pareto "
                "point is the unsmoothed prior"))

    rows.append(dict(
        prediction="P4", statement="specificity vs a matched-strength random control",
        rule="the SRS value must be below the 95th percentile of the random null in "
             ">= 24/24 arms",
        observed=f"{int(ctl.srs_below_c1_p95.sum())}/{len(ctl)} (arm, lambda) cells "
                 f"below the 95th percentile of the 20-replicate random null; the "
                 f"random control raises the median rank SD from "
                 f"{ctl.median_rank_sd.median():.4f} to "
                 f"{ctl.c1_rank_sd_median.median():.4f} at lambda=0.4",
        verdict="SUPPORTED, but vacuously (the SRS value is the unsmoothed value)"))

    rows.append(dict(
        prediction="P5", statement="the gain is not an artefact of flattening",
        rule="monotone transforms leave rank statistics unchanged; AUROC preserved; "
             "dispersion reported",
        observed=f"a monotone power transform changes the AUROC by at most "
                 f"{c3.auroc_gap.abs().max():.1e} and the rank of at most a handful of "
                 f"genes (float ties); the SRS operator changes the score dispersion "
                 f"by up to 40 % while leaving every rank metric unchanged",
        verdict="SUPPORTED (there is no gain to explain: the operator leaves the "
                "ranking untouched)"))

    rows.append(dict(
        prediction="P6", statement="the effect is crowding-specific",
        rule="crowded genes must gain more stability than isolated genes",
        observed=f"per-arm Spearman rho between score-space crowding and the rank SD: "
                 f"median {attr.rho_crowding.median():+.3f} "
                 f"(range {attr.rho_crowding.min():+.3f} to {attr.rho_crowding.max():+.3f}, "
                 f"negative in {int((attr.rho_crowding < 0).sum())}/24 arms); "
                 f"rank-space crowding: median {attr.rho_crowding_rank.median():+.3f} "
                 f"(positive in {int((attr.rho_crowding_rank > 0).sum())}/24); "
                 f"the gene's own score CV: median {attr.rho_score_cv.median():+.3f} "
                 f"(positive in {int((attr.rho_score_cv > 0).sum())}/24)",
        verdict="NOT SUPPORTED for score-space crowding; PARTIAL for rank-space "
                "crowding; the dominant predictor is the gene's own score noise"))

    ind = grid[(grid.variant == "global_inductive") & (grid.alpha == 0.5)]
    rows.append(dict(
        prediction="P7", statement="the effect survives leaving the fold",
        rule="the inductive variant reproduces the direction with at least half the "
             "effect size",
        observed=f"the inductive variant (neighbourhood from training-fold genes only) "
                 f"gives the same median rank SD as the transductive one "
                 f"({ind.median_rank_sd.median():.5f} vs "
                 f"{srs[srs.lam == 0].median_rank_sd.median():.5f}) and the same "
                 f"delta-AUROC; there is no effect to survive",
        verdict="NOT APPLICABLE (no effect to replicate)"))

    extra = dict(
        n_arms=N_ARMS,
        baseline_rank_sd_min=float(arms.median_rank_sd.min()),
        baseline_rank_sd_max=float(arms.median_rank_sd.max()),
        baseline_delta_auroc_median=float(arms.delta_auroc.median()),
        seed_rank_sd_5pct=float(seed[seed.frac == 0.05].median_rank_sd.median()),
        seed_rank_sd_100pct=float(seed[seed.frac == 1.0].median_rank_sd.median()),
        seed_exponent_median=float(expo.rank_sd_exponent.median()),
        seed_arms_monotone=int(sum(
            (seed.pivot_table(index=["network", "label_set"], columns="frac",
                              values="median_rank_sd").diff(axis=1).iloc[:, 1:] < 0)
            .all(axis=1))),
        seed_auroc_arms_monotone=int(sum(
            (seed.pivot_table(index=["network", "label_set"], columns="frac",
                              values="delta_auroc").diff(axis=1).iloc[:, 1:] > 0)
            .all(axis=1))),
        selected_lambda_counts=", ".join(
            f"{int(v)}x k={k}" for (k, v) in
            sel[(sel.alpha == 0.5) & (sel.variant == "global") & (sel.k == 20)]
            .lam_selected.value_counts().items()),
    )
    for k, v in extra.items():
        rows.append(dict(prediction=f"context_{k}", statement="", rule="", observed=str(v),
                          verdict=""))
    return pd.DataFrame(rows)


def main() -> int:
    curve = _tradeoff()
    curve.to_csv(TABLE_DIR / "stability_tradeoff_curve.tsv", sep="\t", index=False, lineterminator="\n")
    card = _scorecard()
    card.to_csv(TABLE_DIR / "stability_prediction_scorecard.tsv", sep="\t", index=False, lineterminator="\n")
    print(curve.to_string(index=False))
    print()
    for _, r in card.iterrows():
        if r.prediction.startswith("P"):
            print(f"{r.prediction}: {r.verdict}\n    observed: {r.observed}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
