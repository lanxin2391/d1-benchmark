#!/usr/bin/env python3
"""Regenerate the corrected section 4.6 statistical tables.

Fixes the three defects listed in `d1/engine/stats.py` and writes:

  results/tables/stats_bootstrap_ci.tsv        per-cell CI, three units
  results/tables/stats_equivalence.tsv         per-cell TOST, three units
  results/tables/stats_ci_all_units.tsv        one row per cell, all units
  results/tables/stats_h1_verdict.tsv          the H1 decision, per unit

The primary unit is the REPEAT (the coverage simulation of work plan
A1/A2 showed it is the calibrated construction); the gene-cluster
bootstrap of pre-registration section 9.3 is reported alongside as the
conservative bound (decision of 2026-10-29).

Usage
-----
    python scripts/run_statistics.py \\
        --combined results/runs_v40/a4_combined.tsv \\
        --gene-ci results/tables/stats_bootstrap_ci_gene.tsv \\
        --suffix _v40
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from d1.engine.intervals import ci_bca, ci_percentile, ci_t  # noqa: E402
from d1.engine.stats import load_primary, tost, unit_values  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--combined", default=str(REPO / "results" / "runs_v40" /
                                             "a4_combined.tsv"))
    ap.add_argument("--gene-ci", default=str(REPO / "results" / "tables" /
                                            "stats_bootstrap_ci_gene.tsv"))
    ap.add_argument("--tables", default=str(REPO / "results" / "tables"))
    ap.add_argument("--suffix", default="_v40")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261029)
    ap.add_argument("--margin", type=float, default=0.05)
    args = ap.parse_args(argv)

    tables = Path(args.tables)
    tables.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    piv = load_primary(args.combined, alpha=0.5)
    cells = sorted({(n, l) for n, l in zip(piv.network, piv.label_set)})

    gene = None
    gp = Path(args.gene_ci)
    if gp.exists():
        gene = pd.read_csv(gp, sep="\t")
        gene = gene[gene.alpha == 0.5].set_index(["network", "label"])

    ci_rows, eq_rows, all_rows = [], [], []
    for (net, lab) in cells:
        g = piv[(piv.network == net) & (piv.label_set == lab)]
        folds = unit_values(g, "fold")
        reps = unit_values(g, "repeat")

        r = {"network": net, "label": lab, "alpha": 0.5,
             "n_folds": len(folds), "n_repeats": len(reps),
             "mean": float(folds.mean()),
             "sd_fold": float(folds.std(ddof=1)),
             "sd_repeat": float(reps.std(ddof=1))}

        for unit, vals in (("fold", folds), ("repeat", reps)):
            lo_p, hi_p = ci_percentile(vals, rng, args.n_boot)
            lo_t, hi_t = ci_t(vals)
            lo_b, hi_b = ci_bca(vals, rng, args.n_boot)
            r[f"ci_{unit}_pct_lo"], r[f"ci_{unit}_pct_hi"] = lo_p, hi_p
            r[f"ci_{unit}_t_lo"], r[f"ci_{unit}_t_hi"] = lo_t, hi_t
            r[f"ci_{unit}_bca_lo"], r[f"ci_{unit}_bca_hi"] = lo_b, hi_b
            r[f"excl0_{unit}_pct"] = bool(not (lo_p <= 0 <= hi_p))
            r[f"excl0_{unit}_t"] = bool(not (lo_t <= 0 <= hi_t))
            r[f"excl0_{unit}_bca"] = bool(not (lo_b <= 0 <= hi_b))
            r[f"width_{unit}_pct"] = hi_p - lo_p

        if gene is not None and (net, lab) in gene.index:
            gr = gene.loc[(net, lab)]
            r["ci_gene_lo"], r["ci_gene_hi"] = float(gr.ci_low), float(gr.ci_high)
            r["excl0_gene"] = bool(gr.excl0)
            r["width_gene"] = float(gr.width)

        # ---- TOST at each unit ----
        for unit, vals in (("fold", folds), ("repeat", reps)):
            t = tost(vals, margin=args.margin)
            for k in ("mean", "sd", "n", "t_low", "p_low", "t_high", "p_high",
                      "ci90_low", "ci90_high", "ci95_low", "ci95_high",
                      "equivalent", "verdict"):
                r[f"tost_{unit}_{k}"] = t[k]
        all_rows.append(r)

        # ---- long-format CI table ----
        ci_rows.append({
            "network": net, "label": lab, "alpha": 0.5,
            "n_iter": args.n_boot,
            "mean": r["mean"],
            "ci_low": r["ci_repeat_pct_lo"], "ci_high": r["ci_repeat_pct_hi"],
            "primary_interval": "repeat_percentile",
            "ci_low_fold_pct": r["ci_fold_pct_lo"],
            "ci_high_fold_pct": r["ci_fold_pct_hi"],
            "ci_low_repeat_pct": r["ci_repeat_pct_lo"],
            "ci_high_repeat_pct": r["ci_repeat_pct_hi"],
            "ci_low_repeat_t": r["ci_repeat_t_lo"],
            "ci_high_repeat_t": r["ci_repeat_t_hi"],
            "ci_low_repeat_bca": r["ci_repeat_bca_lo"],
            "ci_high_repeat_bca": r["ci_repeat_bca_hi"],
            "ci_low_gene": r.get("ci_gene_lo", np.nan),
            "ci_high_gene": r.get("ci_gene_hi", np.nan),
            "excl0_fold_pct": r["excl0_fold_pct"],
            "excl0_repeat_pct": r["excl0_repeat_pct"],
            "excl0_repeat_t": r["excl0_repeat_t"],
            "excl0_repeat_bca": r["excl0_repeat_bca"],
            "excl0_gene": r.get("excl0_gene", np.nan),
        })

        eq = tost(reps, margin=args.margin)
        eq_rows.append({
            "network": net, "label": lab, "alpha": 0.5,
            "margin": args.margin, "unit": "repeat",
            "mean": eq["mean"], "sd": eq["sd"], "n": eq["n"],
            "t_low": eq["t_low"], "p_low": eq["p_low"],
            "t_high": eq["t_high"], "p_high": eq["p_high"],
            "ci90_low": eq["ci90_low"], "ci90_high": eq["ci90_high"],
            "ci95_low": eq["ci95_low"], "ci95_high": eq["ci95_high"],
            "equivalent_to_zero": eq["equivalent"],
            "verdict": eq["verdict"],
        })

    sfx = args.suffix
    pd.DataFrame(ci_rows).to_csv(tables / f"stats_bootstrap_ci{sfx}.tsv",
                                  sep="\t", index=False)
    pd.DataFrame(eq_rows).to_csv(tables / f"stats_equivalence{sfx}.tsv",
                                 sep="\t", index=False)
    all_df = pd.DataFrame(all_rows).sort_values(["network", "label"])
    all_df.to_csv(tables / f"stats_ci_all_units{sfx}.tsv", sep="\t", index=False)

    # ---- H1 verdict per unit ----
    verdict_rows = []
    for unit in ("fold_pct", "repeat_pct", "repeat_t", "repeat_bca", "gene"):
        col = f"excl0_{unit}"
        if col not in all_df.columns:
            continue
        n = int(all_df[col].sum())
        verdict_rows.append({
            "interval": unit, "n_arms": len(all_df), "n_excl0": n,
            "rule_met": bool(n == len(all_df)),
            "rule": "pre-registered H1: the 95% CI excludes 0 in EVERY arm",
            "failing_arms": "; ".join(
                f"{r.network} x {r.label}" for r in all_df[~all_df[col]].itertuples()),
        })
    vdf = pd.DataFrame(verdict_rows)
    vdf.to_csv(tables / f"stats_h1_verdict{sfx}.tsv", sep="\t", index=False)

    print("=== H1 verdict per interval construction (40-repeat design) ===")
    for r in verdict_rows:
        mark = "MET" if r["rule_met"] else "NOT met"
        print(f"  {r['interval']:<12} {r['n_excl0']:>2}/{r['n_arms']}  {mark}"
              + (f"   failing: {r['failing_arms']}" if r["failing_arms"] else ""))
    print()
    print("=== TOST (repeat unit, margin 0.05) ===")
    for r in eq_rows:
        if not r["equivalent_to_zero"]:
            continue
        print(f"  {r['network']} x {r['label']}: mean {r['mean']:+.4f}  "
              f"90% CI [{r['ci90_low']:+.4f}, {r['ci90_high']:+.4f}]  "
              f"p_low {r['p_low']:.2e}  p_high {r['p_high']:.2e}  "
              f"-> EQUIVALENT TO ZERO (null cell)")
    print()
    print(f"wrote {tables}/stats_bootstrap_ci{sfx}.tsv")
    print(f"wrote {tables}/stats_equivalence{sfx}.tsv")
    print(f"wrote {tables}/stats_ci_all_units{sfx}.tsv")
    print(f"wrote {tables}/stats_h1_verdict{sfx}.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
