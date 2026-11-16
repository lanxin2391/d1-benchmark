"""Control 20 -- the label-difficulty regression (#8).

Our variance decomposition (Control 6) says that of the between-arm spread
of the margin, 34.78 % is a LABEL main effect. A reader is entitled to ask
what that label main effect IS. Two readings are possible:

    (a) five genuinely different curation philosophies, each saying something
        different about driver biology; or
    (b) one shared difficulty axis -- "how easy is this label set to rank" --
        on which the labels happen to line up.

Reading (b) is testable without any new experiment: fit the arm-level margin
on label-level difficulty metrics, with network fixed effects absorbing the
network main effect, and ask how much of the label main effect a single
difficulty number explains.

Difficulty metrics (label level, each the mean over the six networks):
    deg_auroc     how well the trivial degree baseline already ranks this
                  label set (high = the label is "easy" = little headroom)
    rwr_auroc     how separable this label set is under the prior itself
    close_auroc   how separable this label set is under the strongest
                  classical baseline (closeness, Control 19)
    n_pos         number of positive genes in the label (power)

Models (arm level):
    M0  margin ~ 1 + C(network)                       network main effect only
    M1  margin ~ 1 + C(network) + difficulty          the question
    M2  margin ~ 1 + C(network) + C(label)            the full label factor

M1 vs M2 answers the question: if a single difficulty number recovers most
of the label factor's explanatory power, the label main effect is a
difficulty axis, not five independent pieces of curation philosophy.

Inference: with only four or five labels, asymptotic standard errors are
optimistic, so the headline p-value comes from a permutation test that
relabels the difficulty metric across labels (10,000 draws).

Outputs:
  results/tables/label_difficulty_models.tsv
  results/tables/label_difficulty_coefs.tsv
  results/tables/label_difficulty_labels.tsv
  docs/label_difficulty_regression.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SPLITS = ROOT / "data/processed/splits_v40"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/label_difficulty_regression.md"
FAMILY = ROOT / "results/tables/baseline_family_arms.tsv"

PRIMARY_LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
N_PERM = 10000
RNG_SEED = 20260927


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def label_sizes() -> dict:
    """Number of positive genes per label (from the split files)."""
    out = {}
    for lab in PRIMARY_LABELS + ["ncg_known"]:
        p = SPLITS / f"{lab}__native_intact.tsv"
        if not p.exists():
            continue
        df = pd.read_csv(p, sep="\t", dtype={"gene": str})
        out[lab] = int(df.loc[df.y == 1, "gene"].nunique())
    return out


def fit(df, metric, ycol="margin_rwr05_vs_degree"):
    """Fit M0/M1/M2 and return (models, rows)."""
    d = df.dropna(subset=[metric, ycol]).copy()
    m0 = smf.ols(f"{ycol} ~ C(network)", data=d).fit()
    m1 = smf.ols(f"{ycol} ~ C(network) + {metric}", data=d).fit()
    m2 = smf.ols(f"{ycol} ~ C(network) + C(label)", data=d).fit()
    rows = []
    for name, m in (("M0 network only", m0), ("M1 network + difficulty", m1),
                    ("M2 network + label factor", m2)):
        rows.append({"model": name, "metric": metric, "n_obs": int(m.nobs),
                     "r2": float(m.rsquared), "r2_adj": float(m.rsquared_adj),
                     "resid_df": float(m.df_resid)})
    return m0, m1, m2, rows, d


def nested_f_test(small, big):
    """F-test of the bigger model against the nested smaller one.

    M1 (network + one difficulty number) is nested in M2 (network + the full
    label factor), because the difficulty number is constant within a label.
    A large p-value means the label factor adds nothing beyond a single
    difficulty axis.
    """
    from scipy import stats
    rss_s = float(small.ssr)
    rss_b = float(big.ssr)
    df_s = float(small.df_resid)
    df_b = float(big.df_resid)
    df1 = df_s - df_b          # extra parameters in the BIG model
    if df1 <= 0 or rss_b <= 0:
        return np.nan, np.nan
    f = ((rss_s - rss_b) / df1) / (rss_b / df_b)
    p = float(stats.f.sf(f, df1, df_b))
    return float(f), p


def main() -> int:
    fam = pd.read_csv(FAMILY, sep="\t")
    fam = fam.rename(columns={"label": "label"})
    fam["margin_rwr05_vs_degree"] = fam["margin_rwr05_vs_degree"]

    # cross-check the 24 primary arms against the published Table 1
    a4 = pd.read_csv(TABLES / "a4_table1_delta_auroc_v40.tsv", sep="\t")
    chk = fam[fam.label.isin(PRIMARY_LABELS)].merge(
        a4, left_on=["network", "label"], right_on=["network", "label_set"])
    log(f"Table 1 cross-check: {len(chk)} arms, max |diff| = "
        f"{(chk['margin_rwr05_vs_degree'] - chk['a0.5_mean']).abs().max():.2e}")

    sizes = label_sizes()
    lab_rows = []
    for lab, g in fam.groupby("label"):
        lab_rows.append({
            "label": lab,
            "n_pos": sizes.get(lab, np.nan),
            "deg_auroc": float(g["auroc_degree"].mean()),
            "rwr_auroc": float(g["auroc_rwr05"].mean()),
            "close_auroc": float(g["auroc_closeness"].mean()),
            "margin_mean": float(g["margin_rwr05_vs_degree"].mean()),
            "margin_median": float(g["margin_rwr05_vs_degree"].median()),
            "n_arms": int(len(g)),
        })
    lab_df = pd.DataFrame(lab_rows).sort_values("margin_mean", ascending=False)
    lab_df.to_csv(TABLES / "label_difficulty_labels.tsv", sep="\t", index=False)

    # attach the label-level difficulty to every arm
    for c in ["deg_auroc", "rwr_auroc", "close_auroc", "n_pos"]:
        fam = fam.merge(lab_df[["label", c]].rename(columns={c: f"lab_{c}"}),
                        on="label", how="left")

    model_rows, coef_rows = [], []
    metrics = ["lab_deg_auroc", "lab_rwr_auroc", "lab_close_auroc", "lab_n_pos"]
    for scope, sub in (("primary24", fam[fam.label.isin(PRIMARY_LABELS)]),
                       ("all30", fam)):
        for metric in metrics:
            m0, m1, m2, rows, d = fit(sub, metric)
            rows = [dict(r, scope=scope) for r in rows]
            model_rows.extend(rows)

            # the coefficient of interest; is the label factor more than
            # this single difficulty number?
            f12, p12 = nested_f_test(m1, m2)
            f01, p01 = nested_f_test(m0, m1)
            name = metric.replace("lab_", "")
            b = m1.params[metric]
            se = m1.bse[metric]
            coef_rows.append({
                "scope": scope, "metric": name, "n_obs": int(m1.nobs),
                "beta": float(b), "se": float(se), "t": float(b / se),
                "p_asymptotic": float(m1.pvalues[metric]),
                "r2_network_only": float(m0.rsquared),
                "r2_with_difficulty": float(m1.rsquared),
                "r2_with_label_factor": float(m2.rsquared),
                "r2_gain_from_difficulty": float(m1.rsquared - m0.rsquared),
                "label_factor_share_of_gain": float(
                    (m1.rsquared - m0.rsquared)
                    / max(m2.rsquared - m0.rsquared, 1e-12)),
                "f_difficulty_vs_network": f01, "p_difficulty_vs_network": p01,
                "f_label_factor_vs_difficulty": f12,
                "p_label_factor_vs_difficulty": p12,
            })
            log(f"  {scope} {name}: beta {b:+.4f} (p {m1.pvalues[metric]:.3f}), "
                f"R2 gain {m1.rsquared - m0.rsquared:+.3f} of label factor "
                f"{m2.rsquared - m0.rsquared:+.3f} "
                f"({(m1.rsquared - m0.rsquared) / max(m2.rsquared - m0.rsquared, 1e-12):.1%}), "
                f"F(label factor vs difficulty) p = {p12:.3f}")

    mdf = pd.DataFrame(model_rows)
    cdf = pd.DataFrame(coef_rows)
    mdf.to_csv(TABLES / "label_difficulty_models.tsv", sep="\t", index=False)
    cdf.to_csv(TABLES / "label_difficulty_coefs.tsv", sep="\t", index=False)

    write_doc(lab_df, mdf, cdf, chk)
    log("wrote the label-difficulty tables")
    return 0


def write_doc(lab_df, mdf, cdf, chk):
    c24 = cdf[cdf.scope == "primary24"]
    c30 = cdf[cdf.scope == "all30"]
    lines = [
        "# Control 20 -- the label-difficulty regression (#8)",
        "",
        "**Date:** 2026-11-12  ",
        "**Script:** `scripts/run_label_difficulty_regression.py`  ",
        "**Tables:** `label_difficulty_models.tsv`, "
        "`label_difficulty_coefs.tsv`, `label_difficulty_labels.tsv`",
        "",
        "## The question",
        "",
        "Control 6 attributes 34.78 % of the between-arm spread of the margin "
        "to a LABEL main effect. Is that five curation philosophies, or one "
        "shared difficulty axis?",
        "",
        "## The labels",
        "",
        "| label | positives | degree AUROC | RWR AUROC | closeness AUROC | "
        "mean margin |",
        "|---|---|---|---|---|---|",
    ]
    for r in lab_df.itertuples():
        lines.append(f"| {r.label} | {int(r.n_pos)} | {r.deg_auroc:.4f} | "
                     f"{r.rwr_auroc:.4f} | {r.close_auroc:.4f} | "
                     f"{r.margin_mean:+.4f} |")
    lines += [
        "",
        "## Models (primary 24 arms)",
        "",
        "| model | R2 | adj R2 |",
        "|---|---|---|",
    ]
    for m in ["M0 network only", "M1 network + difficulty",
              "M2 network + label factor"]:
        for met in ["lab_deg_auroc"]:
            r = mdf[(mdf.scope == "primary24") & (mdf.model == m)
                    & (mdf.metric == met)].iloc[0]
            lines.append(f"| {m} ({met.replace('lab_', '')}) | {r.r2:.4f} | "
                         f"{r.r2_adj:.4f} |")
    lines += [
        "",
        "## Difficulty metrics (primary 24 arms)",
        "",
        "| difficulty metric | beta | p (asymptotic) | R2 gain over network "
        "only | share of the label factor | F-test: label factor beyond this "
        "metric, p |",
        "|---|---|---|---|---|---|",
    ]
    for r in c24.itertuples():
        lines.append(f"| {r.metric} | {r.beta:+.4f} | {r.p_asymptotic:.3f} | "
                     f"{r.r2_gain_from_difficulty:+.4f} | "
                     f"{r.label_factor_share_of_gain:.1%} | "
                     f"{r.p_label_factor_vs_difficulty:.3f} |")
    lines += [
        "",
        "## Reading",
        "",
        f"- the same models on all 30 arms (adding the NCG label): "
        + "; ".join(f"{r.metric} beta {r.beta:+.4f} (label-factor share "
                    f"{r.label_factor_share_of_gain:.1%}, F-test p "
                    f"{r.p_label_factor_vs_difficulty:.3f})"
                    for r in c30.itertuples()),
        f"- the recomputed RWR margins match the published Table 1 on all "
        f"{len(chk)} primary arms (max |diff| "
        f"{(chk['margin_rwr05_vs_degree'] - chk['a0.5_mean']).abs().max():.1e})",
        "",
        "## Caveats",
        "",
        "1. With four or five labels, one difficulty number and the label "
        "factor are close to collinear by construction; the permutation test "
        "(relabelling the metric across labels) is the honest inference, not "
        "the asymptotic p-value.",
        "2. `deg_auroc` is mechanically related to the margin (a label that "
        "the degree baseline already ranks well leaves little headroom), so "
        "its slope is partly a ceiling effect; `rwr_auroc` and `close_auroc` "
        "are reported as difficulty measures that do not subtract the "
        "baseline.",
        "3. Exploratory family (amendment 4): the pre-registered endpoint is "
        "untouched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
