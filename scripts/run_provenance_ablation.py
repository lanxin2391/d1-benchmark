"""Provenance ablation: partition the 6 primary networks into 3
provenance classes (RNA-derived, protein-derived, curated/literature)
and test whether the contribution to cancer-driver identification is
class-invariant.

This addresses H2 of the pre-registration:
    "provenance explains a small share of ΔAUROC variance"
    if true: contributions are similar across classes (provenance-invariant)
    if false: some molecular layer genuinely encodes more driver biology

Protocol §1.2 (per-source classification):
    RNA-derived:    rna_coexp
    Protein-derived: funmap, string_phys700, intact
    Curated/literature: string_full700, reactome

Inputs:
    results/runs/a4_combined.tsv          (4800 rows: 6 nets × 4 labels ×
                                           3 alphas × 50 folds × 2 methods)
    results/tables/results_stats.tsv       (observed + statistical tests)

Outputs:
    results/tables/provenance_ablation.tsv
        per (class, label, alpha): mean ΔAUROC across the networks in
        the class + SD + n_networks.
    results/tables/provenance_summary.tsv
        per (label, alpha): class means + ANOVA / Kruskal effect size.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

RUN_A4 = ROOT / "results" / "runs" / "a4_combined.tsv"
TABLE_DIR = ROOT / "results" / "tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

# Per protocol §1.2 provenance stratification
PROVENANCE = {
    "funmap":          "protein-derived",
    "string_phys700":  "protein-derived",
    "intact":          "protein-derived",
    "rna_coexp":       "rna-derived",
    "string_full700":  "curated-literature",
    "reactome":        "curated-literature",
}
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]


def main():
    df = pd.read_csv(RUN_A4, sep="\t")
    df["repeat"] = df["repeat"].astype(int)
    df["fold"] = df["fold"].astype(int)
    # Compute delta per fold = RWR(auroc) - degree(auroc).
    # We can't set_index on alpha because degree rows have alpha=NaN.
    # Instead, drop alpha, then re-attach after merging.
    rwr = df[df.method == "rwr"].drop(columns=["alpha", "method", "n_iter",
                                              "auprc", "p_at_50", "p_at_100",
                                              "p_at_500", "n_test_pos", "n_test"])
    deg = df[df.method == "degree"].drop(columns=["alpha", "method", "n_iter",
                                              "auprc", "p_at_50", "p_at_100",
                                              "p_at_500", "n_test_pos", "n_test"])
    rwr = rwr.rename(columns={"auroc": "auroc_rwr"}).set_index(
        ["network", "label_set", "repeat", "fold"])
    deg = deg.rename(columns={"auroc": "auroc_deg"}).set_index(
        ["network", "label_set", "repeat", "fold"])
    delta = (rwr["auroc_rwr"] - deg["auroc_deg"]).dropna().reset_index(name="delta")
    # re-attach alpha from the RWR rows (one row per (network, label, alpha, rep, fold))
    alpha_lookup = (df[df.method == "rwr"]
                    [["network", "label_set", "repeat", "fold", "alpha"]]
                    .drop_duplicates())
    delta = delta.merge(alpha_lookup, on=["network", "label_set", "repeat", "fold"])

    delta["provence_class"] = delta["network"].map(PROVENANCE)

    # ---- per (class, label, alpha): mean across networks in that class ----
    print("=" * 70)
    print("Provenance ablation")
    print("=" * 70)
    grouped = (delta.groupby(["provence_class", "label_set", "alpha", "network"])
               .delta.mean().reset_index(name="delta_network_mean"))
    grouped_class = (grouped.groupby(["provence_class", "label_set", "alpha"])
                     .agg(delta_class_mean=("delta_network_mean", "mean"),
                          delta_class_sd=("delta_network_mean", "std"),
                          n_networks=("network", "count"))
                     .reset_index())
    grouped_class = grouped_class.rename(columns={"label_set": "label"})
    grouped_class.to_csv(TABLE_DIR / "provenance_ablation.tsv",
                          sep="\t", index=False)

    print("\nPer-(class, label, alpha=0.5) mean delta_AUROC:")
    pivot = grouped_class[grouped_class.alpha == 0.5].pivot_table(
        index=["provence_class"], columns="label", values="delta_class_mean")
    print(pivot.round(4).to_string())

    # ---- per (label, alpha): class means + Kruskal between classes ----
    # Test H2: do the three classes differ in their mean delta_AUROC?
    # Kruskal-Wallis on per-network class means (1 observation per network).
    summary = []
    for (lab, alpha), sub in grouped.groupby(["label_set", "alpha"]):
        if sub["provence_class"].nunique() < 2:
            summary.append({
                "label": lab, "alpha": alpha,
                "H_kruskal": np.nan, "p_kruskal": np.nan,
                "eta_squared_class": np.nan,
                "n_classes": 1,
                "n_networks": len(sub),
                **{f"{c}_mean": np.nan for c in PROVENANCE.values()},
            })
            continue
        classes = sorted(sub["provence_class"].unique())
        groups = [sub[sub["provence_class"] == c].delta_network_mean.values
                  for c in classes]
        try:
            H, p_kw = scipy_stats.kruskal(*groups)
        except Exception:
            H, p_kw = np.nan, np.nan
        # eta-squared (between-class variance / total variance)
        all_vals = sub.delta_network_mean.values
        if all_vals.var(ddof=1) > 0:
            class_means = sub.groupby("provence_class").delta_network_mean.mean()
            between_var = sum(sub.groupby("provence_class").size() *
                              (class_means - all_vals.mean())**2) / len(sub)
            eta2 = between_var / all_vals.var(ddof=1)
        else:
            eta2 = np.nan
        # per-class means for the table
        class_means = sub.groupby("provence_class").delta_network_mean.mean()
        row = {
            "label": lab, "alpha": alpha,
            "H_kruskal": float(H) if not np.isnan(H) else np.nan,
            "p_kruskal": float(p_kw) if not np.isnan(p_kw) else np.nan,
            "eta_squared_class": float(eta2) if not np.isnan(eta2) else np.nan,
            "n_classes": len(classes),
            "n_networks": len(sub),
        }
        for c in PROVENANCE.values():
            row[f"{c}_mean"] = float(class_means.get(c, np.nan))
        summary.append(row)
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(TABLE_DIR / "provenance_summary.tsv",
                       sep="\t", index=False)
    print("\nPer-(label, alpha=0.5) provenance class means + Kruskal effect size:")
    cls_cols = [f"{c}_mean" for c in PROVENANCE.values()]
    keep = ["label", "alpha", "n_classes", "n_networks",
            "H_kruskal", "p_kruskal", "eta_squared_class"] + cls_cols
    keep = [c for c in keep if c in summary_df.columns]
    print(summary_df[summary_df.alpha == 0.5][keep].round(4).to_string(
        index=False))


if __name__ == "__main__":
    main()