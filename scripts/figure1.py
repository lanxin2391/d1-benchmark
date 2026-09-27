"""Figure 1 — the primary visualisation for the D1 pilot.

Two panels:

* (left) Forest plot of `delta_AUROC = AUROC(RWR) - AUROC(degree)` per arm,
  with 95% CI on each arm and a fixed-effects meta mean. Arms are grouped by
  regime: the original 4-arm grid (top), D-09 temporal-drift (middle), and
  D-11 common-size-cap (bottom).
* (right) Heatmap of mean delta_AUROC for the 4 labels × 2 networks grid,
  with rows = labels and columns = networks. Two heatmaps side by side:
  uncapped vs D-11-capped.

Saves PNG to ``figures/figure1_pilot_overview.png`` and the same data to
``results/tables/figure1_data.tsv`` (so the figure is fully reproducible).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
TABLES = REPO_ROOT / "results" / "tables"
FIG_DIR = REPO_ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main() -> int:
    df = pd.read_csv(TABLES / "results_stats.tsv", sep="\t")
    df_d09 = pd.read_csv(TABLES / "results_stats_d09.tsv", sep="\t")
    df_d11 = pd.read_csv(TABLES / "results_stats_d11.tsv", sep="\t")

    # -------- left panel: forest plot of mean delta_AUROC with 95% CI --------
    fig, axes = plt.subplots(1, 2, figsize=(14, 8.5),
                             gridspec_kw={"width_ratios": [1.8, 1.0]})

    # Collect rows with labels for the forest plot
    rows = []
    for _, r in df.iterrows():
        rows.append(("original", r.network, r.label, r.mean_delta_AUROC,
                     r.ci95_low, r.ci95_high, r.p_one_sided))
    for _, r in df_d09.iterrows():
        rows.append(("D-09 temporal", r.network, r.label, r.mean_delta_AUROC,
                     r.ci95_low, r.ci95_high, r.p_one_sided))
    for _, r in df_d11.iterrows():
        rows.append(("D-11 cap", r.network, r.label, r.mean_delta_AUROC,
                     r.ci95_low, r.ci95_high, r.p_one_sided))
    forest = pd.DataFrame(rows, columns=["regime", "network", "label",
                                         "mean", "lo", "hi", "p"])
    # order: original (4), d09 (2), d11 (8) — but actually d11 already covers
    # intogen2024/ot_all again, so the left panel would double-count.
    # Plot only the regime-specific subset.
    forest_o = forest[forest.regime == "original"].copy()
    forest_d09 = forest[forest.regime == "D-09 temporal"].copy()
    forest_d11_new = forest[(forest.regime == "D-11 cap") &
                            (~forest.label.isin(["intogen2024", "ot_all"]))].copy()
    forest_d11_dup = forest[(forest.regime == "D-11 cap") &
                            (forest.label.isin(["intogen2024", "ot_all"]))].copy()
    # combined forest plot (y axis = row index from bottom to top)
    blocks = [
        ("(a) Original 4-arm grid (uncapped)", forest_o),
        ("(b) D-09 temporal drift", forest_d09),
        ("(c) D-11 common-size cap (new arms)", forest_d11_new),
        ("(d) D-11 cap on original labels (vs uncap)", forest_d11_dup),
    ]
    ax = axes[0]
    y = 0
    yticks, ylabels = [], []
    palette = {"funmap": "#1f77b4", "string_full700": "#d62728"}
    section_y = []   # (y_center, title) for section labels
    for title, blk in blocks:
        ax.axhline(y - 0.5, color="lightgrey", linewidth=0.7, linestyle="--")
        section_y.append((y + (len(blk) - 1) / 2.0, title))
        for _, r in blk.iterrows():
            color = palette.get(r.network, "black")
            ax.errorbar(r["mean"], y, xerr=[[r["mean"] - r["lo"]],
                                            [r["hi"] - r["mean"]]],
                        fmt="o", color=color, ecolor=color, capsize=3, markersize=6)
            ax.text(r["hi"] + 0.005, y, f"{r['mean']:+.4f}",
                    va="center", fontsize=8.5, color=color)
            yticks.append(y)
            ylabels.append(f"{r.network} x {r.label}")
            y += 1
        y += 0.5

    ax.axvline(0, color="black", linewidth=1)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel(r"$\Delta$AUC ROC  (= AUROC_RWR − AUROC_degree)", fontsize=10)
    ax.set_title("Forest plot — pilot arm grid + D-09 + D-11", fontsize=11)
    ax.set_xlim(-0.02, 0.18)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle=":", linewidth=0.5)
    # section labels on the LEFT margin (well outside the plot area)
    for y_center, t in section_y:
        ax.text(-0.18, y_center, t, transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=9.5, fontweight="bold")

    # legend (placed in the upper-left of the forest panel, out of the way
    # of long y-tick labels)
    ax.scatter([], [], color=palette["funmap"], marker="o",
               label="funmap", s=42)
    ax.scatter([], [], color=palette["string_full700"], marker="o",
               label="string_full700", s=42)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)

    # -------- right panel: heatmap (uncapped vs D-11) --------
    net_order = ["funmap", "string_full700"]
    label_order = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]

    def make_mat(df_in: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        M = np.full((len(label_order), len(net_order)), np.nan)
        N = np.full((len(label_order), len(net_order)), np.nan, dtype=int)
        for i, lbl in enumerate(label_order):
            for j, net in enumerate(net_order):
                row = df_in[(df_in.label == lbl) & (df_in.network == net)]
                if len(row):
                    M[i, j] = row.mean_delta_AUROC.iloc[0]
                    N[i, j] = int(row.n_folds.iloc[0])
        return M, N

    M_unc, N_unc = make_mat(df)
    # for D-11 we have intogen_temporal_new only in D-11; others come from df_d11
    M_cap, N_cap = make_mat(df_d11)

    gs = axes[1].get_subplotspec().subgridspec(2, 1, hspace=0.45)
    ax_h1 = fig.add_subplot(gs[0])
    ax_h2 = fig.add_subplot(gs[1])

    for ax_h, M, ttl in [(ax_h1, M_unc, "(uncapped)"),
                         (ax_h2, M_cap, "(D-11 cap)")]:
        im = ax_h.imshow(M, cmap="RdYlGn", vmin=0, vmax=0.16, aspect="auto")
        ax_h.set_xticks(range(len(net_order)))
        ax_h.set_xticklabels(net_order, fontsize=8, rotation=20, ha="right")
        ax_h.set_yticks(range(len(label_order)))
        ax_h.set_yticklabels(label_order, fontsize=8)
        ax_h.set_title(f"ΔAUROC heatmap {ttl}", fontsize=10)
        for i in range(len(label_order)):
            for j in range(len(net_order)):
                v = M[i, j]
                if not np.isnan(v):
                    ax_h.text(j, i, f"{v:+.3f}", ha="center", va="center",
                              fontsize=8, color="black")
                else:
                    # explicit blank cell shading so missing data is visible
                    ax_h.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                                 facecolor="#dddddd",
                                                 edgecolor="white"))
                    ax_h.text(j, i, "n/a", ha="center", va="center",
                              fontsize=8, color="grey")
    # share colour bar with both heatmaps (use im from the last loop)
    fig.colorbar(im, ax=axes[1], fraction=0.05, pad=0.06, shrink=0.85,
                 label=r"$\Delta$AUC ROC")

    # -------- save and report --------
    out_png = FIG_DIR / "figure1_pilot_overview.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_png}")

    # companion data file
    fig_data = forest.rename(columns={"mean": "delta_AUROC_mean",
                                      "lo": "ci95_low",
                                      "hi": "ci95_high",
                                      "p": "p_one_sided"})
    fig_data.to_csv(TABLES / "figure1_data.tsv", sep="\t", index=False)
    print(f"Saved {TABLES / 'figure1_data.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
