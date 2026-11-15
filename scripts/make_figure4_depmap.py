"""Figure 4 (deck version): the DepMap validation, three ways -- the
prior, the degree baseline and random candidates.

The point of the three-way comparison: on a cell-line functional task
the degree baseline carries almost all the signal (65-98 % of the cell
lines enriched against 3-48 % for the prior and 2-6 % for random
candidates), while on the curated driver labels the prior beats degree
in 23 of 24 arms. The prior's contribution is task-specific.
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
STYLE = pathlib.Path(r"C:\Users\PC\AppData\Roaming\com.ai4s.workbench"
                     r"\runtime\xdg-config\opencode\skills"
                     r"\publication-figures\openscience.mplstyle")
plt.style.use(str(STYLE) if STYLE.exists() else "default")

OUT = ROOT / "figures/figure4_depmap_validation.png"
T = ROOT / "results/tables"
INK, MUTED = "#2a2723", "#8c877d"
BLUE, GREY, ORANGE = "#2a78d6", "#9aa0a6", "#eb6834"


def main() -> int:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 14, "axes.titlesize": 16, "axes.labelsize": 14,
        "xtick.labelsize": 13, "ytick.labelsize": 13,
    })
    c = pd.read_csv(T / "depmap_prior_vs_degree.tsv", sep="\t")
    ks = sorted(c["K"].unique())
    labels = ["top 1 %", "top 30", "top 100"]
    prior = [c[c["K"] == k].prior_frac.median() for k in ks]
    degree = [c[c["K"] == k].degree_frac.median() for k in ks]
    random = [c[c["K"] == k].random_frac.median() for k in ks]

    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    fig.subplots_adjust(bottom=0.30, top=0.90, left=0.11, right=0.98)
    x = np.arange(len(ks))
    w = 0.26
    b1 = ax.bar(x - w, prior, w, color=BLUE, label="the prior (RWR)",
                zorder=3)
    b2 = ax.bar(x, degree, w, color=ORANGE, label="the degree baseline",
                zorder=3)
    b3 = ax.bar(x + w, random, w, color=GREY, label="random candidates",
                zorder=3)
    for bars in (b1, b2, b3):
        for r in bars:
            ax.text(r.get_x() + r.get_width() / 2, r.get_height() + 0.015,
                    f"{100 * r.get_height():.0f} %", ha="center",
                    va="bottom", fontsize=12, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("cell lines enriched (p < 0.05)")
    ax.set_ylim(0, 1.18)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.grid(axis="y", zorder=0)
    ax.legend(fontsize=13, frameon=False, loc="upper left")
    ax.set_title("On a cell-line functional task the degree baseline "
                 "carries almost all the signal", loc="left", fontsize=17)
    fig.text(0.012, 0.015,
             "DepMap: 1,361 cancer cell lines; essential = CERES < -0.5 and "
             "log2(TPM+1) > 1; hypergeometric test per cell line; the median "
             "over the 24 arms; the same folds, universe, K and cell lines "
             "for all three rankings.\nSource: "
             "results/tables/depmap_prior_vs_degree.tsv. The prior is better "
             "than random at the top 100 (47.6 % vs 6.1 %) and clearly "
             "worse than degree at every scale.",
             fontsize=11, color=MUTED, ha="left", va="bottom")
    fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
