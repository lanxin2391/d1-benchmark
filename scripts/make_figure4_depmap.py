"""Figure 4 -- the DepMap external validation, sized for the course deck.

The canvas is the final on-slide size, so the in-figure text stays
16-18 pt when the figure is placed on a 13.33 x 7.5 in slide.

Message: for each (network, label) arm, the fraction of the 1,361 DepMap
cell lines in which the prior's top-K candidates are enriched for genes
essential in that cell line (CERES < -0.5 and log2(TPM+1) > 1,
hypergeometric test per cell line). At K = 30 -- about the number of
genuinely essential genes in one cell line -- there is no enrichment
beyond chance; at K = 100 and at 1 % of the network there is. The prior
is informative about coarse, set-level composition, not about individual
gene-level calls.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")   # headless: the Qt backend resizes the figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
STYLE = Path(r"C:\Users\PC\AppData\Roaming\com.ai4s.workbench\runtime"
              r"\xdg-config\opencode\skills\publication-figures"
              r"\openscience.mplstyle")
plt.style.use(str(STYLE) if STYLE.exists() else "default")

OUT = ROOT / "figures" / "figure4_depmap_validation.png"
TABLE = ROOT / "results" / "tables" / "depmap_validation.tsv"

BLUE_SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf",
            "#184f95", "#104281"]
CMAP = LinearSegmentedColormap.from_list("os_blue", BLUE_SEQ)
INK, MUTED, CHANCE = "#2a2723", "#8c877d", "#e34948"

NETS = ["funmap", "rna_coexp", "string_phys700", "intact",
        "string_full700", "reactome"]
# NOTE: in a serif face the "rn" of "rna" reads as an "m", so the
# network ids keep their underscores ("rna_coexp") on the slides.
NETS_SHORT = ["funmap", "rna_coexp", "string_phys", "intact",
               "string_full", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
KS = [30.0, 100.0, 0.01]
KS_TICK = ["K = 30", "K = 100", "K = 1 %"]
XTICK = ["IntOGen", "OT all", "OT nolit", "temporal"]


def main() -> int:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 14,
        "axes.titlesize": 16,
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    })
    d = pd.read_csv(TABLE, sep="\t")
    d = d[d.threshold == -0.5]

    fig = plt.figure(figsize=(9.6, 5.2))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.0, 1, 1, 1], wspace=0.95,
                          left=0.105, right=0.925, top=0.795, bottom=0.215)

    # ---- panel A: the scale dependence ---------------------------
    ax = fig.add_subplot(gs[0, 0])
    meds, los, his = [], [], []
    for k in KS:
        v = d[d.K == k].frac_celllines_p_lt_005
        meds.append(v.median())
        los.append(v.quantile(0.25))
        his.append(v.quantile(0.75))
    x = np.arange(len(KS))
    ax.bar(x, meds, width=0.58, color=BLUE_SEQ[3], zorder=3)
    for i, (lo, hi) in enumerate(zip(los, his)):
        ax.plot([i, i], [lo, hi], color=INK, lw=2.0, zorder=4,
                solid_capstyle="butt")
        ax.text(i, hi + 0.035, f"{100 * meds[i]:.0f} %", ha="center",
                va="bottom", fontsize=14, color=INK)
    ax.axhline(0.05, color=CHANCE, lw=1.8, ls="--", zorder=5)
    ax.text(-0.46, 0.70, "- - chance 5 %", color=CHANCE, fontsize=12,
            ha="left", va="center")
    ax.set_xticks(x)
    ax.set_xticklabels(["30", "100", "1 %"], fontsize=14)
    ax.set_xlabel("candidates per arm (K)", fontsize=14)
    ax.set_ylabel("cell lines enriched (p < 0.05)", fontsize=14)
    ax.set_ylim(0, 0.78)
    ax.set_yticks([0, 0.2, 0.4, 0.6])
    ax.set_title("A  A coarse candidate set\n    is needed", loc="left",
                 fontsize=16)
    ax.grid(axis="y", zorder=0)

    # ---- panels B-D: the heatmaps --------------------------------
    hm = None
    for j, k in enumerate(KS):
        axh = fig.add_subplot(gs[0, j + 1])
        piv = (d[d.K == k]
               .pivot_table(index="network", columns="label",
                            values="frac_celllines_p_lt_005")
               .reindex(index=NETS, columns=LABELS))
        hm = axh.imshow(piv.values, cmap=CMAP, vmin=0, vmax=1, aspect="auto")
        axh.set_xticks(range(len(LABELS)))
        axh.set_xticklabels(XTICK, fontsize=11, rotation=30,
                           ha="right", rotation_mode="anchor")
        axh.set_yticks(range(len(NETS)))
        axh.set_yticklabels(NETS_SHORT if j == 0 else [""] * len(NETS),
                             fontsize=12)
        axh.tick_params(length=0)
        axh.grid(False)
        for s in axh.spines.values():
            s.set_visible(False)
        axh.set_title(f"{'BCD'[j]}  {KS_TICK[j]}", loc="left", fontsize=16)
        vals = piv.values
        for idx in (np.unravel_index(np.nanargmax(vals), vals.shape),
                    np.unravel_index(np.nanargmin(vals), vals.shape)):
            i, k2 = idx
            val = vals[i, k2]
            axh.text(k2, i, f"{100 * val:.0f} %", ha="center", va="center",
                     fontsize=12,
                     color="#ffffff" if val > 0.45 else INK,
                     fontweight="bold")

    cax = fig.add_axes([0.945, 0.215, 0.013, 0.58])
    cb = fig.colorbar(hm, cax=cax)
    cb.set_label("cell lines enriched", fontsize=13)
    cb.ax.tick_params(labelsize=11)
    cb.outline.set_visible(False)

    fig.text(0.012, 0.965,
             "DepMap external validation: set-level signal, "
             "not gene-level calls",
             fontsize=19, color=INK, ha="left", va="top")
    fig.text(0.012, 0.035,
             "Essential = CERES < -0.5 and log2(TPM+1) > 1; hypergeometric "
             "test per cell line; 1,361 DepMap cell lines; median over the "
             "200 (repeat, fold)\nevaluations per arm (panel A: median and "
             "interquartile range over the 24 arms). Source: DepMap "
             "CRISPRGeneEffect, public release.",
             fontsize=11.5, color=MUTED, ha="left", va="bottom")
    fig.savefig(OUT, dpi=200, facecolor="white")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
