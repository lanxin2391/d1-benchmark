"""Figure 5 -- the split-design comparison (Control 23), final layout.

Panel A: for every arm, the range its margin takes over the six designs,
  with the pre-registered design's value marked. Twenty-four tight
  segments is the whole message: the margin does not move.
Panel B: the one quantity that does move -- the median width of the
  conservative (gene-cluster) interval -- expressed relative to the
  pre-registered design, so no truncated axis is needed; the absolute
  widths are printed on the bars.
The three design-invariant verdicts sit in a caption strip under the
figure, where nothing can collide with them.
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
STYLE = pathlib.Path(
    r"C:\Users\PC\AppData\Roaming\com.ai4s.workbench\runtime\xdg-config"
    r"\opencode\skills\publication-figures\openscience.mplstyle")
if STYLE.exists():
    plt.style.use(str(STYLE))

BLUE = "#2a78d6"
AQUA = "#1baf7a"
VIOLET = "#4a3aa7"
GREY = "#8c877d"
INK = "#2a2723"

DESIGNS = ["k5_r10", "k5_r20", "k5_r40", "k10_r10", "k10_r20", "k20_r10"]
LABELS = {"k5_r10": "5 x 10", "k5_r20": "5 x 20", "k5_r40": "5 x 40",
          "k10_r10": "10 x 10", "k10_r20": "10 x 20", "k20_r10": "20 x 10"}
SHORT = {"funmap": "funmap", "rna_coexp": "rna", "string_phys700": "str-phys",
         "intact": "intact", "string_full700": "str-full", "reactome": "reactome"}
SHORT_L = {"intogen2024": "IntOGen", "ot_all": "OT-all",
           "ot_nolit": "OT-nolit", "intogen_temporal_new": "temporal"}


def main() -> int:
    T = ROOT / "results/tables"
    arms = pd.read_csv(T / "c23_split_design_arms.tsv", sep="\t")
    summary = pd.read_csv(T / "c23_split_design_summary.tsv",
                          sep="\t").set_index("design").loc[DESIGNS]

    piv = arms.pivot_table(index=["network", "label"], columns="design",
                            values="margin_mean")
    piv = piv.sort_values("k5_r40", ascending=True)
    y = np.arange(len(piv))

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(10.6, 6.0),
        gridspec_kw={"width_ratios": [2.45, 1]})

    # ---------------- panel A ----------------------------------------
    lo = piv.min(axis=1).to_numpy()
    hi = piv.max(axis=1).to_numpy()
    pre = piv["k5_r40"].to_numpy()
    axA.hlines(y, lo, hi, color=BLUE, lw=1.5, alpha=0.5, zorder=2)
    axA.plot(piv["k10_r10"].to_numpy(), y, "o", ms=4.4, color=AQUA,
             mec="white", mew=0.5, zorder=3, label="the other five designs")
    axA.plot(pre, y, "o", ms=5.2, color=VIOLET, mec="white", mew=0.6,
             zorder=4, label="the pre-registered design (5 x 40)")
    axA.axvline(0, color=GREY, lw=1.0, zorder=1)
    axA.set_yticks(y)
    axA.set_yticklabels([f"{SHORT[n]}  x  {SHORT_L[l]}" for n, l in piv.index],
                        fontsize=7.6)
    axA.set_ylim(len(piv) - 0.4, -1.0)
    axA.set_xlim(-0.030, 0.195)
    axA.set_xlabel("margin over the degree baseline (delta-AUROC)",
                    fontsize=9.5, labelpad=6)
    axA.set_title("A  Every arm, every design: the margin does not move",
                  loc="left", fontsize=11, pad=8)
    axA.grid(axis="y", lw=0.5, color="#e8e5df")
    axA.legend(loc="lower left", bbox_to_anchor=(0.015, 0.02), frameon=False,
               fontsize=8.0, handletextpad=0.4, labelspacing=0.5)

    neg_i = int(np.where(pre < 0)[0][0])
    axA.annotate("the only negative arm\n(Reactome x the newest drivers)",
                 xy=(float(pre[neg_i]), neg_i), xytext=(0.062, 3.2),
                 fontsize=8.0, color=INK, va="center",
                 arrowprops=dict(arrowstyle="->", color=GREY, lw=0.9,
                                 shrinkA=0, shrinkB=4))

    # ---------------- panel B ----------------------------------------
    widths = summary["gene_ci_width_median"].to_numpy()
    base = widths[DESIGNS.index("k5_r40")]
    rel = 100.0 * (widths - base) / base
    x = np.arange(len(DESIGNS))
    cols = [VIOLET if d == "k5_r40" else BLUE for d in DESIGNS]
    axB.bar(x, rel, color=cols, width=0.6, zorder=3)
    axB.axhline(0, color=GREY, lw=1.0, zorder=4)
    for xi, (r, w) in enumerate(zip(rel, widths)):
        axB.text(xi, r + (0.16 if r >= 0 else -0.30), f"{w:.4f}",
                 ha="center", fontsize=7.8, color=INK)
    axB.set_xticks(x)
    axB.set_xticklabels([LABELS[d] for d in DESIGNS], fontsize=8.0,
                         rotation=35, ha="right")
    axB.set_ylim(-1.15, 3.4)
    axB.set_ylabel("interval width relative to the\npre-registered design (%)",
                    fontsize=9.0)
    axB.set_title("B  What the design does move", loc="left", fontsize=11,
                  pad=8)
    axB.grid(axis="y", lw=0.5, color="#e8e5df")
    axB.annotate("the pre-registered design", xy=(2, 0.0), xytext=(0.02, 2.3),
                 fontsize=8.0, color=VIOLET,
                 arrowprops=dict(arrowstyle="->", color=VIOLET, lw=0.9,
                                 shrinkA=0, shrinkB=3))

    fig.suptitle("The conclusion does not depend on the cross-validation "
                 "design", fontsize=12.5, x=0.012, ha="left", y=0.985)
    fig.text(0.012, 0.940,
             "6 networks x 4 catalogues, alpha = 0.5; between designs the "
             "median margin moves by less than 5 % and no arm moves by more "
             "than 0.0075 AUROC",
             fontsize=8.6, color=INK, ha="left")
    fig.text(0.012, 0.022,
             "In every design: 23 of 24 arms positive, 24 of 24 certified by "
             "the repeat-level interval, 20 of 24 by the gene-cluster "
             "bootstrap; the design moves the width of the conservative "
             "interval, by at most 2.3 % (the ten-fold designs are the "
             "widest: +1.6 % and +2.3 %).",
             fontsize=8.2, color=INK, ha="left")

    fig.tight_layout(rect=[0, 0.055, 1, 0.925])
    out = ROOT / "figures/figure5_split_design.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    print("wrote", out, out.stat().st_size // 1024, "KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
