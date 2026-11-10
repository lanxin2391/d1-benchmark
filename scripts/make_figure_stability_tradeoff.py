"""Figure: what a network prior's stability costs, and what buys it.

Three panels. A and B share the (stability, discrimination) plane, so the two
ways of moving a prior can be compared directly.

A. **Post-hoc smoothing does not buy stability.** The pre-registered
   stability-regularised scorer (shrink a gene's score toward the median of its
   score neighbourhood) and its two matched-strength controls, as
   (median per-gene rank SD, delta-AUROC) trajectories over the shrinkage
   strength. Curves are the median over the 24 arms, faint points are the arms.
   The pre-registered operator is order-preserving, so its whole trajectory
   collapses to one point; every operator that *does* move the ranking moves it
   to worse stability *and* worse discrimination.

B. **Stability is bought with seeds, and it costs no discrimination.** The same
   plane, traversed by truncating the training seed set to a fraction of the
   available positives. 24 arms in light grey, the median in blue.

C. **What does predict the instability.** Per-arm Spearman correlation between
   the per-gene rank SD and four candidate predictors, one dot per arm. The
   gene's own score coefficient of variation dominates; score-space crowding
   (the pre-registered premise) is absent or slightly negative; rank-space
   crowding is weakly positive in 24/24 arms.

Usage:
    python scripts/make_figure_stability_tradeoff.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

STYLE = Path(r"C:\Users\Lenovo\AppData\Roaming\com.ai4s.workbench\runtime"
             r"\xdg-config\opencode\skills\publication-figures\openscience.mplstyle")
if STYLE.exists():
    plt.style.use(str(STYLE))

TABLE_DIR = REPO_ROOT / "results" / "tables"
FIG_DIR = REPO_ROOT / "figures"
OUT = FIG_DIR / "figure8_stability_tradeoff.png"

INK = "#2a2723"
GREY = "#cbc6bb"
Y_LIM = (-0.20, 0.13)


def _panel_a(ax, grid):
    """Smoothing axis: the operator families in the (rank SD, delta-AUROC) plane."""
    sub = grid[(grid.alpha == 0.5) & (grid.k == 20) & (grid.variant == "global")]
    series = [
        ("srs", "#2a78d6", "SRS"),
        ("c1", "#1baf7a", "C1 random"),
        ("c2", "#eda100", "C2 global"),
    ]
    for ctl, colour, label in series:
        s = sub[sub.control == ctl]
        med = s.groupby("lam").agg(rank_sd=("median_rank_sd", "median"),
                                    delta=("delta_auroc", "median")).reset_index()
        ax.plot(med.rank_sd, med.delta, "-o", color=colour, label=label,
                markersize=7, markeredgecolor="white", markeredgewidth=0.8,
                zorder=5 if ctl == "srs" else 4)
    base_sd = sub[(sub.control == "srs") & (sub.lam == 0.0)].median_rank_sd.median()
    base_d = sub[(sub.control == "srs") & (sub.lam == 0.0)].delta_auroc.median()
    ax.scatter([base_sd], [base_d], s=170, facecolor="none", edgecolor=INK,
               linewidth=1.3, zorder=6)
    ax.annotate("unsmoothed\nprior", xy=(base_sd, base_d),
                xytext=(base_sd + 0.012, base_d + 0.022), color=INK, fontsize=8.5,
                arrowprops=dict(arrowstyle="-", color=GREY, linewidth=1))
    ax.annotate("all $\\lambda$ land here:\norder-preserving",
                xy=(base_sd, base_d), xytext=(base_sd + 0.090, base_d + 0.016),
                color="#2a78d6", fontsize=8.5,
                arrowprops=dict(arrowstyle="-", color="#2a78d6", linewidth=0.9))
    ax.set_title("A  Smoothing cannot buy stability", loc="left", fontsize=11.5)
    ax.set_xlabel("median per-gene rank SD  (left = more stable)", fontsize=9.5)
    ax.set_ylabel("$\\Delta$AUROC vs degree", fontsize=9.5)
    ax.set_ylim(*Y_LIM)
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.02), fontsize=8,
              handlelength=1.4, borderaxespad=0.3)


def _panel_b(ax, seed):
    """Seed axis: the same plane, traversed by the number of training seeds."""
    sub = seed[seed.alpha == 0.5]
    for _, d in sub.groupby(["network", "label_set"]):
        d = d.sort_values("mean_n_seeds")
        ax.plot(d.median_rank_sd, d.delta_auroc, "-", color=GREY, linewidth=0.9,
                alpha=0.9, zorder=2)
    med = (sub.groupby("frac")
           .agg(rank_sd=("median_rank_sd", "median"),
                delta=("delta_auroc", "median"),
                n=("mean_n_seeds", "median"))
           .sort_values("n"))
    ax.plot(med.rank_sd, med.delta, "-o", color="#2a78d6", linewidth=2.4,
            zorder=5, label="median of 24 arms (15 $\\rightarrow$ 304 seeds)",
            markersize=7, markeredgecolor="white", markeredgewidth=0.8)
    ax.annotate("more seeds: more stable\nand more discriminative",
                xy=(med.rank_sd.iloc[-1], med.delta.iloc[-1]),
                xytext=(med.rank_sd.iloc[-1] - 0.020, med.delta.iloc[-1] - 0.075),
                color="#2a78d6", fontsize=8.5,
                arrowprops=dict(arrowstyle="-", color="#2a78d6", linewidth=0.9))
    ax.set_title("B  Seeds buy stability, and it is free", loc="left", fontsize=11.5)
    ax.set_xlabel("median per-gene rank SD  (left = more stable)", fontsize=9.5)
    ax.set_ylabel("$\\Delta$AUROC vs degree", fontsize=9.5)
    ax.set_ylim(*Y_LIM)
    ax.legend(loc="lower left", fontsize=8.5, handlelength=1.6)


def _panel_c(ax, attrib):
    """Per-arm Spearman correlation between the rank SD and each predictor."""
    labels = [("rho_score_cv", "gene's own score CV", "#1baf7a"),
              ("rho_crowding_rank", "rank-space crowding", "#eda100"),
              ("rho_crowding", "score-space crowding", "#2a78d6"),
              ("rho_degree", "node degree", "#4a3aa7")]
    for i, (col, label, colour) in enumerate(labels):
        v = attrib[col].to_numpy()
        y = np.full(v.size, i) + np.linspace(-0.16, 0.16, v.size)
        ax.plot(v, y, "o", color=colour, alpha=0.75, markersize=6,
                markeredgecolor="white", markeredgewidth=0.6, zorder=3)
        ax.plot([np.median(v)], [i], "D", color=colour, markersize=9,
                markeredgecolor="white", markeredgewidth=0.8, zorder=4)
    ax.axvline(0.0, color=GREY, linewidth=1.0, linestyle="--", zorder=1)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels([lab for _, lab, _ in labels], fontsize=9.5)
    ax.set_xlabel("Spearman $\\rho$ with the per-gene rank SD (one dot per arm)",
                  fontsize=9.5)
    ax.set_title("C  What predicts the instability", loc="left", fontsize=11.5)
    ax.set_ylim(-0.6, len(labels) - 0.4)
    ax.grid(axis="y", visible=False)


def main() -> int:
    grid = pd.read_csv(TABLE_DIR / "stability_regularization_grid.tsv", sep="\t")
    seed = pd.read_csv(TABLE_DIR / "stability_seed_size_scaling.tsv", sep="\t")
    attrib = pd.read_csv(TABLE_DIR / "stability_crowding_attribution.tsv", sep="\t")

    fig = plt.figure(figsize=(7.6, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.78], hspace=0.58, wspace=0.30)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    _panel_a(ax_a, grid)
    _panel_b(ax_b, seed)
    _panel_c(ax_c, attrib)
    fig.suptitle("A network prior's ranking stability: what it costs and what buys it",
                 x=0.01, ha="left", fontsize=12.5, color=INK, y=1.0)
    fig.text(0.01, -0.03,
             "Medians over the 24 (network x label) arms; alpha = 0.5; 40 repeats x 5 folds. "
             "SRS = shrink toward the score-neighbourhood median; C1 = matched-strength random "
             "neighbourhood; C2 = shrink toward the global median. Grey lines in B are the 24 "
             "arms; diamonds in C are the median.",
             fontsize=8, color="#8c877d", ha="left")
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
