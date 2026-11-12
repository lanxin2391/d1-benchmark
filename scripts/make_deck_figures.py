"""The three deck-specific figures for the 1118 defence.

They are drawn at the size they will occupy on the slide (16:9, 13.33 x
7.5 in), in the house palette, with Times New Roman and with the
in-figure text at 13-19 pt, so nothing has to be re-sized afterwards.

  deck_research_map.png   slide 2  -- the logic chain of the project
  deck_main_heatmap.png   slide 11 -- the 30-cell primary result
  deck_stability.png      slide 12 -- where the stability comes from
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
STYLE = Path(r"C:\Users\PC\AppData\Roaming\com.ai4s.workbench\runtime"
              r"\xdg-config\opencode\skills\publication-figures"
              r"\openscience.mplstyle")
plt.style.use(str(STYLE) if STYLE.exists() else "default")
T = ROOT / "results/tables"
FIG = ROOT / "figures"

BLUE = "#2a78d6"
AQUA = "#1baf7a"
YELLOW = "#eda100"
GREEN = "#008300"
VIOLET = "#4a3aa7"
RED = "#e34948"
ORANGE = "#eb6834"
INK, MUTED = "#2a2723", "#8c877d"
DIVERGING = LinearSegmentedColormap.from_list(
    "os_div", ["#104281", "#3987e5", "#cde2fb", "#e7e3da",
               "#f6c3c2", "#e34948", "#8c1f1e"])


def style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 14,
        "axes.titlesize": 16,
        "axes.labelsize": 14,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
    })


# ------------------------------------------------------------------ 1
def research_map():
    """The logic chain, as five horizontal bands (a column of boxes would
    not leave room for the text at a readable size)."""
    style()
    fig = plt.figure(figsize=(11.6, 4.9))
    steps = [
        ("1", "What the field reports", BLUE,
         "RWR beats a degree-only ranking, and the network choice matters "
         "(FunMap > STRING > Reactome, a 0.25-AUROC gap) -- the basis of the "
         "'best network' question."),
        ("2", "What was never measured", VIOLET,
         "How much of that gain is topology, how much is node degree, and "
         "how much is study visibility (how often a gene has been published "
         "about)? No published study separates the three."),
        ("3", "Our design (public data only)", AQUA,
         "6 networks x 5 label sets; RWR frozen at alpha = 0.5; 5 folds x 40 "
         "repeats; 11 leakage controls; pre-registered and frozen "
         "(Zenodo 10.5281/zenodo.23121762)."),
        ("4", "What we found", ORANGE,
         "The margin over degree is +0.081 (median) and positive in 23/24 "
         "arms -- but 10-22 % of our 'negatives' are curated drivers, and "
         "the prior transfers to unseen drivers in 24/24 arms."),
        ("5", "What it means", GREEN,
         "A topology prior adds a real but modest and label-dependent "
         "signal; the headline margins reported for this family of methods "
         "are largely a publication-volume artefact."),
    ]
    top, height, gap = 0.855, 0.128, 0.022
    for i, (num, title, colour, body) in enumerate(steps):
        y = top - i * (height + gap) - height
        ax = fig.add_axes([0.012, y, 0.976, height])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.add_patch(FancyBboxPatch(
            (0.0, 0.0), 1.0, 1.0, boxstyle="round,pad=0.004,rounding_size=0.12",
            linewidth=1.4, edgecolor=colour, facecolor="white"))
        ax.add_patch(FancyBboxPatch(
            (0.006, 0.06), 0.075, 0.88,
            boxstyle="round,pad=0.004,rounding_size=0.12",
            linewidth=0, facecolor=colour))
        ax.text(0.043, 0.50, num, ha="center", va="center", fontsize=19,
                color="white", fontweight="bold")
        ax.text(0.10, 0.66, title, ha="left", va="center", fontsize=15.5,
                color=colour, fontweight="bold")
        ax.text(0.10, 0.28, body, ha="left", va="center", fontsize=13,
                color=INK)
    fig.text(0.012, 0.955,
             "From published claims to a measurable question: how much of a "
             "network prior's advantage is topology?",
             fontsize=18, color=INK, ha="left", va="top")
    fig.text(0.012, 0.055,
             "Three falsifiable hypotheses: H-vis (much of the margin is "
             "study visibility) | H-stable (the prior's stability is a "
             "network property that no learner can remove) | H-lag (curated "
             "networks fail on recently discovered drivers).\nEverything "
             "below is a re-analysis of public data; the three source "
             "papers are Bailey 2018 (Cell), Shi 2025 (Nat Cancer) and "
             "Repana 2019 (Genome Biol).",
             fontsize=12, color=MUTED, ha="left", va="bottom")
    fig.savefig(FIG / "deck_research_map.png", dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote deck_research_map.png")


# ------------------------------------------------------------------ 2
def main_heatmap():
    style()
    t = pd.read_csv(T / "a4_table1_delta_auroc_v2.tsv", sep="\t")
    nets = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
    labs = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new",
            "ncg_known"]
    piv = t.pivot_table(index="network", columns="label_set",
                        values="a0.5_mean").reindex(index=nets, columns=labs)
    vals = piv.values
    fig = plt.figure(figsize=(9.8, 4.6))
    ax = fig.add_axes([0.135, 0.30, 0.72, 0.55])
    norm = TwoSlopeNorm(vmin=min(-0.02, vals.min()), vcenter=0.0,
                        vmax=0.18)
    im = ax.imshow(vals, cmap=DIVERGING, norm=norm, aspect="auto")
    ax.set_xticks(range(len(labs)))
    ax.set_xticklabels(["IntOGen\n2024", "OT all", "OT\nno-lit",
                         "IntOGen\ntemporal", "NCG\n2019"], fontsize=13)
    ax.set_yticks(range(len(nets)))
    ax.set_yticklabels(nets, fontsize=13)
    ax.tick_params(length=0)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            v = vals[i, j]
            dark = abs(v) > 0.09
            ax.text(j, i, f"{v:+.3f}", ha="center", va="center", fontsize=12.5,
                    color="white" if dark else INK)
    # mark the single negative cell
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            if vals[i, j] < 0:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           fill=False, edgecolor=INK,
                                           linewidth=2.2))
    cax = fig.add_axes([0.885, 0.30, 0.018, 0.55])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("delta-AUROC (RWR - degree)", fontsize=13)
    cb.ax.tick_params(labelsize=12)
    cb.outline.set_visible(False)
    fig.text(0.012, 0.015,
             "alpha = 0.5, native universe, 5 folds x 40 repeats; the fifth "
             "column (NCG 2019) is a pre-planned exploratory extension. The "
             "outlined cell is the only negative arm (Reactome x the newest "
             "IntOGen drivers).\nThe pre-registered H1 rule (95 % CI "
             "excludes 0 in every arm) is met under four of five interval "
             "constructions (24/24) and not under the conservative "
             "gene-cluster interval (20/24) -- the four exceptions are all "
             "the newest-driver label.\nSource: our own run, "
             "results/tables/a4_table1_delta_auroc_v2.tsv.",
             fontsize=11.5, color=MUTED, ha="left", va="top")
    fig.savefig(FIG / "deck_main_heatmap.png", dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote deck_main_heatmap.png")


# ------------------------------------------------------------------ 3
def stability():
    style()
    ns = pd.read_csv(T / "network_stability.tsv", sep="\t")
    ns = ns.sort_values("net_stability_sd")
    prov_colour = {"protein-derived": BLUE, "rna-derived": AQUA,
                   "curated-literature": ORANGE}
    fig = plt.figure(figsize=(9.8, 4.6))
    ax = fig.add_axes([0.075, 0.28, 0.40, 0.58])
    cols = [prov_colour[p] for p in ns.provenance]
    ax.barh(range(len(ns)), ns.net_stability_sd, color=cols, height=0.62,
            zorder=3)
    ax.set_yticks(range(len(ns)))
    ax.set_yticklabels(ns.network, fontsize=12)
    ax.set_xlabel("arm-level SD of delta-AUROC (stability)", fontsize=13)
    ax.set_title("A  Which networks are stable?", loc="left", fontsize=16)
    for i, (v, cov) in enumerate(zip(ns.net_stability_sd,
                                     ns.mean_label_coverage)):
        ax.text(v + 0.00025, i, f"{v:.4f}  (cov {cov:.2f})",
                va="center", fontsize=11, color=INK)
    ax.set_xlim(0, 0.0115)
    ax.grid(axis="x", zorder=0)
    handles = [plt.Line2D([], [], marker="s", ls="", ms=8, color=c, label=k)
               for k, c in prov_colour.items()]
    ax.legend(handles=handles, fontsize=11, loc="lower right",
              bbox_to_anchor=(0.995, 0.015), ncol=1, frameon=False,
              borderaxespad=0.0, labelspacing=0.35,
              handletextpad=0.5)

    ax2 = fig.add_axes([0.655, 0.28, 0.315, 0.58])
    att = pd.read_csv(T / "stability_attribution.tsv", sep="\t")
    s = att[att.outcome == "score CV"].copy()
    s["abs"] = s.beta.abs()
    s = s.sort_values("abs", ascending=False).head(6).iloc[::-1]
    y = np.arange(len(s))
    ax2.barh(y, s.beta, color=[BLUE if b > 0 else RED for b in s.beta],
             height=0.62, zorder=3)
    ax2.errorbar(s.beta, y, xerr=[s.beta - s.ci_low, s.ci_high - s.beta],
                 fmt="none", ecolor=INK, elinewidth=1.6, capsize=3, zorder=4)
    ax2.axvline(0, color=MUTED, lw=1.2, zorder=2)
    ax2.set_yticks(y)
    short = {"log_score_cv": "score CV",
             "log_score_mean": "score mean",
             "log_crowding_scott": "crowding",
             "log_local_score_density": "local density",
             "log_degree": "degree",
             "k_core": "k-core"}
    ax2.set_yticklabels([short.get(f, f.replace("log_", "log "))
                         for f in s.feature], fontsize=12)
    ax2.set_xlabel("beta in the gene-level model of the score CV",
                   fontsize=13)
    ax2.set_title("B  What drives gene-level instability?", loc="left",
                  fontsize=16)
    for i, (b, q, r) in enumerate(zip(s.beta, s.q_bh,
                                    s.marginal_spearman_median)):
        ax2.text(0.99, i, f"{b:+.2f}   rho {r:+.2f}", va="center",
                 ha="right", fontsize=11.5, color=INK,
                 transform=ax2.get_yaxis_transform())
    ax2.set_xlim(-0.30, 1.20)
    ax2.grid(axis="x", zorder=0)

    fig.text(0.012, 0.015,
             "A: the arm-level SD over 200 evaluations per arm, 8 networks "
             "(the two STRING sensitivity arms included). B: the "
             "gene-level model of the score CV over 24 arms x 40 repeats; "
             "the gene's own score CV explains it (rho = 1.00, "
             "q = 2e-137),\nno classical topology descriptor does. "
             "Source: results/tables/network_stability.tsv, "
             "stability_attribution.tsv.",
             fontsize=12, color=MUTED, ha="left", va="top")
    fig.savefig(FIG / "deck_stability.png", dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote deck_stability.png")


if __name__ == "__main__":
    research_map()
    main_heatmap()
    stability()
