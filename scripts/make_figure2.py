"""Figure 2: Heatmap of delta-AUROC across networks × labels (the paper's main figure).

Inputs:
    results/tables/a4_table1_delta_auroc.tsv  (Table 1 from A4)
    results/tables/d11_table.tsv             (D-11 cap, alpha=0.5 only)
    results/tables/d12_table.tsv             (D-12 shared, alpha=0.5 only)
    results/tables/d09_table.tsv             (D-09 temporal drift, alpha=0.5)

Output:
    figures/figure2_main_heatmap.png        (primary figure: native, alpha=0.5)
    figures/figure2_d09_d11_d12.png        (secondary: D-09/D-11/D-12 vs native)

The primary heatmap visualises Table 1 at alpha=0.5 across the 6 networks
(rows) and 4 labels (columns). A second figure places the D-09/D-11/D-12
results side by side with the native (uncap) values so the reader can see
how robust each arm is to the controls.

Both figures match the style of B's figure1_pilot_overview.png (which it
eventually replaces in the paper).
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

NET_ORDER = ["funmap", "rna_coexp", "string_phys700",
             "intact", "string_full700", "reactome"]
LABEL_ORDER = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]


def _to_matrix(df: pd.DataFrame, value: str) -> pd.DataFrame:
    """Pivot a long df (rows = network, label_set, ...) to a network x label matrix."""
    if value in df.columns:
        m = df.pivot_table(index="network", columns="label_set", values=value)
    else:
        m = df.pivot_table(index="network", columns="label_set", values=value)
    return m.reindex(index=NET_ORDER, columns=LABEL_ORDER)


def _make_heatmap(ax, mat, title, *, cmap="RdYlGn", vmin=None, vmax=None):
    """Render one annotated heatmap on the given axes."""
    arr = mat.values.astype(float)
    im = ax.imshow(arr, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(len(LABEL_ORDER)))
    ax.set_xticklabels(LABEL_ORDER, rotation=30, ha="right", fontsize=9)
    ax.set_yticks(np.arange(len(NET_ORDER)))
    ax.set_yticklabels(NET_ORDER, fontsize=9)
    ax.set_title(title, fontsize=11)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = arr[i, j]
            if np.isnan(v):
                txt, color = "n/a", "0.6"
            else:
                txt = f"{v:+.3f}"
                # White text for cells with extreme values when vmin/vmax are set
                if vmin is not None and vmax is not None:
                    mid = (vmin + vmax) / 2
                    span = vmax - vmin
                    color = "white" if abs(v - mid) > span * 0.30 else "black"
                else:
                    color = "black"
            ax.text(j, i, txt, ha="center", va="center",
                     fontsize=8, color=color)
    return im


def main():
    # ----- Figure 2: primary (Table 1 at alpha=0.5, native, uncap) -----
    a4 = pd.read_csv(TABLES / "a4_table1_delta_auroc.tsv", sep="\t")
    mat_main = _to_matrix(a4, "a0.5_mean")

    fig, ax = plt.subplots(figsize=(7, 5))
    im = _make_heatmap(ax, mat_main,
                        title="Figure 2: delta-AUROC (alpha=0.5, native, uncap)",
                        cmap="RdYlGn", vmin=0, vmax=0.18)
    fig.colorbar(im, ax=ax, label="delta-AUROC (RWR - degree)")
    fig.tight_layout()
    out1 = FIGURES / "figure2_main_heatmap.png"
    fig.savefig(out1, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out1}")

    # ----- Figure 2 (companion): D-09 / D-11 / D-12 vs native -----
    # Build a 4-row stack of heatmaps.
    a4_native = _to_matrix(a4, "a0.5_mean")

    d09 = pd.read_csv(TABLES / "d09_table.tsv", sep="\t")
    # d09 only has a single "label" column (the test set is always
    # intogen_temporal_new; training is intogen2024 minus temporal_new).
    # Treat it as a one-column result for the strip plot.
    d09_05 = d09[d09.alpha == 0.5].set_index("network")

    d11 = pd.read_csv(TABLES / "d11_table.tsv", sep="\t")
    # d11 has columns: label= "<label>_cap"
    # Strip the "_cap" suffix to align with the canonical label names.
    d11["label_canon"] = d11["label"].str.replace("_cap$", "", regex=True)
    mat_d11 = d11.pivot_table(index="network", columns="label_canon",
                                values="delta_auroc").reindex(
        index=NET_ORDER, columns=LABEL_ORDER)

    d12 = pd.read_csv(TABLES / "d12_table.tsv", sep="\t")
    d12["label_canon"] = d12["label"].str.replace("_shared$", "", regex=True)
    mat_d12 = d12.pivot_table(index="network", columns="label_canon",
                                values="delta_auroc").reindex(
        index=NET_ORDER, columns=LABEL_ORDER)

    # D-09 only covers intogen_temporal_new (one column). Show as a vector
    # rather than a full heatmap to be honest about scope.
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    _make_heatmap(axes[0], a4_native,
                   title="Native (uncap) — primary",
                   cmap="RdYlGn", vmin=0, vmax=0.18)
    _make_heatmap(axes[1], mat_d11,
                   title="D-11 per-network cap",
                   cmap="RdYlGn", vmin=0, vmax=0.18)
    _make_heatmap(axes[2], mat_d12,
                   title="D-12 shared universe",
                   cmap="RdYlGn", vmin=0, vmax=0.18)
    fig.suptitle("Figure 2: delta-AUROC across controls (alpha=0.5)",
                  fontsize=13)
    fig.tight_layout()
    out2 = FIGURES / "figure2_d09_d11_d12.png"
    fig.savefig(out2, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out2}")

    # ----- Figure 2 (extra): D-09 temporal drift as 1-row strip -----
    fig, ax = plt.subplots(figsize=(8, 2.5))
    d09_strip = d09_05[["delta_auroc_mean"]].reindex(NET_ORDER)
    d09_strip.columns = ["delta_auroc_mean"]
    im = ax.imshow(d09_strip.values, cmap="RdYlGn", vmin=-0.02, vmax=0.10,
                    aspect="auto")
    ax.set_xticks([0])
    ax.set_xticklabels(["intogen_temporal_new (test)"], fontsize=10)
    ax.set_yticks(np.arange(len(NET_ORDER)))
    ax.set_yticklabels(NET_ORDER, fontsize=9)
    ax.set_title("D-09 temporal drift: train on intogen2024 \\ temporal_new, "
                  "test on temporal_new (alpha=0.5)", fontsize=11)
    arr = d09_strip.values.astype(float)
    for i in range(arr.shape[0]):
        v = arr[i, 0]
        ax.text(0, i, f"{v:+.3f}" if not np.isnan(v) else "n/a",
                 ha="center", va="center", fontsize=9)
    fig.colorbar(im, ax=ax, label="delta-AUROC")
    fig.tight_layout()
    out3 = FIGURES / "figure2_d09_temporal_strip.png"
    fig.savefig(out3, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out3}")


if __name__ == "__main__":
    main()