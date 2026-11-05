"""S3 -- network-level stability and what it depends on.

S1/S2 worked at the gene level inside one arm. S3 moves up one level:
which NETWORKS are stable, and what topological properties of a network
explain that?

Inputs:
  results/runs_v40/a4_combined.tsv          6 primary networks
  results/runs_v40_sens8/a4_combined.tsv     2 STRING sensitivity arms
                                              (string_full400, string_full900)
  data/processed/networks/*.tsv              the network files
  results/tables/network_features.tsv        cached node-level features
  results/tables/label_network_overlap.tsv   coverage per (label, network)
  results/tables/stability_attribution_loo.csv  the S2 transfer test

Per network we compute
  * descriptors: nodes, edges, density, mean degree, mean clustering,
    modularity (Louvain), degree assortativity, approximate mean shortest
    path length, max k-core, fraction of nodes in the 2-core, the number
    of shared genes, and the mean coverage of the four label sets
  * stability: the SD of the 40 repeat means of delta-AUROC (averaged
    over the four labels), the spread of the four label means (how
    consistent the network's effect is across labels), and the mean
    effect

With eight networks a regression is descriptive only (work plan risk
R18): we therefore report Spearman correlations, a leave-one-network-out
predictive check, and the S2 transfer R² of each network, and we label
everything as hypothesis-generating.

Outputs:
  results/tables/network_descriptors.tsv
  results/tables/network_stability.tsv
  results/tables/network_stability_vs_topology.tsv
  figures/figure6_network_stability.png
  docs/network_stability.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

TABLES = ROOT / "results" / "tables"
FIG_DIR = ROOT / "figures"
NET_DIR = ROOT / "data" / "processed" / "networks"

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome",
            "string_full400", "string_full900"]
PRIMARY = NETWORKS[:6]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
PROVENANCE = {
    "funmap": "protein-derived", "string_phys700": "protein-derived",
    "intact": "protein-derived", "rna_coexp": "rna-derived",
    "string_full700": "curated-literature", "reactome": "curated-literature",
    "string_full400": "curated-literature", "string_full900": "curated-literature",
}
N_PATH_PAIRS = 200
SEED = 20261105


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def descriptors(net_id: str) -> dict:
    import networkx as nx

    df = pd.read_csv(NET_DIR / f"{net_id}.tsv", sep="\t")
    g = nx.Graph()
    for a, b in zip(df["gene_a"].astype(str), df["gene_b"].astype(str)):
        g.add_edge(a, b)
    n, m = g.number_of_nodes(), g.number_of_edges()
    log(f"  {net_id}: {n} nodes, {m} edges")

    t0 = time.time()
    deg = dict(g.degree())
    core = nx.core_number(g)
    clust_mean = float(np.mean(list(nx.clustering(g).values())))
    assort = float(nx.degree_assortativity_coefficient(g))
    comms = nx.community.louvain_communities(g, seed=SEED)
    modularity = float(nx.community.modularity(g, comms))
    # approximate mean shortest path length: 20 random sources, each
    # contributing its distances to every other node (20 BFS runs, not
    # one per pair)
    rng = np.random.default_rng(SEED)
    nodes = list(g.nodes())
    n = len(nodes)
    srcs = rng.choice(n, size=min(20, n - 1), replace=False)
    lengths = []
    for s in srcs:
        d = nx.single_source_shortest_path_length(g, nodes[int(s)])
        lengths.extend(v for k, v in d.items() if k != nodes[int(s)])
    path_len = float(np.mean(lengths)) if lengths else np.nan
    log(f"    descriptors in {time.time()-t0:.1f}s "
        f"(clustering {clust_mean:.4f}, modularity {modularity:.4f}, "
        f"assort {assort:+.3f}, path len {path_len:.2f})")

    return {
        "network": net_id, "n_nodes": n, "n_edges": m,
        "density": 2 * m / (n * (n - 1)),
        "mean_degree": 2 * m / n,
        "mean_clustering": clust_mean,
        "modularity": modularity,
        "degree_assortativity": assort,
        "approx_mean_path_length": path_len,
        "max_k_core": int(max(core.values())),
        "frac_in_2core": float(np.mean([c >= 2 for c in core.values()])),
        "frac_in_10core": float(np.mean([c >= 10 for c in core.values()])),
    }


def arm_stability() -> pd.DataFrame:
    """Per-arm repeat-level stability, from the two combined run files."""
    frames = [pd.read_csv(ROOT / "results" / "runs_v40" / "a4_combined.tsv",
                          sep="\t"),
              pd.read_csv(ROOT / "results" / "runs_v40_sens8" / "a4_combined.tsv",
                          sep="\t")]
    a4 = pd.concat(frames, ignore_index=True)
    a4 = a4[a4["alpha"].isna() | (a4["alpha"] == 0.5)]
    sub = a4[a4["method"].isin(["rwr", "degree"])]
    piv = (sub.pivot_table(index=["network", "label_set", "repeat", "fold"],
                           columns="method", values="auroc")
           .reset_index().dropna(subset=["rwr", "degree"]))
    piv["delta"] = piv["rwr"] - piv["degree"]
    rows = []
    for (net, lab), g in piv.groupby(["network", "label_set"]):
        rep = g.groupby("repeat")["delta"].mean()
        grp = g.groupby("repeat")["delta"]
        k = int(grp.size().iloc[0])
        grand = g["delta"].mean()
        msb = (grp.mean() - grand).pow(2).sum() * k / (len(grp) - 1)
        msw = grp.apply(lambda s: ((s - s.mean()) ** 2).sum()).sum() / (len(g) - len(grp))
        rows.append({
            "network": net, "label": lab,
            "arm_delta_mean": float(rep.mean()),
            "arm_delta_sd_repeat": float(rep.std(ddof=1)),
            "arm_delta_se_repeat": float(rep.std(ddof=1) / np.sqrt(len(rep))),
            "arm_repeat_range": float(rep.max() - rep.min()),
            "arm_icc_intra_repeat": float((msb - msw) / (msb + (k - 1) * msw)),
            "n_folds": int(len(g)),
        })
    return pd.DataFrame(rows)


def main() -> int:
    # ---- descriptors -------------------------------------------------
    desc_path = TABLES / "network_descriptors.tsv"
    if desc_path.exists():
        log(f"using cached descriptors {desc_path.name}")
        desc = pd.read_csv(desc_path, sep="\t")
    else:
        desc = pd.DataFrame([descriptors(n) for n in NETWORKS])
        desc.to_csv(desc_path, sep="\t", index=False)
        log(f"wrote {desc_path.name}")

    # ---- coverage ---------------------------------------------------
    ov = pd.read_csv(TABLES / "label_network_overlap.tsv", sep="\t")
    ov = ov[ov["network"].isin(NETWORKS)]
    cov = ov.groupby("network").agg(
        mean_label_coverage=("frac_label_in_network", "mean"),
        min_label_coverage=("frac_label_in_network", "min")).reset_index()

    # ---- per-arm and per-network stability --------------------------
    arms = arm_stability()
    net = (arms.groupby("network")
           .agg(net_delta_mean=("arm_delta_mean", "mean"),
                net_stability_sd=("arm_delta_sd_repeat", "mean"),
                net_stability_sd_max=("arm_delta_sd_repeat", "max"),
                net_cross_label_spread=("arm_delta_mean", "std"),
                net_icc=("arm_icc_intra_repeat", "mean"))
           .reset_index())
    # the two STRING sensitivity arms also belong to the curated class
    net = net.merge(desc, on="network", how="left")
    net = net.merge(cov, on="network", how="left")
    net["provenance"] = net["network"].map(PROVENANCE)
    net["is_sensitivity_arm"] = ~net["network"].isin(PRIMARY)

    # the S2 transfer result, per network
    loo = pd.read_csv(TABLES / "stability_attribution_loo.csv", sep="\t")
    loo = loo[(loo["outcome"] == "score CV") & (loo["kind"] == "network")]
    net = net.merge(loo[["held_out", "r2_out"]].rename(
        columns={"held_out": "network", "r2_out": "s2_transfer_r2_scorecv"}),
        on="network", how="left")

    net.to_csv(TABLES / "network_stability.tsv", sep="\t", index=False)
    arms.to_csv(TABLES / "network_stability_by_arm.tsv", sep="\t", index=False)
    log(f"wrote network_stability.tsv ({len(net)} networks)")

    # ---- what does the stability depend on? (descriptive, n = 8) ----
    feats = ["n_nodes", "n_edges", "density", "mean_degree",
             "mean_clustering", "modularity", "degree_assortativity",
             "approx_mean_path_length", "max_k_core", "frac_in_2core",
             "frac_in_10core", "mean_label_coverage"]
    rows = []
    for target, tlabel in (("net_stability_sd", "arm-level SD of delta-AUROC"),
                           ("net_cross_label_spread",
                            "spread of the four label means"),
                           ("net_delta_mean", "mean delta-AUROC")):
        for f in feats:
            x = net[f].to_numpy(dtype=float)
            y = net[target].to_numpy(dtype=float)
            ok = np.isfinite(x) & np.isfinite(y)
            rho = float(sps.spearmanr(x[ok], y[ok]).statistic)
            p = float(sps.spearmanr(x[ok], y[ok]).pvalue)
            rows.append({"target": target, "target_label": tlabel,
                         "feature": f, "spearman": rho, "p": p,
                         "n_networks": int(ok.sum())})
    res = pd.DataFrame(rows)
    res["q_bh"] = np.nan
    for t in res["target"].unique():
        m = res["target"] == t
        p = res.loc[m, "p"].to_numpy()
        order = np.argsort(p)
        q = p[order] * len(p) / (np.arange(len(p)) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        out = np.empty_like(q)
        out[order] = np.clip(q, 0, 1)
        res.loc[m, "q_bh"] = out

    # leave-one-network-out predictive check (descriptive)
    loo_pred = []
    for target in ["net_stability_sd", "net_cross_label_spread"]:
        for hold in net["network"]:
            tr = net[net["network"] != hold]
            te = net[net["network"] == hold]
            best = None
            for f in feats:
                x = tr[f].to_numpy(dtype=float)
                y = tr[target].to_numpy(dtype=float)
                ok = np.isfinite(x) & np.isfinite(y)
                if ok.sum() < 5:
                    continue
                b = np.polyfit(x[ok], y[ok], 1)
                pred = np.polyval(b, te[f].to_numpy(dtype=float))
                err = abs(pred - te[target].to_numpy(dtype=float))[0]
                if best is None or err < best[0]:
                    best = (err, f, float(pred[0]))
            loo_pred.append({"target": target, "held_out": hold,
                             "best_feature": best[1] if best else None,
                             "abs_error": best[0] if best else np.nan,
                             "predicted": best[2] if best else np.nan,
                             "observed": float(te[target].iloc[0])})
    res.to_csv(TABLES / "network_stability_vs_topology.tsv", sep="\t",
               index=False)
    pd.DataFrame(loo_pred).to_csv(
        TABLES / "network_stability_loo.tsv", sep="\t", index=False)

    make_figure(net, res)
    print(net.round(4).to_string(index=False))
    return 0


def make_figure(net: pd.DataFrame, res: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    prov_color = {"protein-derived": "#1f77b4", "rna-derived": "#2ca02c",
                  "curated-literature": "#d62728"}
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    ax = axes[0]
    x = net["mean_delta_mean" if False else "net_delta_mean"]
    y = net["net_stability_sd"]
    for i, (_, r) in enumerate(net.iterrows()):
        ax.scatter(r["net_delta_mean"], r["net_stability_sd"],
                   s=70, color=prov_color[r["provenance"]],
                   marker="o" if not r["is_sensitivity_arm"] else "^",
                   edgecolor="black", linewidth=0.6, zorder=3)
        ax.annotate(r["network"], (r["net_delta_mean"], r["net_stability_sd"]),
                    fontsize=7.5, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("mean delta-AUROC over the four labels")
    ax.set_ylabel("arm-level SD of delta-AUROC (stability)")
    ax.set_title("A  Which networks are stable?", fontsize=11)
    handles = [plt.Line2D([], [], color=prov_color[k], lw=6, marker="o",
                          markersize=6) for k in prov_color]
    handles.append(plt.Line2D([], [], color="black", lw=0, marker="^",
                              markersize=8))
    ax.legend(handles, list(prov_color) + ["sensitivity arm"], fontsize=8)

    ax = axes[1]
    sub = res[res["target"] == "net_stability_sd"].sort_values("spearman")
    ypos = np.arange(len(sub))
    ax.barh(ypos, sub["spearman"], color="#4C72B0")
    ax.set_yticks(ypos)
    ax.set_yticklabels(sub["feature"], fontsize=8.5)
    ax.axvline(0, color="black", lw=1)
    ax.set_xlabel("Spearman with the arm-level stability (n = 8 networks)")
    ax.set_title("B  What does stability depend on?", fontsize=11)
    for i, (r, v) in enumerate(zip(sub.itertuples(), sub["spearman"])):
        ax.text(v + (0.02 if v >= 0 else -0.02), i, f"{v:+.2f}",
                va="center", ha="left" if v >= 0 else "right", fontsize=7.5)

    ax = axes[2]
    for i, (_, r) in enumerate(net.iterrows()):
        ax.scatter(r["mean_label_coverage"], r["net_stability_sd"], s=60,
                   color=prov_color[r["provenance"]], zorder=3)
        ax.annotate(r["network"], (r["mean_label_coverage"],
                                   r["net_stability_sd"]),
                    fontsize=7, xytext=(4, 3), textcoords="offset points")
    ax.set_xlabel("mean fraction of the label set inside the network")
    ax.set_ylabel("arm-level SD of delta-AUROC")
    ax.set_title("C  Coverage does not explain stability", fontsize=11)

    fig.suptitle("Figure 6  Network-level stability of the network prior "
                 "(8 networks, 4 labels, 40 repeats)", fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG_DIR / "figure6_network_stability.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    log(f"wrote figure6_network_stability.png")


if __name__ == "__main__":
    raise SystemExit(main())
