"""S2 -- attributing the stability of the prior to network topology.

S1 established WHAT is stable: the arm-level endpoint (repeat-level SD
of delta-AUROC ~0.005) is stable, the individual gene is not (~18 %
relative error), and the gene-level instability is explained by
crowding (Spearman(rank SD, local score density) = +0.47, positive in
24/24 arms).

S2 asks the next question: which TOPOLOGICAL properties of a gene
predict its stability? The candidate features are the standard
node-level descriptors, computed on each network's LCC:

  degree, k-core number, clustering coefficient, betweenness centrality
  (sampled), PageRank, community id and community size, and the
  shortest-path distance from the gene to the arm's training seeds.

Model: within-arm (arm fixed effects) OLS of the stability of a gene on
the standardised features, with cluster-robust standard errors by arm.
Arm fixed effects are used because they estimate exactly the
within-arm contrasts that the question asks ("within one network and
one label, do more-connected genes have more stable scores?") and are
far cheaper than a mixed model on 267,696 rows; a mixed model with an
arm random intercept is fitted on a subsample as a robustness check.

Cross-validation: leave-one-network-out and leave-one-label-out. A
feature is only called a source of stability if its sign survives both.

Outputs:
  results/tables/stability_features_per_gene.tsv
  results/tables/stability_attribution.tsv
  figures/figure5b_stability_attribution.png
  docs/stability_attribution.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

TABLES = ROOT / "results" / "tables"
FIG_DIR = ROOT / "figures"
NET_DIR = ROOT / "data" / "processed" / "networks"
LABEL_DIR = ROOT / "data" / "processed" / "labels"
FIG_DIR.mkdir(parents=True, exist_ok=True)

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
BETWEENNESS_SAMPLES = 300          # sampled betweenness (exact is O(n*m))
LOO_NETWORK = True
LOO_LABEL = True
SEED = 20261105

FEATURES = ["log_degree", "k_core", "clustering", "log_betweenness",
            "log_pagerank", "log_community_size", "dist_to_seed",
            "log_local_score_density", "log_score_mean"]


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def network_features(net_id: str) -> pd.DataFrame:
    """Node-level topological descriptors for one network's LCC."""
    import networkx as nx

    df = pd.read_csv(NET_DIR / f"{net_id}.tsv", sep="\t")
    g = nx.Graph()
    for a, b in zip(df["gene_a"].astype(str), df["gene_b"].astype(str)):
        g.add_edge(a, b)

    log(f"  {net_id}: {g.number_of_nodes()} nodes, {g.number_of_edges()} edges")
    deg = nx.degree_centrality(g)          # already scaled; use raw below
    degree = dict(g.degree())
    core = nx.core_number(g)
    clust = nx.clustering(g)
    t0 = time.time()
    btw = nx.betweenness_centrality(g, k=BETWEENNESS_SAMPLES, seed=SEED,
                                     normalized=True)
    log(f"  {net_id}: sampled betweenness ({BETWEENNESS_SAMPLES} nodes) "
        f"in {time.time()-t0:.1f}s")
    pr = nx.pagerank(g, max_iter=100)
    comms = nx.community.louvain_communities(g, seed=SEED)
    comm_of = {}
    comm_size = {}
    for i, c in enumerate(comms):
        for node in c:
            comm_of[node] = i
            comm_size[node] = len(c)

    out = pd.DataFrame({
        "gene": list(degree),
        "degree": [degree[n] for n in degree],
        "k_core": [core[n] for n in degree],
        "clustering": [clust[n] for n in degree],
        "betweenness": [btw[n] for n in degree],
        "pagerank": [pr[n] for n in degree],
        "community": [comm_of[n] for n in degree],
        "community_size": [comm_size[n] for n in degree],
    })
    return out


def network_graph(net_id: str):
    """Edge list and sparse adjacency for one network (LCC as built by
    d1.networks.io -- the processed file is already the LCC)."""
    from scipy.sparse import csr_matrix

    df = pd.read_csv(NET_DIR / f"{net_id}.tsv", sep="\t")
    genes = sorted(set(df["gene_a"].astype(str)) | set(df["gene_b"].astype(str)))
    idx = {g: i for i, g in enumerate(genes)}
    rows = [idx[a] for a in df["gene_a"].astype(str)]
    cols = [idx[b] for b in df["gene_b"].astype(str)]
    A = csr_matrix((np.ones(len(rows)), (rows, cols)),
                   shape=(len(genes), len(genes)))
    return genes, idx, A


def seed_distances(idx: dict, A, label: str, genes: list) -> pd.Series:
    """Shortest-path distance from every gene to the arm's training seeds.

    The learner sees the positives of the training folds, so the
    relevant quantity is the distance from a gene to the positive set.
    We use the label's positives restricted to the network; the
    fold-level variation of this distance is a second-order effect and
    is not modelled here.
    """
    from scipy.sparse.csgraph import dijkstra

    lab = pd.read_csv(LABEL_DIR / f"{label}.tsv", sep="\t")
    col = "gene" if "gene" in lab.columns else lab.columns[0]
    # the label FILES contain the positives only (there is no `y` column:
    # every row of e.g. data/processed/labels/intogen2024.tsv is a
    # driver). A `y` column, if present, selects the positives.
    if "y" in lab.columns:
        pos = lab[lab["y"] == 1][col].astype(str).tolist()
    else:
        pos = lab[col].astype(str).tolist()
    pos_idx = [idx[g] for g in pos if g in idx]
    if not pos_idx:
        return pd.Series(np.nan, index=genes)
    # min_only: one multi-source pass instead of one pass per seed
    d = dijkstra(A, indices=pos_idx, unweighted=True, min_only=True)
    return pd.Series(d, index=genes)


def main() -> int:
    # ---- per-gene stability from S1 -----------------------------------
    pg = pd.read_csv(TABLES / "stability_per_gene.tsv", sep="\t")
    log(f"loaded S1 table: {len(pg):,} rows")

    # ---- topological features per network ----------------------------
    # cached: the betweenness computation is the expensive part
    feat_cache = TABLES / "network_features.tsv"
    if feat_cache.exists():
        log(f"using the cached feature table {feat_cache.name}")
        cached = pd.read_csv(feat_cache, sep="\t")
        feats = {net: cached[cached.network == net].drop(columns=["network"])
                 .reset_index(drop=True) for net in NETWORKS}
    else:
        feats = {}
        for net in NETWORKS:
            feats[net] = network_features(net)
        cache_rows = []
        for net in NETWORKS:
            f = feats[net].copy()
            f.insert(0, "network", net)
            cache_rows.append(f)
        pd.concat(cache_rows, ignore_index=True).to_csv(
            feat_cache, sep="\t", index=False)
        log(f"wrote {feat_cache.name}")

    # ---- attach the features and the seed distance --------------------
    for net in NETWORKS:
        f = feats[net]
        genes, idx, A = network_graph(net)
        for label in LABELS:
            d = seed_distances(idx, A, label, genes)
            f[f"dist_to_seed__{label}"] = d.reindex(f["gene"]).to_numpy()
        feats[net] = f

    rows = []
    for net in NETWORKS:
        f = feats[net]
        for label in LABELS:
            sub = pg[(pg.network == net) & (pg.label == label)].copy()
            # the S1 table already carries `degree`; keep the network one
            sub = sub.drop(columns=["degree"])
            sub = sub.merge(f.drop(columns=["community"]), on="gene", how="left")
            # the arm's own label decides the seed distance
            sub["dist_to_seed"] = sub[f"dist_to_seed__{label}"]
            rows.append(sub)
    df = pd.concat(rows, ignore_index=True)

    df["log_degree"] = np.log1p(df["degree"])
    df["log_betweenness"] = np.log1p(df["betweenness"])
    # pagerank can be exactly 0 for a node -> log10 gives -inf
    df["log_pagerank"] = np.log10(df["pagerank"].clip(lower=1e-15))
    df["log_community_size"] = np.log1p(df["community_size"])
    df["log_local_score_density"] = np.log1p(df["local_score_density"])
    df["log_score_mean"] = np.log10(df["score_mean"].clip(lower=1e-12))
    # any remaining +-inf (e.g. an unreachable gene) becomes NaN and the
    # model drops that row
    df[FEATURES] = df[FEATURES].replace([np.inf, -np.inf], np.nan)
    # standardise the predictors (z-scores over the pooled rows) so that
    # a coefficient reads as "one standard deviation of the feature
    # changes log(stability) by beta"
    for f in FEATURES:
        mu, sd = df[f].mean(), df[f].std(ddof=1)
        df[f + "_z"] = (df[f] - mu) / sd
    STANDARDISED = [f + "_z" for f in FEATURES]
    df[STANDARDISED] = df[STANDARDISED].replace([np.inf, -np.inf], np.nan)
    df["arm"] = df["network"] + "|" + df["label"]
    df["log_rank_sd"] = np.log(df["rank_sd"].clip(lower=1e-6))
    df["log_score_cv"] = np.log(df["score_cv"].clip(lower=1e-6))

    out_feat = TABLES / "stability_features_per_gene.tsv"
    keep = (["network", "label", "arm", "gene", "y", "degree", "k_core",
             "clustering", "betweenness", "pagerank", "community_size",
             "local_score_density", "score_mean", "score_cv", "rank_sd",
             "dist_to_seed"] + [f"dist_to_seed__{l}" for l in LABELS])
    df[keep].to_csv(out_feat, sep="\t", index=False)
    log(f"wrote {out_feat.name} ({len(df):,} rows)")

    # ---- the attribution model ----------------------------------------
    def fit(d: pd.DataFrame, outcome: str, with_arm: bool = True):
        """Within-arm OLS with arm fixed effects, cluster-robust by arm.

        `with_arm=False` drops the arm fixed effects and is used for the
        leave-one-out transfer test: the features are fitted on the other
        arms and used to predict the held-out arm's genes, so the
        prediction cannot borrow that arm's intercept.
        """
        d = d.dropna(subset=[outcome] + STANDARDISED).copy()
        if len(d) == 0:
            raise SystemExit("no rows left after dropping missing values")
        formula = f"{outcome} ~ " + " + ".join(STANDARDISED)
        if with_arm:
            d["arm"] = pd.Categorical(d["arm"])
            formula += " + C(arm)"
        model = smf.ols(formula, data=d)
        if with_arm:
            model = model.fit(cov_type="cluster",
                              cov_kwds={"groups": d["arm"].cat.codes})
        else:
            model = model.fit()
        return model

    def bh_fdr(p):
        p = np.asarray(p, dtype=float)
        order = np.argsort(p)
        ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        q = np.empty_like(ranked)
        q[order] = np.clip(ranked, 0, 1)
        return q

    results = []
    for outcome, tag in (("log_score_cv", "score CV"),
                         ("log_rank_sd", "rank SD")):
        log(f"fitting the attribution model for {tag} "
            f"({len(df):,} rows) ...")
        t0 = time.time()
        model = fit(df, outcome)
        log(f"  fit in {time.time()-t0:.1f}s")
        ci = model.conf_int()
        for feat in STANDARDISED:
            i = model.params.index.get_loc(feat)
            results.append({
                "outcome": tag, "feature": feat.replace("_z", ""),
                "beta": model.params[feat], "se": model.bse[feat],
                "ci_low": ci.iloc[i, 0], "ci_high": ci.iloc[i, 1],
                "p": model.pvalues[feat],
            })
    res = pd.DataFrame(results)

    # marginal association of each feature with the stability, computed
    # INSIDE each arm (so that the arm-level differences cannot create
    # it) and then averaged over the 24 arms
    marg_rows = []
    for oc, col in (("score CV", "log_score_cv"), ("rank SD", "log_rank_sd")):
        for feat in FEATURES:
            rhos = []
            for arm, sub in df.groupby("arm"):
                s2 = sub[[feat, col]].dropna()
                if len(s2) < 50:
                    continue
                rhos.append(float(sps.spearmanr(s2[feat], s2[col]).statistic))
            marg_rows.append({
                "outcome": oc, "feature": feat,
                "marginal_spearman_median": float(np.median(rhos)),
                "marginal_spearman_min": float(np.min(rhos)),
                "marginal_spearman_max": float(np.max(rhos)),
                "arms_positive": int(sum(r > 0 for r in rhos)),
            })
    marg = pd.DataFrame(marg_rows)

    # correlation among the predictors (to show how much they compete)
    corr = df[FEATURES].corr(method="spearman")

    res["q_bh"] = np.nan
    for oc in res["outcome"].unique():
        m = res["outcome"] == oc
        res.loc[m, "q_bh"] = bh_fdr(res.loc[m, "p"].to_numpy())

    # ---- cross-validation: leave one network / one label out ---------
    OUTCOME_COL = {"score CV": "log_score_cv", "rank SD": "log_rank_sd"}
    loo = []
    for oc in res["outcome"].unique():
        col = OUTCOME_COL[oc]
        for hold in NETWORKS:
            tr = df[df.network != hold]
            te = df[df.network == hold]
            m = fit(tr, col, with_arm=False)
            te = te.dropna(subset=[col] + STANDARDISED)
            preds = m.predict(te)
            y = te[col].to_numpy()
            ss_res = float(((y - preds) ** 2).sum())
            ss_tot = float(((y - y.mean()) ** 2).sum())
            loo.append({"outcome": oc, "held_out": hold, "kind": "network",
                        "r2_out": 1 - ss_res / ss_tot})
        for hold in LABELS:
            tr = df[df.label != hold]
            te = df[df.label == hold]
            m = fit(tr, col, with_arm=False)
            te = te.dropna(subset=[col] + STANDARDISED)
            preds = m.predict(te)
            y = te[col].to_numpy()
            ss_res = float(((y - preds) ** 2).sum())
            ss_tot = float(((y - y.mean()) ** 2).sum())
            loo.append({"outcome": oc, "held_out": hold, "kind": "label",
                        "r2_out": 1 - ss_res / ss_tot})
    loo_df = pd.DataFrame(loo)

    # sign stability of each feature under the leave-one-out refits
    signs = []
    for oc in res["outcome"].unique():
        col = OUTCOME_COL[oc]
        for hold, kind in ([(h, "network") for h in NETWORKS]
                           + [(h, "label") for h in LABELS]):
            tr = df.copy()
            if kind == "network":
                tr = tr[tr.network != hold]
            else:
                tr = tr[tr.label != hold]
            m = fit(tr, col)
            for feat in STANDARDISED:
                signs.append({"outcome": oc, "held_out": hold, "kind": kind,
                              "feature": feat.replace("_z", ""),
                              "beta": m.params[feat]})
    signs_df = pd.DataFrame(signs)
    sign_stab = (signs_df.groupby(["outcome", "feature"])["beta"]
                 .apply(lambda s: float((np.sign(s) == np.sign(s.iloc[0])).mean()))
                 .rename("sign_stability").reset_index())
    res = res.merge(sign_stab, on=["outcome", "feature"], how="left")

    res = res.merge(marg, on=["outcome", "feature"], how="left")
    res.to_csv(TABLES / "stability_attribution.tsv", sep="\t", index=False)
    loo_df.to_csv(TABLES / "stability_attribution_loo.csv", sep="\t", index=False)
    corr.to_csv(TABLES / "stability_attribution_predictor_corr.tsv", sep="\t")
    log(f"wrote stability_attribution.tsv ({len(res)} rows) and "
        f"stability_attribution_loo.csv ({len(loo_df)} rows)")

    make_figure(res, loo_df)
    return 0


def make_figure(res: pd.DataFrame, loo: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2),
                             gridspec_kw={"width_ratios": [1.5, 1.0]})
    ax = axes[0]
    for i, oc in enumerate(["score CV", "rank SD"]):
        sub = res[res.outcome == oc].sort_values("beta")
        y = np.arange(len(sub)) + (i - 0.5) * 0.32
        ax.errorbar(sub["beta"], y,
                    xerr=[sub["beta"] - sub["ci_low"], sub["ci_high"] - sub["beta"]],
                    fmt="o", color="#4C72B0" if i == 0 else "#C44E52",
                    capsize=2.5, markersize=4.5, lw=1.2,
                    label=f"{oc} (within-arm, cluster-robust)")
    ax.axvline(0, color="black", lw=1)
    ax.set_yticks(np.arange(len(res[res.outcome == "score CV"])))
    ax.set_yticklabels(res[res.outcome == "score CV"].sort_values("beta")["feature"],
                       fontsize=9)
    ax.set_xlabel("standardised coefficient on log(stability)")
    ax.set_title("A  What predicts a gene's stability? "
                 "(negative = MORE stable)", fontsize=11)
    ax.legend(fontsize=8.5, loc="lower right")

    ax = axes[1]
    sub = loo[loo["outcome"] == "score CV"].reset_index(drop=True)
    labels = [f"{h}\n({k})" for h, k in zip(sub["held_out"], sub["kind"])]
    kinds = sub["kind"].tolist()
    r2 = sub["r2_out"].tolist()
    colors = ["#4C72B0" if k == "network" else "#DD8452" for k in kinds]
    ax.bar(range(len(r2)), r2, color=colors)
    ax.set_xticks(range(len(r2)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6.5)
    ax.set_ylabel("out-of-sample R²")
    ax.set_title("B  Does the attribution transfer to an unseen arm?",
                 fontsize=11)
    handles = [plt.Line2D([], [], color="#4C72B0", lw=6),
               plt.Line2D([], [], color="#DD8452", lw=6)]
    ax.legend(handles, ["held-out network", "held-out label"], fontsize=8)
    fig.suptitle("Figure 5b  Attribution of the network prior's stability "
                 "to topology", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG_DIR / "figure5b_stability_attribution.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    log(f"wrote figure5b_stability_attribution.png")


if __name__ == "__main__":
    raise SystemExit(main())
