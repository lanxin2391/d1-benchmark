"""S1 -- stability atlas: how stable is the network prior, gene by gene?

The question (teacher, 2026-11-04): what is the QUANTITATIVE SOURCE of
the network's stability? The proteomics paper this came from
(Klingelhuber et al., Nat Metab 2024) treats the stability of the
measurement layer as a first-class quantity: it reports the per-protein
coefficient of variation, the replicate correlation, and the
distribution of both, and then asks which biological features explain
the stable core.

Our analogue: the "measurement" here is the RWR prior's score for a
gene. For every (network, label) arm we have 40 repeats x 5 folds of
per-gene scores (results/scores_v40/*.npz), so we can ask how much a
gene's score moves when the universe is re-partitioned -- and then
(S2) which topological features predict that movement.

Definitions (per arm, alpha = 0.5):

  evaluated score   for every gene and every repeat, the score of the
                    ONE fold in which that gene was in the test set.
                    This is the score that actually entered the AUROC.
  all-fold score    the mean over the 5 folds of a repeat (the
                    repeat-level estimate used by the primary CI).
  score SD / CV     across the 40 repeats (sampling error of the score)
  percentile rank   within each repeat, the rank of the gene among all
                    genes of the universe, scaled to [0, 1]
  rank SD           across the 40 repeats -- scale-free, and the
                    quantity that says "does this gene keep its place
                    in the ranking?"
  top-k hits        in how many of the 40 repeats the gene is inside
                    the top 10 / 50 / 1 % of the universe

Outputs:
  results/tables/stability_per_gene.tsv    (24 arms x n_genes rows)
  results/tables/stability_by_arm.tsv      (24 rows of summaries)
  figures/figure5_stability_atlas.png
  docs/stability_atlas.md

Usage:
  E:\\anaconda\\python.exe scripts/run_stability_atlas.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SCORES = ROOT / "results" / "scores_v40"
TABLES = ROOT / "results" / "tables"
FIG_DIR = ROOT / "figures"
DOC = ROOT / "docs" / "stability_atlas.md"
TABLES.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

ALPHA = 0.5
NET_ORDER = ["funmap", "rna_coexp", "string_phys700", "intact",
             "string_full700", "reactome"]
LABEL_ORDER = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
PROVENANCE = {
    "funmap": "protein-derived", "string_phys700": "protein-derived",
    "intact": "protein-derived", "rna_coexp": "rna-derived",
    "string_full700": "curated-literature", "reactome": "curated-literature",
}


def load_arm(net: str, label: str):
    fp = SCORES / f"{net}__{label}.npz"
    d = np.load(fp, allow_pickle=True)
    alphas = list(d["alphas"])
    ai = int(np.argmin(np.abs(np.asarray(alphas, dtype=float) - ALPHA)))
    rwr = np.asarray(d["rwr"][ai], dtype=np.float64)      # (n_genes, 200)
    repeat_idx = np.asarray(d["repeat_idx"]).astype(int)
    fold_idx = np.asarray(d["fold_idx"]).astype(int)
    in_test = np.asarray(d["in_test"], dtype=bool)         # (200, n_genes)
    genes = np.asarray(d["genes"]).astype(str)
    y = np.asarray(d["y"]).astype(int)
    degree = np.asarray(d["degree"], dtype=float)
    n_rep = int(repeat_idx.max()) + 1

    # evaluated score: the fold in which the gene was in the test set
    eval_scores = np.full((len(genes), n_rep), np.nan)
    for r in range(n_rep):
        for j in np.where(repeat_idx == r)[0]:
            m = in_test[j]
            eval_scores[m, r] = rwr[m, j]
    # all-fold mean per repeat (the repeat-level estimate)
    allfold = np.empty_like(eval_scores)
    for r in range(n_rep):
        js = np.where(repeat_idx == r)[0]
        allfold[:, r] = rwr[:, js].mean(axis=1)

    return dict(genes=genes, y=y, degree=degree, eval_scores=eval_scores,
                allfold=allfold, n_rep=n_rep)


def arm_metrics(arm: dict) -> pd.DataFrame:
    """Per-gene stability metrics.

    The rank metrics use the SHARED engine (d1/engine/stability.py), so
    that the D-07 tie-break (score descending, then symbol ascending) is
    the one used everywhere else in the project. Before the merge this
    script broke score ties by index order, which made the median rank
    SD differ from the one B's engine produces (0.0267 vs 0.027715).
    """
    from d1.engine.stability import crowding_counts, percentile_ranks
    S = arm["eval_scores"]
    n_genes, n_rep = S.shape
    mean = S.mean(axis=1)
    sd = S.std(axis=1, ddof=1)
    genes = arm["genes"]
    ranks = np.empty_like(S, dtype=float)
    for r in range(n_rep):
        ranks[:, r] = percentile_ranks(S[:, r], genes)
    rank_mean = ranks.mean(axis=1)
    rank_sd = ranks.std(axis=1, ddof=1)
    rank_range = ranks.max(axis=1) - ranks.min(axis=1)
    # crowding, two definitions: B's Scott-bandwidth count and the
    # +-10 % score band. They disagree in SIGN (B's attribution finds
    # the Scott-bandwidth crowding negatively correlated with rank
    # instability in 21 of 24 arms), so both are reported.
    crowd = crowding_counts(mean)
    s_sorted = np.sort(mean)
    lo = np.searchsorted(s_sorted, mean * 0.9, side="left")
    hi = np.searchsorted(s_sorted, mean * 1.1, side="right")
    density = hi - lo

    top10 = (ranks <= 10 / (n_genes - 1)).sum(axis=1)
    top50 = (ranks <= 50 / (n_genes - 1)).sum(axis=1)
    top1pct = (ranks <= 0.01).sum(axis=1)
    top5pct = (ranks <= 0.05).sum(axis=1)

    cv = np.divide(sd, mean, out=np.full_like(sd, np.nan),
                   where=np.abs(mean) > 1e-12)

    return pd.DataFrame({
        "gene": arm["genes"], "y": arm["y"], "degree": arm["degree"],
        "score_mean": mean, "score_sd": sd, "score_cv": cv,
        "rank_mean": rank_mean, "rank_sd": rank_sd, "rank_range": rank_range,
        "top10_hits": top10, "top50_hits": top50,
        "top1pct_hits": top1pct, "top5pct_hits": top5pct,
        "crowding_scott": crowd,
        "local_score_density": density,
    })


def arm_level_stability() -> pd.DataFrame:
    """Arm-level stability: how much does the arm's own endpoint move
    when the universe is re-partitioned?

    This is the number that maps most directly onto "how stable is this
    network's contribution": the SD of the 40 repeat means of
    delta-AUROC, and the intra-repeat correlation of the fold-level
    deltas (negative correlation between the folds of one repeat means
    the five folds of a repeat are not five independent observations).
    """
    a4 = pd.read_csv(ROOT / "results" / "runs_v40" / "a4_combined.tsv",
                     sep="\t")
    sub = a4[(a4["method"] == "degree") | (a4["alpha"] == ALPHA)]
    sub = sub[sub["method"].isin(["rwr", "degree"])]
    piv = (sub.pivot_table(index=["network", "label_set", "repeat", "fold"],
                           columns="method", values="auroc")
           .reset_index().dropna(subset=["rwr", "degree"]))
    piv["delta"] = piv["rwr"] - piv["degree"]
    out = []
    for (net, lab), g in piv.groupby(["network", "label_set"]):
        rep = g.groupby("repeat")["delta"].mean()
        # one-way ICC(1) of the fold-level deltas
        grp = g.groupby("repeat")["delta"]
        k = grp.size().iloc[0]
        grand = g["delta"].mean()
        msb = (grp.mean() - grand).pow(2).sum() * k / (len(grp) - 1)
        msw = grp.apply(lambda s: ((s - s.mean()) ** 2).sum()).sum() / (len(g) - len(grp))
        icc = (msb - msw) / (msb + (k - 1) * msw)
        out.append({
            "network": net, "label": lab,
            "arm_delta_mean": float(rep.mean()),
            "arm_delta_sd_repeat": float(rep.std(ddof=1)),
            "arm_delta_se_repeat": float(rep.std(ddof=1) / np.sqrt(len(rep))),
            "arm_repeat_range": float(rep.max() - rep.min()),
            "arm_icc_intra_repeat": float(icc),
        })
    return pd.DataFrame(out)


def make_figure(per_gene: pd.DataFrame, by_arm: pd.DataFrame,
                 prefix: str = "") -> Path:
    """Figure 5 -- the stability atlas."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    C_BG, C_POS, C_LINE, C_ACC = "#4C72B0", "#C44E52", "#999999", "#DD8452"
    prov_color = {"protein-derived": "#1f77b4", "rna-derived": "#2ca02c",
                  "curated-literature": "#d62728"}

    pg = per_gene.copy()
    # score decile within each arm. pd.qcut on the rank is used (not
    # int(pct*10), which lumps the bottom 20 % into one bin).
    pg["score_dec"] = pg.groupby(["network", "label"])["score_mean"] \
        .transform(lambda s: pd.qcut(s.rank(method="first"), 10,
                                     labels=False, duplicates="drop") + 1)

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.5))

    # ---- A: U-shape, pooled, background vs drivers --------------------
    ax = axes[0, 0]
    for y_val, colour, lab in ((0, C_BG, "background genes"),
                               (1, C_POS, "driver genes")):
        sub = pg[pg.y == y_val]
        med = sub.groupby("score_dec")["rank_sd"].median()
        ax.plot(med.index, med.values, "o-", color=colour, label=lab, lw=2)
    ax.set_xlabel("score-magnitude decile (within arm)")
    ax.set_ylabel("median rank SD (across 40 repeats)")
    ax.set_title("A  Rank instability peaks in the crowded middle", fontsize=11)
    ax.legend(fontsize=9)
    ax.text(0.02, 0.95, "rank SD is bounded at the top and the bottom\n"
                        "(a gene pinned near an edge cannot move much)",
            transform=ax.transAxes, fontsize=8, va="top", color="grey")

    # ---- B: the U-shape holds in every arm ----------------------------
    ax = axes[0, 1]
    for (net, lab), sub in pg[pg.y == 0].groupby(["network", "label"]):
        med = sub.groupby("score_dec")["rank_sd"].median()
        ax.plot(med.index, med.values, "-", color=C_LINE, lw=0.8, alpha=0.7)
    sub = pg[pg.y == 0]
    med = sub.groupby("score_dec")["rank_sd"].median()
    ax.plot(med.index, med.values, "k-", lw=2.2, label="pooled (24 arms)")
    ax.set_xlabel("score-magnitude decile (within arm)")
    ax.set_ylabel("median rank SD")
    ax.set_title("B  The shape is the same in all 24 arms", fontsize=11)
    ax.legend(fontsize=9)

    # ---- C: arm-level stability ---------------------------------------
    ax = axes[1, 0]
    st = by_arm.copy()
    st["label_short"] = st["label"].replace({
        "intogen2024": "intogen2024", "ot_all": "ot_all",
        "ot_nolit": "ot_nolit",
        "intogen_temporal_new": "temporal_new"})
    order = [f"{n} | {l}" for n in NET_ORDER for l in LABEL_ORDER]
    st["key"] = st["network"] + " | " + st["label"]
    st = st.set_index("key").reindex([o for o in order if o in set(st["key"])])
    ax.bar(range(len(st)), st["arm_delta_sd_repeat"].to_numpy(),
           color=[prov_color[p] for p in st["provenance"]])
    ax.set_xticks(range(len(st)))
    ax.set_xticklabels(st.index, rotation=90, fontsize=7.5)
    ax.set_ylabel("SD of delta-AUROC across the 40 repeats")
    ax.set_title("C  Arm-level stability: how much the endpoint itself moves",
                 fontsize=11)
    handles = [plt.Line2D([], [], color=prov_color[k], lw=6)
               for k in prov_color]
    ax.legend(handles, list(prov_color), fontsize=8, loc="upper left",
              title="provenance", title_fontsize=8)

    # ---- D: the mechanism -- local score density ----------------------
    ax = axes[1, 1]
    sub = pg[(pg.y == 0) & (pg.local_score_density > 0)]
    ax.scatter(np.log10(sub.local_score_density), sub.rank_sd, s=2,
               alpha=0.05, color=C_BG, rasterized=True)
    rho = by_arm["spearman_ranked_sd_vs_local_density"]
    rho_cv = by_arm["spearman_cv_vs_local_density"]
    ax.set_xlabel("log10 number of genes within +-10 % of this gene's score")
    ax.set_ylabel("rank SD")
    ax.set_title("D  Crowding depends on the window; the gene's own\n"
                 "     score noise is the robust predictor",
                 fontsize=11)
    ax.text(0.02, 0.95,
            f"per-arm Spearman(rank SD, local density):\n"
            f"  median {rho.median():+.2f}, range {rho.min():+.2f} .. "
            f"{rho.max():+.2f}  (positive in {int((rho > 0).sum())}/24)\n"
            f"per-arm Spearman(score CV, local density):\n"
            f"  median {rho_cv.median():+.2f}, positive in "
            f"{int((rho_cv > 0).sum())}/24",
            transform=ax.transAxes, fontsize=8.5, va="top", family="monospace")

    fig.suptitle("Figure 5  Stability atlas of the network prior "
                 "(per-gene ranking stability across 40 re-partitions)",
                 fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = FIG_DIR / f"figure5_stability_atlas{prefix}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores-dir", default=None,
                    help="directory holding the per-gene score dumps "
                         "(default results/scores_v40)")
    ap.add_argument("--networks", nargs="+", default=None)
    ap.add_argument("--labels", nargs="+", default=None)
    ap.add_argument("--prefix", default="",
                    help="prefix for the output table names, e.g. '_clingen'")
    args = ap.parse_args()

    global SCORES, NET_ORDER, LABEL_ORDER
    if args.scores_dir:
        SCORES = ROOT / args.scores_dir
    if args.networks:
        NET_ORDER = list(args.networks)
    if args.labels:
        LABEL_ORDER = list(args.labels)
    prefix = args.prefix

    per_gene_rows = []
    arm_rows = []

    for net in NET_ORDER:
        for label in LABEL_ORDER:
            fp = SCORES / f"{net}__{label}.npz"
            if not fp.exists():
                print(f"  [skip] {net} x {label}: no score dump")
                continue
            arm = load_arm(net, label)
            df = arm_metrics(arm)
            df.insert(0, "label", label)
            df.insert(0, "network", net)
            df["provenance"] = PROVENANCE[net]
            per_gene_rows.append(df)

            pos = df["y"] == 1
            bg = ~pos
            # stability of the all-fold variant, for the sensitivity check
            S2 = arm["allfold"]
            n_genes = S2.shape[0]
            o2 = np.argsort(np.argsort(S2, axis=0), axis=0) / (n_genes - 1)
            rank_sd2 = o2.std(axis=1, ddof=1)

            def sp(a, b):
                a = np.asarray(a, dtype=float)
                b = np.asarray(b, dtype=float)
                if len(a) < 5 or np.all(a == a[0]) or np.all(b == b[0]):
                    return np.nan
                return float(sps.spearmanr(a, b).statistic)

            arm_rows.append({
                "network": net, "label": label,
                "provenance": PROVENANCE[net],
                "n_genes": int(len(df)), "n_positives": int(pos.sum()),
                "n_repeats": int(arm["n_rep"]),
                # background
                "bg_rank_sd_median": float(df.loc[bg, "rank_sd"].median()),
                "bg_rank_sd_iqr": float(df.loc[bg, "rank_sd"].quantile(.75)
                                        - df.loc[bg, "rank_sd"].quantile(.25)),
                "bg_score_cv_median": float(df.loc[bg, "score_cv"].median()),
                "bg_top1pct_hits_mean": float(df.loc[bg, "top1pct_hits"].mean()),
                "bg_top5pct_hits_mean": float(df.loc[bg, "top5pct_hits"].mean()),
                "bg_frac_top1pct_ever": float((df.loc[bg, "top1pct_hits"] > 0).mean()),
                # positives
                "pos_rank_sd_median": float(df.loc[pos, "rank_sd"].median()),
                "pos_rank_sd_iqr": float(df.loc[pos, "rank_sd"].quantile(.75)
                                         - df.loc[pos, "rank_sd"].quantile(.25)),
                "pos_score_cv_median": float(df.loc[pos, "score_cv"].median()),
                "pos_top1pct_hits_mean": float(df.loc[pos, "top1pct_hits"].mean()),
                "pos_top5pct_hits_mean": float(df.loc[pos, "top5pct_hits"].mean()),
                "pos_frac_top1pct_ever": float((df.loc[pos, "top1pct_hits"] > 0).mean()),
                # how much more stable are the drivers than the background?
                "rank_sd_advantage": float(df.loc[bg, "rank_sd"].median()
                                           - df.loc[pos, "rank_sd"].median()),
                "top1pct_hits_ratio": float(df.loc[pos, "top1pct_hits"].mean()
                                            / max(df.loc[bg, "top1pct_hits"].mean(), 1e-9)),
                # what correlates with stability (background genes)
                "spearman_ranked_sd_vs_degree": sp(df.loc[bg, "rank_sd"],
                                                   df.loc[bg, "degree"]),
                "spearman_ranked_sd_vs_scoremean": sp(df.loc[bg, "rank_sd"],
                                                      df.loc[bg, "score_mean"]),
                "spearman_ranked_sd_vs_degree_pos": sp(df.loc[pos, "rank_sd"],
                                                       df.loc[pos, "degree"]),
                # the crowding signal, TWO definitions: the relative
                # +-10 % score band and the Scott bandwidth. They
                # correlate with the rank instability with OPPOSITE
                # signs, which is why "crowding" is not a mechanism --
                # the sign depends on the window.
                "spearman_ranked_sd_vs_local_density": sp(df.loc[bg, "rank_sd"],
                                                          df.loc[bg, "local_score_density"]),
                "spearman_ranked_sd_vs_crowding_scott": sp(
                    df.loc[bg, "rank_sd"], df.loc[bg, "crowding_scott"]),
                # the robust predictor: the gene's OWN score noise
                "spearman_ranked_sd_vs_score_cv": sp(df.loc[bg, "rank_sd"],
                                                      df.loc[bg, "score_cv"]),
                "spearman_cv_vs_local_density": sp(df.loc[bg, "score_cv"],
                                                   df.loc[bg, "local_score_density"]),
                "spearman_cv_vs_crowding_scott": sp(df.loc[bg, "score_cv"],
                                                     df.loc[bg, "crowding_scott"]),
                "spearman_cv_vs_degree": sp(df.loc[bg, "score_cv"],
                                            df.loc[bg, "degree"]),
                # sensitivity: evaluated score vs all-fold mean
                "spearman_ranked_sd_evaluated_vs_allfold": sp(
                    df["rank_sd"], rank_sd2),
            })
            print(f"  [done] {net} x {label}", flush=True)

    per_gene = pd.concat(per_gene_rows, ignore_index=True)
    by_arm = pd.DataFrame(arm_rows)
    # arm-level stability (repeat-level movement of the endpoint itself)
    arm_stab = arm_level_stability()
    by_arm = by_arm.merge(arm_stab, on=["network", "label"], how="left")
    per_gene.to_csv(TABLES / f"stability_per_gene{prefix}.tsv", sep="\t",
                    index=False, lineterminator="\n")
    by_arm.to_csv(TABLES / f"stability_by_arm{prefix}.tsv", sep="\t",
                  index=False, lineterminator="\n")
    print(f"\nwrote {TABLES/f'stability_per_gene{prefix}.tsv'} "
          f"({len(per_gene)} rows)")
    print(f"wrote {TABLES/f'stability_by_arm{prefix}.tsv'} ({len(by_arm)} rows)")
    fig = make_figure(per_gene, by_arm)
    print(f"wrote {fig}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
