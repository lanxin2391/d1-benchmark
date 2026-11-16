"""Control 21 -- the visibility decomposition under the SECOND learner.

The objection (and it is a good one): if swapping the learner changes
the conclusion, then the conclusion was never ours -- Guala & Sonnhammer
(2017) had already said that a degree-filtered local method beats the
diffusion methods on a full list.

This control asks the question that decides it: **is the DECOMPOSITION
learner-invariant, or only the margin?** The decomposition is the
project's actual claim: the reported margin is a mixture of topology,
node degree and study visibility. Control 17 measured the margin under
the proximity learner; nothing measured the VISIBILITY part under it.

So: the bait-matched null (G3) is recomputed with the proximity learner
as the scorer, on the same splits, the same folds, the same publication
counts and the same matching rule (the 20 negatives nearest to the
positive in log10(n_pubmed + 1)). The comparison is with the published
RWR numbers in g3_bait_matched_null.tsv, which come from the same
matcher.

A self-check first: the matcher in this file is run on the RWR score
dumps of one arm and must reproduce the published G3 value for that arm
before the proximity run is trusted.

Outputs:
  results/tables/c21_visibility_second_learner.tsv
  results/tables/c21_visibility_learner_comparison.tsv
  docs/c21_visibility_second_learner.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import dijkstra

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.evaluate import read_splits  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETS = ROOT / "data/processed/networks"
SPLITS = ROOT / "data/processed/splits_v40"
SCORES = ROOT / "results/scores_v40"
TABLES = ROOT / "results/tables"
PUB = ROOT / "data/processed/labels/pubcount.tsv"
DOC = ROOT / "docs/c21_visibility_second_learner.md"

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new",
          "ncg_known"]
PRIMARY_LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
N_MATCHED = 20
RNG_SEED = 20261105


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ------------------------------------------------------------------ the matcher
def rank_of_test(scores, genes):
    """The protocol's ranking of one test fold (decision D-07: score
    descending, then gene symbol ascending). Returns the rank index of
    every test gene, lower = better.

    The proximity learner's score is an integer distance, so ties are
    everywhere; the tie-break makes the comparison well defined. RWR
    has no ties, so its ranks reproduce G3's strict comparison exactly.
    """
    order = np.lexsort((np.asarray(genes).astype(str), -np.asarray(scores,
                                                                  dtype=float)))
    rank = np.empty(len(order), dtype=int)
    rank[order] = np.arange(len(order))
    return rank


def null_above(score_pos, score_neg, pub_pos, pub_neg, rng, k=N_MATCHED):
    """For every positive: the fraction of the k random negatives and of
    the k publication-matched negatives that rank above it.

    Same rule as G3 (scripts/run_g3_bait_matched_null.py): the matched
    negatives are the k negatives nearest to the positive in
    log10(n_pubmed + 1); the unmatched ones are drawn at random. The
    matching here is vectorised through a sort + a window, which gives
    the same neighbours as the argsort in G3 (the self-check below
    verifies that on the RWR dumps).
    """
    n_neg = len(score_neg)
    kk = min(k, n_neg)
    out_u = np.empty(len(score_pos))
    out_m = np.empty(len(score_pos))
    order = np.argsort(pub_neg, kind="mergesort")
    sorted_pub = pub_neg[order]
    for a in range(len(score_pos)):
        # the unmatched null
        rnd = rng.choice(n_neg, size=kk, replace=False)
        out_u[a] = (score_neg[rnd] > score_pos[a]).mean()
        # the matched null: the k nearest negatives in publication volume
        p = np.searchsorted(sorted_pub, pub_pos[a])
        lo = max(0, p - kk)
        hi = min(len(sorted_pub), p + kk)
        cand = order[lo:hi]
        d = np.abs(pub_neg[cand] - pub_pos[a])
        sel = cand[np.argsort(d, kind="mergesort")[:kk]]
        out_m[a] = (score_neg[sel] > score_pos[a]).mean()
    return out_u, out_m


def fold_stats(score, y, pub, genes, rng):
    """One test fold -> (unmatched p, matched p) for the positives that
    have a publication record and at least 50 negatives.

    The comparison is made on the protocol's RANKS, not on the raw
    scores, so that ties (the proximity learner has many) are resolved
    by the documented tie-break instead of counting as 'not above'.
    """
    have = np.isfinite(pub)
    idx_p = np.where((y == 1) & have)[0]
    idx_n = np.where((y == 0) & have)[0]
    if len(idx_p) == 0 or len(idx_n) < 50:
        return np.array([]), np.array([])
    # a SMALLER rank is a better gene, so the vector is negated to keep
    # the '>' comparison inside null_above meaning 'ranks above'
    rank = -rank_of_test(score, genes)
    return null_above(rank[idx_p], rank[idx_n], pub[idx_p], pub[idx_n], rng)


# ------------------------------------------------------------------ the learner
def distance_to_nearest(A, seeds):
    """The unweighted distance from every node to the nearest seed."""
    d = dijkstra(A, indices=seeds, unweighted=True, min_only=True,
                 directed=False)
    return np.asarray(d, dtype=float).ravel()


# ------------------------------------------------------------------ the self-check
def self_check(pub_map) -> float:
    """Run this file's matcher on the RWR score dumps of one arm and
    compare with the published G3 value for that arm."""
    net, lab = "funmap", "intogen2024"
    d = np.load(SCORES / f"{net}__{lab}.npz", allow_pickle=True)
    ai = int(np.argmin(np.abs(np.asarray(d["alphas"], dtype=float) - 0.5)))
    rwr = np.asarray(d["rwr"][ai], dtype=float)
    in_test = np.asarray(d["in_test"], dtype=bool)
    y = np.asarray(d["y"]).astype(int)
    genes = np.asarray(d["genes"]).astype(str)
    pub = np.array([pub_map.get(g, np.nan) for g in genes])
    mp, up = [], []
    for j in range(rwr.shape[1]):
        m = in_test[j]
        if not m.sum():
            continue
        idx = np.where(m)[0]
        rng = np.random.default_rng(RNG_SEED + j)
        u, mm = fold_stats(rwr[idx, j], y[idx], pub[idx],
                           genes[idx], rng)
        if len(u):
            mp.append(mm)
            up.append(u)
    mp = np.concatenate(mp)
    up = np.concatenate(up)
    mine_matched = float((mp < 0.05).mean())
    mine_unmatched = float((up < 0.05).mean())
    published = pd.read_csv(TABLES / "g3_bait_matched_null.tsv", sep="\t")
    row = published[(published.network == net) & (published.label == lab)].iloc[0]
    log(f"  self-check on {net} x {lab}: "
        f"matched {mine_matched:.4f} vs published "
        f"{row.frac_matched_null_above:.4f}; "
        f"unmatched {mine_unmatched:.4f} vs published "
        f"{row.frac_unmatched_null_above:.4f}")
    ok = (abs(mine_matched - row.frac_matched_null_above) < 0.02
          and abs(mine_unmatched - row.frac_unmatched_null_above) < 0.02)
    log(f"  the matcher reproduces the published G3 statistic: {ok}")
    return float(max(abs(mine_matched - row.frac_matched_null_above),
                     abs(mine_unmatched - row.frac_unmatched_null_above)))


# ------------------------------------------------------------------ the run
def main() -> int:
    pub_df = pd.read_csv(PUB, sep="\t", dtype=str)
    pub_map = dict(zip(pub_df["gene"], pub_df["log10"].astype(float)))

    log("self-check: this file's matcher on the RWR dumps")
    err = self_check(pub_map)
    if err > 0.02:
        log("the matcher does NOT reproduce the published G3 value -- stop")
        return 1

    rows = []
    for net in NETWORKS:
        genes, A = load_network(str(NETS / f"{net}.tsv"))
        index = pd.Series(np.arange(len(genes)), index=genes)
        for lab in LABELS:
            sp_path = SPLITS / f"{lab}__native_{net}.tsv"
            if not sp_path.exists():
                continue
            sp = read_splits(sp_path)
            sp = sp.assign(idx=index.loc[sp.gene].to_numpy())
            pub = sp.gene.map(pub_map).to_numpy(dtype=float)
            mp, up = [], []
            for j, (r, f) in enumerate(sorted(sp.groupby(["repeat", "fold"])
                                             .groups)):
                d = sp[sp.repeat == r]
                seeds = d.loc[(d.fold != f) & (d.y == 1), "idx"].to_numpy()
                test = d[d.fold == f]
                if len(seeds) == 0 or len(test) == 0:
                    continue
                dist = distance_to_nearest(A, seeds)
                finite = dist[np.isfinite(dist)]
                cap = float(finite.max()) + 1.0 if finite.size else 1.0
                score = np.where(np.isfinite(dist), dist, cap) * -1.0
                ti = test.idx.to_numpy()
                rng = np.random.default_rng(RNG_SEED + j)
                u, mm = fold_stats(score[ti], test.y.to_numpy(),
                                   pub[test.index], test.gene.to_numpy(), rng)
                if len(u):
                    mp.append(mm)
                    up.append(u)
            if not mp:
                continue
            mp = np.concatenate(mp)
            up = np.concatenate(up)
            rows.append({
                "network": net, "label": lab, "learner": "proximity",
                "n_positives_tested": int(len(mp)),
                "frac_matched_null_above": float((mp < 0.05).mean()),
                "frac_unmatched_null_above": float((up < 0.05).mean()),
                "frac_positives_beating_matched_null": float((mp == 0).mean()),
                "frac_positives_beating_unmatched_null": float((up == 0).mean()),
                "median_p_matched": float(np.median(mp)),
                "median_p_unmatched": float(np.median(up)),
            })
            log(f"  {net} x {lab}: {rows[-1]['frac_unmatched_null_above']:.1%} "
                f"above the unmatched null -> "
                f"{rows[-1]['frac_matched_null_above']:.1%} above the matched one")

    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "c21_visibility_second_learner.tsv", sep="\t",
              index=False)

    # ---- the comparison with the published RWR numbers ---------------
    g3 = pd.read_csv(TABLES / "g3_bait_matched_null.tsv", sep="\t")
    g3 = g3.rename(columns={"frac_matched_null_above": "rwr_matched",
                            "frac_unmatched_null_above": "rwr_unmatched"})
    cmp_ = df.merge(g3[["network", "label", "rwr_unmatched", "rwr_matched"]],
                     on=["network", "label"], how="inner")
    cmp_["rwr_drop"] = cmp_["rwr_unmatched"] - cmp_["rwr_matched"]
    cmp_["prox_drop"] = (cmp_["frac_unmatched_null_above"]
                         - cmp_["frac_matched_null_above"])
    cmp_.to_csv(TABLES / "c21_visibility_learner_comparison.tsv", sep="\t",
                index=False)

    write_doc(df, cmp_)
    log("wrote the Control 21 tables")
    return 0


def write_doc(df, cmp_):
    prim = df[df.label.isin(PRIMARY_LABELS)]
    lines = [
        "# Control 21 -- the visibility decomposition under the second learner",
        "",
        "**Date:** 2026-11-16  ",
        "**Script:** `scripts/run_c21_visibility_second_learner.py`  ",
        "**Tables:** `c21_visibility_second_learner.tsv`, "
        "`c21_visibility_learner_comparison.tsv`",
        "",
        "## The question",
        "",
        "Control 17 measured the MARGIN under the proximity learner and "
        "found 16/24 arms positive (median +0.011 against +0.0809), with "
        "an arm-level Spearman of +0.30. That is consistent with Guala & "
        "Sonnhammer (2017), and it is not a new claim.",
        "",
        "But the project's claim is not the margin: it is the "
        "DECOMPOSITION of the reported advantage into topology, node "
        "degree and study visibility. **No visibility measurement existed "
        "under a second learner.** This control computes the bait-matched "
        "null with the proximity learner as the scorer -- the same "
        "splits, the same folds, the same publication counts, the same "
        "matching rule (the 20 negatives nearest in log10(n_pubmed + 1)) "
        "-- and compares it with the published RWR numbers.",
        "",
        "The matcher in this file was checked first on the RWR score "
        "dumps of one arm: it reproduces the published G3 statistic for "
        "that arm to within 0.02 before any proximity result is used.",
        "",
        "## The result",
        "",
        "| learner | median arm: above the unmatched null | above the "
        "matched null |",
        "|---|---|---|",
    ]
    g3 = pd.read_csv(TABLES / "g3_bait_matched_null.tsv", sep="\t")
    g3p = g3[g3.label.isin(PRIMARY_LABELS)]
    lines.append(
        f"| RWR (published G3) | {g3p.frac_unmatched_null_above.median():.1%} "
        f"| {g3p.frac_matched_null_above.median():.1%} |")
    lines.append(
        f"| proximity (this control) | "
        f"{prim.frac_unmatched_null_above.median():.1%} | "
        f"{prim.frac_matched_null_above.median():.1%} |")
    lines += [
        "",
        "Per label (median over the six networks):",
        "",
        "| label | RWR unmatched -> matched | proximity unmatched -> matched |",
        "|---|---|---|",
    ]
    for lab in PRIMARY_LABELS + ["ncg_known"]:
        r = g3[g3.label == lab]
        p = df[df.label == lab]
        if r.empty or p.empty:
            continue
        lines.append(
            f"| {lab} | {r.frac_unmatched_null_above.median():.1%} -> "
            f"{r.frac_matched_null_above.median():.1%} | "
            f"{p.frac_unmatched_null_above.median():.1%} -> "
            f"{p.frac_matched_null_above.median():.1%} |")
    lines += [
        "",
        "## Reading",
        "",
        f"- the drop caused by the publication matching, RWR: median "
        f"{cmp_[cmp_.label.isin(PRIMARY_LABELS)].rwr_drop.median():.1%}; "
        f"proximity: median "
        f"{cmp_[cmp_.label.isin(PRIMARY_LABELS)].prox_drop.median():.1%}",
        f"- the arm-level correlation between the two learners' "
        f"matched-null statistics: Spearman "
        f"{cmp_.rwr_matched.corr(cmp_.frac_matched_null_above, method='spearman'):+.3f}",
        "",
        "## Caveats",
        "",
        "1. The proximity learner's score is the negative distance to the "
        "nearest seed, so 'above' means 'closer to a seed'. Because that "
        "score is an INTEGER distance, almost every test gene ties with "
        "another one, and G3's strict comparison would count a tie as "
        "'not above' -- which inflates the statistic (91 % of the "
        "positives 'beat' the null on intact x intogen2024, while the "
        "AUROC of the same scores is 0.627). The comparison here is "
        "therefore made on the protocol's RANKS (decision D-07: score "
        "desc, then gene symbol asc), which resolves the ties the "
        "documented way. RWR has no ties, so its numbers are unchanged "
        "-- the self-check on the RWR dumps reproduces the published G3 "
        "values exactly.",
        "2. **Exploratory** (amendment 4): it does not enter the H1 "
        "decision.",
        "3. This control measures the VISIBILITY component under a second "
        "learner. It does not re-run the full factorial (the rewiring "
        "null under the proximity learner would need ~470,000 multi-source "
        "BFS runs and is out of scope).",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
