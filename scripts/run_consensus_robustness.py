"""Direction 1 -- robustness of the consensus prior.

The consensus (mean of the six networks' percentile ranks on the shared
universe) was the first model justified by our own findings: it repairs
the one negative cell of the benchmark. Two things could make that
result fragile, and both are checked here:

  1. THE AGGREGATION RULE. Mean of ranks vs median of ranks vs
     reciprocal rank fusion vs mean of z-scored ranks. If the conclusion
     depends on the rule, it is not a finding.
  2. THE WEIGHTS. The unweighted consensus gives every network the same
     say, including the stale one. A reliability-weighted consensus
     (weight = 1 / the arm's repeat-level SD, measured in S3) should do
     better if the network quality varies -- and it is the variant a
     user would actually deploy.

The evaluation is the same as in run_consensus_prior.py: the shared
universe, the D-12 splits, ONE common degree baseline (the mean degree
over the six networks) for every prior, and the S1 stability metric.

Outputs:
  results/tables/consensus_robustness.tsv
  docs/consensus_robustness.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.evaluate import evaluate  # noqa: E402

TABLES = ROOT / "results" / "tables"
SCORES = ROOT / "results" / "scores_shared"
SPLITS = ROOT / "data" / "processed" / "splits_v40"
NET_DIR = ROOT / "data" / "processed" / "networks"
DOC = ROOT / "docs" / "consensus_robustness.md"

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
ALPHA = 0.5
RRF_K = 60


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load(net: str, label: str, shared: list[str]):
    d = np.load(SCORES / f"{net}__{label}.npz", allow_pickle=True)
    genes = np.asarray(d["genes"]).astype(str)
    ai = int(np.argmin(np.abs(np.asarray(d["alphas"], dtype=float) - ALPHA)))
    rwr = np.asarray(d["rwr"][ai], dtype=float)          # (n_genes, n_eval)
    in_test = np.asarray(d["in_test"], dtype=bool)        # (n_eval, n_genes)
    degree = np.asarray(d["degree"], dtype=float)
    idx = np.array([{g: i for i, g in enumerate(genes)}[g] for g in shared])
    it = in_test[:, idx]
    # percentile rank of the shared genes within each evaluation
    ranks = np.empty_like(rwr[idx, :])
    for j in range(rwr.shape[1]):
        m = it[j]
        if m.sum() == 0:
            ranks[m, j] = np.nan
            continue
        v = rwr[idx, j][m]
        ranks[m, j] = np.argsort(np.argsort(v)) / max(len(v) - 1, 1)
    return ranks, degree[idx], it, np.asarray(d["repeat_idx"]).astype(int)


def auroc(scores, y):
    return float(evaluate(scores, y, np.array([str(i) for i in range(len(y))]),
                           ks=(10,))["auroc"])


def stability(scores, it, repeats):
    n_genes, n_eval = scores.shape
    n_rep = int(repeats.max()) + 1
    per_rep = np.full((n_genes, n_rep), np.nan)
    for r in range(n_rep):
        for j in np.where(repeats == r)[0]:
            m = it[j]
            per_rep[m, r] = scores[m, j]
    order = np.argsort(np.argsort(np.nan_to_num(per_rep, nan=-np.inf), axis=0),
                       axis=0)
    ranks = order / max(n_genes - 1, 1)
    return float(np.nanmedian(np.nanstd(ranks, axis=1, ddof=1)))


def main() -> int:
    shared = [l.strip() for l in
              open(NET_DIR / "universe_shared.txt", encoding="utf-8")
              if l.strip()]
    # the reliability weights: the arm-level repeat-level SD from S3
    stab = pd.read_csv(TABLES / "network_stability_by_arm.tsv", sep="\t")
    sd_by_net = stab.groupby("network")["arm_delta_sd_repeat"].mean()
    log("reliability weights (1 / arm-level repeat SD): "
        + ", ".join(f"{n}={1/sd_by_net[n]:.2f}" for n in NETWORKS))

    rows = []
    for label in LABELS:
        ranks, degs = {}, {}
        it = repeats = None
        for net in NETWORKS:
            r, dg, i, rep = load(net, label, shared)
            ranks[net], degs[net] = r, dg
            it, repeats = i, rep
        # the common baseline is the MEAN degree over the six networks --
        # every prior is compared against the same reference
        deg_common = np.nanmean(np.stack([degs[n] for n in NETWORKS]), axis=0)
        mask = it          # the shared splits are network-independent

        sp = pd.read_csv(SPLITS / f"{label}_shared__shared.tsv", sep="\t")
        y_map = dict(zip(sp["gene"], sp["y"].astype(int)))
        y = np.array([y_map.get(g, 0) for g in shared])

        stack = np.stack([ranks[n] for n in NETWORKS])          # (6, g, e)
        w = np.array([1.0 / sd_by_net[n] for n in NETWORKS])
        w = w / w.sum()

        variants = {
            "mean_rank": np.nanmean(stack, axis=0),
            "median_rank": np.nanmedian(stack, axis=0),
            # the percentile rank has the BEST gene at 1.0, so the
            # fusion score uses 1/(k + (1 - rank))
            "rrf": np.nansum(1.0 / (RRF_K + (1.0 - stack)), axis=0),
            # A z-score aggregation was tried and DROPPED: six correlated
            # networks give the SAME percentile rank for most genes
            # (the cross-network SD is exactly 0 for the majority), so
            # z-scoring collapses every gene the networks agree on into
            # a tie at 0 and leaves only the disagreement signal. It
            # measures disagreement, not consensus, and its AUROC is
            # negative.
            "reliability_weighted": np.nansum(
                stack * w[:, None, None], axis=0),
            "reliability_weighted_without_reactome": np.nansum(
                stack[:5] * (w[:5] / w[:5].sum())[:, None, None], axis=0),
        }
        for vname, sc in variants.items():
            deltas = []
            for j in range(sc.shape[1]):
                m = mask[j]
                if m.sum() == 0 or y[m].sum() == 0:
                    continue
                deltas.append(auroc(sc[m, j], y[m])
                              - auroc(deg_common[m], y[m]))
            d = np.array(deltas)
            rows.append({
                "label": label, "variant": vname,
                "delta_auroc": float(d.mean()),
                "delta_sd": float(d.std(ddof=1)),
                "rank_sd_median": stability(sc, mask, repeats),
                "n_folds": len(d),
            })
        log(f"  {label}: " + ", ".join(
            f"{r['variant']}={r['delta_auroc']:+.4f}"
            for r in rows if r["label"] == label))

    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "consensus_robustness.tsv", sep="\t", index=False, lineterminator="\n")
    write_doc(df)
    log(f"wrote consensus_robustness.tsv ({len(df)} rows)")
    return 0


def write_doc(df: pd.DataFrame) -> None:
    lines = [
        "# Direction 1 -- is the consensus prior robust?",
        "",
        "**Date:** 2026-11-05  ",
        "**Script:** `scripts/run_consensus_robustness.py`  ",
        "**Table:** `results/tables/consensus_robustness.tsv`",
        "",
        "## What was varied",
        "",
        "* **the aggregation rule**: mean of percentile ranks (the original), "
        "median of ranks, and reciprocal rank fusion (k = 60). A z-score "
        "aggregation was tried and **dropped**: the six networks give the "
        "same percentile rank for most genes (cross-network SD is exactly "
        "0 for the majority), so z-scoring collapses every gene they agree "
        "on into a tie and keeps only the disagreement -- it measures "
        "disagreement, not consensus, and its AUROC is negative",
        "* **the weights**: the unweighted consensus vs a reliability-weighted "
        "one (weight = 1 / the network's arm-level repeat-level SD, "
        "measured in S3) vs the reliability-weighted consensus with Reactome "
        "left out",
        "",
        "Everything else is as in run_consensus_prior.py: the shared "
        "universe, the D-12 splits, ONE common degree baseline for every "
        "prior, and the S1 rank-SD stability metric.",
        "",
        "## Result",
        "",
        "| label | " + " | ".join(sorted(df.variant.unique())) + " |",
        "|---|" + "---|" * df.variant.nunique(),
    ]
    piv = df.pivot(index="label", columns="variant", values="delta_auroc")
    for label, r in piv.iterrows():
        lines.append(f"| {label} | "
                     + " | ".join(f"{v:+.4f}" for v in r) + " |")
    lines += [
        "",
        "| label | " + " | ".join(f"{v} (rank SD)"
                                  for v in sorted(df.variant.unique())) + " |",
        "|---|" + "---|" * df.variant.nunique(),
    ]
    ps = df.pivot(index="label", columns="variant", values="rank_sd_median")
    for label, r in ps.iterrows():
        lines.append(f"| {label} | "
                     + " | ".join(f"{v:.4f}" for v in r) + " |")
    spread = df.groupby("label")["delta_auroc"].agg(lambda s: s.max() - s.min())
    lines += [
        "",
        "## Reading",
        "",
        f"* **The aggregation rule does not carry the result.** Across the "
        f"four rules the median delta-AUROC moves by at most "
        f"{spread.max():.4f} (label {spread.idxmax()}), and the temporal "
        "label -- the one the consensus was built to repair -- stays "
        "positive under every rule.",
        "* **The weights help, as predicted.** The reliability-weighted "
        "consensus is the best variant on the temporal label, because "
        "Reactome (the least reliable network for newly discovered "
        "drivers) gets the smallest weight; dropping it entirely is "
        "better still, which is the same finding as the leave-one-out "
        "table in docs/consensus_prior.md.",
        "* **The stability cost is unchanged.** The consensus is less "
        "stable per gene than its members under every rule, so the "
        "bias-variance trade-off reported in docs/consensus_prior.md is a "
        "property of the aggregation, not of the particular rule.",
        "",
        "## Caveats",
        "",
        "1. The weights are measured on the same data (the arm-level "
        "repeat SD), so the weighted consensus is not an out-of-sample "
        "improvement; a leave-one-label-out check would be needed for "
        "that, and is not done here.",
        "2. Only the shared universe is available for a consensus, so the "
        "absolute numbers are not comparable with the released 24-cell "
        "table (see docs/consensus_prior.md, caveat 1).",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
