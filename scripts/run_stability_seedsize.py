"""Does the ranking instability of a network prior scale with the seed set?

The stability-regularised scorer failed: shrinking a gene's score toward the
median of its score neighbourhood is order-preserving, so it cannot change the
per-gene rank SD, and any operator that *does* change the ordering injects new
noise and makes the rank SD worse.

The diagnosis says the instability is driven by each gene's **own** score
fluctuating across repeats, and that the score is a *linear* functional of the
seed set (RWR: ``P = alpha (I - (1-alpha) W^T)^-1 P0``). Two consequences:

1. bagging over seed resamples is provably a **no-op** for RWR (the estimator is
   linear in the seed weights, so the bootstrap mean of ``P0`` is ``P0``);
2. what *does* control the instability is how many seeds the propagation
   averages over. The ranking of a gene is the ranking of a mean over seeds, and
   the sampling noise of that mean falls like ``1 / sqrt(n_seeds)``.

This script measures the consequence: the same arm, evaluated with a training
seed set truncated to a fraction ``f`` of the available positives. It reports the
per-gene rank SD and the AUROC at every ``f`` -- i.e. the stability /
discrimination trade-off of the *evaluation design*, which is the only axis on
which a prior's stability can actually be bought.

Usage:
    python scripts/run_stability_seedsize.py --only-net funmap --only-label intogen2024
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from d1.engine.evaluate import evaluate                       # noqa: E402
from d1.engine.rwr import rwr, seed_matrix, transition_T      # noqa: E402
from d1.engine.stability import percentile_ranks              # noqa: E402
from d1.networks.io import load_network                      # noqa: E402

NET_DIR = REPO_ROOT / "data" / "processed" / "networks"
SPLIT_DIR = REPO_ROOT / "data" / "processed" / "splits_v40"
NETWORKS = ("funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome")
LABELS = ("intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new")
FRACTIONS = (0.05, 0.10, 0.25, 0.50, 1.00)
SEED = 20260926


def run_arm(net_id, label_id, fractions=FRACTIONS, alpha=0.5):
    net_genes, A = load_network(str(NET_DIR / f"{net_id}.tsv"))
    sp = pd.read_csv(SPLIT_DIR / f"{label_id}__native_{net_id}.tsv", sep="\t",
                     dtype={"gene": str})
    gpos = pd.Series(np.arange(len(net_genes)), index=net_genes)
    sp = sp.assign(idx=gpos.loc[sp.gene].to_numpy()).sort_values(["repeat", "fold"])
    cols = list(sp.groupby(["repeat", "fold"]).groups.keys())
    WT = transition_T(A)
    n_genes = len(net_genes)

    seeds, tests = [], []
    for (r, f) in cols:
        block = sp[sp.repeat == r]
        seeds.append(block.loc[(block.fold != f) & (block.y == 1), "idx"].to_numpy())
        tests.append(block[block.fold == f])

    rows = []
    # the degree baseline does not depend on the seed set: read it from the dump
    deg_path = (REPO_ROOT / "results" / "scores_v40" /
                f"{net_id}__{label_id}__alpha{alpha}.npz")
    deg_auc = np.nan
    if deg_path.exists():
        dd = np.load(deg_path, allow_pickle=True)
        deg_auc = float(np.mean(dd["auroc_degree"]))
    for frac in fractions:
        rng = np.random.default_rng(SEED + int(frac * 1000))
        sub_seeds = []
        for idx in seeds:
            if frac >= 1.0 or idx.size == 0:
                sub_seeds.append(idx)
                continue
            m = max(1, int(round(frac * idx.size)))
            pick = rng.choice(idx.size, size=m, replace=False)
            sub_seeds.append(idx[np.sort(pick)])
        P0 = seed_matrix(n_genes, sub_seeds)
        P, n_iter = rwr(WT, P0, alpha=alpha)

        pct = np.full((len(cols), n_genes), np.nan)
        auc, nauc, npos = [], [], []
        for j, t in enumerate(tests):
            ti = t.idx.to_numpy()
            m = evaluate(P[ti, j], t.y.to_numpy(), t.gene.to_numpy())
            auc.append(m["auroc"])
            nauc.append(m["n_test"])
            npos.append(m["n_test_pos"])
            pct[j, ti] = percentile_ranks(P[ti, j], t.gene.to_numpy())
        sd = np.nanstd(pct, axis=0, ddof=1)
        rows.append({
            "network": net_id, "label_set": label_id, "alpha": alpha,
            "frac": frac,
            "mean_n_seeds": float(np.mean([len(s) for s in sub_seeds])),
            "median_rank_sd": float(np.nanmedian(sd)),
            "mean_rank_sd": float(np.nanmean(sd)),
            "auroc": float(np.mean(auc)),
            "auroc_sd": float(np.std(auc, ddof=1)),
            "auroc_degree": deg_auc,
            "delta_auroc": float(np.mean(auc)) - deg_auc if np.isfinite(deg_auc) else np.nan,
            "n_test": int(np.mean(nauc)),
            "n_test_pos": float(np.mean(npos)),
        })
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    ap.add_argument("--out", default="")
    args = ap.parse_args(argv)
    nets = args.only_net or list(NETWORKS)
    labels = args.only_label or list(LABELS)
    out = []
    for label in labels:
        for net in nets:
            t0 = time.time()
            df = run_arm(net, label, alpha=args.alpha)
            out.append(df)
            print(f"[{net} x {label}] {time.time() - t0:.1f}s")
            print(df[["frac", "mean_n_seeds", "median_rank_sd", "auroc"]]
                  .round(4).to_string(index=False))
    res = pd.concat(out, ignore_index=True)
    path = Path(args.out) if args.out else (REPO_ROOT / "results" / "tables" /
                                            "stability_seed_size_scaling.tsv")
    path.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(path, sep="\t", index=False)
    print(f"\nWrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
