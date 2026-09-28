"""W3 Control 5 — Cross-catalogue label transfer, single-shot.

Loads the network and the transfer split (produced by
``scripts/build_transfer_splits.py``). The transfer split has:

    gene, fold, y, label_id, universe_id, transfer_id, role
        role == "target"    (label_b positives ∩ LCC)
        role == "background" (LCC \\ (label_a ∪ label_b))
        fold ∈ [0..n_folds-1]  — for transfer evaluation we do NOT
        need to split the target into folds; the target set is the
        target set. We therefore aggregate across folds.

Seeds = label_a positives ∩ LCC \\ (label_b positives ∩ LCC)
       — i.e. we use the asymmetric label_a \\ label_b as the seed set,
       to avoid leaking label_b into the training signal.

Output: a single arm result with N bootstrap folds (random background
draws) for RWR and degree baselines. C5 schema matches A's ``run_arm``
(network, label_set, universe, alpha, repeat, fold, method, n_iter,
auroc, ...). The single fixed-target set is recorded once per fold;
statistics (mean, SD, paired t-test) are computed by
``scripts/results_stats_transfer.py``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from d1.engine.rwr import ALPHA, rwr, seed_matrix, transition_T  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

LABELS_DIR = REPO_ROOT / "data" / "processed" / "labels"


def _auroc_subset(scores: np.ndarray, target_mask: np.ndarray,
                  bg_indices: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score
    idx = np.concatenate([np.flatnonzero(target_mask), bg_indices])
    y = np.zeros(len(scores), dtype=int)
    y[target_mask] = 1
    return float(roc_auc_score(y[idx], scores[idx]))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--network", required=True, type=Path)
    ap.add_argument("--transfer-split", required=True, type=Path)
    ap.add_argument("--seeds-label", required=True, type=str,
                    help="Label id whose positives are the training seeds.")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--n-bootstrap", type=int, default=100,
                    help="Number of background draws to compute per-fold SD.")
    ap.add_argument("--background-size", type=int, default=1000)
    ap.add_argument("--alpha", type=float, default=ALPHA)
    ap.add_argument("--seed", type=int, default=20260926)
    args = ap.parse_args(argv)

    transfer_df = pd.read_csv(args.transfer_split, sep="\t", dtype=str)
    seeds_df = pd.read_csv(LABELS_DIR / f"{args.seeds_label}.tsv",
                           sep="\t", dtype=str)[["gene"]].dropna()

    # load network
    genes, A = load_network(args.network)
    index = pd.Series(np.arange(len(genes)), index=genes)

    target_genes = sorted(set(transfer_df[transfer_df.role == "target"].gene))
    background_genes = sorted(set(transfer_df[transfer_df.role == "background"].gene))
    # CRITICAL: use the NETWORK's full LCC for the seed pool, not the
    # transfer file's sampled LCC. The transfer file samples only 1000
    # background genes for storage reasons, but the seed pool must be
    # over the whole LCC.
    network_lcc = set(genes)
    seed_pool = sorted((set(seeds_df.gene) & network_lcc) - set(target_genes))
    if not seed_pool:
        raise ValueError(f"empty seed pool for {args.seeds_label} on "
                         f"{args.network.name}; check LCC overlap "
                         f"(label ∩ LCC = {len(set(seeds_df.gene) & network_lcc)})")

    target_idx = index.loc[target_genes].to_numpy()
    bg_idx_pool = index.loc[background_genes].to_numpy()
    seed_idx = index.loc[seed_pool].to_numpy()
    label_b = transfer_df.label_id.iloc[0]
    print(f"  transfer arm {args.seeds_label} -> {label_b} on {args.network.stem}")
    print(f"  n_target={len(target_idx)}  n_background={len(bg_idx_pool)}  "
          f"n_seeds={len(seed_idx)}  alpha={args.alpha}")

    # RWR
    WT = transition_T(A)
    P0 = seed_matrix(len(genes), [seed_idx])
    P, n_iter = rwr(WT, P0, alpha=args.alpha)
    rwr_scores = np.asarray(P[:, 0]).ravel()

    # degree baseline
    deg = np.asarray(A.sum(axis=1)).ravel()

    # bootstrap
    rng = np.random.default_rng(args.seed)
    target_mask = np.zeros(len(genes), dtype=bool)
    target_mask[target_idx] = True
    bg_pool_size = len(bg_idx_pool)
    if bg_pool_size < args.background_size:
        bg_size = bg_pool_size
        print(f"  [warn] bg pool {bg_pool_size} < requested {args.background_size}, "
              f"using {bg_size}")
    else:
        bg_size = args.background_size

    rows = []
    base = {
        "network": args.network.stem,
        "label_set": f"{args.seeds_label}->{label_b}",
        "universe": transfer_df.universe_id.iloc[0],
        "alpha": args.alpha,
        "n_iter": n_iter,
    }
    for fold in range(args.n_bootstrap):
        sel = rng.choice(bg_idx_pool, size=bg_size, replace=False)
        auroc_rwr = _auroc_subset(rwr_scores, target_mask, sel)
        auroc_deg = _auroc_subset(deg, target_mask, sel)
        rows.append({**base, "repeat": 0, "fold": fold, "method": "rwr",
                     "auroc": auroc_rwr,
                     "n_test_pos": int(target_mask.sum()),
                     "n_test": int(target_mask.sum() + bg_size)})
        rows.append({**base, "repeat": 0, "fold": fold, "method": "degree",
                     "auroc": auroc_deg,
                     "n_test_pos": int(target_mask.sum()),
                     "n_test": int(target_mask.sum() + bg_size)})

    # full-set AUROC (single deterministic eval against the whole bg)
    from sklearn.metrics import roc_auc_score
    y_full = np.zeros(len(genes), dtype=int)
    y_full[target_mask] = 1
    auroc_rwr_full = float(roc_auc_score(y_full, rwr_scores))
    auroc_deg_full = float(roc_auc_score(y_full, deg))
    rows.append({**base, "repeat": 0, "fold": -1, "method": "rwr",
                 "auroc": auroc_rwr_full,
                 "n_test_pos": int(target_mask.sum()),
                 "n_test": int(len(genes))})
    rows.append({**base, "repeat": 0, "fold": -1, "method": "degree",
                 "auroc": auroc_deg_full,
                 "n_test_pos": int(target_mask.sum()),
                 "n_test": int(len(genes))})

    out_df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {args.out} ({len(out_df)} rows)")

    # quick bootstrap summary
    boot = out_df[out_df.fold >= 0]
    pvt = boot.pivot_table(index="method", values="auroc",
                           aggfunc=["mean", "std", "count"])
    print("\nBootstrap summary:")
    print(pvt.to_string())
    delta = (boot[boot.method == "rwr"].auroc.values
             - boot[boot.method == "degree"].auroc.values)
    from scipy import stats
    t = stats.ttest_1samp(delta, popmean=0, alternative="greater")
    print(f"\nmean ΔAUROC = {delta.mean():.4f} ± {delta.std():.4f}  "
          f"t = {t.statistic:.2f}  p = {t.pvalue:.3g}")
    print(f"\nFull-set AUROC: rwr = {auroc_rwr_full:.4f}, "
          f"degree = {auroc_deg_full:.4f}, "
          f"delta = {auroc_rwr_full - auroc_deg_full:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())