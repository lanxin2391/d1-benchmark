"""D-09 temporal drift control.

Per B's proposal (locked at S2):

    D-09 temporal split label = intogen_temporal_new (152 genes):
    IntOGen drivers added in the 2020 -> 2024 window (present in 2024 release
    but absent from 2020 release). Use this set as the *training-set temporal
    holdout*: train ranking on intogen2024 \\ temporal_new, evaluate
    AUROC uplift on temporal_new alone.

Algorithm
---------
1. Load all 6 primary networks.
2. Load intogen2024 labels (633 positives) and intogen_temporal_new
   labels (152 positives). Training seeds = intogen2024 - temporal_new.
3. For each repeat (we do n_rep = 50, more than the 10 the main
   pipeline uses, to stabilise the AUROC estimate for this small label):
      a. Seed matrix over the network's LCC, using only training seeds
         that map to nodes in the network.
      b. RWR(alpha) and degree scores over all genes.
      c. AUROC on (positive = temporal_new, negative = all other genes
         in the network's LCC that are NOT in temporal_new).

Output:
    * data/processed/labels/intogen2024_minus_temporal.tsv   (training seeds)
    * data/processed/labels/intogen_temporal_new.tsv        (B already created this)
    * results/runs/d09_<network>.tsv                       (per-network C5)
    * results/tables/d09_table.tsv                         (Table for S2)
"""
import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.evaluate import evaluate  # noqa: E402
from d1.engine.rwr import transition_T, seed_matrix, rwr  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETWORKS = ["funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome"]
LABELS_DIR = ROOT / "data" / "processed" / "labels"
PROCESSED_NET = ROOT / "data" / "processed" / "networks"
RUNS_DIR = ROOT / "results" / "runs"
TABLE_DIR = ROOT / "results" / "tables"
RUNS_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

N_REP = 50                                    # more repeats for the small label


def _label_genes(label_id: str) -> set[str]:
    return set(pd.read_csv(LABELS_DIR / f"{label_id}.tsv", sep="\t",
                            dtype=str).gene)


def main():
    intogen = _label_genes("intogen2024")
    temporal_new = _label_genes("intogen_temporal_new")
    train_seeds_global = intogen - temporal_new        # 481 genes
    test_pos = temporal_new                          # 152 genes

    print(f"D-09 temporal drift: {len(intogen)} intogen positives, "
          f"{len(test_pos)} temporal_new positives, "
          f"{len(train_seeds_global)} training seeds "
          f"(intogen2024 minus temporal_new)")

    # Write training-seeds set as a label file for posterity
    pd.DataFrame(
        {"gene": sorted(train_seeds_global), "rank": range(1, len(train_seeds_global) + 1),
         "source": "intogen2024", "release": "minus_temporal_new"}
    ).to_csv(LABELS_DIR / "intogen2024_minus_temporal.tsv", sep="\t", index=False)
    print(f"  wrote {LABELS_DIR / 'intogen2024_minus_temporal.tsv'}")

    summary_rows = []
    for net_id in NETWORKS:
        net_path = PROCESSED_NET / f"{net_id}.tsv"
        print(f"\n[{net_id}] loading network from {net_path.name}")
        genes, A = load_network(str(net_path))
        WT = transition_T(A)
        gene_to_idx = pd.Series(np.arange(len(genes)), index=genes)

        # Map seeds to network indices
        train_in_net = gene_to_idx.reindex(sorted(train_seeds_global)).dropna()
        test_in_net = gene_to_idx.reindex(sorted(test_pos)).dropna()
        if len(test_in_net) == 0:
            print(f"  SKIP: no temporal_new positives in network")
            continue
        train_idx = train_in_net.to_numpy(dtype=int)
        test_idx = test_in_net.to_numpy(dtype=int)
        # Test labels (1 if gene in temporal_new, else 0)
        y_full = np.zeros(len(genes), dtype=int)
        y_full[test_idx] = 1
        # Degree baseline
        deg = A.sum(axis=1).A.ravel().astype(float)

        # Test set = ALL genes in the network (positives = temporal_new,
        # negatives = everyone else). AUROC requires both classes.
        neg_pool = np.setdiff1d(np.arange(len(genes)), test_idx)
        test_all_idx = np.concatenate([test_idx, neg_pool])
        test_y_all = np.zeros(len(genes), dtype=int)
        test_y_all[test_idx] = 1
        test_y = test_y_all[test_all_idx]
        test_g = [genes[i] for i in test_all_idx]

        # For each repeat, build a fresh seed matrix.
        rwr_rows = []
        deg_rows = []
        for rep in range(N_REP):
            rng = np.random.default_rng(rep)
            seed_idx_shuf = train_idx.copy()
            rng.shuffle(seed_idx_shuf)
            P0 = seed_matrix(len(genes), [seed_idx_shuf.tolist()])
            test_d = deg[test_all_idx]
            deg_metrics = evaluate(test_d, test_y, test_g)
            deg_rows.append({"network": net_id, "repeat": rep, "method": "degree",
                             "n_test_pos": int(test_y.sum()),
                             "n_test": int(len(test_y)),
                             **deg_metrics})
            for alpha in (0.3, 0.5, 0.7):
                P, _ = rwr(WT, P0, alpha=alpha)
                test_rwr = P[test_all_idx, 0]
                rwr_metrics = evaluate(test_rwr, test_y, test_g)
                rwr_rows.append({"network": net_id, "repeat": rep, "method": "rwr",
                                 "alpha": alpha,
                                 "n_test_pos": int(test_y.sum()),
                                 "n_test": int(len(test_y)),
                                 **rwr_metrics})
        deg_df = pd.DataFrame(deg_rows)
        rwr_df = pd.DataFrame(rwr_rows)
        # Join RWR rows to degree rows by (network, repeat) to compute delta
        joined = rwr_df.merge(deg_df[["network", "repeat", "auroc"]]
                              .rename(columns={"auroc": "degree_auroc"}),
                              on=["network", "repeat"])
        joined["delta_auroc"] = joined["auroc"] - joined["degree_auroc"]
        # Per-arm C5 file
        out_path = RUNS_DIR / f"d09_{net_id}.tsv"
        joined.to_csv(out_path, sep="\t", index=False)

        # Summary row
        for alpha in (0.3, 0.5, 0.7):
            sub = joined[joined.alpha == alpha]
            summary_rows.append(
                {"network": net_id, "alpha": alpha,
                 "n_train_seeds": len(train_idx),
                 "n_test_pos": int(sub.n_test_pos.iloc[0]) if len(sub) else 0,
                 "n_test": int(sub.n_test.iloc[0]) if len(sub) else 0,
                 "rwr_auroc_mean": sub.auroc.mean(),
                 "degree_auroc_mean": sub.degree_auroc.mean(),
                 "delta_auroc_mean": sub.delta_auroc.mean(),
                 "delta_auroc_sd": sub.delta_auroc.std(),
                 "n_rep": len(sub)}
            )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TABLE_DIR / "d09_table.tsv", sep="\t", index=False)
    print("\nD-09 Table:")
    pivot = summary.pivot_table(
        index="network", columns="alpha", values="delta_auroc_mean"
    ).round(4)
    print(pivot.to_string())


if __name__ == "__main__":
    main()