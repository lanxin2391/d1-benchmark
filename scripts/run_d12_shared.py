"""Shared-universe control (D-12).

Per the protocol:

    D-12 shared-universe control: RWR propagates on each network's full LCC;
    seeds and evaluation genes are restricted to the shared gene set;
    degree is the degree in the full LCC.

Algorithm
---------
1. Load all 6 primary networks. The "shared" universe = intersection
   of their gene sets (genes present in ALL six networks).
2. For each of the 4 labels, restrict positives to those in the
   shared universe.
3. Build splits using the shared universe (using build_splits.py).
4. For each (network, label) combination, run run_arm with the
   shared-universe splits. RWR propagates on the full network LCC
   (run_arm does this automatically); evaluation is restricted to
   the shared-universe genes.

Output:
    * data/processed/labels/<label>_shared.tsv         (shared-universe capped labels)
    * data/processed/splits/<label>_shared__native_<net>.tsv
    * results/runs/d12_<network>__<label>_shared.tsv
    * results/tables/d12_table.tsv
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.networks.io import load_network  # noqa: E402

LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
NETWORKS = ["funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome"]
LABELS_DIR = ROOT / "data" / "processed" / "labels"
PROCESSED_NET = ROOT / "data" / "processed" / "networks"
PROCESSED_SPLIT = ROOT / "data" / "processed" / "splits"
RUNS_DIR = ROOT / "results" / "runs"
TABLE_DIR = ROOT / "results" / "tables"
RUNS_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)


def main():
    # Step 1: shared universe
    print("Step 1: compute shared universe = intersection of 6 networks' gene sets")
    net_genes_list = []
    for net_id in NETWORKS:
        genes, _ = load_network(str(PROCESSED_NET / f"{net_id}.tsv"))
        net_genes_list.append(set(genes))
    shared = net_genes_list[0]
    for s in net_genes_list[1:]:
        shared &= s
    print(f"  shared universe: {len(shared):,} genes "
          f"(in all 6 networks' LCCs)")
    for s, net_id in zip(net_genes_list, NETWORKS):
        print(f"  {net_id}: {len(s):,} genes; shared fraction = {len(shared)/len(s):.1%}")
    # Save shared universe
    pd.Series(sorted(shared), name="gene").to_csv(
        PROCESSED_NET / "universe_shared.txt", index=False, header=False
    )

    # Step 2: write shared-universe-capped label files
    print("\nStep 2: cap labels to shared universe")
    for label in LABELS:
        df = pd.read_csv(LABELS_DIR / f"{label}.tsv", sep="\t", dtype=str)
        df_shared = df[df.gene.isin(shared)].copy()
        df_shared["source"] = f"{df_shared['source'].iloc[0]}_shared"
        df_shared["release"] = f"{df_shared['release'].iloc[0]}_shared"
        df_shared.to_csv(LABELS_DIR / f"{label}_shared.tsv", sep="\t", index=False)
        print(f"  {label}_shared: {len(df_shared)} positives")

    # Step 3: regenerate splits for shared universe
    print("\nStep 3: regenerate shared-universe splits")
    for net_id in NETWORKS:
        for label in LABELS:
            shared_label = f"{label}_shared"
            split_path = PROCESSED_SPLIT / f"{shared_label}__native_{net_id}.tsv"
            if split_path.exists():
                print(f"  [SKIP] {split_path.name}")
                continue
            print(f"  building {split_path.name}")
            subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "build_splits.py"),
                 "--label", shared_label,
                 "--universe", f"native_{net_id}"],
                check=True, cwd=str(ROOT))

    # Step 4: run arms at alpha=0.5
    print("\nStep 4: run 24 shared-universe arms at alpha=0.5")
    from d1.engine.evaluate import run_arm, margins
    summary_rows = []
    for net_id in NETWORKS:
        for label in LABELS:
            shared_label = f"{label}_shared"
            net_path = PROCESSED_NET / f"{net_id}.tsv"
            split_path = PROCESSED_SPLIT / f"{shared_label}__native_{net_id}.tsv"
            out_path = RUNS_DIR / f"d12_{net_id}__{shared_label}.tsv"
            res = run_arm(network=str(net_path), splits=str(split_path),
                          network_id=net_id, label_id=shared_label,
                          universe_id=f"native_{net_id}",
                          alphas=(0.5,))
            res.to_csv(out_path, sep="\t", index=False)
            m = margins(res)
            d = m[m.alpha == 0.5]
            n_test_pos = int(res[(res.alpha == 0.5) & (res.method == "rwr")].n_test_pos.iloc[0])
            n_test = int(res[(res.alpha == 0.5) & (res.method == "rwr")].n_test.iloc[0])
            if len(d):
                summary_rows.append({
                    "network": net_id, "label": shared_label, "alpha": 0.5,
                    "shared_genes": len(shared),
                    "n_test_pos": n_test_pos, "n_test": n_test,
                    "rwr_auroc": d.rwr_auroc.mean(),
                    "degree_auroc": d.degree_auroc.mean(),
                    "delta_auroc": d.delta_auroc.mean(),
                    "delta_auroc_sd": d.delta_auroc.std(),
                })
            print(f"  [{net_id} x {shared_label}] {len(res)} rows -> delta={d.delta_auroc.mean():.4f}")

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TABLE_DIR / "d12_table.tsv", sep="\t", index=False)
    print("\nD-12 Table (delta_AUROC at alpha=0.5):")
    pivot = summary.pivot_table(index="network", columns="label",
                                 values="delta_auroc").round(4)
    print(pivot.to_string())


if __name__ == "__main__":
    main()