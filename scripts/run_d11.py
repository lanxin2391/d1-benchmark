"""D-11 per-network common-size cap.

Per B's locked decision (S2):

    D-11 common size N = per-network min. For each (label x network) arm
    separately, N = min(positives_in_LCC) across the four labels for that
    network. Result (computed dynamically below):
        funmap N = N_f, string_full400 N = N_400, string_full700 N = N_700,
        string_full900 N = N_900, string_phys700 N = N_p,
        intact N = N_i, reactome N = N_r, rna_coexp N = N_rn.

Algorithm
---------
For each network:
    1. For each of the 4 labels, count positives in the network's LCC.
    2. N = min(counts) over the 4 labels.
    3. For each label, truncate the rank-1..N subset (keep the highest-ranked
       positives) so every label has exactly N positives.
    4. Hand the truncated label sets to scripts/build_splits.py to
       regenerate the C4 splits with the new positives.
    5. Run run_a4.py on the new (network, label_cap, native_<net>)
       combinations.

Output:
    * data/processed/labels/<label>_cap_<N>.tsv  (per-network capped labels)
    * data/processed/splits/<label>_cap_<N>__native_<net>.tsv
    * results/runs/d11_<network>__<label>.tsv   (per-arm C5)
    * results/tables/d11_table.tsv             (Table 1 for D-11)
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

NETWORKS = ["funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
LABELS_DIR = ROOT / "data" / "processed" / "labels"
PROCESSED_NET = ROOT / "data" / "processed" / "networks"
PROCESSED_SPLIT = ROOT / "data" / "processed" / "splits"
RUNS_DIR = ROOT / "results" / "runs"
TABLE_DIR = ROOT / "results" / "tables"
RUNS_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)
SCRATCH_DIR = ROOT / "scratch"
SCRATCH_DIR.mkdir(exist_ok=True)


def _label_genes(label_id: str) -> set[str]:
    df = pd.read_csv(LABELS_DIR / f"{label_id}.tsv", sep="\t", dtype=str)
    return set(df.gene)


def main():
    # Step 1: compute N per network from current label x network intersections.
    print("Step 1: compute per-network N (the minimum positives over labels)")
    counts = {}                                            # counts[net] = {label: int}
    for net_id in NETWORKS:
        genes, _ = load_network(str(PROCESSED_NET / f"{net_id}.tsv"))
        net_genes = set(genes)
        counts[net_id] = {}
        for label in LABELS:
            label_genes = _label_genes(label)
            counts[net_id][label] = len(label_genes & net_genes)
        n_min = min(counts[net_id].values())
        print(f"  {net_id}: " + ", ".join(f"{l}={c}" for l, c in counts[net_id].items())
              + f"  -> N = {n_min}")

    # Step 2: for each (network, label) pair, write a capped label file
    # that contains only the top-N (rank-1..N) positives.
    print("\nStep 2: write capped label files (top-N positives only)")
    caps = {}                                               # caps[(net, label)] = N
    for net_id in NETWORKS:
        n_min = min(counts[net_id].values())
        for label in LABELS:
            caps[(net_id, label)] = n_min
            df = pd.read_csv(LABELS_DIR / f"{label}.tsv", sep="\t", dtype=str)
            df = df.sort_values("rank", kind="stable")
            df_cap = df.head(n_min).copy()
            df_cap["source"] = f"{df_cap['source'].iloc[0]}_cap_{n_min}"
            df_cap["release"] = f"{df_cap['release'].iloc[0]}_cap"
            cap_label = f"{label}_cap"
            df_cap.to_csv(LABELS_DIR / f"{cap_label}.tsv", sep="\t", index=False)

    # Step 3: regenerate splits for every (capped_label, native_<net>) pair.
    print("\nStep 3: regenerate capped splits (one per (label, network))")
    for net_id in NETWORKS:
        for label in LABELS:
            cap_label = f"{label}_cap"
            split_path = PROCESSED_SPLIT / f"{cap_label}__native_{net_id}.tsv"
            if split_path.exists():
                print(f"  [SKIP] {split_path.name} already exists")
                continue
            print(f"  building {split_path.name}")
            subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "build_splits.py"),
                 "--label", cap_label,
                 "--universe", f"native_{net_id}"],
                check=True, cwd=str(ROOT))

    # Step 4: run each arm (6 networks x 4 capped labels) at alpha=0.5.
    print(f"\nStep 4: run 24 capped arms at alpha=0.5")
    summary_rows = []
    for net_id in NETWORKS:
        for label in LABELS:
            cap_label = f"{label}_cap"
            net_path = PROCESSED_NET / f"{net_id}.tsv"
            split_path = PROCESSED_SPLIT / f"{cap_label}__native_{net_id}.tsv"
            out_path = RUNS_DIR / f"d11_{net_id}__{cap_label}.tsv"
            t0 = pd.Timestamp.now()
            # Use a small inline helper to avoid re-running the whole A4.
            from d1.engine.evaluate import run_arm
            res = run_arm(network=str(net_path),
                          splits=str(split_path),
                          network_id=net_id,
                          label_id=cap_label,
                          universe_id=f"native_{net_id}",
                          alphas=(0.5,))
            res.to_csv(out_path, sep="\t", index=False)
            # Extract delta_AUROC for this arm (alpha=0.5 row)
            from d1.engine.evaluate import margins
            m = margins(res)
            d = m[m.alpha == 0.5]
            n_test_pos = int(res[(res.alpha == 0.5) & (res.method == "rwr")].n_test_pos.iloc[0])
            n_test = int(res[(res.alpha == 0.5) & (res.method == "rwr")].n_test.iloc[0])
            if len(d):
                summary_rows.append({
                    "network": net_id, "label": cap_label, "alpha": 0.5,
                    "N": caps[(net_id, label)],
                    "n_test_pos": n_test_pos,
                    "n_test": n_test,
                    "rwr_auroc": d.rwr_auroc.mean(),
                    "degree_auroc": d.degree_auroc.mean(),
                    "delta_auroc": d.delta_auroc.mean(),
                    "delta_auroc_sd": d.delta_auroc.std(),
                })
            print(f"  [{net_id} x {cap_label}] {len(res)} rows -> delta={d.delta_auroc.mean():.4f} +- {d.delta_auroc.std():.4f}")

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TABLE_DIR / "d11_table.tsv", sep="\t", index=False)
    print("\nD-11 Table (delta_AUROC at alpha=0.5):")
    pivot = summary.pivot_table(index="network", columns="label",
                                 values="delta_auroc").round(4)
    print(pivot.to_string())


if __name__ == "__main__":
    main()