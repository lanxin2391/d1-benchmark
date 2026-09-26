"""End-to-end check of the harness on synthetic networks with known answers (W0, task A1).

Runs the full pipeline exactly as it will run on real data:
    toy edges -> finalize_edges (C2 file) -> split file (C4) -> run_arm -> result table (C5)
for three scenarios (see d1/toy.py) and three alphas.

Usage:  python scripts/run_toy_benchmark.py
Output: results/toy/networks/, results/toy/splits/, results/toy/toy_runs.tsv,
        results/toy/toy_summary.tsv, and a summary table printed to the screen.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd  # noqa: E402

from d1.engine.evaluate import run_arm, summarize  # noqa: E402
from d1.networks.io import check_network_file, finalize_edges, load_network  # noqa: E402
from d1.toy import make_toy_splits, toy_network  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "toy")
EXPECT = {
    "random": "RWR ~0.5, degree ~0.5, margin ~0  (no signal anywhere)",
    "hubs":   "degree very high, margin ~0 or negative  (prior = degree proxy)",
    "module": "degree ~0.5, RWR high, margin clearly > 0  (real topology signal)",
}


def main():
    t0 = time.time()
    all_res = []
    for scenario in ("random", "hubs", "module"):
        edges, pos = toy_network(scenario, n=3000, n_pos=150, seed=0)
        net_id = f"toy_{scenario}"
        _, meta = finalize_edges(edges, net_id, {"source_file": "synthetic"},
                                 out_dir=os.path.join(OUT, "networks"))
        net_path = os.path.join(OUT, "networks", f"{net_id}.tsv")
        check_network_file(net_path)
        genes, A = load_network(net_path)

        splits = make_toy_splits(genes, pos & set(genes))
        os.makedirs(os.path.join(OUT, "splits"), exist_ok=True)
        split_path = os.path.join(OUT, "splits", f"toy__native_{net_id}.tsv")
        splits.to_csv(split_path, sep="\t", index=False)

        res = run_arm(net_path, split_path, net_id, "toy_labels", f"native_{net_id}",
                      alphas=(0.3, 0.5, 0.7))
        all_res.append(res)
        print(f"{net_id}: {meta['n_nodes_lcc']} genes, {meta['n_edges_lcc']} edges, "
              f"{len(pos)} positives")

    res = pd.concat(all_res, ignore_index=True)
    res.to_csv(os.path.join(OUT, "toy_runs.tsv"), sep="\t", index=False)
    summ = summarize(res)
    summ.to_csv(os.path.join(OUT, "toy_summary.tsv"), sep="\t", index=False)

    pd.set_option("display.width", 160)
    print("\nMean over 5 folds x 10 repeats:")
    print(summ[["network", "alpha", "rwr_auroc", "degree_auroc", "delta_auroc",
                "delta_auroc_sd"]].round(3).to_string(index=False))
    print("\nWhat you should see:")
    for k, v in EXPECT.items():
        print(f"  toy_{k:7s} {v}")
    print(f"\nDone in {time.time() - t0:.1f} s. Results in results/toy/")


if __name__ == "__main__":
    main()
