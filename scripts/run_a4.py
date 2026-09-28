"""A4: Wire all 6 networks into the harness end-to-end (full factorial).

This is the W3 task that produces **Table 1** of the paper: the main
benchmark result showing delta-AUROC for every (network, label) pair
under the primary analysis settings.

What this script does
---------------------
1. Generate the C4 split file for every (label, native_<net_id>) pair we
   do not already have on disk. We split on a per-call basis because
   build_splits.py is reproducible given the (label, universe) marker.
3. Run ``run_arm`` for every (network, label) combination:
      * 6 primary networks: funmap, rna_coexp, string_phys700, intact,
        string_full700, reactome
      * 4 primary labels:    intogen2024, ot_all, ot_nolit, intogen_temporal_new
      * 3 alpha values:      0.3, 0.5, 0.7 (sensitivity on RWR restart prob)
   -> 6 x 4 = 24 arms x 3 alphas x 50 folds x 2 methods = 7,200 rows.
4. Concatenate all per-arm C5 result tables into one big file.
5. Compute mean margins for the three alpha values and write Table 1
   (delta_AUROC +- SD) to ``results/tables/a4_table1_delta_auroc.tsv``.

Why we keep the per-arm C5 files as well
-----------------------------------------
The protocol asks for one C5 result file per arm (contract C5). This
script writes them under ``results/runs/`` and *also* aggregates them
into a single combined table for convenience. The aggregate does not
replace the per-arm files; both are committed.

Outputs
-------
* results/runs/a4_<network>__<label>__alpha<0.5>.tsv  (one per arm + alpha)
* results/runs/a4_combined.tsv                        (all rows in one)
* results/tables/a4_table1_delta_auroc.tsv            (Table 1, paper)
* results/tables/a4_table1_rwr_auroc.tsv             (raw RWR AUROC)
* results/tables/a4_table1_degree_auroc.tsv           (raw degree AUROC)
"""
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.evaluate import run_arm, summarize, margins  # noqa: E402

NETWORKS = [
    "funmap", "rna_coexp", "string_phys700",
    "intact", "string_full700", "reactome",
]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
ALPHAS = (0.3, 0.5, 0.7)                       # primary is 0.5
N_REP = 10
K = 5

PROCESSED_NET = ROOT / "data" / "processed" / "networks"
PROCESSED_SPLIT = ROOT / "data" / "processed" / "splits"
RUNS_DIR = ROOT / "results" / "runs"
TABLE_DIR = ROOT / "results" / "tables"
RUNS_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)


def _split_path(label: str, net_id: str) -> Path:
    return PROCESSED_SPLIT / f"{label}__native_{net_id}.tsv"


def _ensure_split(label: str, net_id: str) -> Path:
    """Run build_splits.py if the split file is missing."""
    p = _split_path(label, net_id)
    if p.exists():
        return p
    script = ROOT / "scripts" / "build_splits.py"
    print(f"  generating {p.name} via {script.name} ...")
    subprocess.run(
        [sys.executable, str(script),
         "--label", label,
         "--universe", f"native_{net_id}"],
        check=True, cwd=str(ROOT),
    )
    return p


def _run_arm(net_id: str, label: str, alphas: tuple[float, ...]) -> Path:
    """Run one network x label x alpha-set arm; return the output path."""
    net_path = PROCESSED_NET / f"{net_id}.tsv"
    split_path = _ensure_split(label, net_id)
    out_path = RUNS_DIR / f"a4_{net_id}__{label}.tsv"
    t0 = time.time()
    res = run_arm(
        network=str(net_path),
        splits=str(split_path),
        network_id=net_id,
        label_id=label,
        universe_id=f"native_{net_id}",
        alphas=alphas,
    )
    res.to_csv(out_path, sep="\t", index=False)
    dt = time.time() - t0
    n = len(res)
    print(f"  [{net_id} x {label}] {n} rows in {dt:.1f}s -> {out_path.name}")
    return out_path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only-net", action="append", default=None,
                    help="restrict to one or more net_ids (default: all 6)")
    ap.add_argument("--only-label", action="append", default=None,
                    help="restrict to one or more label_ids (default: all 4)")
    ap.add_argument("--alpha", type=float, nargs="+", default=list(ALPHAS))
    ap.add_argument("--n-rep", type=int, default=N_REP)
    ap.add_argument("--k", type=int, default=K)
    args = ap.parse_args(argv)

    nets = args.only_net or NETWORKS
    labels = args.only_label or LABELS
    alphas = tuple(args.alpha)
    print(f"A4 factorial: {len(nets)} networks x {len(labels)} labels x "
          f"{len(alphas)} alphas = {len(nets) * len(labels)} arms")

    # Make sure every needed split file exists.
    print("\nStep 1: ensure split files exist (generate if missing)")
    for net in nets:
        for label in labels:
            _ensure_split(label, net)

    # Run each arm.
    print(f"\nStep 2: run {len(nets) * len(labels)} arms (one C5 file per arm)")
    arm_paths = []
    for net in nets:
        for label in labels:
            arm_paths.append(_run_arm(net, label, alphas))

    # Concatenate every C5 row into one combined file.
    print(f"\nStep 3: concatenate {len(arm_paths)} arm tables into results/runs/a4_combined.tsv")
    combined = pd.concat([pd.read_csv(p, sep="\t") for p in arm_paths],
                         ignore_index=True)
    combined.to_csv(RUNS_DIR / "a4_combined.tsv", sep="\t", index=False)
    print(f"  combined rows: {len(combined):,}")

    # Build Table 1 (paper Table 1): mean +- SD of delta_AUROC per (net, label, alpha).
    # IMPORTANT: pivot on `margins()` not `summarize()`. summarize() has one
    # row per (network, label, alpha) so its std is undefined; margins()
    # has one row per fold so the std is meaningful.
    print(f"\nStep 4: aggregate to Table 1 -> results/tables/a4_table1_delta_auroc.tsv")
    m = margins(combined)
    g = m.groupby(["network", "label_set", "alpha"])["delta_auroc"]
    table1_long = pd.DataFrame({
        "mean": g.mean(),
        "std": g.std(),
        "n": g.size(),
    }).reset_index()
    table1 = table1_long.pivot_table(
        index=["network", "label_set"],
        columns="alpha",
        values=["mean", "std", "n"],
    ).round(4)
    table1.columns = [f"a{alpha}_{stat}" for stat, alpha in table1.columns]
    table1 = table1.reindex(columns=[
        "a0.3_mean", "a0.3_std", "a0.3_n",
        "a0.5_mean", "a0.5_std", "a0.5_n",
        "a0.7_mean", "a0.7_std", "a0.7_n",
    ])
    table1.to_csv(TABLE_DIR / "a4_table1_delta_auroc.tsv", sep="\t")
    print("  Table 1 (delta_AUROC, mean +- SD across 50 folds):")
    print(table1.to_string())

    # Raw RWR and degree AUROC tables (for the supplementary).
    rwr_df = combined[combined["method"] == "rwr"]
    table1_rwr = rwr_df.groupby(["network", "label_set", "alpha"])["auroc"].agg(
        ["mean", "std"]
    ).round(4).reset_index().pivot_table(
        index=["network", "label_set"], columns="alpha", values=["mean", "std"]
    )
    table1_rwr.columns = [f"a{alpha}_{stat}" for stat, alpha in table1_rwr.columns]
    table1_rwr.to_csv(TABLE_DIR / "a4_table1_rwr_auroc.tsv", sep="\t")

    deg_df = combined[combined["method"] == "degree"]
    table1_deg = deg_df.groupby(["network", "label_set"])["auroc"].agg(
        ["mean", "std"]
    ).round(4)
    table1_deg.columns = [f"degree_{stat}" for stat in table1_deg.columns]
    table1_deg.to_csv(TABLE_DIR / "a4_table1_degree_auroc.tsv", sep="\t")

    print(f"\nDone. {len(arm_paths)} arms; combined file {len(combined):,} rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())