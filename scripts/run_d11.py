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
    * data/processed/labels/<label>_cap_<net>.tsv  (per-(network, label)
      capped labels -- the file name MUST carry the network, see the fix
      below)
    * data/processed/splits/<label>_cap_<net>__native_<net>.tsv
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
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-rep", type=int, default=40,
                    help="repeats for the capped splits; must match the "
                         "primary design (40, work plan A2)")
    ap.add_argument("--splits-dir", default="data/processed/splits_v40",
                    help="directory holding the primary splits; the capped "
                         "splits are written here so the control runs on the "
                         "same design as the primary factorial")
    args = ap.parse_args()
    splits_dir = (ROOT / args.splits_dir) if args.splits_dir else PROCESSED_SPLIT
    splits_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: compute N per network from current label x network intersections.
    print("Step 1: compute per-network N (the minimum positives over labels)")
    counts = {}                                            # counts[net] = {label: int}
    net_gene_sets = {}                           # net_gene_sets[net] = set of genes
    for net_id in NETWORKS:
        genes, _ = load_network(str(PROCESSED_NET / f"{net_id}.tsv"))
        net_genes = set(genes)
        net_gene_sets[net_id] = net_genes
        counts[net_id] = {}
        for label in LABELS:
            label_genes = _label_genes(label)
            counts[net_id][label] = len(label_genes & net_genes)
        n_min = min(counts[net_id].values())
        print(f"  {net_id}: " + ", ".join(f"{l}={c}" for l, c in counts[net_id].items())
              + f"  -> N = {n_min}")

    # Step 2: for each (network, label) pair, write a capped label file
    # that contains only the top-N (rank-1..N) positives.
    #
    # FIX (P0-4): the file name used to be `{label}_cap` for EVERY
    # network, so each network OVERWROTE the previous network's cap file,
    # and the splits built in step 3 then used whichever cap happened to
    # be on disk at the time. The released d11_table.tsv claimed one N per
    # network while the splits actually used a different one (e.g. funmap
    # claimed N = 104 but its split held 42 positives).
    # The name must carry the network.
    print("\nStep 2: write capped label files (top-N positives only, "
          "one file per (network, label))")
    caps = {}                                               # caps[(net, label)] = N
    cap_labels = {}                                 # cap_labels[(net, label)] = id
    for net_id in NETWORKS:
        n_min = min(counts[net_id].values())
        for label in LABELS:
            caps[(net_id, label)] = n_min
            df = pd.read_csv(LABELS_DIR / f"{label}.tsv", sep="\t", dtype=str)
            # FIX (P0-4, second defect): the cap N is defined as
            #   N = min over labels of (positives in the network's LCC),
            # so the label must be restricted to the universe BEFORE the
            # top-N is taken. Taking the top-N of the FULL label file
            # lets genes that are not in the LCC into the cap, and the
            # split then holds fewer than N positives (funmap: 73 instead
            # of 104).
            df = df[df.gene.isin(net_gene_sets[net_id])]
            df = df.sort_values("rank", kind="stable")
            df_cap = df.head(n_min).copy()
            if len(df_cap) < n_min:
                raise AssertionError(
                    f"[{net_id} x {label}] only {len(df_cap)} of the label's "
                    f"positives are in the universe, cannot reach the "
                    f"cap N = {n_min}"
                )
            df_cap["source"] = f"{df_cap['source'].iloc[0]}_cap_{n_min}"
            df_cap["release"] = f"{df_cap['release'].iloc[0]}_cap_{net_id}"
            cap_label = f"{label}_cap_{net_id}"
            cap_labels[(net_id, label)] = cap_label
            df_cap.to_csv(LABELS_DIR / f"{cap_label}.tsv", sep="\t", index=False, lineterminator="\n")

    # Step 3: regenerate splits for every (capped_label, native_<net>) pair.
    #
    # FIX (P0-4): the old code skipped the split whenever the file
    # already existed. The stale splits (built from the colliding cap
    # files) were therefore never rebuilt, and the `[SKIP]` made that
    # invisible. The splits are now always rebuilt from the matching
    # per-network capped label, and their content is verified afterwards.
    print("\nStep 3: regenerate capped splits (one per (label, network), "
          f"n_rep={args.n_rep})")
    for net_id in NETWORKS:
        for label in LABELS:
            cap_label = cap_labels[(net_id, label)]
            split_path = splits_dir / f"{cap_label}__native_{net_id}.tsv"
            print(f"  building {split_path.name}")
            subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "build_splits.py"),
                 "--label", cap_label,
                 "--universe", f"native_{net_id}",
                 "--n-rep", str(args.n_rep),
                 "--out-dir", str(splits_dir)],
                check=True, cwd=str(ROOT))

            # VERIFY: the split must contain exactly the intended number
            # of positives. This is the check that would have caught the
            # original defect.
            sp = pd.read_csv(split_path, sep="\t")
            n_actual = int(sp.loc[sp.y.astype(int) == 1, "gene"].nunique())
            n_want = caps[(net_id, label)]
            if n_actual != n_want:
                raise AssertionError(
                    f"[{net_id} x {cap_label}] split holds {n_actual} "
                    f"positives, expected the capped N = {n_want}. The "
                    f"capped label file or the split is wrong."
                )

    # Step 4: run each arm (6 networks x 4 capped labels) at alpha=0.5.
    print(f"\nStep 4: run 24 capped arms at alpha=0.5")
    summary_rows = []
    for net_id in NETWORKS:
        for label in LABELS:
            cap_label = cap_labels[(net_id, label)]
            net_path = PROCESSED_NET / f"{net_id}.tsv"
            split_path = splits_dir / f"{cap_label}__native_{net_id}.tsv"
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
            res.to_csv(out_path, sep="\t", index=False, lineterminator="\n")
            # Extract delta_AUROC for this arm (alpha=0.5 row)
            from d1.engine.evaluate import margins
            m = margins(res)
            d = m[m.alpha == 0.5]
            n_test_pos = int(res[(res.alpha == 0.5) & (res.method == "rwr")].n_test_pos.iloc[0])
            n_test = int(res[(res.alpha == 0.5) & (res.method == "rwr")].n_test.iloc[0])
            if len(d):
                sp = pd.read_csv(split_path, sep="\t")
                n_actual_pos = int(
                    sp.loc[sp.y.astype(int) == 1, "gene"].nunique())
                summary_rows.append({
                    "network": net_id, "label": cap_label,
                    # the base label, so the table can be pivoted into a
                    # label x network grid without string surgery
                    "label_base": label,
                    "alpha": 0.5,
                    "N": caps[(net_id, label)],
                    "N_in_split": n_actual_pos,
                    "n_test_pos": n_test_pos,
                    "n_test": n_test,
                    "rwr_auroc": d.rwr_auroc.mean(),
                    "degree_auroc": d.degree_auroc.mean(),
                    "delta_auroc": d.delta_auroc.mean(),
                    "delta_auroc_sd": d.delta_auroc.std(),
                })
            print(f"  [{net_id} x {cap_label}] {len(res)} rows -> delta={d.delta_auroc.mean():.4f} +- {d.delta_auroc.std():.4f}")

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TABLE_DIR / "d11_table.tsv", sep="\t", index=False, lineterminator="\n")
    print("\nD-11 Table (delta_AUROC at alpha=0.5):")
    pivot = summary.pivot_table(index="network", columns="label_base",
                                 values="delta_auroc").round(4)
    print(pivot.to_string())
    print("\nCap N used per network (must equal the number of positives in "
          "the split):")
    print(summary.pivot_table(index="network", columns="label_base",
                               values="N").to_string())
    bad = summary[summary.N != summary.N_in_split]
    if len(bad):
        raise AssertionError(
            f"{len(bad)} arms where the split does not hold the intended N:\n"
            f"{bad[['network', 'label', 'N', 'N_in_split']].to_string()}"
        )
    print("\nAll 24 arms: the split holds exactly the intended N.")


if __name__ == "__main__":
    main()