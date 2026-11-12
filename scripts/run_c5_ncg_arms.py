"""C5 -- put the NCG label into the main factorial.

What this does
--------------
1. **The primary convention (uncapped).** The released 24-arm primary
   factorial runs each label with its own positive set (461 / 394 /
   374 / 104 positives per repeat in funmap); the D-11 common-size cap
   is a *control* arm, not the primary design (see `scripts/run_d11.py`
   and `results/tables/d11_table.tsv`). The NCG arm is therefore built
   the same way as the other 24: the label's full positive set, the
   same 5 folds x 40 repeats, the same three alphas.

2. **The D-11 control convention (capped).** The control exists to
   equalise the positive-set size across labels within an arm, so the
   NCG label gets a capped arm as well. The cap N per network is the
   minimum over the **four canonical labels** (it does not change when
   NCG is added: NCG has 407-694 positives per network, always above
   the temporal label's 63-150), and the cap is applied *after* the
   universe restriction (the P0-4 / A2.3i fix). Every capped split is
   verified to hold exactly N positives.

Outputs:
  data/processed/labels/ncg_known_cap_<net>.tsv
  data/processed/splits_v40/ncg_known__native_<net>.tsv         (primary)
  data/processed/splits_v40/ncg_known_cap_<net>__native_<net>.tsv (capped)
  results/runs_v40_ncg/a4_<net>__ncg_known.tsv, a4_combined.tsv
  results/runs/d11/d11_<net>__ncg_known_cap_<net>.tsv
  results/tables/a4_table1_delta_auroc_v2.tsv        (30 combos)
  results/tables/d11_table_v2.tsv                    (30 rows)
  results/runs_v40/a4_combined_v2.tsv                (30 combos, the
                                                      input for C6)
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.evaluate import margins, run_arm  # noqa: E402

NETWORKS = ["funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome"]
CANON = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
NCG_LABEL = "ncg_known"
ALPHAS = (0.3, 0.5, 0.7)

LABELS_DIR = ROOT / "data/processed/labels"
NET_DIR = ROOT / "data/processed/networks"
SPLITS = ROOT / "data/processed/splits_v40"
RUNS_NCG = ROOT / "results/runs_v40_ncg"
RUNS_D11 = ROOT / "results/runs/d11"
TABLES = ROOT / "results/tables"
RUNS_V40 = ROOT / "results/runs_v40"


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def net_genes(net_id: str) -> set:
    df = pd.read_csv(NET_DIR / f"{net_id}.tsv", sep="\t", dtype=str)
    return set(df.gene_a) | set(df.gene_b)


def label_genes(label_id: str) -> set:
    df = pd.read_csv(LABELS_DIR / f"{label_id}.tsv", sep="\t", dtype=str)
    return set(df.gene)


def build_splits(label_id: str, net_id: str) -> Path:
    subprocess.run(
        [sys.executable, str(ROOT / "scripts/build_splits.py"),
         "--label", label_id,
         "--universe", f"native_{net_id}",
         "--n-rep", "40",
         "--out-dir", str(SPLITS)],
        check=True, cwd=str(ROOT),
        stdout=subprocess.DEVNULL)
    return SPLITS / f"{label_id}__native_{net_id}.tsv"


def main() -> int:
    RUNS_NCG.mkdir(parents=True, exist_ok=True)
    RUNS_D11.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    # ---- the cap N per network (over the four canonical labels) ------
    n_cap = {}
    for net in NETWORKS:
        genes = net_genes(net)
        counts = {lab: len(label_genes(lab) & genes) for lab in CANON}
        n_cap[net] = min(counts.values())
        log(f"  {net}: canonical counts {counts} -> N = {n_cap[net]}")

    # ---- 1. the primary (uncapped) NCG arms --------------------------
    log("Step 1: primary NCG splits (uncapped, 5 folds x 40 repeats)")
    for net in NETWORKS:
        p = build_splits(NCG_LABEL, net)
        df = pd.read_csv(p, sep="\t")
        n_pos = int(df.loc[df.y.astype(int) == 1, "gene"].nunique())
        log(f"  {p.name}: {n_pos} positives in the universe")

    log("Step 2: run the 6 primary NCG arms (3 alphas)")
    for net in NETWORKS:
        split_path = SPLITS / f"{NCG_LABEL}__native_{net}.tsv"
        out_path = RUNS_NCG / f"a4_{net}__{NCG_LABEL}.tsv"
        t0 = time.time()
        res = run_arm(network=str(NET_DIR / f"{net}.tsv"),
                      splits=str(split_path),
                      network_id=net,
                      label_id=NCG_LABEL,
                      universe_id=f"native_{net}",
                      alphas=ALPHAS,
                      dump_scores=str(ROOT / f"results/scores_v40/"
                                           f"{net}__{NCG_LABEL}.npz"))
        res.to_csv(out_path, sep="\t", index=False, lineterminator="\n")
        log(f"  {net}: {len(res)} rows in {time.time() - t0:.1f}s")

    # ---- 3. the D-11 capped NCG arms ---------------------------------
    log("Step 3: D-11 capped NCG arms")
    d11_rows = []
    for net in NETWORKS:
        genes = net_genes(net)
        n = n_cap[net]
        df = pd.read_csv(LABELS_DIR / f"{NCG_LABEL}.tsv", sep="\t", dtype=str)
        df = df[df.gene.isin(genes)]          # universe FIRST (A2.3i)
        df = df.sort_values("rank", kind="stable")
        df_cap = df.head(n).copy()
        if len(df_cap) < n:
            raise AssertionError(
                f"[{net} x {NCG_LABEL}] only {len(df_cap)} of the label's "
                f"positives are in the universe, cannot reach N = {n}")
        df_cap["source"] = f"{df_cap['source'].iloc[0]}_cap_{n}"
        df_cap["release"] = f"{df_cap['release'].iloc[0]}_cap_{net}"
        cap_label = f"{NCG_LABEL}_cap_{net}"
        df_cap.to_csv(LABELS_DIR / f"{cap_label}.tsv", sep="\t", index=False,
                      lineterminator="\n")
        split_path = build_splits(cap_label, net)
        sp = pd.read_csv(split_path, sep="\t")
        n_actual = int(sp.loc[sp.y.astype(int) == 1, "gene"].nunique())
        if n_actual != n:
            raise AssertionError(
                f"[{net} x {cap_label}] split holds {n_actual} positives, "
                f"expected N = {n}")
        out_path = RUNS_D11 / f"d11_{net}__{cap_label}.tsv"
        res = run_arm(network=str(NET_DIR / f"{net}.tsv"),
                      splits=str(split_path),
                      network_id=net,
                      label_id=cap_label,
                      universe_id=f"native_{net}",
                      alphas=(0.5,))
        res.to_csv(out_path, sep="\t", index=False, lineterminator="\n")
        m = margins(res)
        d = m[m.alpha == 0.5]
        rwr = res[(res.alpha == 0.5) & (res.method == "rwr")]
        d11_rows.append({
            "network": net, "label": cap_label, "label_base": NCG_LABEL,
            "alpha": 0.5, "N": n, "N_in_split": n_actual,
            "n_test_pos": int(rwr.n_test_pos.iloc[0]),
            "n_test": int(rwr.n_test.iloc[0]),
            "rwr_auroc": float(d.rwr_auroc.mean()),
            "degree_auroc": float(d.degree_auroc.mean()),
            "delta_auroc": float(d.delta_auroc.mean()),
            "delta_auroc_sd": float(d.delta_auroc.std()),
        })
        log(f"  {net}: N = {n}, delta = {d11_rows[-1]['delta_auroc']:+.4f}")

    # ---- 4. the v2 tables --------------------------------------------
    # 4a. the 30-combination main table
    v1 = pd.read_csv(TABLES / "a4_table1_delta_auroc_v40.tsv", sep="\t")
    log(f"the released v40 Table 1 has {len(v1)} rows")
    ncg_rows = []
    for net in NETWORKS:
        p = RUNS_NCG / f"a4_{net}__{NCG_LABEL}.tsv"
        if not p.exists():
            continue
        res = pd.read_csv(p, sep="\t")
        m = margins(res)
        row = {"network": net, "label_set": NCG_LABEL}
        for a in ALPHAS:
            d = m[m.alpha == a]
            if d.empty:
                continue
            row[f"a{a}_mean"] = float(d.delta_auroc.mean())
            row[f"a{a}_std"] = float(d.delta_auroc.std())
            row[f"a{a}_n"] = float(len(d))
        ncg_rows.append(row)
    ncg_tab = pd.DataFrame(ncg_rows)
    log(f"the NCG arms add {len(ncg_tab)} rows")

    # the combined per-fold table (the input for the variance decomposition)
    comb_v1 = pd.read_csv(RUNS_V40 / "a4_combined.tsv", sep="\t")
    comb_ncg = pd.concat([pd.read_csv(RUNS_NCG / f"a4_{n}__{NCG_LABEL}.tsv",
                                       sep="\t")
                           for n in NETWORKS
                           if (RUNS_NCG / f"a4_{n}__{NCG_LABEL}.tsv").exists()],
                          ignore_index=True)
    comb_v2 = pd.concat([comb_v1, comb_ncg], ignore_index=True)
    comb_v2.to_csv(RUNS_V40 / "a4_combined_v2.tsv", sep="\t", index=False,
                   lineterminator="\n")
    log(f"wrote a4_combined_v2.tsv: {len(comb_v1):,} + {len(comb_ncg):,} "
        f"= {len(comb_v2):,} rows, "
        f"{comb_v2.groupby(['network', 'label_set']).ngroups} combos")

    # 4b. the D-11 table v2 (the 24 released rows + the 6 NCG rows)
    d11_v1 = pd.read_csv(TABLES / "d11_table.tsv", sep="\t")
    d11_v2 = pd.concat([d11_v1, pd.DataFrame(d11_rows)], ignore_index=True)
    d11_v2.to_csv(TABLES / "d11_table_v2.tsv", sep="\t", index=False,
                  lineterminator="\n")
    log(f"wrote d11_table_v2.tsv: {len(d11_v1)} + {len(d11_rows)} rows")

    # 4c. the main table v2, in the same shape as the released one
    t1_v1 = pd.read_csv(TABLES / "a4_table1_delta_auroc_v40.tsv", sep="\t")
    log(f"the released Table-1 columns: {list(t1_v1.columns)}")
    t1_v2 = pd.concat([t1_v1, ncg_tab], ignore_index=True)
    t1_v2.to_csv(TABLES / "a4_table1_delta_auroc_v2.tsv", sep="\t", index=False,
                 lineterminator="\n")
    log(f"wrote a4_table1_delta_auroc_v2.tsv: {len(t1_v2)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
