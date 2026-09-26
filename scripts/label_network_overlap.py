"""B6: label-network overlap report and common-size capping.

For each (label, network) pair, count how many of the label's positive genes
land inside the network's LCC (the universe that RWR will rank on). This
**alone** moves absolute AUROC per protocol §4.2, so it is reported as a table
before any scoring is run. Then cap every label set to the same number N of
positives (the smallest count observed in the shared universe), so the four
arms differ in composition, not in positive-class prevalence.

Outputs
-------
* ``results/tables/label_network_overlap.tsv`` (one row per label x network)
* ``data/processed/labels/{label_id}.capped.tsv`` (capped to N positives)
* A short summary printed to stdout.

Decisions referenced
--------------------
* D-OT-CAP initial top-N = 600. We recompute N at the end as the minimum
  (label_count_in_shared) across all four label sets.
* D-11 placeholder for the agreed N.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
LABELS_DIR = REPO_ROOT / "data" / "processed" / "labels"
NETWORKS_DIR = REPO_ROOT / "data" / "processed" / "networks"
RESULTS_DIR = REPO_ROOT / "results" / "tables"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

LABEL_FILES = {
    "intogen2024": "intogen2024.tsv",
    "ot_all": "ot_all.tsv",
    "ot_nolit": "ot_nolit.tsv",
    "intogen_temporal_new": "intogen_temporal_new.tsv",
}


def _read_network_genes(net_id: str) -> set[str]:
    p = NETWORKS_DIR / f"{net_id}.tsv"
    if not p.exists():
        raise FileNotFoundError(f"network file {p} not found; "
                                "wait for Person A's network panel")
    df = pd.read_csv(p, sep="\t", usecols=["gene_a", "gene_b"], dtype=str)
    return set(df.gene_a) | set(df.gene_b)


def _read_universe_shared() -> set[str] | None:
    """Read the shared universe if Person A has produced it."""
    p = NETWORKS_DIR / "universe_shared.txt"
    if not p.exists():
        return None
    return {g.strip() for g in open(p) if g.strip()}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--network", action="append",
                    help="restrict to one or more net_ids (default: all found in data/processed/networks)")
    args = ap.parse_args(argv)

    # discover networks
    if args.network:
        net_ids = args.network
    else:
        net_ids = sorted(p.stem for p in NETWORKS_DIR.glob("*.tsv")
                         if not p.stem.endswith("universe"))
    if not net_ids:
        print("!! no network files found under", NETWORKS_DIR)
        return 1
    net_genes = {n: _read_network_genes(n) for n in net_ids}
    shared = _read_universe_shared()

    # discover labels
    label_sets: dict[str, set[str]] = {}
    for label_id, fn in LABEL_FILES.items():
        p = LABELS_DIR / fn
        if not p.exists():
            print(f"  [skip missing] {label_id} -> {p}")
            continue
        df = pd.read_csv(p, sep="\t", dtype=str)
        label_sets[label_id] = set(df.gene.dropna().astype(str))
    if not label_sets:
        print("!! no label files found under", LABELS_DIR)
        return 1

    # build the overlap table
    rows = []
    for lid, lgenes in label_sets.items():
        for nid, ngenes in net_genes.items():
            overlap = lgenes & ngenes
            rows.append({
                "label_set": lid,
                "network": nid,
                "n_label_positives": len(lgenes),
                "n_network_lcc": len(ngenes),
                "n_overlap": len(overlap),
                "frac_label_in_network": round(len(overlap) / max(1, len(lgenes)), 4),
            })
    overlap_df = pd.DataFrame(rows)
    out_csv = RESULTS_DIR / "label_network_overlap.tsv"
    overlap_df.to_csv(out_csv, sep="\t", index=False)
    print(f"Wrote {out_csv}")

    # common-size cap using the shared universe, if available
    if shared is not None:
        in_shared = {lid: sorted(genes & shared) for lid, genes in label_sets.items()}
        counts = {lid: len(g) for lid, g in in_shared.items()}
        n = min(counts.values()) if counts else 0
        print(f"\nShared universe: {len(shared):,} genes")
        print(f"Per-label positives in shared universe:")
        for lid, c in sorted(counts.items()):
            print(f"  {lid:25s} {c:5d}")
        print(f"\nCommon cap N = {n} (D-11 candidate)")
        # write capped files
        for lid, sorted_genes in in_shared.items():
            head = sorted_genes[:n]
            df_full = pd.read_csv(LABELS_DIR / LABEL_FILES[lid], sep="\t", dtype=str)
            head_set = set(head)
            df_capped = df_full[df_full.gene.isin(head_set)].copy()
            df_capped["rank"] = range(1, len(df_capped) + 1)
            out = LABELS_DIR / f"{lid}.capped.tsv"
            df_capped[["gene", "rank", "source", "release"]].to_csv(out, sep="\t", index=False)
            print(f"  {lid:25s} -> {out.relative_to(REPO_ROOT)} ({len(df_capped)} genes)")
    else:
        print("\nNo universe_shared.txt yet; overlap table written without capping.")
        print("Re-run after Person A writes data/processed/networks/universe_shared.txt")

    # also: print a pivot for quick review
    print("\nLabel-network overlap (rows = label, cols = network):")
    pivot = overlap_df.pivot(index="label_set", columns="network", values="n_overlap")
    print(pivot.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
