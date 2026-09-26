"""B5 end-to-end: produce C4 split files for one or more (label, universe) pairs.

Given a label file (C3 schema) and a universe (typically ``native_<net_id>``
or ``shared``), writes ``data/processed/splits/{label_id}__{universe_id}.tsv``
under contract C4. One invocation can produce several split files at once.

Usage
-----
    python scripts/build_splits.py --label intogen2024 --universe native_funmap
    python scripts/build_splits.py --label intogen2024 --label ot_all \\
            --universe native_funmap --universe shared

Network LCCs are read from ``data/processed/networks/{net_id}.tsv``; for the
``shared`` universe the script reads ``data/processed/networks/universe_shared.txt``
when it exists, otherwise it falls back to the LCC of the first network found
(clearly logged).

The split files are deterministic given ``base_seed`` (D-BASE = 20260926) and
the (label, universe) pair, and pass ``d1.engine.evaluate.check_splits``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
LABELS_DIR = REPO_ROOT / "data" / "processed" / "labels"
NETWORKS_DIR = REPO_ROOT / "data" / "processed" / "networks"
SPLITS_DIR = REPO_ROOT / "data" / "processed" / "splits"
SPLITS_DIR.mkdir(parents=True, exist_ok=True)


def _read_label_genes(label_id: str) -> set[str]:
    p = LABELS_DIR / f"{label_id}.tsv"
    if not p.exists():
        raise FileNotFoundError(p)
    return set(pd.read_csv(p, sep="\t", dtype=str).gene.dropna().astype(str))


def _read_native_universe(net_id: str) -> set[str]:
    p = NETWORKS_DIR / f"{net_id}.tsv"
    if not p.exists():
        raise FileNotFoundError(p)
    df = pd.read_csv(p, sep="\t", dtype=str)
    return set(df.gene_a) | set(df.gene_b)


def _read_shared_universe() -> set[str]:
    p = NETWORKS_DIR / "universe_shared.txt"
    if p.exists():
        return {g.strip() for g in open(p) if g.strip()}
    raise FileNotFoundError(
        "data/processed/networks/universe_shared.txt not found; "
        "Person A's network panel must land before the shared universe exists"
    )


def _resolve_universe(spec: str) -> tuple[str, set[str]]:
    """Return (universe_id, gene_set) for a universe spec like ``native_funmap`` or ``shared``."""
    if spec == "shared":
        return "shared", _read_shared_universe()
    if spec.startswith("native_"):
        net_id = spec[len("native_"):]
        return spec, _read_native_universe(net_id)
    raise ValueError(f"unknown universe spec {spec!r}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--label", action="append", required=True,
                    help="label_id (e.g. intogen2024, ot_all, ot_nolit, intogen_temporal_new); "
                         "may be repeated")
    ap.add_argument("--universe", action="append", required=True,
                    help="universe spec (e.g. shared, native_funmap); may be repeated")
    ap.add_argument("--base-seed", type=int, default=20260926,
                    help="D-BASE; default 20260926")
    ap.add_argument("--n-rep", type=int, default=10)
    ap.add_argument("--k", type=int, default=5)
    args = ap.parse_args(argv)

    sys.path.insert(0, str(REPO_ROOT))
    from d1.splits import make_splits, write_splits, summarise_splits

    written = []
    for label_id in args.label:
        try:
            positives = _read_label_genes(label_id)
        except FileNotFoundError as e:
            print(f"!! {e}")
            continue
        for u_spec in args.universe:
            try:
                u_id, universe = _resolve_universe(u_spec)
            except FileNotFoundError as e:
                print(f"  [skip] {label_id} x {u_spec}: {e}")
                continue
            df = make_splits(universe, positives,
                             base_seed=args.base_seed, n_rep=args.n_rep, k=args.k)
            path = write_splits(df, label_id, u_id, out_dir=SPLITS_DIR)
            summ = summarise_splits(df)
            print(f"[{label_id} x {u_id}] {len(universe):,} genes, {len(positives):,} positives, "
                  f"{summ['n_positives']} positives in universe -> {path.relative_to(REPO_ROOT)} "
                  f"({len(df):,} rows)")
            written.append(path)

    if not written:
        print("!! no split files written")
        return 1
    print(f"\nWrote {len(written)} split files to {SPLITS_DIR.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
