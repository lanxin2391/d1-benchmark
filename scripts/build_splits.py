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

D-11 common-size cap
--------------------
With ``--apply-d11-cap`` the script first computes, for each requested universe,
the minimum number of positives across all four canonical labels
(intogen2024, ot_all, ot_nolit, intogen_temporal_new) — call that N_u — and
caps every label's positives to that N_u *before* generating splits. This
makes positive-class prevalence identical across labels within the same arm,
which is what the protocol's §4.2 "common-size cap" means. The cap is
**per-universe**, not global, so sparser networks still keep more positives
than a global cap would allow.

A companion ``data/processed/splits/N_per_universe.json`` records the chosen
N_u for downstream audit.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
LABELS_DIR = REPO_ROOT / "data" / "processed" / "labels"
NETWORKS_DIR = REPO_ROOT / "data" / "processed" / "networks"
SPLITS_DIR = REPO_ROOT / "data" / "processed" / "splits"
SPLITS_DIR.mkdir(parents=True, exist_ok=True)

CANONICAL_D11_LABELS = ("intogen2024", "ot_all", "ot_nolit",
                        "intogen_temporal_new")


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
    """Return (universe_id, gene_set) for a universe spec like ``native_<net>`` or ``shared``."""
    if spec == "shared":
        return "shared", _read_shared_universe()
    if spec.startswith("native_"):
        net_id = spec[len("native_"):]
        return spec, _read_native_universe(net_id)
    raise ValueError(f"unknown universe spec {spec!r}")


def _cap_label_genes(label_pos: set[str], n_cap: int,
                     label_file: Path | None = None) -> set[str]:
    """Cap a positive set to its top-N by the label file's ``rank`` column.

    If the label file has more than ``n_cap`` rows, take the first ``n_cap``
    by ``rank`` (C3 contract: rank starts at 1, ascending). If it has fewer
    than ``n_cap``, the set is unchanged.
    """
    if label_file is None or not label_file.exists():
        # fall back: take any n_cap elements (stable, sort)
        return set(sorted(label_pos)[:n_cap])
    df = pd.read_csv(label_file, sep="\t", dtype=str)
    df = df.sort_values("rank", kind="mergesort")
    keep = set(df.gene.head(n_cap))
    return label_pos & keep


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
    ap.add_argument("--apply-d11-cap", action="store_true",
                    help="cap each label's positives to the per-universe minimum N_u "
                         "(D-11 proposal A, B-side preferred)")
    args = ap.parse_args(argv)

    sys.path.insert(0, str(REPO_ROOT))
    from d1.splits import make_splits, write_splits, summarise_splits

    # ---- D-11 precompute: per-universe N_u = min(positives_in_universe across canonical labels)
    d11_N: dict[str, int] = {}
    if args.apply_d11_cap:
        # build (label -> positives) once for the four canonical labels
        canon_pos = {}
        for lab in CANONICAL_D11_LABELS:
            try:
                canon_pos[lab] = _read_label_genes(lab)
            except FileNotFoundError:
                print(f"  [D-11] skip {lab}: file not found")
        for u_spec in args.universe:
            try:
                _, universe = _resolve_universe(u_spec)
            except FileNotFoundError as e:
                print(f"  [D-11] skip {u_spec}: {e}")
                continue
            counts = {lab: len(p & universe) for lab, p in canon_pos.items()}
            n_u = min(counts.values()) if counts else 0
            d11_N[u_spec] = n_u
            print(f"[D-11] {u_spec}: per-label counts = {counts}, N_u = {n_u}")

    written = []
    n_log = []
    for label_id in args.label:
        try:
            positives_full = _read_label_genes(label_id)
        except FileNotFoundError as e:
            print(f"!! {e}")
            continue
        label_file = LABELS_DIR / f"{label_id}.tsv"
        for u_spec in args.universe:
            try:
                u_id, universe = _resolve_universe(u_spec)
            except FileNotFoundError as e:
                print(f"  [skip] {label_id} x {u_spec}: {e}")
                continue

            positives = positives_full
            cap_note = ""
            if args.apply_d11_cap and u_spec in d11_N:
                n_u = d11_N[u_spec]
                positives = _cap_label_genes(positives_full, n_u, label_file)
                cap_note = f" [D-11 cap → {len(positives)}]"

            df = make_splits(universe, positives,
                             base_seed=args.base_seed, n_rep=args.n_rep, k=args.k)
            path = write_splits(df, label_id, u_id, out_dir=SPLITS_DIR)
            summ = summarise_splits(df)
            print(f"[{label_id} x {u_id}] {len(universe):,} genes, "
                  f"{len(positives):,} positives (of {len(positives_full):,} full), "
                  f"{summ['n_positives']} positives in universe -> "
                  f"{path.relative_to(REPO_ROOT)} ({len(df):,} rows){cap_note}")
            written.append(path)
            n_log.append({"label": label_id, "universe": u_id,
                          "n_universe": len(universe),
                          "n_positives_full": len(positives_full),
                          "n_positives_used": len(positives),
                          "rows": len(df)})

    if not written:
        print("!! no split files written")
        return 1

    # record D-11 audit
    if args.apply_d11_cap:
        n_json = SPLITS_DIR / "N_per_universe.json"
        n_json.write_text(json.dumps({
            "scheme": "per_universe_min_across_4_canonical_labels",
            "canonical_labels": list(CANONICAL_D11_LABELS),
            "N_per_universe": d11_N,
        }, indent=2))
        print(f"\nWrote D-11 audit to {n_json.relative_to(REPO_ROOT)}")

    log_df = pd.DataFrame(n_log)
    log_path = SPLITS_DIR / "build_log.tsv"
    log_df.to_csv(log_path, sep="\t", index=False)
    print(f"\nWrote {len(written)} split files to "
          f"{SPLITS_DIR.relative_to(REPO_ROOT)} + audit log")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
