"""Build the 40-repeat split panel (``data/processed/splits_v40/``) and verify it.

Why
---
The released split panel (``data/processed/splits/``) uses 10 repeats x 5 folds.
Work-plan v3.1 task A2 raised the repeat count to 40, and that is the design the
stability analysis needs: a per-gene stability statistic is estimated from the
across-repeat variability of that gene's rank, so the number of repeats is the
sample size of the *dependent variable* of this analysis.

The new panel is written to a **separate directory** and the released panel is
never touched (work-plan v3.1 risk R3: "all new splits use a new directory +
version number, v1 is never modified").

Determinism
-----------
``d1.splits.make_splits`` draws ``StratifiedKFold(random_state=base_seed + r)``
with ``base_seed = 20260926`` (decision D-BASE). The repeat index ``r`` is the
*only* thing that changes between the 10-repeat and the 40-repeat panels, so
repeats 0-9 of the new panel must be identical to the released panel. This
script **checks** that instead of assuming it: every row of repeats 0-9 is
compared against the released file, row by row, and the script fails loudly on
any mismatch.

The label positive sets are the plain label files (no D-11 cap): the cap is a
*control* arm (D-11) and lives in ``data/processed/splits_d11/``. The primary
factorial uses the uncapped positive set of each (label, network) pair, which is
what the released primary split files contain.

Usage
-----
    python scripts/build_splits_v40.py                 # build + verify all 24
    python scripts/build_splits_v40.py --verify-only   # only re-check vs v1
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from d1.splits import make_splits, summarise_splits, write_splits  # noqa: E402

LABELS_DIR = REPO_ROOT / "data" / "processed" / "labels"
NETWORKS_DIR = REPO_ROOT / "data" / "processed" / "networks"
V1_DIR = REPO_ROOT / "data" / "processed" / "splits"
V40_DIR = REPO_ROOT / "data" / "processed" / "splits_v40"

NETWORKS = ("funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome")
LABELS = ("intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new")
BASE_SEED = 20260926
K = 5


def _label_genes(label_id: str) -> set[str]:
    df = pd.read_csv(LABELS_DIR / f"{label_id}.tsv", sep="\t", dtype=str)
    return set(df.gene.dropna().astype(str))


def _network_genes(net_id: str) -> set[str]:
    df = pd.read_csv(NETWORKS_DIR / f"{net_id}.tsv", sep="\t", dtype=str)
    return set(df.gene_a) | set(df.gene_b)


def _verify_prefix(v40_path: Path, v1_path: Path, n_shared_rep: int) -> None:
    """The repeats 0..n_shared_rep-1 of the v40 panel must equal the v1 panel."""
    if not v1_path.exists():
        raise FileNotFoundError(f"released panel missing: {v1_path}")
    v40 = pd.read_csv(v40_path, sep="\t", dtype={"gene": str})
    v1 = pd.read_csv(v1_path, sep="\t", dtype={"gene": str})
    a = v40[v40.repeat < n_shared_rep].reset_index(drop=True)
    b = v1[v1.repeat < n_shared_rep].reset_index(drop=True)
    if len(a) != len(b):
        raise AssertionError(f"{v40_path.name}: row count {len(a)} != v1 {len(b)}")
    if not a.equals(b):
        bad = a.compare(b)
        raise AssertionError(
            f"{v40_path.name}: repeats 0-{n_shared_rep - 1} differ from the released "
            f"panel ({len(bad)} differing cells); first diff:\n{bad.head()}"
        )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-rep", type=int, default=40)
    ap.add_argument("--k", type=int, default=K)
    ap.add_argument("--base-seed", type=int, default=BASE_SEED)
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    args = ap.parse_args(argv)

    nets = args.only_net or list(NETWORKS)
    labels = args.only_label or list(LABELS)
    V40_DIR.mkdir(parents=True, exist_ok=True)

    n_shared = 10  # repeats 0..9 must be shared with the released 10-repeat panel
    rows = []
    for label in labels:
        pos = _label_genes(label)
        for net in nets:
            universe = _network_genes(net)
            name = f"{label}__native_{net}"
            v40_path = V40_DIR / f"{name}.tsv"
            v1_path = V1_DIR / f"{name}.tsv"
            t0 = time.time()
            if not args.verify_only:
                df = make_splits(universe, pos, base_seed=args.base_seed,
                                 n_rep=args.n_rep, k=args.k)
                write_splits(df, label, f"native_{net}", out_dir=V40_DIR)
            _verify_prefix(v40_path, v1_path, n_shared)
            summ = summarise_splits(pd.read_csv(v40_path, sep="\t", dtype={"gene": str}))
            dt = time.time() - t0
            print(f"[{name}] ok  genes={summ['n_genes']:,} pos={summ['n_positives']} "
                  f"rep={summ['n_repeats']} fold={summ['n_folds']} "
                  f"pos/fold {summ['pos_per_fold_min']}-{summ['pos_per_fold_max']} "
                  f"({dt:.1f}s)")
            rows.append({"label": label, "universe": f"native_{net}",
                         "n_universe": summ["n_genes"], "n_positives": summ["n_positives"],
                         "n_repeats": summ["n_repeats"], "n_folds": summ["n_folds"],
                         "rows": summ["n_genes"] * summ["n_repeats"],
                         "sha256": _sha256(v40_path)})

    log = pd.DataFrame(rows)
    log_path = V40_DIR / "build_log_v40.tsv"
    log.to_csv(log_path, sep="\t", index=False, lineterminator="\n")
    print(f"\nWrote {len(rows)} split files + {log_path.relative_to(REPO_ROOT)}")
    return 0


def _sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


if __name__ == "__main__":
    raise SystemExit(main())
