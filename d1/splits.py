"""Stratified CV split generator (B5, contract C4).

Produces one ``data/processed/splits/{label_id}__{universe_id}.tsv`` file per
(label, universe) pair. Every split file is fully determined by:

* the universe (the list of genes to score)
* the positive set (capped to the universe; everything else is negative)
* a base random seed, written once into decisions.md (D-BASE)
* the number of repeats (default 10) and folds (default 5)

Re-running this module with identical inputs produces byte-identical files. The
split files are released as artefacts rather than regenerated from seeds, per
protocol section 4.7.

Decision D-06 says: ``seeds = training-fold positives; they are not in the
ranked test set.`` The harness enforces that on read. Here we only generate the
gene/repeat/fold/y table.

Decision D-BASE (to be confirmed at S0 with Person A) is the integer seed used
for ``StratifiedKFold(random_state=BASE + r)``; we keep it here as a module
constant for now and write the chosen value into docs/decisions.md.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

SPLIT_COLUMNS = ["gene", "repeat", "fold", "y"]
DEFAULT_BASE_SEED = 20260926   # B's candidate; mark D-BASE in decisions.md
DEFAULT_N_REP = 10
DEFAULT_K = 5


def make_splits(universe, positives, *, base_seed: int = DEFAULT_BASE_SEED,
                n_rep: int = DEFAULT_N_REP, k: int = DEFAULT_K) -> pd.DataFrame:
    """Build a C4 split table.

    Parameters
    ----------
    universe : iterable of str
        Every gene that may appear in any test fold. Genes absent here are
        excluded from the CV entirely.
    positives : iterable of str
        Set of positive genes. A gene is positive iff it appears in *both*
        ``universe`` and ``positives``. Everything else in ``universe`` is a
        negative.
    base_seed, n_rep, k : see DEFAULT_*.

    Returns
    -------
    DataFrame with columns ``gene repeat fold y`` and a stable row order.
    """
    universe = np.array(sorted(set(universe)), dtype=object)
    pos = set(positives) & set(universe)
    y = np.array([g in pos for g in universe], dtype=int)

    rows = []
    for r in range(n_rep):
        skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=base_seed + r)
        for f, (_, test_idx) in enumerate(skf.split(universe, y)):
            test_idx = np.asarray(test_idx)
            for i in test_idx:
                rows.append((str(universe[i]), r, f, int(y[i])))
    df = pd.DataFrame(rows, columns=SPLIT_COLUMNS)
    df = df.sort_values(["repeat", "fold", "gene"], kind="mergesort").reset_index(drop=True)
    return df


def write_splits(df: pd.DataFrame, label_id: str, universe_id: str,
                 out_dir: str | Path = "data/processed/splits") -> Path:
    """Write ``df`` as a C4 file under ``out_dir``.

    Returns the written path.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{label_id}__{universe_id}.tsv"
    df[SPLIT_COLUMNS].to_csv(path, sep="\t", index=False)
    return path


def load_splits(path: str | Path) -> pd.DataFrame:
    """Read a C4 file; returns the DataFrame unchanged."""
    return pd.read_csv(path, sep="\t", dtype=str)


def summarise_splits(df: pd.DataFrame) -> dict:
    """Small per-arm diagnostic summary, useful at S2 to compare the panels."""
    df = df.copy()
    df["y"] = df["y"].astype(int)
    n_total = df.drop_duplicates("gene").shape[0]
    n_pos = int(df.loc[df.y == 1, "gene"].drop_duplicates().shape[0])
    by_fold = df.groupby(["repeat", "fold"]).y.agg(["sum", "size"])
    pos_per_fold = by_fold["sum"]
    return {
        "n_genes": int(n_total),
        "n_positives": n_pos,
        "n_repeats": int(df.repeat.nunique()),
        "n_folds": int(df.fold.nunique()),
        "pos_per_fold_min": int(pos_per_fold.min()),
        "pos_per_fold_max": int(pos_per_fold.max()),
        "pos_per_fold_mean": float(pos_per_fold.mean()),
    }
