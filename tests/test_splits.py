"""Tests for the stratified CV split generator (B5, contract C4)."""
from __future__ import annotations

import pandas as pd
import pytest

from d1.engine.evaluate import check_splits
from d1.splits import (DEFAULT_BASE_SEED, make_splits, SPLIT_COLUMNS,
                       summarise_splits, write_splits)


def _toy_universe(n=200, n_pos=20):
    genes = [f"G{i:04d}" for i in range(n)]
    pos = set(genes[:n_pos])
    return genes, pos


def test_split_columns_and_shape():
    genes, pos = _toy_universe()
    df = make_splits(genes, pos, n_rep=3, k=5)
    assert list(df.columns) == SPLIT_COLUMNS
    assert set(df.repeat.unique()) == {0, 1, 2}
    assert set(df.fold.unique()) == {0, 1, 2, 3, 4}
    assert len(df) == 3 * len(genes)


def test_splits_are_deterministic():
    genes, pos = _toy_universe()
    a = make_splits(genes, pos, n_rep=2, k=5)
    b = make_splits(genes, pos, n_rep=2, k=5)
    pd.testing.assert_frame_equal(a, b)


def test_passes_check_splits():
    genes, pos = _toy_universe()
    df = make_splits(genes, pos, n_rep=2, k=5)
    check_splits(df)                                       # contract C4


def test_gene_outside_universe_is_dropped():
    genes, pos = _toy_universe()
    pos_with_extra = pos | {"G9999", "FOO"}
    df = make_splits(genes, pos_with_extra, n_rep=1, k=5)
    assert df.gene.isin({"G9999", "FOO"}).sum() == 0


def test_every_gene_once_per_repeat():
    genes, pos = _toy_universe()
    df = make_splits(genes, pos, n_rep=4, k=5)
    grp = df.groupby(["repeat", "gene"]).size()
    assert (grp == 1).all()


def test_stratification_keeps_pos_rate_per_fold():
    genes, pos = _toy_universe(n=1000, n_pos=100)
    df = make_splits(genes, pos, n_rep=3, k=5)
    rates = df.assign(y=df.y.astype(int)).groupby(["repeat", "fold"]).y.mean()
    overall = 100 / 1000
    assert (rates - overall).abs().max() < 0.01


def test_summarise_splits_keys():
    genes, pos = _toy_universe()
    df = make_splits(genes, pos, n_rep=2, k=5)
    s = summarise_splits(df)
    assert s["n_genes"] == 200
    assert s["n_positives"] == 20
    assert s["n_repeats"] == 2 and s["n_folds"] == 5


def test_write_splits_round_trip(tmp_path):
    genes, pos = _toy_universe(n=80, n_pos=8)
    df = make_splits(genes, pos, n_rep=1, k=4)
    path = write_splits(df, "toy_label", "toy_universe", out_dir=tmp_path)
    assert path.name == "toy_label__toy_universe.tsv"
    loaded = pd.read_csv(path, sep="\t", dtype=str)
    assert list(loaded.columns) == SPLIT_COLUMNS
    assert len(loaded) == len(df)
    assert (loaded.y.astype(int).sum()) == 8


def test_default_base_seed_is_a_known_value():
    # If this changes, also update docs/decisions.md (D-BASE) and the protocol.
    assert DEFAULT_BASE_SEED == 20260926
