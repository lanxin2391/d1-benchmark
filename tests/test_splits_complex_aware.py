"""Tests for the complex-aware (group-preserving) CV split generator.

Phase D, task D2. The leakage channel being closed (D1): a held-out
positive has a complex-linked gene in the TRAINING fold, and RWR
propagates that partner's signal into the held-out gene. The fix keeps
each complex TOGETHER in one fold.
"""
from __future__ import annotations

import pandas as pd
import pytest

from d1.splits import (DEFAULT_BASE_SEED, complex_groups, make_splits,
                       make_splits_complex_aware, SPLIT_COLUMNS,
                       summarise_splits)


def _toy_universe(n=200, n_pos=20):
    genes = [f"G{i:04d}" for i in range(n)]
    pos = set(genes[:n_pos])
    return genes, pos


def _synthetic_complex_map(genes, n_complex=5, size=12, overlap=False):
    """A map with a few big complexes; optionally two overlapping ones."""
    cmap = {}
    n = len(genes)
    for c in range(n_complex):
        members = genes[c * size:(c + 1) * size]
        for g in members:
            cmap.setdefault(g, []).append(f"CPX-{c}")
    if overlap:
        members = genes[0:size]
        for g in members:
            cmap.setdefault(g, []).append("CPX-overlap")
    return cmap


# --------------------------------------------------------------------------
# 1. a complex is never split across folds
# --------------------------------------------------------------------------
def test_complex_group_is_never_split_across_folds():
    """Every complex group lies entirely in one fold, except for the
    documented repair (one positive member moved out to give an empty
    fold a positive)."""
    genes, pos = _toy_universe(n=200, n_pos=20)
    cmap = _synthetic_complex_map(genes, n_complex=5, size=12)
    exceptions: list[dict] = []
    df = make_splits_complex_aware(genes, pos, complex_map=cmap,
                                   n_rep=3, k=5, exceptions_out=exceptions)
    split_genes = {e["gene"] for e in exceptions}
    for c, members in {c: [g for g, cs in cmap.items() if c in cs]
                       for c in {x for cs in cmap.values() for x in cs}}.items():
        for r, g in df.groupby("repeat"):
            in_cpx = g[g.gene.isin(members)]
            if in_cpx["fold"].nunique() > 1:
                # the only genes that may differ are the repaired ones
                others = [x for x in members if x not in split_genes]
                kept = g[g.gene.isin(others)]
                assert kept["fold"].nunique() == 1, (
                    f"complex {c} was split across folds without a "
                    f"documented repair")


def test_overlapping_complexes_merge_into_one_group():
    genes, pos = _toy_universe(n=200, n_pos=20)
    cmap = _synthetic_complex_map(genes, n_complex=5, size=12, overlap=True)
    groups = complex_groups(genes, cmap)
    # the overlapping complex links C0 and CPX-overlap -> they merge
    big = [g for g in groups if len(g) >= 12]
    assert len(big) == 5, "overlapping complexes must merge into one group"


def test_sparsity_guard_keeps_the_largest_complexes():
    """With max_links_per_gene=1 a gene keeps only its LARGEST complex,
    so the smaller complex's links disappear."""
    genes = [f"G{i:03d}" for i in range(60)]
    cmap = {g: ["CPX-a"] for g in genes[:10]}          # 10 members
    for g in genes[10:50]:
        cmap.setdefault(g, []).append("CPX-b")         # 40 members
    for g in genes[:10]:
        cmap[g].append("CPX-b")
    for g in genes[50:]:
        cmap.setdefault(g, [])
    # CPX-b has 40 members, CPX-a has 10 -> with the guard every gene
    # keeps only CPX-b, so CPX-a's members are no longer linked to it
    groups = complex_groups(genes, cmap, max_links_per_gene=1)
    assert max(len(g) for g in groups) == 50
    assert len([g for g in groups if len(g) == 10]) == 0


# --------------------------------------------------------------------------
# 2. determinism
# --------------------------------------------------------------------------
def test_complex_aware_splits_are_deterministic():
    genes, pos = _toy_universe()
    cmap = _synthetic_complex_map(genes)
    a = make_splits_complex_aware(genes, pos, complex_map=cmap, n_rep=3, k=5)
    b = make_splits_complex_aware(genes, pos, complex_map=cmap, n_rep=3, k=5)
    pd.testing.assert_frame_equal(a, b)


def test_repeats_differ_from_each_other():
    """The repeats must be different partitions, not 40 copies of one."""
    genes, pos = _toy_universe(n=400, n_pos=120)
    cmap = _synthetic_complex_map(genes, n_complex=8, size=20)
    df = make_splits_complex_aware(genes, pos, complex_map=cmap, n_rep=4, k=5)
    assignments = df.groupby("repeat").apply(
        lambda g: dict(zip(g.gene, g.fold)), include_groups=False)
    distinct = {tuple(sorted(a.items())) for a in assignments}
    assert len(distinct) == 4, "all repeats produced the same partition"


# --------------------------------------------------------------------------
# 3./5. schema, coverage, positives per fold
# --------------------------------------------------------------------------
def test_schema_and_coverage():
    genes, pos = _toy_universe(n=200, n_pos=20)
    cmap = _synthetic_complex_map(genes)
    df = make_splits_complex_aware(genes, pos, complex_map=cmap, n_rep=3, k=5)
    assert list(df.columns) == SPLIT_COLUMNS
    assert len(df) == 3 * len(genes)
    for r in range(3):
        assert df[df.repeat == r].gene.nunique() == len(genes)
    per_fold = df.groupby(["repeat", "fold"]).y.sum()
    assert (per_fold > 0).all(), "a fold has no positives"


def test_positives_per_fold_are_balanced_on_a_realistic_panel():
    """A group-preserving split must still balance the positives well
    when the groups are small relative to the positive set."""
    genes = [f"G{i:04d}" for i in range(1000)]
    pos = set(genes[:300])
    cmap = _synthetic_complex_map(genes, n_complex=20, size=12)
    df = make_splits_complex_aware(genes, pos, complex_map=cmap, n_rep=2, k=5)
    per_fold = df.groupby(["repeat", "fold"]).y.sum()
    target = 300 / 5
    dev = (per_fold - target).abs() / target
    assert dev.max() < 0.35, f"worst deviation {dev.max():.2f}"


# --------------------------------------------------------------------------
# 4. v1 parity with an empty map
# --------------------------------------------------------------------------
def test_empty_complex_map_reproduces_v1_exactly():
    genes, pos = _toy_universe(n=300, n_pos=40)
    v1 = make_splits(genes, pos, base_seed=DEFAULT_BASE_SEED, n_rep=3, k=5)
    cx = make_splits_complex_aware(genes, pos, complex_map={},
                                    base_seed=DEFAULT_BASE_SEED, n_rep=3, k=5)
    pd.testing.assert_frame_equal(v1, cx)


def test_summarise_splits_works_on_complex_aware_output():
    genes, pos = _toy_universe()
    cmap = _synthetic_complex_map(genes)
    df = make_splits_complex_aware(genes, pos, complex_map=cmap, n_rep=2, k=5)
    s = summarise_splits(df)
    assert s["n_genes"] == len(genes)
    assert s["n_positives"] == len(pos)
    assert s["n_repeats"] == 2
