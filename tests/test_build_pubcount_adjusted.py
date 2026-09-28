"""Tests for ``scripts/build_pubcount_adjusted.py`` (W3-B Control 7).

Covers:
  * the 4 canonical adjusted labels are generated
  * C3 schema
  * the trim reduces the size to exactly 50 % of positives
  * all adjusted-label genes are in the original positives
  * β is positive (cancer drivers are better studied than non-drivers)
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))


_spec = importlib.util.spec_from_file_location(
    "build_pubcount_adjusted",
    REPO_ROOT / "scripts" / "build_pubcount_adjusted.py",
)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


CANONICAL = ("intogen2024", "ot_all", "ot_nolit", "clingen_label")
LABELS_DIR = REPO_ROOT / "data" / "processed" / "labels"
PUBCOUNT_PATH = LABELS_DIR / "pubcount.tsv"


def _adjusted_path(label: str) -> Path:
    return LABELS_DIR / f"{label}_pubcount_adjusted.tsv"


# ------------------------------------------------------------------ schema
@pytest.mark.parametrize("label", CANONICAL)
def test_adjusted_file_exists(label: str) -> None:
    p = _adjusted_path(label)
    assert p.exists(), f"missing adjusted file: {p}"


@pytest.mark.parametrize("label", CANONICAL)
def test_adjusted_schema(label: str) -> None:
    df = pd.read_csv(_adjusted_path(label), sep="\t")
    expected = {"gene", "rank", "source", "release"}
    assert expected.issubset(set(df.columns)), \
        f"missing columns: {expected - set(df.columns)}"


@pytest.mark.parametrize("label", CANONICAL)
def test_adjusted_subset_of_original(label: str) -> None:
    """Every adjusted-label gene must appear in the original positives."""
    adj = pd.read_csv(_adjusted_path(label), sep="\t")
    orig = pd.read_csv(LABELS_DIR / f"{label}.tsv", sep="\t")
    assert set(adj.gene).issubset(set(orig.gene))


@pytest.mark.parametrize("label", CANONICAL)
def test_trim_is_50_percent(label: str) -> None:
    """At q=0.5 we keep ~50 % of the positives (the implementation uses
    a `>=` threshold against the median quantile, which keeps a few
    extra genes when the median value is tied). Allow a ±5 % window."""
    adj = pd.read_csv(_adjusted_path(label), sep="\t")
    orig = pd.read_csv(LABELS_DIR / f"{label}.tsv", sep="\t")
    expected_low = int(len(orig) * 0.45)
    expected_high = int(len(orig) * 0.55) + 1
    assert expected_low <= len(adj) <= expected_high, \
        f"{label}: kept {len(adj)}, expected 45-55% of {len(orig)} = [{expected_low}, {expected_high}]"


@pytest.mark.parametrize("label", CANONICAL)
def test_rank_is_dense_and_starts_at_1(label: str) -> None:
    df = pd.read_csv(_adjusted_path(label), sep="\t")
    assert df["rank"].iloc[0] == 1
    assert df["rank"].is_monotonic_increasing
    assert set(df["rank"]) == set(range(1, len(df) + 1))


# --------------------------------------------------------------- regression
def test_beta_is_positive() -> None:
    """β = pos_mean - neg_mean should be positive (cancer drivers are
    over-studied). We check the index file written by build_pubcount_adjusted.
    """
    idx_path = LABELS_DIR / "pubcount_adjusted_index.tsv"
    assert idx_path.exists(), "missing pubcount_adjusted_index.tsv"
    idx = pd.read_csv(idx_path, sep="\t")
    assert (idx.beta > 0).all(), "all β should be positive"


def test_smoke_build_one_label() -> None:
    """Direct call to build_adjusted works on intogen2024."""
    label_path = LABELS_DIR / "intogen2024.tsv"
    out_path = LABELS_DIR / "intogen2024_pubcount_adjusted_REDO.tsv"
    s = _mod.build_adjusted("intogen2024", PUBCOUNT_PATH, label_path,
                            out_path, residual_top_quantile=0.5)
    # 45-55% of 633 should pass
    assert 285 <= s["n_kept"] <= 350
    assert s["n_positives_total"] == 633
    assert s["beta"] > 0
    out_path.unlink(missing_ok=True)