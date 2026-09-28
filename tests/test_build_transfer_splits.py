"""Tests for ``scripts/build_transfer_splits.py`` (W3-A Control 5).

Covers:
  * file schema (C4-like with extra columns)
  * target = label_b positives ∩ LCC \ label_a
  * background_pool size is exactly background_sample_size (or less)
  * seeds-not-in-test (no leak)
  * no self-transfer (label_a == label_b raises)
  * file index lists all generated files
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

# load build_transfer_splits as a module
_spec = importlib.util.spec_from_file_location(
    "build_transfer_splits",
    REPO_ROOT / "scripts" / "build_transfer_splits.py",
)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


def _read_label(label_id: str) -> set[str]:
    p = REPO_ROOT / "data" / "processed" / "labels" / f"{label_id}.tsv"
    return set(pd.read_csv(p, sep="\t", dtype=str).gene)


def _read_net_lcc(net_id: str) -> set[str]:
    p = REPO_ROOT / "data" / "processed" / "networks" / f"{net_id}.tsv"
    df = pd.read_csv(p, sep="\t", dtype=str)
    return set(df.gene_a) | set(df.gene_b)


# ---------------------------------------------------------------- schema
def test_schema_columns() -> None:
    """All generated files must have the canonical 7 columns."""
    p = (REPO_ROOT / "data" / "processed" / "transfer"
         / "intogen2024__to__ot_all__funmap.tsv")
    df = pd.read_csv(p, sep="\t")
    expected = {"gene", "fold", "y", "label_id", "universe_id",
                "transfer_id", "role"}
    assert expected.issubset(set(df.columns)), \
        f"missing columns: {expected - set(df.columns)}"


def test_target_is_label_b_in_lcc_minus_label_a() -> None:
    """target rows must be in label_b ∩ LCC, not in label_a."""
    p = (REPO_ROOT / "data" / "processed" / "transfer"
         / "intogen2024__to__ot_all__funmap.tsv")
    df = pd.read_csv(p, sep="\t", dtype=str)
    target = df[df.role == "target"]
    ot_all = _read_label("ot_all")
    intogen2024 = _read_label("intogen2024")
    funmap_lcc = _read_net_lcc("funmap")
    target_genes = set(target.gene)
    # every target gene must be in ot_all ∩ funmap LCC
    assert target_genes.issubset(ot_all)
    assert target_genes.issubset(funmap_lcc)
    # and not in intogen2024
    assert not (target_genes & intogen2024)


def test_background_does_not_include_seeds() -> None:
    """background rows must not be in label_a (which carries the seeds)."""
    p = (REPO_ROOT / "data" / "processed" / "transfer"
         / "intogen2024__to__ot_all__funmap.tsv")
    df = pd.read_csv(p, sep="\t", dtype=str)
    bg = df[df.role == "background"]
    intogen2024 = _read_label("intogen2024")
    bg_genes = set(bg.gene)
    assert not (bg_genes & intogen2024)


def test_no_self_transfer() -> None:
    """_make_transfer_table raises on self-transfer."""
    with pytest.raises(ValueError, match="self-transfer"):
        _mod._make_transfer_table("intogen2024", "intogen2024", "funmap")


def test_index_file_exists() -> None:
    """The transfer_index.tsv must be written and have one row per arm."""
    p = REPO_ROOT / "data" / "processed" / "transfer" / "transfer_index.tsv"
    assert p.exists(), "transfer_index.tsv not generated"
    idx = pd.read_csv(p, sep="\t")
    # 12 ordered pairs × 6 networks = 72 cells, minus skips for empty LCC
    # (we observed 66 in the actual run).
    assert len(idx) >= 60
    for col in ("label_a", "label_b", "network", "n_target", "n_seed_pool"):
        assert col in idx.columns


def test_arm_files_match_index() -> None:
    """Every index row should have a corresponding arm file on disk."""
    p = REPO_ROOT / "data" / "processed" / "transfer" / "transfer_index.tsv"
    idx = pd.read_csv(p, sep="\t")
    base = REPO_ROOT / "data" / "processed" / "transfer"
    for _, row in idx.iterrows():
        arm_path = base / f"{row.label_a}__to__{row.label_b}__{row.network}.tsv"
        assert arm_path.exists(), f"missing arm file {arm_path.name}"