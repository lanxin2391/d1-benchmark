"""Tests for the end-to-end split-file builder (B5).

Uses an in-memory fake label set + a fake network panel written via
A's d1.networks.io.finalize_edges, so the test is independent of the real
funmap / STRING / IntAct downloads.
"""
from __future__ import annotations

import importlib.util
import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_splits.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("build_splits", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture()
def fake_layout(tmp_path, monkeypatch):
    """Set up two tiny networks via finalize_edges + two tiny label files."""
    net_dir = tmp_path / "data" / "processed" / "networks"
    label_dir = tmp_path / "data" / "processed" / "labels"
    splits_dir = tmp_path / "data" / "processed" / "splits"
    for d in (net_dir, label_dir, splits_dir):
        d.mkdir(parents=True, exist_ok=True)

    sys_mod_path = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(sys_mod_path))
    from d1.networks.io import finalize_edges

    # Two toy networks sharing most of their nodes.
    # 30 nodes is enough for a 5-fold StratifiedKFold; all gene names appear
    # in both networks' LCCs to avoid an "empty positive fold" failure.
    genes = [f"G_{i:02d}" for i in range(28)]
    # toy_a = 32 nodes: 28 G_*, TP53, KRAS, BRCA1, EGFR; edges form a cycle
    a_edges_a = genes + ["TP53", "KRAS", "BRCA1", "EGFR"]
    # toy_b = 30 nodes: 28 G_* + TP53 + EGFR; edges form a cycle
    b_edges_a = genes + ["TP53", "EGFR"]
    a_edges_b = [a_edges_a[(i + 1) % len(a_edges_a)] for i in range(len(a_edges_a))]
    b_edges_b = [b_edges_a[(i + 1) % len(b_edges_a)] for i in range(len(b_edges_a))]
    pd.DataFrame({
        "gene_a": a_edges_a,
        "gene_b": a_edges_b,
        "weight": [1] * len(a_edges_a),
    }).to_csv(tmp_path / "toy_a_edges.tsv", sep="\t", index=False)
    df_a = pd.read_csv(tmp_path / "toy_a_edges.tsv", sep="\t")
    finalize_edges(df_a, "toy_a", out_dir=net_dir)

    pd.DataFrame({
        "gene_a": b_edges_a,
        "gene_b": b_edges_b,
        "weight": [1] * len(b_edges_a),
    }).to_csv(tmp_path / "toy_b_edges.tsv", sep="\t", index=False)
    df_b = pd.read_csv(tmp_path / "toy_b_edges.tsv", sep="\t")
    finalize_edges(df_b, "toy_b", out_dir=net_dir)

    # Label: pick positives that are in both LCCs
    pd.DataFrame({
        "gene": ["TP53", "KRAS", "BRCA1", "EGFR", "G_00", "G_05", "G_10", "G_15"],
        "rank": [1, 2, 3, 4, 5, 6, 7, 8],
        "source": "toy",
        "release": "test",
    }).to_csv(label_dir / "toy_label.tsv", sep="\t", index=False)

    # shared universe = intersection of toy_a LCC + toy_b LCC
    def _lcc(net_id):
        df = pd.read_csv(net_dir / f"{net_id}.tsv", sep="\t", dtype=str)
        return set(df.gene_a) | set(df.gene_b)
    shared = sorted(_lcc("toy_a") & _lcc("toy_b"))
    (net_dir / "universe_shared.txt").write_text("\n".join(shared) + "\n")

    mod = _load_module()
    monkeypatch.setattr(mod, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(mod, "LABELS_DIR", label_dir)
    monkeypatch.setattr(mod, "NETWORKS_DIR", net_dir)
    monkeypatch.setattr(mod, "SPLITS_DIR", splits_dir)
    return tmp_path, mod


def test_split_file_written_for_native_universe(fake_layout):
    _, mod = fake_layout
    rc = mod.main(["--label", "toy_label", "--universe", "native_toy_a",
                   "--n-rep", "2"])                     # 2 reps for speed
    assert rc == 0
    out = mod.SPLITS_DIR / "toy_label__native_toy_a.tsv"
    assert out.exists()
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == ["gene", "repeat", "fold", "y"]
    # rows = universe_size * n_rep (each gene appears in exactly one fold per repeat,
    #                              fold size = universe_size / k = 32/5 ≈ 6)
    n_universe = len(pd.read_csv(mod.NETWORKS_DIR / "toy_a.tsv", sep="\t"))
    assert len(df) == n_universe * 2
    assert df.groupby("repeat").gene.nunique().min() == n_universe      # every gene present
    assert df.groupby("repeat").fold.nunique().min() == 5                # 5 folds
    assert set(df.y.unique()) <= {0, 1}


def test_split_file_for_shared_universe(fake_layout):
    _, mod = fake_layout
    rc = mod.main(["--label", "toy_label", "--universe", "shared",
                   "--n-rep", "2"])
    assert rc == 0
    out = mod.SPLITS_DIR / "toy_label__shared.tsv"
    assert out.exists()
    df = pd.read_csv(out, sep="\t")
    shared = (mod.NETWORKS_DIR / "universe_shared.txt").read_text().split()
    assert len(df) == len(shared) * 2
    assert set(df.gene.unique()) == set(shared)
    assert df.groupby("repeat").fold.nunique().min() == 5


def test_unknown_universe_spec_is_skipped(fake_layout):
    _, mod = fake_layout
    rc = mod.main(["--label", "toy_label", "--universe", "native_does_not_exist"])
    # nothing written, returns 1
    assert rc == 1


def test_missing_label_exits_one(tmp_path, monkeypatch):
    mod = _load_module()
    monkeypatch.setattr(mod, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(mod, "LABELS_DIR", tmp_path)
    monkeypatch.setattr(mod, "NETWORKS_DIR", tmp_path)
    monkeypatch.setattr(mod, "SPLITS_DIR", tmp_path)
    rc = mod.main(["--label", "missing", "--universe", "shared"])
    assert rc == 1
