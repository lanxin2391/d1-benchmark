"""Tests for the label-network overlap builder (B6).

Uses an in-memory fake network panel and the four label files produced by
the other B-side builders, so the test is independent of A's actual network
files landing in data/processed/networks/.
"""
from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path

import pandas as pd
import pytest

# Locate scripts/label_network_overlap.py and import it as a module.
SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "label_network_overlap.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("label_network_overlap", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture()
def fake_panel(tmp_path, monkeypatch):
    """Write a tiny network panel of two networks and patch the SCRIPT's paths."""
    net_dir = tmp_path / "data" / "processed" / "networks"
    net_dir.mkdir(parents=True)
    label_dir = tmp_path / "data" / "processed" / "labels"
    label_dir.mkdir(parents=True)
    results_dir = tmp_path / "results" / "tables"
    results_dir.mkdir(parents=True)

    # Network A: {TP53, KRAS, BRCA1} as a triangle
    pd.DataFrame({
        "gene_a": ["TP53", "TP53", "KRAS"],
        "gene_b": ["KRAS", "BRCA1", "BRCA1"],
        "weight": [1, 1, 1],
    }).to_csv(net_dir / "tiny_a.tsv", sep="\t", index=False)
    # Network B: {TP53, EGFR} only
    pd.DataFrame({
        "gene_a": ["TP53"],
        "gene_b": ["EGFR"],
        "weight": [1],
    }).to_csv(net_dir / "tiny_b.tsv", sep="\t", index=False)
    # Shared universe: {TP53, KRAS, BRCA1, EGFR}
    (net_dir / "universe_shared.txt").write_text("TP53\nKRAS\nBRCA1\nEGFR\n")

    # Two label sets
    pd.DataFrame({
        "gene": ["TP53", "KRAS", "BRCA1", "EGFR", "FAKE"],
        "rank": [1, 2, 3, 4, 5],
        "source": "toy",
        "release": "test",
    }).to_csv(label_dir / "label_a.tsv", sep="\t", index=False)
    pd.DataFrame({
        "gene": ["TP53", "EGFR"],
        "rank": [1, 2],
        "source": "toy",
        "release": "test",
    }).to_csv(label_dir / "label_b.tsv", sep="\t", index=False)

    # Patch SCRIPT constants to point at our temp dirs
    mod = _load_module()
    monkeypatch.setattr(mod, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(mod, "LABELS_DIR", label_dir)
    monkeypatch.setattr(mod, "NETWORKS_DIR", net_dir)
    monkeypatch.setattr(mod, "RESULTS_DIR", results_dir)
    monkeypatch.setattr(mod, "LABEL_FILES", {"label_a": "label_a.tsv", "label_b": "label_b.tsv"})
    return tmp_path, mod


def test_overlap_table_shape(fake_panel):
    _, mod = fake_panel
    rc = mod.main([])
    assert rc == 0
    out = Path(mod.RESULTS_DIR) / "label_network_overlap.tsv"
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 2 * 2                                  # 2 labels x 2 networks
    assert set(df.columns) >= {"label_set", "network", "n_overlap",
                                "n_label_positives", "n_network_lcc",
                                "frac_label_in_network"}
    pivot = df.pivot(index="label_set", columns="network", values="n_overlap")
    assert pivot.loc["label_a", "tiny_a"] == 3                # TP53, KRAS, BRCA1
    assert pivot.loc["label_a", "tiny_b"] == 2                # TP53, EGFR
    assert pivot.loc["label_b", "tiny_a"] == 1                # TP53
    assert pivot.loc["label_b", "tiny_b"] == 2                # TP53, EGFR


def test_common_cap_uses_shared_universe(fake_panel):
    _, mod = fake_panel
    mod.main([])
    # label_a in shared = {TP53, KRAS, BRCA1, EGFR, FAKE} ∩ shared = 4
    # label_b in shared = {TP53, EGFR} ∩ shared = 2
    # N = min = 2
    cap_a = pd.read_csv(Path(mod.LABELS_DIR) / "label_a.capped.tsv", sep="\t")
    cap_b = pd.read_csv(Path(mod.LABELS_DIR) / "label_b.capped.tsv", sep="\t")
    assert len(cap_a) == 2
    assert len(cap_b) == 2
    assert set(cap_b.gene) == {"TP53", "EGFR"}


def test_missing_network_dir_exits_one(tmp_path, monkeypatch):
    """If no networks are found under data/processed/networks, main() exits 1."""
    net_dir = tmp_path / "empty" / "networks"
    net_dir.mkdir(parents=True)
    mod = _load_module()
    monkeypatch.setattr(mod, "NETWORKS_DIR", net_dir)
    rc = mod.main([])
    assert rc == 1
