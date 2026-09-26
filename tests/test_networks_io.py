"""Tests for the shared edge cleaner and network reader (contract C2)."""
import json

import numpy as np
import pandas as pd
import pytest

from d1.networks.io import check_network_file, finalize_edges, load_network


def messy_edges():
    return pd.DataFrame({
        "gene_a": ["TP53", "MDM2", "TP53", "EGFR", "KRAS", "BRAF", None, "X1"],
        "gene_b": ["MDM2", "TP53", "TP53", "KRAS", "EGFR", "KRAS", "TP53", "X2"],
        "weight": [0.9, 0.8, 1.0, 0.5, 0.7, 0.6, 0.3, 0.2],
    })
    # TP53-MDM2 twice (reversed), a self-loop, EGFR-KRAS twice, an unmapped gene,
    # and a separate small component X1-X2; plus a bridge-less component {EGFR,KRAS,BRAF}


def test_finalize_edges_cleans_everything(tmp_path):
    df, meta = finalize_edges(messy_edges(), "toy", out_dir=tmp_path)
    # components: {TP53, MDM2}, {EGFR, KRAS, BRAF}, {X1, X2}; LCC = EGFR-KRAS-BRAF
    assert set(df.gene_a) | set(df.gene_b) == {"BRAF", "EGFR", "KRAS"}
    assert (df.gene_a < df.gene_b).all()
    assert not df.duplicated(["gene_a", "gene_b"]).any()
    assert df.loc[(df.gene_a == "EGFR") & (df.gene_b == "KRAS"), "weight"].item() == 0.7
    assert meta["n_input_edges"] == 8
    assert meta["n_mapped_edges"] == 7
    assert meta["n_self_loops_removed"] == 1
    assert meta["n_components"] == 3
    assert meta["n_nodes_lcc"] == 3 and meta["n_edges_lcc"] == 2
    assert json.loads((tmp_path / "toy.meta.json").read_text())["net_id"] == "toy"
    check_network_file(tmp_path / "toy.tsv")


def test_load_network_symmetric_binary(tmp_path):
    finalize_edges(messy_edges(), "toy", out_dir=tmp_path)
    genes, A = load_network(tmp_path / "toy.tsv")
    assert list(genes) == ["BRAF", "EGFR", "KRAS"]
    assert (A != A.T).nnz == 0
    assert set(np.unique(A.data)) == {1.0}
    assert A.sum() == 4                      # 2 edges, stored in both directions


def test_check_network_file_catches_bad_rows(tmp_path):
    p = tmp_path / "bad.tsv"
    pd.DataFrame({"gene_a": ["B", "A"], "gene_b": ["A", "B"], "weight": [1, 1]}).to_csv(
        p, sep="\t", index=False)
    with pytest.raises(ValueError, match="contract C2"):
        check_network_file(p)
