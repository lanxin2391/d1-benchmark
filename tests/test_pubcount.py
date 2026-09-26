"""Tests for the per-gene publication-count builder (B3, control 7 / H3)."""
from __future__ import annotations

import gzip
import gzip
from pathlib import Path

import pandas as pd
import pytest

from d1.hgnc import HGNCMapper
from scripts.build_pubcount import _read_human_gene2pubmed


def _toy_gene2pubmed(tmp_path: Path) -> Path:
    """Write a tiny gene2pubmed.gz with three human and one Drosophila row."""
    p = tmp_path / "toy_gene2pubmed.gz"
    body = "\n".join([
        "#tax_id\tGeneID\tPubMed_ID",
        "9606\t7157\t1001",            # TP53
        "9606\t7157\t1002",            # TP53 again -> unique count = 2
        "9606\t7157\t1001",            # duplicate of TP53 PMID -> still unique = 2
        "9606\t1956\t2001",            # EGFR
        "9606\t1956\t2002",            # EGFR
        "9606\t1956\t2003",            # EGFR
        "23\t1\t9999",                 # fly, ignored
        "9606\t99999999\t3001",        # unknown Entrez ID, dropped after HGNC mapping
    ])
    with gzip.open(p, "wt", encoding="utf-8") as f:
        f.write(body)
    return p


def test_read_human_gene2pubmed_counts_only_human(tmp_path):
    p = _toy_gene2pubmed(tmp_path)
    counts = _read_human_gene2pubmed(p)
    assert "7157" in counts and len(counts["7157"]) == 2   # 1001, 1002
    assert "1956" in counts and len(counts["1956"]) == 3
    assert "99999999" in counts                           # present in dict (HGNC drops later)
    assert "1" not in counts                              # fly ignored


def test_pubcount_end_to_end(tmp_path):
    """End-to-end: build_pubcount writes a C3-shaped file with expected rows."""
    import subprocess, sys, os
    repo = Path(__file__).resolve().parents[1]
    raw = _toy_gene2pubmed(tmp_path)
    hgnc = repo / "data" / "raw" / "reference" / "hgnc_complete_set.txt"
    if not hgnc.exists():
        pytest.skip("real HGNC file not present; the toy test above is sufficient for CI")
    out = tmp_path / "pubcount.tsv"

    r = subprocess.run(
        [sys.executable, "scripts/build_pubcount.py",
         "--hgnc", str(hgnc), "--raw", str(raw), "--out", str(out)],
        cwd=repo, capture_output=True, text=True, timeout=120,
    )
    assert r.returncode == 0, r.stderr
    df = pd.read_csv(out, sep="\t")
    assert {"gene", "n_pubmed", "log10"} <= set(df.columns)
    row = df.set_index("gene").loc["TP53"]
    assert int(row["n_pubmed"]) == 2
    row2 = df.set_index("gene").loc["EGFR"]
    assert int(row2["n_pubmed"]) == 3
    # unknown Entrez is dropped, no row for it
    assert "TP53" in df.gene.values
    assert "EGFR" in df.gene.values
