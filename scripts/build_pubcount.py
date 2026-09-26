"""B3: per-gene publication counts (for control 7 / H3).

Reads ``gene2pubmed.gz`` (NCBI gene2pubmed), keeps only human rows
(tax_id == 9606), counts unique PubMed IDs per GeneID, maps GeneID -> HGNC
symbol, and writes ``data/processed/reference/pubcount.tsv``.

Output schema (handed to Person A's control-7 pipeline later):
    gene      HGNC approved symbol
    n_pubmed  integer count
    log10     log10(n_pubmed + 1)   (the -1 keeps n=0 genes at 0, not log-shifted)

Decision D-LOG10: use log10(n + 1) (not log10(n)) so genes with no publications
are not artificially given a negative value; this matches what the protocol's
control-7 model regresses out.

Notes
-----
* Reading 287 MB compressed (~13 GB raw, tens of millions of rows) is done with
  pandas' chunked CSV reader to keep memory flat at a few hundred MB.
* ``keep_default_na=False`` so that empty / ``-`` tokens stay as empty strings.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW = REPO_ROOT / "data" / "raw" / "reference" / "gene2pubmed.gz"
HGNC = REPO_ROOT / "data" / "raw" / "reference" / "hgnc_complete_set.txt"
OUT = REPO_ROOT / "data" / "processed" / "reference" / "pubcount.tsv"
OUT.parent.mkdir(parents=True, exist_ok=True)

HUMAN_TAX_ID = "9606"
CHUNK_ROWS = 200_000


def _read_human_gene2pubmed(path: Path) -> pd.DataFrame:
    """Return a Series ``{GeneID: set_of_PubMed_IDs}`` over human rows only.

    We do not materialise the whole file; instead we build a dict of
    ``{GeneID: set(PubMed_ID)}`` chunk-by-chunk to keep memory bounded by the
    number of distinct human GeneIDs (~20k protein-coding + ~30k non-coding).
    """
    counts: dict[str, set[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as f:
        # Skip the "#tax_id" header
        header = f.readline()
        if not header.startswith("#"):
            raise ValueError(f"unexpected header: {header!r}")
        for chunk in pd.read_csv(
            f, sep="\t", header=None, names=["tax_id", "gene_id", "pmid"],
            dtype=str, chunksize=CHUNK_ROWS,
        ):
            chunk = chunk[chunk["tax_id"] == HUMAN_TAX_ID]
            if chunk.empty:
                continue
            for gid, pmid in zip(chunk["gene_id"].to_numpy(),
                                 chunk["pmid"].to_numpy()):
                counts.setdefault(gid, set()).add(pmid)
    return counts


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hgnc", type=Path, default=HGNC)
    ap.add_argument("--raw", type=Path, default=RAW)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    if not args.raw.exists():
        print(f"!! raw gene2pubmed.gz not found at {args.raw}", file=sys.stderr)
        return 1
    if not args.hgnc.exists():
        print(f"!! HGNC complete set not found at {args.hgnc}", file=sys.stderr)
        return 1

    print(f"Reading {args.raw.name} (this takes ~1-2 min)...")
    counts = _read_human_gene2pubmed(args.raw)
    print(f"  distinct human GeneIDs: {len(counts):,}")

    # Map GeneID -> HGNC symbol via the HGNC mapper (id_type='entrez').
    sys.path.insert(0, str(REPO_ROOT))
    from d1.hgnc import HGNCMapper
    mapper = HGNCMapper(args.hgnc)
    gene_ids = pd.Series(list(counts.keys()))
    symbols = mapper.map(gene_ids, "entrez")
    n_pubmed = pd.Series([len(counts[g]) for g in gene_ids], dtype=int)

    df = pd.DataFrame({"gene": symbols, "n_pubmed": n_pubmed})
    df = df.dropna(subset=["gene"])
    df = df.drop_duplicates(subset=["gene"], keep="first")     # multiple GeneIDs may collapse onto the same symbol
    df = df.sort_values("gene", kind="mergesort").reset_index(drop=True)
    df["log10"] = np.log10(df["n_pubmed"].to_numpy() + 1)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"\nWrote {args.out} ({len(df):,} approved-symbol rows)")
    print(f"  median publications: {int(df['n_pubmed'].median())}")
    print(f"  top 10 by publications:")
    for _, row in df.sort_values("n_pubmed", ascending=False).head(10).iterrows():
        print(f"    {row['gene']:10s}  {row['n_pubmed']:>6d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
