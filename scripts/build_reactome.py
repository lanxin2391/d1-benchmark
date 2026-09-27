"""Build the Reactome network (A3 task).

Reads:
    data/raw/networks/reactome/reactome.homo_sapiens.interactions.tab-delimited.txt
        tab-separated, 9 columns, header line starting with "#".
        Columns (after #):
          0: Interactor 1 UniProt id          (e.g. "uniprotkb:Q9Y287")
          1: Interactor 1 Ensembl gene id     (pipe-separated, e.g. "ENSG00000136156|...")
          2: Interactor 1 Entrez Gene id      ("-" if absent)
          3: Interactor 2 UniProt id
          4: Interactor 2 Ensembl gene id
          5: Interactor 2 Entrez Gene id
          6: Interaction type                  ("physical association", ...)
          7: Interaction context                (e.g. "reactome:R-HSA-976871")
          8: Pubmed references                  (pipe-separated)

Writes:
    data/processed/networks/reactome.tsv        (contract C2)
    data/processed/networks/reactome.meta.json  (contract C2b)

Decisions referenced:
    D-04 human-only is implicit (file is named homo_sapiens.*).
    D-02 HGNC symbol mapping via Person B's d1.hgnc.HGNCMapper (UniProt -> HGNC).

Notes
-----
* Reactome has many self-interactions (a protein interacting with itself via
  different domains). finalize_edges drops these in the self-loop step.
* Many columns are missing Ensembl/Entrez entries ("-"). We use UniProt (col 0
  and col 3) because it is always present and Person B's HGNCMapper handles
  the mapping robustly.
"""
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from d1.networks.io import finalize_edges, check_network_file  # noqa: E402
from d1.hgnc import HGNCMapper  # noqa: E402

RAW_PATH = os.path.join(ROOT, "data", "raw", "networks", "reactome",
                        "reactome.homo_sapiens.interactions.tab-delimited.txt")
MAPPING_PATH = os.path.join(ROOT, "data", "raw", "reference", "hgnc_complete_set.txt")
OUT_DIR = os.path.join(ROOT, "data", "processed", "networks")
NET_ID = "reactome"


def _strip_prefix(s: str) -> str:
    if s.startswith("uniprotkb:"):
        return s[len("uniprotkb:"):]
    return s


def main():
    print("Reactome network (A3)")
    print("Step 1: loading HGNC mapper (UniProt -> HGNC symbol)")
    mapper = HGNCMapper(MAPPING_PATH)
    print(f"  mapper covers {len(mapper.maps['uniprot']):,} UniProt -> HGNC entries")

    print("\nStep 2: reading Reactome interactions")
    df = pd.read_csv(RAW_PATH, sep="\t", comment="#",
                     names=["u1", "e1", "g1", "u2", "e2", "g2", "itype", "icontext", "pubs"],
                     dtype=str, keep_default_na=False, na_values=[""])
    print(f"  raw rows: {len(df):,}")
    df = df.assign(gene_a=df["u1"].map(_strip_prefix), gene_b=df["u2"].map(_strip_prefix))
    df = df[["gene_a", "gene_b", "itype"]]

    print("\nStep 3: map UniProt -> HGNC via Person B's HGNCMapper")
    before = len(df)
    df = df.assign(
        gene_a=mapper.map(df["gene_a"], "uniprot"),
        gene_b=mapper.map(df["gene_b"], "uniprot"),
    )
    df = df.dropna(subset=["gene_a", "gene_b"])
    print(f"  after mapping: {len(df):,}/{before:,} ({100*len(df)/before:.1f}%)")
    itypes = df["itype"].value_counts().head(5).to_dict()
    print(f"  top interaction types: {itypes}")

    print("\nStep 4: finalize_edges (drop self-loops, dedup, LCC)")
    out_df, meta = finalize_edges(
        df[["gene_a", "gene_b"]], net_id=NET_ID,
        meta={"source_file": "reactome.homo_sapiens.interactions.tab-delimited.txt",
              "threshold": None,
              "description": "Reactome curated pathway interactions (human only)"},
        out_dir=OUT_DIR,
    )
    print(f"  edges: {len(out_df):,}  LCC nodes: {meta['n_nodes_lcc']:,}")
    print(f"  self-loops removed: {meta['n_self_loops_removed']}")
    print(f"  duplicates merged: {meta['n_edges_dedup'] - len(out_df)}")

    check_network_file(os.path.join(OUT_DIR, f"{NET_ID}.tsv"))
    print(f"  OK: {NET_ID}.tsv passes contract C2")


if __name__ == "__main__":
    main()