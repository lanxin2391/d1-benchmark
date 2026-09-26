"""Build STRING networks (A2 task).

Reads:
    data/raw/networks/string/9606.protein.info.v12.0.txt.gz
        ENSP -> gene symbol mapping (preferred_name column).
    data/raw/networks/string/9606.protein.physical.links.v12.0.txt.gz
        physical-only interaction subnetwork, 3 columns: protein1, protein2, combined_score
    data/raw/networks/string/9606.protein.links.detailed.v12.0.txt.gz
        full network with per-channel scores, combined_score is the aggregate

Writes:
    data/processed/networks/string_phys700.tsv + meta.json
    data/processed/networks/string_full400.tsv + meta.json
    data/processed/networks/string_full700.tsv + meta.json
    data/processed/networks/string_full900.tsv + meta.json

Decision references:
    D-03 STRING physical >= 700 -> string_phys700
    D-01 unweighted: we keep the score column for the threshold cut, then drop it
    D-02 HGNC symbol mapping via STRING's preferred_name (sometimes a protein has no
         symbol, e.g. an uncharacterized ENSP; these edges are dropped silently)

Why three full-network thresholds (400/700/900):
    STRING combined_score is a Bayesian integration of seven evidence channels; the
    conventional "medium" cut is 400, the "high" cut is 700, and 900 is "highest".
    We pre-build all three so we can measure sensitivity of ΔAUROC to threshold
    choice. 400 gives ~6x more edges than 900; 700 is the protocol's primary cut.
"""
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from d1.networks.io import finalize_edges, check_network_file  # noqa: E402

RAW_DIR = os.path.join(ROOT, "data", "raw", "networks", "string")
OUT_DIR = os.path.join(ROOT, "data", "processed", "networks")
MAPPING_PATH = os.path.join(RAW_DIR, "9606.protein.info.v12.0.txt.gz")
PHYSICAL_PATH = os.path.join(RAW_DIR, "9606.protein.physical.links.v12.0.txt.gz")
DETAILED_PATH = os.path.join(RAW_DIR, "9606.protein.links.detailed.v12.0.txt.gz")


def load_mapping():
    """Return {ENSP: gene_symbol} (with the 9606. prefix stripped)."""
    info = pd.read_csv(MAPPING_PATH, sep="\t")
    info["ensp"] = info["#string_protein_id"].str.replace("9606.", "", regex=False)
    return dict(zip(info["ensp"], info["preferred_name"]))


def load_and_map(src_path):
    """Read a STRING link file and return a gene_a / gene_b / weight frame in HGNC symbols."""
    df = pd.read_csv(src_path, sep=" ")
    df["ensp1"] = df["protein1"].str.replace("9606.", "", regex=False)
    df["ensp2"] = df["protein2"].str.replace("9606.", "", regex=False)
    mapping = load_mapping()
    df["gene_a"] = df["ensp1"].map(mapping)
    df["gene_b"] = df["ensp2"].map(mapping)
    # Keep only edges where both ends mapped to an HGNC symbol
    before = len(df)
    df = df.dropna(subset=["gene_a", "gene_b"])
    mapped = len(df)
    print(f"  ENSP -> symbol: kept {mapped}/{before} ({100*mapped/before:.1f}%)")
    return df.rename(columns={"combined_score": "weight"})[["gene_a", "gene_b", "weight"]]


def build(net_id, src_path, threshold):
    print(f"\n--- {net_id} (threshold={threshold}) ---")
    df = load_and_map(src_path)
    before = len(df)
    df = df[df["weight"] >= threshold]
    print(f"  after threshold cut: {len(df)}/{before} edges")
    out_df, meta = finalize_edges(
        df, net_id=net_id,
        meta={"source_file": os.path.basename(src_path), "threshold": threshold,
              "score_column": "combined_score"},
        out_dir=OUT_DIR,
    )
    print(f"  after finalize_edges: {len(out_df)} edges, {meta['n_nodes_lcc']} nodes in LCC")
    check_network_file(os.path.join(OUT_DIR, f"{net_id}.tsv"))
    print(f"  OK: {net_id}.tsv passes contract C2")


def main():
    print("STRING networks (A2)")
    build("string_phys700", PHYSICAL_PATH, threshold=700)
    build("string_full400", DETAILED_PATH, threshold=400)
    build("string_full700", DETAILED_PATH, threshold=700)
    build("string_full900", DETAILED_PATH, threshold=900)


if __name__ == "__main__":
    main()