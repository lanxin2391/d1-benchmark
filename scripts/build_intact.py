"""Build IntAct networks (A3 task).

Reads:
    data/raw/networks/intact/intact.zip
        PSI-MITAB 2.7, 42 columns, all species, ~1.3 GB.
        The main file inside is ``intact.txt`` (~11 GB uncompressed).
Writes:
    data/processed/networks/intact.tsv        (contract C2)
    data/processed/networks/intact.meta.json  (contract C2b)
    data/processed/networks/intact_miscore045.tsv  (sensitivity arm)

Decisions referenced:
    D-04 human-human pairs (both taxid 9606); UniProt on both sides;
        isoform suffix dropped; primary arm has no MI threshold;
        MI-score >= 0.45 is the sensitivity arm.
    D-02 HGNC symbol mapping via Person B's d1.hgnc.HGNCMapper.

Why a sensitivity arm:
    IntAct provides a confidence score (`intact-miscore` between 0 and 1, computed
    from the number and type of supporting experiments). Cutting at 0.45 keeps
    ~50% of edges while removing low-confidence ones; lets us measure how
    sensitive the result is to threshold choice.

Why we filter taxid at parse time (not inside finalize_edges):
    finalize_edges does not know about taxid and would happily include
    yeast-mouse pairs. The biological filter has to happen before
    network-building.
"""
import os
import sys
import zipfile

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from d1.networks.io import finalize_edges, check_network_file  # noqa: E402
from d1.hgnc import HGNCMapper  # noqa: E402

RAW_ZIP = os.path.join(ROOT, "data", "raw", "networks", "intact", "intact.zip")
MAPPING_PATH = os.path.join(ROOT, "data", "raw", "reference", "hgnc_complete_set.txt")
OUT_DIR = os.path.join(ROOT, "data", "processed", "networks")
MI_THRESHOLD = 0.45
PRIMARY_NET = "intact"
SENSITIVITY_NET = "intact_miscore045"


def _strip_prefix(s: str) -> str:
    """IntAct IDs are 'uniprotkb:P37840' or 'uniprotkb:P37840-2'; strip the prefix.
    Isoform suffixes are stripped automatically by HGNCMapper."""
    if s.startswith("uniprotkb:"):
        return s[len("uniprotkb:"):]
    return s


def _read_intact():
    """Stream the large intact.txt from inside intact.zip; filter taxid; extract pairs."""
    rows = []
    with zipfile.ZipFile(RAW_ZIP, "r") as z:
        with z.open("intact.txt") as f:
            header = f.readline().decode("utf-8").rstrip("\n").split("\t")
            # Column indices (PSI-MITAB 2.7, 42 columns):
            # 0  #ID(s) interactor A
            # 1  ID(s) interactor B
            # 9  Taxid interactor A
            # 10 Taxid interactor B
            # 14 Confidence value(s)
            for raw in f:
                line = raw.decode("utf-8", errors="replace").rstrip("\n")
                if not line:
                    continue
                cells = line.split("\t")
                if len(cells) < 15:
                    continue
                a_id = _strip_prefix(cells[0])
                b_id = _strip_prefix(cells[1])
                tax_a = cells[9]
                tax_b = cells[10]
                # Only human-human
                if "9606" not in tax_a or "9606" not in tax_b:
                    continue
                # Confidence: parse "intact-miscore:0.67"
                mi = 0.0
                for tok in cells[14].split("|"):
                    if tok.startswith("intact-miscore:"):
                        try:
                            mi = float(tok.split(":", 1)[1])
                        except ValueError:
                            mi = 0.0
                        break
                rows.append((a_id, b_id, mi))
    print(f"  read {len(rows):,} human-human pairs from intact.txt")
    return pd.DataFrame(rows, columns=["gene_a", "gene_b", "mi"])


def _map_and_finalize(df, net_id, mi_threshold, mapper):
    """Map UniProt -> HGNC symbol via B's HGNCMapper, then finalize_edges."""
    if mi_threshold is not None:
        before = len(df)
        df = df[df["mi"] >= mi_threshold]
        print(f"  after MI >= {mi_threshold}: {len(df):,}/{before:,}")

    # Map both columns. The mapper already strips isoform suffixes.
    df = df.assign(
        gene_a=mapper.map(df["gene_a"], "uniprot"),
        gene_b=mapper.map(df["gene_b"], "uniprot"),
    )
    before = len(df)
    df = df.dropna(subset=["gene_a", "gene_b"])
    print(f"  after UniProt -> HGNC: {len(df):,}/{before:,} ({100*len(df)/before:.1f}%)")

    # finalize_edges handles self-loops, dedup, alphabetise, LCC, etc.
    out_df, meta = finalize_edges(
        df[["gene_a", "gene_b"]], net_id=net_id,
        meta={"source_file": "intact.txt (in intact.zip)", "threshold": mi_threshold,
              "taxid_filter": "9606 (human-human only)",
              "score_column": "intact-miscore"},
        out_dir=OUT_DIR,
    )
    print(f"  after finalize_edges: {len(out_df):,} edges, "
          f"{meta['n_nodes_lcc']:,} nodes in LCC, "
          f"removed {meta['n_self_loops_removed']} self-loops, "
          f"{meta['n_edges_dedup'] - len(out_df)} duplicates merged")

    check_network_file(os.path.join(OUT_DIR, f"{net_id}.tsv"))
    print(f"  OK: {net_id}.tsv passes contract C2")


def main():
    print("IntAct networks (A3)")
    print("Step 1: building HGNC mapper (UniProt -> HGNC symbol)")
    mapper = HGNCMapper(MAPPING_PATH)
    print(f"  mapper covers {len(mapper.maps['uniprot']):,} UniProt -> HGNC entries")

    print("\nStep 2: reading IntAct (streaming, filter human-human)")
    raw = _read_intact()

    print("\nStep 3: primary arm (no MI threshold)")
    _map_and_finalize(raw.copy(), PRIMARY_NET, mi_threshold=None, mapper=mapper)

    print("\nStep 4: sensitivity arm (MI >= 0.45)")
    _map_and_finalize(raw.copy(), SENSITIVITY_NET, mi_threshold=MI_THRESHOLD, mapper=mapper)


if __name__ == "__main__":
    main()