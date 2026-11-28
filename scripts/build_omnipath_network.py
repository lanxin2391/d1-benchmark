"""Item 17 (part 3) -- the fourth curated PPI, at FunMap's scale.

The class effect cannot be separated from the network's size, because every
provenance class happens to sit at its own scale. The fix is a **fourth
curated PPI at FunMap's scale** plus a third co-expression variant, so that
every class has members inside a comparable size band.

The source is OmniPath's curated interaction aggregate (reachable from this
machine, no credential). Its default download gives UniProt accessions, so
the symbols come from the STRING alias table we already have.
"""
from __future__ import annotations

import gzip
import ssl
import sys
import time
import urllib.request
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(r"D:\Grade3\swxxx\final\D1\d1-benchmark")
sys.path.insert(0, str(ROOT))

import pandas as pd  # noqa: E402

from d1.networks.io import finalize_edges  # noqa: E402

NETS_DIR = ROOT / "data/processed/networks"
RAW_DIR = ROOT / "data/raw/networks/omnipath"
ALIASES = ROOT / "data/raw/networks/string/9606.protein.aliases.v12.0.txt.gz"
NET_ID = "omnipath_curated"
# the curated interaction sources (protein-protein, not the pathways)
CURATED = {"reactome", "intact", "signor", "corum", "matrixdb", "proteinfull",
            "hpmr", "mentions", "phosphosite_plus", "biogrid", "intact_complex",
            "ncbi_complex", "complexportal", "cello", "hpa", "mips", "psicquic",
            "reactome_physical", "corum_3d", "hpa_secretome"}
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


DATASETS = ["reactome", "intact", "signor", "corum", "matrixdb",
             "proteinfull", "hpmr", "mentions", "biogrid", "complexportal",
             "ncbi_complex", "intact_complex", "phosphosite_plus", "hpa",
             "irefindex", "pdb"]


def download() -> Path:
    """One request per curated database (the default download does not
    carry the database names)."""
    import urllib.parse
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    out = RAW_DIR / "omnipath_interactions.tsv"
    if out.exists() and out.stat().st_size > 1000:
        return out
    import urllib.request
    frames = []
    for ds in DATASETS:
        url = ("https://omnipathdb.org/interactions?"
               + urllib.parse.urlencode({"datasets": ds}))
        try:
            req = urllib.request.Request(url,
                                         headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=180, context=CTX) as r:
                body = r.read()
        except Exception as exc:
            log(f"  {ds}: {type(exc).__name__}")
            continue
        if body.startswith(b"Something is not entirely good"):
            log(f"  {ds}: the API refused the request")
            continue
        frames.append(body.decode("utf-8", errors="replace"))
        log(f"  {ds}: {len(body):,} bytes")
    if not frames:
        raise SystemExit("no OmniPath database could be downloaded")
    header = frames[0].splitlines()[0]
    rows = []
    for f in frames[1:]:
        lines = f.splitlines()
        if lines and lines[0] == header:
            lines = lines[1:]
        rows.extend(lines)
    out.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
    log(f"  {out.stat().st_size:,} bytes in total")
    return out


def uniprot_to_symbol() -> dict[str, str]:
    """The UniProt accession -> the gene symbol (the STRING alias table)."""
    out: dict[str, str] = {}
    with gzip.open(ALIASES, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_src = header.index("source")
        i_alias = header.index("alias")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(i_src, i_alias):
                continue
            if p[i_src].startswith("uniprot"):
                out.setdefault(p[i_alias], p[i_alias])
    # the alias table maps the accession to the alias string; the canonical
    # symbol is the one in the "protein.info" table, so the accession's own
    # symbol is taken from the info table's preferred name
    return out


def main() -> int:
    raw_path = download()
    df = pd.read_csv(raw_path, sep="\t", dtype=str)
    log(f"the raw table: {df.shape}")
    # the UniProt accession -> the gene symbol, from the STRING info table
    import gzip
    info = ROOT / ("data/raw/networks/string/9606.protein.info.v12.0.txt.gz")
    # the STRING protein id -> the preferred name (the symbol)
    id2sym: dict[str, str] = {}
    with gzip.open(info, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_id = header.index("#string_protein_id")
        i_pref = header.index("preferred_name")
        for line in fh:
            q = line.rstrip("\n").split("\t")
            if len(q) > max(i_id, i_pref):
                id2sym.setdefault(q[i_id], q[i_pref])
    # the alias table: the UniProt accession -> the STRING protein id
    aliases = ROOT / ("data/raw/networks/string/9606.protein.aliases.v12.0.txt.gz")
    acc2sym: dict[str, str] = {}
    with gzip.open(aliases, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_src = header.index("source")
        i_alias = header.index("alias")
        for line in fh:
            q = line.rstrip("\n").split("\t")
            if len(q) <= max(i_src, i_alias):
                continue
            if ("uniprot" in q[i_src].lower()
                    and q[i_alias] in id2sym):
                acc2sym.setdefault(q[i_src], id2sym[q[i_alias]])
    log(f"  the UniProt mappings: {len(acc2sym):,}")
    log(f"  the curated rows: {len(df):,}")
    edges = pd.DataFrame({
        "gene_a": df.source.map(acc2sym),
        "gene_b": df.target.map(acc2sym)})
    edges = edges.dropna()
    edges = edges[edges.gene_a != edges.gene_b]
    edges = edges.drop_duplicates()
    log(f"  the harmonised edges: {len(edges):,}")
    meta = {"source_file": raw_path.name,
            "description": "OmniPath curated protein-protein interactions",
            "net_id": NET_ID, "date": time.strftime("%Y-%m-%d"),
            "script": "scripts/build_omnipath_network.py"}
    out = finalize_edges(edges, NET_ID, meta=meta, out_dir=str(NETS_DIR),
                         write=True)
    log(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
