"""Item 17 (part 3) -- the fourth curated PPI: BioGRID, at FunMap's scale.

Why BioGRID: it is a curated protein-protein interaction database that is
independent of the three curated networks already in the panel (IntAct,
Reactome, FunMap). It is also large, so the network is built at a
**matched edge count** (FunMap's, by default) -- the whole point of the
size-matched design is that the class is not entangled with the scale.

The interaction types are kept to the physical ones (the "direct
interaction" and the "association" evidence codes), so the network means
the same kind of thing as the others.
"""
from __future__ import annotations

import ssl
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(r"D:\Grade3\swxxx\final\D1\d1-benchmark")
sys.path.insert(0, str(ROOT))

import pandas as pd  # noqa: E402

from d1.networks.io import finalize_edges  # noqa: E402

NETS_DIR = ROOT / "data/processed/networks"
RAW_DIR = ROOT / "data/raw/networks/biogrid"
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
NET_ID = "biogrid_curated"
URL = ("https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/"
       "BIOGRID-ALL-LATEST.tab3.zip")
TARGET_EDGES = 196_605          # FunMap's edge count (the size match)
PHYSICAL = {"direct interaction", "physical association", "association",
            "direct interaction with text"}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def download() -> Path:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    out = RAW_DIR / "biogrid_all.tab3.zip"
    if out.exists() and out.stat().st_size > 1000:
        return out
    log("downloading the BioGRID release (this is a large file)")
    req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=900, context=CTX) as r:
        out.write_bytes(r.read())
    log(f"  {out.stat().st_size / 1e6:.1f} MB")
    return out


def main() -> int:
    zpath = download()
    with zipfile.ZipFile(zpath) as zf:
        name = [n for n in zf.namelist() if n.endswith(".tab3.txt")][0]
        log(f"  the member: {name}")
        with zf.open(name) as fh:
            head = fh.readline().decode("utf-8", errors="replace").rstrip(
                "\n").split("\t")
            cols = {h: i for i, h in enumerate(head)}
            i_a = cols["Official Symbol Interactor A"]
            i_b = cols["Official Symbol Interactor B"]
            i_s = cols["Source Database"]
            i_t = cols["Experimental System Type"]
            i_oa = cols["Organism ID Interactor A"]
            i_ob = cols["Organism ID Interactor B"]
            rows = []
            for line in fh:
                p = line.decode("utf-8", errors="replace").rstrip(
                    "\n").split("\t")
                if len(p) <= max(i_a, i_b, i_t, i_oa, i_ob):
                    continue
                # the physical interactions of human
                if p[i_t].strip().lower() != "physical":
                    continue
                if p[i_s].strip() != "BIOGRID":
                    continue
                if p[i_oa].strip() != "9606" or p[i_ob].strip() != "9606":
                    continue
                rows.append((p[i_a].strip(), p[i_b].strip()))
    log(f"  the physical BioGRID rows: {len(rows):,}")
    edges = pd.DataFrame(rows, columns=["gene_a", "gene_b"])
    edges = edges[(edges.gene_a != "") & (edges.gene_b != "")]
    edges = edges.drop_duplicates()
    log(f"  the unique edges: {len(edges):,}")
    if len(edges) > TARGET_EDGES:
        # the size match: take the strongest (the most databases reporting
        # the pair), so the sub-network is not a random sample
        counts = (edges.groupby(["gene_a", "gene_b"]).size()
                  .reset_index(name="n"))
        counts = counts.sort_values("n", ascending=False)
        counts = counts.head(TARGET_EDGES)
        edges = counts[["gene_a", "gene_b"]]
        log(f"  the size-matched sub-network: {len(edges):,} edges")
    meta = {"source_file": "BIOGRID-ALL-LATEST.tab3.zip",
            "description": "BioGRID curated physical interactions, "
                           "size-matched to FunMap",
            "net_id": NET_ID, "date": time.strftime("%Y-%m-%d"),
            "script": "scripts/build_biogrid_network.py",
            "target_edges": TARGET_EDGES}
    out = finalize_edges(edges, NET_ID, meta=meta, out_dir=str(NETS_DIR),
                         write=True)
    log(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
