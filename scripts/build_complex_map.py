"""Build the gene -> protein-complex map for the complex-aware CV (Phase D, task C1).

Source and why
--------------
The canonical ComplexPortal dump is NOT reachable from this machine
(no external FTP or website access; ftp.ebi.ac.uk returns 404 for the
ComplexPortal paths and www.ebi.ac.uk/complexportal is a JavaScript
single-page app). UniProt's REST API IS reachable and carries the
ComplexPortal cross-references of every reviewed protein, so the
complex identifiers here are ComplexPortal identifiers (CPX-...), taken
from UniProt's curation of them.

Two properties make this sufficient for the leakage audit and the
complex-aware splitter:

  * two genes are complex-linked iff their CPX-ID lists intersect, so
    the member lists are not needed -- only gene -> CPX-IDs;
  * gene -> UniProt accession is an EXACT mapping from the HGNC
    complete set already in the repo (data/raw/reference/
    hgnc_complete_set.txt), so no symbol has to be guessed.

Provenance (release, date, query, SHA-256 of the raw responses) is
written to data/processed/complexes/provenance.md.

Usage:
    E:\\anaconda\\python.exe scripts/build_complex_map.py
"""
from __future__ import annotations

import hashlib
import subprocess
import sys
import time
import urllib.parse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "complexes"
OUT = ROOT / "data" / "processed" / "complexes"
NET_DIR = ROOT / "data" / "processed" / "networks"
HGNC = ROOT / "data" / "raw" / "reference" / "hgnc_complete_set.txt"
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
BATCH = 60
SEARCH = "https://rest.uniprot.org/uniprotkb/search"


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


CANONICAL = __import__("re").compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")


def curl(url: str, tries: int = 2):
    """(body, headers) with retries; None when the query is rejected.

    A single invalid term (an obsolete or isoform accession) makes
    UniProt reject the WHOLE query with HTTP 400, so the caller splits
    the batch and isolates the offending term.
    """
    last = None
    for i in range(tries):
        r = subprocess.run(
            ["curl.exe", "-s", "-D", "-", "--max-time", "120", "-L", url],
            capture_output=True, text=True)
        txt = r.stdout
        # split headers from body
        parts = txt.split("\r\n\r\n", 1)
        if len(parts) == 2:
            head, body = parts[0], parts[1]
        else:
            head, body = "", txt
        if "Entry\t" in body or body.strip().startswith("Entry"):
            return body, head
        last = body[:200].replace("\n", " ")
        time.sleep(1 + 2 * i)
    log(f"    query rejected ({last[:90]}); the batch will be split")
    return None


def panel_genes() -> list[str]:
    genes = set()
    for net in NETWORKS:
        df = pd_read_tsv(NET_DIR / f"{net}.tsv")
        genes |= set(df["gene_a"].astype(str)) | set(df["gene_b"].astype(str))
    return sorted(genes)


def pd_read_tsv(path: Path):
    import pandas as pd
    return pd.read_csv(path, sep="\t", dtype=str)


def symbol_to_uniprot(symbols: list[str]) -> dict[str, str]:
    """symbol -> UniProt accession (exact mapping from the HGNC set)."""
    import pandas as pd
    hgnc = pd.read_csv(HGNC, sep="\t", dtype=str)
    hgnc = hgnc[hgnc["symbol"].isin(set(symbols))]
    out = {}
    for _, r in hgnc.iterrows():
        acc = r.get("uniprot_ids")
        if not isinstance(acc, str) or not acc.strip():
            continue
        accs = [a.strip() for a in acc.split(";") if a.strip()]
        if not accs:
            continue
        # prefer the reviewed (Swiss-Prot) accession
        chosen = None
        for a in accs:
            if "-" not in a and len(a) in (6, 10):
                chosen = a
                break
        if chosen is None:
            chosen = accs[0]
        out[r["symbol"]] = chosen
    return out


def fetch_chunk(chunk: list[str], acc_to_syms: dict[str, list[str]],
               headers_seen: list[str]) -> list[tuple]:
    """Query UniProt for one chunk of accessions; split the query when
    UniProt rejects it (it refuses queries with too many OR terms)."""
    q = " OR ".join(f"accession:{a}" for a in chunk)
    params = urllib.parse.urlencode({
        "query": f"({q}) AND organism_id:9606",
        "format": "tsv",
        "fields": "accession,gene_primary,xref_complexportal",
        "size": "500"})
    res = curl(f"{SEARCH}?{params}")
    if res is None:
        body = head = None
    else:
        body, head = res
    if body is None:
        half = len(chunk) // 2
        if half < 5:
            log(f"    giving up on {len(chunk)} accessions")
            return []
        out = fetch_chunk(chunk[:half], acc_to_syms, headers_seen)
        out += fetch_chunk(chunk[half:], acc_to_syms, headers_seen)
        return out
    for h in head.splitlines():
        if h.lower().startswith("x-uniprot-release") and h not in headers_seen:
            headers_seen.append(h)
    out = []
    for line in body.splitlines():
        if line.startswith("Entry\t") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        acc = parts[0].strip()
        cpx = [c.strip() for c in parts[2].split(";") if c.strip()]
        if not cpx:
            continue
        for c in cpx:
            for sym in acc_to_syms.get(acc, []):
                out.append((c, sym, acc))
    return out


def main() -> int:
    genes = panel_genes()
    log(f"panel genes (union of the six LCCs): {len(genes):,}")
    mapping = symbol_to_uniprot(genes)
    log(f"mapped to a UniProt accession: {len(mapping):,} "
        f"({len(genes)-len(mapping):,} unmapped)")

    acc_to_syms: dict[str, list[str]] = {}
    for sym, acc in mapping.items():
        acc_to_syms.setdefault(acc, []).append(sym)
    # keep only canonical UniProtKB accessions: an isoform suffix
    # ("P04637-2") or an obsolete id makes UniProt reject the whole
    # query with HTTP 400
    dropped = {a: syms for a, syms in acc_to_syms.items()
               if not CANONICAL.match(a)}
    for a in dropped:
        del acc_to_syms[a]
    accs = sorted(acc_to_syms)
    log(f"distinct accessions to query: {len(accs):,} "
        f"({len(dropped):,} non-canonical accessions dropped)")

    raw_path = RAW / "uniprot_complexportal.tsv"
    rows = []
    headers_seen = []
    n_batches = (len(accs) + BATCH - 1) // BATCH
    for bi in range(n_batches):
        chunk = accs[bi * BATCH:(bi + 1) * BATCH]
        rows.extend(fetch_chunk(chunk, acc_to_syms, headers_seen))
        if bi % 20 == 0 or bi == n_batches - 1:
            log(f"  batch {bi+1}/{n_batches}: {len(rows):,} (gene, complex) rows")
    # de-duplicate
    rows = sorted(set(rows))
    log(f"raw (gene, complex) pairs: {len(rows):,}")
    raw_path.write_text("\n".join("\t".join(r) for r in rows) + "\n",
                        encoding="utf-8")
    sha = hashlib.sha256(raw_path.read_bytes()).hexdigest()

    import pandas as pd
    df = pd.DataFrame(rows, columns=["complex_id", "gene", "uniprot"])
    df["source"] = "ComplexPortal (via UniProtKB cross-reference)"
    df["n_members"] = df.groupby("complex_id")["gene"].transform("size")
    out_path = OUT / "complex_map.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    n_complex = df["complex_id"].nunique()
    n_genes = df["gene"].nunique()
    log(f"wrote {out_path.name}: {n_complex:,} complexes, {n_genes:,} genes")

    prov = OUT / "provenance.md"
    rel = next((h for h in headers_seen if "release" in h.lower()), "unknown")
    prov.write_text(f"""# Provenance of the complex map

**Built:** {time.strftime('%Y-%m-%d %H:%M')} by `scripts/build_complex_map.py`

## Source

ComplexPortal complex identifiers (CPX-...) as curated in **UniProtKB**,
retrieved through the UniProt REST search API
(`https://rest.uniprot.org/uniprotkb/search`). UniProt release header:
`{rel}`.

**Why not the canonical ComplexPortal dump:** this machine has no
external FTP or website access. `ftp.ebi.ac.uk/pub/databases/complexportal/`
returns 404 and `www.ebi.ac.uk/complexportal` is a JavaScript
single-page application with no static data files; CORUM
(`mips.helmholtz-muenchen.de`) is unreachable. UniProt's REST API is
reachable and carries the ComplexPortal cross-references, so the
identifiers are ComplexPortal identifiers and the underlying curation is
ComplexPortal's. The difference from the canonical dump is that only the
complexes annotated to a UniProt entry are seen, and the complex
*member list* is not used at all (see below).

## Why the member list is not needed

Two genes are complex-linked iff their CPX-ID lists intersect. The
leakage audit (D1) and the complex-aware splitter (D2) both work on
gene-level complex links, so only `gene -> CPX-IDs` is required.

## Gene identifier mapping

`data/raw/reference/hgnc_complete_set.txt` (the HGNC complete set,
already in the repo) supplies the exact `symbol -> UniProt accession`
mapping from its `uniprot_ids` column, preferring the reviewed
(Swiss-Prot) accession when several are listed. No symbol was guessed.

## Counts

| quantity | value |
|---|---|
| panel genes (union of the six networks' LCCs) | {len(genes):,} |
| mapped to a UniProt accession | {len(mapping):,} |
| unmapped | {len(genes)-len(mapping):,} |
| accessions queried | {len(accs):,} |
| complexes found | {n_complex:,} |
| distinct genes in at least one complex | {n_genes:,} |
| (gene, complex) pairs | {len(rows):,} |

## Raw artefact

`data/raw/complexes/uniprot_complexportal.tsv` -- SHA-256
`{sha}`

## Caveat to state in the paper

The complex set is UniProt's view of ComplexPortal, not ComplexPortal's
complete dump, and the complexes are the *curated* ones. It is used
only to decide which genes may share a fold; a gene in no complex is
treated as unconstrained, exactly as the protocol specifies.
""", encoding="utf-8")
    log(f"wrote {prov.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
