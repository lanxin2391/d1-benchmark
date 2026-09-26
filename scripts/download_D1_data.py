#!/usr/bin/env python3
"""
D1 data downloader - fetches every public dataset the D1 protocol needs.

Usage (PowerShell, from the folder where you want the data):
    python download_D1_data.py              # download everything
    python download_D1_data.py --only string funmap   # just some groups
    python download_D1_data.py --list       # show what would be downloaded

- Standard library only (no pip install needed). Python 3.8+.
- Resumable: re-run the script after a failure and it continues where it stopped.
- Skips files that already exist with the right size / MD5.
- Detects fake downloads (HTML error pages saved as .gz/.tgz, like the 1 KB
  funmap file you got from Zenodo's 403 page) and deletes them.
"""
import argparse, hashlib, os, re, sys, time, urllib.request, urllib.error
from html.parser import HTMLParser

ROOT = os.path.join(os.getcwd(), "D1_data")
UA = "Mozilla/5.0 (D1-benchmark-downloader; academic use)"
OT = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.06/output"

# (group, relative path, url, md5 or None, note)
FILES = [
    # ---------- Networks ----------
    ("funmap", "networks/funmap/funmap.tsv",
     "https://funmap.linkedomics.org/data/download/funmap.tsv", None,
     "FunMap edge list (protein-derived network, ~196,800 edges)"),
    ("funmap", "networks/funmap/funmap_input_expression_data.tgz",
     "https://zenodo.org/records/7948944/files/funmap_input_expression_data.tgz?download=1",
     "f490643dc8bec13dac72ed818d8b707d",
     "FunMap input mRNA + protein matrices, 228.6 MB (for the RNA network)"),
    ("string", "networks/string/9606.protein.links.detailed.v12.0.txt.gz",
     "https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz",
     None, "STRING full, per-channel scores (you already have this, 133 MB)"),
    ("string", "networks/string/9606.protein.physical.links.v12.0.txt.gz",
     "https://stringdb-downloads.org/download/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz",
     None, "STRING physical subnetwork"),
    ("string", "networks/string/9606.protein.info.v12.0.txt.gz",
     "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz",
     None, "ENSP -> gene symbol mapping"),
    ("string", "networks/string/9606.protein.aliases.v12.0.txt.gz",
     "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz",
     None, "ENSP -> HGNC/Entrez aliases"),
    ("reactome", "networks/reactome/reactome.homo_sapiens.interactions.tab-delimited.txt",
     "https://reactome.org/download/current/interactors/reactome.homo_sapiens.interactions.tab-delimited.txt",
     None, "Reactome human interactions"),
    ("intact", "networks/intact/intact.zip",
     "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip",
     None, "IntAct PSI-MITAB, all species (large, ~1+ GB); filter taxid:9606 later"),
    # ---------- Labels ----------
    ("intogen", "labels/intogen/IntOGen-Drivers-20240920.zip",
     "https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip", None,
     "IntOGen 2024.09.20 driver catalogue (label set 1)"),
    ("intogen", "labels/intogen/IntOGen-Drivers-20230531.zip",
     "https://www.intogen.org/download?file=IntOGen-Drivers-20230531.zip", None,
     "IntOGen 2023.05.31 (temporal split, control 6)"),
    ("intogen", "labels/intogen/IntOGen-Drivers-20200201.zip",
     "https://www.intogen.org/download?file=IntOGen-Drivers-20200201.zip", None,
     "IntOGen 2020.02.01 (temporal split, control 6)"),
    ("clingen", "labels/clingen/clingen_gene_validity.csv",
     "https://search.clinicalgenome.org/kb/gene-validity/download", None,
     "ClinGen gene-disease validity"),
    # ---------- Cohorts ----------
    ("tcga", "cohorts/mc3.v0.2.8.PUBLIC.maf.gz",
     "https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc", None,
     "TCGA PanCanAtlas MC3 somatic mutations, ALL 33 cancer types (~750 MB)"),
    # ---------- Reference ----------
    ("reference", "reference/hgnc_complete_set.txt",
     "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
     None, "HGNC symbols (harmonise all gene IDs)"),
    ("reference", "reference/gene2pubmed.gz",
     "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz", None,
     "Per-gene publication counts (control 7 / H3)"),
]

# Open Targets 26.06 datasets are folders of parquet files -> crawled
OT_DATASETS = [
    ("association_by_datasource_direct", "labels 2+3: lets you drop the literature (europepmc) channel"),
    ("association_by_datatype_direct", "per-datatype scores"),
    ("association_overall_direct", "all-evidence score (label set 2)"),
    ("target", "target annotations (Ensembl ID -> symbol)"),
    ("disease", "disease ontology (find oncology / neoplasm EFO_0000616 descendants)"),
]


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def looks_like_html(path):
    try:
        with open(path, "rb") as f:
            head = f.read(512).lstrip().lower()
        return head.startswith(b"<!doctype html") or head.startswith(b"<html")
    except OSError:
        return False


def fetch(url, dest, md5=None, retries=6, allow_html=False):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest):
        if not allow_html and looks_like_html(dest):
            print(f"  ! {os.path.basename(dest)} is an HTML error page, deleting")
            os.remove(dest)
        elif md5 and md5sum(dest) != md5:
            print(f"  ! MD5 mismatch on existing file, re-downloading")
            os.remove(dest)
        else:
            print(f"  = already have {os.path.basename(dest)} "
                  f"({os.path.getsize(dest)/1e6:.1f} MB)")
            return True
    part = dest + ".part"
    for attempt in range(1, retries + 1):
        try:
            have = os.path.getsize(part) if os.path.exists(part) else 0
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            if have:
                req.add_header("Range", f"bytes={have}-")
            with urllib.request.urlopen(req, timeout=60) as r:
                if have and r.status != 206:      # server ignored Range
                    have = 0
                total = r.headers.get("Content-Length")
                total = int(total) + have if total else None
                mode = "ab" if have else "wb"
                done, t0 = have, time.time()
                with open(part, mode) as f:
                    while True:
                        buf = r.read(1 << 20)
                        if not buf:
                            break
                        f.write(buf)
                        done += len(buf)
                        if total:
                            pct = 100 * done / total
                            spd = (done - have) / max(time.time() - t0, 1e-6) / 1e6
                            print(f"\r    {done/1e6:8.1f}/{total/1e6:.1f} MB "
                                  f"({pct:5.1f}%)  {spd:5.2f} MB/s", end="")
                print()
            if not allow_html and looks_like_html(part):
                os.remove(part)
                raise RuntimeError("server returned an HTML page instead of data "
                                   "(blocked, login form, or wrong link)")
            if md5 and md5sum(part) != md5:
                os.remove(part)
                raise RuntimeError("MD5 mismatch")
            os.replace(part, dest)
            print(f"  + saved {dest}")
            return True
        except (urllib.error.URLError, OSError, RuntimeError) as e:
            wait = min(2 ** attempt, 60)
            print(f"\n  x attempt {attempt}/{retries} failed: {e}")
            if attempt < retries:
                print(f"    retrying in {wait}s ...")
                time.sleep(wait)
    return False


class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for k, v in attrs:
                if k == "href" and v:
                    self.links.append(v)


def list_dir(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        p = LinkParser()
        p.feed(r.read().decode("utf-8", "replace"))
    out = []
    for h in p.links:
        if h.startswith(("?", "/", "..", "http")):
            continue
        out.append(h)
    return out


def crawl(url, dest_dir, depth=0):
    """Recursively download a directory listing (Open Targets parquet folders)."""
    ok = True
    for name in list_dir(url):
        sub = url.rstrip("/") + "/" + name
        if name.endswith("/"):
            if depth < 3:
                ok &= crawl(sub, os.path.join(dest_dir, name.rstrip("/")), depth + 1)
        else:
            ok &= fetch(sub, os.path.join(dest_dir, name))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", help="groups: funmap string reactome intact "
                    "intogen clingen opentargets tcga reference")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    want = set(a.only) if a.only else None

    if a.list:
        for g, p, u, _, n in FILES:
            print(f"[{g:10}] {p}\n             {n}")
        for d, n in OT_DATASETS:
            print(f"[opentargets] labels/opentargets_26.06/{d}/  ({n})")
        return

    print(f"Downloading into {ROOT}\n")
    failed = []
    for g, p, u, m, n in FILES:
        if want and g not in want:
            continue
        print(f"[{g}] {n}")
        if not fetch(u, os.path.join(ROOT, p), m):
            failed.append((p, u))
    if not want or "opentargets" in want:
        for d, n in OT_DATASETS:
            print(f"[opentargets] {d}: {n}")
            try:
                if not crawl(f"{OT}/{d}/", os.path.join(ROOT, "labels", "opentargets_26.06", d)):
                    failed.append((d, f"{OT}/{d}/"))
            except Exception as e:
                print(f"  x could not list {OT}/{d}/: {e}")
                failed.append((d, f"{OT}/{d}/"))

    print("\n" + "=" * 70)
    if failed:
        print("FAILED - rerun the script, or download these by hand in a browser:")
        for p, u in failed:
            print(f"  {p}\n     {u}")
        print("\nIntOGen tip: if the .zip links fail, open https://www.intogen.org/download "
              "in a browser, download the three releases, and put the zips in "
              "D1_data/labels/intogen/")
        sys.exit(1)
    print("All D1 data downloaded.")


if __name__ == "__main__":
    main()
