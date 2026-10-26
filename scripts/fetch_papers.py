#!/usr/bin/env python3
"""Fetch the PDFs of every paper the D1 benchmark cites.

Sources, in priority order per reference:
  1. Europe PMC REST API (https://www.ebi.ac.uk/europepmc/webservices/rest)
     - OA articles: PDF pulled from the Europe PMC render endpoint.
  2. arXiv (for conference papers behind a paywall).
  3. A hard-coded direct URL (author-hosted copy of a classic paper).
  4. Otherwise: reported as MANUAL so the author can use institutional access.

All DOIs in REFERENCES were verified against Crossref on 2026-10-26.
The script also looks for NEWER database-issue updates of the resources we
cite (the "Q" series), because the foundational paper of a database is often
several years older than its most recent NAR update.

Usage (PowerShell, from the repo root):
    python scripts/fetch_papers.py --list          # show the reference list
    python scripts/fetch_papers.py                 # download everything
    python scripts/fetch_papers.py --only P1 P2    # just these ids
    python scripts/fetch_papers.py --only-manual   # list what still needs a human

Output:
    data/papers/<id>_<slug>.pdf
    data/papers/manifest.tsv   (status per reference: OK / MANUAL / ERROR)
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

# Windows consoles default to GBK here; author names contain non-ASCII (e.g.
# "del Toro" is fine but some surnames carry accents). Force UTF-8 output.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "data" / "papers"

EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EPMC_PDF = "https://europepmc.org/articles/{pmcid}?pdf=render"
ARXIV_PDF = "https://arxiv.org/pdf/{arxiv_id}"

UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36")
CURL_HEADERS = [
    "-A", UA,
    "-H", "Accept: text/html,application/xhtml+xml,application/pdf,*/*;q=0.8",
    "-H", "Accept-Language: en-US,en;q=0.9",
    "-H", "Sec-Fetch-Dest: document",
    "-H", "Sec-Fetch-Mode: navigate",
    "-H", "Sec-Fetch-Site: none",
]

# id, slug, citation, lookup spec
#   spec keys: doi | title | arxiv | direct (list of fallback URLs)
REFERENCES = [
    # ---------- the 16 papers the manuscript cites ----------
    ("P1", "funmap_natcancer_2025",
     "Shi Z, Lei JT, Elizarraras JM, Zhang B. 2025. Mapping the functional network of human cancer through machine learning and pan-cancer proteogenomics. Nat Cancer 6:205-222. doi:10.1038/s43018-024-00869-z",
     {"doi": "10.1038/s43018-024-00869-z"}),
    ("P2", "string_v12_nar_2023",
     "Szklarczyk D, Kirsch R, Koutrouli M, et al. 2023. The STRING database in 2023. Nucleic Acids Res 51:D638-D646. doi:10.1093/nar/gkac1000",
     {"doi": "10.1093/nar/gkac1000"}),
    ("P3", "intact_nar_2022",
     "del Toro N, Shrivastava A, Ragueneau E, et al. 2022. The IntAct database: efficient access to fine-grained molecular interaction data. Nucleic Acids Res 50:D648-D653. doi:10.1093/nar/gkab1006",
     {"doi": "10.1093/nar/gkab1006"}),
    ("P4", "reactome_nar_2024",
     "Milacic M, Beavers D, Conley P, et al. 2024. The Reactome Pathway Knowledgebase 2024. Nucleic Acids Res 52:D672-D678. doi:10.1093/nar/gkad1025",
     {"doi": "10.1093/nar/gkad1025"}),
    ("P5", "intogen_natrevcancer_2020",
     "Martinez-Jimenez F, Muiños F, Sentís I, et al. 2020. A compendium of mutational cancer driver genes. Nat Rev Cancer 20:555-572. doi:10.1038/s41568-020-0290-x",
     {"doi": "10.1038/s41568-020-0290-x"}),
    ("P6", "opentargets_nar_2023",
     "Ochoa D, Hercules A, Carmona M, et al. 2023. The next-generation Open Targets Platform: reimagined, redesigned, rebuilt. Nucleic Acids Res 51:D1353-D1359. doi:10.1093/nar/gkac1046",
     {"doi": "10.1093/nar/gkac1046"}),
    ("P7", "clingen_nejm_2015",
     "Rehm HL, Berg JS, Brooks LD, et al. 2015. ClinGen - the Clinical Genome Resource. N Engl J Med 372:2235-2242. doi:10.1056/NEJMsr1406261",
     {"doi": "10.1056/NEJMsr1406261"}),
    ("P8", "hgnc_nar_2023",
     "Seal RL, Braschi B, Gray K, et al. 2023. Genenames.org: the HGNC resources in 2023. Nucleic Acids Res 51:D1003-D1009. doi:10.1093/nar/gkac888",
     {"doi": "10.1093/nar/gkac888"}),
    ("P9", "refseq_nar_2016",
     "Maglott DR, Kim J, Tipton ML, et al. 2016. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res 44:D733-D745. doi:10.1093/nar/gkv1189",
     {"doi": "10.1093/nar/gkv1189"}),
    ("P10", "node2vec_kdd_2016",
     "Grover A, Leskovec J. 2016. node2vec: Scalable Feature Learning for Networks. ACM KDD. arXiv:1403.6652",
     {"arxiv": "1403.6652"}),
    ("P11", "bh_fdr_1995",
     "Benjamini Y, Hochberg Y. 1995. Controlling the False Discovery Rate: a Practical and Powerful Approach to Multiple Testing. J R Stat Soc Series B 57:289-300. doi:10.1111/j.2517-6161.1995.tb02031.x",
     {"doi": "10.1111/j.2517-6161.1995.tb02031.x"}),
    ("P12", "laird_ware_1982",
     "Laird NM, Ware JH. 1982. Random-effects models for longitudinal data. Biometrics 38:963-974. doi:10.2307/2529876",
     {"doi": "10.2307/2529876",
      "direct": ["https://www.jstor.org/stable/pdf/2529876.pdf"]}),
    ("P13", "efron_bootstrap_1979",
     "Efron B. 1979. Bootstrap methods: another look at the jackknife. Ann Statist 7:1-26. doi:10.1214/aos/1176344552",
     {"doi": "10.1214/aos/1176344552",
      "direct": ["https://www-stat.stanford.edu/~johns/papers/annstat_bootstrap.pdf"]}),
    ("P14", "schuirmann_tost_1987",
     "Schuirmann DJ. 1987. A comparison of the Two One-Sided Tests Procedure and the Power Approach for assessing the equivalence of average bioavailability. J Pharmacokinet Biopharm 15:657-680. doi:10.1007/BF01068419",
     {"doi": "10.1007/BF01068419",
      "direct": ["https://link.springer.com/content/pdf/10.1007/BF01068419.pdf"]}),
    ("P15", "pearson_1895",
     "Pearson K. 1895. Notes on Regression and Inheritance. Proc R Soc Lond 58:1-30. doi:10.1098/rspl/1895.0060",
     {"doi": "10.1098/rspl/1895.0060",
      "direct": ["https://royalsocietypublishing.org/doi/pdf/10.1098/rspl/1895.0060"]}),
    ("P16", "fisher_z_1915",
     "Fisher RA. 1915. Frequency distribution of the values of the correlation coefficients in samples from an indefinitely large population. Biometrika 10:507-521. doi:10.1093/biomet/10.4.507",
     {"doi": "10.1093/biomet/10.4.507",
      "direct": ["https://academic.oup.com/biomet/article-pdf/10/4/507/1119416/10.4.507.pdf"]}),
    # ---------- newer database-issue updates (work_plan Phase 1) ----------
    ("Q1", "string_nar_2025",
     "Szklarczyk D, Nastou K, Koutrouli M, et al. 2025. The STRING database in 2025: protein networks with directionality of regulation. Nucleic Acids Res 53:D730-D737. doi:10.1093/nar/gkae1113",
     {"doi": "10.1093/nar/gkae1113"}),
    ("Q2", "hgnc_nar_2026",
     "Seal RL, Braschi B, Gray K, et al. 2026. Genenames.org: the HGNC and PGNC resources in 2026. Nucleic Acids Res 54:D1098-D1107. doi:10.1093/nar/gkaf1229",
     {"doi": "10.1093/nar/gkaf1229"}),
    # ---------- methods for the planned Phase 2 / 6 / 8 work ----------
    ("M1", "gcn_iclr_2017",
     "Kipf TN, Welling M. 2017. Semi-Supervised Classification with Graph Convolutional Networks. ICLR. arXiv:1609.02907",
     {"arxiv": "1609.02907"}),
    ("M2", "graphsage_neurips_2017",
     "Hamilton WL, Ying R, Leskovec J. 2017. Inductive Representation Learning on Large Graphs. NeurIPS. arXiv:1706.02216",
     {"arxiv": "1706.02216"}),
    ("M3", "hotnet2_natgenet_2015",
     "Leiserson MD, Vandin F, Wu HT, et al. 2015. Pan-cancer network analysis identifies combinations of rare somatic mutations across pathways and protein complexes. Nat Genet 47:106-114. doi:10.1038/ng.3168",
     {"doi": "10.1038/ng.3168"}),
]


# ----------------------------------------------------------------- helpers

def curl_download(url: str, dest: Path, timeout: int = 60) -> tuple[bool, str]:
    """Download with curl (survives the TLS/403 quirks that urllib hits)."""
    tmp = dest.with_suffix(".part")
    cmd = ["curl.exe", "-sL", "-m", str(timeout), *CURL_HEADERS,
           "-w", "%{http_code}", "-o", str(tmp), url]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout + 15)
    except subprocess.TimeoutExpired:
        return False, "curl timeout"
    code = (p.stdout or "").strip()[-3:]
    if not tmp.exists():
        return False, f"http={code} no body"
    body = tmp.read_bytes()
    tmp.unlink(missing_ok=True)
    if not body.startswith(b"%PDF"):
        return False, f"http={code} non-PDF ({len(body)} bytes)"
    dest.write_bytes(body)
    return True, f"http={code} {len(body)} bytes"


def epmc_search(spec: dict, tries: int = 5) -> dict | None:
    """Return the best Europe PMC core record for the spec."""
    if "arxiv" in spec:
        return None
    if "doi" in spec:
        query = f'DOI:"{spec["doi"]}"'
    else:
        title = spec.get("title", "")
        query = f'TITLE:"{title}"'
        if spec.get("journal"):
            query += f' AND JOURNAL:"{spec["journal"]}"'
        if spec.get("year_min"):
            query += f" AND (FIRST_PDATE:[{spec['year_min']}-01-01 TO 2026-12-31])"
    url = (f"{EPMC_SEARCH}?query={urllib.parse.quote(query)}"
           f"&resultType=core&format=json&pageSize=5")
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=45) as r:
                data = json.loads(r.read())
            results = data.get("resultList", {}).get("result", [])
            if results:
                return results[0]
            return None
        except Exception as e:
            if attempt == tries - 1:
                print(f"    EPMC lookup failed: {e!r}")
            time.sleep(2.5 * (attempt + 1))
    return None


def pdf_url_from_record(rec: dict) -> str | None:
    for u in rec.get("fullTextUrlList", {}).get("fullTextUrl", []):
        if u.get("availability") == "Open access" and u.get("documentStyle") == "pdf":
            return u.get("url")
    pmcid = rec.get("pmcid")
    if pmcid and rec.get("isOpenAccess") == "Y":
        return EPMC_PDF.format(pmcid=pmcid)
    return None


def fetch_one(ref_id: str, slug: str, citation: str, spec: dict) -> tuple[str, str]:
    dest = OUT_DIR / f"{ref_id}_{slug}.pdf"

    if dest.exists() and dest.stat().st_size > 10_000:
        with dest.open("rb") as fh:
            if fh.read(4) == b"%PDF":
                return "OK", f"already present ({dest.stat().st_size} bytes)"

    notes: list[str] = []

    # --- 1. arXiv ---
    if "arxiv" in spec:
        url = ARXIV_PDF.format(arxiv_id=spec["arxiv"])
        for attempt in range(3):
            ok, note = curl_download(url, dest)
            if ok:
                return "OK", f"arXiv {spec['arxiv']}: {note}"
            notes.append(f"arXiv {note}")
            time.sleep(2 + attempt * 2)
        return "ERROR", "; ".join(notes)

    # --- 2. Europe PMC ---
    rec = epmc_search(spec)
    pmcid = rec.get("pmcid") if rec else None
    year = rec.get("pubYear") if rec else None
    oa = rec.get("isOpenAccess") if rec else None
    if rec:
        pdf_url = pdf_url_from_record(rec)
        if pdf_url:
            for attempt in range(4):
                ok, note = curl_download(pdf_url, dest, timeout=90)
                if ok:
                    return "OK", f"EuropePMC {pmcid} ({year}, OA={oa}): {note}"
                notes.append(f"EuropePMC {note}")
                if "429" in note:
                    time.sleep(15 + attempt * 12)
                else:
                    time.sleep(3 + attempt * 3)
        else:
            notes.append(f"no OA PDF in EPMC (pmcid={pmcid}, year={year}, OA={oa})")
    else:
        notes.append("not found in Europe PMC")

    # --- 3. direct author/publisher URLs ---
    for url in spec.get("direct", []):
        ok, note = curl_download(url, dest, timeout=90)
        if ok:
            return "OK", f"direct: {note} <- {url}"
        notes.append(f"direct {note}")
        time.sleep(2)

    if rec is not None and rec.get("isOpenAccess") != "Y":
        return "MANUAL", (f"paywalled (pmcid={pmcid}, year={year}); "
                          f"needs institutional access - " + " | ".join(notes))
    return "MANUAL", " | ".join(notes)


# ----------------------------------------------------------------- main

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--only-manual", action="store_true",
                    help="list references that still need manual acquisition")
    args = ap.parse_args()

    manifest = OUT_DIR / "manifest.tsv"
    if args.only_manual:
        if not manifest.exists():
            print("no manifest yet - run the script first")
            return 1
        for line in manifest.read_text(encoding="utf-8").splitlines()[1:]:
            f = line.split("\t")
            if len(f) >= 4 and f[3] != "OK":
                print(f"  {f[0]:<4} {f[1]:<30} {f[4] if len(f) > 4 else ''}")
        return 0

    print("Reference list")
    print("=" * 78)
    for ref_id, slug, citation, _ in REFERENCES:
        print(f"  {ref_id:<4} {citation}")
    print()
    if args.list:
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = REFERENCES
    if args.only:
        wanted = {x.upper() for x in args.only}
        targets = [r for r in REFERENCES if r[0] in wanted]

    rows: list[tuple] = []
    for ref_id, slug, citation, spec in targets:
        print(f"[{ref_id}] {citation[:68]}")
        status, note = fetch_one(ref_id, slug, citation, spec)
        print(f"    -> {status}: {note}")
        rows.append((ref_id, slug, citation, status, note))
        time.sleep(2.0)

    lines = ["id\tslug\tcitation\tstatus\tnote"]
    for r in rows:
        lines.append("\t".join(str(x).replace("\t", " ") for x in r))
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")

    n_ok = sum(1 for r in rows if r[3] == "OK")
    n_man = sum(1 for r in rows if r[3] == "MANUAL")
    n_err = sum(1 for r in rows if r[3] == "ERROR")
    print()
    print(f"OK={n_ok}  MANUAL={n_man}  ERROR={n_err}  (total {len(rows)})")
    print(f"PDFs     -> {OUT_DIR}")
    print(f"Manifest -> {manifest}")
    if n_man or n_err:
        print()
        print("Still needed (institutional access / publisher):")
        for r in rows:
            if r[3] != "OK":
                print(f"  {r[0]:<4} {r[1]:<30} {r[4][:110]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
