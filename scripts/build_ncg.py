"""Build the NCG driver-gene label from the Genome Biology 2019
supplementary files (Repana et al., doi:10.1186/s13059-018-1612-0).

Why this label exists: it is an INDEPENDENT arbitration catalogue.
None of the project's four label sets contains it, so it can be used
(a) to measure how much of our "negative" set is actually curated
    driver genes we happen not to label (Bailey et al. 2018; Repana et
    al. 2019) and (b) as an external check on the ranking (Boutros et
    al. 2014). These two uses are ours: the literature review
    (prior_benchmark_literature.md) lists Repana's NCG as the concrete
    case of the "the label sets overlap and leak into each other"
    worry, and the work plan's Phase C item C2 is what first put the
    NCG catalogue on the table.

It also supplies a literature count per gene -- the number of
publications NCG cites for it -- which is the closest thing we have to
Schaefer et al.'s "bait" count, i.e. the study-visibility covariate
that their degree-matched null model is built on.

Outputs:
  data/processed/labels/ncg.tsv
  results/tables/ncg_build_summary.tsv
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

RAW = ROOT / "data" / "raw" / "labels" / "ncg"
LIT_FILE = RAW / "13059_2018_1612_MOESM3_ESM.xlsx"
HGNC = ROOT / "data" / "raw" / "reference" / "hgnc_complete_set.txt"
OUT = ROOT / "data" / "processed" / "labels" / "ncg.tsv"
SUMMARY = ROOT / "results" / "tables" / "ncg_build_summary.tsv"

SHEET_LIT = "2372_cancergenes_literature"
SHEET_SYS = "2372_systemslevelproperties"
SHEET_KNOWN = "711_known_tsgoncogene"
SHEET_CAND = "104_candidatemultiplescreens"


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> int:
    xl = pd.ExcelFile(LIT_FILE)
    lit = xl.parse(SHEET_LIT)
    sysp = xl.parse(SHEET_SYS)
    known = xl.parse(SHEET_KNOWN)
    cand = xl.parse(SHEET_CAND)
    log(f"NCG {len(lit)} literature rows, {len(sysp)} gene rows, "
        f"{len(known)} known TS/oncogenes, {len(cand)} candidates")

    # ---- the per-gene literature count (the bait-count proxy) -------
    lit["pubmed_id"] = lit["pubmed_id"].astype(str)
    n_pub = (lit.groupby("symbol")["pubmed_id"]
             .agg(lambda s: s.nunique()).rename("n_publications"))
    # The cancer type of a gene is spread over ALL of its literature
    # rows: the "Known Cancer" rows carry primary_site "multiple" and
    # an empty cancer_type for exactly the genes whose type we want,
    # so restricting to them loses the information. Aggregate over all
    # rows instead, and treat "multiple" as "pan-cancer" only when no
    # specific site is available for that gene.
    def _join(values) -> str:
        keep = {str(x) for x in values if str(x) not in ("-", "nan", "")}
        return ";".join(sorted(keep))

    def _site(values) -> str:
        keep = {str(x) for x in values if str(x) not in ("-", "nan", "")}
        specific = sorted(keep - {"multiple"})
        return ";".join(specific) if specific else ("multiple" if keep else "")

    site = lit.groupby("symbol")["primary_site"].agg(_site)
    ctype = lit.groupby("symbol")["cancer_type"].agg(_join)

    # ---- the tiers and the roles -------------------------------------
    known_set = set(known["symbol"])
    cand_set = set(cand["symbol"])
    role = known.set_index("symbol")[["cgc_annotation", "vogelstein_annotation",
                                      "NCG6_oncogene", "NCG6_tsg"]]
    role = role.rename(columns={"NCG6_oncogene": "is_oncogene",
                                "NCG6_tsg": "is_tsg"})

    df = sysp.copy()
    df = df.rename(columns={"essential_celllines_percentage":
                             "essential_celllines_pct"})
    df = df.merge(n_pub, left_on="symbol", right_index=True, how="left")
    df = df.merge(site.rename("primary_site"), left_on="symbol",
                   right_index=True, how="left")
    df = df.merge(ctype.rename("cancer_type"), left_on="symbol",
                   right_index=True, how="left")
    df["cancer_type"] = df["cancer_type"].fillna("")
    df["primary_site"] = df["primary_site"].fillna("")
    df = df.merge(role, left_on="symbol", right_index=True, how="left")
    df["is_oncogene"] = df["is_oncogene"].fillna(0).astype(int)
    df["is_tsg"] = df["is_tsg"].fillna(0).astype(int)
    df["tier"] = np.where(df["symbol"].isin(known_set), "known",
                          np.where(df["symbol"].isin(cand_set),
                                   "candidate", "screening"))
    df["n_publications"] = df["n_publications"].fillna(0).astype(int)

    # ---- the HGNC check ----------------------------------------------
    hgnc = pd.read_csv(HGNC, sep="\t", dtype=str)
    hgnc = hgnc[["hgnc_id", "symbol", "entrez_id"]].copy()
    hgnc["entrez_id"] = hgnc["entrez_id"].fillna("")
    id_map = {}
    for _, r in hgnc.iterrows():
        for e in str(r["entrez_id"]).split(";"):
            e = e.strip()
            if e and e.isdigit():
                id_map.setdefault(e, r["hgnc_id"])
    df["hgnc_id"] = df["entrez"].astype(str).map(id_map).fillna("")
    sym_map = dict(zip(hgnc["entrez_id"].fillna(""), hgnc["symbol"]))
    df["symbol_hgnc"] = df["entrez"].astype(str).map(sym_map).fillna("")
    df["symbol_matches_hgnc"] = (df["symbol"] == df["symbol_hgnc"])

    cols = ["symbol", "entrez", "hgnc_id", "symbol_hgnc",
            "symbol_matches_hgnc", "tier", "is_oncogene", "is_tsg",
            "cgc_annotation", "vogelstein_annotation", "primary_site",
            "cancer_type", "n_publications", "essential_celllines_pct",
            "ppin_degree", "ppin_betweenness", "ppin_clustering",
            "duplicability", "origin", "complexes", "mirna",
            "expressed_tissues_rnaseq", "expressed_tissues_protein"]
    df[[c for c in cols if c in df.columns]].to_csv(OUT, sep="\t", index=False)
    log(f"wrote {OUT.relative_to(ROOT)}: {len(df)} genes")

    summary = pd.DataFrame([
        {"quantity": "genes in the NCG catalogue", "value": len(df)},
        {"quantity": "tier=known (711 TS/oncogenes)", "value":
            int((df.tier == "known").sum())},
        {"quantity": "tier=candidate (104 multi-screen)", "value":
            int((df.tier == "candidate").sum())},
        {"quantity": "tier=screening (the rest)", "value":
            int((df.tier == "screening").sum())},
        {"quantity": "oncogenes (NCG6_oncogene=1)", "value":
            int(df.is_oncogene.sum())},
        {"quantity": "tumour suppressors (NCG6_tsg=1)", "value":
            int(df.is_tsg.sum())},
        {"quantity": "genes with a literature count > 0", "value":
            int((df.n_publications > 0).sum())},
        {"quantity": "symbols matching the HGNC record", "value":
            int(df.symbol_matches_hgnc.sum())},
        {"quantity": "symbols NOT matching the HGNC record", "value":
            int((~df.symbol_matches_hgnc).sum())},
        {"quantity": "source", "value":
            "Repana et al. 2019 Genome Biology 20:1, "
            "doi:10.1186/s13059-018-1612-0, Additional file 3"},
    ])
    summary.to_csv(SUMMARY, sep="\t", index=False)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
