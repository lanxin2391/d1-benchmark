"""C7, part 1 -- the quality tier of every positive gene in every label.

Each catalogue carries its own evidence, and the evidence is not equally
strong inside a catalogue. This script gives every positive gene a tier
(high / mid / low) using each catalogue's own notion of confidence:

* **ClinGen** -- the review status (star level). The gene's best level
  over its disease entries wins: Definitive/Strong = high, Moderate =
  mid, Limited/Disputed = low.
* **IntOGen** -- the number of *distinct detecting methods* that flagged
  the gene (IntOGen combines seven complementary methods and votes).
  Four or more votes = high, two or three = mid, one = low.
* **Open Targets** -- how many independent datasources support the
  gene's cancer associations. A gene whose associations are mostly
  carried by a single datasource (and, at worst, by the literature
  channel alone) is a low-confidence positive.
* **NCG** -- the number of publications the catalogue cites for the
  gene (NCG 6.0's own confidence level is unavailable).

Output: results/tables/label_quality_tiers.tsv
"""
from __future__ import annotations

import json
import sys
import time
import zipfile
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

RAW_LABELS = ROOT / "data/raw/labels"
PROC_LABELS = ROOT / "data/processed/labels"
OUT = ROOT / "results/tables/label_quality_tiers.tsv"
CANCER_IDS = PROC_LABELS / "_cancer_disease_ids.json"

LABELS = {
    "intogen2024": "intogen",
    "intogen_temporal_new": "intogen",
    "ot_all": "opentargets",
    "ot_nolit": "opentargets",
    "clingen_label": "clingen",
    "ncg_known": "ncg",
}
CLINGEN_STARS = {"Definitive": 3, "Strong": 2, "Moderate": 1,
                  "Limited": 0, "Disputed": 0}


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def clingen_tiers() -> dict:
    """gene -> (star, tier) from the ClinGen gene-validity CSV."""
    csv = RAW_LABELS / "clingen/clingen_gene_validity.csv"
    # the file has a 6-line preamble and 10 columns; read it by position
    df = pd.read_csv(csv, header=None, skiprows=6).iloc[:, :8]
    df.columns = ["gene", "hgnc", "disease", "mondo", "mode", "sop",
                  "classification", "url"]
    df = df[df.gene.notna() & df.classification.notna()].copy()
    df["stars"] = df.classification.map(CLINGEN_STARS).fillna(0)
    best = df.sort_values("stars").groupby("gene").tail(1)
    out = {}
    for r in best.itertuples():
        tier = "high" if r.stars >= 2 else ("mid" if r.stars == 1 else "low")
        out[str(r.gene)] = (str(r.classification), tier, int(r.stars))
    log(f"ClinGen: {len(out)} genes, "
        f"{sum(1 for v in out.values() if v[1] == 'high')} high")
    return out


def intogen_tiers(release_zip: Path) -> dict:
    """gene -> (n_methods, tier) from the IntOGen compendium."""
    z = zipfile.ZipFile(release_zip)
    name = [n for n in z.namelist()
            if n.endswith("Compendium_Cancer_Genes.tsv")][0]
    df = pd.read_csv(z.open(name), sep="\t", low_memory=False)
    df = df[df.IS_DRIVER.astype(str).str.lower() == "true"].copy()
    methods = {}
    for r in df.itertuples():
        g = str(r.SYMBOL)
        ms = methods.setdefault(g, set())
        for m in str(r.METHODS).split(","):
            m = m.strip()
            if m and m.lower() != "nan":
                ms.add(m)
    out = {}
    for g, ms in methods.items():
        n = len(ms)
        tier = "high" if n >= 4 else ("mid" if n >= 2 else "low")
        out[g] = (n, tier, sorted(ms))
    log(f"IntOGen ({release_zip.name}): {len(out)} genes, "
        f"{sum(1 for v in out.values() if v[1] == 'high')} high")
    return out


def opentargets_tiers() -> dict:
    """gene -> (frac_literature_only, tier, note) from the datasource table.

    The datasource id lives in ``aggregationValue`` (the constant
    ``aggregationType`` is ``datasourceId`` on every row). Almost every
    (target, cancer disease) pair carries exactly ONE datasource -- 88,510
    of 89,111 in the first part -- so the number of sources per
    association cannot separate the labels. What does separate them is
    **whether the evidence is literature**: a gene whose cancer
    associations are all carried by ``europepmc`` /
    ``uniprot_literature`` is a literature-only positive, and that is the
    tier that matters here (it is the label form of the bait-matched
    null in §F.12).

    Tier: high = at least half of the gene's cancer associations have a
    non-literature source, mid = some but less than half, low = none
    (literature only).
    """
    folder = RAW_LABELS / "opentargets_26.06/association_by_datasource_direct"
    cancer_ids = set(json.loads(CANCER_IDS.read_text()))
    log(f"OT: {len(cancer_ids)} cancer diseases")
    lit_sources = {"europepmc", "uniprot_literature"}
    per_target = {}
    for p in sorted(folder.glob("part-*.snappy.parquet")):
        try:
            d = pq.read_table(str(p), columns=["targetId", "diseaseId",
                                                "aggregationValue",
                                                "associationScore"]).to_pandas()
        except Exception as e:
            log(f"  skipped {p.name}: {e}")
            continue
        d = d[d.diseaseId.isin(cancer_ids)]
        d = d[d.associationScore.notna()]
        if d.empty:
            continue
        d["src"] = d.aggregationValue.astype(str).str.strip().str.lower()
        g = d.groupby(["targetId", "diseaseId"])["src"].agg(
            lambda s: tuple(sorted(set(s))))
        for (t, _dis), srcs in g.items():
            rec = per_target.setdefault(t, {"n": 0, "lit_only": 0,
                                            "nonlit": 0, "n_sources": 0})
            rec["n"] += 1
            rec["n_sources"] += len(srcs)
            if all(s in lit_sources for s in srcs):
                rec["lit_only"] += 1
            else:
                rec["nonlit"] += 1
    log(f"OT: {len(per_target)} targets with cancer associations")

    # the target id -> symbol map
    tfolder = RAW_LABELS / "opentargets_26.06/target"
    tdf = pq.read_table([str(p) for p in sorted(tfolder.glob("part-*.parquet"))],
                        columns=["id", "approvedSymbol"]).to_pandas()
    tdf = tdf.dropna(subset=["id"]).drop_duplicates("id", keep="first")
    tdf["ensembl"] = tdf["id"].astype(str).str.split(".").str[0]
    tmap = dict(zip(tdf.ensembl, tdf.approvedSymbol.astype(str)))

    out = {}
    for t, rec in per_target.items():
        sym = tmap.get(t)
        if sym is None:
            continue
        frac_single = rec["lit_only"] / rec["n"] if rec["n"] else 0.0
        out[sym] = (frac_single, rec)
    log(f"OT: {len(out)} targets mapped to a symbol")
    return out


def ncg_tiers() -> dict:
    df = pd.read_csv(PROC_LABELS / "ncg.tsv", sep="\t", dtype=str).fillna("")
    df = df[df.tier == "known"].copy()
    df["n_pub"] = pd.to_numeric(df.n_publications, errors="coerce").fillna(0)
    out = {}
    for r in df.itertuples():
        n = int(r.n_pub)
        tier = "high" if n >= 10 else ("mid" if n >= 2 else "low")
        out[str(r.symbol_hgnc or r.symbol)] = (n, tier)
    log(f"NCG: {len(out)} genes, "
        f"{sum(1 for v in out.values() if v[1] == 'high')} high")
    return out


def main() -> int:
    rows = []
    clin = clingen_tiers()
    into = intogen_tiers(RAW_LABELS / "intogen/IntOGen-Drivers-20240920.zip")
    ncg = ncg_tiers()
    ot = opentargets_tiers()

    for label, kind in LABELS.items():
        p = PROC_LABELS / f"{label}.tsv"
        if not p.exists():
            log(f"skip {label}: no label file")
            continue
        lab = pd.read_csv(p, sep="\t", dtype=str)
        for r in lab.itertuples():
            g = str(r.gene)
            if kind == "clingen":
                info = clin.get(g)
                tier = info[1] if info else "untiered"
                detail = info[0] if info else ""
                extra = info[2] if info else ""
            elif kind == "intogen":
                info = into.get(g)
                tier = info[1] if info else "untiered"
                detail = f"{info[0]} methods" if info else ""
                extra = ";".join(info[2]) if info else ""
            elif kind == "ncg":
                info = ncg.get(g)
                tier = info[1] if info else "untiered"
                detail = f"{info[0]} publications" if info else ""
                extra = ""
            else:
                info = ot.get(g)
                if info is None:
                    tier, detail, extra = "untiered", "", ""
                else:
                    frac, rec = info
                    tier = ("high" if rec["nonlit"] >= 0.5 * rec["n"] else
                            "mid" if rec["nonlit"] > 0 else "low")
                    detail = f"{frac:.0%} literature-only"
                    extra = (f"nonlit {rec['nonlit']}/{rec['n']}, "
                             f"mean {rec['n_sources'] / rec['n']:.2f} src")
            rows.append({"label": label, "gene": g, "tier": tier,
                         "detail": detail, "evidence": extra})
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False, lineterminator="\n")
    log(f"wrote {OUT.name}: {len(df)} rows")
    print(df.groupby(["label", "tier"]).size().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
