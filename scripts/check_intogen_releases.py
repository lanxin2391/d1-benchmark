#!/usr/bin/env python3
"""Check the IntOGen release zips and report what is needed for the
per-year temporal split (plan item 19).

What it does
-------------
1. Lists every ``IntOGen-Drivers-*.zip`` in ``data/raw/labels/intogen/``.
2. For each release, opens the per-cohort driver table
   (``Compendium_Cancer_Genes.tsv`` -- the file the label builder reads)
   and reports the schema, the gene count and the cancer-type count.
3. Compares the gene set of each release with the ones before it, so the
   year-over-year additions are visible (this is what the per-year split
   needs).
4. Checks the releases we already have against the labels we already
   built (``data/processed/labels/intogen20*.tsv``) so a new download can
   be trusted the same way.

Usage
-----
    python scripts/check_intogen_releases.py
"""
from __future__ import annotations

import io
import sys
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "labels" / "intogen"
LABELS = ROOT / "data" / "processed" / "labels"
COMPENDIUM = "Compendium_Cancer_Genes.tsv"


def release_key(name: str) -> str:
    """The date inside IntOGen-Drivers-YYYYMMDD.zip."""
    stem = name.replace("IntOGen-Drivers-", "").replace(".zip", "")
    if len(stem) == 8 and stem.isdigit():
        return f"{stem[:4]}-{stem[4:6]}-{stem[6:]}"
    return stem


def main() -> int:
    zips = sorted(RAW.glob("IntOGen-Drivers-*.zip"),
                  key=lambda p: release_key(p.name))
    if not zips:
        print(f"no IntOGen zips in {RAW}")
        return 1
    print(f"{len(zips)} IntOGen release(s) in {RAW}")
    print()
    sets: dict[str, set] = {}
    rows = []
    for z in zips:
        rel = release_key(z.name)
        with zipfile.ZipFile(z) as zf:
            names = zf.namelist()
            comp = [n for n in names if n.endswith(COMPENDIUM)]
            if not comp:
                print(f"  {z.name}: {COMPENDIUM} is missing -- the zip looks "
                      "wrong (an HTML error page?)")
                continue
            raw = pd.read_csv(io.TextIOWrapper(zf.open(comp[0]),
                                                encoding="utf-8"),
                               sep="\t", dtype=str, low_memory=False)
        sym = raw["SYMBOL"].dropna().unique()
        sets[rel] = set(sym)
        n_ct = raw.groupby("SYMBOL")["CANCER_TYPE"].nunique()
        rows.append({
            "release": rel, "file": z.name,
            "kb": z.stat().st_size // 1024,
            "rows": len(raw), "genes": len(sets[rel]),
            "cancer_types": int(raw["CANCER_TYPE"].nunique()),
            "cohorts": int(raw["COHORT"].nunique()),
            "top_rank": int(n_ct.max()),
        })
        print(f"  {z.name}")
        print(f"      the table: {len(raw):,} rows, {len(raw.columns)} columns")
        print(f"      the genes: {len(sets[rel]):,}; the cancer types: "
              f"{raw['CANCER_TYPE'].nunique()}; the cohorts: "
              f"{raw['COHORT'].nunique()}")
        print(f"      the columns: {', '.join(list(raw.columns)[:10])}"
              f"{' ...' if len(raw.columns) > 10 else ''}")
    df = pd.DataFrame(rows)
    print()
    print("the year-over-year additions (the input of the per-year split):")
    rels = sorted(sets)
    for a, b in zip(rels, rels[1:]):
        added = sets[b] - sets[a]
        lost = sets[a] - sets[b]
        print(f"  {a} -> {b}: +{len(added):,} genes, -{len(lost):,} genes")
        if added:
            sample = sorted(added)[:8]
            print(f"      e.g. {', '.join(sample)}")
    if len(rels) < 4:
        print()
        print("  NOTE: the per-year split wants one release per year "
              "(2020, 2021, 2022, 2023, 2024). The releases that are "
              "missing have to be downloaded from "
              "https://www.intogen.org/download and put into this folder.")
    print()
    print("the labels we have already built (the consistency check):")
    for rel in rels:
        lab = LABELS / f"intogen{rel[:4]}.tsv"
        if not lab.exists():
            continue
        genes = set(pd.read_csv(lab, sep="\t")["gene"])
        inter = len(genes & sets[rel])
        print(f"  {lab.name}: {len(genes):,} genes; the overlap with the "
              f"release: {inter:,} "
              f"({inter / max(len(genes), 1):.1%} of the label)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
