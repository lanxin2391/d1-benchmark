"""B-side build for the ClinGen cancer label set (clingen_label.tsv).

Source: ``data/raw/labels/clingen/clingen_gene_validity.csv`` — ClinGen's
official Gene Validity Curation dump, 3680 (gene, disease) rows. Columns:
GENE SYMBOL, GENE ID (HGNC), DISEASE LABEL, DISEASE ID (MONDO), CLASSIFICATION, ...

Filter:
* Cancer diseases only: disease ID must be in ``_cancer_disease_ids.json``
  (3,661 IDs cached from ``build_disease.py``; the format there is
  ``MONDO_xxxx`` with underscore, while ClinGen uses ``MONDO:xxxx`` with
  colon — we normalise).
* Classification in {Definitive, Strong, Moderate}: excludes Limited,
  Disputed, Refuted, No Known Disease Relationship. "Definitive"/"Strong"
  are the two ClinGen tiers that other curators treat as solid evidence;
  "Moderate" we include to keep the label non-trivial in size.
* Gene must map through HGNCMapper; previous/alias symbols are accepted
  only when not ambiguous.

Rank definition (D-RANK, replicated from build_intogen.py):
* score = (# cancer types in which the gene is a ClinGen-Definitive-or-Strong
  driver). For moderate we still count, but as a fraction 0.5 — equivalent
  to saying "half a vote".
* ties broken by ascending gene symbol (D-07).
* rank starts at 1.

Output: ``data/processed/labels/clingen_label.tsv`` (C3 schema).
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_CLINGEN = REPO_ROOT / "data" / "raw" / "labels" / "clingen" / "clingen_gene_validity.csv"
HGNC_PATH = REPO_ROOT / "data" / "raw" / "reference" / "hgnc_complete_set.txt"
CANCER_IDS = REPO_ROOT / "data" / "processed" / "labels" / "_cancer_disease_ids.json"
OUT_PATH = REPO_ROOT / "data" / "processed" / "labels" / "clingen_label.tsv"

CLASSIFICATION_WEIGHT = {
    "Definitive": 1.0,
    "Strong": 1.0,
    "Moderate": 0.5,
    "Limited": 0.0,
    "Disputed": 0.0,
    "Refuted": 0.0,
    "No Known Disease Relationship": 0.0,
}


def _normalise_disease_id(s: str) -> str:
    """ClinGen uses ``MONDO:0013212``; OT/our cache uses ``MONDO_0013212``."""
    return s.replace(":", "_")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clingen", type=Path, default=RAW_CLINGEN)
    ap.add_argument("--hgnc", type=Path, default=HGNC_PATH)
    ap.add_argument("--cancer-ids", type=Path, default=CANCER_IDS)
    ap.add_argument("--out", type=Path, default=OUT_PATH)
    args = ap.parse_args(argv)

    sys.path.insert(0, str(REPO_ROOT))
    from d1.hgnc import HGNCMapper

    # ---- 1. load ClinGen
    print(f"Reading {args.clingen.name} ...")
    raw = pd.read_csv(args.clingen, skiprows=4)
    print(f"  raw rows: {len(raw):,}")
    raw = raw[raw["GENE SYMBOL"].notna() & (raw["GENE SYMBOL"] != "+++++++++++")]
    print(f"  after removing header divider: {len(raw):,}")

    raw["DISEASE_NORM"] = raw["DISEASE ID (MONDO)"].astype(str).map(_normalise_disease_id)

    # ---- 2. cancer filter
    with open(args.cancer_ids) as f:
        cancer = set(json.load(f))
    print(f"  cancer disease IDs: {len(cancer):,}")

    cancer_rows = raw[raw["DISEASE_NORM"].isin(cancer)].copy()
    print(f"  rows in cancer diseases: {len(cancer_rows):,}")

    # ---- 3. classification filter
    cancer_rows = cancer_rows.assign(
        weight=cancer_rows["CLASSIFICATION"].map(CLASSIFICATION_WEIGHT).fillna(0.0)
    )
    kept = cancer_rows[cancer_rows["weight"] > 0].copy()
    print(f"  rows after dropping Limited/Disputed/Refuted/NoDisease: {len(kept):,}")

    # ---- 4. HGNC map symbol -> approved symbol
    print(f"Loading HGNC mapper from {args.hgnc.name} ...")
    mapper = HGNCMapper(args.hgnc)
    mapped = mapper.map(kept["GENE SYMBOL"].astype(str), id_type="symbol")
    kept = kept.assign(approved=mapped)
    kept = kept.dropna(subset=["approved"])
    print(f"  rows after HGNC mapping: {len(kept):,}")
    print(f"  unique approved symbols: {kept['approved'].nunique():,}")

    # ---- 5. aggregate score (# weighted cancer types per approved symbol)
    score_per_gene = defaultdict(float)
    cancer_types_per_gene = defaultdict(set)
    for _, r in kept.iterrows():
        g = r["approved"]
        score_per_gene[g] += r["weight"]
        cancer_types_per_gene[g].add(r["DISEASE_NORM"])

    df = pd.DataFrame({
        "gene": list(score_per_gene.keys()),
        "score": [score_per_gene[g] for g in score_per_gene],
        "n_cancer_types": [len(cancer_types_per_gene[g]) for g in score_per_gene],
    })
    df = df.sort_values(
        ["score", "gene"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)
    df["rank"] = df.index + 1
    df["source"] = "clingen_cancer_validity"
    df["release"] = "2026-09-26 (ClinGen + MONDO_0045024 cancer subset)"

    out_df = df[["gene", "rank", "source", "release"]]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"\nWrote {args.out.relative_to(REPO_ROOT)} ({len(out_df):,} rows)")
    print(f"  top 10 by score:")
    print(out_df.head(10).to_string(index=False))
    print(f"  median score: {df['score'].median():.2f}")
    print(f"  n with score >= 2: {(df['score'] >= 2).sum()}")
    print(f"  n with score >= 1: {(df['score'] >= 1).sum()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
