"""Cache the Open Targets cancer disease IDs (MONDO_0045024 descendants).

Used by ``scripts/build_opentargets.py`` to filter association tables to
oncology-relevant rows.

The cancer ontology root is ``MONDO_0045024`` ("cancer or benign tumor"). We
include every disease whose ``therapeuticAreas`` list contains this ID, plus
every disease whose ``ancestors`` list contains it (catches diseases whose
direct therapeuticAreas is empty but whose lineage passes through cancer).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pyarrow.parquet as pq

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW = REPO_ROOT / "data" / "raw" / "labels" / "opentargets_26.06"
OUT = REPO_ROOT / "data" / "processed" / "labels" / "_cancer_disease_ids.json"

ONC_ID = "MONDO_0045024"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    disease_file = RAW / "disease" / "disease.parquet"
    if not disease_file.exists():
        print(f"!! {disease_file} not found", file=sys.stderr)
        return 1

    print(f"Reading {disease_file.name} ...")
    df = pq.read_table(str(disease_file),
                       columns=["id", "name", "therapeuticAreas", "ancestors"]).to_pandas()
    print(f"  total diseases: {len(df):,}")

    areas = df["therapeuticAreas"].apply(lambda x: list(x) if x is not None else [])
    anc = df["ancestors"].apply(lambda x: list(x) if x is not None else [])
    is_cancer = (areas.apply(lambda x: ONC_ID in x) |
                 anc.apply(lambda x: ONC_ID in x))
    n = int(is_cancer.sum())
    cancer_ids = sorted(df.loc[is_cancer, "id"].tolist())
    print(f"  cancer diseases (MONDO_0045024): {n:,}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(cancer_ids, f)
    print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
