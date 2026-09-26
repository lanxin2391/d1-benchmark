"""B4: Open Targets 26.06 label sets (ot_all, ot_nolit).

Two label sets are produced, both in the contract-C3 schema ``gene rank source
release``:

* **ot_all**: per (target, disease) take ``association_overall_direct``'s
  ``associationScore``, restrict to cancer diseases (MONDO_0045024 descendants),
  per target take the max score across cancer diseases, rank descending, keep
  the top N (D-11; default 600 to match the protocol's "~600 genes").
* **ot_nolit**: same but using ``association_by_datasource_direct``, dropping
  the ``europepmc`` (literature-only) datasource, and recombining the remaining
  scores per (target, disease) with the harmonic-sum formula documented below.
  This is control 7 / H3 in label form: the network prior's contribution net
  of literature evidence channels.

Decisions recorded at S1 / W2 sync point:
* D-OT-AREA: cancer ontology = ``MONDO_0045024`` ("cancer or benign tumor"),
  expanded by ``ancestors`` so descendant cancers are included.
* D-OT-NOLIT-SCORE: harmonic sum ``H(s) = (s_1/1^2 + s_2/2^2 + ...) / (pi^2 / 6)``
  over the per-(target, disease, datasource) scores after dropping europepmc.
  The denominator is the Basel problem value (Harmonic series on squares).
* D-OT-CAP: top N = 600 (protocol §1.1); refine at S2 once the shared-universe
  size N is known.

Notes
-----
* ``target`` parquet provides ``approvedSymbol`` directly. We cross-check it
  with the HGNC mapper over ``ensemblId`` and keep the HGNC result when it
  disagrees, because protocol C1 says HGNC is authoritative.
* Cancer disease IDs are cached in ``data/processed/labels/_cancer_disease_ids.json``
  (built by ``build_disease.py``).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW = REPO_ROOT / "data" / "raw" / "labels" / "opentargets_26.06"
PROCESSED = REPO_ROOT / "data" / "processed" / "labels"
PROCESSED.mkdir(parents=True, exist_ok=True)
CANCER_IDS_JSON = PROCESSED / "_cancer_disease_ids.json"

ONC_ID = "MONDO_0045024"
LITERATURE_DATASOURCE = "europepmc"
TOP_N = 600
SOURCE_TAG = "open_targets_26.06"


# --------------------------------------------------------------- helpers
def _load_cancer_ids() -> set[str]:
    if not CANCER_IDS_JSON.exists():
        raise FileNotFoundError(
            f"{CANCER_IDS_JSON} not found; run scripts/build_disease.py first"
        )
    with open(CANCER_IDS_JSON) as f:
        return set(json.load(f))


def _load_target_ensembl_to_symbol(target_dir: Path, mapper) -> pd.Series:
    """Return a Series ``{ensembl_id: approved_symbol}`` for every target row.

    The HGNC mapper is the source of truth for symbol spelling; we keep the
    Open Targets ``approvedSymbol`` as a fallback only when the HGNC mapper
    returns NaN (rare).
    """
    parts = sorted(target_dir.glob("part-*.snappy.parquet"))
    if not parts:
        raise FileNotFoundError(f"no target parquet parts in {target_dir}")
    df = pq.read_table([str(p) for p in parts],
                       columns=["id", "approvedSymbol"]).to_pandas()
    df = df.dropna(subset=["id"]).drop_duplicates("id", keep="first")
    df["ensembl_clean"] = df["id"].str.split(".").str[0]
    mapped = mapper.map(df["ensembl_clean"], "ensembl")
    symbols = mapped.where(mapped.notna(), df["approvedSymbol"])
    return pd.Series(symbols.values, index=df["ensembl_clean"].values, name="symbol")


def _read_overall(folder: Path, cancer_ids: set[str]) -> pd.DataFrame:
    """Read all parts of ``association_overall_direct``, restrict to cancer diseases.

    Skips parts that are not yet complete (parquet magic bytes fail to read);
    this lets the pipeline produce a partial result while downloads are still
    in progress. ``scripts/verify_downloads.py`` records the remaining size.

    Returns a DataFrame ``[targetId (Ensembl), associationScore]``.
    """
    parts = sorted(folder.glob("part-*.snappy.parquet"))
    if not parts:
        raise FileNotFoundError(f"no overall parquet parts in {folder}")
    dfs = []
    for p in parts:
        try:
            d = pq.read_table(str(p),
                              columns=["targetId", "diseaseId", "associationScore"]).to_pandas()
        except Exception as e:
            print(f"  [skip partial] {p.name}: {type(e).__name__}")
            continue
        d = d[d.diseaseId.isin(cancer_ids)]
        if not d.empty:
            dfs.append(d)
    return pd.concat(dfs, ignore_index=True)


def _harmonise_target_ids(df: pd.DataFrame, target_ids: pd.Series) -> pd.DataFrame:
    """Map ``targetId`` (Ensembl) -> HGNC symbol using the precomputed map."""
    out = df.copy()
    out["ensembl_clean"] = out["targetId"].astype(str).str.split(".").str[0]
    out["gene"] = out["ensembl_clean"].map(target_ids)
    out = out.dropna(subset=["gene"])
    return out[["gene", "associationScore"]]


def _to_label_table(scored: pd.DataFrame, top_n: int, label_id: str,
                    release: str, source_tag: str = SOURCE_TAG) -> pd.DataFrame:
    """Take a DataFrame ``[gene, score]`` and return a C3 label table."""
    per_target = (scored.groupby("gene", as_index=False)["associationScore"]
                  .max().sort_values(["associationScore", "gene"],
                                     ascending=[False, True], kind="mergesort"))
    per_target = per_target.head(top_n).reset_index(drop=True)
    per_target["rank"] = range(1, len(per_target) + 1)
    return pd.DataFrame({
        "gene": per_target["gene"].astype(str),
        "rank": per_target["rank"].astype(int),
        "source": source_tag,
        "release": release,
    })


# --------------------------------------------------------------- build_ot_all
def build_ot_all(target_ids: pd.Series, cancer_ids: set[str], top_n: int = TOP_N) -> Path:
    folder = RAW / "association_overall_direct"
    print(f"[ot_all] reading {folder.name} ...")
    df = _read_overall(folder, cancer_ids)
    print(f"[ot_all] {len(df):,} (target, cancer-disease) rows")
    df = _harmonise_target_ids(df, target_ids)
    print(f"[ot_all] {df.gene.nunique():,} unique target genes after HGNC mapping")
    label = _to_label_table(df, top_n=top_n, label_id="ot_all",
                            release="26.06 (overall)")
    out = PROCESSED / "ot_all.tsv"
    label.to_csv(out, sep="\t", index=False)
    print(f"[ot_all] wrote {out.relative_to(REPO_ROOT)} (top {top_n} genes)")
    return out


# --------------------------------------------------------------- build_ot_nolit
def _harmonic_sum(scores: np.ndarray) -> float:
    """Normalised harmonic sum ``H(s) = sum_i s_i / i^2  /  (pi^2 / 6)``.

    For a single score this returns the score unchanged (because
    ``s / 1^2 / (pi^2/6) * (pi^2/6)`` is ill-defined; the formula is defined
    for two or more scores; the caller is expected to handle n=1 separately).
    For one score we return the score directly.
    """
    s = np.sort(np.asarray(scores, dtype=float))[::-1]
    if s.size == 1:
        return float(s[0])
    ranks = np.arange(1, s.size + 1, dtype=float)
    return float((s / ranks ** 2).sum() / (np.pi ** 2 / 6.0))


def build_ot_nolit(target_ids: pd.Series, cancer_ids: set[str],
                   top_n: int = TOP_N) -> Path | None:
    folder = RAW / "association_by_datasource_direct"
    if not folder.exists() or not list(folder.glob("part-*.snappy.parquet")):
        print(f"[ot_nolit] SKIP: {folder} has no parquet parts (run download_D1_data.py)")
        return None
    print(f"[ot_nolit] reading {folder.name} ...")
    parts = sorted(folder.glob("part-*.snappy.parquet"))
    rows = []
    for p in parts:
        d = pq.read_table(str(p),
                          columns=["targetId", "diseaseId", "datasourceId", "score"]).to_pandas()
        d = d[d.diseaseId.isin(cancer_ids) & (d.datasourceId != LITERATURE_DATASOURCE)]
        if not d.empty:
            rows.append(d)
    df = pd.concat(rows, ignore_index=True)
    print(f"[ot_nolit] {len(df):,} (target, cancer-disease, datasource) rows after dropping europepmc")
    df = df.assign(ensembl_clean=df["targetId"].astype(str).str.split(".").str[0])
    df["gene"] = df["ensembl_clean"].map(target_ids)
    df = df.dropna(subset=["gene"])
    pair = (df.groupby(["gene", "diseaseId"])["score"]
            .apply(_harmonic_sum).rename("associationScore").reset_index())
    out_df = _to_label_table(pair, top_n=top_n, label_id="ot_nolit",
                             release="26.06 (no-literature harmonic sum)",
                             source_tag="open_targets_26.06_no_lit")
    out = PROCESSED / "ot_nolit.tsv"
    out_df.to_csv(out, sep="\t", index=False)
    print(f"[ot_nolit] wrote {out.relative_to(REPO_ROOT)} (top {top_n} genes)")
    return out


# --------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hgnc", type=Path,
                    default=REPO_ROOT / "data" / "raw" / "reference" / "hgnc_complete_set.txt")
    ap.add_argument("--top-n", type=int, default=TOP_N)
    ap.add_argument("--only", choices=["ot_all", "ot_nolit", "both"], default="both")
    args = ap.parse_args(argv)

    sys.path.insert(0, str(REPO_ROOT))
    from d1.hgnc import HGNCMapper
    print(f"Loading HGNC mapper from {args.hgnc.name} ...")
    mapper = HGNCMapper(args.hgnc)

    print("Loading target parquet parts ...")
    target_ids = _load_target_ensembl_to_symbol(RAW / "target", mapper)
    print(f"  {len(target_ids):,} target rows")

    cancer_ids = _load_cancer_ids()
    print(f"  {len(cancer_ids):,} cancer disease IDs")

    if args.only in ("ot_all", "both"):
        build_ot_all(target_ids, cancer_ids, top_n=args.top_n)
    if args.only in ("ot_nolit", "both"):
        build_ot_nolit(target_ids, cancer_ids, top_n=args.top_n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
