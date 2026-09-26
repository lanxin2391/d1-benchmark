"""B2: parse IntOGen driver catalogues into the C3 label-file format.

Inputs (under data/raw/labels/intogen/):
    IntOGen-Drivers-20240920.zip   primary label set
    IntOGen-Drivers-20230531.zip   temporal split (mid)
    IntOGen-Drivers-20200201.zip   temporal split (early)

Outputs (under data/processed/labels/):
    intogen2024.tsv
    intogen2023.tsv
    intogen2020.tsv
    intogen_temporal_new.tsv       positives = genes added in 2024 vs 2020

Each output follows contract C3: header `gene rank source release`. Rank = 1 for
the strongest driver evidence (used by B6 for common-size capping), ties broken
by symbol.

Rank definition (decision D-RANK, recorded at S1):
    rank is by descending count of distinct CANCER_TYPEs in which the gene
    appears as a driver; ties broken by ascending SYMBOL. This is the simplest
    publication-friendly proxy for "how widely supported as a driver".
"""
from __future__ import annotations

import argparse
import io
import os
import zipfile
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_INTOGEN = REPO_ROOT / "data" / "raw" / "labels" / "intogen"
OUT_DIR = REPO_ROOT / "data" / "processed" / "labels"

RELEASE_FILES = {
    "2024-09-20": "IntOGen-Drivers-20240920.zip",
    "2023-05-31": "IntOGen-Drivers-20230531.zip",
    "2020-02-01": "IntOGen-Drivers-20200201.zip",
}
RELEASE_LABEL_ID = {
    "2024-09-20": "intogen2024",
    "2023-05-31": "intogen2023",
    "2020-02-01": "intogen2020",
}
COMPENDIUM_NAME = "Compendium_Cancer_Genes.tsv"


def _read_compendium(zip_path: Path) -> pd.DataFrame:
    """Read IntOGen's per-cohort driver table (one row per gene-cohort)."""
    with zipfile.ZipFile(zip_path) as z:
        names = [n for n in z.namelist() if n.endswith(COMPENDIUM_NAME)]
        if not names:
            raise FileNotFoundError(
                f"{COMPENDIUM_NAME} not found inside {zip_path.name}"
            )
        with z.open(names[0]) as fh:
            return pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8"),
                               sep="\t", dtype=str, low_memory=False)


def build_label_table(zip_path: Path, release_date: str, mapper,
                      source_tag: str = "intogen") -> pd.DataFrame:
    """Return a C3 label table (gene rank source release) for one release."""
    raw = _read_compendium(zip_path)
    if "SYMBOL" not in raw.columns or "CANCER_TYPE" not in raw.columns:
        raise ValueError(
            f"{zip_path.name} is missing expected columns SYMBOL/CANCER_TYPE; "
            f"got {list(raw.columns)[:6]}"
        )
    # Map outdated symbols through HGNC; drop unmapped / NaN rows.
    sym = mapper.map(raw["SYMBOL"], "symbol")
    raw = raw.assign(symbol=sym).dropna(subset=["symbol"])

    # Rank = number of distinct cancer types the gene is a driver in;
    # ties broken by symbol ascending (so the rank order is stable).
    grp = raw.groupby("symbol")
    counts = grp["CANCER_TYPE"].nunique().rename("n_cancer_types")
    symbols = grp.size().rename("_n")  # placeholder
    df = pd.concat([counts, symbols], axis=1).reset_index()
    df = df.sort_values(
        ["n_cancer_types", "symbol"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)
    out = pd.DataFrame({
        "gene": df["symbol"].astype(str),
        "rank": df["rank"].astype(int),
        "source": source_tag,
        "release": release_date,
    })
    return out


def write_label(df: pd.DataFrame, label_id: str, out_dir: Path = OUT_DIR) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{label_id}.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def build_temporal_new(genes_2024: set[str], genes_2020: set[str]) -> pd.DataFrame:
    """Control 6: training set = 2020 drivers; new positives = 2024 \\ 2020."""
    new_genes = sorted(genes_2024 - genes_2020)
    df = pd.DataFrame({
        "gene": new_genes,
        "rank": range(1, len(new_genes) + 1),
        "source": "intogen_temporal_new",
        "release": "2024-09-20 minus 2020-02-01",
    })
    return df


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hgnc", type=Path,
                    default=REPO_ROOT / "data" / "raw" / "reference" / "hgnc_complete_set.txt",
                    help="HGNC complete set file used for symbol harmonisation")
    args = ap.parse_args(argv)

    from d1.hgnc import HGNCMapper
    mapper = HGNCMapper(args.hgnc)

    summary = []
    per_release_genes: dict[str, set[str]] = {}
    for release_date, zip_name in RELEASE_FILES.items():
        zip_path = RAW_INTOGEN / zip_name
        if not zip_path.exists():
            print(f"!! {zip_path} missing; skipping release {release_date}")
            continue
        label_id = RELEASE_LABEL_ID[release_date]
        df = build_label_table(zip_path, release_date, mapper)
        path = write_label(df, label_id)
        n = len(df)
        per_release_genes[release_date] = set(df.gene)
        unmapped = max(0, int(round(n * 0)))  # placeholder
        summary.append((release_date, label_id, n, str(path.relative_to(REPO_ROOT))))
        print(f"[intogen {release_date}] wrote {path.relative_to(REPO_ROOT)} "
              f"({n} unique genes after HGNC harmonisation)")

    if {"2024-09-20", "2020-02-01"}.issubset(per_release_genes):
        new_df = build_temporal_new(per_release_genes["2024-09-20"],
                                    per_release_genes["2020-02-01"])
        new_path = write_label(new_df, "intogen_temporal_new")
        print(f"[temporal] wrote {new_path.relative_to(REPO_ROOT)} "
              f"({len(new_df)} genes added between 2020 and 2024)")

    print("\nSummary:")
    for r, lid, n, p in summary:
        print(f"  {r}  {lid:18s}  n={n:5d}  {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
