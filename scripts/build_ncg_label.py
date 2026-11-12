"""Build the NCG label in the C3 schema, so it can be used as a label set.

The NCG catalogue itself is `data/processed/labels/ncg.tsv` (2,372
genes, built by `scripts/build_ncg.py`). This script turns its
**known** tier (the high-confidence tumour suppressors and oncogenes)
into a label file with the C3 contract (gene, rank, source, release),
the same shape as the other four labels, so the split builder, the
runner and the statistics can treat it like any other label.

Decisions
---------
* **Which tier is the label?** The `known` tier (708 genes). The 104
  `candidate` genes are multi-screen candidates and are kept aside for
  the label-quality tiering (plan item C7); the 1,560 `screening` genes
  are not drivers at all and are never used as positives.
* **Symbols (D-02).** The HGNC approved symbol is used, with the
  previous symbol mapped to its current one where the catalogue carries
  one (9 of the 708: CARS -> CARS1, MKL1 -> MRTFA, FGFR1OP -> CEP43,
  ...). `build_ncg.py` already resolved those against the HGNC
  complete set, so the mapping is taken from there and not redone.
* **Rank (D-13, adapted).** NCG 6.0's own confidence level is not
  available (the portal is down for maintenance), so the rank is the
  number of publications NCG cites for the gene, descending, ties broken
  by the symbol ascending. This matters beyond cosmetics: the D-11 cap
  keeps the top-N by rank, so the capped set is the most-studied
  (and by the C2 check the *least* degree-loaded of the NCG genes).

Outputs:
  data/processed/labels/ncg_known.tsv
  results/tables/ncg_label_build_summary.tsv
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CATALOGUE = ROOT / "data/processed/labels/ncg.tsv"
OUT = ROOT / "data/processed/labels/ncg_known.tsv"
SUMMARY = ROOT / "results/tables/ncg_label_build_summary.tsv"
RELEASE = "2019-08-01"


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> int:
    ncg = pd.read_csv(CATALOGUE, sep="\t", dtype=str).fillna("")
    known = ncg[ncg["tier"] == "known"].copy()
    known["n_publications"] = pd.to_numeric(known["n_publications"],
                                             errors="coerce").fillna(0)
    # the HGNC approved symbol (the previous symbol is mapped to it)
    known["gene"] = known["symbol_hgnc"].where(known["symbol_hgnc"] != "",
                                                known["symbol"])
    known = known.drop_duplicates(subset="gene")
    # D-13 (adapted): the literature count, descending; ties by symbol
    known = known.sort_values(["n_publications", "gene"],
                               ascending=[False, True], kind="mergesort")
    known["rank"] = range(1, len(known) + 1)
    out = known[["gene", "rank", "n_publications"]].copy()
    out["source"] = "ncg"
    out["release"] = RELEASE
    out = out[["gene", "rank", "source", "release", "n_publications"]]
    OUT.write_text(out.to_csv(sep="\t", index=False, lineterminator="\n"),
                   encoding="utf-8")
    log(f"wrote {OUT.name}: {len(out)} positives "
        f"(tier 'known', {len(ncg) - len(known)} other rows kept aside)")

    renamed = known[known["symbol"] != known["gene"]]
    summary = pd.DataFrame([
        {"item": "n_known", "value": len(known)},
        {"item": "n_candidate_tier", "value": int((ncg.tier == "candidate").sum())},
        {"item": "n_screening_tier", "value": int((ncg.tier == "screening").sum())},
        {"item": "n_symbol_renamed_by_hgnc", "value": len(renamed)},
        {"item": "n_publications_min", "value": int(out.n_publications.min())},
        {"item": "n_publications_median",
         "value": float(out.n_publications.median())},
        {"item": "n_publications_max", "value": int(out.n_publications.max())},
        {"item": "rank_rule", "value": "n_publications desc, gene asc"},
        {"item": "release", "value": RELEASE},
        {"item": "source", "value": "ncg"},
    ])
    summary.to_csv(SUMMARY, sep="\t", index=False)
    log("wrote the build summary")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
