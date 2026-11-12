"""GR-9 -- the arbitration test, stratified by cancer type.

B's review (section 7.2, GR-9) asks whether the prior's advantage
depends on the cancer type. The NCG catalogue annotates each known
driver with the cancer types it was curated for (426 of the 708 known
drivers carry at least one), so the GR-7 arbitration test can be run
separately per cancer type: the arbitration positives are the NCG-known
drivers of that type that our label does not contain, and the background
is the rest of the test fold.

The background is sub-sampled to BG genes with a fixed seed (the
background is only there to define the ROC curve, and 20,000 points
already pin the AUROC to about +/-0.001). The full-background numbers
are in GR-7; `background_effect.tsv` re-computes one arm both ways so
the difference is on record.

Outputs:
  results/tables/gr9_cancer_type_arbitration.tsv
  results/tables/gr9_cancer_type_summary.tsv
  results/tables/gr9_background_effect.tsv
  docs/gr9_cancer_type_arbitration.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_pauc_top10 import auroc  # noqa: E402

SCORES = ROOT / "results" / "scores_v40"
TABLES = ROOT / "results" / "tables"
NCG = ROOT / "data/processed/labels/ncg.tsv"
DOC = ROOT / "docs/gr9_cancer_type_arbitration.md"
ALPHA = 0.5
BG = 20000
SEED = 20261105
MIN_DRIVERS = 25

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> int:
    ncg = pd.read_csv(NCG, sep="\t", dtype=str).fillna("")
    known = ncg[ncg["tier"] == "known"]
    types: dict[str, set[str]] = {}
    for r in known.itertuples():
        for t in str(r.cancer_type).split(";"):
            t = t.strip()
            if t and t not in ("pan-cancer_adult", "pan-cancer_paediatric"):
                types.setdefault(t, set()).add(r.symbol)
    # keep the types with enough drivers to be worth a test
    keep = sorted([(t, g) for t, g in types.items()
                   if len(g) >= MIN_DRIVERS], key=lambda x: -len(x[1]))
    log(f"{len(types)} cancer types annotated; {len(keep)} with >= "
        f"{MIN_DRIVERS} known drivers")
    for t, g in keep:
        log(f"   {t}: {len(g)} drivers")

    rng = np.random.default_rng(SEED)
    rows = []
    for net in NETWORKS:
        for label in LABELS:
            fp = SCORES / f"{net}__{label}.npz"
            if not fp.exists():
                continue
            d = np.load(fp, allow_pickle=True)
            ai = int(np.argmin(np.abs(np.asarray(d["alphas"], dtype=float)
                                      - ALPHA)))
            rwr = np.asarray(d["rwr"][ai], dtype=float)
            deg = np.asarray(d["degree"], dtype=float)
            in_test = np.asarray(d["in_test"], dtype=bool)
            y = np.asarray(d["y"]).astype(int)
            genes = np.asarray(d["genes"]).astype(str)
            gset = {g: i for i, g in enumerate(genes)}
            ours = y == 1
            for ctype, drivers in keep:
                sel_all = np.array([gset[g] for g in drivers
                                    if g in gset and not ours[gset[g]]],
                                   dtype=int)
                if len(sel_all) < 5:
                    continue
                deltas = []
                for j in range(rwr.shape[1]):
                    m = in_test[j]
                    idx = np.where(m)[0]
                    if len(idx) == 0:
                        continue
                    pos = np.intersect1d(sel_all, idx, assume_unique=False)
                    if len(pos) < 5:
                        continue
                    rest = np.setdiff1d(idx, pos, assume_unique=True)
                    if len(rest) > BG:
                        rest = rng.choice(rest, size=BG, replace=False)
                    yy = np.concatenate([np.ones(len(pos), dtype=int),
                                         np.zeros(len(rest), dtype=int)])
                    sc = np.concatenate([rwr[pos, j], rwr[rest, j]])
                    dc = np.concatenate([deg[pos], deg[rest]])
                    deltas.append(auroc(sc, yy) - auroc(dc, yy))
                if not deltas:
                    continue
                rows.append({
                    "network": net, "label": label, "cancer_type": ctype,
                    "n_drivers": int(len(sel_all)),
                    "n_drivers_in_test": int(len(sel_all)),
                    "delta_auroc": float(np.mean(deltas)),
                    "delta_auroc_sd": float(np.std(deltas, ddof=1)),
                    "n_folds": len(deltas),
                })
            log(f"  {net} x {label}: {len(keep)} cancer types done")
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "gr9_cancer_type_arbitration.tsv", sep="\t",
              index=False)

    summ = (df.groupby("cancer_type")
            .agg(n_arms=("delta_auroc", "size"),
                 median_delta_auroc=("delta_auroc", "median"),
                 min_delta_auroc=("delta_auroc", "min"),
                 max_delta_auroc=("delta_auroc", "max"),
                 n_positive=("delta_auroc", lambda s: int((s > 0).sum())),
                 drivers=("n_drivers", "median"))
            .sort_values("median_delta_auroc", ascending=False)
            .reset_index())
    summ.to_csv(TABLES / "gr9_cancer_type_summary.tsv", sep="\t", index=False)

    # the background-subsampling effect, on one arm
    eff = background_effect()
    eff.to_csv(TABLES / "gr9_background_effect.tsv", sep="\t", index=False)

    write_doc(df, summ, eff)
    log(f"wrote the GR-9 tables ({len(df)} arm x cancer-type rows)")
    return 0


def background_effect() -> pd.DataFrame:
    """Re-compute one arm with the full background and with the
    sub-sampled background, to show that the sub-sampling is harmless."""
    rng = np.random.default_rng(SEED)
    net, label = "string_phys700", "intogen2024"
    ncg = pd.read_csv(NCG, sep="\t", dtype=str).fillna("")
    known = set(ncg[(ncg["tier"] == "known")]["symbol"])
    d = np.load(SCORES / f"{net}__{label}.npz", allow_pickle=True)
    ai = int(np.argmin(np.abs(np.asarray(d["alphas"], dtype=float) - ALPHA)))
    rwr = np.asarray(d["rwr"][ai], dtype=float)
    deg = np.asarray(d["degree"], dtype=float)
    in_test = np.asarray(d["in_test"], dtype=bool)
    y = np.asarray(d["y"]).astype(int)
    genes = np.asarray(d["genes"]).astype(str)
    is_known = np.array([g in known for g in genes])
    ours = y == 1
    rows = []
    for bg in (BG, 200000, 10 ** 9):
        ds = []
        for j in range(rwr.shape[1]):
            idx = np.where(in_test[j])[0]
            pos = idx[is_known[idx] & ~ours[idx]]
            rest = idx[~is_known[idx]]
            if len(pos) < 5 or len(rest) < 50:
                continue
            if len(rest) > bg:
                rest = rng.choice(rest, size=bg, replace=False)
            yy = np.concatenate([np.ones(len(pos), dtype=int),
                                 np.zeros(len(rest), dtype=int)])
            ds.append(auroc(np.concatenate([rwr[pos, j], rwr[rest, j]]), yy)
                       - auroc(np.concatenate([deg[pos], deg[rest]]), yy))
        rows.append({"background": ("full" if bg > 10 ** 8 else str(bg)),
                     "delta_auroc": float(np.mean(ds)),
                     "n_folds": len(ds)})
    return pd.DataFrame(rows)


def write_doc(df: pd.DataFrame, summ: pd.DataFrame, eff: pd.DataFrame) -> None:
    lines = [
        "# GR-9 -- the arbitration test, stratified by cancer type",
        "",
        "**Date:** 2026-11-05  ",
        "**Script:** `scripts/run_gr9_cancer_type_arbitration.py`  ",
        "**Tables:** `results/tables/gr9_cancer_type_arbitration.tsv`, "
        "`gr9_cancer_type_summary.tsv`, `gr9_background_effect.tsv`",
        "",
        "## The question",
        "",
        "Does the prior's advantage over the degree baseline depend on "
        "the cancer type? The arbitration positives are the NCG-known "
        "drivers of one cancer type that our label does not contain; "
        "the background is the rest of the test fold.",
        "",
        "## Result",
        "",
        "| cancer type | drivers | median dAUROC | range | arms positive |",
        "|---|---|---|---|---|",
    ]
    for r in summ.itertuples():
        lines.append(f"| {r.cancer_type} | {r.drivers:.0f} | "
                     f"{r.median_delta_auroc:+.4f} | "
                     f"{r.min_delta_auroc:+.4f} .. {r.max_delta_auroc:+.4f} | "
                     f"{r.n_positive}/{r.n_arms} |")
    lines += [
        "",
        "## The background sub-sampling",
        "",
        "The background is sub-sampled to 20,000 genes (fixed seed) to "
        "keep the 24 arms x cancer types tractable. The table below "
        "re-computes one arm with the full background, so the effect is "
        "on record.",
        "",
        "| background | dAUROC |",
        "|---|---|",
    ]
    for r in eff.itertuples():
        lines.append(f"| {r.background} | {r.delta_auroc:+.4f} |")
    lines += [
        "",
        "## Caveats",
        "",
        "1. The NCG cancer-type annotation is itself curated and "
        "incomplete (426 of the 708 known drivers carry a type), so the "
        "strata are not exhaustive.",
        "2. The pan-cancer entries are excluded from the stratification "
        "on purpose: they would dominate every stratum.",
        "3. **Exploratory** (amendment 4); H1 untouched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
