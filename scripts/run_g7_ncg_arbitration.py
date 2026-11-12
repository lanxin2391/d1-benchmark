"""G7 -- the NCG catalogue as an untouched arbitration set.

Boutros et al. (2014): the gold standard dominates any ranking, so the
evaluator should be blind to the algorithm. The NCG catalogue is
independent of the project's four label sets (it shares no curation
with them), which makes it the arbitration set we do not control.

The question: with the project's own label as the seed set, does the
prior rank the drivers we DID NOT label above the rest of the universe?
If it does, the prior's ranking is not an artefact of the label we
seeded it with.

Per (network, label) arm and per fold:
    P = the NCG-known drivers (708) that are NOT in this label's
        positive set and that are in the network;
    B = the rest of the test fold;
    the AUROC of the RWR score and of the degree baseline on (P, B), and
    their difference.

The seeds are the label's own training-fold positives, so this measures
transfer: the prior is asked to rank drivers it was never trained on.

Outputs:
  results/tables/g7_ncg_arbitration.tsv
  docs/g7_ncg_arbitration.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_pauc_top10 import auroc, pauc_area, pauc_excess  # noqa: E402

SCORES = ROOT / "results" / "scores_v40"
TABLES = ROOT / "results" / "tables"
NCG = ROOT / "data/processed/labels/ncg.tsv"
DOC = ROOT / "docs/g7_ncg_arbitration.md"
ALPHA = 0.5

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> int:
    ncg = pd.read_csv(NCG, sep="\t", dtype=str)
    known = set(ncg[ncg["tier"] == "known"]["symbol"])
    log(f"the NCG-known tier: {len(known)} genes")

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
            is_known = np.array([g in known for g in genes])
            is_ours = y == 1

            auroc_d, pauc_d, np_pos = [], [], []
            for j in range(rwr.shape[1]):
                m = in_test[j]
                if m.sum() == 0:
                    continue
                idx = np.where(m)[0]
                # the arbitration positives: NCG-known, not in our label
                sel = idx[is_known[idx] & ~is_ours[idx]]
                rest = idx[~is_known[idx]]
                if len(sel) < 5 or len(rest) < 50:
                    continue
                yy = np.concatenate([np.ones(len(sel), dtype=int),
                                      np.zeros(len(rest), dtype=int)])
                sc = np.concatenate([rwr[sel, j], rwr[rest, j]])
                dc = np.concatenate([deg[sel], deg[rest]])
                auroc_d.append(auroc(sc, yy) - auroc(dc, yy))
                pauc_d.append(pauc_excess(pauc_area(sc, yy),
                                          pauc_area(dc, yy)))
                np_pos.append(len(sel))
            if not auroc_d:
                continue
            rows.append({
                "network": net, "label": label,
                "n_arbitration_positives_median": int(np.median(np_pos)),
                "delta_auroc": float(np.mean(auroc_d)),
                "delta_auroc_sd": float(np.std(auroc_d, ddof=1)),
                "delta_pauc": float(np.mean(pauc_d)),
                "n_folds": len(auroc_d),
            })
            log(f"  {net} x {label}: dAUROC {rows[-1]['delta_auroc']:+.4f} "
                f"(NCG-known, not ours: {rows[-1]['n_arbitration_positives_median']} "
                f"per fold)")
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "g7_ncg_arbitration.tsv", sep="\t", index=False)
    write_doc(df)
    log(f"wrote g7_ncg_arbitration.tsv ({len(df)} arms)")
    return 0


def write_doc(df: pd.DataFrame) -> None:
    lines = [
        "# G7 -- the NCG catalogue as an untouched arbitration set",
        "",
        "**Date:** 2026-11-05  ",
        "**Script:** `scripts/run_g7_ncg_arbitration.py`  ",
        "**Table:** `results/tables/g7_ncg_arbitration.tsv`",
        "",
        "## The question",
        "",
        "The seeds are the label's own training-fold positives, so this "
        "is a **transfer** test: can the prior rank the drivers it was "
        "never seeded with? The arbitration positives are the NCG-known "
        "drivers (708) that are **not** in this label's positive set; the "
        "background is the rest of the test fold. The margin is again "
        "the prior over the degree baseline.",
        "",
        "## Result",
        "",
        "| network | label | arbitration positives per fold | dAUROC | "
        "dpAUC |",
        "|---|---|---|---|---|",
    ]
    for r in df.itertuples():
        lines.append(f"| {r.network} | {r.label} | "
                     f"{r.n_arbitration_positives_median} | "
                     f"{r.delta_auroc:+.4f} | {r.delta_pauc:+.4f} |")
    lines += [
        "",
        "## Caveats",
        "",
        "1. The arbitration positives are a *subset* of the genes that are "
        "not in our label, so they are not random: they are curated "
        "drivers that our label misses. The test asks whether the prior "
        "finds those, which is exactly what an arbitration catalogue is "
        "for (Boutros et al. 2014: the gold standard dominates the "
        "ranking).",
        "2. The seeds are our label's positives, so this is not a "
        "label-free evaluation: it measures transfer from our label to "
        "an independent catalogue, not the quality of a label-free "
        "prior.",
        "3. **Exploratory** (amendment 4): it does not enter the H1 "
        "decision.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
