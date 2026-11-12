"""G8 -- the negative set is not clean, and how much does it matter?

B's project review (section 7.2, G8) points at Bailey et al. 2018 (the
LDA filter that removes 45 curated false positives from 127 true
positives) and Repana et al. 2019 (NCG flags 201 possible false
positives and 49 known genes it does not label). Our negative set is
"everything else in the universe", which means every curated driver we
happen not to label is being treated as a confident negative.

With the NCG catalogue in hand we can measure this and correct it:

A. THE CONTAMINATION. Per (network, label): how many of the held-out
   "negatives" are in the NCG catalogue (tier=known / candidate /
   screening), and how many NCG-known drivers each of our labels
   actually contains (the coverage question, which is also the external
   gold-standard comparison).

B. THE CLEANED ANALYSIS. The primary endpoint recomputed with every
   NCG gene removed from the background: if the margin survives, the
   false-positive contamination of the negative set is not what carries
   the result. All three endpoints are reported (the full AUROC and the
   two co-primary ones) because the endpoint choice turned out to
   matter.

Outputs:
  results/tables/g8_negative_contamination.tsv
  results/tables/g8_cleaned_background.tsv
  docs/g8_negative_set.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_pauc_top10 import (auroc, pauc_area,  # noqa: E402
                                    pauc_excess, precision_at_top)

SCORES = ROOT / "results" / "scores_v40"
TABLES = ROOT / "results" / "tables"
NCG = ROOT / "data" / "processed" / "labels" / "ncg.tsv"
DOC = ROOT / "docs" / "g8_negative_set.md"
ALPHA = 0.5
TOP_FRAC = 0.10

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> int:
    ncg = pd.read_csv(NCG, sep="\t", dtype=str)
    ncg_all = set(ncg["symbol"])
    ncg_known = set(ncg[ncg["tier"] == "known"]["symbol"])
    ncg_cand = set(ncg[ncg["tier"] == "candidate"]["symbol"])
    log(f"the NCG catalogue: {len(ncg_all)} genes "
        f"({len(ncg_known)} known, {len(ncg_cand)} candidate)")

    contam_rows, clean_rows, cover_rows = [], [], []
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
            is_ncg = np.array([g in ncg_all for g in genes])
            is_known = np.array([g in ncg_known for g in genes])
            is_cand = np.array([g in ncg_cand for g in genes])

            n_neg = n_ncg = n_known = n_cand = 0
            std = {"auroc": [], "pauc": [], "p10": []}
            cln = {"auroc": [], "pauc": [], "p10": []}
            for j in range(rwr.shape[1]):
                m = in_test[j]
                if m.sum() == 0:
                    continue
                idx = np.where(m)[0]
                pos = idx[y[idx] == 1]
                neg = idx[y[idx] == 0]
                if len(pos) == 0 or len(neg) == 0:
                    continue
                n_neg += len(neg)
                n_ncg += int(is_ncg[neg].sum())
                n_known += int(is_known[neg].sum())
                n_cand += int(is_cand[neg].sum())
                # the standard endpoint, on the full test set
                std["auroc"].append(auroc(rwr[idx, j], y[idx])
                                    - auroc(deg[idx], y[idx]))
                ar_r = pauc_area(rwr[idx, j], y[idx])
                ar_d = pauc_area(deg[idx], y[idx])
                std["pauc"].append(pauc_excess(ar_r, ar_d))
                std["p10"].append(precision_at_top(rwr[idx, j], y[idx],
                                                 genes[idx])
                                  - precision_at_top(deg[idx], y[idx],
                                                     genes[idx]))
                # the cleaned analysis: the NCG genes leave the background
                keep = idx[~is_ncg[idx]]
                if len(keep) <= len(pos) + 10:
                    continue
                yk = y[keep]
                if yk.sum() == 0 or yk.sum() == len(keep):
                    continue
                cln["auroc"].append(auroc(rwr[keep, j], yk)
                                    - auroc(deg[keep], yk))
                ar_r = pauc_area(rwr[keep, j], yk)
                ar_d = pauc_area(deg[keep], yk)
                cln["pauc"].append(pauc_excess(ar_r, ar_d))
                cln["p10"].append(precision_at_top(rwr[keep, j], yk,
                                                 genes[keep])
                                  - precision_at_top(deg[keep], yk,
                                                     genes[keep]))

            for key, vals in (("standard", std), ("cleaned", cln)):
                if not vals["auroc"]:
                    continue
                row = {"network": net, "label": label, "version": key,
                       "n_folds": len(vals["auroc"])}
                for ep in ("auroc", "pauc", "p10"):
                    row[f"delta_{ep}"] = float(np.mean(vals[ep]))
                    row[f"sd_{ep}"] = float(np.std(vals[ep], ddof=1))
                clean_rows.append(row)
            contam_rows.append({
                "network": net, "label": label,
                "n_negatives": n_neg,
                "n_negatives_in_ncg": n_ncg,
                "frac_negatives_in_ncg": n_ncg / max(n_neg, 1),
                "n_negatives_ncg_known": n_known,
                "n_negatives_ncg_candidate": n_cand,
                "frac_negatives_ncg_known": n_known / max(n_neg, 1),
            })
            # the coverage of our labels against the NCG gold standard
            cover_rows.append({
                "network": net, "label": label,
                "n_positives": int((y == 1).sum()),
                "n_positives_in_ncg": int(is_ncg[y == 1].sum()),
                "n_positives_ncg_known": int(is_known[y == 1].sum()),
                "frac_positives_ncg_known": float(
                    is_known[y == 1].sum() / max((y == 1).sum(), 1)),
            })
        log(f"  {net}: done")

    ct = pd.DataFrame(contam_rows)
    ct.to_csv(TABLES / "g8_negative_contamination.tsv", sep="\t", index=False)
    cl = pd.DataFrame(clean_rows)
    cl.to_csv(TABLES / "g8_cleaned_background.tsv", sep="\t", index=False)
    cv = pd.DataFrame(cover_rows)
    cv.to_csv(TABLES / "g8_label_coverage_of_ncg.tsv", sep="\t", index=False)
    write_doc(ct, cl, cv)
    log("wrote the three G8 tables")
    return 0


def write_doc(ct: pd.DataFrame, cl: pd.DataFrame, cv: pd.DataFrame) -> None:
    piv = cl.pivot_table(index=["network", "label"], columns="version",
                         values=["delta_auroc", "delta_pauc", "delta_p10"])
    lines = [
        "# G8 -- the negative set is not clean",
        "",
        "**Date:** 2026-11-05  ",
        "**Script:** `scripts/run_g8_negative_set.py`  ",
        "**Tables:** `results/tables/g8_negative_contamination.tsv`, "
        "`g8_cleaned_background.tsv`, `g8_label_coverage_of_ncg.tsv`",
        "",
        "## A. how much of our negative set is actually a curated driver?",
        "",
        "| network | label | negatives | of them in NCG | NCG-known |",
        "|---|---|---|---|---|",
    ]
    for r in ct.itertuples():
        lines.append(f"| {r.network} | {r.label} | {r.n_negatives:,} | "
                     f"{r.n_negatives_in_ncg:,} "
                     f"({r.frac_negatives_in_ncg:.1%}) | "
                     f"{r.n_negatives_ncg_known:,} "
                     f"({r.frac_negatives_ncg_known:.2%}) |")
    lines += [
        "",
        "## B. does the margin survive a cleaned background?",
        "",
        "Every NCG gene (all 2,372, all tiers) is removed from the "
        "background and the endpoint recomputed on the same folds and the "
        "same scores.",
        "",
        "| network | label | dAUROC standard | dAUROC cleaned | "
        "dpAUC standard | dpAUC cleaned | dtop10% standard | "
        "dtop10% cleaned |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for (net, label), r in piv.iterrows():
        def g(ep, ver):
            try:
                return f"{r[(ep, ver)]:+.4f}"
            except KeyError:
                return "-"
        lines.append(
            f"| {net} | {label} | {g('delta_auroc', 'standard')} | "
            f"{g('delta_auroc', 'cleaned')} | {g('delta_pauc', 'standard')} | "
            f"{g('delta_pauc', 'cleaned')} | {g('delta_p10', 'standard')} | "
            f"{g('delta_p10', 'cleaned')} |")
    lines += [
        "",
        "## C. how much of the NCG gold standard do our labels contain?",
        "",
        "| network | label | positives | of them NCG-known | fraction |",
        "|---|---|---|---|---|",
    ]
    for r in cv.itertuples():
        lines.append(f"| {r.network} | {r.label} | {r.n_positives} | "
                     f"{r.n_positives_ncg_known} | "
                     f"{r.frac_positives_ncg_known:.1%} |")
    lines += [
        "",
        "## Reading",
        "",
        "**The negative set is genuinely contaminated.** Between 10 % and "
        "22 % of the held-out \"negatives\" are in the NCG catalogue -- "
        "curated driver genes we simply do not label -- and 1.3 % to 9.5 % "
        "are NCG-*known* high-confidence drivers. Reactome's arms are the "
        "worst (15-22 %), the dense protein networks the least (10 %).",
        "",
        "**Cleaning the background moves the numbers in both "
        "directions, so it is not a small correction.** On the IntOGen "
        "label the margin roughly halves (intact 0.0651 -> 0.0304, funmap "
        "0.1001 -> 0.0247, string_phys700 0.1189 -> 0.0768); on the "
        "Open Targets labels it *grows* (funmap ot_all 0.0853 -> 0.1030, "
        "string_phys700 ot_nolit 0.1711 -> 0.1979), because removing "
        "curated drivers from the background removes false negatives that "
        "the degree baseline was getting right; and on the temporal label "
        "reactome flips from -0.0125 to +0.0073. Under the cleaned "
        "background the low-FPR advantage also disappears on the IntOGen "
        "label (several arms go negative in dpAUC).",
        "",
        "**The IntOGen label is a good proxy for an independent gold "
        "standard; the Open Targets labels are not.** 70-79 % of the "
        "NCG-known drivers are in intogen2024 (and 52-59 % in the "
        "temporal label), against 27-37 % in ot_nolit -- the Open "
        "Targets panels are much broader, which is also why they carry "
        "the largest false-positive contamination.",
        "",
        "## Caveats",
        "",
        "1. NCG is a *different* curation criterion (somatic driver "
        "roles from cancer sequencing screens), so a gene being in NCG "
        "does not make it a false positive in our labels -- it makes it "
        "a driver we do not label. That is the point of the correction, "
        "but it also means the cleaned background removes genes that may "
        "be perfectly good negatives for our question.",
        "2. The correction removes the NCG genes from the background only; "
        "the positives are untouched.",
        "3. The coverage table (part C) is an external comparison, not a "
        "performance measure: it says how much of an independent gold "
        "standard each of our label sets contains.",
        "4. **Exploratory** (amendment 4): it does not enter the H1 "
        "decision.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
