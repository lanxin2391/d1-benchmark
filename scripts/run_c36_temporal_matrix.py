"""Control 36 -- the temporal transfer matrix (plan item 19, redefined).

IntOGen has no 2021 or 2022 release, so the plan's "per-year" split
cannot be built from it. What the data do allow is a **four-generation
chain** (2019-11-12, 2020-02-01, 2023-05-31, 2024-09-20), which is the
finest temporal resolution the catalogue actually has. The question this
control answers is the one the defence asked about: *does the network
prior still help when the positives it is scored against were curated
years after the seeds it was seeded from?*

For every (train year, test year) pair with train < test:

    seeds        the drivers known in the train-year release
    test set     the drivers present in the test-year release but absent
                 from the train-year release (the newly curated genes)
    negatives    the network's genes that are not in the test-year release
    endpoint     delta-AUROC of the prior (RWR, alpha = 0.5) against the
                 degree baseline, on the same test genes

Each cell is repeated with a bootstrap draw of the seed set (the same
device as the D-09 temporal control, which fixes the "the repeats were
identical" defect).

Outputs:
  results/tables/c36_temporal_matrix.tsv
  results/tables/c36_temporal_summary.tsv
  docs/c36_temporal_matrix.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.baselines import degree_scores  # noqa: E402
from d1.engine.evaluate import evaluate  # noqa: E402
from d1.engine.rwr import rwr, seed_matrix, transition_T  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETS = ROOT / "data/processed/networks"
LABELS_DIR = ROOT / "data/processed/labels"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/c36_temporal_matrix.md"
NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
RELEASES = {"2019": "intogen2019", "2020": "intogen2020",
            "2023": "intogen2023", "2024": "intogen2024"}
# the release dates, for the gap in years
DATES = {"2019": (2019, 11, 12), "2020": (2020, 2, 1),
         "2023": (2023, 5, 31), "2024": (2024, 9, 20)}
PAIRS = [("2019", "2020"), ("2019", "2023"), ("2019", "2024"),
         ("2020", "2023"), ("2020", "2024"), ("2023", "2024")]
ALPHA = 0.5
N_REP = 10
SEED = 20261126


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def label_genes(label_id: str) -> set[str]:
    p = LABELS_DIR / f"{label_id}.tsv"
    if not p.exists():
        return set()
    return set(pd.read_csv(p, sep="\t", dtype=str)["gene"])


def gap_years(a: str, b: str) -> float:
    ya, ma, da = DATES[a]
    yb, mb, db = DATES[b]
    return (yb - ya) + (mb - ma) / 12.0 + (db - da) / 365.0


def main() -> int:
    genes_by_year = {y: label_genes(l) for y, l in RELEASES.items()}
    for y, g in genes_by_year.items():
        log(f"  {y}: {len(g)} drivers")
    rows = []
    for net in NETWORKS:
        net_file = NETS / f"{net}.tsv"
        if not net_file.exists():
            continue
        genes, A = load_network(str(net_file))
        deg = degree_scores(A)
        WT = transition_T(A)
        gset = set(genes)
        index = {g: i for i, g in enumerate(genes)}
        for tr, te in PAIRS:
            train = genes_by_year[tr] & gset
            test_new = (genes_by_year[te] - genes_by_year[tr]) & gset
            negatives = gset - genes_by_year[te]
            if len(train) < 20 or len(test_new) < 10 or len(negatives) < 100:
                log(f"  {net} {tr}->{te}: too small, skipped "
                    f"(the seeds {len(train)}, the new {len(test_new)})")
                continue
            seed_idx = np.array([index[g] for g in sorted(train)])
            pos_idx = np.array([index[g] for g in sorted(test_new)])
            neg_idx = np.array([index[g] for g in sorted(negatives)])
            rng = np.random.default_rng(SEED + hash((net, tr, te)) % 9973)
            t0 = time.time()
            for rep in range(N_REP):
                draw = rng.choice(seed_idx, size=len(seed_idx), replace=True)
                P0 = seed_matrix(len(genes), [draw])
                P, _ = rwr(WT, P0, alpha=ALPHA)
                sc = P[:, 0]
                pos = sc[pos_idx]
                neg = sc[neg_idx]
                labels = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
                values = np.concatenate([pos, neg])
                ev_p = evaluate(values, labels, np.arange(len(values)))
                d_pos = deg[pos_idx]
                d_neg = deg[neg_idx]
                d_values = np.concatenate([d_pos, d_neg])
                ev_d = evaluate(d_values, labels, np.arange(len(d_values)))
                rows.append({
                    "network": net, "train": tr, "test": te,
                    "gap_years": round(gap_years(tr, te), 3),
                    "repeat": rep, "n_seeds": len(seed_idx),
                    "n_new": len(pos_idx),
                    "prior_auroc": ev_p["auroc"],
                    "degree_auroc": ev_d["auroc"],
                    "delta_auroc": ev_p["auroc"] - ev_d["auroc"],
                })
            log(f"  {net} {tr}->{te}: {N_REP} bootstraps "
                f"({time.time() - t0:.0f}s)")
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "c36_temporal_matrix.tsv", sep="\t", index=False)
    cell = df.groupby(["network", "train", "test", "gap_years"]).agg(
        delta=("delta_auroc", "mean"), sd=("delta_auroc", "std"),
        prior=("prior_auroc", "mean"), degree=("degree_auroc", "mean"),
        n_new=("n_new", "first"), n=("delta_auroc", "size")).reset_index()
    cell.to_csv(TABLES / "c36_temporal_summary.tsv", sep="\t", index=False)
    write_doc(cell)
    log(f"wrote the Control 36 tables ({len(df)} rows)")
    return 0


def write_doc(cell: pd.DataFrame) -> None:
    piv = cell.pivot_table(index="network", columns=["train", "test"],
                            values="delta")
    lines = [
        "# Control 36 — the temporal transfer matrix",
        "",
        "**Date:** 2026-11-26  ",
        "**Script:** `scripts/run_c36_temporal_matrix.py`  ",
        "**Tables:** `c36_temporal_matrix.tsv`, `c36_temporal_summary.tsv`",
        "",
        "## The design",
        "",
        "IntOGen has no 2021 or 2022 release, so the plan's per-year split "
        "cannot be built from it; what the catalogue offers is a "
        "four-generation chain (2019-11-12, 2020-02-01, 2023-05-31, "
        "2024-09-20). For every (train year, test year) pair the seeds are "
        "the drivers known in the train-year release, the test positives "
        "are the drivers present in the test-year release but absent from "
        "the train-year release, and the endpoint is the delta-AUROC of "
        "the prior against the degree baseline on the same test genes "
        f"({N_REP} bootstrap draws of the seed set per cell).",
        "",
        "## The mean delta-AUROC against degree, per network and pair",
        "",
        "| network | " + " | ".join(f"{a}->{b}" for a, b in PAIRS) + " |",
        "|---|" + "---|" * len(PAIRS),
    ]
    for net in NETWORKS:
        if net not in piv.index.get_level_values(0):
            continue
        cells = []
        for a, b in PAIRS:
            try:
                v = piv.loc[net, (a, b)]
                v = float(v.iloc[0]) if hasattr(v, "iloc") else float(v)
            except KeyError:
                v = np.nan
            cells.append(f"{v:+.4f}" if pd.notna(v) else "-")
        lines.append(f"| {net} | " + " | ".join(cells) + " |")
    lines += [
        "",
        "## The decay curve (the median over the six networks)",
        "",
        "| the pair | the gap in years | the median margin | the arms "
        "positive |",
        "|---|---|---|---|",
    ]
    for a, b in PAIRS:
        sub = cell[(cell.train == a) & (cell.test == b)]
        if not len(sub):
            continue
        lines.append(f"| {a}->{b} | {gap_years(a, b):.2f} | "
                     f"{sub.delta.median():+.4f} | "
                     f"{int((sub.delta > 0).sum())}/{len(sub)} |")
    lines += [
        "",
        "## Reading",
        "",
        "## Caveats",
        "",
        "1. The 2019->2020 pair is 83 days apart, so it is a consistency "
        "check rather than a decay step (the two catalogues are almost the "
        "same curation).",
        "2. The test positives are the newly curated genes, which are few "
        "in each cell (the count is in the table); the AUROC on a few "
        "positives is noisier than the main endpoint, which is why each "
        "cell is repeated with bootstrap seed draws.",
        "3. The degree baseline is the same baseline as the primary "
        "endpoint, so the numbers are comparable in sign and scale with "
        "Table 1.",
        "4. **Exploratory** (amendment 4): no pre-registered endpoint is "
        "touched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
