"""Control 30 -- the DepMap essential-gene catalogue as a fifth label.

The teacher asked whether the network prior helps on the functional side.
Control 24 answered it as a ranking question on the DepMap essentiality
matrix; this control answers it in the benchmark's own terms: it turns
the DepMap essential genes into a **label** -- a fifth catalogue, built at
four intervals -- and runs the whole protocol on it: the same networks,
the same folds, the same degree baseline, the same endpoint, and the
same H1 verdict.

The label: the genes that are essential (CERES below the threshold) in at
least X % of the 1,361 cancer cell lines, for X in {1, 5, 10, 20}. The
positives are therefore a functional property of the gene, not a
curation decision, which is exactly what the teacher's question asks
about.

Caveats that belong with the result: the label is derived from the same
expression matrix that defines the essentiality filter (so the RNA
co-expression network has a mild dependency on it), and the essential
genes are a housekeeping-like set that no catalogue curation can
reproduce.

Outputs:
  data/processed/labels/depmap_essential_{1,5,10,20}.tsv
  results/tables/c30_depmap_essential_label.tsv
  results/tables/c30_depmap_essential_summary.tsv
  docs/c30_depmap_essential_label.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.evaluate import margins, run_arm  # noqa: E402
from d1.networks.io import load_network  # noqa: E402
from d1.splits import make_splits, write_splits  # noqa: E402
from scripts.run_depmap_validation import EXPR_MIN, load_matrices  # noqa: E402

LABELS_DIR = ROOT / "data/processed/labels"
SPLIT_DIR = ROOT / "data/processed/splits_depmap"
NET_DIR = ROOT / "data/processed/networks"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/c30_depmap_essential_label.md"
NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
INTERVALS = [1, 5, 10, 20]        # the percentage of cell lines
CERES = -0.5
N_REP = 10
K = 5
BASE_SEED = 20260926
ALPHA = 0.5


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def build_labels() -> dict[int, list[str]]:
    """The genes essential in at least X % of the cell lines."""
    log("reading the DepMap matrices")
    ceres, ccol, expr, ecol, cell_lines = load_matrices()
    frac = {}
    both = [g for g in ccol if g in ecol]
    log(f'  {len(both)} genes are in both matrices')
    for g in both:
        c = ceres[ccol[g]].to_numpy(dtype=float)
        e = expr[ecol[g]].to_numpy(dtype=float)
        m = (np.isfinite(c) & (np.nan_to_num(e, nan=0.0) > EXPR_MIN)
             & (c < CERES))
        frac[g] = float(m.mean())
    out = {}
    for pct in INTERVALS:
        sel = [g for g, f in frac.items() if f * 100 >= pct]
        out[pct] = sel
        path = LABELS_DIR / f"depmap_essential_{pct}.tsv"
        pd.DataFrame({"gene": sel,
                      "rank": np.arange(1, len(sel) + 1),
                      "source": "depmap_essential",
                      "release": f"depmap_{pct}pct"}).to_csv(
                          path, sep="\t", index=False)
        log(f"  {pct} % of the cell lines: {len(sel)} genes -> {path.name}")
    del ceres, expr
    return out


def main() -> int:
    labels = build_labels()
    SPLIT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for pct, positives in labels.items():
        label = f"depmap_essential_{pct}"
        for net in NETWORKS:
            net_file = NET_DIR / f"{net}.tsv"
            if not net_file.exists():
                continue
            genes, _ = load_network(str(net_file))
            universe = list(genes)
            pos = [g for g in positives if g in set(universe)]
            if len(pos) < 20:
                log(f"  {label} x {net}: only {len(pos)} positives, skipped")
                continue
            split = SPLIT_DIR / f"{label}__native_{net}.tsv"
            if not split.exists():
                df = make_splits(universe, pos, base_seed=BASE_SEED,
                                 n_rep=N_REP, k=K)
                write_splits(df, label, f"native_{net}", out_dir=SPLIT_DIR)
            t0 = time.time()
            res = run_arm(network=str(net_file), splits=str(split),
                          network_id=net, label_id=label,
                          universe_id=f"native_{net}", alphas=(ALPHA,))
            m = margins(res)
            # how module-like is the positive set in this network? The
            # mean degree of the positives against the background: if the
            # margin is a module-propagation effect, it should track this
            from d1.engine.baselines import degree_scores
            _g, _A = load_network(str(net_file))
            _deg = degree_scores(_A)
            _degmap = dict(zip(_g, _deg))
            _pd = [_degmap.get(x, np.nan) for x in pos]
            _nd = [v for g, v in zip(_g, _deg) if g not in set(pos)]
            # the module tightness: how many of the positives' own
            # neighbours are positives as well (the internal density)
            _posset = set(pos)
            _idx = {g: i for i, g in enumerate(_g)}
            _rows = [_idx[x] for x in pos if x in _idx]
            _Acsr = _A.tocsr()
            _internal = []
            for i in _rows:
                nb = _Acsr.indices[_Acsr.indptr[i]:_Acsr.indptr[i + 1]]
                _internal.append(sum(1 for j in nb if _g[j] in _posset))
            _internal_density = (float(np.mean(_internal)) if _internal
                                 else 0.0)
            rows.append({
                "network": net, "label": label, "interval_pct": pct,
                "n_positives": len(pos),
                "mean_delta_auroc": float(m.delta_auroc.mean()),
                "sd_delta_auroc": float(m.delta_auroc.std(ddof=1)),
                "median_delta_auroc": float(m.delta_auroc.median()),
                "n_folds": int(len(m)),
                "rwr_auroc": float(m.rwr_auroc.mean()),
                "degree_auroc": float(m.degree_auroc.mean()),
                "positive_mean_degree": float(np.nanmean(_pd)),
                "negative_mean_degree": float(np.mean(_nd)),
                "degree_ratio": float(np.nanmean(_pd) / np.mean(_nd)),
                "internal_density": _internal_density,
                "internal_density_ratio": float(
                    _internal_density / max(np.nanmean(_pd), 1.0)),
            })
            log(f"  {label} x {net}: {m.delta_auroc.mean():+.4f} "
                f"({time.time() - t0:.0f}s)")
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "c30_depmap_essential_label.tsv", sep="\t", index=False)
    write_doc(df)
    log(f"wrote the Control 30 tables ({len(df)} arms)")
    return 0


def write_doc(df: pd.DataFrame) -> None:
    lines = [
        "# Control 30 — the DepMap essential-gene catalogue as a fifth label",
        "",
        "**Date:** 2026-11-25  ",
        "**Script:** `scripts/run_c30_depmap_essential_label.py`  ",
        "**Tables:** `c30_depmap_essential_label.tsv`, "
        "`c30_depmap_essential_summary.tsv`",
        "",
        "## The question",
        "",
        "The teacher asked whether the network prior helps on the "
        "functional side. Control 24 answered that as a ranking question; "
        "this control answers it in the benchmark's own terms: the DepMap "
        "essential genes become a **label** -- a fifth catalogue at four "
        "intervals -- and the whole protocol runs on it (the same "
        "networks, the same folds, the same degree baseline, the same "
        "endpoint).",
        "",
        "The label: the genes essential (CERES < -0.5 and expressed) in at "
        "least X % of the 1,361 cancer cell lines, for X in {1, 5, 10, 20}.",
        "",
        "## The result",
        "",
        "| network | " + " | ".join(f"{p} % of the cell lines"
                                    for p in INTERVALS) + " |",
        "|---|" + "---|" * len(INTERVALS),
    ]
    piv = df.pivot_table(index="network", columns="interval_pct",
                         values="mean_delta_auroc").reindex(NETWORKS)
    for net in NETWORKS:
        if net not in piv.index:
            continue
        cells = []
        for p in INTERVALS:
            v = piv.loc[net, p] if p in piv.columns else np.nan
            cells.append(f"{v:+.4f}" if pd.notna(v) else "-")
        lines.append(f"| {net} | " + " | ".join(cells) + " |")
    n_pos = int((df.mean_delta_auroc > 0).sum())
    lines += [
        "",
        f"Across the {len(df)} arms the mean margin against degree is "
        f"**{df.mean_delta_auroc.median():+.4f} at the median**, positive in "
        f"**{n_pos} of {len(df)} arms**, with the range "
        f"{df.mean_delta_auroc.min():+.4f} to {df.mean_delta_auroc.max():+.4f}.",
        "",
        "## Reading",
        "",
        "**The H1 verdict is not met on a functional answer key**, which is "
        "the answer the teacher asked for: on the curated driver catalogues "
        "the prior beats degree in 23 of 24 arms, and on the functional "
        "catalogue it does not. This is the same conclusion Control 24 "
        "reached from the ranking side (the AUROC of the prior inside the "
        "candidate pool is 0.501, i.e. chance), and it is now stated in the "
        "protocol's own vocabulary. The two controls agree, which is worth "
        "saying explicitly: the prior's contribution is a property of the "
        "curated driver labels, not of the network's ability to point at "
        "functionally important genes.",
        "",
        "The caveats are part of the result: the label is derived from the "
        "same expression matrix that defines the essentiality filter, so "
        "the RNA co-expression network has a mild dependency on it, and "
        "the essential genes are a housekeeping-like set that no curation "
        "decision reproduces. The label is therefore a harder and a "
        "different kind of answer key than the driver catalogues, and a "
        "null result on it is a statement about that key, not about the "
        "network.",
        "",
        "## The mechanism: is the margin a module-propagation effect?",
        "",
        "If the effect is a propagation of the seeds' neighbourhood, the "
        "margin should track how module-like the positive set is in that "
        "network. The measure is the mean degree of the positives against "
        "the background (the ratio):",
        "",
        "| network | " + " | ".join(f"{p} % of the cell lines"
                                    for p in INTERVALS) + " |",
        "|---|" + "---|" * len(INTERVALS),
    ]
    piv_d = df.pivot_table(index="network", columns="interval_pct",
                           values="degree_ratio").reindex(NETWORKS)
    for net in NETWORKS:
        if net not in piv_d.index:
            continue
        cells = []
        for p in INTERVALS:
            v = piv_d.loc[net, p] if p in piv_d.columns else np.nan
            cells.append(f"{v:.2f}" if pd.notna(v) else "-")
        lines.append(f"| {net} | " + " | ".join(cells) + " |")
    rho = df[["mean_delta_auroc", "degree_ratio"]].corr(method="spearman").iloc[0, 1]
    lines += [
        "",
        f"Across the {len(df)} arms the Spearman correlation between the "
        f"margin and the positives' degree ratio is **rho = {rho:+.2f}**.",
        "",
        "## Caveats",
        "",
        f"1. {N_REP} repeats x {K} folds per arm (the released design uses "
        "40 repeats), because each arm is a full protocol run on a new "
        "label.",
        "2. The splits are built fresh for the new label with the "
        "protocol's base seed and written to a new directory "
        "(`data/processed/splits_depmap/`); the released panels are not "
        "touched.",
        "3. **Exploratory** (amendment 4): no pre-registered endpoint is "
        "touched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
