"""Control 24 -- does the network add anything on the FUNCTIONAL task,
once hubness is taken out?

Q1 (the gene level): for each arm and fold, and for each cell line,
the AUROC of three scores for discriminating the functionally essential
genes inside the candidate pool:
    prior     the RWR score,
    degree    the degree baseline,
    residual  the prior's score after regressing out log1p(degree).
Q2 (the stratum level): split the pool into degree strata and compare
the prior's top-K picks, the degree top-K picks, the residual's top-K
picks and a RANDOM pick from the SAME stratum against the rest of the
pool.

Both are vectorised over the 1,361 cell lines.

Outputs:
  results/tables/c24_depmap_gene_level.tsv
  results/tables/c24_depmap_within_degree.tsv
  docs/c24_depmap_within_degree.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_depmap_validation import (  # noqa: E402
    ALPHA, EXPR_MIN, LABELS, NETWORKS, SCORES, THRESHOLDS, load_matrices)

TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/c24_depmap_within_degree.md"
SEED = 20261116
K_VALUES = [30, 100]
N_STRATA = 5
Q2_EVERY = 5          # one fold in five for the stratified test


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def tie_ranks(x: np.ndarray) -> np.ndarray:
    """The tie-averaged 1-based ranks of a score vector."""
    order = np.argsort(x, kind="mergesort")
    s = x[order]
    ranks = np.empty(len(s), dtype=float)
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        ranks[i:j + 1] = 0.5 * (i + j) + 1.0
        i = j + 1
    out = np.empty(len(x), dtype=float)
    out[order] = ranks
    return out


def auc_per_column(scores: np.ndarray, ess: np.ndarray):
    """The AUROC of one score for every column (cell line) at once."""
    ranks = tie_ranks(scores)
    E = ess[np.argsort(scores, kind="mergesort")]
    n_pos = E.sum(axis=0).astype(float)
    n_neg = E.shape[0] - n_pos
    sum_pos = (ranks[:, None] * E).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        auc = (sum_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    ok = (n_pos >= 3) & (n_neg >= 3)
    return auc, ok


def residualise(score: np.ndarray, degree: np.ndarray) -> np.ndarray:
    """The part of the score that log-degree does not explain."""
    x = np.log1p(degree)
    ok = np.isfinite(score) & np.isfinite(x)
    out = np.zeros_like(score, dtype=float)
    if ok.sum() < 10:
        return out
    A = np.vstack([x[ok], np.ones(ok.sum())]).T
    coef, *_ = np.linalg.lstsq(A, score[ok], rcond=None)
    out[ok] = score[ok] - (coef[0] * x[ok] + coef[1])
    return out


def main() -> int:
    from scipy import stats as sps

    q1_only_gene = "--only-q1" in sys.argv
    q2_only = "--only-q2" in sys.argv

    arms = {}
    need = set()
    for net in NETWORKS:
        for label in LABELS:
            fp = SCORES / f"{net}__{label}.npz"
            if not fp.exists():
                continue
            d = np.load(fp, allow_pickle=True)
            ai = int(np.argmin(np.abs(np.asarray(d["alphas"], dtype=float)
                                      - ALPHA)))
            arms[(net, label)] = (
                np.asarray(d["rwr"][ai], dtype=float),
                np.asarray(d["degree"], dtype=float),
                np.asarray(d["in_test"], dtype=bool),
                np.asarray(d["y"]).astype(int),
                np.asarray(d["genes"]).astype(str))
            for j in range(arms[(net, label)][2].shape[0]):
                idx = np.where(arms[(net, label)][2][j])[0]
                neg = idx[arms[(net, label)][3][idx] == 0]
                if len(neg) == 0:
                    continue
                kmax = max(100, int(np.ceil(0.01 * len(idx))))
                rwr = arms[(net, label)][0]
                order = neg[np.argsort(-rwr[neg, j], kind="mergesort")]
                need.update(arms[(net, label)][4][order[:kmax]].tolist())
    log(f"the candidate pool touches {len(need)} genes")

    ceres, ccol, expr, ecol, cell_lines = load_matrices()
    wanted = [g for g in sorted(need) if g in ccol and g in ecol]
    log(f"  {len(wanted)} of them are in both matrices")
    index = {g: i for i, g in enumerate(wanted)}
    score = np.empty((len(wanted), len(cell_lines)), dtype=np.float32)
    expressed = np.empty((len(wanted), len(cell_lines)), dtype=bool)
    for i, g in enumerate(wanted):
        score[i] = ceres[ccol[g]].to_numpy(dtype=np.float32)
        e = expr[ecol[g]].to_numpy(dtype=float)
        expressed[i] = np.nan_to_num(e, nan=0.0) > EXPR_MIN
    del ceres, expr

    gene_rows = []
    strat_rows = []
    for threshold in THRESHOLDS:
        ess = (np.isfinite(score) & expressed & (score < threshold))
        log(f"=== threshold CERES < {threshold} ===")
        for (net, label), (rwr, deg, in_test, y, genes) in arms.items():
            for j in range(in_test.shape[0]):
                idx = np.where(in_test[j])[0]
                neg = idx[y[idx] == 0]
                if len(neg) < 50:
                    continue
                in_mat = np.array([g in index for g in genes[neg]],
                                  dtype=bool)
                pool = np.array([index[g] for g in genes[neg][in_mat]],
                                dtype=int)
                if len(pool) < 50:
                    continue
                # the scores, restricted to the same genes as the pool
                s_prior = rwr[neg, j][in_mat]
                s_degree = deg[neg][in_mat]
                s_resid = residualise(s_prior, s_degree)
                e_pool = ess[pool]                     # (genes, cell lines)

                # ---- Q1: the gene level ----
                if not q2_only:
                    for tag, sc in (("prior", s_prior),
                                    ("degree", s_degree),
                                    ("residual", s_resid)):
                        auc, ok = auc_per_column(sc, e_pool)
                        if ok.sum() == 0:
                            continue
                        gene_rows.append({
                            "network": net, "label": label,
                            "threshold": threshold, "scorer": tag,
                            "fold": j, "n_pool": len(pool),
                            "n_cell_lines": int(ok.sum()),
                            "auc_median": float(np.median(auc[ok])),
                        })

                # ---- Q2: inside degree strata ----
                if q1_only_gene or (j % Q2_EVERY):
                    continue
                order_deg = np.argsort(s_degree, kind="mergesort")
                strata = np.empty(len(pool), dtype=int)
                strata[order_deg] = np.repeat(
                    np.arange(N_STRATA),
                    int(np.ceil(len(pool) / N_STRATA)))[:len(pool)]
                for st in range(N_STRATA):
                    in_st = np.where(strata == st)[0]
                    if len(in_st) < 10:
                        continue
                    rest = np.setdiff1d(np.arange(len(pool)), in_st)
                    rng = np.random.default_rng(SEED + 1000 * j + st)
                    for tag, sc in (("prior", s_prior),
                                    ("degree", s_degree),
                                    ("residual", s_resid)):
                        for K in K_VALUES:
                            k = min(K, len(in_st))
                            if len(in_st) - k < 5:
                                continue
                            cand = in_st[np.argsort(-sc[in_st],
                                                    kind="mergesort")[:k]]
                            picks = {"top": cand}
                            picks["random"] = rng.permutation(in_st)[:k]
                            for pick_tag, pick in picks.items():
                                rows = pool[pick]
                                bg = pool[rest]
                                if len(bg) < 20:
                                    continue
                                k_cand = ess[rows].sum(axis=0).astype(int)
                                k_bg = ess[bg].sum(axis=0).astype(int)
                                M = len(rows) + len(bg)
                                with np.errstate(divide="ignore",
                                                 invalid="ignore"):
                                    pv = sps.hypergeom.sf(k_cand - 1, M,
                                                          len(rows), k_bg)
                                good = (k_bg > 0) & np.isfinite(pv)
                                if good.sum() == 0:
                                    continue
                                strat_rows.append({
                                    "network": net, "label": label,
                                    "threshold": threshold, "stratum": st,
                                    "scorer": tag, "pick": pick_tag, "K": K,
                                    "n_pool": len(pool),
                                    "n_stratum": len(in_st),
                                    "n_cells": int(good.sum()),
                                    "frac_enriched": float(
                                        (pv[good] < 0.05).mean()),
                                    "mean_ess_cand": float(k_cand.mean()),
                                    "mean_ess_bg": float(k_bg.mean()),
                                })
            log(f"  {net} x {label} done")
    if not q2_only:
        g = pd.DataFrame(gene_rows)
        g.to_csv(TABLES / "c24_depmap_gene_level.tsv", sep="\t",
                 index=False)
    s = pd.DataFrame(strat_rows)
    s.to_csv(TABLES / "c24_depmap_within_degree.tsv", sep="\t", index=False)
    write_doc(g if not q2_only else
              pd.read_csv(TABLES / "c24_depmap_gene_level.tsv", sep="\t"),
              s)
    log("wrote the Control 24 tables")
    return 0


def write_doc(g: pd.DataFrame, s: pd.DataFrame) -> None:
    lines = [
        "# Control 24 — does the network add anything on the functional "
        "task, once hubness is taken out?",
        "",
        "**Date:** 2026-11-18  ",
        "**Script:** `scripts/run_c24_depmap_within_degree.py`  ",
        "**Tables:** `c24_depmap_gene_level.tsv`, "
        "`c24_depmap_within_degree.tsv`",
        "",
        "## The question",
        "",
        "The released DepMap validation says the prior does not beat the "
        "degree baseline on functional essentiality (at K = 30: prior "
        "5.8 %, degree 69.2 %, random 3.1 %). Two things were open. "
        "**Q1**: at the gene level, does the prior discriminate the "
        "essential genes better than degree does, and does it still do so "
        "after the degree effect is regressed out of its score? **Q2**: "
        "inside a degree stratum, does the prior's pick find essential "
        "genes more often than a random pick from the same stratum?",
        "",
        "## Q1: the gene level (the median AUROC for functional "
        "essentiality, over the cell lines)",
        "",
        "| threshold | the prior | degree | the residual (the prior minus "
        "its degree part) |",
        "|---|---|---|---|",
    ]
    for th in sorted(g.threshold.unique()):
        sub = g[g.threshold == th]
        vals = []
        for tag in ("prior", "degree", "residual"):
            v = sub[sub.scorer == tag].auc_median
            vals.append(f"{v.median():.4f}" if len(v) else "-")
        lines.append(f"| {th} | " + " | ".join(vals) + " |")

    lines += ["", "### Q1 per arm (the threshold -0.5 rows, the median over "
              "the folds)", "",
              "| network | label | the prior | degree | the residual | the "
              "size of the pool |", "|---|---|---|---|---|---|"]
    sub = g[g.threshold == -0.5]
    for (net, label), grp in sub.groupby(["network", "label"]):
        vals = []
        for tag in ("prior", "degree", "residual"):
            v = grp[grp.scorer == tag].auc_median
            vals.append(f"{v.median():.4f}" if len(v) else "-")
        npool = grp["n_pool"].median() if "n_pool" in grp.columns else np.nan
        lines.append(f"| {net} | {label} | " + " | ".join(vals)
                     + f" | {npool:.0f} |")

    lines += [
        "",
        "## Q2: inside a degree stratum (the share of cell lines enriched, "
        "K = 30)",
        "",
        "| threshold | stratum (by degree, 1 = lowest) | the prior's picks "
        "| degree's picks | the residual's picks | a random pick from the "
        "same stratum |",
        "|---|---|---|---|---|---|",
    ]
    for th in sorted(s.threshold.unique()):
        for st in sorted(s[s.threshold == th].stratum.unique()):
            sub2 = s[(s.threshold == th) & (s.stratum == st)
                     & (s.K == 30)]
            if not len(sub2):
                continue
            cells = []

            def med(tag, pick):
                v = sub2[(sub2.scorer == tag) & (sub2.pick == pick)]
                return f"{v.frac_enriched.median():.1%}" if len(v) else "-"
            cells.append(med("prior", "top"))
            cells.append(med("degree", "top"))
            cells.append(med("residual", "top"))
            cells.append(med("prior", "random"))
            lines.append(f"| {th} | {st} | " + " | ".join(cells) + " |")

    lines += [
        "",
        "## Caveats",
        "",
        "1. The pools are the held-out negatives of each fold, the same "
        "pool the released DepMap analysis used, so the numbers are "
        "comparable with it.",
        "2. The residual is the ordinary-least-squares residual of the "
        "prior's score on log1p(degree) inside the pool.",
        "3. Q2 uses one fold in five (the per-cell-line medians are stable "
        "across folds; Q1 uses all of them).",
        "4. The pools differ in size between the labels (the label with "
        "many positives leaves few negatives), so Q2 can only be run where "
        "a stratum holds enough genes for a top-K pick.",
        "5. **Exploratory** (amendment 4): nothing in the pre-registered "
        "endpoint is touched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
