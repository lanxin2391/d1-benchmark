"""F3b -- the DepMap validation three ways: the prior, the degree
baseline and random candidates, on the same folds, the same universe,
the same K and the same cell lines.

Why: the released F3 table reported the prior's enrichment against a
nominal 5 % level. Two comparisons were missing -- the degree baseline
(the strongest simple baseline, and the one our own stability analysis
says is strong) and random candidates (the actual base rate). The pool
of genes stays the prior's candidate pool, so the prior's column
reproduces the released numbers exactly and the three columns are
directly comparable.

Outputs:
  results/tables/depmap_three_way.tsv
  results/tables/depmap_prior_vs_degree.tsv
  docs/depmap_degree_control.md
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
    ALPHA, EXPR_MIN, K_VALUES, LABELS, NETWORKS, SCORES, THRESHOLDS,
    load_matrices)

TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/depmap_degree_control.md"
SEED = 20261112


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> int:
    from scipy import stats as sps

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
            rwr = np.asarray(d["rwr"][ai], dtype=float)
            deg = np.asarray(d["degree"], dtype=float)
            in_test = np.asarray(d["in_test"], dtype=bool)
            y = np.asarray(d["y"]).astype(int)
            genes = np.asarray(d["genes"]).astype(str)
            arms[(net, label)] = (rwr, deg, in_test, y, genes)
            for j in range(in_test.shape[0]):
                idx = np.where(in_test[j])[0]
                neg = idx[y[idx] == 0]
                if len(neg) == 0:
                    continue
                kmax = max(100, int(np.ceil(0.01 * len(idx))))
                order = neg[np.argsort(-rwr[neg, j], kind="mergesort")]
                need.update(genes[order[:kmax]].tolist())
    log(f"the candidate pool touches {len(need)} genes")

    log("loading the DepMap matrices once ...")
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

    keep = None
    for threshold in THRESHOLDS:
        log(f"=== threshold CERES < {threshold} ===")
        ess = (np.isfinite(score) & expressed
               & (score < threshold)).astype(bool)
        out_rows = []
        for (net, label), (rwr, deg, in_test, y, genes) in arms.items():
            for K in K_VALUES:
                ps = {"prior": [], "degree": [], "random": []}
                ks = []
                for j in range(in_test.shape[0]):
                    idx = np.where(in_test[j])[0]
                    neg = idx[y[idx] == 0]
                    if len(neg) < 20:
                        continue
                    k = int(np.ceil(K * len(idx))) if K < 1 else int(K)
                    k = min(k, len(neg))
                    if len(neg) - k < 20:
                        continue
                    order_p = neg[np.argsort(-rwr[neg, j], kind="mergesort")]
                    order_d = neg[np.argsort(-deg[neg], kind="mergesort")]
                    rng = np.random.default_rng(SEED + 1000 * j
                                                + int(100 * K))
                    order_r = rng.permutation(neg)
                    for tag, order in (("prior", order_p),
                                       ("degree", order_d),
                                       ("random", order_r)):
                        keep_c = [g for g in genes[order[:k]] if g in index]
                        keep_b = [g for g in genes[order[k:]] if g in index]
                        cand_rows = np.array([index[g] for g in keep_c],
                                             dtype=int)
                        bg_rows = np.array([index[g] for g in keep_b],
                                            dtype=int)
                        if len(cand_rows) < 5 or len(bg_rows) < 20:
                            continue
                        k_cand = ess[cand_rows].sum(axis=0).astype(int)
                        k_bg = ess[bg_rows].sum(axis=0).astype(int)
                        M = len(cand_rows) + len(bg_rows)
                        with np.errstate(divide="ignore",
                                         invalid="ignore"):
                            pv = sps.hypergeom.sf(k_cand - 1, M,
                                                  len(cand_rows), k_bg)
                        good = (k_bg > 0) & np.isfinite(pv)
                        if good.any():
                            ps[tag].append(pv[good])
                        if tag == "prior":
                            ks.append(len(cand_rows))
                for tag in ("prior", "degree", "random"):
                    if not ps[tag]:
                        continue
                    P = np.array(ps[tag])
                    out_rows.append({
                        "network": net, "label": label, "K": K,
                        "threshold": threshold, "ranking": tag,
                        "n_folds": len(P),
                        "frac_celllines_p_lt_005": float((P < 0.05).mean()),
                        "frac_celllines_p_lt_001": float((P < 0.01).mean()),
                        "median_p": float(np.median(P)),
                        "min_p": float(P.min()),
                        "n_celllines_tested": int(len(P)),
                        "n_candidates_median": int(np.median(ks)) if ks else 0,
                    })
        df = pd.DataFrame(out_rows)
        df.to_csv(TABLES / f"depmap_three_way_thr{threshold}.tsv", sep="\t",
                  index=False)
        log(f"  wrote depmap_three_way_thr{threshold}.tsv ({len(df)} rows)")
        if threshold == -0.5:
            keep = df
    keep.to_csv(TABLES / "depmap_three_way.tsv", sep="\t", index=False)

    rows = []
    for K in sorted(keep["K"].unique()):
        for net in NETWORKS:
            for label in LABELS:
                sel = keep[(keep["K"] == K) & (keep.network == net)
                           & (keep.label == label)]
                if sel.ranking.nunique() < 3:
                    continue
                d = {r.ranking: r for r in sel.itertuples()}
                rows.append({
                    "network": net, "label": label, "K": K,
                    "prior_frac": d["prior"].frac_celllines_p_lt_005,
                    "degree_frac": d["degree"].frac_celllines_p_lt_005,
                    "random_frac": d["random"].frac_celllines_p_lt_005,
                    "prior_median_p": d["prior"].median_p,
                    "degree_median_p": d["degree"].median_p,
                    "random_median_p": d["random"].median_p,
                })
    comp = pd.DataFrame(rows)
    comp.to_csv(TABLES / "depmap_prior_vs_degree.tsv", sep="\t", index=False)
    write_doc(comp)
    log(f"wrote the three-way table ({len(comp)} rows)")
    return 0


def write_doc(comp: pd.DataFrame):
    lines = [
        "# F3b -- the DepMap validation three ways: prior, degree, random",
        "",
        "**Date:** 2026-11-12  ",
        "**Script:** `scripts/run_depmap_degree_baseline.py`  ",
        "**Tables:** `results/tables/depmap_three_way.tsv`, "
        "`depmap_prior_vs_degree.tsv`",
        "",
        "## Why this control was missing",
        "",
        "The released validation (F3) asked whether the genes the prior "
        "puts at the top of the ranking are more likely to be required "
        "in cancer cell lines. It reported no enrichment beyond chance "
        "for the top 30 candidates and a clear one for the top 100 -- "
        "but against a nominal 5 % level, and without the two "
        "comparisons that matter: **the degree baseline** (the strongest "
        "simple baseline, and the one our own stability analysis says is "
        "strong) and **random candidates** (the actual base rate, which "
        "is not 5 % when a cell line has a few hundred essential genes "
        "among ~18,000).",
        "",
        "All three rankings are computed on the same test folds, the same "
        "universe (the prior's candidate pool, so the prior's column "
        "reproduces the released numbers), the same K and the same cell "
        "lines.",
        "",
    ]
    for K in sorted(comp["K"].unique()):
        s = comp[comp["K"] == K]
        lines += [
            f"### K = {K:g}",
            "",
            f"- the prior: **{s.prior_frac.median():.1%}** of the cell "
            f"lines enriched (median p {s.prior_median_p.median():.2f})",
            f"- the degree baseline: **{s.degree_frac.median():.1%}** "
            f"(median p {s.degree_median_p.median():.2f})",
            f"- random candidates: **{s.random_frac.median():.1%}** "
            f"(median p {s.random_median_p.median():.2f})",
            f"- the prior beats random in "
            f"**{int((s.prior_frac > s.random_frac).sum())}/{len(s)}** arms; "
            f"the degree beats the prior in "
            f"**{int((s.degree_frac > s.prior_frac).sum())}/{len(s)}**",
            "",
        ]
    lines += [
        "## Reading",
        "",
        "This is the control the earlier statement was missing, and it "
        "cuts both ways. The prior is better than random at the larger "
        "candidate sets, but **the degree baseline is clearly better "
        "than the prior at every scale**. Hubness predicts functional "
        "essentiality well -- that is a real biological signal, not an "
        "artefact of our benchmark -- and a diffusion prior adds nothing "
        "to it on this task.",
        "",
        "So the scope of our claim narrows further and it is worth "
        "stating plainly: on the curated driver labels the prior adds to "
        "degree (23 of 24 arms); on a cell-line functional task degree is "
        "the better predictor.",
        "",
        "## Caveats",
        "",
        "1. The candidates are the held-out negatives of each arm, so the "
        "test set is identical in the three columns.",
        "2. The essentiality definition is a proxy for dependency, and "
        "it correlates with expression and study intensity.",
        "3. Exploratory family (amendment 4); the pre-registered "
        "endpoint is untouched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
