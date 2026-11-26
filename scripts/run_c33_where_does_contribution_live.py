"""Control 33 -- where does the network's contribution live?

This is the repositioned item 7 of the post-defence plan. The original
item (the second paper's H6) was a **conditional repair**: "if the RWR
advantage disappears once the hubs are removed, we need a
degree-orthogonalised model". The trigger was tested on 2026-11-05 and it
did not fire -- the degree-penalised walk is *larger* (+0.0979), not
smaller -- so as a repair the idea was retired. The plan listed it again
without that context, and a first note (after Control 24) said the
residual had "no signal", which over-generalised: Control 24 measured
the residual on the **functional** task, not on the driver task this
result is about.

What was never done, and is the useful version of the item, is the
**decomposition of the network's contribution on the primary task**:

  1. the per-gene residual of the prior's score after regressing out the
     degree, after regressing out the publication count, and after
     regressing out both -- is the network's own component still above
     chance for the driver label?
  2. the **self-written model** (the other half of H6): a logistic model
     with the degree and the publication count as features, with and
     without the network score -- how much does the network add to a
     degree-and-visibility model?

Everything is computed on the same folds, the same test genes and the
same endpoint as the primary analysis, from the released score dumps.

Outputs:
  results/tables/c33_residual_decomposition.tsv
  results/tables/c33_self_written_model.tsv
  results/tables/c33_summary.tsv
  docs/c33_where_does_the_contribution_live.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SCORES = ROOT / "results/scores_v40"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/c33_where_does_the_contribution_live.md"
NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
ALPHA = 0.5


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def tie_ranks(x: np.ndarray) -> np.ndarray:
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


def auc(scores: np.ndarray, labels: np.ndarray) -> float:
    pos = scores[labels > 0]
    neg = scores[labels <= 0]
    if pos.size < 2 or neg.size < 2:
        return np.nan
    r = tie_ranks(scores)
    n_pos = float((labels > 0).sum())
    n_neg = float((labels <= 0).sum())
    return float((r[labels > 0].sum() - n_pos * (n_pos + 1) / 2)
                 / (n_pos * n_neg))


def residualise(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """The part of y that the columns of X do not explain."""
    A = np.column_stack([np.ones(len(y)), X])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ coef


def logistic(X: np.ndarray, y: np.ndarray, ridge: float = 1e-6,
              iters: int = 60) -> np.ndarray:
    """A small ridge-stabilised IRLS logistic fit."""
    A = np.column_stack([np.ones(len(y)), X])
    beta = np.zeros(A.shape[1])
    for _ in range(iters):
        eta = A @ beta
        p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        W = np.clip(p * (1 - p), 1e-9, None)
        z = eta + (y - p) / W
        WA = A * W[:, None]
        try:
            new = np.linalg.solve(WA.T @ A + ridge * np.eye(A.shape[1]),
                                   WA.T @ z)
        except np.linalg.LinAlgError:
            break
        if np.max(np.abs(new - beta)) < 1e-8:
            beta = new
            break
        beta = new
    return beta


def main() -> int:
    pc = pd.read_csv(ROOT / "data/processed/labels/pubcount.tsv", sep="\t")
    pub = {g: float(np.log10(n + 1.0)) for g, n in
           zip(pc["gene"], pc["n_pubmed"])}
    res_rows = []
    mod_rows = []
    for net in NETWORKS:
        for label in LABELS:
            fp = SCORES / f"{net}__{label}.npz"
            if not fp.exists():
                continue
            d = np.load(fp, allow_pickle=True)
            ai = int(np.argmin(np.abs(np.asarray(d["alphas"], dtype=float)
                                      - ALPHA)))
            genes = np.asarray(d["genes"]).astype(str)
            y = np.asarray(d["y"]).astype(int)
            degree = np.asarray(d["degree"], dtype=float)
            rwr = np.asarray(d["rwr"][ai], dtype=float)
            in_test = np.asarray(d["in_test"], dtype=bool)
            vis = np.array([pub.get(g, 0.0) for g in genes])
            for j in range(in_test.shape[0]):
                idx = np.where(in_test[j])[0]
                lab = y[idx]
                if (lab > 0).sum() < 5 or (lab <= 0).sum() < 20:
                    continue
                s = rwr[idx, j]
                x1 = np.log1p(degree[idx])
                x2 = vis[idx]
                # standardise for the model fits
                def z(v: np.ndarray) -> np.ndarray:
                    sd = v.std()
                    return (v - v.mean()) / (sd if sd > 0 else 1.0)
                z1, z2, zs = z(x1), z(x2), z(s)
                variants = {
                    "the prior": s,
                    "the degree": x1,
                    "the publication count": x2,
                    "the prior, degree removed":
                        residualise(s, x1[:, None]),
                    "the prior, the publication count removed":
                        residualise(s, x2[:, None]),
                    "the prior, the degree and the publication count removed":
                        residualise(s, np.column_stack([x1, x2])),
                    "the degree, the publication count removed":
                        residualise(x1, x2[:, None]),
                }
                for name, sc in variants.items():
                    a = auc(sc, lab)
                    if np.isfinite(a):
                        res_rows.append({
                            "network": net, "label": label, "fold": j,
                            "scorer": name, "auroc": a,
                            "auroc_minus_half": a - 0.5})
                # ---- the self-written model ----
                b_base = logistic(np.column_stack([z1, z2]), lab.astype(float))
                b_full = logistic(np.column_stack([z1, z2, zs]),
                                  lab.astype(float))
                a_base = auc(b_base[0] + b_base[1] * z1 + b_base[2] * z2, lab)
                a_full = auc(b_full[0] + b_full[1] * z1 + b_full[2] * z2
                             + b_full[3] * zs, lab)
                mod_rows.append({
                    "network": net, "label": label, "fold": j,
                    "beta_network": float(b_full[3]),
                    "beta_degree": float(b_full[1]),
                    "beta_visibility": float(b_full[2]),
                    "auc_degree_visibility": a_base,
                    "auc_with_network": a_full,
                    "delta_auc": a_full - a_base,
                })
            log(f"  {net} x {label}: done")
    r = pd.DataFrame(res_rows)
    r.to_csv(TABLES / "c33_residual_decomposition.tsv", sep="\t", index=False)
    m = pd.DataFrame(mod_rows)
    m.to_csv(TABLES / "c33_self_written_model.tsv", sep="\t", index=False)
    write_doc(r, m)
    log(f"wrote the Control 33 tables ({len(r)} fold rows, {len(m)} model "
        "rows)")
    return 0


def write_doc(r: pd.DataFrame, m: pd.DataFrame) -> None:
    arm = r.groupby(["network", "label", "scorer"]).auroc.median().unstack()
    mod = m.groupby(["network", "label"]).agg(
        beta_network=("beta_network", "median"),
        delta_auc=("delta_auc", "mean")).reset_index()
    lines = [
        "# Control 33 — where does the network's contribution live?",
        "",
        "**Date:** 2026-11-26  ",
        "**Script:** `scripts/run_c33_where_does_contribution_live.py`  ",
        "**Tables:** `c33_residual_decomposition.tsv`, "
        "`c33_self_written_model.tsv`, `c33_summary.tsv`",
        "",
        "## What this item is, after the repositioning",
        "",
        "The original item (the second paper's H6) was a **conditional "
        "repair**: write a degree-orthogonalised model *if* the RWR "
        "advantage disappeared once the hubs were removed. That trigger was "
        "tested on 2026-11-05 and it did not fire — the degree-penalised "
        "walk is larger (+0.0979), not smaller — so as a repair the idea was "
        "retired. The useful version of the item is the **decomposition of "
        "the network's contribution on the primary task**:",
        "",
        "1. the per-gene residual of the prior's score after regressing out "
        "the degree, after regressing out the publication count, and after "
        "regressing out both — is the network's own component still above "
        "chance for the driver label?",
        "2. the **self-written model** (the other half of H6): a logistic "
        "model with the degree and the publication count as features, with "
        "and without the network score.",
        "",
        "Both are computed on the same folds, the same test genes and the "
        "same endpoint as the primary analysis.",
        "",
        "## 1. The AUROC of each variant (the median over the folds, per "
        "arm)",
        "",
        "| network | label | the prior | the degree | the publication count "
        "| the prior, degree removed | the prior, visibility removed | the "
        "prior, both removed |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for _, row in mod.iterrows():
        net, lab = row["network"], row["label"]
        if net not in arm.index.get_level_values(0):
            continue
        a = arm.loc[net]
        cells = []
        for c in ("the prior", "the degree", "the publication count",
                  "the prior, degree removed",
                  "the prior, the publication count removed",
                  "the prior, the degree and the publication count removed"):
            v = a.loc[lab, c] if (lab in a.index and c in a.columns) else np.nan
            cells.append(f"{v:.4f}" if pd.notna(v) else "-")
        lines.append(f"| {net} | {lab} | " + " | ".join(cells) + " |")
    lines += [
        "",
        "## 2. The self-written model",
        "",
        "| network | label | the coefficient on the network score | the AUROC "
        "gain of adding the network to the degree+visibility model |",
        "|---|---|---|---|",
    ]
    for _, row in mod.iterrows():
        lines.append(f"| {row['network']} | {row['label']} | "
                     f"{row['beta_network']:+.4f} | "
                     f"{row['delta_auc']:+.4f} |")
    lines += [
        "",
        "## Reading",
        "",
        "## Caveats",
        "",
        "1. The regressions are fitted on the test genes of each fold (the "
        "same population the AUROC is computed on), with no label "
        "information, so the residual is 'what the degree and the "
        "publication count do not explain *here*'.",
        "2. The publication count is a per-gene covariate; the "
        "publication-matched null of F.12 is a different (and stronger) "
        "control, because it matches the publication count at the level of "
        "the bait rather than regressing it out of the score.",
        "3. The logistic model is a plain ridge-stabilised IRLS fit with "
        "standardised features; the coefficient on the network score is "
        "therefore comparable across arms.",
        "4. **Exploratory** (amendment 4): no pre-registered endpoint is "
        "touched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
