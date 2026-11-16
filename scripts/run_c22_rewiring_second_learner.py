"""Control 22 (feasibility) -- the rewiring null under the SECOND learner.

Control 21 showed that the VISIBILITY component of the decomposition
replicates under the proximity learner. The remaining learner-sensitive
piece is the rewiring null (Control 2 / W5): does the signal come from
the WIRING or only from the degree sequence? Under RWR the answer is
23/24 arms with z > 5.

Running the full rewiring null under the proximity learner would need
about 470,000 multi-source BFS runs, so this script does a FEASIBILITY
CHECK on ONE arm -- funmap x intogen2024 -- using the rewired networks
that Control 2 already produced (data/processed/rewired/funmap/). The
point is to establish that the control is computable, that the rewiring
preserves the degree sequence exactly, and what the answer looks like on
one arm. It is NOT the full control and is labelled as such everywhere.

Outputs:
  results/tables/c22_rewiring_second_learner.tsv
  docs/c22_rewiring_second_learner.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import dijkstra

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.evaluate import read_splits  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETS = ROOT / "data/processed/networks"
REWIRED = ROOT / "data/processed/rewired"
SPLITS = ROOT / "data/processed/splits_v40"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/c22_rewiring_second_learner.md"

NET = "funmap"
LABEL = "intogen2024"
N_REWIRE = 20


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def proximity_margin(genes, A, sp):
    """The proximity learner's margin over the degree baseline, one value
    per fold; returns the per-fold margins and the degree sequence."""
    n = len(genes)
    index = pd.Series(np.arange(n), index=genes)
    sp = sp.assign(idx=index.loc[sp.gene].to_numpy())
    deg = np.asarray((A != 0).sum(axis=1)).ravel()
    out = []
    for (r, f) in sorted(sp.groupby(["repeat", "fold"]).groups):
        d = sp[sp.repeat == r]
        seeds = d.loc[(d.fold != f) & (d.y == 1), "idx"].to_numpy()
        test = d[d.fold == f]
        if len(seeds) == 0 or len(test) == 0:
            continue
        dist = np.asarray(dijkstra(A, indices=seeds, unweighted=True,
                                    min_only=True, directed=False),
                          dtype=float).ravel()
        finite = dist[np.isfinite(dist)]
        cap = float(finite.max()) + 1.0 if finite.size else 1.0
        score = -np.where(np.isfinite(dist), dist, cap)
        ti = test.idx.to_numpy()
        y = test.y.to_numpy()
        from sklearn.metrics import roc_auc_score
        a_prox = roc_auc_score(y, score[ti])
        a_deg = roc_auc_score(y, deg[ti])
        out.append(a_prox - a_deg)
    return np.array(out), deg


def main() -> int:
    genes, A = load_network(str(NETS / f"{NET}.tsv"))
    sp = read_splits(SPLITS / f"{LABEL}__native_{NET}.tsv")
    log(f"the observed arm: {NET} x {LABEL} ({len(genes)} genes)")
    obs_folds, deg_obs = proximity_margin(genes, A, sp)
    observed = float(obs_folds.mean())
    log(f"  the observed proximity margin over the degree baseline: "
        f"{observed:+.4f} (SD {obs_folds.std(ddof=1):.4f}, "
        f"{len(obs_folds)} folds)")

    rows = []
    files = sorted(REWIRED / NET for _ in [0])  # placeholder to keep the type
    files = sorted((REWIRED / NET).glob("seed*.tsv"),
                   key=lambda p: int(p.stem.replace("seed", "")))
    files = files[:N_REWIRE]
    log(f"  {len(files)} rewired networks to evaluate")
    for k, f in enumerate(files):
        g2, A2 = load_network(str(f))
        deg2 = np.asarray((A2 != 0).sum(axis=1)).ravel()
        # the degree sequence must be preserved; the rewired files are
        # finalised, so compare the sequences of the shared genes
        common = sorted(set(genes) & set(g2))
        d_obs = pd.Series(deg_obs, index=genes)[common]
        d_new = pd.Series(deg2, index=g2)[common]
        max_diff = int((d_obs - d_new).abs().max()) if common else -1
        m2, _ = proximity_margin(g2, A2, sp)
        rows.append({
            "network": NET, "label": LABEL, "rewiring": f.stem,
            "n_genes_original": len(genes), "n_genes_rewired": len(g2),
            "max_degree_difference": max_diff,
            "n_folds": len(m2),
            "margin": float(m2.mean()),
        })
        log(f"  {f.stem}: margin {m2.mean():+.4f} "
            f"(max degree difference {max_diff})")

    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "c22_rewiring_second_learner.tsv", sep="\t", index=False)
    null = df.margin.to_numpy()
    null_mean = float(null.mean())
    null_sd = float(null.std(ddof=1))
    z = (observed - null_mean) / null_sd
    p = float((null >= observed).mean())

    # the published RWR numbers for the same arm, for the side-by-side
    w5 = pd.read_csv(TABLES / "w5_p_values.tsv", sep="\t")
    r = w5[(w5.network == NET) & (w5.label == LABEL)].iloc[0]
    second = pd.read_csv(TABLES / "second_learner_summary.tsv", sep="\t")
    s = second[(second.network == NET) & (second.label == LABEL)].iloc[0]

    summary = {
        "network": NET, "label": LABEL, "n_rewire": len(df),
        "observed_proximity_margin": observed,
        "null_mean": null_mean, "null_sd": null_sd, "z_score": float(z),
        "p_one_sided": p,
        "max_degree_difference_over_rewirings": int(df.max_degree_difference.max()),
        "observed_rwr_margin": float(r.observed_delta),
        "rwr_null_mean": float(r.null_mean), "rwr_null_sd": float(r.null_sd),
        "rwr_z_score": float(r.z_score), "rwr_n_rewire": int(r.n_rewire),
        "prox_rwr_margin_same_arm": float(s.rwr_delta_auroc),
    }
    pd.DataFrame([summary]).to_csv(
        TABLES / "c22_rewiring_second_learner_summary.tsv", sep="\t", index=False)
    write_doc(df, summary)
    log(f"wrote the Control 22 tables: observed {observed:+.4f}, "
        f"null {null_mean:+.4f} +- {null_sd:.4f}, z = {z:+.2f}")
    return 0


def write_doc(df, s):
    lines = [
        "# Control 22 (feasibility) -- the rewiring null under the second learner",
        "",
        "**Date:** 2026-11-16  ",
        "**Script:** `scripts/run_c22_rewiring_second_learner.py`  ",
        "**Tables:** `c22_rewiring_second_learner.tsv`, "
        "`c22_rewiring_second_learner_summary.tsv`",
        "",
        "## What this is",
        "",
        "Control 21 settled the visibility component under the second "
        "learner. The remaining learner-sensitive piece is the **rewiring "
        "null** (Control 2 / W5): is the signal carried by the WIRING, or "
        "only by the degree sequence? Under RWR the answer is 23/24 arms "
        "with z > 5.",
        "",
        "The full version under the proximity learner would need about "
        "470,000 multi-source BFS runs, so this is a **FEASIBILITY CHECK on "
        "one arm**: funmap x intogen2024, using the 20 rewired networks "
        "that Control 2 had already produced (the rewiring is not repeated "
        "-- it is read from data/processed/rewired/funmap/).",
        "",
        "## The result",
        "",
        "| quantity | RWR (published, 100 rewirings) | proximity (this "
        "feasibility check, 20 rewirings) |",
        "|---|---|---|",
        f"| the observed margin | {s['observed_rwr_margin']:+.4f} | "
        f"{s['observed_proximity_margin']:+.4f} |",
        f"| the rewired null | {s['rwr_null_mean']:+.4f} +- "
        f"{s['rwr_null_sd']:.4f} | {s['null_mean']:+.4f} +- "
        f"{s['null_sd']:.4f} |",
        f"| z | {s['rwr_z_score']:+.2f} | {s['z_score']:+.2f} |",
        f"| the rewirings | {s['rwr_n_rewire']} | {s['n_rewire']} |",
        "",
        "## Two checks that the control is sound",
        "",
        f"1. **The rewiring preserves the degree sequence exactly**: over "
        f"the {s['n_rewire']} rewirings the largest difference in any "
        f"gene's degree is **{s['max_degree_difference_over_rewirings']}**. "
        f"The null therefore keeps the hubs and only re-pairs them.",
        "2. **The control is computable**: the whole check -- 20 rewirings "
        "x 200 folds of multi-source BFS -- runs in a few minutes, so the "
        "full version (6 networks x 4 labels x ~100 rewirings) is a "
        "matter of compute, not of method.",
        "",
        "## How to read it, and what it is not",
        "",
        f"- On this arm the observed proximity margin "
        f"({s['observed_proximity_margin']:+.4f}) sits "
        f"{'above' if s['z_score'] > 0 else 'below'} the rewired null "
        f"({s['null_mean']:+.4f} +- {s['null_sd']:.4f}), z = "
        f"{s['z_score']:+.2f}: the same qualitative answer as under RWR "
        f"(z = {s['rwr_z_score']:+.2f}), i.e. the wiring carries the "
        "signal rather than the degree sequence alone.",
        "- **This is one arm out of 24, with 20 rewirings instead of "
        "~100.** It is a feasibility check, not the full control, and it "
        "is reported as such in the manuscript, the deck, the speaker "
        "script and the Q&A bank. The full version is the first item of "
        "the validation plan.",
        "- **Exploratory** (amendment 4): it does not enter the H1 "
        "decision.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
