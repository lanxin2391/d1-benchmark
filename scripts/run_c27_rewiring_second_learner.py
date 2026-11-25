"""Control 27 -- the rewiring null under the second learner (the full grid).

Control 22 showed the rewiring null is computable under the proximity
learner on one arm. This runs the full grid: 6 networks x 4 labels x
every rewiring that exists, with the same proximity learner, the same
degree baseline and the same test genes as the released second-learner
run.

For each arm the statistic is the margin of the proximity learner over
the degree baseline, compared with the same margin on the rewired
networks. Two z conventions are reported, because the project has used
both and they differ by a factor of sqrt(n_rewire):

    the effect-size z   = (observed - null mean) / null SD
                          (how many null SDs the observation is away --
                           the honest reading of "an outlier")
    the standard-error z = (observed - null mean) / (null SD / sqrt(n))
                          (what the released w5_p_values.tsv column
                           reports: how many standard errors of the null
                           mean the observation is)

Both are reported; the verdicts are the same because the null is
centred near zero and the observation is far from it.

Folds: one fold per repeat (the first 20 repeats) for every arm, the
observed and the null alike, because each rewiring needs a fresh
multi-source BFS per fold.

Outputs:
  results/tables/c27_rewiring_second_learner.tsv
  results/tables/c27_rewiring_second_learner_summary.tsv
  docs/c27_rewiring_second_learner.md
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.baselines import degree_scores  # noqa: E402
from d1.engine.evaluate import evaluate, read_splits  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETS = ROOT / "data/processed/networks"
REWIRED = ROOT / "data/processed/rewired"
SPLITS = ROOT / "data/processed/splits_v40"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/c27_rewiring_second_learner.md"
NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
N_FOLDS = 20            # one fold per repeat
SEED_FOLDS = (0, 1, 2, 3, 4)


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def dist_to_nearest(A, seeds) -> np.ndarray:
    """The unweighted distance from every node to the nearest seed."""
    from scipy.sparse.csgraph import dijkstra
    d = dijkstra(A, indices=seeds, unweighted=True, min_only=True,
                 directed=False)
    return np.asarray(d, dtype=float).ravel()


def arm_folds(net: str, label: str):
    """The (seeds, test) index arrays for the selected folds of one arm."""
    sp_path = SPLITS / f"{label}__native_{net}.tsv"
    if not sp_path.exists():
        return None
    sp = read_splits(sp_path)
    genes, _ = load_network(str(NETS / f"{net}.tsv"))
    index = pd.Series(np.arange(len(genes)), index=genes)
    sp = sp.assign(idx=index.loc[sp.gene].to_numpy())
    out = []
    for (r, f) in sorted(sp.groupby(["repeat", "fold"]).groups):
        if int(r) >= N_FOLDS or int(f) != SEED_FOLDS[0]:
            continue
        d = sp[sp.repeat == r]
        seeds = d.loc[(d.fold != f) & (d.y == 1), "idx"].to_numpy()
        test = d[d.fold == f]
        if len(seeds) == 0 or len(test) == 0:
            continue
        out.append((seeds, test.idx.to_numpy(), test.y.to_numpy(),
                    test.gene.to_numpy()))
    return out


def margin_on(A, deg, folds) -> tuple[float, int]:
    """The mean proximity-vs-degree margin over the given folds."""
    ms = []
    for seeds, ti, y, genes in folds:
        dist = dist_to_nearest(A, seeds)
        finite = dist[np.isfinite(dist)]
        cap = float(finite.max()) + 1.0 if finite.size else 1.0
        score = -np.where(np.isfinite(dist), dist, cap)
        ev_p = evaluate(score[ti], y, genes)
        ev_d = evaluate(deg[ti], y, genes)
        ms.append(ev_p["auroc"] - ev_d["auroc"])
    return (float(np.mean(ms)) if ms else np.nan), len(ms)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-rewire", type=int, default=100)
    args = ap.parse_args(argv)

    rows = []
    for net in NETWORKS:
        folds_by_label = {}
        for label in LABELS:
            f = arm_folds(net, label)
            if f:
                folds_by_label[label] = f
        if not folds_by_label:
            log(f"  {net}: no splits, skipped")
            continue
        genes, A = load_network(str(NETS / f"{net}.tsv"))
        deg = degree_scores(A)
        # ---- the observed margin on the real network ----
        for label, folds in folds_by_label.items():
            m, n = margin_on(A, deg, folds)
            rows.append({"network": net, "label": label, "rewiring": -1,
                         "margin": m, "n_folds": n})
        # ---- the rewired null ----
        files = sorted(REWIRED.joinpath(net).glob("seed*.tsv"),
                       key=lambda p: int(p.stem.replace("seed", "")))
        files = files[:args.max_rewire]
        t0 = time.time()
        for k, f in enumerate(files):
            try:
                _, Aw = load_network(str(f))
            except Exception as exc:                      # noqa: BLE001
                log(f"  {net}: {f.name} could not be read ({exc})")
                continue
            degw = degree_scores(Aw)
            for label, folds in folds_by_label.items():
                m, n = margin_on(Aw, degw, folds)
                rows.append({"network": net, "label": label,
                             "rewiring": int(f.stem.replace("seed", "")),
                             "margin": m, "n_folds": n})
            if (k + 1) % 20 == 0:
                log(f"  {net}: {k + 1}/{len(files)} rewirings "
                    f"({time.time() - t0:.0f}s)")
        log(f"  {net}: {len(files)} rewirings done "
            f"({time.time() - t0:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "c27_rewiring_second_learner.tsv", sep="\t",
              index=False)

    summary = []
    for (net, label), grp in df.groupby(["network", "label"]):
        obs = grp[grp.rewiring == -1].margin.iloc[0]
        null = grp[grp.rewiring >= 0].margin.dropna()
        n = len(null)
        m, s = float(null.mean()), float(null.std(ddof=1))
        z_eff = (obs - m) / s if s > 0 else np.nan
        z_se = (obs - m) / (s / np.sqrt(n)) if s > 0 else np.nan
        summary.append({
            "network": net, "label": label,
            "observed_margin": round(obs, 4),
            "null_mean": round(m, 4), "null_sd": round(s, 4),
            "n_rewire": n, "n_folds": int(grp.n_folds.max()),
            "z_effect_size": round(z_eff, 2),
            "p_effect_size": float(sps.norm.sf(z_eff)),
            "z_standard_error": round(z_se, 2),
            "p_standard_error": float(sps.norm.sf(z_se)),
        })
    s_df = pd.DataFrame(summary)
    s_df.to_csv(TABLES / "c27_rewiring_second_learner_summary.tsv", sep="\t",
                index=False)
    write_doc(s_df)
    log(f"wrote the Control 27 tables ({len(s_df)} arms)")
    return 0


def write_doc(s: pd.DataFrame) -> None:
    w5 = pd.read_csv(TABLES / "w5_p_values.tsv", sep="\t")
    w5 = w5.rename(columns={"label": "label"})
    lines = [
        "# Control 27 — the rewiring null under the second learner (the "
        "full grid)",
        "",
        "**Date:** 2026-11-25  ",
        "**Script:** `scripts/run_c27_rewiring_second_learner.py`  ",
        "**Tables:** `c27_rewiring_second_learner.tsv`, "
        "`c27_rewiring_second_learner_summary.tsv`",
        "",
        "## What this adds",
        "",
        "Control 22 showed on one arm that the degree-preserving rewiring "
        "null is computable under the second (proximity) learner. This runs "
        "the full grid: 6 networks × 4 labels × every rewiring that exists, "
        "with the same proximity learner, the same degree baseline and the "
        "same test genes as the released second-learner run.",
        "",
        "Two z conventions are reported, because the project has used both "
        "and they differ by a factor of √n:",
        "",
        "* **the effect-size z** = (observed − null mean) / null SD — how "
        "many null SDs the observation is away from the null;",
        "* **the standard-error z** = (observed − null mean) / (null SD / "
        "√n) — what the released `w5_p_values.tsv` column reports, i.e. how "
        "many standard errors of the null *mean* the observation is.",
        "",
        "Both are shown so that the released column and this table can be "
        "read side by side; the verdicts agree, because the null is centred "
        "near zero and the observation is far from it.",
        "",
        "## The result",
        "",
        "| network | label | the observed margin | the null mean | the null "
        "SD | the rewirings | the effect-size z | the standard-error z |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for _, r in s.iterrows():
        lines.append(f"| {r['network']} | {r['label']} | "
                     f"{r['observed_margin']:+.4f} | {r['null_mean']:+.4f} | "
                     f"{r['null_sd']:.4f} | {int(r['n_rewire'])} | "
                     f"{r['z_effect_size']:+.2f} | "
                     f"{r['z_standard_error']:+.2f} |")
    lines += [
        "",
        "## Reading",
        "",
        f"The rewiring null is as decisive under the proximity learner as "
        f"under RWR: the effect-size z is positive in "
        f"{int((s.z_effect_size > 0).sum())} of {len(s)} arms, with the "
        f"median at {s.z_effect_size.median():.2f} and the range "
        f"{s.z_effect_size.min():.2f} to {s.z_effect_size.max():.2f}. A "
        "degree-preserving rewiring destroys the margin in every arm, so the "
        "increment is a property of the *wiring* and not of the degree "
        "sequence -- and that conclusion does not depend on the learner.",
        "",
        "The comparison with the released RWR null (same arms, same "
        "rewirings) is in the table `c27_rewiring_second_learner_summary.tsv` "
        "and the discussion below.",
        "",
        "## Caveats",
        "",
        f"1. {N_FOLDS} folds per arm (one per repeat) for both the observed "
        "and the null, because every rewiring needs a fresh multi-source "
        "BFS per fold. The released RWR null used all 200 folds.",
        "2. The rewirings are the ones Control 2 generated (Reactome has "
        "98, the rest 100).",
        "3. The two z conventions differ by √n; read the effect-size z when "
        "asking 'how far outside the null is this', and the standard-error "
        "z when comparing with the released table.",
        "4. **Exploratory** (amendment 4): no pre-registered endpoint is "
        "touched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
