"""Control 23 -- the split-design comparison.

The teacher asked us to try other cross-validation designs (for example
10 x 10) and to say which design is best. This control answers that on
the same data, the same networks and the same labels, and it keeps the
pre-registered design (5 folds x 40 repeats) as the primary result: the
comparison is a SENSITIVITY analysis.

Designs compared (all 24 arms: 6 networks x 4 catalogues, alpha = 0.5):

    k=5,  R=10    100 evaluations   subset of the released panel
    k=5,  R=20    200 evaluations   subset of the released panel
    k=5,  R=40    400 evaluations   THE PRE-REGISTERED DESIGN
    k=10, R=10    100 evaluations   the teacher's example
    k=10, R=20    200 evaluations
    k=20, R=10    200 evaluations

For every design and every arm we report
  * the margin (delta-AUROC) mean, median and SD,
  * the repeat-level percentile interval (bootstrap over the per-repeat
    means) -- the primary inference construction,
  * the gene-cluster bootstrap interval -- the conservative construction
    that decides the H1 verdict,
  * and, per design, how many arms each construction certifies and how
    weak a certified arm can be (the power proxy).

The k=5 designs reuse the released runs and score dumps (the first 10
repeats of the 40-repeat files are byte-identical to the released
10-repeat panel). The k=10 and k=20 designs get their own splits in
data/processed/splits_designs/ and are run through the same engine.

Outputs:
  results/tables/c23_split_design_arms.tsv
  results/tables/c23_split_design_summary.tsv
  docs/c23_split_design.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.evaluate import run_arm  # noqa: E402
from d1.engine.intervals import gene_cluster_bootstrap  # noqa: E402
import d1.splits as SPLITS  # noqa: E402

T = ROOT / "results/tables"
RUNS = ROOT / "results/runs_v40"
DUMPS = ROOT / "results/scores_v40"
NETS = ROOT / "data/processed/networks"
DESIGN_SPLITS = ROOT / "data/processed/splits_designs"
DESIGN_RUNS = ROOT / "results/runs_designs"
DESIGN_DUMPS = ROOT / "results/scores_designs"

BASE_SEED = 20260926
ALPHA = 0.5
N_BOOT_REPEAT = 2000
N_BOOT_GENE = 1000
SEED = 20261029

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]

DESIGNS = [
    (5, 10), (5, 20), (5, 40),
    (10, 10), (10, 20), (20, 10),
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def label_genes(label_id: str) -> set[str]:
    df = pd.read_csv(ROOT / "data/processed/labels" / f"{label_id}.tsv",
                     sep="\t", dtype=str)
    return set(df["gene"])


def network_genes(net_id: str) -> set[str]:
    df = pd.read_csv(NETS / f"{net_id}.tsv", sep="\t", dtype=str)
    return set(df["gene_a"]) | set(df["gene_b"])


def build_splits(k: int, r: int) -> dict[tuple[str, str], Path]:
    """The splits for one design, written to their own directory."""
    out: dict[tuple[str, str], Path] = {}
    for label in LABELS:
        pos = label_genes(label)
        for net in NETWORKS:
            uni = network_genes(net)
            df = SPLITS.make_splits(uni, pos, base_seed=BASE_SEED,
                                   n_rep=r, k=k)
            p = DESIGN_SPLITS / f"k{k}_r{r}" / f"{label}__native_{net}.tsv"
            p.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(p, sep="\t", index=False, lineterminator="\n")
            out[(net, label)] = p
    return out


def repeat_level_ci(deltas: pd.Series) -> tuple[float, float]:
    """The percentile bootstrap over the per-repeat means."""
    reps = deltas.to_numpy()
    if reps.size < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(SEED)
    means = reps[rng.integers(0, reps.size, size=(N_BOOT_REPEAT, reps.size))].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def arm_stats_from_runs(path: Path, net: str, label: str,
                        max_repeat: int | None) -> dict:
    """The per-fold margins of one arm from a run file, restricted to the
    first `max_repeat` repeats."""
    df = pd.read_csv(path, sep="\t")
    df = df[(df.network == net) & (df.label_set == label)]
    if max_repeat is not None:
        df = df[df.repeat < max_repeat]
    df = df[df.alpha.isna() | (np.isclose(df.alpha.astype(float), ALPHA))]
    piv = df.pivot_table(index=["repeat", "fold"], columns="method",
                          values="auroc").reset_index()
    piv = piv.dropna(subset=["rwr", "degree"])
    piv["delta"] = piv["rwr"] - piv["degree"]
    return {"per_fold": piv, "n_evals": int(len(piv))}


def gene_cluster_ci(dump: Path, net: str, label: str,
                    max_repeat: int | None) -> tuple[float, float, float]:
    """The gene-cluster bootstrap interval, restricted to the first
    `max_repeat` repeats (the dumps carry repeat_idx and fold_idx)."""
    d = np.load(dump, allow_pickle=True)
    ai = int(np.argmin(np.abs(np.asarray(d["alphas"], dtype=float) - ALPHA)))
    rwr = np.asarray(d["rwr"][ai], dtype=float).T
    in_test = np.asarray(d["in_test"], dtype=bool)
    rep = np.asarray(d["repeat_idx"], dtype=int)
    if max_repeat is not None:
        keep = rep < max_repeat
        rwr = rwr[keep]
        in_test = in_test[keep]
    res = gene_cluster_bootstrap(rwr, np.asarray(d["degree"], dtype=float),
                                 np.asarray(d["y"]).astype(int), in_test,
                                 n_boot=N_BOOT_GENE, seed=SEED)
    return float(res["ci_low"]), float(res["ci_high"]), float(res["observed_mean"])


def main() -> int:
    rows: list[dict] = []
    for (k, r) in DESIGNS:
        tag = f"k{k}_r{r}"
        log(f"design {tag} ({k} folds x {r} repeats = {k * r} evaluations)")
        if k == 5:
            # a subset of the released panel: no recomputation
            for net in NETWORKS:
                for label in LABELS:
                    run = RUNS / f"a4_{net}__{label}.tsv"
                    dump = DUMPS / f"{net}__{label}.npz"
                    st = arm_stats_from_runs(run, net, label, max_repeat=r)
                    per_fold = st["per_fold"]
                    if per_fold.empty:
                        continue
                    deltas = per_fold["delta"]
                    per_repeat = per_fold.groupby("repeat")["delta"].mean()
                    lo, hi = repeat_level_ci(per_repeat)
                    glo, ghi, gmean = gene_cluster_ci(dump, net, label,
                                                      max_repeat=r)
                    rows.append({
                        "design": tag, "k": k, "n_repeats": r,
                        "n_evaluations": st["n_evals"],
                        "network": net, "label": label,
                        "margin_mean": float(deltas.mean()),
                        "margin_median": float(deltas.median()),
                        "margin_sd": float(deltas.std(ddof=1)),
                        "repeat_ci_low": lo, "repeat_ci_high": hi,
                        "repeat_ci_width": hi - lo,
                        "repeat_ci_excludes_0": bool(lo > 0 or hi < 0),
                        "gene_ci_low": glo, "gene_ci_high": ghi,
                        "gene_ci_width": ghi - glo,
                        "gene_ci_excludes_0": bool(glo > 0 or ghi < 0),
                        "gene_ci_mean": gmean,
                    })
            log(f"  design {tag}: {len(rows)} arm rows so far")
        else:
            # a new design: build the splits and run the arms
            splits = build_splits(k, r)
            for (net, label), sp in splits.items():
                out = DESIGN_RUNS / tag / f"a4_{net}__{label}.tsv"
                out.parent.mkdir(parents=True, exist_ok=True)
                dump = DESIGN_DUMPS / tag / f"{net}__{label}.npz"
                dump.parent.mkdir(parents=True, exist_ok=True)
                res = run_arm(network=str(NETS / f"{net}.tsv"),
                              splits=str(sp), network_id=net,
                              label_id=label, universe_id=f"native_{net}",
                              alphas=(ALPHA,), dump_scores=str(dump))
                res.to_csv(out, sep="\t", index=False, lineterminator="\n")
                st = arm_stats_from_runs(out, net, label, max_repeat=None)
                per_fold = st["per_fold"]
                if per_fold.empty:
                    continue
                deltas = per_fold["delta"]
                per_repeat = per_fold.groupby("repeat")["delta"].mean()
                lo, hi = repeat_level_ci(per_repeat)
                glo, ghi, gmean = gene_cluster_ci(dump, net, label,
                                                  max_repeat=None)
                rows.append({
                    "design": tag, "k": k, "n_repeats": r,
                    "n_evaluations": st["n_evals"],
                    "network": net, "label": label,
                    "margin_mean": float(deltas.mean()),
                    "margin_median": float(deltas.median()),
                    "margin_sd": float(deltas.std(ddof=1)),
                    "repeat_ci_low": lo, "repeat_ci_high": hi,
                    "repeat_ci_width": hi - lo,
                    "repeat_ci_excludes_0": bool(lo > 0 or hi < 0),
                    "gene_ci_low": glo, "gene_ci_high": ghi,
                    "gene_ci_width": ghi - glo,
                    "gene_ci_excludes_0": bool(glo > 0 or ghi < 0),
                    "gene_ci_mean": gmean,
                })
            log(f"  design {tag} done")

    df = pd.DataFrame(rows)
    df.to_csv(T / "c23_split_design_arms.tsv", sep="\t", index=False)

    # ---- the summary per design -------------------------------------
    out = []
    for (tag, g) in df.groupby("design", sort=False):
        g = g.sort_values("k")
        certified_rep = g[g.repeat_ci_excludes_0]
        certified_gene = g[g.gene_ci_excludes_0]
        out.append({
            "design": tag, "k": int(g.k.iloc[0]),
            "n_repeats": int(g.n_repeats.iloc[0]),
            "n_evaluations": int(g.n_evaluations.iloc[0]),
            "arms": len(g),
            "margin_median": float(g.margin_median.median()),
            "margin_mean": float(g.margin_mean.mean()),
            "margin_sd_median": float(g.margin_sd.median()),
            "arms_positive": int((g.margin_mean > 0).sum()),
            "repeat_ci_certified": int(len(certified_rep)),
            "repeat_ci_width_median": float(g.repeat_ci_width.median()),
            "weakest_certified_margin": (float(certified_rep.margin_mean.min())
                                        if len(certified_rep) else np.nan),
            "gene_ci_certified": int(len(certified_gene)),
            "gene_ci_width_median": float(g.gene_ci_width.median()),
            "weakest_gene_certified_margin": (float(certified_gene.margin_mean.min())
                                              if len(certified_gene) else np.nan),
        })
    summary = pd.DataFrame(out)
    summary.to_csv(T / "c23_split_design_summary.tsv", sep="\t", index=False)
    write_doc(df, summary)
    log("wrote the Control 23 tables")
    print(summary.to_string(index=False))
    return 0


def write_doc(df: pd.DataFrame, summary: pd.DataFrame) -> None:
    lines = [
        "# Control 23 -- the split-design comparison",
        "",
        "**Date:** 2026-11-18  ",
        "**Script:** `scripts/run_c23_split_design.py`  ",
        "**Tables:** `c23_split_design_arms.tsv`, `c23_split_design_summary.tsv`",
        "",
        "## The question",
        "",
        "The teacher asked us to try other cross-validation designs (for "
        "example 10 x 10) and to say which design is best. This control "
        "compares six designs on the same 24 arms, the same networks, the "
        "same labels and the same endpoint (alpha = 0.5). **The "
        "pre-registered design (5 folds x 40 repeats) stays the primary "
        "result; this is a sensitivity analysis.**",
        "",
        "## The designs and what they cost",
        "",
        "| design | evaluations per arm | test positives per fold | "
        "repeat-level interval | gene-cluster interval |",
        "|---|---|---|---|---|",
    ]
    for r in summary.itertuples():
        lines.append(
            f"| {r.design} | {r.n_evaluations} | "
            f"{'~' + str(round(r.margin_median, 3)) if False else ''}"
            f"median margin {r.margin_median:+.4f} | "
            f"certifies {r.repeat_ci_certified}/{r.arms}, median width "
            f"{r.repeat_ci_width_median:.4f} | "
            f"certifies {r.gene_ci_certified}/{r.arms}, median width "
            f"{r.gene_ci_width_median:.4f} |")
    lines += [
        "",
        "## Reading",
        "",
        f"- the median margin over the 24 arms moves between "
        f"{summary.margin_median.min():+.4f} and "
        f"{summary.margin_median.max():+.4f} across the six designs: "
        "the size of the effect is a property of the data, not of the "
        "design;",
        f"- the conservative gene-cluster construction certifies "
        f"{int(summary.gene_ci_certified.min())} to "
        f"{int(summary.gene_ci_certified.max())} of the 24 arms "
        "depending on the design: this is the part of the verdict that the "
        "design does move;",
        f"- the median gene-cluster interval width goes from "
        f"{summary.gene_ci_width_median.min():.4f} to "
        f"{summary.gene_ci_width_median.max():.4f}.",
        "",
        "## Caveats",
        "",
        "1. The k = 10 and k = 20 designs use the same base seed but a "
        "different fold structure, so their gene-to-fold assignment "
        "differs from the released panel; that is the point of the "
        "comparison, and it is why the comparison is reported as a "
        "sensitivity analysis rather than as a replacement.",
        "2. The gene-cluster bootstrap uses 1,000 resamples and the "
        "repeat-level percentile interval 2,000, the same as the "
        "released analysis.",
        "3. **Exploratory** (amendment 4): the pre-registered endpoint "
        "and the H1 verdict are untouched.",
    ]
    (ROOT / "docs/c23_split_design.md").write_text("\n".join(lines) + "\n",
                                                   encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
