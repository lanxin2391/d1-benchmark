"""Control 38 -- rescuing H2: does the provenance class explain the margin
once the network's size and coverage are taken out?

H2 ("the networks of one provenance class behave alike") was stated as
**not testable** on the released panel, because provenance is confounded
with size and coverage and six networks cannot separate the two. The
collection of item 17 adds four networks, and this control asks the
question in a form that the confounding cannot spoil:

    1. the released framing: a Kruskal-Wallis test of the per-network
       margin across the three provenance classes;
    2. the residualised framing: the margin is regressed on the network's
       size (log nodes) and its coverage (the share of the label's
       positives inside the network), and the class effect is tested on
       the **residuals** -- the same trick Control 33 uses for the score;
    3. the variance shares: how much of the margin's variance the class
       explains, how much the size and coverage explain, and how much is
       left over.

Outputs:
  results/tables/c38_h2_rescue.tsv
  results/tables/c38_h2_summary.tsv
  docs/c38_h2_rescue.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

TABLES = ROOT / "results/tables"
NETS_DIR = ROOT / "data/processed/networks"
LABELS_DIR = ROOT / "data/processed/labels"
DOC = ROOT / "docs/c38_h2_rescue.md"
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
PROVENANCE = {
    "funmap": "curated-literature", "reactome": "curated-literature",
    "intact": "curated-literature",
    "string_phys700": "protein", "string_full700": "protein",
    "string_full400": "protein", "string_full900": "protein",
    "rna_coexp": "rna", "rna_coexp_d150k": "rna",
    "rna_coexp_d250k": "rna",
    # the fourth curated member, size-matched to FunMap (item 17)
    "biogrid_curated": "curated-literature",
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def network_genes(net: str) -> set[str]:
    p = NETS_DIR / f"{net}.tsv"
    genes: set[str] = set()
    with p.open(encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genes.add(parts[0])
            genes.add(parts[1])
    return genes


def load_margins() -> pd.DataFrame:
    rel = pd.read_csv(TABLES / "a4_table1_delta_auroc_v40.tsv", sep="\t")
    rel = rel[["network", "label_set", "a0.5_mean"]].rename(
        columns={"label_set": "label", "a0.5_mean": "margin"})
    rel["panel"] = "the released panel"
    extra_path = TABLES / "a4_table1_delta_auroc_extra.tsv"
    if extra_path.exists():
        ex = pd.read_csv(extra_path, sep="\t")
        ex = ex[["network", "label_set", "a0.5_mean"]].rename(
            columns={"label_set": "label", "a0.5_mean": "margin"})
        ex["panel"] = "the extra networks"
        rel = pd.concat([rel, ex], ignore_index=True)
    return rel


def main() -> int:
    df = load_margins()
    df = df[df.margin.notna()].copy()
    # the network-level covariates
    cov_rows = []
    for net in sorted(df.network.unique()):
        genes = network_genes(net)
        edges = sum(1 for _ in (NETS_DIR / f"{net}.tsv").open(encoding="utf-8")) - 1
        cover = []
        for lab in LABELS:
            p = LABELS_DIR / f"{lab}.tsv"
            if not p.exists():
                continue
            pos = set(pd.read_csv(p, sep="\t", dtype=str)["gene"])
            if pos:
                cover.append(len(pos & genes) / len(pos))
        cov_rows.append({
            "network": net, "provenance": PROVENANCE.get(net, "unknown"),
            "n_nodes": len(genes), "n_edges": edges,
            "coverage": float(np.mean(cover)) if cover else np.nan})
    cov = pd.DataFrame(cov_rows)
    df = df.merge(cov, on="network", how="left")
    df.to_csv(TABLES / "c38_h2_rescue.tsv", sep="\t", index=False)

    # the per-network score: the median margin over the labels
    per = (df.groupby(["network", "provenance", "panel"])
           .agg(margin=("margin", "median"), n_nodes=("n_nodes", "first"),
                n_edges=("n_edges", "first"), coverage=("coverage", "first"),
                n_labels=("margin", "size")).reset_index())

    def kw(sub: pd.DataFrame) -> tuple[float, float]:
        groups = [g.margin.values for _, g in sub.groupby("provenance")
                  if len(g)]
        if len(groups) < 2:
            return np.nan, np.nan
        h, p = sps.kruskal(*groups)
        return float(h), float(p)

    rows = []
    for panel, sub in list(per.groupby("panel")) + [("the extended panel", per)]:
        h_raw, p_raw = kw(sub)
        # the residualised margin: the size and the coverage out
        X = np.column_stack([np.ones(len(sub)), np.log10(sub.n_nodes.values),
                             sub.coverage.values])
        y = sub.margin.values
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        sub = sub.copy()
        sub["residual"] = resid
        h_res, p_res = kw(sub)
        # the variance shares
        ss_tot = float(((y - y.mean()) ** 2).sum())
        ss_class = float(sum(len(g) * (g.margin.mean() - y.mean()) ** 2
                             for _, g in sub.groupby("provenance")))
        r2_size = 1 - float(((y - X @ beta) ** 2).sum()) / ss_tot
        rows.append({
            "panel": panel, "n_networks": len(sub),
            "kruskal_h_raw": h_raw, "kruskal_p_raw": p_raw,
            "kruskal_h_residual": h_res, "kruskal_p_residual": p_res,
            "variance_share_class": ss_class / ss_tot if ss_tot else np.nan,
            "r2_size_and_coverage": r2_size,
            "residual_sd": float(resid.std(ddof=1)) if len(resid) > 1 else np.nan,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(TABLES / "c38_h2_summary.tsv", sep="\t", index=False)

    # the per-class medians on the residuals
    lines = []
    for panel, sub in list(per.groupby("panel")) + [("the extended panel", per)]:
        X = np.column_stack([np.ones(len(sub)), np.log10(sub.n_nodes.values),
                             sub.coverage.values])
        beta, *_ = np.linalg.lstsq(X, sub.margin.values, rcond=None)
        sub = sub.copy()
        sub["residual"] = sub.margin.values - X @ beta
        for prov, g in sub.groupby("provenance"):
            lines.append((panel, prov, len(g), g.margin.median(),
                          g.residual.median()))
    pd.DataFrame(lines, columns=["panel", "provenance", "n_networks",
                                 "the median margin",
                                 "the median residual"]
                 ).to_csv(TABLES / "c38_h2_class_medians.tsv", sep="\t",
                          index=False)
    write_doc(per, summary, pd.DataFrame(
        lines, columns=["panel", "provenance", "n_networks",
                        "the median margin", "the median residual"]))
    log(f"wrote the Control 38 tables ({len(per)} networks)")
    return 0


def write_doc(per: pd.DataFrame, summary: pd.DataFrame,
              cls: pd.DataFrame) -> None:
    lines = [
        "# Control 38 — rescuing H2: the class effect once the size and the "
        "coverage are out",
        "",
        "**Date:** 2026-11-27  ",
        "**Script:** `scripts/run_c38_h2_rescue.py`  ",
        "**Tables:** `c38_h2_rescue.tsv`, `c38_h2_summary.tsv`, "
        "`c38_h2_class_medians.tsv`",
        "",
        "## The question",
        "",
        "H2 (\"the networks of one provenance class behave alike\") was "
        "stated as **not testable** on the released panel: provenance is "
        "confounded with size and coverage, and six networks cannot "
        "separate the two. The collection of item 17 adds four networks, "
        "which does not remove the confounding, but it does allow the "
        "confound to be removed **arithmetically**: the margin can be "
        "regressed on the network's size and its coverage, and the class "
        "effect can then be tested on what is left.",
        "",
        "## The result",
        "",
        "| panel | the networks | the Kruskal-Wallis p on the raw margins | "
        "the p on the size- and coverage-adjusted margins | the variance "
        "share of the class | the R2 of the size and the coverage |",
        "|---|---|---|---|---|---|",
    ]
    for _, r in summary.iterrows():
        lines.append(f"| {r['panel']} | {int(r['n_networks'])} | "
                     f"{r['kruskal_p_raw']:.3f} | "
                     f"{r['kruskal_p_residual']:.3f} | "
                     f"{r['variance_share_class']:.3f} | "
                     f"{r['r2_size_and_coverage']:.3f} |")
    lines += [
        "",
        "## The per-class medians",
        "",
        "| panel | the class | the networks | the median margin | the median "
        "adjusted margin |",
        "|---|---|---|---|---|",
    ]
    for _, r in cls.iterrows():
        lines.append(f"| {r['panel']} | {r['provenance']} | "
                     f"{int(r['n_networks'])} | {r['the median margin']:+.4f} "
                     f"| {r['the median residual']:+.4f} |")
    lines += [
        "",
        "## Reading",
        "",
        "## Caveats",
        "",
        "1. The adjustment is a linear regression on log10(nodes) and the "
        "coverage with ten (or six) networks: three parameters on ten "
        "observations leaves little room, and the adjusted test is a "
        "diagnostic, not a confirmation.",
        "2. The extra networks share their sources with the released ones "
        "(the same STRING download, the same TCGA expression), so the "
        "within-class variation is *source* variation, not independent "
        "replication.",
        "3. The class labels are three categories of one judgement; the "
        "test is a Kruskal-Wallis on ten numbers.",
        "4. **Exploratory** (amendment 4): no pre-registered endpoint is "
        "touched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
