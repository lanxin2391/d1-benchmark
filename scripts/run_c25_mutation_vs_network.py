"""Control 25 -- mutation-first versus network-first.

The teacher's question after the defence was whether the network prior
adds anything to the *mutation* evidence. Two experiments, both on data
already in the repository (no new engine run).

25a  The head-to-head as a prioritisation strategy.
     The answer key is the NCG-known catalogue, which is curated
     independently of IntOGen. On the same test genes and the same
     folds we compare four scorers:
         the mutation score  -- the IntOGen 2024 mutation evidence
                                (the breadth: the number of cohorts in
                                which the gene is a driver; and the depth:
                                the number of samples),
         the network prior   -- RWR at alpha = 0.5, seeded with the
                                training positives,
         degree,
         the combination     -- the rank sum of the prior and the mutation
                                score.
     Because IntOGen and the NCG overlap (Jaccard 0.522), the same
     comparison is run with the IntOGen driver genes removed from the
     test universe, which removes the direct circularity.

25b  Does the network know which drivers are recurrently mutated?
     Among the held-out IntOGen 2024 positives, the prior's score is
     compared for the recurrent drivers (seen in at least three samples
     in some cohort) and the single-mutation ones, with the degree and
     the publication count as the comparators.

Outputs:
  results/tables/c25_mutation_vs_network.tsv
  results/tables/c25_recurrence.tsv
  results/tables/c25_summary.tsv
  docs/c25_mutation_vs_network.md
"""
from __future__ import annotations

import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/c25_mutation_vs_network.md"
SCORES = ROOT / "results/scores_v40"
NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
ALPHA = 0.5
INTOGEN_ZIP = (ROOT / "data/raw/labels/intogen/"
               "IntOGen-Drivers-20240920.zip")
INTOGEN_MEMBER = "2024-06-18_IntOGen-Drivers/Unfiltered_drivers.tsv"


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


def auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """The AUROC of a score for the positive class (ties count as 0.5)."""
    pos = scores[labels > 0]
    neg = scores[labels <= 0]
    if pos.size < 2 or neg.size < 2:
        return np.nan
    r = tie_ranks(scores)
    n_pos = float((labels > 0).sum())
    n_neg = float((labels <= 0).sum())
    sum_pos = r[labels > 0].sum()
    return float((sum_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def unit_ranks(x: np.ndarray) -> np.ndarray:
    """The ranks scaled to [0, 1] (for the rank-sum combination)."""
    r = tie_ranks(x)
    return (r - 1.0) / max(len(x) - 1.0, 1.0)


def mutation_evidence() -> pd.DataFrame:
    """The per-gene IntOGen 2024 mutation evidence."""
    z = zipfile.ZipFile(INTOGEN_ZIP)
    raw = z.read(INTOGEN_MEMBER).decode("utf-8", errors="replace")
    it = pd.read_csv(pd.io.common.StringIO(raw), sep="\t", low_memory=False)
    drv = it[it["DRIVER"] == True]          # noqa: E712
    agg = drv.groupby("SYMBOL").agg(
        n_cohorts=("COHORT", "nunique"),
        n_samples=("SAMPLES_COHORT", "sum"),
        n_recurrent_samples=("SAMPLES_3MUTS", "sum"),
        best_tier=("TIER", "min"),
        n_rows=("SYMBOL", "size"),
    ).reset_index()
    rec = (it.groupby("SYMBOL")["SAMPLES_3MUTS"].max() > 0)
    agg = agg.merge(rec.rename("recurrent").reset_index(), on="SYMBOL",
                    how="left")
    agg["recurrent"] = agg["recurrent"].fillna(False)
    return agg


def main() -> int:
    log("reading the IntOGen 2024 mutation evidence")
    ev = mutation_evidence()
    log(f"  {len(ev)} driver genes in the IntOGen 2024 release "
        f"({int(ev['recurrent'].sum())} of them recurrent)")

    pc = pd.read_csv(ROOT / "data/processed/labels/pubcount.tsv", sep="\t")
    pub = dict(zip(pc["gene"], pc["log10"]))

    mut_rows = []
    rec_rows = []
    for net in NETWORKS:
        fp = SCORES / f"{net}__ncg_known.npz"
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
        gidx = {g: i for i, g in enumerate(genes)}
        mut_coh = np.zeros(len(genes))
        mut_samp = np.zeros(len(genes))
        for g, c, s in zip(ev["SYMBOL"], ev["n_cohorts"], ev["n_samples"]):
            i = gidx.get(g)
            if i is not None:
                mut_coh[i] = c
                mut_samp[i] = s
        # the overlap set: the IntOGen 2024 driver label as the project
        # builds it (the same 623 genes everywhere else)
        label = pd.read_csv(ROOT / "data/processed/labels/intogen2024.tsv",
                            sep="\t")
        in_intogen = np.array([g in set(label["gene"]) for g in genes])
        log(f"{net}: {len(genes)} genes, {int(y.sum())} positives, "
            f"{int(in_intogen.sum())} in the IntOGen driver set")
        for j in range(in_test.shape[0]):
            idx = np.where(in_test[j])[0]
            lab = y[idx]
            if (lab > 0).sum() < 5 or (lab <= 0).sum() < 20:
                continue
            s_prior = rwr[idx, j]
            s_deg = degree[idx]
            s_m1 = mut_coh[idx]
            s_m2 = mut_samp[idx]
            s_comb = unit_ranks(s_prior) + unit_ranks(s_m1)
            s_comb2 = unit_ranks(s_prior) + unit_ranks(s_m2)
            for universe in ("all", "minus_intogen"):
                keep = np.ones(len(idx), dtype=bool)
                if universe == "minus_intogen":
                    keep = ~in_intogen[idx]
                    if (lab[keep] > 0).sum() < 5:
                        continue
                for tag, sc in (("the prior", s_prior),
                                ("degree", s_deg),
                                ("the mutation score (the cohorts)",
                                 s_m1),
                                ("the mutation score (the samples)",
                                 s_m2),
                                ("the prior + the mutation score",
                                 s_comb),
                                ("the prior + the samples score", s_comb2)):
                    a = auc(sc[keep], lab[keep])
                    if np.isfinite(a):
                        mut_rows.append({
                            "network": net, "fold": j,
                            "universe": universe, "scorer": tag,
                            "auc": a,
                            "n_test": int(keep.sum()),
                            "n_pos": int((lab[keep] > 0).sum()),
                        })
        log(f"  {net}: the head-to-head is done")

    # ---- 25b: recurrent versus single-mutation drivers -------------------
    for net in NETWORKS:
        fp = SCORES / f"{net}__intogen2024.npz"
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
        gidx = {g: i for i, g in enumerate(genes)}
        recurrent = np.zeros(len(genes), dtype=bool)
        for g, r in zip(ev["SYMBOL"], ev["recurrent"]):
            i = gidx.get(g)
            if i is not None:
                recurrent[i] = bool(r)
        pubcount = np.array([pub.get(g, 0.0) for g in genes])
        for j in range(in_test.shape[0]):
            idx = np.where(in_test[j])[0]
            pos = idx[y[idx] == 1]
            if len(pos) < 10:
                continue
            lab = recurrent[pos].astype(int)
            if (lab == 1).sum() < 5 or (lab == 0).sum() < 5:
                continue
            for tag, sc in (("the prior", rwr[pos, j]),
                            ("degree", degree[pos]),
                            ("the publication count", pubcount[pos])):
                a = auc(sc, lab)
                if np.isfinite(a):
                    rec_rows.append({
                        "network": net, "fold": j, "scorer": tag,
                        "auc": a, "n_positives": len(pos),
                        "n_recurrent": int((lab == 1).sum()),
                    })
        log(f"  {net}: the recurrence comparison is done")

    m = pd.DataFrame(mut_rows)
    m.to_csv(TABLES / "c25_mutation_vs_network.tsv", sep="\t", index=False)
    r = pd.DataFrame(rec_rows)
    r.to_csv(TABLES / "c25_recurrence.tsv", sep="\t", index=False)
    summary = summarise(m, r)
    summary.to_csv(TABLES / "c25_summary.tsv", sep="\t", index=False)
    write_doc(m, r, summary)
    log("wrote the Control 25 tables")
    return 0


def summarise(m: pd.DataFrame, r: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for universe in sorted(m.universe.unique()):
        sub = m[m.universe == universe]
        for net, grp in sub.groupby("network"):
            med = grp.groupby("scorer").auc.median()
            rows.append({"part": "25a", "universe": universe, "network": net,
                         **{f"auc_{k}": round(v, 4) for k, v in med.items()},
                         "n_folds": grp.fold.nunique()})
            rows.append({"part": "25a", "universe": universe, "network": net,
                         "prior_minus_mutation_cohorts":
                             round(med.get("the prior", np.nan)
                                   - med.get("the mutation score (the cohorts)",
                                             np.nan), 4),
                         "combination_minus_prior":
                             round(med.get("the prior + the mutation score",
                                           np.nan) - med.get("the prior",
                                                             np.nan), 4),
                         "n_folds": grp.fold.nunique()})
    for net, grp in r.groupby("network"):
        med = grp.groupby("scorer").auc.median()
        rows.append({"part": "25b", "universe": "intogen2024",
                     "network": net,
                     **{f"auc_{k}": round(v, 4) for k, v in med.items()},
                     "n_folds": grp.fold.nunique()})
    return pd.DataFrame(rows)


def write_doc(m: pd.DataFrame, r: pd.DataFrame, summary: pd.DataFrame) -> None:
    lines = [
        "# Control 25 — mutation-first versus network-first",
        "",
        "**Date:** 2026-11-18  ",
        "**Script:** `scripts/run_c25_mutation_vs_network.py`  ",
        "**Tables:** `c25_mutation_vs_network.tsv`, `c25_recurrence.tsv`, "
        "`c25_summary.tsv`",
        "",
        "## The question",
        "",
        "After the defence the question left open was whether the network "
        "prior adds anything to the **mutation** evidence. Two experiments, "
        "both on data already in the repository.",
        "",
        "**25a, the head-to-head.** The answer key is the NCG-known "
        "catalogue, curated independently of IntOGen. On the same test "
        "genes and the same folds, four scorers are compared: the IntOGen "
        "2024 mutation evidence (the number of cohorts in which the gene "
        "is a driver, and the number of samples), the network prior (RWR "
        "at alpha = 0.5, seeded with the training positives), degree, and "
        "the rank sum of the prior and the mutation score. Because "
        "IntOGen and the NCG overlap, the same comparison is run with the "
        "IntOGen driver genes removed from the test universe.",
        "",
        "**25b, the recurrence comparison.** Among the held-out IntOGen "
        "2024 positives, the prior's score is compared for the recurrently "
        "mutated drivers (seen in at least three samples in some cohort) "
        "and the single-mutation ones, with degree and the publication "
        "count as the comparators.",
        "",
        "## 25a: the AUROC on the NCG-known label (the median over the "
        "folds)",
        "",
        "| network | universe | the prior | degree | the mutation score "
        "(the cohorts) | the mutation score (the samples) | the prior + the "
        "mutation score |",
        "|---|---|---|---|---|---|---|",
    ]
    for universe in sorted(m.universe.unique()):
        for net in sorted(m[m.universe == universe].network.unique()):
            sub = m[(m.universe == universe) & (m.network == net)]
            med = sub.groupby("scorer").auc.median()
            lines.append(
                f"| {net} | {universe} | "
                f"{med.get('the prior', float('nan')):.4f} | "
                f"{med.get('degree', float('nan')):.4f} | "
                f"{med.get('the mutation score (the cohorts)', float('nan')):.4f} | "
                f"{med.get('the mutation score (the samples)', float('nan')):.4f} | "
                f"{med.get('the prior + the mutation score', float('nan')):.4f} |")
    lines += [
        "",
        "## 25b: does the network know which drivers are recurrently "
        "mutated? (the AUROC for 'is recurrent' among the held-out "
        "positives)",
        "",
        "| network | the prior | degree | the publication count |",
        "|---|---|---|---|",
    ]
    for net in sorted(r.network.unique()):
        sub = r[r.network == net]
        med = sub.groupby("scorer").auc.median()
        lines.append(
            f"| {net} | {med.get('the prior', float('nan')):.4f} | "
            f"{med.get('degree', float('nan')):.4f} | "
            f"{med.get('the publication count', float('nan')):.4f} |")
    lines += [
        "",
        "## Caveats",
        "",
        "1. The mutation score is a static per-gene quantity; the prior is "
        "fold-dependent. Both are evaluated on the same test genes in the "
        "same folds.",
        "2. The NCG catalogue and IntOGen overlap (the Jaccard is 0.522), "
        "so the `minus_intogen` universe is the one to read.",
        "3. The recurrence status is the IntOGen annotation of the gene; a "
        "gene being recurrent is also a gene that has been looked at more, "
        "which is why the publication count is the comparator to read "
        "alongside.",
        "4. **Exploratory** (amendment 4): nothing in the pre-registered "
        "endpoint is touched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
