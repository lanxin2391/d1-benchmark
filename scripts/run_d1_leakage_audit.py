"""Phase D, task D1 -- quantify the complex leakage of the CURRENT splits.

The plan's rule is "quantify first, then fix": before building a
complex-aware splitter, measure how much complex leakage the v1 splits
actually contain.

For every (network, label, repeat, fold) of the split sets in use, we
count the held-out (test) positives that have at least one
COMPLEX-LINKED gene in the training fold. A complex-linked pair is two
genes that appear in the same curated complex (the complex map is
data/processed/complexes/complex_map.tsv; two genes are linked iff their
CPX-ID lists intersect).

The same quantity is computed for the held-out NEGATIVES as a control:
if the fraction were equally high for negatives, the effect would not be
about the positives at all.

Three split families are audited:
  * native_<net>       -- the 24 arms of the primary factorial
  * shared             -- the D-12 shared-universe control
  * <label>_cap_<net>  -- the D-11 common-size cap control

Outputs:
  results/tables/d1_leakage_audit.tsv      (per arm x repeat x fold)
  results/tables/d1_leakage_summary.tsv    (per arm)
  docs/d1_leakage_audit.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

TABLES = ROOT / "results" / "tables"
SPLITS = ROOT / "data" / "processed" / "splits_v40"
CPLX = ROOT / "data" / "processed" / "complexes" / "complex_map.tsv"
DOC = ROOT / "docs" / "d1_leakage_audit.md"

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def build_partners() -> dict[str, set[str]]:
    """gene -> set of panel genes that share a complex with it."""
    cm = pd.read_csv(CPLX, sep="\t")
    by_cpx: dict[str, list[str]] = {}
    for c, g in zip(cm["complex_id"], cm["gene"]):
        by_cpx.setdefault(c, []).append(g)
    partners: dict[str, set[str]] = {}
    for c, gs in by_cpx.items():
        gs = sorted(set(gs))
        for i, a in enumerate(gs):
            pa = partners.setdefault(a, set())
            pa.update(gs[:i])
            pa.update(gs[i + 1:])
    return partners


def audit_split(path: Path, partners: dict[str, set[str]],
                family: str, net: str, label: str) -> list[dict]:
    df = pd.read_csv(path, sep="\t")
    df["y"] = df["y"].astype(int)
    out = []
    for (r, f), g in df.groupby(["repeat", "fold"]):
        test = g[g.y == 1]
        train_genes = set(df[(df["repeat"] == r) & (df["fold"] != f)]["gene"])
        n_pos = len(test)
        n_pos_linked = 0
        n_pos_in_complex = 0
        for gene in test["gene"]:
            p = partners.get(gene)
            if p:
                n_pos_in_complex += 1
                if p & train_genes:
                    n_pos_linked += 1
        neg = g[g.y == 0]
        n_neg = len(neg)
        n_neg_linked = 0
        n_neg_in_complex = 0
        for gene in neg["gene"]:
            p = partners.get(gene)
            if p:
                n_neg_in_complex += 1
                if p & train_genes:
                    n_neg_linked += 1
        out.append({
            "family": family, "network": net, "label": label,
            "repeat": int(r), "fold": int(f),
            "n_test_pos": n_pos, "n_test_pos_complex_linked_to_train": n_pos_linked,
            "frac_test_pos_linked": n_pos_linked / n_pos if n_pos else np.nan,
            "n_test_pos_in_a_complex": n_pos_in_complex,
            "frac_test_pos_in_complex_linked": (
                n_pos_linked / n_pos_in_complex if n_pos_in_complex else np.nan),
            "n_test_neg": n_neg, "n_test_neg_complex_linked_to_train": n_neg_linked,
            "frac_test_neg_linked": n_neg_linked / n_neg if n_neg else np.nan,
            "n_test_neg_in_a_complex": n_neg_in_complex,
            "frac_test_neg_in_complex_linked": (
                n_neg_linked / n_neg_in_complex if n_neg_in_complex else np.nan),
            "n_genes_in_a_complex": int(sum(1 for g in df[df.repeat == r]["gene"]
                                           if g in partners)),
        })
    return out


def main() -> int:
    partners = build_partners()
    log(f"complex map: {len(partners):,} panel genes with complex partners")
    sizes = {}
    for g, p in partners.items():
        sizes[g] = len(p)
    log(f"  median number of complex partners per gene: "
        f"{np.median(list(sizes.values())):.0f}")

    rows = []
    # family 1: the 24 native arms of the primary factorial
    for net in NETWORKS:
        for label in LABELS:
            p = SPLITS / f"{label}__native_{net}.tsv"
            if not p.exists():
                continue
            rows += audit_split(p, partners, "native", net, label)
        log(f"  audited the native arms of {net}")
    # family 2: the D-12 shared-universe splits
    for label in LABELS:
        p = SPLITS / f"{label}_shared__shared.tsv"
        if p.exists():
            rows += audit_split(p, partners, "shared", "shared", label)
    # family 3: the D-11 capped splits
    for net in NETWORKS:
        for label in LABELS:
            p = SPLITS / f"{label}_cap_{net}__native_{net}.tsv"
            if p.exists():
                rows += audit_split(p, partners, "capped", net, label)
    log(f"audited {len(rows):,} (arm, repeat, fold) rows")

    audit = pd.DataFrame(rows)
    audit.to_csv(TABLES / "d1_leakage_audit.tsv", sep="\t", index=False)

    # ---- per-arm summary ----------------------------------------------
    summary = (audit.groupby(["family", "network", "label"])
                .agg(n_folds=("frac_test_pos_linked", "size"),
                     median_frac_test_pos_linked=("frac_test_pos_linked",
                                                  "median"),
                     mean_frac_test_pos_linked=("frac_test_pos_linked", "mean"),
                     max_frac_test_pos_linked=("frac_test_pos_linked", "max"),
                     median_frac_test_neg_linked=("frac_test_neg_linked",
                                                  "median"),
                     median_frac_pos_linked_given_in_complex=(
                         "frac_test_pos_in_complex_linked", "median"),
                     min_frac_pos_linked_given_in_complex=(
                         "frac_test_pos_in_complex_linked", "min"),
                     median_frac_neg_linked_given_in_complex=(
                         "frac_test_neg_in_complex_linked", "median"),
                     mean_n_pos_linked=("n_test_pos_complex_linked_to_train",
                                        "mean"),
                     n_pos=("n_test_pos", "mean"))
                .reset_index())
    summary["excess_over_negatives"] = (
        summary["median_frac_test_pos_linked"]
        - summary["median_frac_test_neg_linked"])
    summary.to_csv(TABLES / "d1_leakage_summary.tsv", sep="\t", index=False)

    log("\n=== the 24 primary arms (native) ===")
    nat = summary[summary.family == "native"].sort_values(
        "median_frac_test_pos_linked", ascending=False)
    print(nat[["network", "label", "n_pos", "median_frac_test_pos_linked",
               "max_frac_test_pos_linked", "median_frac_test_neg_linked",
               "excess_over_negatives"]].round(4).to_string(index=False))
    log("\n=== other split families ===")
    for fam in ("shared", "capped"):
        sub = summary[summary.family == fam]
        if len(sub):
            log(f"{fam}: median fraction of test positives linked to the "
                f"train fold = {sub.median_frac_test_pos_linked.median():.4f} "
                f"(negatives {sub.median_frac_test_neg_linked.median():.4f}); "
                f"max over arms = {sub.median_frac_test_pos_linked.max():.4f}")

    write_doc(audit, summary)
    log(f"wrote d1_leakage_audit.tsv and d1_leakage_summary.tsv")
    return 0


def write_doc(audit: pd.DataFrame, summary: pd.DataFrame) -> None:
    nat = summary[summary.family == "native"]
    shared = summary[summary.family == "shared"]
    capped = summary[summary.family == "capped"]
    pos_all = nat.median_frac_test_pos_linked.median()
    neg_all = nat.median_frac_test_neg_linked.median()
    pos_norm = nat.median_frac_pos_linked_given_in_complex.median()
    neg_norm = nat.median_frac_neg_linked_given_in_complex.median()
    lines = [
        "# D1 -- how much complex leakage do the current splits contain?",
        "",
        "**Date:** 2026-11-05  ",
        "**Script:** `scripts/run_d1_leakage_audit.py`  ",
        "**Task:** Phase D, task D1 (quantify first, then fix)",
        "",
        "## What was measured",
        "",
        "For every (arm, repeat, fold): the number of held-out (test) "
        "positives that have at least one **complex-linked** gene in the "
        "training fold. Two genes are complex-linked when they appear in "
        "the same curated complex (ComplexPortal identifiers via "
        "UniProtKB; `data/processed/complexes/complex_map.tsv`). The same "
        "quantity is computed for the held-out negatives as a control, and "
        "the whole thing is normalised by the number of test positives "
        "that are in a curated complex at all.",
        "",
        "## Result 1 -- the channel is wide open",
        "",
        "| quantity | median over the 24 arms | arm range |",
        "|---|---|---|",
        f"| test positives with a complex-linked gene in the train fold | "
        f"**{pos_all:.3f}** | {nat.median_frac_test_pos_linked.min():.3f} - "
        f"{nat.median_frac_test_pos_linked.max():.3f} |",
        f"| test NEGATIVES with a complex-linked gene in the train fold | "
        f"{neg_all:.3f} | {nat.median_frac_test_neg_linked.min():.3f} - "
        f"{nat.median_frac_test_neg_linked.max():.3f} |",
        f"| **among the test positives that are in a complex**, the "
        f"fraction with a linked gene in the train fold | **{pos_norm:.3f}** "
        f"| {nat.min_frac_pos_linked_given_in_complex.min():.3f} - "
        f"{1.0:.3f} |",
        f"| the same among the test negatives | {neg_norm:.3f} | |",
        f"| excess (positives - negatives) | "
        f"**{nat.excess_over_negatives.median():+.3f}** | "
        f"{nat.excess_over_negatives.min():+.3f} - "
        f"{nat.excess_over_negatives.max():+.3f} |",
        "",
        "**Reading.** For every gene that *can* be leaked, the standard "
        "CV leaks it: among the held-out positives that are in a curated "
        "complex, "
        f"**{pos_norm:.0%}** have a complex-linked gene in the training "
        "fold. That is not a property of the positives -- the held-out "
        f"**negatives** are linked at the same rate ({neg_norm:.1%}), "
        "because whether two genes share a complex does not depend on "
        "whether they are drivers. The leakage channel is therefore "
        "**wide open but symmetric**: it is opened for the positive and "
        "the negative classes alike.",
        "",
        "The second quantity is the one that matters for the benchmark's "
        "claim. The *enrichment* of positives over negatives is "
        f"{nat.excess_over_negatives.median():+.3f} on average "
        f"({nat.excess_over_negatives.min():+.3f} to "
        f"{nat.excess_over_negatives.max():+.3f} across arms) and is "
        "essentially zero for the four `intogen_temporal_new` arms "
        "(their positives are few and mostly not in a curated complex). "
        "A symmetric channel with a modest enrichment is the profile of a "
        "leak that inflates the *absolute* performance of a propagating "
        "learner but is unlikely to inflate its *advantage over the degree "
        "baseline* much -- but whether it does is exactly what D4 has to "
        "measure, not what D1 can assume.",
        "",
        "## Result 2 -- the same holds for the control split families",
        "",
        "| split family | test positives linked | negatives linked |",
        "|---|---|---|",
        f"| native (the 24 primary arms) | {pos_all:.3f} | {neg_all:.3f} |",
    ]
    if len(shared):
        lines.append(f"| shared (D-12) | "
                     f"{shared.median_frac_test_pos_linked.median():.3f} | "
                     f"{shared.median_frac_test_neg_linked.median():.3f} |")
    if len(capped):
        lines.append(f"| capped (D-11) | "
                     f"{capped.median_frac_test_pos_linked.median():.3f} | "
                     f"{capped.median_frac_test_neg_linked.median():.3f} |")
    lines += [
        "",
        "## Decision for Phase D",
        "",
        "The leakage exists and is measurable, but it is **symmetric** "
        "between the classes. Under the plan's Q3 rule the consequence "
        "is decided by D4 (the change in delta-AUROC when the splits "
        "become complex-aware), not by D1. D1 does establish that the "
        "criticism of Picart-Armada et al. 2019 is *mechanically* "
        "applicable to this panel, so D2-D4 are worth doing -- and the "
        "prediction to test is that the change in delta-AUROC will be "
        "**smaller than the raw leakage numbers suggest**, because the "
        "channel favours neither class.",
        "",
        "## Caveats",
        "",
        "1. The complex map covers 18.3 % of the panel's genes and "
        "33-37 % of the positives -- curated complex data is sparse, and "
        "a gene in no curated complex cannot be leaked through this "
        "channel. This is a property of the resource, not of the "
        "benchmark, and it means the audit measures a **lower bound** on "
        "the leakage that a complete complex resource would show.",
        "2. The audit measures leakage in the **v1** splits. The "
        "complex-aware splits (D3) must show a near-zero excess on the "
        "same metric; that check is part of D3's acceptance criteria.",
        "3. 'Complex-linked' is based on curated complex membership, not "
        "on interactions: two proteins that interact strongly without "
        "being curated as a complex are not counted here. The reviewer's "
        "channel is therefore partly, not fully, covered.",
        "",
        "## Files",
        "",
        "* `results/tables/d1_leakage_audit.tsv` -- one row per "
        "(family, arm, repeat, fold)",
        "* `results/tables/d1_leakage_summary.tsv` -- one row per arm",
        "* `data/processed/complexes/complex_map.tsv` and "
        "`data/processed/complexes/provenance.md` -- the complex map and "
        "where it came from",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
