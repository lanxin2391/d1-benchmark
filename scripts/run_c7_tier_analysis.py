"""C7, part 2 -- what the label quality does to the margin.

Part 1 (`run_c7_build_tiers.py`) gave every positive gene a tier. This
script asks what that tier is worth, in three steps.

1. **The tier-restricted margin.** Using the per-gene score dumps of the
   40-repeat design, the margin is recomputed with only one tier's
   positives counted as positives (the other tiers' positives become
   background). The seeds -- and therefore the prior itself -- do not
   change, so the comparison isolates the *evaluation* noise: if the
   low tier carries the margin, the margin is a label-noise artefact.

2. **The noise anchor.** ClinGen's expert panels disagree on 17 % of
   multi-laboratory variants (the figure the literature review nominates
   as the external anchor for label noise). That rate is simulated
   directly: in each fold a random 17 % of the positives is demoted to
   background, and the margin is recomputed. The theoretical ceiling of
   a perfect scorer under the same misassignment is 1 - f/2 = 0.915.

3. **The interval width.** The plan's B3 question: is a wide interval
   explained by a low-quality label and a small positive set?

Outputs:
  results/tables/label_tier_delta_auroc.tsv
  results/tables/clingen_noise_anchor.tsv
  results/tables/label_quality_interval_width.tsv
  docs/label_quality.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_pauc_top10 import auroc  # noqa: E402

SCORES = ROOT / "results/scores_v40"
SCORES_CLINGEN = ROOT / "results/scores_v40_clingen"
TIERS = ROOT / "results/tables/label_quality_tiers.tsv"
TABLES = ROOT / "results/tables"
CI = TABLES / "stats_bootstrap_ci.tsv"
DOC = ROOT / "docs/label_quality.md"
ALPHA = 0.5
TIER_ORDER = ["high", "mid", "low"]
NOISE_RATE = 0.17          # the ClinGen inter-laboratory disagreement rate
SEED = 20261105


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_dump(path: Path):
    d = np.load(path, allow_pickle=True)
    alphas = np.asarray(d["alphas"], dtype=float)
    ai = int(np.argmin(np.abs(alphas - ALPHA)))
    return (np.asarray(d["genes"]).astype(str),
            np.asarray(d["y"]).astype(int),
            np.asarray(d["degree"], dtype=float),
            np.asarray(d["rwr"][ai], dtype=float),
            np.asarray(d["in_test"], dtype=bool))


def tier_margins(path: Path, tier_of: dict, rng=None,
                 demote_frac: float = 0.0):
    """The margin per tier, and (optionally) with a random fraction of
    the positives demoted to background (the noise simulation)."""
    genes, y, degree, rwr, in_test = load_dump(path)
    is_pos = y == 1
    tiers = np.array([tier_of.get(g, "untiered") for g in genes])
    out = {}
    for tier in TIER_ORDER:
        sel = is_pos & (tiers == tier)
        if sel.sum() == 0:
            continue
        ds = []
        for j in range(rwr.shape[1]):
            idx = np.where(in_test[j])[0]
            if len(idx) == 0:
                continue
            pos = idx[sel[idx]]
            if len(pos) < 3:
                continue
            if demote_frac > 0:
                k = int(round(demote_frac * len(pos)))
                if k > 0:
                    drop = rng.choice(pos, size=k, replace=False)
                    pos = np.setdiff1d(pos, drop)
                    if len(pos) < 3:
                        continue
            rest = np.setdiff1d(idx, pos, assume_unique=True)
            yy = np.concatenate([np.ones(len(pos), dtype=int),
                                 np.zeros(len(rest), dtype=int)])
            ds.append(auroc(np.concatenate([rwr[pos, j], rwr[rest, j]]), yy)
                       - auroc(np.concatenate([degree[pos], degree[rest]]), yy))
        if ds:
            out[tier] = (float(np.mean(ds)), float(np.std(ds, ddof=1)),
                         len(ds), int(sel.sum()))
    return out, int(is_pos.sum())


def main() -> int:
    tiers = pd.read_csv(TIERS, sep="\t")
    tier_maps = {lab: dict(zip(g.gene, g.tier))
                 for lab, g in tiers.groupby("label")}

    # ---- 1. the tier-restricted margins -------------------------------
    rows = []
    arms = []
    for p in sorted(SCORES.glob("*.npz")):
        net, lab = p.stem.split("__")
        arms.append((net, lab, p))
    for p in sorted(SCORES_CLINGEN.glob("*.npz")):
        net, lab = p.stem.split("__")
        arms.append((net, lab, p))
    for net, lab, p in arms:
        tm = tier_maps.get(lab)
        if not tm:
            log(f"  {net} x {lab}: no tier table")
            continue
        res, n_pos = tier_margins(p, tm)
        full = res.get("high")
        for tier, (mean, sd, nfold, nsel) in res.items():
            rows.append({"network": net, "label": lab, "tier": tier,
                         "n_positives_in_label": nsel,
                         "n_folds": nfold,
                         "delta_auroc": mean, "delta_auroc_sd": sd})
        # the "all positives" margin from the same dump, for reference
        allres, _ = tier_margins(p, {g: "high" for g in
                                      pd.read_csv(TIERS, sep="\t")
                                      .query("label == @lab").gene})
        log(f"  {net} x {lab}: " + ", ".join(
            f"{t}={res[t][0]:+.4f}" for t in TIER_ORDER if t in res))
    tier_tbl = pd.DataFrame(rows)
    tier_tbl.to_csv(TABLES / "label_tier_delta_auroc.tsv", sep="\t",
                    index=False, lineterminator="\n")
    log(f"wrote label_tier_delta_auroc.tsv ({len(tier_tbl)} rows)")

    # ---- 2. the ClinGen noise anchor ---------------------------------
    rng = np.random.default_rng(SEED)
    anchor_rows = []
    for p in sorted(SCORES_CLINGEN.glob("*.npz")):
        net, lab = p.stem.split("__")
        tm = tier_maps.get(lab, {})
        clean, n_pos = tier_margins(p, {g: "high" for g in tm})
        noisy, _ = tier_margins(p, {g: "high" for g in tm}, rng=rng,
                                demote_frac=NOISE_RATE)
        if not clean or not noisy:
            continue
        c = clean["high"][0]
        n = noisy["high"][0]
        anchor_rows.append({
            "network": net, "label": lab,
            "n_positives": n_pos,
            "delta_auroc_clean": c,
            "delta_auroc_17pct_misassigned": n,
            "change": n - c,
            "relative_change": (n - c) / c if c else np.nan,
            "theoretical_auroc_ceiling": 1 - NOISE_RATE / 2,
        })
        log(f"  {net}: clean {c:+.4f} -> 17 % misassigned {n:+.4f} "
            f"({n - c:+.4f})")
    anchor = pd.DataFrame(anchor_rows)
    anchor.to_csv(TABLES / "clingen_noise_anchor.tsv", sep="\t", index=False,
                  lineterminator="\n")
    log(f"wrote clingen_noise_anchor.tsv ({len(anchor)} rows)")

    # ---- 3. the interval width ----------------------------------------
    ci = pd.read_csv(CI, sep="\t")
    ci = ci[ci.alpha == ALPHA].copy()
    ci["ci_width_primary"] = ci.ci_high - ci.ci_low
    ci["ci_width_gene"] = ci.ci_high_gene - ci.ci_low_gene
    # the positive-set size and the high-tier share per arm
    sizes = {}
    for net, lab, p in arms:
        genes, y, _d, _r, _t = load_dump(p)
        tm = tier_maps.get(lab, {})
        is_pos = y == 1
        tiers_arr = np.array([tm.get(g, "untiered") for g in genes])
        sizes[(net, lab)] = (int(is_pos.sum()),
                             float((tiers_arr[is_pos] == "high").mean()))
    ci["n_positives"] = [sizes.get((r.network, r.label), (np.nan, np.nan))[0]
                         for r in ci.itertuples()]
    ci["frac_high_tier"] = [sizes.get((r.network, r.label), (np.nan, np.nan))[1]
                            for r in ci.itertuples()]
    ci.to_csv(TABLES / "label_quality_interval_width.tsv", sep="\t",
              index=False, lineterminator="\n")
    ok = ci.dropna(subset=["n_positives"])
    rho_size = stats.spearmanr(ok.n_positives, ok.ci_width_gene)
    rho_qual = stats.spearmanr(ok.frac_high_tier, ok.ci_width_gene)
    log(f"the gene-cluster CI width vs the positive-set size: "
        f"rho = {rho_size.statistic:+.3f} (p = {rho_size.pvalue:.3g})")
    log(f"the gene-cluster CI width vs the high-tier share: "
        f"rho = {rho_qual.statistic:+.3f} (p = {rho_qual.pvalue:.3g})")

    write_doc(tier_tbl, anchor, ci, rho_size, rho_qual)
    log("wrote the doc")
    return 0


def write_doc(tier_tbl, anchor, ci, rho_size, rho_qual) -> None:
    piv = tier_tbl.pivot_table(index=["network", "label"], columns="tier",
                               values="delta_auroc")
    lines = [
        "# C7 -- label quality: which positives carry the margin?",
        "",
        "**Date:** 2026-11-05  ",
        "**Scripts:** `scripts/run_c7_build_tiers.py`, "
        "`scripts/run_c7_tier_analysis.py`  ",
        "**Tables:** `label_quality_tiers.tsv`, "
        "`label_tier_delta_auroc.tsv`, `clingen_noise_anchor.tsv`, "
        "`label_quality_interval_width.tsv`",
        "",
        "## The tiers",
        "",
        "Every positive gene gets a tier from its own catalogue's "
        "evidence: ClinGen's review status (Definitive/Strong = high), "
        "IntOGen's number of distinct detecting methods (>= 4 = high), "
        "the NCG literature count (>= 10 = high), and for Open Targets "
        "whether the cancer associations are carried by a non-literature "
        "source.",
        "",
        "| label | high | mid | low |",
        "|---|---|---|---|",
    ]
    t = pd.read_csv(TIERS, sep="\t")
    for lab, g in t.groupby("label"):
        vc = g.tier.value_counts()
        lines.append(f"| {lab} | {vc.get('high', 0)} | {vc.get('mid', 0)} | "
                     f"{vc.get('low', 0)} |")
    lines += [
        "",
        "## 1. The tier-restricted margin",
        "",
        "The margin is recomputed with only one tier's positives counted "
        "as positives; the other tiers' positives fall into the "
        "background. The seeds, and therefore the prior, do not change, "
        "so this isolates the evaluation noise.",
        "",
        "| network | label | high | mid | low |",
        "|---|---|---|---|---|",
    ]
    for (net, lab), r in piv.iterrows():
        cells = []
        for tier in TIER_ORDER:
            v = r.get(tier)
            cells.append(f"{v:+.4f}" if pd.notna(v) else "-")
        lines.append(f"| {net} | {lab} | " + " | ".join(cells) + " |")
    lines += [
        "",
        "## 2. The noise anchor (ClinGen's 17 % disagreement rate)",
        "",
        "ClinGen's expert panels disagree on 17 % of multi-laboratory "
        "variants. Simulating exactly that rate -- a random 17 % of the "
        "positives demoted to background in every fold -- gives:",
        "",
        "| network | n positives | clean | 17 % misassigned | change |",
        "|---|---|---|---|---|",
    ]
    for r in anchor.itertuples():
        lines.append(f"| {r.network} | {r.n_positives} | "
                     f"{r.delta_auroc_clean:+.4f} | "
                     f"{r.delta_auroc_17pct_misassigned:+.4f} | "
                     f"{r.change:+.4f} |")
    lines += [
        "",
        "The theoretical ceiling of a *perfect* scorer under the same "
        "misassignment is 1 - f/2 = 0.915, so a 17 % misassignment rate "
        "costs at most 8.5 points of AUROC however well the prior ranks.",
        "",
        "## 3. Why some intervals are wide",
        "",
        "The gene-cluster bootstrap interval width against the arm's "
        f"positive-set size: Spearman rho = {rho_size.statistic:+.3f} "
        f"(p = {rho_size.pvalue:.3g}); against the share of high-tier "
        f"positives: rho = {rho_qual.statistic:+.3f} "
        f"(p = {rho_qual.pvalue:.3g}).",
        "",
        "## Caveats",
        "",
        "1. The tier-restricted margin moves the other tiers' positives "
        "into the background, so it is not a label's own margin; it is a "
        "diagnostic of where the margin comes from.",
        "2. The Open Targets tier uses the literature/non-literature "
        "split because almost every (target, cancer disease) pair in the "
        "datasource table carries exactly one source (88,510 of 89,111 in "
        "the first part), so a source-count tier would be constant.",
        "3. **Exploratory** (amendment 4); the H1 decision is untouched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
