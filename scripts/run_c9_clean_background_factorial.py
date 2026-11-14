"""C9 -- the primary factorial re-run with the negative set cleaned.

G8 measured how much the contaminated background matters, but it did
so POST HOC on the evaluation set of the existing runs. The stronger
question is: **does the headline survive if the benchmark is re-run with
the curated drivers removed from the universe in the first place?**

This script does exactly that, end to end:

1. **Clean the splits.** For every (network, label) split file of the
   40-repeat design, the rows of the genes that are in the NCG catalogue
   (all three tiers) **and are negatives of that label** are removed.
   The positives are untouched, so the comparison isolates one thing:
   the background. (G8's post-hoc version was slightly stronger -- it
   also dropped the NCG genes that happen to be positives of the label,
   from the evaluation set -- so the numbers here are the conservative
   version.)
2. **Re-run the 24 arms** on the cleaned splits, with the same locked
   learner, the same 5 folds x 40 repeats and the same three alphas,
   writing the per-gene score dumps again.
3. **Recompute the statistics**: the gene-cluster bootstrap CI, the
   interval variants (the five constructions) and the H1 verdict.
4. **Compare** the standard and the cleaned design, arm by arm.

Outputs:
  data/processed/splits_v40_cleanbg/            the cleaned splits
  results/runs_v40_cleanbg/                      the per-arm runs
  results/scores_v40_cleanbg/                   the per-gene dumps
  results/tables/stats_bootstrap_ci_gene_cleanbg.tsv
  results/tables/stats_ci_all_units_cleanbg.tsv
  results/tables/stats_*_cleanbg_v40.tsv        the statistics
  results/tables/c9_clean_background_comparison.tsv
  docs/c9_clean_background_factorial.md
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SPLITS = ROOT / "data/processed/splits_v40"
CLEAN_SPLITS = ROOT / "data/processed/splits_v40_cleanbg"
RUNS = ROOT / "results/runs_v40_cleanbg"
SCORES = ROOT / "results/scores_v40_cleanbg"
TABLES = ROOT / "results/tables"
NCG = ROOT / "data/processed/labels/ncg.tsv"
DOC = ROOT / "docs/c9_clean_background_factorial.md"

PY = r"E:\anaconda\python.exe"


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def run(args, label):
    log(label)
    r = subprocess.run([PY] + args, capture_output=True, text=True,
                       cwd=str(ROOT))
    if r.returncode != 0:
        print(r.stdout[-2000:])
        print(r.stderr[-2000:])
        raise SystemExit(f"{label} failed")
    tail = (r.stdout or r.stderr).strip().splitlines()
    for line in tail[-4:]:
        print("   ", line)
    return r.stdout


def main() -> int:
    ncg = pd.read_csv(NCG, sep="\t", dtype=str)
    ncg_all = set(ncg["symbol"])
    log(f"the NCG catalogue: {len(ncg_all)} genes")

    # ---- 1. the cleaned splits ------------------------------------
    CLEAN_SPLITS.mkdir(parents=True, exist_ok=True)
    n_files = 0
    removed_rows = 0
    removed_genes = set()
    for src in sorted(SPLITS.glob("*.tsv")):
        if src.name in ("N_per_universe.json",):
            continue
        if "__native_" not in src.name:
            continue
        d = pd.read_csv(src, sep="\t")
        if not {"gene", "y"} <= set(d.columns):
            continue
        is_ncg = d.gene.isin(ncg_all).to_numpy()
        is_neg = (d.y.to_numpy() == 0)
        drop = is_ncg & is_neg
        if not drop.any():
            out = d
        else:
            removed_rows += int(drop.sum())
            removed_genes |= set(d.gene.to_numpy()[drop])
            out = d[~drop]
        dst = CLEAN_SPLITS / src.name
        out.to_csv(dst, sep="\t", index=False)
        n_files += 1
    log(f"the cleaned splits: {n_files} files, {removed_rows} negative rows "
        f"removed ({len(removed_genes)} distinct genes)")

    # ---- 2. the 24 arms -------------------------------------------
    RUNS.mkdir(parents=True, exist_ok=True)
    SCORES.mkdir(parents=True, exist_ok=True)
    if (RUNS / "a4_combined.tsv").exists():
        log("the cleaned arms have already been run -- skipping")
    else:
        run(["scripts/run_a4.py",
             "--splits-dir", str(CLEAN_SPLITS),
             "--runs-dir", str(RUNS),
             "--scores-dir", str(SCORES),
             "--table-suffix", "_cleanbg_v40",
             "--n-rep", "40"],
            "re-running the 24 arms on the cleaned splits")

    # ---- 3. the statistics ---------------------------------------
    run(["scripts/run_gene_cluster_ci.py",
         "--scores-dir", str(SCORES),
         "--out", str(TABLES / "stats_bootstrap_ci_gene_cleanbg.tsv")],
        "the gene-cluster bootstrap CI")
    run(["scripts/run_ci_variants.py",
         "--combined", str(RUNS / "a4_combined.tsv"),
         "--out", str(TABLES / "stats_ci_all_units_cleanbg.tsv")],
        "the interval variants")
    run(["scripts/run_statistics.py",
         "--combined", str(RUNS / "a4_combined.tsv"),
         "--gene-ci", str(TABLES / "stats_bootstrap_ci_gene_cleanbg.tsv"),
         "--tables", str(TABLES),
         "--suffix", "_cleanbg_v40"], "the CI, the TOST and the H1 verdict")

    # ---- 4. the comparison ----------------------------------------
    std = pd.read_csv(TABLES / "a4_table1_delta_auroc_v40.tsv", sep="\t")
    cln = pd.read_csv(TABLES / "a4_table1_delta_auroc_cleanbg_v40.tsv",
                      sep="\t")
    rows = []
    std_idx = {(r["network"], r["label_set"]): r
               for r in std.to_dict("records")}
    for c in cln.to_dict("records"):
        key = (c["network"], c["label_set"])
        s = std_idx.get(key)
        if s is None:
            continue
        rows.append({
            "network": c["network"], "label": c["label_set"],
            "standard_a05": float(s["a0.5_mean"]),
            "cleaned_a05": float(c["a0.5_mean"]),
            "change": float(c["a0.5_mean"]) - float(s["a0.5_mean"]),
            "standard_sd": float(s["a0.5_std"]),
            "cleaned_sd": float(c["a0.5_std"]),
        })
    comp = pd.DataFrame(rows)
    comp.to_csv(TABLES / "c9_clean_background_comparison.tsv", sep="\t",
                index=False)

    # the H1 verdicts, standard and cleaned
    v_std = pd.read_csv(TABLES / "stats_h1_verdict_v40.tsv", sep="\t")
    v_cln = pd.read_csv(TABLES / "stats_h1_verdict_cleanbg_v40.tsv", sep="\t")

    # the five interval constructions, both designs
    units_std = pd.read_csv(TABLES / "stats_ci_all_units_v40.tsv", sep="\t")
    units_cln = pd.read_csv(TABLES / "stats_ci_all_units_cleanbg.tsv",
                            sep="\t")

    write_doc(comp, v_std, v_cln, units_std, units_cln, removed_genes,
              removed_rows)
    log("wrote the comparison table and the doc")
    return 0


def write_doc(comp, v_std, v_cln, units_std, units_cln, removed_genes,
              removed_rows):
    lines = [
        "# C9 -- the factorial re-run with a cleaned negative set",
        "",
        "**Date:** 2026-11-12  ",
        "**Script:** `scripts/run_c9_clean_background_factorial.py`  ",
        "**Table:** `results/tables/c9_clean_background_comparison.tsv`",
        "",
        "## What was done",
        "",
        "G8 measured the contamination of our negative set post hoc, on the "
        "evaluation set of the existing runs. This is the stronger version: "
        "**the benchmark is re-run from the splits with the curated drivers "
        "removed from the universe.** The NCG genes (all three tiers) that "
        "are negatives of a label leave the split files; the positives are "
        "untouched, so the comparison isolates the background. The learner, "
        "the folds, the repeats and the alphas are unchanged.",
        "",
        f"- negative rows removed: **{removed_rows:,}** "
        f"({len(removed_genes)} distinct genes)",
        "",
        "## The margin, arm by arm",
        "",
        "| network | label | standard | cleaned | change |",
        "|---|---|---|---|---|",
    ]
    for r in comp.sort_values("standard_a05", ascending=False).itertuples():
        lines.append(f"| {r.network} | {r.label} | {r.standard_a05:+.4f} | "
                     f"{r.cleaned_a05:+.4f} | {r.change:+.4f} |")
    pos_std = int((comp.standard_a05 > 0).sum())
    pos_cln = int((comp.cleaned_a05 > 0).sum())
    flips = comp[((comp.standard_a05 > 0) & (comp.cleaned_a05 < 0))
                 | ((comp.standard_a05 < 0) & (comp.cleaned_a05 > 0))]
    lines += [
        "",
        "## The verdict",
        "",
        f"- positive arms: **{pos_std}/24** before, **{pos_cln}/24** after; "
        f"sign flips: **{len(flips)}**",
        f"- the median margin: **{comp.standard_a05.median():+.4f}** before, "
        f"**{comp.cleaned_a05.median():+.4f}** after",
        "",
        "| interval | arms excluding 0 (standard) | arms excluding 0 "
        "(cleaned) |",
        "|---|---|---|",
    ]
    for u in ("fold_pct", "repeat_pct", "repeat_t", "repeat_bca", "gene"):
        s = units_std[units_std.interval == u] if "interval" in units_std \
            else None
        c = units_cln[units_cln.interval == u] if "interval" in units_cln \
            else None
        if s is None or c is None:
            continue
        lines.append(f"| {u} | {int(s.n_excl0.iloc[0])}/{len(s)} | "
                     f"{int(c.n_excl0.iloc[0])}/{len(c)} |")
    lines += [
        "",
        "## Reading",
        "",
        "This is the sensitivity analysis the negative-set audit demands: if "
        "the headline survives with the curated drivers removed from the "
        "background, then the contamination of the negative set is not what "
        "carries the result.",
        "",
        "## Caveats",
        "",
        "1. The exploration family (amendment 4): the pre-registered "
        "endpoint is the 24-arm standard design, untouched.",
        "2. Removing the NCG genes from the universe changes the walk's "
        "normalisation and the fold sizes slightly, so the cleaned margins "
        "are not the standard margins with a few negatives deleted; they "
        "are a separate, equally valid run.",
        "3. The cleaning uses one catalogue (NCG). A driver curated only in "
        "another source (for example the Reactome-curated positives of the "
        "open-targets labels) is not removed.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
