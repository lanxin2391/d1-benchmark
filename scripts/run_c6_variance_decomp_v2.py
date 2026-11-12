"""C6 -- the variance decomposition with the fifth label, and the
before/after comparison the plan asks for.

The plan's C6 hypothesis was: "add an INDEPENDENT label and the label
main effect + interaction should take a LARGER share of the
between-arm spread, with the fold noise taking a larger share". The C2
circularity check showed that NCG is NOT independent of IntOGen
(Jaccard 0.522), so the hypothesis has to be re-read as: **what does
adding a label that is largely a restatement of an existing one do to
the decomposition?** That is the comparison this script writes.

It runs the corrected decomposition (`scripts/run_variance_decomp.py`)
on the 24-arm primary table and on the 30-arm extended table, and puts
the two side by side.

Outputs:
  results/tables/stats_variance_decomp_v24.tsv
  results/tables/stats_variance_components_v24.tsv
  results/tables/stats_variance_decomp_v2.tsv      (30 arms)
  results/tables/stats_variance_components_v2.tsv
  results/tables/stats_variance_decomp_comparison.tsv
  docs/variance_decomposition_v2.md
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/variance_decomposition_v2.md"
SCRIPT = ROOT / "scripts/run_variance_decomp.py"

GRID24 = ROOT / "results/runs_v40/a4_combined.tsv"
GRID30 = ROOT / "results/runs_v40/a4_combined_v2.tsv"


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def run_decomp(a4: Path, tag: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the decomposition and move its outputs to the tagged names."""
    log(f"running the decomposition on {a4.name} ({tag})")
    subprocess.run([sys.executable, str(SCRIPT), "--a4", str(a4)],
                   check=True, cwd=str(ROOT))
    comp = pd.read_csv(TABLES / "stats_variance_components.tsv", sep="\t")
    decomp = pd.read_csv(TABLES / "stats_variance_decomp.tsv", sep="\t")
    comp.to_csv(TABLES / f"stats_variance_components_{tag}.tsv", sep="\t",
                index=False, lineterminator="\n")
    decomp.to_csv(TABLES / f"stats_variance_decomp_{tag}.tsv", sep="\t",
                  index=False, lineterminator="\n")
    return decomp, comp


def main() -> int:
    for p in (GRID24, GRID30):
        if not p.exists():
            raise FileNotFoundError(p)

    decomp24, comp24 = run_decomp(GRID24, "v24")
    decomp30, comp30 = run_decomp(GRID30, "v2")

    rows = []
    for (name, a), (_, b) in zip(decomp24.iterrows(), decomp30.iterrows()):
        rows.append({
            "variance_source": a["variance_source"],
            "sigma2_24arms": a["sigma2"],
            "sigma2_30arms": b["sigma2"],
            "frac_between_24": a["fraction_of_between_arm_spread"],
            "frac_between_30": b["fraction_of_between_arm_spread"],
            "delta_frac_between": (b["fraction_of_between_arm_spread"]
                                    - a["fraction_of_between_arm_spread"]),
            "frac_total_24": a["fraction_of_total_variance"],
            "frac_total_30": b["fraction_of_total_variance"],
        })
    comp_tbl = pd.DataFrame(rows)
    comp_tbl.to_csv(TABLES / "stats_variance_decomp_comparison.tsv", sep="\t",
                    index=False, lineterminator="\n")
    log("wrote the comparison table")

    # restore the released (24-arm) tables as the canonical ones, since
    # the pre-registered endpoint is the 24-arm grid
    subprocess.run([sys.executable, str(SCRIPT), "--a4", str(GRID24)],
                   check=True, cwd=str(ROOT), stdout=subprocess.DEVNULL)
    log("restored the released 24-arm variance tables")

    write_doc(comp_tbl, comp24, comp30)
    log("wrote the doc")
    return 0


def write_doc(comp_tbl: pd.DataFrame, v24: pd.DataFrame,
              v30: pd.DataFrame) -> None:
    lines = [
        "# Variance decomposition with the fifth label (C6)",
        "",
        "**Date:** 2026-11-05  ",
        "**Script:** `scripts/run_c6_variance_decomp_v2.py`  ",
        "**Tables:** `stats_variance_decomp_comparison.tsv`, "
        "`stats_variance_decomp_v24.tsv`, `stats_variance_decomp_v2.tsv`",
        "",
        "## What is being compared",
        "",
        "The same corrected decomposition (`scripts/run_variance_decomp.py`) "
        "on the released 24-arm grid (6 networks x 4 labels) and on the "
        "30-arm grid (6 networks x 5 labels, NCG added).",
        "",
        "**The plan's hypothesis, re-read.** The plan expected that adding "
        "an *independent* label would let the label main effect and the "
        "interaction take a larger share of the between-arm spread. The "
        "C2 check showed NCG is not independent of IntOGen (Jaccard "
        "0.522), so the question is: what does adding a label that is "
        "largely a restatement of an existing one do to the "
        "decomposition?",
        "",
        "## Result",
        "",
        "| variance source | sigma2 (24) | sigma2 (30) | % between-arm "
        "(24) | % between-arm (30) | change |",
        "|---|---|---|---|---|---|",
    ]
    for r in comp_tbl.itertuples():
        lines.append(
            f"| {r.variance_source} | {r.sigma2_24arms:.4g} | "
            f"{r.sigma2_30arms:.4g} | {100*r.frac_between_24:.2f} % | "
            f"{100*r.frac_between_30:.2f} % | "
            f"{100*r.delta_frac_between:+.2f} pp |")
    lines += [
        "",
        "| grid | sigma2 network | sigma2 label | sigma2 interaction | "
        "sigma2 within-fold |",
        "|---|---|---|---|---|",
    ]
    for tag, v in (("24 arms", v24), ("30 arms", v30)):
        row = v.iloc[0]
        lines.append(f"| {tag} | {row.sigma2_network:.4g} | "
                     f"{row.sigma2_label:.4g} | {row.sigma2_interaction:.4g} "
                     f"| {row.sigma2_within_fold:.4g} |")
    lines += [
        "",
        "## Caveats",
        "",
        "1. The 30-arm grid mixes a label that is 65 % identical to "
        "IntOGen with four others; the decomposition therefore measures "
        "the effect of adding a *redundant* label, not of adding an "
        "independent one.",
        "2. The pre-registered endpoint stays the 24-arm grid; the "
        "30-arm decomposition is a pre-planned exploratory family (plan "
        "item C5, which pre-dates the analysis).",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
