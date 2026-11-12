# Variance decomposition with the fifth label (C6)

**Date:** 2026-11-05  
**Script:** `scripts/run_c6_variance_decomp_v2.py`  
**Tables:** `stats_variance_decomp_comparison.tsv`, `stats_variance_decomp_v24.tsv`, `stats_variance_decomp_v2.tsv`

## What is being compared

The same corrected decomposition (`scripts/run_variance_decomp.py`) on the released 24-arm grid (6 networks x 4 labels) and on the 30-arm grid (6 networks x 5 labels, NCG added).

**The plan's hypothesis, re-read.** The plan expected that adding an *independent* label would let the label main effect and the interaction take a larger share of the between-arm spread. The C2 check showed NCG is not independent of IntOGen (Jaccard 0.522), so the question is: what does adding a label that is largely a restatement of an existing one do to the decomposition?

## Result

| variance source | sigma2 (24) | sigma2 (30) | % between-arm (24) | % between-arm (30) | change |
|---|---|---|---|---|---|
| network (between) | 0.0001918 | 0.0001788 | 24.82 % | 26.95 % | +2.13 pp |
| label (between) | 0.0002687 | 0.0002041 | 34.78 % | 30.75 % | -4.02 pp |
| network x label (interaction) | 0.0003095 | 0.0002784 | 40.07 % | 41.95 % | +1.88 pp |
| fold (within cell) | 2.609e-06 | 2.351e-06 | 0.34 % | 0.35 % | +0.02 pp |

| grid | sigma2 network | sigma2 label | sigma2 interaction | sigma2 within-fold |
|---|---|---|---|---|
| 24 arms | 0.0001918 | 0.0002687 | 0.0003095 | 0.0005218 |
| 30 arms | 0.0001788 | 0.0002041 | 0.0002784 | 0.0004702 |

## Caveats

1. The 30-arm grid mixes a label that is 65 % identical to IntOGen with four others; the decomposition therefore measures the effect of adding a *redundant* label, not of adding an independent one.
2. The pre-registered endpoint stays the 24-arm grid; the 30-arm decomposition is a pre-planned exploratory family (plan item C5, which pre-dates the analysis).
