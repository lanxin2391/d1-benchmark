# Control 22 (feasibility) -- the rewiring null under the second learner

**Date:** 2026-11-16  
**Script:** `scripts/run_c22_rewiring_second_learner.py`  
**Tables:** `c22_rewiring_second_learner.tsv`, `c22_rewiring_second_learner_summary.tsv`

## What this is

Control 21 settled the visibility component under the second learner. The remaining learner-sensitive piece is the **rewiring null** (Control 2 / W5): is the signal carried by the WIRING, or only by the degree sequence? Under RWR the answer is 23/24 arms with z > 5.

The full version under the proximity learner would need about 470,000 multi-source BFS runs, so this is a **FEASIBILITY CHECK on one arm**: funmap x intogen2024, using the 20 rewired networks that Control 2 had already produced (the rewiring is not repeated -- it is read from data/processed/rewired/funmap/).

## The result

| quantity | RWR (published, 100 rewirings) | proximity (this feasibility check, 20 rewirings) |
|---|---|---|
| the observed margin | +0.0998 | +0.0804 |
| the rewired null | +0.0024 +- 0.0103 | -0.0124 +- 0.0072 |
| z | +94.16 | +12.92 |
| the rewirings | 100 | 20 |

## Two checks that the control is sound

1. **The rewiring preserves the degree sequence exactly**: over the 20 rewirings the largest difference in any gene's degree is **0**. The null therefore keeps the hubs and only re-pairs them.
2. **The control is computable**: the whole check -- 20 rewirings x 200 folds of multi-source BFS -- runs in a few minutes, so the full version (6 networks x 4 labels x ~100 rewirings) is a matter of compute, not of method.

## How to read it, and what it is not

- On this arm the observed proximity margin (+0.0804) sits above the rewired null (-0.0124 +- 0.0072), z = +12.92: the same qualitative answer as under RWR (z = +94.16), i.e. the wiring carries the signal rather than the degree sequence alone.
- **This is one arm out of 24, with 20 rewirings instead of ~100.** It is a feasibility check, not the full control, and it is reported as such in the manuscript, the deck, the speaker script and the Q&A bank. The full version is the first item of the validation plan.
- **Exploratory** (amendment 4): it does not enter the H1 decision.
