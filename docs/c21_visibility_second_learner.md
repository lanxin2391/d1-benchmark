# Control 21 -- the visibility decomposition under the second learner

**Date:** 2026-11-16  
**Script:** `scripts/run_c21_visibility_second_learner.py`  
**Tables:** `c21_visibility_second_learner.tsv`, `c21_visibility_learner_comparison.tsv`

## The question

Control 17 measured the MARGIN under the proximity learner and found 16/24 arms positive (median +0.011 against +0.0809), with an arm-level Spearman of +0.30. That is consistent with Guala & Sonnhammer (2017), and it is not a new claim.

But the project's claim is not the margin: it is the DECOMPOSITION of the reported advantage into topology, node degree and study visibility. **No visibility measurement existed under a second learner.** This control computes the bait-matched null with the proximity learner as the scorer -- the same splits, the same folds, the same publication counts, the same matching rule (the 20 negatives nearest in log10(n_pubmed + 1)) -- and compares it with the published RWR numbers.

The matcher in this file was checked first on the RWR score dumps of one arm: it reproduces the published G3 statistic for that arm to within 0.02 before any proximity result is used.

## The result

| learner | median arm: above the unmatched null | above the matched null |
|---|---|---|
| RWR (published G3) | 17.6% | 11.6% |
| proximity (this control) | 7.7% | 4.8% |

Per label (median over the six networks):

| label | RWR unmatched -> matched | proximity unmatched -> matched |
|---|---|---|
| intogen2024 | 20.7% -> 10.9% | 8.4% -> 5.6% |
| ot_all | 31.0% -> 20.2% | 8.3% -> 5.0% |
| ot_nolit | 27.1% -> 18.9% | 7.2% -> 4.7% |
| intogen_temporal_new | 6.8% -> 4.7% | 7.5% -> 4.8% |

## Reading

- the drop caused by the publication matching, RWR: median 5.9%; proximity: median 2.3%
- the arm-level correlation between the two learners' matched-null statistics: Spearman +0.283

## Caveats

1. The proximity learner's score is the negative distance to the nearest seed, so 'above' means 'closer to a seed'. Because that score is an INTEGER distance, almost every test gene ties with another one, and G3's strict comparison would count a tie as 'not above' -- which inflates the statistic (91 % of the positives 'beat' the null on intact x intogen2024, while the AUROC of the same scores is 0.627). The comparison here is therefore made on the protocol's RANKS (decision D-07: score desc, then gene symbol asc), which resolves the ties the documented way. RWR has no ties, so its numbers are unchanged -- the self-check on the RWR dumps reproduces the published G3 values exactly.
2. **Exploratory** (amendment 4): it does not enter the H1 decision.
3. This control measures the VISIBILITY component under a second learner. It does not re-run the full factorial (the rewiring null under the proximity learner would need ~470,000 multi-source BFS runs and is out of scope).
