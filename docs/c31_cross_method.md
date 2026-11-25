# Control 31 — the cross-method complementarity

**Date:** 2026-11-25  
**Script:** `scripts/run_c31_cross_method.py`  
**Tables:** `c31_cross_method_overlap.tsv`, `c31_cross_method_union.tsv`

## The question

Control 19 compared the methods by their accuracy against degree. This control compares them by their **candidates**: what does a user get by running a second method on top of the first? Six scorers on the same folds, the same test genes and the same top 100:

* the RWR prior (the headline learner),
* degree (the baseline),
* the betweenness centrality (a path-based topology score),
* the 1-hop seed count and the 2-hop seed count,
* the proximity learner (the distance to the nearest seed).

## The Jaccard overlap of the top-100 picks (the mean over the folds)

| | vs the RWR prior | vs degree | vs betweenness | vs the 1-hop seed count | vs the 2-hop seed count | vs the proximity learner |
|---|---|---|---|---|---|---|
| **the RWR prior** | - | 0.277 | 0.215 | 0.462 | 0.299 | 0.103 |
| **degree** | 0.277 | - | 0.235 | 0.392 | 0.401 | 0.070 |
| **betweenness** | 0.215 | 0.235 | - | 0.193 | 0.242 | 0.059 |
| **the 1-hop seed count** | 0.462 | 0.392 | 0.193 | - | 0.445 | 0.141 |
| **the 2-hop seed count** | 0.299 | 0.401 | 0.242 | 0.445 | - | 0.078 |
| **the proximity learner** | 0.103 | 0.070 | 0.059 | 0.141 | 0.078 | - |

## How much is shared and how much is unique

| scorer | the picks | the genes every method picks | the genes only this method picks |
|---|---|---|---|
| the RWR prior | 100 | 3.1 | 24.0 |
| degree | 100 | 3.1 | 24.6 |
| betweenness | 100 | 3.1 | 43.7 |
| the 1-hop seed count | 100 | 3.1 | 9.2 |
| the 2-hop seed count | 100 | 3.1 | 18.3 |
| the proximity learner | 100 | 3.1 | 65.7 |

## Reading

## Caveats

1. 20 folds per arm (one per repeat), the same folds for every scorer.
2. The closeness centrality is not in this table: its landmark estimate costs 256 breadth-first searches per fold, and Control 19 already reports it. Its accuracy against degree (+0.0312 in 24 of 24 arms) is in that section.
3. The overlap is a property of the top-K sets, not of the scores; two scorers can agree on the picks and disagree on the order.
4. **Exploratory** (amendment 4): no pre-registered endpoint is touched.
