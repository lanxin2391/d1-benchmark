# Direction 1 -- is the consensus prior robust?

**Date:** 2026-11-05  
**Script:** `scripts/run_consensus_robustness.py`  
**Table:** `results/tables/consensus_robustness.tsv`

## What was varied

* **the aggregation rule**: mean of percentile ranks (the original), median of ranks, and reciprocal rank fusion (k = 60). A z-score aggregation was tried and **dropped**: the six networks give the same percentile rank for most genes (cross-network SD is exactly 0 for the majority), so z-scoring collapses every gene they agree on into a tie and keeps only the disagreement -- it measures disagreement, not consensus, and its AUROC is negative
* **the weights**: the unweighted consensus vs a reliability-weighted one (weight = 1 / the network's arm-level repeat-level SD, measured in S3) vs the reliability-weighted consensus with Reactome left out

Everything else is as in run_consensus_prior.py: the shared universe, the D-12 splits, ONE common degree baseline for every prior, and the S1 rank-SD stability metric.

## Result

| label | mean_rank | median_rank | reliability_weighted | reliability_weighted_without_reactome | rrf |
|---|---|---|---|---|---|
| intogen2024 | +0.1128 | +0.1092 | +0.1213 | +0.1201 | +0.1128 |
| intogen_temporal_new | +0.0358 | +0.0321 | +0.0388 | +0.0512 | +0.0356 |
| ot_all | +0.1421 | +0.1533 | +0.1541 | +0.1412 | +0.1423 |
| ot_nolit | +0.1384 | +0.1528 | +0.1527 | +0.1416 | +0.1387 |

| label | mean_rank (rank SD) | median_rank (rank SD) | reliability_weighted (rank SD) | reliability_weighted_without_reactome (rank SD) | rrf (rank SD) |
|---|---|---|---|---|---|
| intogen2024 | 0.0409 | 0.0465 | 0.0397 | 0.0398 | 0.0410 |
| intogen_temporal_new | 0.0528 | 0.0560 | 0.0500 | 0.0497 | 0.0528 |
| ot_all | 0.0364 | 0.0409 | 0.0355 | 0.0363 | 0.0363 |
| ot_nolit | 0.0379 | 0.0421 | 0.0366 | 0.0371 | 0.0380 |

## Reading

* **The aggregation rule does not carry the result.** Across the four rules the median delta-AUROC moves by at most 0.0191 (label intogen_temporal_new), and the temporal label -- the one the consensus was built to repair -- stays positive under every rule.
* **The weights help, as predicted.** The reliability-weighted consensus is the best variant on the temporal label, because Reactome (the least reliable network for newly discovered drivers) gets the smallest weight; dropping it entirely is better still, which is the same finding as the leave-one-out table in docs/consensus_prior.md.
* **The stability cost is unchanged.** The consensus is less stable per gene than its members under every rule, so the bias-variance trade-off reported in docs/consensus_prior.md is a property of the aggregation, not of the particular rule.

## Caveats

1. The weights are measured on the same data (the arm-level repeat SD), so the weighted consensus is not an out-of-sample improvement; a leave-one-label-out check would be needed for that, and is not done here.
2. Only the shared universe is available for a consensus, so the absolute numbers are not comparable with the released 24-cell table (see docs/consensus_prior.md, caveat 1).
