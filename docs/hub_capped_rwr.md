# E3 -- hub-capped RWR: is the prior's advantage a hub effect?

**Date:** 2026-11-05  
**Script:** `scripts/run_hub_capped_rwr.py`  
**Table:** `results/tables/hub_capped_rwr.tsv`

## What was run

Three variants on the SAME v1 splits as the primary factorial, with the degree baseline untouched:

* `trunc128` / `trunc256` -- the nodes with degree > K are removed from the PROPAGATION (they stay in the universe and score 0); `*_exclude_hubs` additionally drops them from the evaluation set, because a tie block at the bottom changes the AUROC denominator
* `degpen0.5` -- the degree-penalised transition `W[i,j] = A[i,j] / (deg(i)^0.5 * deg(j)^0.5)`

## Result

| variant | median delta-AUROC | change vs baseline (median) | arms within 5 % of the baseline |
|---|---|---|---|
| baseline | +0.0809 | +0.0000 | 24/24 |
| trunc128 | +0.0252 | -0.0499 | 0/24 |
| trunc128_exclude_hubs | +0.0849 | -0.0001 | 15/24 |
| trunc256 | +0.0648 | -0.0054 | 10/24 |
| trunc256_exclude_hubs | +0.0853 | +0.0007 | 15/24 |
| degpen0.5 | +0.0979 | +0.0063 | 7/24 |

## Per arm

| network | label | baseline | trunc128 | trunc128 excl-hubs | degpen0.5 |
|---|---|---|---|---|---|
| funmap | intogen2024 | +0.1001 | +0.0388 | +0.1196 | +0.1323 |
| funmap | intogen_temporal_new | +0.0196 | -0.0360 | +0.0150 | +0.0338 |
| funmap | ot_all | +0.0853 | +0.0183 | +0.0894 | +0.1050 |
| funmap | ot_nolit | +0.0780 | +0.0190 | +0.0744 | +0.0950 |
| intact | intogen2024 | +0.0651 | -0.1928 | +0.0404 | +0.0714 |
| intact | intogen_temporal_new | +0.0071 | -0.1190 | -0.0253 | -0.0088 |
| intact | ot_all | +0.0842 | -0.1816 | +0.0802 | +0.1009 |
| intact | ot_nolit | +0.0784 | -0.1568 | +0.0882 | +0.0944 |
| reactome | intogen2024 | +0.0429 | +0.0332 | +0.0411 | +0.0365 |
| reactome | intogen_temporal_new | -0.0125 | -0.0113 | -0.0126 | -0.0177 |
| reactome | ot_all | +0.1103 | +0.0950 | +0.1099 | +0.1143 |
| reactome | ot_nolit | +0.0996 | +0.0878 | +0.1001 | +0.1048 |
| rna_coexp | intogen2024 | +0.0521 | -0.0695 | +0.0719 | +0.0774 |
| rna_coexp | intogen_temporal_new | +0.0127 | +0.0037 | +0.0323 | +0.0245 |
| rna_coexp | ot_all | +0.0767 | +0.0337 | +0.0843 | +0.0946 |
| rna_coexp | ot_nolit | +0.0834 | +0.0392 | +0.0855 | +0.1067 |
| string_full700 | intogen2024 | +0.1133 | -0.0099 | +0.1170 | +0.1238 |
| string_full700 | intogen_temporal_new | +0.0604 | +0.0519 | +0.0657 | +0.0601 |
| string_full700 | ot_all | +0.1201 | -0.0106 | +0.1200 | +0.1264 |
| string_full700 | ot_nolit | +0.1267 | +0.0314 | +0.1258 | +0.1329 |
| string_phys700 | intogen2024 | +0.1189 | +0.1063 | +0.1198 | +0.1241 |
| string_phys700 | intogen_temporal_new | +0.0378 | +0.0434 | +0.0377 | +0.0385 |
| string_phys700 | ot_all | +0.1582 | +0.1391 | +0.1542 | +0.1625 |
| string_phys700 | ot_nolit | +0.1711 | +0.1550 | +0.1670 | +0.1752 |

## Reading

**The two truncation conventions answer two different questions, and both must be reported.** The trigger for the hub-debiased model of the second paper (work plan H6) is: *if the RWR advantage DISAPPEARS once the hubs are removed, the orthogonalised model becomes necessary*. Under the clean convention (hubs excluded from the evaluation) it does not disappear; under the harsh one (hubs silenced, i.e. kept in the evaluation with score 0) it shrinks to a third of its size and flips sign in a third of the arms. The trigger is therefore **not met, but not trivially so**:

* **The harshest convention** (the hubs are silenced -- they stay in the universe and score 0, so they tie at the bottom) leaves a median delta-AUROC of +0.0252 against the baseline's +0.0809: the advantage shrinks to 31% of its size but stays clearly positive, and only 8/24 arms change sign.
* **The clean convention** (the hubs are removed from the evaluation set as well, so the AUROC denominator is unchanged) leaves the effect essentially untouched (+0.0849 vs +0.0809, median change -0.0001, 15/24 arms within 5 %).
* **The degree-penalised walk** (`W[i,j] = A[i,j]/(deg(i)^0.5 deg(j)^0.5)`, which down-weights the hubs WITHOUT removing them) *increases* the effect (+0.0979, median change +0.0063).

So the benchmark's headline claim -- that the network prior contributes beyond node degree -- is **not a hub artefact**. That is the same conclusion the GraphSAGE arm had to reach the opposite way round: it had to down-sample edges to degree <= 128 to survive its own degree distribution, while the locked RWR learner survives hub removal and even improves under degree penalisation.

**What this says about the paper's claim, precisely.**

* The claim *the prior contributes beyond node degree* survives hub removal: on the part of the network without the hubs the prior still beats the degree baseline by the same amount (+0.0849 vs +0.0809). The headline is **not** a hub artefact.
* But a large part of the raw margin IS hub-mediated: silencing the hubs instead of removing them takes the median effect to +0.0252 (31% of its size) and flips the sign in 8/24 arms. The two numbers differ because silencing pushes the hubs to the bottom of the ranking and the hubs carry a large share of the positives in the Open Targets labels.
* Therefore the paper must state both numbers and quote the smaller one. The honest sentence is: *the prior's advantage is not entirely hub-mediated, but roughly two thirds of the raw margin sits in the hubs' own ranking position.*

**Consequence for the second paper.** The orthogonalisation idea that H6 proposes is, in effect, what `degpen0.5` already tests -- a walk that down-weights high-degree nodes. It **increases** the effect (median +0.0979 vs +0.0809, median change +0.0063), so there is no hub-driven failure left for it to fix. If the second paper wants a method, the degree-penalised walk is the candidate that this project has already measured; H6 as a *corrective* model is not needed.

A useful by-product: the comparison quantifies how much of the raw margin is hub-mediated -- roughly two thirds of it -- which is a number the field does not have, and it is the same quantity the S1 stability atlas measures from the other side (the top 1 % of the ranking contains almost no drivers).

## Caveats

1. The truncation variants change the evaluation set (the hubs tie at the bottom, or are dropped), so the comparison with the baseline is not a pure propagation change; both conventions are reported for that reason.
2. `degpen0.5` changes the walk, not the graph: it down-weights the contribution of high-degree nodes without removing them.
3. This is an exploratory control (work plan family: exploratory; it does not enter the H1 verdict).
