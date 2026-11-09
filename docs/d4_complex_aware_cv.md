# D4 -- the complex-aware (v2) result against v1

**Date:** 2026-11-05  
**Script:** `scripts/run_d4_compare_v1_v2.py`  
**Tables:** `results/tables/d4_v1_vs_v2_comparison.tsv`, `results/tables/d4_v2_stats.tsv`

## Decision (the plan's Q3 rule)

* median |relative change| in delta-AUROC: **0.142** (max 0.455)
* arms whose sign flips: **0/24**
* H1 (all four interval families exclude 0): v1 **24/24**, v2 **24/24**
* cells equivalent to zero (TOST): v1 **6/24**, v2 **8/24**
* median repeat-level interval width, v2 / v1: **1.176**

**=> > 10 % -- the complex-aware result is the SECOND paper's main result and the first paper lists it as a limitation**

## Per arm

| network | label | v1 delta | v2 delta | rel. change | v1 95 % CI (repeat pct) | v2 95 % CI (repeat pct) |
|---|---|---|---|---|---|---|
| funmap | intogen2024 | +0.1001 | +0.0895 | -10.6% | [+0.0987, +0.1014] | [+0.0883, +0.0906] |
| funmap | intogen_temporal_new | +0.0196 | +0.0240 | +22.2% | [+0.0164, +0.0226] | [+0.0213, +0.0268] |
| funmap | ot_all | +0.0853 | +0.0504 | -41.0% | [+0.0835, +0.0873] | [+0.0486, +0.0522] |
| funmap | ot_nolit | +0.0780 | +0.0425 | -45.5% | [+0.0763, +0.0796] | [+0.0406, +0.0442] |
| intact | intogen2024 | +0.0651 | +0.0636 | -2.3% | [+0.0645, +0.0656] | [+0.0629, +0.0642] |
| intact | intogen_temporal_new | +0.0071 | +0.0066 | -7.0% | [+0.0049, +0.0090] | [+0.0048, +0.0084] |
| intact | ot_all | +0.0842 | +0.0672 | -20.2% | [+0.0834, +0.0850] | [+0.0664, +0.0679] |
| intact | ot_nolit | +0.0784 | +0.0580 | -26.0% | [+0.0775, +0.0793] | [+0.0573, +0.0588] |
| reactome | intogen2024 | +0.0429 | +0.0234 | -45.3% | [+0.0418, +0.0439] | [+0.0223, +0.0247] |
| reactome | intogen_temporal_new | -0.0125 | -0.0178 | -42.3% | [-0.0153, -0.0096] | [-0.0213, -0.0145] |
| reactome | ot_all | +0.1103 | +0.0875 | -20.6% | [+0.1087, +0.1118] | [+0.0854, +0.0894] |
| reactome | ot_nolit | +0.0996 | +0.0757 | -24.0% | [+0.0986, +0.1005] | [+0.0741, +0.0772] |
| rna_coexp | intogen2024 | +0.0521 | +0.0459 | -11.9% | [+0.0505, +0.0537] | [+0.0446, +0.0474] |
| rna_coexp | intogen_temporal_new | +0.0127 | +0.0146 | +15.2% | [+0.0087, +0.0167] | [+0.0107, +0.0185] |
| rna_coexp | ot_all | +0.0767 | +0.0576 | -24.9% | [+0.0748, +0.0784] | [+0.0556, +0.0594] |
| rna_coexp | ot_nolit | +0.0834 | +0.0536 | -35.7% | [+0.0814, +0.0854] | [+0.0516, +0.0553] |
| string_full700 | intogen2024 | +0.1133 | +0.1110 | -2.0% | [+0.1126, +0.1140] | [+0.1102, +0.1117] |
| string_full700 | intogen_temporal_new | +0.0604 | +0.0612 | +1.3% | [+0.0583, +0.0625] | [+0.0590, +0.0633] |
| string_full700 | ot_all | +0.1201 | +0.1115 | -7.2% | [+0.1193, +0.1209] | [+0.1108, +0.1121] |
| string_full700 | ot_nolit | +0.1267 | +0.1205 | -4.9% | [+0.1260, +0.1273] | [+0.1198, +0.1212] |
| string_phys700 | intogen2024 | +0.1189 | +0.1162 | -2.3% | [+0.1177, +0.1202] | [+0.1153, +0.1169] |
| string_phys700 | intogen_temporal_new | +0.0378 | +0.0418 | +10.5% | [+0.0356, +0.0400] | [+0.0397, +0.0441] |
| string_phys700 | ot_all | +0.1582 | +0.1372 | -13.3% | [+0.1572, +0.1592] | [+0.1356, +0.1387] |
| string_phys700 | ot_nolit | +0.1711 | +0.1495 | -12.6% | [+0.1702, +0.1720] | [+0.1475, +0.1513] |

## How to read this

1. **The interval widths are not comparable.** The complex-aware partitions are less diverse than the stratified ones (the group constraint removes freedom), so the v2 repeat-level intervals are 1.18x the v1 width for the same data. A narrower interval is not a stronger result; the comparison that matters is the point estimate and the sign.
2. **The estimator is not identical.** A v2 fold holds whole complexes, so the test set is a different kind of set; the comparison is between two valid estimators of the same quantity on the same data, not between the same estimator twice.
3. **What was removed.** The D1 audit measured that 83-95 % of the held-out positives that are in a curated complex had a complex-linked gene in the training fold. On the v2 splits that number is 0.02-0.18. The v2 result is therefore the one that does NOT benefit from the complex-partner leak.

## Caveats

* The complex map covers 18.3 % of the panel's genes, so the control only constrains the genes that have curated complex membership.
* The residual leakage comes from the sparsity guard (a gene keeps only its 10 largest complexes) and from the complexes skipped by the group-size cap; it is measured, not assumed away.
* The v2 analysis is an additional, exploratory analysis: the pre-registration fixes RWR on the stratified splits, so the v2 numbers belong to the exploratory family (amendment 4).
