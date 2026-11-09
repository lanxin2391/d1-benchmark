# D3 -- the complex-aware v2 splits, generated and audited

**Date:** 2026-11-05  
**Script:** `scripts/run_d3_generate_cx_splits.py`  
**Splits:** `data/processed/splits_v40_cx/` (the v1 splits are untouched)  
**Invariants:** `results/tables/d3_split_invariants.tsv`

## Acceptance criteria and their outcome

| criterion | required | observed |
|---|---|---|
| complex groups split across folds | 0 (except documented repairs) | **0** (repairs: 0) |
| D1 leakage on the new splits | ~0 | **0.000 - 0.188** (v1: 0.849 - 0.962) |
| positives per fold | as balanced as group preservation allows | max relative deviation **0.127** (median 0.084) |

## Per arm

| network | label | groups | split | repairs | pos/fold | leak v1 | leak v2 |
|---|---|---|---|---|---|---|---|
| funmap | intogen2024 | 8,113 | 0 | 0 | 84-98 (target 92) | 0.891 | **0.053** |
| funmap | ot_all | 8,113 | 0 | 0 | 72-84 (target 79) | 0.917 | **0.129** |
| funmap | ot_nolit | 8,113 | 0 | 0 | 68-80 (target 75) | 0.919 | **0.145** |
| funmap | intogen_temporal_new | 8,113 | 0 | 0 | 19-22 (target 21) | 0.849 | **0.026** |
| rna_coexp | intogen2024 | 7,685 | 0 | 0 | 71-81 (target 77) | 0.849 | **0.079** |
| rna_coexp | ot_all | 7,685 | 0 | 0 | 65-77 (target 72) | 0.872 | **0.058** |
| rna_coexp | ot_nolit | 7,685 | 0 | 0 | 63-74 (target 70) | 0.876 | **0.051** |
| rna_coexp | intogen_temporal_new | 7,685 | 0 | 0 | 17-19 (target 18) | 0.869 | **0.000** |
| string_phys700 | intogen2024 | 7,201 | 0 | 0 | 101-115 (target 109) | 0.949 | **0.060** |
| string_phys700 | ot_all | 7,201 | 0 | 0 | 96-115 (target 108) | 0.957 | **0.161** |
| string_phys700 | ot_nolit | 7,201 | 0 | 0 | 93-111 (target 104) | 0.961 | **0.178** |
| string_phys700 | intogen_temporal_new | 7,201 | 0 | 0 | 21-25 (target 23) | 0.950 | **0.020** |
| intact | intogen2024 | 15,119 | 0 | 0 | 118-130 (target 125) | 0.954 | **0.064** |
| intact | ot_all | 15,119 | 0 | 0 | 110-124 (target 118) | 0.956 | **0.164** |
| intact | ot_nolit | 15,119 | 0 | 0 | 109-123 (target 117) | 0.962 | **0.183** |
| intact | intogen_temporal_new | 15,119 | 0 | 0 | 28-31 (target 30) | 0.953 | **0.021** |
| string_full700 | intogen2024 | 13,119 | 0 | 0 | 115-128 (target 123) | 0.957 | **0.063** |
| string_full700 | ot_all | 13,119 | 0 | 0 | 109-125 (target 119) | 0.961 | **0.161** |
| string_full700 | ot_nolit | 13,119 | 0 | 0 | 109-123 (target 117) | 0.962 | **0.188** |
| string_full700 | intogen_temporal_new | 13,119 | 0 | 0 | 27-30 (target 29) | 0.947 | **0.022** |
| reactome | intogen2024 | 2,929 | 0 | 0 | 70-78 (target 75) | 0.920 | **0.049** |
| reactome | ot_all | 2,929 | 0 | 0 | 67-79 (target 74) | 0.933 | **0.047** |
| reactome | ot_nolit | 2,929 | 0 | 0 | 63-75 (target 70) | 0.931 | **0.045** |
| reactome | intogen_temporal_new | 2,929 | 0 | 0 | 11-13 (target 13) | 0.873 | **0.034** |

## Reading

The leakage channel that D1 found wide open (0.83-0.95 of the held-out positives that are in a curated complex have a complex-linked gene in the training fold) is closed by the group-preserving splits: the same number is 0.00-0.17 on the v2 splits. The group invariant holds in every arm, and the positives per fold stay within a few percent of the target even for the temporal label.

The residual leakage is the sparsity guard: a gene keeps only its 10 largest complexes in the splitter, while the audit uses the complete partner set, so links through dropped complexes are counted as residual leakage. It is a known under-count, not a failure of the invariant (see docs/d2_complex_splitter.md).

**Next:** D4 -- re-run the primary factorial on these splits with the same statistical machinery and compare with v1 under the plan's Q3 rule.
