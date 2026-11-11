# G4 -- the endpoints the published evidence actually uses

**Date:** 2026-11-05  
**Script:** `scripts/run_pauc_top10.py`  
**Tables:** `results/tables/pauc_top10.tsv`, `results/tables/pauc_top10_summary.tsv`  
**Figure:** `figures/figure9_pauc_top10.png`

## Why

Guala & Sonnhammer 2017 compare AUCs only up to an FPR of 0.02 (partial AUCs) and report that the diffusion methods win only in the **top 10 %** of the ranking. Our benchmark reports the full AUROC, which is the metric on which propagation methods look worst -- so the published effect sizes are invisible to us. These two endpoints are therefore computed as **co-primary** endpoints, on the same folds, the same seeds and the same degree baseline as the primary one.

## Result (median over the 24 arms)

| endpoint | median margin over degree | arms positive |
|---|---|---|
| full AUROC (primary) | +0.0809 | 23/24 |
| partial AUC, FPR <= 0.02 | +0.0308 | 20/24 |
| precision in the top 10 % | +0.0412 | 20/24 |

## The per-arm table

| network | label | dAUROC | dpAUC@0.02 | d top-10 % |
|---|---|---|---|---|
| funmap | intogen2024 | +0.1001 | -0.0001 | +0.0276 |
| funmap | intogen_temporal_new | +0.0196 | -0.0027 | -0.0020 |
| funmap | ot_all | +0.0853 | +0.0741 | +0.0470 |
| funmap | ot_nolit | +0.0780 | +0.0765 | +0.0427 |
| intact | intogen2024 | +0.0651 | +0.0311 | +0.0363 |
| intact | intogen_temporal_new | +0.0071 | +0.0001 | -0.0005 |
| intact | ot_all | +0.0842 | +0.0534 | +0.0511 |
| intact | ot_nolit | +0.0784 | +0.0305 | +0.0402 |
| reactome | intogen2024 | +0.0429 | +0.0032 | +0.0422 |
| reactome | intogen_temporal_new | -0.0125 | -0.0031 | -0.0091 |
| reactome | ot_all | +0.1103 | +0.0551 | +0.1223 |
| reactome | ot_nolit | +0.0996 | +0.0335 | +0.1329 |
| rna_coexp | intogen2024 | +0.0521 | -0.0072 | +0.0166 |
| rna_coexp | intogen_temporal_new | +0.0127 | +0.0170 | -0.0020 |
| rna_coexp | ot_all | +0.0767 | +0.0059 | +0.0373 |
| rna_coexp | ot_nolit | +0.0834 | +0.0062 | +0.0349 |
| string_full700 | intogen2024 | +0.1133 | +0.1140 | +0.0803 |
| string_full700 | intogen_temporal_new | +0.0604 | +0.0211 | +0.0067 |
| string_full700 | ot_all | +0.1201 | +0.2246 | +0.1143 |
| string_full700 | ot_nolit | +0.1267 | +0.1759 | +0.1214 |
| string_phys700 | intogen2024 | +0.1189 | +0.0700 | +0.1193 |
| string_phys700 | intogen_temporal_new | +0.0378 | +0.0087 | +0.0055 |
| string_phys700 | ot_all | +0.1582 | +0.1908 | +0.1945 |
| string_phys700 | ot_nolit | +0.1711 | +0.1721 | +0.2056 |

## Caveats

1. **The low-FPR partial AUC is noisy by construction**: at FPR <= 0.02 a fold holds ~30-120 positives, so the truncated ROC is estimated from very few points. The intervals in `pauc_top10_summary.tsv` are the repeat-level ones; the arm-level noise is larger than for the full AUROC, and the paper should say so rather than quote the point estimates alone.
2. **The top-10 % precision is a coverage metric, not a ranking metric**: with a small positive set, the top 10 % of the test set contains 1-2 positives, so the metric is coarse and sensitive to the tie-break. It is reported because the published comparisons use it, not because it is precise.
3. These are **exploratory** endpoints (amendment 4). The pre-registered endpoint stays the full AUROC; the co-primary endpoints are reported alongside it and do not change the H1 verdict.
4. They answer a question the paper should now ask openly: *is the prior's advantage concentrated where the literature measures it?* If the answer is no, that is a finding about the method; if yes, the full-AUROC numbers we have been quoting are the wrong lens.
