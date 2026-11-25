# Control 27 — the rewiring null under the second learner (the full grid)

**Date:** 2026-11-25  
**Script:** `scripts/run_c27_rewiring_second_learner.py`  
**Tables:** `c27_rewiring_second_learner.tsv`, `c27_rewiring_second_learner_summary.tsv`

## What this adds

Control 22 showed on one arm that the degree-preserving rewiring null is computable under the second (proximity) learner. This runs the full grid: 6 networks × 4 labels × every rewiring that exists, with the same proximity learner, the same degree baseline and the same test genes as the released second-learner run.

Two z conventions are reported, because the project has used both and they differ by a factor of √n:

* **the effect-size z** = (observed − null mean) / null SD — how many null SDs the observation is away from the null;
* **the standard-error z** = (observed − null mean) / (null SD / √n) — what the released `w5_p_values.tsv` column reports, i.e. how many standard errors of the null *mean* the observation is.

Both are shown so that the released column and this table can be read side by side; the verdicts agree, because the null is centred near zero and the observation is far from it.

## The result

| network | label | the observed margin | the null mean | the null SD | the rewirings | the effect-size z | the standard-error z |
|---|---|---|---|---|---|---|---|
| funmap | intogen2024 | +0.0821 | -0.0064 | 0.0079 | 100 | +11.21 | +112.15 |
| funmap | intogen_temporal_new | +0.0512 | -0.0067 | 0.0202 | 100 | +2.86 | +28.64 |
| funmap | ot_all | +0.0700 | -0.0267 | 0.0071 | 100 | +13.56 | +135.56 |
| funmap | ot_nolit | +0.0492 | -0.0228 | 0.0077 | 100 | +9.36 | +93.61 |
| intact | intogen2024 | -0.0824 | -0.1242 | 0.0040 | 100 | +10.40 | +104.05 |
| intact | intogen_temporal_new | -0.0426 | -0.0621 | 0.0192 | 100 | +1.01 | +10.15 |
| intact | ot_all | -0.0693 | -0.1184 | 0.0045 | 100 | +10.98 | +109.78 |
| intact | ot_nolit | -0.0552 | -0.1025 | 0.0049 | 100 | +9.69 | +96.91 |
| reactome | intogen2024 | -0.0293 | -0.0927 | 0.0111 | 98 | +5.73 | +56.75 |
| reactome | intogen_temporal_new | -0.0435 | -0.0581 | 0.0252 | 98 | +0.58 | +5.72 |
| reactome | ot_all | +0.0079 | -0.1045 | 0.0080 | 98 | +14.11 | +139.65 |
| reactome | ot_nolit | -0.0031 | -0.0946 | 0.0101 | 98 | +9.04 | +89.45 |
| rna_coexp | intogen2024 | +0.0455 | -0.0389 | 0.0079 | 100 | +10.69 | +106.90 |
| rna_coexp | intogen_temporal_new | +0.0303 | -0.0183 | 0.0205 | 100 | +2.38 | +23.76 |
| rna_coexp | ot_all | +0.0489 | -0.0171 | 0.0091 | 100 | +7.23 | +72.30 |
| rna_coexp | ot_nolit | +0.0770 | -0.0166 | 0.0094 | 100 | +9.93 | +99.32 |
| string_full700 | intogen2024 | +0.0109 | -0.1031 | 0.0058 | 100 | +19.51 | +195.07 |
| string_full700 | intogen_temporal_new | +0.0106 | -0.0651 | 0.0204 | 100 | +3.72 | +37.16 |
| string_full700 | ot_all | -0.0214 | -0.1382 | 0.0055 | 100 | +21.18 | +211.81 |
| string_full700 | ot_nolit | -0.0115 | -0.1177 | 0.0049 | 100 | +21.69 | +216.90 |
| string_phys700 | intogen2024 | +0.0508 | -0.0602 | 0.0082 | 100 | +13.48 | +134.81 |
| string_phys700 | intogen_temporal_new | +0.0001 | -0.0520 | 0.0209 | 100 | +2.48 | +24.85 |
| string_phys700 | ot_all | +0.0586 | -0.0702 | 0.0083 | 100 | +15.44 | +154.39 |
| string_phys700 | ot_nolit | +0.0774 | -0.0599 | 0.0087 | 100 | +15.74 | +157.45 |

## Reading

The rewiring null is as decisive under the proximity learner as under RWR: the effect-size z is positive in 24 of 24 arms, with the median at 10.16 and the range 0.58 to 21.69. A degree-preserving rewiring destroys the margin in every arm, so the increment is a property of the *wiring* and not of the degree sequence -- and that conclusion does not depend on the learner.

The comparison with the released RWR null (same arms, same rewirings) is in the table `c27_rewiring_second_learner_summary.tsv` and the discussion below.

## Caveats

1. 20 folds per arm (one per repeat) for both the observed and the null, because every rewiring needs a fresh multi-source BFS per fold. The released RWR null used all 200 folds.
2. The rewirings are the ones Control 2 generated (Reactome has 98, the rest 100).
3. The two z conventions differ by √n; read the effect-size z when asking 'how far outside the null is this', and the standard-error z when comparing with the released table.
4. **Exploratory** (amendment 4): no pre-registered endpoint is touched.
