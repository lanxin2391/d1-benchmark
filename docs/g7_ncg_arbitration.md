# G7 -- the NCG catalogue as an untouched arbitration set

**Date:** 2026-11-05  
**Script:** `scripts/run_g7_ncg_arbitration.py`  
**Table:** `results/tables/g7_ncg_arbitration.tsv`

## The question

The seeds are the label's own training-fold positives, so this is a **transfer** test: can the prior rank the drivers it was never seeded with? The arbitration positives are the NCG-known drivers (708) that are **not** in this label's positive set; the background is the rest of the test fold. The margin is again the prior over the degree baseline.

## Result

| network | label | arbitration positives per fold | dAUROC | dpAUC |
|---|---|---|---|---|
| funmap | intogen2024 | 33 | +0.0286 | -0.0042 |
| funmap | ot_all | 64 | +0.0477 | -0.0159 |
| funmap | ot_nolit | 73 | +0.0453 | -0.0158 |
| funmap | intogen_temporal_new | 87 | +0.0207 | -0.0153 |
| rna_coexp | intogen2024 | 26 | +0.0123 | -0.0018 |
| rna_coexp | ot_all | 54 | +0.0093 | -0.0298 |
| rna_coexp | ot_nolit | 61 | +0.0064 | -0.0252 |
| rna_coexp | intogen_temporal_new | 72 | +0.0190 | -0.0177 |
| string_phys700 | intogen2024 | 37 | +0.1269 | +0.0353 |
| string_phys700 | ot_all | 76 | +0.0722 | +0.0119 |
| string_phys700 | ot_nolit | 89 | +0.0590 | +0.0253 |
| string_phys700 | intogen_temporal_new | 107 | +0.0763 | +0.0455 |
| intact | intogen2024 | 46 | +0.0521 | +0.0199 |
| intact | ot_all | 92 | +0.0374 | +0.0028 |
| intact | ot_nolit | 105 | +0.0293 | +0.0069 |
| intact | intogen_temporal_new | 121 | +0.0193 | +0.0152 |
| string_full700 | intogen2024 | 47 | +0.1183 | +0.0451 |
| string_full700 | ot_all | 91 | +0.0758 | +0.0252 |
| string_full700 | ot_nolit | 106 | +0.0521 | +0.0352 |
| string_full700 | intogen_temporal_new | 121 | +0.0679 | +0.0632 |
| reactome | intogen2024 | 22 | +0.0470 | +0.0044 |
| reactome | ot_all | 45 | +0.0344 | -0.0073 |
| reactome | ot_nolit | 56 | +0.0361 | -0.0135 |
| reactome | intogen_temporal_new | 74 | +0.0089 | -0.0346 |

## Caveats

1. The arbitration positives are a *subset* of the genes that are not in our label, so they are not random: they are curated drivers that our label misses. The test asks whether the prior finds those, which is exactly what an arbitration catalogue is for (Boutros et al. 2014: the gold standard dominates the ranking).
2. The seeds are our label's positives, so this is not a label-free evaluation: it measures transfer from our label to an independent catalogue, not the quality of a label-free prior.
3. **Exploratory** (amendment 4): it does not enter the H1 decision.
