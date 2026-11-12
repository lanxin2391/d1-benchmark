# GR-9 -- the arbitration test, stratified by cancer type

**Date:** 2026-11-05  
**Script:** `scripts/run_gr9_cancer_type_arbitration.py`  
**Tables:** `results/tables/gr9_cancer_type_arbitration.tsv`, `gr9_cancer_type_summary.tsv`, `gr9_background_effect.tsv`

## The question

Does the prior's advantage over the degree baseline depend on the cancer type? The arbitration positives are the NCG-known drivers of one cancer type that our label does not contain; the background is the rest of the test fold.

## Result

| cancer type | drivers | median dAUROC | range | arms positive |
|---|---|---|---|---|
| ampullary_adenocarcinoma | 16 | +0.0708 | +0.0085 .. +0.1894 | 11/11 |
| colorectal_adenocarcinoma | 20 | +0.0568 | -0.0278 .. +0.1782 | 17/18 |
| breast_cancer | 35 | +0.0546 | -0.0540 .. +0.1232 | 20/23 |
| nasopharyngeal_carcinoma | 16 | +0.0538 | -0.0202 .. +0.1408 | 16/18 |
| diffuse_large_B-cell_lymphoma | 21 | +0.0493 | -0.0441 .. +0.1507 | 16/18 |
| chronic_myelomonocytic_leukemia | 16 | +0.0421 | -0.0420 .. +0.1479 | 14/16 |
| triple_negative_breast_cancer | 28 | +0.0415 | -0.0213 .. +0.1103 | 14/16 |
| lung_cancer_(all_histologies) | 24 | +0.0379 | -0.0362 .. +0.1064 | 10/12 |
| bladder_cancer | 22 | +0.0366 | -0.0233 .. +0.1117 | 14/18 |
| melanoma | 20 | +0.0328 | -0.0442 .. +0.1199 | 14/17 |
| prostate_cancer | 20 | +0.0323 | -0.0072 .. +0.2160 | 16/17 |
| squamous_head_and_neck_cancer | 16 | +0.0284 | -0.0507 .. +0.0904 | 9/13 |
| pancreatic_cancer_(all_histologies) | 13 | +0.0276 | -0.0134 .. +0.1710 | 13/14 |
| acute_lymphoblastic_leukemia | 24 | +0.0266 | -0.0594 .. +0.0834 | 14/18 |
| cholangiocarcinoma | 14 | +0.0265 | -0.0913 .. +0.1767 | 12/13 |
| lung_adenocarcinoma | 18 | +0.0236 | -0.0733 .. +0.1209 | 11/16 |
| chronic_lymphocytic_leukemia | 16 | +0.0229 | -0.0202 .. +0.1177 | 13/17 |
| anaplastic_thyroid_carcinoma | 11 | +0.0220 | -0.0302 .. +0.2505 | 9/15 |
| acute_myeloid_leukemia | 17 | +0.0211 | -0.1747 .. +0.0981 | 12/16 |
| endometrial_cancer | 20 | +0.0172 | -0.1013 .. +0.1021 | 13/16 |
| pancreatic_ductal_adenocarcinoma | 14 | +0.0170 | -0.2580 .. +0.1552 | 14/18 |
| hepatocellular_carcinoma | 14 | +0.0151 | -0.0335 .. +0.1052 | 13/18 |
| glioblastoma | 14 | +0.0143 | -0.0188 .. +0.0851 | 8/10 |
| pan-glioma | 17 | +0.0131 | -0.1250 .. +0.0934 | 11/15 |
| T-cell_leukemia/lymphoma | 18 | +0.0101 | -0.0976 .. +0.0904 | 11/18 |
| biliary_tract_cancer | 14 | +0.0097 | -0.0480 .. +0.1193 | 9/12 |

## The background sub-sampling

The background is sub-sampled to 20,000 genes (fixed seed) to keep the 24 arms x cancer types tractable. The table below re-computes one arm with the full background, so the effect is on record.

| background | dAUROC |
|---|---|
| 20000 | +0.1269 |
| 200000 | +0.1269 |
| full | +0.1269 |

## Caveats

1. The NCG cancer-type annotation is itself curated and incomplete (426 of the 708 known drivers carry a type), so the strata are not exhaustive.
2. The pan-cancer entries are excluded from the stratification on purpose: they would dominate every stratum.
3. **Exploratory** (amendment 4); H1 untouched.
