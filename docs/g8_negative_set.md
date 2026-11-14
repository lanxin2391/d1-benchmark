# G8 -- the negative set is not clean

**Date:** 2026-11-05  
**Script:** `scripts/run_g8_negative_set.py`  
**Tables:** `results/tables/g8_negative_contamination.tsv`, `g8_cleaned_background.tsv`, `g8_label_coverage_of_ncg.tsv`

## A. how much of our negative set is actually a curated driver?

| network | label | negatives | of them in NCG | NCG-known |
|---|---|---|---|---|
| funmap | intogen2024 | 389,120 | 39,320 (10.1%) | 6,600 (1.70%) |
| funmap | ot_all | 391,800 | 47,160 (12.0%) | 12,600 (3.22%) |
| funmap | ot_nolit | 392,600 | 48,600 (12.4%) | 14,600 (3.72%) |
| funmap | intogen_temporal_new | 403,400 | 52,560 (13.0%) | 17,360 (4.30%) |
| rna_coexp | intogen2024 | 350,680 | 39,160 (11.2%) | 5,080 (1.45%) |
| rna_coexp | ot_all | 351,600 | 45,600 (13.0%) | 10,600 (3.01%) |
| rna_coexp | ot_nolit | 352,120 | 46,880 (13.3%) | 12,240 (3.48%) |
| rna_coexp | intogen_temporal_new | 362,520 | 50,080 (13.8%) | 14,400 (3.97%) |
| string_phys700 | intogen2024 | 371,440 | 43,200 (11.6%) | 7,480 (2.01%) |
| string_phys700 | ot_all | 371,680 | 52,080 (14.0%) | 15,240 (4.10%) |
| string_phys700 | ot_nolit | 372,440 | 54,080 (14.5%) | 17,800 (4.78%) |
| string_phys700 | intogen_temporal_new | 388,560 | 59,440 (15.3%) | 21,320 (5.49%) |
| intact | intogen2024 | 691,800 | 66,240 (9.6%) | 9,200 (1.33%) |
| intact | ot_all | 693,160 | 76,840 (11.1%) | 18,400 (2.65%) |
| intact | ot_nolit | 693,320 | 78,920 (11.4%) | 21,120 (3.05%) |
| intact | intogen_temporal_new | 710,720 | 83,960 (11.8%) | 24,200 (3.40%) |
| string_full700 | intogen2024 | 610,760 | 62,280 (10.2%) | 9,400 (1.54%) |
| string_full700 | ot_all | 611,520 | 72,640 (11.9%) | 18,360 (3.00%) |
| string_full700 | ot_nolit | 611,800 | 74,720 (12.2%) | 21,080 (3.45%) |
| string_full700 | intogen_temporal_new | 629,480 | 79,920 (12.7%) | 24,320 (3.86%) |
| reactome | intogen2024 | 143,200 | 22,280 (15.6%) | 4,360 (3.04%) |
| reactome | ot_all | 143,280 | 27,280 (19.0%) | 8,840 (6.17%) |
| reactome | ot_nolit | 144,080 | 29,080 (20.2%) | 11,000 (7.63%) |
| reactome | intogen_temporal_new | 155,600 | 34,360 (22.1%) | 14,800 (9.51%) |

## B. does the margin survive a cleaned background?

Every NCG gene (all 2,372, all tiers) is removed from the background and the endpoint recomputed on the same folds and the same scores.

| network | label | dAUROC standard | dAUROC cleaned | dpAUC standard | dpAUC cleaned | dtop10% standard | dtop10% cleaned |
|---|---|---|---|---|---|---|---|
| funmap | intogen2024 | +0.1001 | +0.0247 | -0.0001 | -0.0141 | +0.0276 | -0.0009 |
| funmap | intogen_temporal_new | +0.0196 | +0.0402 | -0.0027 | -0.0153 | -0.0020 | -0.0005 |
| funmap | ot_all | +0.0853 | +0.1030 | +0.0741 | +0.1452 | +0.0470 | +0.0493 |
| funmap | ot_nolit | +0.0780 | +0.0979 | +0.0765 | +0.1318 | +0.0427 | +0.0482 |
| intact | intogen2024 | +0.0651 | +0.0304 | +0.0311 | +0.0040 | +0.0363 | +0.0007 |
| intact | intogen_temporal_new | +0.0071 | +0.0025 | +0.0001 | -0.0252 | -0.0005 | -0.0010 |
| intact | ot_all | +0.0842 | +0.0968 | +0.0534 | +0.0576 | +0.0511 | +0.0315 |
| intact | ot_nolit | +0.0784 | +0.0866 | +0.0305 | +0.0403 | +0.0402 | +0.0317 |
| reactome | intogen2024 | +0.0429 | +0.0140 | +0.0032 | -0.0581 | +0.0422 | -0.0047 |
| reactome | intogen_temporal_new | -0.0125 | +0.0073 | -0.0031 | -0.0527 | -0.0091 | -0.0100 |
| reactome | ot_all | +0.1103 | +0.1575 | +0.0551 | +0.1253 | +0.1223 | +0.0783 |
| reactome | ot_nolit | +0.0996 | +0.1378 | +0.0335 | +0.1010 | +0.1329 | +0.1049 |
| rna_coexp | intogen2024 | +0.0521 | +0.0299 | -0.0072 | +0.0207 | +0.0166 | +0.0003 |
| rna_coexp | intogen_temporal_new | +0.0127 | +0.0134 | +0.0170 | +0.0022 | -0.0020 | -0.0017 |
| rna_coexp | ot_all | +0.0767 | +0.1063 | +0.0059 | +0.0192 | +0.0373 | +0.0365 |
| rna_coexp | ot_nolit | +0.0834 | +0.1110 | +0.0062 | +0.0222 | +0.0349 | +0.0379 |
| string_full700 | intogen2024 | +0.1133 | +0.0780 | +0.1140 | +0.0253 | +0.0803 | +0.0027 |
| string_full700 | intogen_temporal_new | +0.0604 | +0.0327 | +0.0211 | +0.0019 | +0.0067 | -0.0001 |
| string_full700 | ot_all | +0.1201 | +0.1411 | +0.2246 | +0.2410 | +0.1143 | +0.0808 |
| string_full700 | ot_nolit | +0.1267 | +0.1388 | +0.1759 | +0.2164 | +0.1214 | +0.0919 |
| string_phys700 | intogen2024 | +0.1189 | +0.0768 | +0.0700 | -0.0061 | +0.1193 | +0.0065 |
| string_phys700 | intogen_temporal_new | +0.0378 | +0.0384 | +0.0087 | +0.0000 | +0.0055 | +0.0013 |
| string_phys700 | ot_all | +0.1582 | +0.1806 | +0.1908 | +0.2143 | +0.1945 | +0.1345 |
| string_phys700 | ot_nolit | +0.1711 | +0.1979 | +0.1721 | +0.2126 | +0.2056 | +0.1607 |

## C. how much of the NCG gold standard do our labels contain?

| network | label | positives | of them NCG-known | fraction |
|---|---|---|---|---|
| funmap | intogen2024 | 461 | 326 | 70.7% |
| funmap | ot_all | 394 | 176 | 44.7% |
| funmap | ot_nolit | 374 | 126 | 33.7% |
| funmap | intogen_temporal_new | 104 | 57 | 54.8% |
| rna_coexp | intogen2024 | 385 | 279 | 72.5% |
| rna_coexp | ot_all | 362 | 141 | 39.0% |
| rna_coexp | ot_nolit | 349 | 100 | 28.7% |
| rna_coexp | intogen_temporal_new | 89 | 46 | 51.7% |
| string_phys700 | intogen2024 | 544 | 414 | 76.1% |
| string_phys700 | ot_all | 538 | 220 | 40.9% |
| string_phys700 | ot_nolit | 519 | 156 | 30.1% |
| string_phys700 | intogen_temporal_new | 116 | 68 | 58.6% |
| intact | intogen2024 | 623 | 458 | 73.5% |
| intact | ot_all | 589 | 228 | 38.7% |
| intact | ot_nolit | 585 | 160 | 27.4% |
| intact | intogen_temporal_new | 150 | 83 | 55.3% |
| string_full700 | intogen2024 | 613 | 454 | 74.1% |
| string_full700 | ot_all | 594 | 230 | 38.7% |
| string_full700 | ot_nolit | 587 | 162 | 27.6% |
| string_full700 | intogen_temporal_new | 145 | 81 | 55.9% |
| reactome | intogen2024 | 373 | 296 | 79.4% |
| reactome | ot_all | 371 | 184 | 49.6% |
| reactome | ot_nolit | 351 | 130 | 37.0% |
| reactome | intogen_temporal_new | 63 | 35 | 55.6% |

## Reading

**The negative set is genuinely contaminated.** Between 10 % and 22 % of the held-out "negatives" are in the NCG catalogue -- curated driver genes we simply do not label -- and 1.3 % to 9.5 % are NCG-*known* high-confidence drivers. Reactome's arms are the worst (15-22 %), the dense protein networks the least (10 %).

**Read the next section before quoting these numbers.** The cleaned column below is a POST-HOC computation on the evaluation set, and it removes the genes that are positives of the label being evaluated (87 % of the IntOGen positives), so the "the margin roughly halves" reading is mostly a statement about a smaller evaluation set. The clean version -- background only, positives untouched, whole benchmark re-run -- leaves the headline unchanged (23/24 arms, median +0.082).

 independent gold standard; the Open Targets labels are not.** 70-79 % of the NCG-known drivers are in intogen2024 (and 52-59 % in the temporal label), against 27-37 % in ot_nolit -- the Open Targets panels are much broader, which is also why they carry the largest false-positive contamination.

## Correction (2026-11-12): how the cleaned numbers must be read

The cleaned column above is computed **post hoc on the evaluation set**:
every NCG gene leaves the evaluation set, *including the genes that
are positives of the label being evaluated*. Because NCG and IntOGen
overlap heavily (463 of the 708 known NCG drivers are IntOGen
positives), that removes **87 % of the IntOGen positives** (633 -> 80),
and the "the margin roughly halves" reading is mostly a statement about
a much smaller evaluation set.

The clean version of the question is Control 13
(`scripts/run_c9_clean_background_factorial.py`): the curated drivers
leave the **background only**, the positives are untouched, and the
whole benchmark is re-run. The result: **23 of 24 arms positive, no sign
flip, the median margin +0.0809 -> +0.0822**, and the gene-cluster
interval passes in 21 of 24 arms instead of 20. The contamination of the
negative set is real and it does not carry the result.

## Caveats

1. NCG is a *different* curation criterion (somatic driver roles from cancer sequencing screens), so a gene being in NCG does not make it a false positive in our labels -- it makes it a driver we do not label. That is the point of the correction, but it also means the cleaned background removes genes that may be perfectly good negatives for our question.
2. The correction removes the NCG genes from the background only; the positives are untouched.
3. The coverage table (part C) is an external comparison, not a performance measure: it says how much of an independent gold standard each of our label sets contains.
4. **Exploratory** (amendment 4): it does not enter the H1 decision.
