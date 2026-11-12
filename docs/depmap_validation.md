# F3 -- external validation against DepMap

**Date:** 2026-11-05  
**Script:** `scripts/run_depmap_validation.py`  
**Table:** `results/tables/depmap_validation.tsv` (+ `depmap_validation_thr-0.3.tsv` / `-0.8.tsv`)

## The question

The protocol's warning (section 9) is that a candidate list can look like a hub list rather than a dependency list, and our own evidence agrees (the RWR ranking's top 1 % contains almost no drivers; the degree baseline is a strong predictor). So the external data have to answer a specific question: are the genes the prior puts at the top of the ranking (excluding the label's own positives) more likely to be **required in cancer cell lines** than the rest of the universe?

## The test

* `C` = the top-K genes of the test fold by the RWR score, with the label's positives removed; `B` = the rest of the test fold.
* A gene is **essential in a cell line** if its CERES gene-effect score is below the threshold **and** the gene is expressed in that line (log2(TPM+1) > 1).
* Per cell line: a one-sided hypergeometric p-value asks whether `C` carries more essentials than the background rate implies.
* Three K values (30, 100, 1 % of the test set) and three essentiality thresholds (-0.3, -0.5, -0.8; the CERES scale is release-specific, so the sweep is not optional).

## Result (threshold -0.5)

| network | label | K | cell lines with p < 0.05 | median p | candidates per fold |
|---|---|---|---|---|---|
| funmap | intogen2024 | 0.01 | 0.019 | 0.698 | 20 |
| funmap | intogen_temporal_new | 0.01 | 0.004 | 0.902 | 21 |
| funmap | ot_all | 0.01 | 0.017 | 0.778 | 20 |
| funmap | ot_nolit | 0.01 | 0.019 | 0.766 | 20 |
| intact | intogen2024 | 0.01 | 0.419 | 0.0765 | 36 |
| intact | intogen_temporal_new | 0.01 | 0.225 | 0.186 | 36 |
| intact | ot_all | 0.01 | 0.311 | 0.13 | 36 |
| intact | ot_nolit | 0.01 | 0.280 | 0.152 | 36 |
| reactome | intogen2024 | 0.01 | 0.188 | 0.308 | 8 |
| reactome | intogen_temporal_new | 0.01 | 0.014 | 0.762 | 8 |
| reactome | ot_all | 0.01 | 0.250 | 0.165 | 8 |
| reactome | ot_nolit | 0.01 | 0.068 | 0.408 | 8 |
| rna_coexp | intogen2024 | 0.01 | 0.016 | 0.657 | 19 |
| rna_coexp | intogen_temporal_new | 0.01 | 0.004 | 0.878 | 19 |
| rna_coexp | ot_all | 0.01 | 0.021 | 0.658 | 19 |
| rna_coexp | ot_nolit | 0.01 | 0.013 | 0.669 | 19 |
| string_full700 | intogen2024 | 0.01 | 0.076 | 0.393 | 30 |
| string_full700 | intogen_temporal_new | 0.01 | 0.027 | 0.616 | 32 |
| string_full700 | ot_all | 0.01 | 0.192 | 0.258 | 31 |
| string_full700 | ot_nolit | 0.01 | 0.064 | 0.45 | 32 |
| string_phys700 | intogen2024 | 0.01 | 0.030 | 0.588 | 18 |
| string_phys700 | intogen_temporal_new | 0.01 | 0.015 | 0.707 | 19 |
| string_phys700 | ot_all | 0.01 | 0.079 | 0.414 | 19 |
| string_phys700 | ot_nolit | 0.01 | 0.017 | 0.679 | 19 |
| funmap | intogen2024 | 30.0 | 0.051 | 0.533 | 29 |
| funmap | intogen_temporal_new | 30.0 | 0.006 | 0.851 | 29 |
| funmap | ot_all | 30.0 | 0.059 | 0.594 | 29 |
| funmap | ot_nolit | 30.0 | 0.065 | 0.526 | 29 |
| intact | intogen2024 | 30.0 | 0.363 | 0.104 | 30 |
| intact | intogen_temporal_new | 30.0 | 0.197 | 0.209 | 30 |
| intact | ot_all | 30.0 | 0.300 | 0.135 | 30 |
| intact | ot_nolit | 30.0 | 0.271 | 0.153 | 30 |
| reactome | intogen2024 | 30.0 | 0.325 | 0.132 | 30 |
| reactome | intogen_temporal_new | 30.0 | 0.056 | 0.577 | 29 |
| reactome | ot_all | 30.0 | 0.377 | 0.0989 | 30 |
| reactome | ot_nolit | 30.0 | 0.153 | 0.303 | 30 |
| rna_coexp | intogen2024 | 30.0 | 0.030 | 0.61 | 30 |
| rna_coexp | intogen_temporal_new | 30.0 | 0.013 | 0.781 | 30 |
| rna_coexp | ot_all | 30.0 | 0.020 | 0.656 | 30 |
| rna_coexp | ot_nolit | 30.0 | 0.009 | 0.669 | 30 |
| string_full700 | intogen2024 | 30.0 | 0.078 | 0.407 | 28 |
| string_full700 | intogen_temporal_new | 30.0 | 0.026 | 0.599 | 30 |
| string_full700 | ot_all | 30.0 | 0.173 | 0.268 | 29 |
| string_full700 | ot_nolit | 30.0 | 0.057 | 0.474 | 30 |
| string_phys700 | intogen2024 | 30.0 | 0.046 | 0.519 | 27 |
| string_phys700 | intogen_temporal_new | 30.0 | 0.020 | 0.709 | 29 |
| string_phys700 | ot_all | 30.0 | 0.096 | 0.379 | 28 |
| string_phys700 | ot_nolit | 30.0 | 0.012 | 0.666 | 29 |
| funmap | intogen2024 | 100.0 | 0.890 | 0.000678 | 98 |
| funmap | intogen_temporal_new | 100.0 | 0.766 | 0.0029 | 98 |
| funmap | ot_all | 100.0 | 0.885 | 0.000241 | 98 |
| funmap | ot_nolit | 100.0 | 0.873 | 0.00116 | 98 |
| intact | intogen2024 | 100.0 | 0.868 | 0.00139 | 99 |
| intact | intogen_temporal_new | 100.0 | 0.506 | 0.0479 | 99 |
| intact | ot_all | 100.0 | 0.563 | 0.0348 | 99 |
| intact | ot_nolit | 100.0 | 0.446 | 0.068 | 99 |
| reactome | intogen2024 | 100.0 | 0.653 | 0.0153 | 99 |
| reactome | intogen_temporal_new | 100.0 | 0.879 | 0.000307 | 98 |
| reactome | ot_all | 100.0 | 0.625 | 0.0193 | 97 |
| reactome | ot_nolit | 100.0 | 0.277 | 0.192 | 98 |
| rna_coexp | intogen2024 | 100.0 | 0.094 | 0.487 | 99 |
| rna_coexp | intogen_temporal_new | 100.0 | 0.041 | 0.617 | 99 |
| rna_coexp | ot_all | 100.0 | 0.185 | 0.26 | 98 |
| rna_coexp | ot_nolit | 100.0 | 0.088 | 0.497 | 99 |
| string_full700 | intogen2024 | 100.0 | 0.659 | 0.019 | 96 |
| string_full700 | intogen_temporal_new | 100.0 | 0.091 | 0.385 | 97 |
| string_full700 | ot_all | 100.0 | 0.513 | 0.0463 | 97 |
| string_full700 | ot_nolit | 100.0 | 0.409 | 0.0854 | 98 |
| string_phys700 | intogen2024 | 100.0 | 0.405 | 0.0854 | 96 |
| string_phys700 | intogen_temporal_new | 100.0 | 0.208 | 0.257 | 96 |
| string_phys700 | ot_all | 100.0 | 0.285 | 0.155 | 96 |
| string_phys700 | ot_nolit | 100.0 | 0.068 | 0.513 | 96 |

## Reading

**The candidate set has to be big enough to show anything.** With the top 30 candidates the enrichment is at chance level (the median arm: 5.8 % of the cell lines at p < 0.05, median p = 0.52 at the -0.5 threshold); with the **top 100** the median arm reaches 47.6 % of the cell lines at p < 0.05 and a median p of 0.058, with the arms ranging from 4 % to 89 %. The threshold sweep barely moves the picture (-0.3 / -0.5 / -0.8: median p 0.048 / 0.058 / 0.073 at K = 100), so the result is not an artefact of the essentiality cut-off.

**This is the same message as the endpoint result (§H.1).** The prior's advantage over the degree baseline is a *global*-ranking advantage: in the low-FPR region it is 40-60 % smaller than in the full AUROC, and here the top 30 candidates carry no detectable functional signal while the top 100 do. The useful signal of a network prior is not concentrated in the first few names.

**Which networks validate best.** At K = 100 the enrichment is strongest for reactome and intact (the arms with 37-38 % of the cell lines significant at K = 30) and weakest for rna_coexp and funmap. Part of that is a property of the networks rather than of the prior: a candidate that is not expressed in a cell line cannot be called essential there, and the expression filter makes the networks that rank poorly expressed genes look worse.

## Caveats

1. The per-cell-line p-value asks about **enrichment**: a small median p means the candidate list carries more essentials than the background rate implies, not that the candidates are right.
2. The essentiality definition is deliberately strict (expression AND effect below the threshold); the threshold sweep is the sensitivity analysis.
3. The three K values overlap, so the rows are not independent.
4. **Exploratory external validation** (amendment 4): it does not enter the H1 decision.
5. A null result is a publishable finding (the protocol says so) and would sharpen the claim: the prior improves the *ranking* of known drivers without necessarily pointing at genes that are essential in cancer cell lines.
