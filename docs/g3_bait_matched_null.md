# G3 -- the bait-matched null and the degree profile of the labels

**Date:** 2026-11-05  
**Script:** `scripts/run_g3_bait_matched_null.py`  
**Tables:** `results/tables/g3_label_degree_enrichment.tsv`, `results/tables/g3_bait_matched_null.tsv`

## Why

Schaefer et al. 2015: the number of publications that mention a gene (the "bait" count) correlates with its degree (r = 0.520), and a null that ignores this "introduces circularity into network-based disease gene prediction". The degree baseline removes the symptom but not the common cause (how much the gene has been studied).

## A. the degree profile of the four label sets

| network | label | positives mean deg | negatives mean deg | ratio | positives median deg | negatives median deg |
|---|---|---|---|---|---|---|
| funmap | intogen2024 | 49.5 | 38.1 | 1.30 | 10 | 10 |
| funmap | ot_all | 52.0 | 38.1 | 1.37 | 20 | 10 |
| funmap | ot_nolit | 49.6 | 38.2 | 1.30 | 18 | 10 |
| funmap | intogen_temporal_new | 47.4 | 38.5 | 1.23 | 12 | 10 |
| rna_coexp | intogen2024 | 84.3 | 41.0 | 2.05 | 22 | 10 |
| rna_coexp | ot_all | 56.5 | 42.3 | 1.33 | 16 | 10 |
| rna_coexp | ot_nolit | 53.4 | 42.5 | 1.26 | 14 | 10 |
| rna_coexp | intogen_temporal_new | 49.1 | 42.8 | 1.15 | 11 | 10 |
| string_phys700 | intogen2024 | 27.4 | 16.8 | 1.63 | 16 | 6 |
| string_phys700 | ot_all | 27.9 | 16.8 | 1.66 | 16 | 6 |
| string_phys700 | ot_nolit | 25.5 | 17.0 | 1.51 | 14 | 6 |
| string_phys700 | intogen_temporal_new | 17.9 | 17.4 | 1.03 | 11 | 7 |
| intact | intogen2024 | 138.2 | 60.3 | 2.29 | 87 | 32 |
| intact | ot_all | 139.3 | 60.4 | 2.30 | 87 | 33 |
| intact | ot_nolit | 122.7 | 61.0 | 2.01 | 77 | 33 |
| intact | intogen_temporal_new | 96.9 | 62.7 | 1.54 | 72 | 34 |
| string_full700 | intogen2024 | 67.4 | 28.3 | 2.38 | 46 | 13 |
| string_full700 | ot_all | 74.6 | 28.1 | 2.66 | 56 | 13 |
| string_full700 | ot_nolit | 66.7 | 28.4 | 2.35 | 52 | 13 |
| string_full700 | intogen_temporal_new | 42.0 | 29.7 | 1.41 | 29 | 14 |
| reactome | intogen2024 | 19.4 | 9.2 | 2.10 | 11 | 5 |
| reactome | ot_all | 22.6 | 8.9 | 2.54 | 13 | 5 |
| reactome | ot_nolit | 20.0 | 9.2 | 2.16 | 12 | 5 |
| reactome | intogen_temporal_new | 13.9 | 10.1 | 1.37 | 9 | 5 |

## B. the bait-matched null

For every held-out positive: the fraction of random negatives that score above it, once with the negatives drawn at random and once with them **matched to the positive on log10(n_pubmed + 1)** (nearest neighbour, 20 draws). A positive that beats its matched nulls is ranked by the prior and not by how well it has been studied.

| network | label | positives beating the matched null | beating the unmatched null | median p (matched) | median p (unmatched) |
|---|---|---|---|---|---|
| funmap | intogen2024 | 8.5% | 10.8% | 0.450 | 0.300 |
| funmap | ot_all | 14.1% | 17.4% | 0.400 | 0.250 |
| funmap | ot_nolit | 14.1% | 16.7% | 0.400 | 0.300 |
| funmap | intogen_temporal_new | 3.9% | 5.7% | 0.500 | 0.400 |
| rna_coexp | intogen2024 | 13.3% | 15.0% | 0.350 | 0.300 |
| rna_coexp | ot_all | 9.2% | 11.9% | 0.375 | 0.300 |
| rna_coexp | ot_nolit | 8.5% | 11.3% | 0.400 | 0.300 |
| rna_coexp | intogen_temporal_new | 5.9% | 6.3% | 0.500 | 0.450 |
| string_phys700 | intogen2024 | 12.3% | 24.9% | 0.250 | 0.150 |
| string_phys700 | ot_all | 29.9% | 40.3% | 0.100 | 0.050 |
| string_phys700 | ot_nolit | 32.1% | 40.2% | 0.100 | 0.050 |
| string_phys700 | intogen_temporal_new | 4.1% | 7.4% | 0.450 | 0.350 |
| intact | intogen2024 | 8.4% | 23.6% | 0.300 | 0.150 |
| intact | ot_all | 12.5% | 28.3% | 0.250 | 0.100 |
| intact | ot_nolit | 10.9% | 24.8% | 0.300 | 0.150 |
| intact | intogen_temporal_new | 5.0% | 9.2% | 0.450 | 0.300 |
| string_full700 | intogen2024 | 16.2% | 35.4% | 0.200 | 0.050 |
| string_full700 | ot_all | 33.1% | 49.2% | 0.100 | 0.050 |
| string_full700 | ot_nolit | 32.7% | 44.7% | 0.100 | 0.050 |
| string_full700 | intogen_temporal_new | 5.9% | 11.8% | 0.450 | 0.250 |
| reactome | intogen2024 | 9.5% | 17.9% | 0.300 | 0.200 |
| reactome | ot_all | 26.3% | 33.6% | 0.150 | 0.100 |
| reactome | ot_nolit | 23.6% | 29.5% | 0.150 | 0.100 |
| reactome | intogen_temporal_new | 4.4% | 6.0% | 0.400 | 0.400 |

## Reading

**The label sets are degree-enriched** (mean degree of the positives over that of the negatives, inside the same LCC): 1.0-2.7x, and the temporal label is the least enriched.

**Most of the prior's margin over random genes is study visibility.** The median arm has 17.6 % of its held-out positives above 20 random negatives, but only 11.6 % above literature-matched ones (chance = 5 %). Per label: intogen2024 20.7 -> 10.9 %, ot_all 31.0 -> 20.2 %, ot_nolit 27.1 -> 18.9 %, and **intogen_temporal_new 6.8 -> 4.7 %, i.e. exactly chance**.

So on the three curated label systems a real part of the advantage survives (four times chance on ot_all), while on the temporal label the advantage is a study-visibility effect and not a topological one. This is the same conclusion the Control-7 matched null reaches from the other side, and it is the paper's most consequential limitation.

## Caveats

1. The bait count here is the project's publication count per gene (log10(n_pubmed + 1)), which is a study-visibility proxy rather than Schaefer's MeSH "bait" count. The two are related but not identical, so this is the same argument with a slightly different covariate.
2. The matching is one-to-one on the literature score within a fold; genes with no publication record are dropped from the test (reported as the number of positives tested).
3. The test asks about the positives being ranked above literature-matched non-drivers. It is **not** a test of whether the positives are functional (that is the DepMap validation, F3).
4. **Exploratory** (amendment 4): it does not enter the H1 decision.
