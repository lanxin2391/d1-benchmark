# Two questions about the labels: why four arms fail the gene-cluster interval, and how many label opinions there really are

**Date:** 2026-11-12  
**Scripts:** `scripts/run_control15_16.py`  
**Tables:** `c15_size_matched_gene_ci.tsv`, `c16_label_redundancy.tsv`, `c16_label_pca.tsv`

## Control 15: the size-matched diagnostic

The pre-registered gene-cluster interval excludes zero in 20 of 24 arms and the four failures are all `intogen_temporal_new`, the label with the smallest positive set. Two explanations: the temporal label's effect is too small, or the interval is wide *because* the positive set is small.

The D-11 capped panel settles it: there all four labels have the same number of positives inside each arm.

| network | label | positives (standard) | gene CI excludes 0 (standard) | positives (capped) | gene CI excludes 0 (capped) |
|---|---|---|---|---|---|
| funmap | intogen2024 | 92 | yes | 21 | yes |
| funmap | intogen_temporal_new | 21 | no | 21 | no |
| funmap | ot_all | 79 | yes | 21 | no |
| funmap | ot_nolit | 75 | yes | 21 | yes |
| intact | intogen2024 | 125 | yes | 30 | yes |
| intact | intogen_temporal_new | 30 | no | 30 | no |
| intact | ot_all | 118 | yes | 30 | yes |
| intact | ot_nolit | 117 | yes | 30 | yes |
| reactome | intogen2024 | 75 | yes | 13 | no |
| reactome | intogen_temporal_new | 13 | no | 13 | no |
| reactome | ot_all | 74 | yes | 13 | yes |
| reactome | ot_nolit | 70 | yes | 13 | yes |
| rna_coexp | intogen2024 | 77 | yes | 18 | yes |
| rna_coexp | intogen_temporal_new | 18 | no | 18 | no |
| rna_coexp | ot_all | 73 | yes | 18 | yes |
| rna_coexp | ot_nolit | 70 | yes | 18 | yes |
| string_full700 | intogen2024 | 123 | yes | 29 | yes |
| string_full700 | intogen_temporal_new | 29 | yes | 29 | yes |
| string_full700 | ot_all | 119 | yes | 29 | yes |
| string_full700 | ot_nolit | 118 | yes | 29 | yes |
| string_phys700 | intogen2024 | 108 | yes | 23 | yes |
| string_phys700 | intogen_temporal_new | 23 | yes | 23 | yes |
| string_phys700 | ot_all | 107 | yes | 23 | yes |
| string_phys700 | ot_nolit | 103 | yes | 23 | yes |

**Reading.** Standard design: 20/24 arms exclude zero. Equal-positive-size panel: 18/24.

## Control 16: how many independent label opinions are there?

- the participation ratio of the overlap matrix: **1.54**; on the centred overlap: **1.00** (for 5 labels). The panel carries fewer independent curation opinions than its label count suggests.
- the eigenvalues of the overlap matrix: 0.00, 0.00, 0.00, 0.27, 0.93

| label | mean Jaccard with the others | R^2 of its arm means from the other labels | mean arm margin |
|---|---|---|---|
| intogen2024 | 0.272 | 0.84 | +0.0821 |
| ot_all | 0.251 | 1.00 | +0.1058 |
| ot_nolit | 0.212 | 1.00 | +0.1062 |
| intogen_temporal_new | 0.101 | 0.99 | +0.0208 |
| ncg | 0.246 | 0.99 | +0.0876 |

| label | PC1 | PC2 | PC3 | PC4 | PC5 |
|---|---|---|---|---|---|
| intogen2024 | -0.42 | -0.51 | -0.71 | -0.26 | +0.01 |
| ot_all | -0.40 | +0.50 | -0.19 | +0.15 | -0.73 |
| ot_nolit | -0.50 | +0.50 | +0.14 | -0.52 | +0.47 |
| intogen_temporal_new | -0.30 | -0.47 | +0.64 | -0.35 | -0.39 |
| ncg | -0.57 | -0.16 | +0.19 | +0.72 | +0.30 |

## Caveats

1. The capped panel changes the positive sets, so its margins are not the primary design's; it is used here only to hold the size constant.
2. The redundancy measures are descriptive: they say how much independent information the label panel carries, not how much the results depend on it.
3. Exploratory family (amendment 4).
