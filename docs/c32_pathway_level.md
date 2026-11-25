# Control 32 — the pathway level

**Date:** 2026-11-25  
**Script:** `scripts/run_c32_pathway_level.py`  
**Tables:** `c32_pathway_enrichment_drivers.tsv`, `c32_pathway_picks.tsv`, `c32_pathway_permutation.tsv`

## The setup

The pathway-to-gene mapping comes from the raw Reactome release (the `Interaction context` column is the pathway), with the Ensembl ids mapped to symbols through HGNC, and it is restricted to the genes of the Reactome network so that the universe is the network's own gene set.

## 1. Which pathways are the drivers in?

| label | the pathways tested | the pathways enriched at q < 0.05 |
|---|---|---|
| intogen2024 | 4294 | 257 |
| intogen_temporal_new | 4294 | 0 |
| ot_all | 4294 | 992 |
| ot_nolit | 4294 | 887 |

## 2. Do the prior's picks concentrate on those pathways more than the degree baseline's?

| network | label | the prior: the enriched pathways | degree: the enriched pathways | the random baseline |
|---|---|---|---|---|
| funmap | intogen2024 | 76.7 | 93.5 | 91.2 |
| funmap | ot_all | 74.5 | 97.8 | 91.4 |
| funmap | ot_nolit | 81.9 | 94.1 | 93.0 |
| funmap | intogen_temporal_new | 96.2 | 92.4 | 92.8 |
| rna_coexp | intogen2024 | 119.8 | 109.0 | 85.4 |
| rna_coexp | ot_all | 123.6 | 108.8 | 87.9 |
| rna_coexp | ot_nolit | 127.5 | 101.0 | 87.2 |
| rna_coexp | intogen_temporal_new | 112.9 | 95.8 | 89.2 |
| string_phys700 | intogen2024 | 342.9 | 444.2 | 66.5 |
| string_phys700 | ot_all | 385.4 | 435.3 | 66.7 |
| string_phys700 | ot_nolit | 395.5 | 422.2 | 66.9 |
| string_phys700 | intogen_temporal_new | 169.9 | 425.4 | 67.1 |
| intact | intogen2024 | 163.2 | 123.0 | 89.3 |
| intact | ot_all | 175.0 | 126.5 | 89.5 |
| intact | ot_nolit | 165.1 | 121.2 | 91.3 |
| intact | intogen_temporal_new | 128.9 | 114.7 | 91.1 |
| string_full700 | intogen2024 | 335.0 | 286.1 | 86.3 |
| string_full700 | ot_all | 407.1 | 286.7 | 87.0 |
| string_full700 | ot_nolit | 402.4 | 290.9 | 86.4 |
| string_full700 | intogen_temporal_new | 219.2 | 291.7 | 86.9 |
| reactome | intogen2024 | 632.3 | 849.3 | 54.6 |
| reactome | ot_all | 682.9 | 861.0 | 54.9 |
| reactome | ot_nolit | 672.0 | 846.1 | 55.1 |
| reactome | intogen_temporal_new | 490.6 | 855.8 | 54.7 |

## Reading

## Caveats

1. The pathways come from Reactome, so the Reactome arm has a shared source with the pathway annotation; the other five networks do not.
2. The pathway-level null permutes the picks among the test genes, so it tests the enrichment statistic, not the pathway annotation itself.
3. **Exploratory** (amendment 4): no pre-registered endpoint is touched.
