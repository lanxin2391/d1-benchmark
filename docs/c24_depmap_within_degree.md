# Control 24 — does the network add anything on the functional task, once hubness is taken out?

**Date:** 2026-11-18  
**Script:** `scripts/run_c24_depmap_within_degree.py`  
**Tables:** `c24_depmap_gene_level.tsv`, `c24_depmap_within_degree.tsv`

## The question

The released DepMap validation says the prior does not beat the degree baseline on functional essentiality (at K = 30: prior 5.8 %, degree 69.2 %, random 3.1 %). Two things were open. **Q1**: at the gene level, does the prior discriminate the essential genes better than degree does, and does it still do so after the degree effect is regressed out of its score? **Q2**: inside a degree stratum, does the prior's pick find essential genes more often than a random pick from the same stratum?

## Q1: the gene level (the median AUROC for functional essentiality, over the cell lines)

| threshold | the prior | degree | the residual (the prior minus its degree part) |
|---|---|---|---|
| -0.8 | 0.5000 | 0.4927 | 0.5006 |
| -0.5 | 0.5007 | 0.4942 | 0.5005 |
| -0.3 | 0.5010 | 0.4954 | 0.5006 |

### Q1 per arm (the threshold -0.5 rows, the median over the folds)

| network | label | the prior | degree | the residual | the size of the pool |
|---|---|---|---|---|---|
| funmap | intogen2024 | 0.4997 | 0.4910 | 0.4989 | 726 |
| funmap | intogen_temporal_new | 0.4971 | 0.4899 | 0.4993 | 777 |
| funmap | ot_all | 0.5030 | 0.4903 | 0.5007 | 732 |
| funmap | ot_nolit | 0.5025 | 0.4921 | 0.5002 | 736 |
| intact | intogen2024 | 0.5067 | 0.5041 | 0.5064 | 957 |
| intact | intogen_temporal_new | 0.5039 | 0.5024 | 0.5059 | 1022 |
| intact | ot_all | 0.5057 | 0.5073 | 0.4991 | 958 |
| intact | ot_nolit | 0.5040 | 0.5057 | 0.5034 | 964 |
| reactome | intogen2024 | 0.5067 | 0.4998 | 0.5038 | 416 |
| reactome | intogen_temporal_new | 0.4991 | 0.5017 | 0.5048 | 466 |
| reactome | ot_all | 0.5001 | 0.4977 | 0.5013 | 414 |
| reactome | ot_nolit | 0.4987 | 0.4956 | 0.4940 | 417 |
| rna_coexp | intogen2024 | 0.5012 | 0.5076 | 0.5004 | 649 |
| rna_coexp | intogen_temporal_new | 0.5034 | 0.5024 | 0.5003 | 695 |
| rna_coexp | ot_all | 0.5062 | 0.5039 | 0.4989 | 652 |
| rna_coexp | ot_nolit | 0.5017 | 0.5059 | 0.4981 | 656 |
| string_full700 | intogen2024 | 0.4997 | 0.4816 | 0.4946 | 915 |
| string_full700 | intogen_temporal_new | 0.4950 | 0.4780 | 0.5073 | 982 |
| string_full700 | ot_all | 0.4972 | 0.4780 | 0.4937 | 922 |
| string_full700 | ot_nolit | 0.5017 | 0.4770 | 0.4964 | 926 |
| string_phys700 | intogen2024 | 0.5003 | 0.4939 | 0.5011 | 754 |
| string_phys700 | intogen_temporal_new | 0.4972 | 0.4912 | 0.5057 | 816 |
| string_phys700 | ot_all | 0.4994 | 0.4935 | 0.4948 | 753 |
| string_phys700 | ot_nolit | 0.4961 | 0.4884 | 0.5036 | 760 |

## Q2: inside a degree stratum (the share of cell lines enriched, K = 30)

| threshold | stratum (by degree, 1 = lowest) | the prior's picks | degree's picks | the residual's picks | a random pick from the same stratum |
|---|---|---|---|---|---|
| -0.8 | 0 | 0.0% | 0.0% | 0.0% | 0.0% |
| -0.8 | 1 | 0.0% | 0.0% | 0.0% | 0.0% |
| -0.8 | 2 | 0.0% | 0.0% | 0.0% | 0.0% |
| -0.8 | 3 | 0.0% | 3.5% | 0.0% | 0.6% |
| -0.8 | 4 | 39.4% | 99.8% | 34.5% | 94.2% |
| -0.5 | 0 | 0.0% | 0.0% | 0.0% | 0.0% |
| -0.5 | 1 | 0.0% | 0.0% | 0.0% | 0.0% |
| -0.5 | 2 | 0.0% | 0.1% | 0.0% | 0.0% |
| -0.5 | 3 | 0.1% | 7.7% | 0.1% | 2.1% |
| -0.5 | 4 | 50.2% | 99.9% | 45.2% | 96.6% |
| -0.3 | 0 | 0.0% | 0.0% | 0.0% | 0.0% |
| -0.3 | 1 | 0.0% | 0.0% | 0.0% | 0.0% |
| -0.3 | 2 | 0.0% | 0.2% | 0.0% | 0.1% |
| -0.3 | 3 | 0.3% | 11.1% | 0.3% | 4.6% |
| -0.3 | 4 | 50.5% | 99.5% | 47.7% | 93.7% |

## Caveats

1. The pools are the held-out negatives of each fold, the same pool the released DepMap analysis used, so the numbers are comparable with it.
2. The residual is the ordinary-least-squares residual of the prior's score on log1p(degree) inside the pool.
3. Q2 uses one fold in five (the per-cell-line medians are stable across folds; Q1 uses all of them).
4. The pools differ in size between the labels (the label with many positives leaves few negatives), so Q2 can only be run where a stratum holds enough genes for a top-K pick.
5. **Exploratory** (amendment 4): nothing in the pre-registered endpoint is touched.
