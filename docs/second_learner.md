# Control 17 -- a second, fixed learner

**Date:** 2026-11-12  
**Script:** `scripts/run_second_learner.py`  
**Tables:** `second_learner_margins.tsv`, `second_learner_summary.tsv`

## The learner

The **proximity learner**: score(g) = −distance(g, nearest training seed), the unweighted single-source distance on the same network. It is a deliberately different family from diffusion -- local, greedy, no restart, no stationary distribution -- and it is the family Guala & Sonnhammer's benchmark favours. Same splits, same folds, same degree baseline, same endpoint.

## The result

| network | label | proximity margin | RWR margin |
|---|---|---|---|
| string_phys700 | ot_nolit | +0.0814 | +0.1711 |
| funmap | intogen2024 | +0.0804 | +0.1001 |
| rna_coexp | ot_nolit | +0.0706 | +0.0834 |
| funmap | ot_all | +0.0624 | +0.0853 |
| string_phys700 | ot_all | +0.0621 | +0.1582 |
| rna_coexp | ot_all | +0.0605 | +0.0767 |
| funmap | ot_nolit | +0.0561 | +0.0780 |
| string_phys700 | intogen2024 | +0.0483 | +0.1189 |
| funmap | intogen_temporal_new | +0.0369 | +0.0196 |
| rna_coexp | intogen2024 | +0.0367 | +0.0521 |
| rna_coexp | intogen_temporal_new | +0.0272 | +0.0127 |
| string_phys700 | intogen_temporal_new | +0.0167 | +0.0378 |
| string_full700 | intogen_temporal_new | +0.0057 | +0.0604 |
| reactome | ot_all | +0.0041 | +0.1103 |
| string_full700 | intogen2024 | +0.0031 | +0.1133 |
| reactome | ot_nolit | +0.0002 | +0.0996 |
| string_full700 | ot_nolit | -0.0132 | +0.1267 |
| string_full700 | ot_all | -0.0248 | +0.1201 |
| reactome | intogen_temporal_new | -0.0264 | -0.0125 |
| reactome | intogen2024 | -0.0301 | +0.0429 |
| intact | intogen_temporal_new | -0.0422 | +0.0071 |
| intact | ot_nolit | -0.0515 | +0.0784 |
| intact | ot_all | -0.0695 | +0.0842 |
| intact | intogen2024 | -0.0815 | +0.0651 |

## Reading

- the proximity learner is positive in **16/24** arms (median **+0.0112**); the RWR margin on the same arms is +0.0809
- the rank correlation between the two learners' margins: Spearman +0.304

## Caveats

1. The two learners are not equally powerful by construction: diffusion uses the whole network and the seeds together, the proximity learner uses one step. The comparison is therefore about whether the increment over degree survives a change of family, not about which learner is better.
2. Genes unreachable from every seed get the maximum finite distance; their share is reported per arm in the table.
3. Exploratory family (amendment 4): the pre-registered endpoint is untouched.
