# Control 28 — the topology prior suite

**Date:** 2026-11-25  
**Script:** `scripts/run_c28_topology_suite.py`  
**Tables:** `c28_topology_suite.tsv`, `c28_topology_suite_summary.tsv`

## The question

The headline is a network prior against a degree baseline, and Control 19 already showed that closeness beats degree in every arm while degree is never the strongest classical baseline. This control asks the sharper question: **is the prior's increment over degree something any topology-aware score gets, or is it specific to a diffusion prior?**

Three new scorers, on the same folds, the same test genes and the same degree baseline as the primary endpoint: the betweenness centrality (a classical 'how central is this gene in the network's paths' score, approximated for the large networks), the number of training seeds one edge away, and the number within two edges.

## The mean delta-AUROC against degree, per arm

| network | label | the RWR prior | betweenness | the 1-hop seed count | the 2-hop seed count |
|---|---|---|---|---|---|
| funmap | intogen2024 | - | +0.0000 | +0.0852 | +0.0937 |
| funmap | intogen_temporal_new | - | -0.0293 | +0.0133 | +0.0400 |
| funmap | ot_all | - | -0.0114 | +0.0816 | +0.0507 |
| funmap | ot_nolit | - | -0.0099 | +0.0720 | +0.0574 |
| intact | intogen2024 | - | -0.0171 | +0.0804 | +0.0601 |
| intact | intogen_temporal_new | - | -0.0101 | -0.0361 | +0.0285 |
| intact | ot_all | - | -0.0126 | +0.0909 | +0.0367 |
| intact | ot_nolit | - | -0.0142 | +0.0752 | +0.0261 |
| reactome | intogen2024 | - | -0.0179 | +0.0081 | +0.0328 |
| reactome | intogen_temporal_new | - | -0.0125 | -0.0803 | -0.0091 |
| reactome | ot_all | - | -0.0395 | +0.0745 | +0.0594 |
| reactome | ot_nolit | - | -0.0370 | +0.0571 | +0.0659 |
| rna_coexp | intogen2024 | - | -0.0169 | +0.0611 | +0.0622 |
| rna_coexp | intogen_temporal_new | - | -0.0098 | -0.0056 | +0.0299 |
| rna_coexp | ot_all | - | +0.0078 | +0.0716 | +0.0505 |
| rna_coexp | ot_nolit | - | +0.0123 | +0.0686 | +0.0468 |
| string_full700 | intogen2024 | - | -0.0014 | +0.1051 | +0.1038 |
| string_full700 | intogen_temporal_new | - | -0.0039 | -0.0403 | +0.0773 |
| string_full700 | ot_all | - | -0.0068 | +0.0999 | +0.0635 |
| string_full700 | ot_nolit | - | -0.0041 | +0.1029 | +0.0623 |
| string_phys700 | intogen2024 | - | +0.0346 | +0.0922 | +0.1034 |
| string_phys700 | intogen_temporal_new | - | +0.0183 | -0.0290 | +0.0192 |
| string_phys700 | ot_all | - | +0.0300 | +0.1393 | +0.0937 |
| string_phys700 | ot_nolit | - | +0.0219 | +0.1438 | +0.0925 |

## Reading

* **betweenness**: the median margin against degree -0.0098, the range -0.0395 to +0.0346, positive in 7 of 24 arms
* **the 1-hop seed count**: the median margin against degree +0.0732, the range -0.0803 to +0.1438, positive in 19 of 24 arms
* **the 2-hop seed count**: the median margin against degree +0.0584, the range -0.0091 to +0.1038, positive in 23 of 24 arms

## Caveats

1. The betweenness is approximated with 100 sampled sources for the networks with more than 3,000 nodes (exact for the rest), so its values are noisy at the level of the sampling error; the comparison is at the arm level.
2. The seed-count scores are fold-dependent by construction (the seeds change), the betweenness is not (it is a property of the network alone).
3. **Exploratory** (amendment 4): no pre-registered endpoint is touched.
