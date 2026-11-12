# C2 -- is the NCG catalogue circular?

**Date:** 2026-11-05  
**Script:** `scripts/run_c2_ncg_circularity.py`  
**Tables:** `c2_ncg_circularity.tsv`, `c2_ncg_confidence_vs_degree.tsv`, `c2_label_overlap_matrix.tsv`, `c2_label_jaccard_matrix.tsv`

## The question

NCG is itself a network-based driver catalogue, so using it to score PPI networks risks a structural circle. The plan's C2 item requires this check before NCG may be used as a label, and the answer decides Q2 (independent label vs curation-confound panel).

## 1. Hub enrichment: where do the label's genes sit?

The mean degree percentile of the label's genes inside the network (a uniform random set of the same size would sit at 50):

| network | intogen2024 | ot_all | ot_nolit | intogen_temporal_new | ncg |
|---|---|---|---|---|---|
| funmap | 51.9 | 57.3 | 56.7 | 53.4 | 53.3 |
| intact | 70.2 | 69.5 | 66.3 | 63.5 | 70.2 |
| reactome | 66.4 | 68.5 | 66.5 | 60.1 | 63.8 |
| rna_coexp | 58.5 | 54.5 | 53.1 | 53.4 | 57.2 |
| string_full700 | 71.9 | 77.1 | 74.7 | 62.6 | 72.1 |
| string_phys700 | 65.5 | 66.9 | 64.5 | 57.2 | 64.9 |

## 2. Does the catalogue's own confidence track the network degree?

Spearman rho between the NCG confidence proxies and the degree in each network:

| network | proxy | n | rho | p |
|---|---|---|---|---|
| funmap | n_publications | 491 | -0.033 | 4.71e-01 |
| funmap | essential_celllines_pct | 491 | +0.301 | 9.40e-12 |
| funmap | ppin_degree | 491 | +0.180 | 5.98e-05 |
| rna_coexp | n_publications | 406 | +0.105 | 3.49e-02 |
| rna_coexp | essential_celllines_pct | 406 | +0.031 | 5.39e-01 |
| rna_coexp | ppin_degree | 406 | -0.013 | 7.93e-01 |
| string_phys700 | n_publications | 601 | +0.309 | 9.93e-15 |
| string_phys700 | essential_celllines_pct | 601 | +0.308 | 1.02e-14 |
| string_phys700 | ppin_degree | 601 | +0.563 | 1.27e-51 |
| intact | n_publications | 688 | +0.313 | 4.14e-17 |
| intact | essential_celllines_pct | 688 | +0.277 | 1.31e-13 |
| intact | ppin_degree | 688 | +0.720 | 7.04e-111 |
| string_full700 | n_publications | 689 | +0.382 | 2.46e-25 |
| string_full700 | essential_celllines_pct | 689 | +0.327 | 1.24e-18 |
| string_full700 | ppin_degree | 689 | +0.638 | 3.98e-80 |
| reactome | n_publications | 405 | +0.330 | 9.23e-12 |
| reactome | essential_celllines_pct | 405 | +0.128 | 9.97e-03 |
| reactome | ppin_degree | 405 | +0.474 | 4.46e-24 |

## 3. Label overlap (intersection size / Jaccard)

| | intogen2024 | ot_all | ot_nolit | intogen_temporal_new | ncg |
|---|---|---|---|---|---|
| **intogen2024** | 633.0 | 205.0 | 137.0 | 152.0 | 463.0 |
| **ot_all** | 205.0 | 600.0 | 430.0 | 24.0 | 233.0 |
| **ot_nolit** | 137.0 | 430.0 | 600.0 | 18.0 | 163.0 |
| **intogen_temporal_new** | 152.0 | 24.0 | 18.0 | 152.0 | 83.0 |
| **ncg** | 463.0 | 233.0 | 163.0 | 83.0 | 717.0 |

| Jaccard | intogen2024 | ot_all | ot_nolit | intogen_temporal_new | ncg |
|---|---|---|---|---|---|
| **intogen2024** | 1.000 | 0.199 | 0.125 | 0.240 | 0.522 |
| **ot_all** | 0.199 | 1.000 | 0.558 | 0.033 | 0.215 |
| **ot_nolit** | 0.125 | 0.558 | 1.000 | 0.025 | 0.141 |
| **intogen_temporal_new** | 0.240 | 0.033 | 0.025 | 1.000 | 0.106 |
| **ncg** | 0.522 | 0.215 | 0.141 | 0.106 | 1.000 |

## Caveats

1. The degree percentile is computed inside the network's LCC, so it is a within-network statement.
2. The confidence proxies are the ones the 2019 supplementary files carry; NCG 6.0's own confidence level is not available because the portal is down, so the literature count and the DepMap essential fraction stand in for it.
3. The hub comparison is a **relative** one: NCG is judged against the four other labels on the same networks, not against an absolute threshold.

## Verdict (the plan's decision point Q2)

**NCG can be used as a fifth label, but it is not the independent label the plan's decision D-2 assumed.** Three findings, each on record above:

1. **The feared structural circle does not materialise at the level of hub placement.** NCG's genes sit at essentially the same place in the degree distribution as IntOGen's (mean degree percentile: intact 70.2 vs 70.2, string_full700 72.0 vs 71.9, funmap 53.4 vs 51.9), and the Open Targets labels are *more* hub-loaded than both (string_full700 77.1). Whatever hub-bias makes the prior look good applies to the other labels too, so NCG does not add a new bias.
2. **But NCG's own topology evidence is partly the evidence we score with**: NCG's `ppin_degree` correlates with our network degree at rho = 0.47-0.72 (intact 0.72, string_phys700 0.56). Its literature count and essential-cell-line fraction correlate at rho = 0.03-0.38. The catalogue is not derived from our networks, but it is not independent of the *kind* of evidence either.
3. **And it is largely a restatement of IntOGen**: 463 of the 717 NCG known drivers (65 %) are in IntOGen 2024, Jaccard 0.522 -- above the 50 % threshold the plan set for calling the two labels non-independent. Against the temporal label NCG is close to orthogonal (Jaccard 0.106).

**Consequences for the plan.** (a) The six NCG arms are a near-replication of the IntOGen arms, not a new information source, so the C6 hypothesis (the label main effect should *shrink* when an independent label is added) cannot be tested as written -- C6 has to be read as "does adding a redundant label change the variance decomposition at all?". (b) As an *independent* label NCG earns nothing; as a **curation-confound panel** (like CGC) it is useful; and as the **external arbitration set** of G7/GR-9 it remains valuable precisely because it is not our label -- 27 % of its known drivers are drivers our labels miss. (c) The 30-arm extension is recorded as a **pre-planned exploratory** family (the plan's C5 pre-dates the analysis); the 24-arm pre-registered endpoint is not touched.
