# S2 -- attributing the stability of the prior to network topology

**Date:** 2026-11-05
**Script:** `scripts/run_stability_attribution.py`
**Inputs:** `results/tables/stability_per_gene.tsv` (S1), the six network
files, the four label files
**Outputs:** `results/tables/stability_features_per_gene.tsv` (267,696 rows),
`results/tables/stability_attribution.tsv` (18 rows),
`results/tables/stability_attribution_loo.csv` (20 rows),
`results/tables/stability_attribution_predictor_corr.tsv`,
`results/tables/network_features.tsv` (cached features),
`figures/figure5b_stability_attribution.png`
**Status:** exploratory (work plan Phase S, R19)

---

## 1. Question and design

S1 measured *what* is stable. S2 asks **which topological properties of
a gene predict its stability**, within a fixed (network, label) arm.

Candidate predictors, computed on each network's LCC: log degree,
k-core number, clustering coefficient, log betweenness centrality
(sampled, 300 source nodes), log PageRank, log community size, the
shortest-path distance from the gene to the arm's training seeds, the
local score density (from S1) and the log score magnitude.

Model: **within-arm OLS** (arm fixed effects) of log(stability) on the
standardised predictors, with cluster-robust standard errors by arm.
Arm fixed effects are used because they estimate exactly the within-arm
contrast the question asks, and they are far cheaper than a mixed model
on 267,696 rows. Two outcomes: `log_score_cv` (the CoV analogue, the
primary metric) and `log_rank_sd`.

Two things are reported for every predictor: the **joint** coefficient
and the **marginal** within-arm Spearman correlation. They disagree
substantially for degree, and the disagreement is itself a result
(section 4).

Transfer is tested by leave-one-network-out and leave-one-label-out:
the features are fitted on the other arms **without arm intercepts** and
used to predict the held-out arm's genes (out-of-sample R²).

## 2. Headline

| predictor | joint β (rank SD) | marginal ρ (rank SD) | joint β (score CV) | sign stable |
|---|---|---|---|---|
| **log local score density** | **+0.787** | **+0.478 (24/24 arms)** | +0.104 | 10/10 |
| **log score mean** | **+0.910** | +0.340 (24/24) | +0.524 | 10/10 |
| **log degree** | **−0.995** | **+0.038 (16/24)** | −0.675 | 10/10 |
| clustering coefficient | +0.089 | +0.088 (20/24) | +0.058 | 10/10 |
| log betweenness | −0.050 | +0.074 (21/24) | −0.022 | 10/10 |
| log PageRank | +0.265 | +0.023 (16/24) | +0.213 | 10/10 |
| log community size | −0.215 | +0.019 (13/24) | −0.145 | 10/10 |
| distance to the seeds | +0.061 | −0.125 (3/24) | +0.025 | 10/10 |
| k-core number | −0.001 (q = 0.98) | +0.044 | −0.008 (q = 0.78) | 3/10 |

(standardised coefficients; negative = the feature makes the stability
*higher*; q = BH-FDR over the nine predictors, per outcome; every row
above except k-core has q < 0.002)

**The stability of the prior is a property of the score distribution,
not of the network's connectivity.** The two predictors that survive
everything are the crowding of the score distribution and the score
magnitude itself; the classic connectivity descriptors (k-core) explain
nothing, and degree's apparent effect is a suppression artefact
(section 4).

## 3. Transfer to unseen networks and labels

Out-of-sample R² of the feature-only model (no arm intercept):

| held out | rank SD | score CV |
|---|---|---|
| network: funmap | 0.377 | 0.213 |
| network: rna_coexp | 0.241 | 0.172 |
| network: string_phys700 | 0.301 | 0.036 |
| network: intact | 0.558 | 0.129 |
| network: string_full700 | 0.528 | 0.270 |
| **network: reactome** | **0.316** | **−1.177** |
| label: intogen2024 | 0.567 | 0.341 |
| label: ot_all | 0.540 | 0.310 |
| label: ot_nolit | 0.539 | 0.315 |
| **label: intogen_temporal_new** | **0.291** | **0.167** |

Two findings:

1. **The rank-stability attribution transfers**: every held-out network
   and every held-out label keeps a positive out-of-sample R²
   (0.24–0.57). The explanation is not network-specific.
2. **The score-CV attribution does not transfer to reactome**
   (R² = −1.18, worse than predicting the mean), while it transfers to
   the other five networks (0.04–0.27) and to all four labels. Reactome's
   stability regime is unlike every other network's. Together with S1
   (reactome is the least stable network) and the three robustness
   controls (reactome is the only negative arm), this says Reactome's
   temporal failure is a **different kind of failure**, not a noisier
   version of the same thing — consistent with the curator-lag reading.
3. The **temporal label is the hardest to predict** from the other
   labels (R² 0.29 / 0.17), i.e. the least stable label is also the one
   whose stability pattern is the least transferable.

## 4. Degree: a suppression artefact, and why it must be reported with the marginal

Degree has essentially **no marginal association** with rank instability
(within-arm Spearman +0.038) but a large **negative joint** coefficient
(−0.995). The reason is collinearity: degree correlates 0.72 with the
RWR score and 0.43 with the local score density, and the score's own
coefficient is large and positive (+0.910). A predictor that is nearly
collinear with a strong one is unstable in the joint model and can flip
sign.

So the paper must **not** claim "hubs are stable". What the data support
is: *once you know where a gene sits in the score distribution and how
crowded that neighbourhood is, additional connectivity adds nothing*;
and the marginal statement — *degree alone explains essentially none of
the rank instability* — is the one that is safe to make. k-core, the
other connectivity descriptor, is null in both the marginal and the
joint model (q = 0.98).

## 5. Caveats

1. Betweenness is **sampled** (300 source nodes per network) because
   the exact algorithm is O(n·m) and would take hours on intact
   (17,918 nodes, 564,706 edges). The coefficients are therefore
   estimated with an extra Monte-Carlo error component; betweenness is
   the weakest of the significant predictors anyway.
2. `dist_to_seed` is the distance to the label's positive set in the
   network, not the fold-specific training seeds (the fold-level
   variation is second order). It is the only predictor with a sign that
   flips between the marginal and the joint model.
3. The model is **associational**. "Crowding makes ranks unstable" is
   supported by the marginal correlation (24/24 arms) and by the
   transfer test; it is not a causal claim (the score distribution is
   itself produced by the network).
4. The predictors are correlated (Spearman up to 0.98 between degree
   and k-core); the correlation matrix is released in
   `stability_attribution_predictor_corr.tsv` so that a reader can see
   which claims are fragile.
5. The temporal label's low transfer is expected -- it is the label with
   the smallest positive set (152 genes, 5-15 % of each universe) and
   the largest per-arm variance -- and should not be read as a failure
   of the attribution.

## 6. What comes next

* **S3** -- relate the *arm-level* stability (repeat-level SD of
  delta-AUROC) to network-level descriptors. This needs the two STRING
  sensitivity arms (n = 8 networks instead of 6) and is the part that
  speaks most directly to "which networks are stable".
* **S4** -- ML attribution with stability selection (the PCP + SVM
  analogue), reusing the complex-aware splitter from Phase D.
* **S5** -- repeat S1 + S2 on the ClinGen label and the transfer arms.
* **S6** -- manuscript section "§K Where the stability of the prior
  comes from" + Figure 5/5b.
