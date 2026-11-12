# The numbers the 1118 deck quotes

**Generated:** `scripts/dump_report_facts.py` (do not edit by hand).
Every figure on a slide comes from here, so the deck, the Q&A bank
and the manuscript cannot drift apart.

## 1. The primary endpoint (24 arms, alpha = 0.5)

- the margin over the degree baseline is positive in **23/24** arms; median **+0.0809**, range -0.0125 .. +0.1711
- the best arm: **string_phys700 x ot_nolit** (+0.1711); the worst: **reactome x intogen_temporal_new** (-0.0125)
- the single negative cell: [{'network': 'reactome', 'label_set': 'intogen_temporal_new', 'm': -0.0125}]

| network | label | dAUROC | SD |
|---|---|---|---|
| string_phys700 | ot_nolit | +0.1711 | 0.0191 |
| string_phys700 | ot_all | +0.1582 | 0.0150 |
| string_full700 | ot_nolit | +0.1267 | 0.0125 |
| string_full700 | ot_all | +0.1201 | 0.0126 |
| string_phys700 | intogen2024 | +0.1189 | 0.0170 |
| string_full700 | intogen2024 | +0.1133 | 0.0128 |
| reactome | ot_all | +0.1103 | 0.0208 |
| funmap | intogen2024 | +0.1001 | 0.0175 |
| reactome | ot_nolit | +0.0996 | 0.0250 |
| funmap | ot_all | +0.0853 | 0.0210 |
| intact | ot_all | +0.0842 | 0.0109 |
| rna_coexp | ot_nolit | +0.0834 | 0.0209 |
| intact | ot_nolit | +0.0784 | 0.0121 |
| funmap | ot_nolit | +0.0780 | 0.0188 |
| rna_coexp | ot_all | +0.0767 | 0.0205 |
| intact | intogen2024 | +0.0651 | 0.0086 |
| string_full700 | intogen_temporal_new | +0.0604 | 0.0273 |
| rna_coexp | intogen2024 | +0.0521 | 0.0163 |
| reactome | intogen2024 | +0.0429 | 0.0207 |
| string_phys700 | intogen_temporal_new | +0.0378 | 0.0352 |
| funmap | intogen_temporal_new | +0.0196 | 0.0292 |
| rna_coexp | intogen_temporal_new | +0.0127 | 0.0372 |
| intact | intogen_temporal_new | +0.0071 | 0.0195 |
| reactome | intogen_temporal_new | -0.0125 | 0.0501 |

**The H1 decision table** (`stats_h1_verdict_v40.tsv`):

```
  interval  n_arms  n_excl0  rule_met                                                  rule                                                                                                                    failing_arms
  fold_pct      24       24      True pre-registered H1: the 95% CI excludes 0 in EVERY arm                                                                                                                             NaN
repeat_pct      24       24      True pre-registered H1: the 95% CI excludes 0 in EVERY arm                                                                                                                             NaN
  repeat_t      24       24      True pre-registered H1: the 95% CI excludes 0 in EVERY arm                                                                                                                             NaN
repeat_bca      24       24      True pre-registered H1: the 95% CI excludes 0 in EVERY arm                                                                                                                             NaN
      gene      24       20     False pre-registered H1: the 95% CI excludes 0 in EVERY arm funmap x intogen_temporal_new; intact x intogen_temporal_new; reactome x intogen_temporal_new; rna_coexp x intogen_temporal_new
```

## 2. The control stack (the leakage controls)

- **control 2, degree-matched random seeds** (results/runs/a5/a5/ + results/tables/stats_fine_z_summary.tsv): 23/24 arms beat the null at z > 5; cell-level z = 3.09-226.9; the null mean is -0.052 .. +0.005, i.e. degree-matched random seeds are WORSE than the real seeds
- **control 3, rewired network** (results/tables/stats_fine_z_summary.tsv): 23/24 arms beat the rewired null at z > 5
- **D-09, temporal drift (50 seed resamples)** (results/tables/d09_table.tsv): -0.0095 +/- 0.0116
- **D-11, common-size cap per network** (results/tables/d11_table.tsv): -0.0125; reactome x intogen2024 flips from +0.0429 to negative
- **D-12, shared gene universe (1,605 genes)** (results/tables/d12_table.tsv): -0.0360
- **complex-aware cross-validation (Phase D)** (results/tables/d4_v1_v2_comparison.tsv): the leakage drops from 0.87-0.96 to 0.02-0.18; the median |relative change| of the margin is 0.142; H1 still holds in 24/24 and no sign flips
- **hub-capped RWR (E3)** (results/tables/hub_capped_rwr.tsv): the margin barely moves (+0.0849 -> +0.0809 clean; +0.0252 in the harshest cap, 8/24 arms flip) -> the degree-orthogonalised model of the second paper is NOT needed
- **Control 7, publication-count adjusted label** (results/tables/control7_matched_null.tsv): IntOGen +0.0826 -> +0.0804 -> +0.0587 (survives); temporal +0.0162 -> +0.0308 -> +0.0022 (collapses)

## 3. The results added since the manuscript draft (v1.0-rc3)

### G3 -- the bait-matched null (study visibility)

```
       network                label  n_positives_tested  frac_matched_null_above  median_p_matched  frac_unmatched_null_above  median_p_unmatched  frac_positives_beating_matched_null  frac_positives_beating_unmatched_null
        funmap          intogen2024               18440                 0.085412             0.450                   0.108460                0.30                             0.085412                               0.108460
        funmap               ot_all               15760                 0.140990             0.400                   0.173858                0.25                             0.140990                               0.173858
        funmap             ot_nolit               14960                 0.141110             0.400                   0.167246                0.30                             0.141110                               0.167246
        funmap intogen_temporal_new                4160                 0.038942             0.500                   0.057452                0.40                             0.038942                               0.057452
     rna_coexp          intogen2024               15400                 0.133182             0.350                   0.150325                0.30                             0.133182                               0.150325
     rna_coexp               ot_all               14480                 0.091713             0.375                   0.119061                0.30                             0.091713                               0.119061
     rna_coexp             ot_nolit               13960                 0.084885             0.400                   0.112966                0.30                             0.084885                               0.112966
     rna_coexp intogen_temporal_new                3560                 0.059270             0.500                   0.062640                0.45                             0.059270                               0.062640
string_phys700          intogen2024               21760                 0.122840             0.250                   0.248667                0.15                             0.122840                               0.248667
string_phys700               ot_all               21520                 0.298978             0.100                   0.402835                0.05                             0.298978                               0.402835
string_phys700             ot_nolit               20760                 0.320665             0.100                   0.402360                0.05                             0.320665                               0.402360
string_phys700 intogen_temporal_new                4640                 0.040517             0.450                   0.074138                0.35                             0.040517                               0.074138
        intact          intogen2024               24920                 0.083587             0.300                   0.235995                0.15                             0.083587                               0.235995
        intact               ot_all               23560                 0.125042             0.250                   0.282725                0.10                             0.125042                               0.282725
        intact             ot_nolit               23400                 0.108974             0.300                   0.247607                0.15                             0.108974                               0.247607
        intact intogen_temporal_new                6000                 0.049667             0.450                   0.092167                0.30                             0.049667                               0.092167
string_full700          intogen2024               24520                 0.161582             0.200                   0.354078                0.05                             0.161582                               0.354078
string_full700               ot_all               23760                 0.330598             0.100                   0.492172                0.05                             0.330598                               0.492172
string_full700             ot_nolit               23480                 0.326661             0.100                   0.447061                0.05                             0.326661                               0.447061
string_full700 intogen_temporal_new                5800                 0.059483             0.450                   0.118276                0.25                             0.059483                               0.118276
      reactome          intogen2024               14920                 0.094839             0.300                   0.178820                0.20                             0.094839                               0.178820
      reactome               ot_all               14840                 0.262871             0.150                   0.336456                0.10                             0.262871                               0.336456
      reactome             ot_nolit               14040                 0.236182             0.150                   0.294872                0.10                             0.236182                               0.294872
      reactome intogen_temporal_new                2520                 0.044444             0.400                   0.059524                0.40                             0.044444                               0.059524
```

### G7 -- the NCG arbitration test

- dAUROC positive in **24/24** arms, median **+0.0413**, range +0.0064 .. +0.1269
- dpAUC (FPR <= 0.02) positive in **13/24** arms, median **+0.0036**; the split is systematic: the four protein networks are positive in every label, funmap and RNA co-expression negative in all four

### G8 -- the negative set is not clean

- **9.6% to 22.1%** of the held-out negatives are curated NCG drivers; 1.3% to 9.5% are NCG-*known*
- the effect of cleaning the background (delta_AUROC, standard -> cleaned):

| network | label | standard | NCG-cleaned |
|---|---|---|---|
| funmap | intogen2024 | +0.1001 | +0.0247 |
| funmap | intogen_temporal_new | +0.0196 | +0.0402 |
| funmap | ot_all | +0.0853 | +0.1030 |
| funmap | ot_nolit | +0.0780 | +0.0979 |
| intact | intogen2024 | +0.0651 | +0.0304 |
| intact | intogen_temporal_new | +0.0071 | +0.0025 |
| intact | ot_all | +0.0842 | +0.0968 |
| intact | ot_nolit | +0.0784 | +0.0866 |
| reactome | intogen2024 | +0.0429 | +0.0140 |
| reactome | intogen_temporal_new | -0.0125 | +0.0073 |
| reactome | ot_all | +0.1103 | +0.1575 |
| reactome | ot_nolit | +0.0996 | +0.1378 |
| rna_coexp | intogen2024 | +0.0521 | +0.0299 |
| rna_coexp | intogen_temporal_new | +0.0127 | +0.0134 |
| rna_coexp | ot_all | +0.0767 | +0.1063 |
| rna_coexp | ot_nolit | +0.0834 | +0.1110 |
| string_full700 | intogen2024 | +0.1133 | +0.0780 |
| string_full700 | intogen_temporal_new | +0.0604 | +0.0327 |
| string_full700 | ot_all | +0.1201 | +0.1411 |
| string_full700 | ot_nolit | +0.1267 | +0.1388 |
| string_phys700 | intogen2024 | +0.1189 | +0.0768 |
| string_phys700 | intogen_temporal_new | +0.0378 | +0.0384 |
| string_phys700 | ot_all | +0.1582 | +0.1806 |
| string_phys700 | ot_nolit | +0.1711 | +0.1979 |

### GR-9 -- the cancer-type stratification

- **26/26** cancer types have a positive median margin (median over the types +0.0271)

## 4. The fifth label (NCG)

- the six NCG arms (alpha = 0.5): funmap +0.0771, rna_coexp +0.0426, string_phys700 +0.1352, intact +0.0673, string_full700 +0.1387, reactome +0.0647

- on the D-11 capped panel (equal positive-set size):

| network | intogen2024 | intogen_temporal_new | ncg_known | ot_all | ot_nolit |
|---|---|---|---|---|---|
| funmap | +0.1268 | +0.0196 | +0.1127 | +0.0191 | +0.0885 |
| intact | +0.0681 | +0.0071 | +0.0858 | +0.0664 | +0.0973 |
| reactome | -0.0110 | -0.0125 | +0.0893 | +0.0627 | +0.1649 |
| rna_coexp | +0.0829 | +0.0127 | +0.0720 | +0.0403 | +0.0956 |
| string_full700 | +0.1020 | +0.0604 | +0.1177 | +0.0785 | +0.1388 |
| string_phys700 | +0.0955 | +0.0378 | +0.1222 | +0.0938 | +0.1819 |

  **Reactome is positive on NCG (+0.0893) and negative on the two newest catalogues** -> the Reactome failure is label-specific (curator lag), not network-specific.

### The variance decomposition, 24 arms vs 30 arms

| variance source | % between-arm (24) | % between-arm (30) |
|---|---|---|
| network (between) | 24.82 % | 26.95 % |
| label (between) | 34.78 % | 30.75 % |
| network x label (interaction) | 40.07 % | 41.95 % |
| fold (within cell) | 0.34 % | 0.35 % |

## 5. Label quality (C7)

| label | high | mid | low |
|---|---|---|---|
| clingen_label | 80 | 4 | 0 |
| intogen2024 | 252 | 203 | 177 |
| intogen_temporal_new | 11 | 74 | 66 |
| ncg_known | 115 | 313 | 280 |
| ot_all | 245 | 355 | 0 |
| ot_nolit | 236 | 364 | 0 |

The tier-restricted margin (median over the six networks):

| label | high | mid | low |
|---|---|---|---|
| clingen_label | +0.1042 | +0.1207 | - |
| intogen2024 | +0.0892 | +0.0633 | +0.0653 |
| intogen_temporal_new | +0.0743 | +0.0123 | +0.0232 |
| ncg_known | +0.0666 | +0.0874 | +0.0650 |
| ot_all | +0.1120 | +0.0908 | - |
| ot_nolit | +0.1097 | +0.0775 | - |

**The ClinGen 17 % noise anchor:** the margin moves by -0.0018 .. +0.0026 when a random 17 % of the positives is demoted to background in every fold (a perfect scorer would lose 8.5 AUROC points).

## 6. The stability analysis (Phase S)

### stability_per_gene.tsv (267696 rows, columns: network, label, gene, y, degree, score_mean, score_sd, score_cv)

```
network       label  gene  y  degree   score_mean     score_sd  score_cv  rank_mean  rank_sd  rank_range  top10_hits  top50_hits  top1pct_hits  top5pct_hits  crowding_scott  local_score_density      provenance
 funmap intogen2024  A1BG  0   107.0 2.540538e-05 5.268542e-06  0.207379   0.362868 0.042644    0.156851           0           0             0             0             861                  351 protein-derived
 funmap intogen2024   A2M  0    93.0 2.992524e-05 6.479200e-06  0.216513   0.336096 0.041186    0.136533           0           0             0             0             697                  326 protein-derived
 funmap intogen2024 A2ML1  0    22.0 2.560937e-05 4.068137e-06  0.158853   0.359150 0.026784    0.092069           0           0             0             0             858                  357 protein-derived
 funmap intogen2024  AAAS  0    22.0 1.892030e-05 3.098196e-06  0.163750   0.409121 0.029210    0.102572           0           0             0             0            1214                  335 protein-derived
 funmap intogen2024  AACS  0     2.0 9.614867e-07 2.176916e-07  0.226411   0.835687 0.023119    0.079309           0           0             0             0            4283                  183 protein-derived
 funmap intogen2024 AAGAB  0     7.0 1.955915e-05 7.763371e-06  0.396918   0.422809 0.096826    0.256871           0           0             0             0            1178                  339 protein-derived
 funmap intogen2024  AAK1  0     3.0 2.310961e-06 4.333659e-07  0.187526   0.740118 0.025660    0.115037           0           0             0             0            4577                  264 protein-derived
 funmap intogen2024 AAMDC  0     1.0 1.326497e-06 3.474608e-07  0.261939   0.805629 0.034051    0.134766           0           0             0             0            4366                  214 protein-derived
 funmap intogen2024  AAMP  0     6.0 3.520186e-06 6.658168e-07  0.189142   0.681562 0.031911    0.101394           0           0             0             0            4787                  318 protein-derived
 funmap intogen2024  AAR2  0    13.0 1.567730e-05 2.470446e-06  0.157581   0.440931 0.029330    0.119651           0           0             0             0            1457                  366 protein-derived
 funmap intogen2024 AARS1  0    18.0 1.431643e-05 4.540743e-06  0.317170   0.463442 0.064246    0.185807           0           0             0             0            1609                  372 protein-derived
 funmap intogen2024 AARS2  0    16.0 3.052971e-06 5.523966e-07  0.180937   0.702488 0.028456    0.116215           0           0             0             0            4718                  284 protein-derived
```

### network_stability.tsv (8 rows, columns: network, net_delta_mean, net_stability_sd, net_stability_sd_max, net_cross_label_spread, net_icc, n_nodes, n_edges)

```
       network  net_delta_mean  net_stability_sd  net_stability_sd_max  net_cross_label_spread   net_icc  n_nodes  n_edges  density  mean_degree  mean_clustering  modularity  degree_assortativity  approx_mean_path_length  max_k_core  frac_in_2core  frac_in_10core  mean_label_coverage  min_label_coverage         provenance  is_sensitivity_arm  s2_transfer_r2_scorecv
        funmap        0.070759          0.006478              0.010128                0.035296 -0.141277    10189   196605 0.003788    38.591618         0.288042    0.756950              0.435776                 4.742462         115       0.856512        0.431348              0.68374              0.6233    protein-derived               False                0.212947
        intact        0.058695          0.003464              0.006512                0.035292 -0.165044    17918   564706 0.003518    63.032258         0.095770    0.370036             -0.059700                 2.924198          82       0.952450        0.761860              0.98316              0.9750    protein-derived               False                0.128862
      reactome        0.060056          0.005366              0.009861                0.056686 -0.205926     3953    20143 0.002579    10.191247         0.511802    0.759887              0.058277                 5.016409          24       0.847458        0.230458              0.58666              0.4145 curated-literature               False               -1.177396
     rna_coexp        0.056213          0.007679              0.012894                0.031980 -0.121940     9152   196172 0.004685    42.869755         0.384850    0.707356              0.392310                 4.690739         166       0.835555        0.472465              0.59240              0.5817        rna-derived               False                0.171597
string_full400        0.073187          0.002065              0.005113                0.012890 -0.216605    19486   929471 0.004896    95.398850         0.260567    0.496093              0.129033                 2.967229         169       0.995741        0.906446              0.99570              0.9917 curated-literature                True                     NaN
string_full700        0.105138          0.003442              0.006744                0.030298 -0.198267    15882   236712 0.001877    29.808840         0.383762    0.681319              0.279801                 4.090385         112       0.904357        0.538534              0.97812              0.9539 curated-literature               False                0.270122
string_full900        0.130704          0.004317              0.007761                0.039849 -0.193247    11693   100383 0.001469    17.169760         0.436939    0.772446              0.393367                 4.714104          85       0.860515        0.353117              0.91462              0.8289 curated-literature                True                     NaN
string_phys700        0.121509          0.004394              0.007229                0.060035 -0.196359     9830    85576 0.001771    17.411190         0.478059    0.783354              0.336463                 4.669437          87       0.847813        0.317091              0.86734              0.7632    protein-derived               False                0.036164
```

### stability_attribution.tsv (22 rows, columns: outcome, feature, beta, se, ci_low, ci_high, p, q_bh)

```
 outcome                 feature      beta       se    ci_low   ci_high             p          q_bh  sign_stability  marginal_spearman_median  marginal_spearman_min  marginal_spearman_max  arms_positive
score CV              log_degree  0.065954 0.015577  0.035424  0.096484  2.294617e-05  6.310195e-05             1.0                 -0.151816              -0.377281               0.027418              2
score CV                  k_core -0.036520 0.009464 -0.055070 -0.017970  1.139905e-04  2.507792e-04             1.0                 -0.156576              -0.399760               0.031041              2
score CV              clustering -0.002319 0.001916 -0.006074  0.001436  2.260755e-01  2.486831e-01             1.0                 -0.029662              -0.289652               0.088833              7
score CV         log_betweenness -0.000010 0.001835 -0.003608  0.003587  9.954555e-01  9.954555e-01             0.6                 -0.098029              -0.182252               0.021994              1
score CV            log_pagerank -0.011949 0.007614 -0.026872  0.002974  1.165759e-01  1.602918e-01             1.0                 -0.092900              -0.232741               0.030640              5
score CV      log_community_size  0.005324 0.003628 -0.001787  0.012434  1.422875e-01  1.739069e-01             0.9                 -0.060015              -0.449434               0.094802              7
score CV            dist_to_seed -0.015812 0.004495 -0.024622 -0.007002  4.350844e-04  6.837040e-04             1.0                 -0.094773              -0.356007               0.187713              7
score CV log_local_score_density  0.081345 0.018355  0.045370  0.117320  9.344896e-06  5.139693e-05             1.0                  0.097001              -0.078795               0.245356             20
score CV          log_score_mean -0.125563 0.029055 -0.182509 -0.068617  1.548952e-05  5.679489e-05             1.0                  0.144720              -0.043893               0.422651             18
score CV            log_score_cv  0.616390 0.024602  0.568171  0.664610 1.573834e-138 1.731217e-137             1.0                  1.000000               1.000000               1.000000             24
score CV      log_crowding_scott -0.086891 0.023498 -0.132946 -0.040836  2.174751e-04  3.987043e-04             1.0                 -0.125963              -0.353694               0.092102              7
 rank SD              log_degree -0.083405 0.020700 -0.123976 -0.042833  5.597382e-05  1.231424e-04             1.0                  0.038001              -0.098519               0.238072             16
```

## 7. DepMap external validation (figure 4)

- K = 0.01: median over the 24 arms of the fraction of cell lines with p < 0.05 = **2.8%** (chance = 5 %); median p = 0.602
- K = 30.0: median over the 24 arms of the fraction of cell lines with p < 0.05 = **5.8%** (chance = 5 %); median p = 0.523
- K = 100.0: median over the 24 arms of the fraction of cell lines with p < 0.05 = **47.6%** (chance = 5 %); median p = 0.058

## 8. The design (panel and labels)

| network | provenance | LCC genes | LCC edges |
|---|---|---|---|
| funmap | protein (CPTAC proteogenomics) | 10,189 | 196,605 |
| rna_coexp | RNA-derived | 9,152 | 196,172 |
| string_phys700 | protein (STRING physical >= 700) | 9,830 | 85,576 |
| intact | protein (IntAct) | 17,918 | 564,706 |
| string_full700 | curated (STRING all >= 700) | 15,882 | 236,712 |
| reactome | curated (Reactome v86) | 3,953 | 20,143 |

| label | source | positives |
|---|---|---|
| intogen2024 | mutation recurrence (MC3/IntOGen) | 633 |
| ot_all / ot_nolit | Open Targets 26.06, with / without the literature channel | 600 / 600 |
| intogen_temporal_new | the 152 newest IntOGen drivers | 152 |
| clingen_label | ClinGen cancer gene-disease validity | 84 |
| ncg_known | NCG known tier (Repana 2019) | 708 (capped per arm) |

- design: 6 networks x 4 labels x 3 alphas x 5 folds x 40 repeats = 1,200 evaluations per arm; the learner (RWR, alpha = 0.5) is frozen and pre-registered; the pre-registration is frozen at Zenodo DOI 10.5281/zenodo.23121762 and every departure is an amendment.

