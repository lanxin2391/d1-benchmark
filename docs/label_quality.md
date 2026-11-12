# C7 -- label quality: which positives carry the margin?

**Date:** 2026-11-05  
**Scripts:** `scripts/run_c7_build_tiers.py`, `scripts/run_c7_tier_analysis.py`  
**Tables:** `label_quality_tiers.tsv`, `label_tier_delta_auroc.tsv`, `clingen_noise_anchor.tsv`, `label_quality_interval_width.tsv`

## The tiers

Every positive gene gets a tier from its own catalogue's evidence: ClinGen's review status (Definitive/Strong = high), IntOGen's number of distinct detecting methods (>= 4 = high), the NCG literature count (>= 10 = high), and for Open Targets whether the cancer associations are carried by a non-literature source.

| label | high | mid | low |
|---|---|---|---|
| clingen_label | 80 | 4 | 0 |
| intogen2024 | 252 | 203 | 177 |
| intogen_temporal_new | 11 | 74 | 66 |
| ncg_known | 115 | 313 | 280 |
| ot_all | 245 | 355 | 0 |
| ot_nolit | 236 | 364 | 0 |

## 1. The tier-restricted margin

The margin is recomputed with only one tier's positives counted as positives; the other tiers' positives fall into the background. The seeds, and therefore the prior, do not change, so this isolates the evaluation noise.

| network | label | high | mid | low |
|---|---|---|---|---|
| funmap | clingen_label | +0.0801 | -0.0486 | - |
| funmap | intogen2024 | +0.1367 | +0.0733 | +0.0638 |
| funmap | intogen_temporal_new | +0.0976 | +0.0154 | +0.0108 |
| funmap | ncg_known | +0.1158 | +0.0964 | +0.0309 |
| funmap | ot_all | +0.1103 | +0.0631 | - |
| funmap | ot_nolit | +0.0975 | +0.0617 | - |
| intact | clingen_label | +0.1133 | +0.0981 | - |
| intact | intogen2024 | +0.0702 | +0.0533 | +0.0669 |
| intact | intogen_temporal_new | -0.0312 | +0.0092 | +0.0110 |
| intact | ncg_known | +0.0506 | +0.0783 | +0.0577 |
| intact | ot_all | +0.0859 | +0.0809 | - |
| intact | ot_nolit | +0.0847 | +0.0724 | - |
| reactome | clingen_label | +0.1602 | - | - |
| reactome | intogen2024 | +0.0629 | +0.0135 | +0.0246 |
| reactome | intogen_temporal_new | +0.0589 | -0.0193 | -0.0165 |
| reactome | ncg_known | +0.0464 | +0.0596 | +0.0722 |
| reactome | ot_all | +0.1136 | +0.1007 | - |
| reactome | ot_nolit | +0.1219 | +0.0815 | - |
| rna_coexp | clingen_label | +0.0712 | - | - |
| rna_coexp | intogen2024 | +0.0607 | +0.0496 | +0.0362 |
| rna_coexp | intogen_temporal_new | -0.0470 | -0.0064 | +0.0531 |
| rna_coexp | ncg_known | +0.0599 | +0.0386 | +0.0357 |
| rna_coexp | ot_all | +0.0862 | +0.0669 | - |
| rna_coexp | ot_nolit | +0.0944 | +0.0734 | - |
| string_full700 | clingen_label | +0.0951 | +0.1433 | - |
| string_full700 | intogen2024 | +0.1082 | +0.1026 | +0.1222 |
| string_full700 | intogen_temporal_new | +0.0897 | +0.0572 | +0.0583 |
| string_full700 | ncg_known | +0.0733 | +0.1363 | +0.1597 |
| string_full700 | ot_all | +0.1287 | +0.1103 | - |
| string_full700 | ot_nolit | +0.1369 | +0.1164 | - |
| string_phys700 | clingen_label | +0.1635 | +0.2877 | - |
| string_phys700 | intogen2024 | +0.1247 | +0.1034 | +0.1107 |
| string_phys700 | intogen_temporal_new | +0.1819 | +0.0198 | +0.0353 |
| string_phys700 | ncg_known | +0.0989 | +0.1490 | +0.1225 |
| string_phys700 | ot_all | +0.1516 | +0.1557 | - |
| string_phys700 | ot_nolit | +0.1688 | +0.1654 | - |

## 2. The noise anchor (ClinGen's 17 % disagreement rate)

ClinGen's expert panels disagree on 17 % of multi-laboratory variants. Simulating exactly that rate -- a random 17 % of the positives demoted to background in every fold -- gives:

| network | n positives | clean | 17 % misassigned | change |
|---|---|---|---|---|
| funmap | 61 | +0.0735 | +0.0733 | -0.0002 |
| intact | 83 | +0.1124 | +0.1134 | +0.0010 |
| reactome | 61 | +0.1710 | +0.1736 | +0.0026 |
| rna_coexp | 49 | +0.0665 | +0.0647 | -0.0018 |
| string_full700 | 84 | +0.0974 | +0.0970 | -0.0005 |
| string_phys700 | 80 | +0.1686 | +0.1685 | -0.0001 |

The theoretical ceiling of a *perfect* scorer under the same misassignment is 1 - f/2 = 0.915, so a 17 % misassignment rate costs at most 8.5 points of AUROC however well the prior ranks.

## 3. Why some intervals are wide

The gene-cluster bootstrap interval width against the arm's positive-set size: Spearman rho = -0.966 (p = 2.01e-14); against the share of high-tier positives: rho = -0.475 (p = 0.0191).

## Caveats

1. The tier-restricted margin moves the other tiers' positives into the background, so it is not a label's own margin; it is a diagnostic of where the margin comes from.
2. The Open Targets tier uses the literature/non-literature split because almost every (target, cancer disease) pair in the datasource table carries exactly one source (88,510 of 89,111 in the first part), so a source-count tier would be constant.
3. **Exploratory** (amendment 4); the H1 decision is untouched.
