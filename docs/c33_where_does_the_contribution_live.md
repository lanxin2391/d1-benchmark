# Control 33 — where does the network's contribution live?

**Date:** 2026-11-26  
**Script:** `scripts/run_c33_where_does_contribution_live.py`  
**Tables:** `c33_residual_decomposition.tsv`, `c33_self_written_model.tsv`, `c33_summary.tsv`

## What this item is, after the repositioning

The original item (the second paper's H6) was a **conditional repair**: write a degree-orthogonalised model *if* the RWR advantage disappeared once the hubs were removed. That trigger was tested on 2026-11-05 and it did not fire — the degree-penalised walk is larger (+0.0979), not smaller — so as a repair the idea was retired. The useful version of the item is the **decomposition of the network's contribution on the primary task**:

1. the per-gene residual of the prior's score after regressing out the degree, after regressing out the publication count, and after regressing out both — is the network's own component still above chance for the driver label?
2. the **self-written model** (the other half of H6): a logistic model with the degree and the publication count as features, with and without the network score.

Both are computed on the same folds, the same test genes and the same endpoint as the primary analysis.

## 1. The AUROC of each variant (the median over the folds, per arm)

| network | label | the prior | the degree | the publication count | the prior, degree removed | the prior, visibility removed | the prior, both removed |
|---|---|---|---|---|---|---|---|
| funmap | intogen2024 | 0.6201 | 0.5219 | 0.8198 | 0.6408 | 0.4752 | 0.5807 |
| funmap | intogen_temporal_new | 0.5513 | 0.5352 | 0.7315 | 0.5367 | 0.4444 | 0.5181 |
| funmap | ot_all | 0.6631 | 0.5773 | 0.8583 | 0.6429 | 0.5215 | 0.5758 |
| funmap | ot_nolit | 0.6517 | 0.5728 | 0.8388 | 0.6253 | 0.5084 | 0.5680 |
| intact | intogen2024 | 0.7743 | 0.7084 | 0.8351 | 0.6773 | 0.5630 | 0.5928 |
| intact | intogen_temporal_new | 0.6461 | 0.6354 | 0.7562 | 0.4987 | 0.4404 | 0.4508 |
| intact | ot_all | 0.7867 | 0.7023 | 0.8526 | 0.7218 | 0.5860 | 0.6366 |
| intact | ot_nolit | 0.7499 | 0.6684 | 0.8223 | 0.7014 | 0.5737 | 0.6288 |
| reactome | intogen2024 | 0.7229 | 0.6796 | 0.7697 | 0.5892 | 0.5768 | 0.5322 |
| reactome | intogen_temporal_new | 0.5913 | 0.6058 | 0.6508 | 0.4719 | 0.4507 | 0.4490 |
| reactome | ot_all | 0.8128 | 0.7050 | 0.7743 | 0.7302 | 0.6811 | 0.6721 |
| reactome | ot_nolit | 0.7824 | 0.6831 | 0.7459 | 0.7029 | 0.6628 | 0.6596 |
| rna_coexp | intogen2024 | 0.6418 | 0.5882 | 0.8156 | 0.6183 | 0.5165 | 0.5438 |
| rna_coexp | intogen_temporal_new | 0.5380 | 0.5351 | 0.7183 | 0.4829 | 0.4226 | 0.4309 |
| rna_coexp | ot_all | 0.6224 | 0.5464 | 0.8545 | 0.6176 | 0.5028 | 0.5302 |
| rna_coexp | ot_nolit | 0.6157 | 0.5289 | 0.8239 | 0.6202 | 0.4992 | 0.5431 |
| string_full700 | intogen2024 | 0.8423 | 0.7274 | 0.8284 | 0.7773 | 0.7103 | 0.7225 |
| string_full700 | intogen_temporal_new | 0.6873 | 0.6272 | 0.7453 | 0.5553 | 0.4767 | 0.5022 |
| string_full700 | ot_all | 0.9029 | 0.7815 | 0.8314 | 0.8393 | 0.7991 | 0.7958 |
| string_full700 | ot_nolit | 0.8837 | 0.7583 | 0.8050 | 0.8275 | 0.7853 | 0.7879 |
| string_phys700 | intogen2024 | 0.7815 | 0.6634 | 0.8093 | 0.7232 | 0.6389 | 0.6391 |
| string_phys700 | intogen_temporal_new | 0.6084 | 0.5683 | 0.7184 | 0.5173 | 0.4303 | 0.4545 |
| string_phys700 | ot_all | 0.8364 | 0.6776 | 0.8052 | 0.8028 | 0.7178 | 0.7257 |
| string_phys700 | ot_nolit | 0.8242 | 0.6546 | 0.7767 | 0.8099 | 0.7188 | 0.7430 |

## 2. The self-written model

| network | label | the coefficient on the network score | the AUROC gain of adding the network to the degree+visibility model |
|---|---|---|---|
| funmap | intogen2024 | +0.2037 | +0.0051 |
| funmap | intogen_temporal_new | -0.0752 | +0.0018 |
| funmap | ot_all | +0.3774 | +0.0141 |
| funmap | ot_nolit | +0.3608 | +0.0135 |
| intact | intogen2024 | +0.1377 | +0.0038 |
| intact | intogen_temporal_new | -0.1443 | -0.0006 |
| intact | ot_all | +0.2615 | +0.0128 |
| intact | ot_nolit | +0.1756 | +0.0113 |
| reactome | intogen2024 | +0.1746 | +0.0056 |
| reactome | intogen_temporal_new | -0.1736 | +0.0081 |
| reactome | ot_all | +0.7380 | +0.0546 |
| reactome | ot_nolit | +0.4941 | +0.0538 |
| rna_coexp | intogen2024 | +0.2649 | +0.0060 |
| rna_coexp | intogen_temporal_new | -0.0932 | +0.0081 |
| rna_coexp | ot_all | +0.2449 | +0.0053 |
| rna_coexp | ot_nolit | +0.2240 | +0.0061 |
| string_full700 | intogen2024 | +0.4453 | +0.0303 |
| string_full700 | intogen_temporal_new | +0.0286 | +0.0041 |
| string_full700 | ot_all | +0.7337 | +0.0662 |
| string_full700 | ot_nolit | +0.6194 | +0.0774 |
| string_phys700 | intogen2024 | +0.3152 | +0.0183 |
| string_phys700 | intogen_temporal_new | -0.0428 | +0.0036 |
| string_phys700 | ot_all | +0.7399 | +0.0753 |
| string_phys700 | ot_nolit | +0.7061 | +0.1001 |

## Reading

## Caveats

1. The regressions are fitted on the test genes of each fold (the same population the AUROC is computed on), with no label information, so the residual is 'what the degree and the publication count do not explain *here*'.
2. The publication count is a per-gene covariate; the publication-matched null of F.12 is a different (and stronger) control, because it matches the publication count at the level of the bait rather than regressing it out of the score.
3. The logistic model is a plain ridge-stabilised IRLS fit with standardised features; the coefficient on the network score is therefore comparable across arms.
4. **Exploratory** (amendment 4): no pre-registered endpoint is touched.
