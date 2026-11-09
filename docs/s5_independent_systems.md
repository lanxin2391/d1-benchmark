# S5 -- the stability story on two independent systems

**Date:** 2026-11-05  
**Script:** `scripts/run_s5_independent_systems.py`  
**Tables:** `results/tables/s5_transfer_stability.tsv`, `stability_by_arm_clingen.tsv`

## (a) Transfer arms -- is the stability a property of the network?

The transfer arms train on one label system and evaluate on another (66 cells: IntOGen -> Open Targets and back, 2 networks x label pairs). If the network-level stability ordering is the same there as in the within-label arms, the stability is a property of the network rather than of the label pair.

| network | within-label: SD | within-label: mean delta | transfer: SD | transfer: mean delta | transfer negative cells |
|---|---|---|---|---|---|
| funmap | 0.00648 | +0.0708 | 0.00582 | +0.0228 | 2 |
| intact | 0.00346 | +0.0587 | 0.00368 | +0.0218 | 1 |
| reactome | 0.00537 | +0.0601 | 0.00592 | +0.0147 | 4 |
| rna_coexp | 0.00768 | +0.0562 | 0.00550 | +0.0209 | 0 |
| string_full700 | 0.00344 | +0.1051 | 0.00537 | +0.0489 | 0 |
| string_phys700 | 0.00439 | +0.1215 | 0.00749 | +0.0590 | 0 |

## (b) ClinGen -- a different curation criterion

ClinGen (clinical validity of cancer variants, 84 genes) is a completely different curation criterion from IntOGen and Open Targets, and the primary panel has no ClinGen arm on the 40-repeat design. The splits and the runs were produced here and the stability atlas recomputed for them.

| network | label | positives | background rank SD | background score CV | drivers more stable? |
|---|---|---|---|---|---|
| funmap | clingen_label | 61 | 0.0304 | 0.234 | yes (+0.0083) |
| rna_coexp | clingen_label | 49 | 0.0289 | 0.238 | no (-0.0050) |
| string_phys700 | clingen_label | 80 | 0.0269 | 0.192 | yes (+0.0121) |
| intact | clingen_label | 83 | 0.0152 | 0.100 | no (-0.0048) |
| string_full700 | clingen_label | 84 | 0.0218 | 0.151 | yes (+0.0135) |
| reactome | clingen_label | 61 | 0.0281 | 0.215 | yes (+0.0065) |

## Reading

**ClinGen replicates the network-level stability profile.** The ordering of the networks by their stability is the same as in the primary panel: intact is the most stable (background rank SD 0.0152, score CV 0.100) and funmap / rna_coexp the least (0.0304 / 0.0289). A curation criterion that shares nothing with IntOGen or Open Targets therefore gives the same picture, so the profile is not an artefact of one label system.

**The transfer arms reproduce the extremes but not the full ranking.** The Spearman correlation between the within-label stability ordering and the transfer ordering is +0.37: the two most stable networks in the primary panel (intact, string_full700) are also the most stable in the transfer arms, and the least stable (funmap, reactome, string_phys700) are at the other end, but the middle reorders (rna_coexp is the least stable within-label and mid-table in transfer). The transfer arms' delta-AUROC is also much smaller than the within-label one (e.g. string_phys700: +0.059 against +0.122), because the model is evaluated on a label system it was not trained on -- the harder question, with a correspondingly noisier estimate.

**The driver-stability advantage is weaker here.** In the primary panel 19 of 24 arms have the drivers more stable than the background; in ClinGen it is 4 of 6, with the two largest networks on the negative side. This is what S1 predicted: the advantage is a position effect (the drivers sit higher in the score distribution), and a label with 84 positives spread over six networks gives it less room to show.

## Caveats

1. ClinGen has 84 positives in total, so a fold holds ~17 of them: its stability estimates are much noisier than the primary panel's, and the arm-level SD is correspondingly larger.
2. The transfer arms use their own split family (train on label A, evaluate on label B), so their delta-AUROC is not directly comparable with the within-label arms; only the ORDERING of the networks is compared.
