# Control 7 -- the publication-count control, completed

**Date:** 2026-11-05  
**Script:** `scripts/run_control7_matched_null.py`  
**Table:** `results/tables/control7_matched_null.tsv`

## What was missing and what is now there

Control 7 (protocol §4.4) adjusts each driver label for the genes' literature volume and keeps the half of the positives whose residual literature score is above the median of the union (decision D-PUB-MARGIN). It covered **23 of the 24** primary arms: there was no pubcount-adjusted version of `intogen_temporal_new`. It is built now (152 → 76 positives, label effect β = +1.045), and **all five labels are re-run on the 40-repeat design** so the control is comparable with the result it is a control for (the earlier arms were on the 10-repeat design; that table is kept as `pubcount_adj_comparison_10rep.tsv`).

## The matched null

Trimming the positives is not enough on its own: the *background* still contains every popular gene, so a residual advantage could still be a popularity effect. The matched null replaces the background by a set of genes matched to the positives one-to-one on log10(n_pubmed + 1) (nearest-neighbour matching, without replacement, seeded), so the background has the same size and the same literature distribution as the positives.

Three numbers per arm: the pre-registered endpoint (the original label, standard background), the Control-7 quantity (the publication-adjusted label, standard background) and the matched null (the adjusted label, literature-matched background).

| network | label | original label | adjusted label | adjusted + matched background |
|---|---|---|---|---|
| funmap | intogen2024 | +0.1001 | +0.1036 | +0.0803 |
| funmap | ot_all | +0.0853 | +0.0715 | +0.0463 |
| funmap | ot_nolit | +0.0780 | +0.0718 | +0.0472 |
| funmap | intogen_temporal_new | +0.0196 | +0.0205 | +0.0009 |
| funmap | clingen_label | +nan | +0.1177 | +0.1230 |
| rna_coexp | intogen2024 | +0.0521 | +0.0782 | +0.0518 |
| rna_coexp | ot_all | +0.0767 | +0.0636 | +0.0530 |
| rna_coexp | ot_nolit | +0.0834 | +0.0608 | +0.0433 |
| rna_coexp | intogen_temporal_new | +0.0127 | +0.0601 | +0.0247 |
| rna_coexp | clingen_label | +nan | +0.0280 | +0.0079 |
| string_phys700 | intogen2024 | +0.1189 | +0.1087 | +0.0540 |
| string_phys700 | ot_all | +0.1582 | +0.1261 | +0.1095 |
| string_phys700 | ot_nolit | +0.1711 | +0.1213 | +0.1199 |
| string_phys700 | intogen_temporal_new | +0.0378 | +0.0413 | +0.0035 |
| string_phys700 | clingen_label | +nan | +0.1157 | +0.1154 |
| intact | intogen2024 | +0.0651 | +0.0826 | +0.0634 |
| intact | ot_all | +0.0842 | +0.0847 | +0.0719 |
| intact | ot_nolit | +0.0784 | +0.0726 | +0.0651 |
| intact | intogen_temporal_new | +0.0071 | +0.0014 | -0.0109 |
| intact | clingen_label | +nan | +0.0751 | +0.0787 |
| string_full700 | intogen2024 | +0.1133 | +0.0687 | +0.1127 |
| string_full700 | ot_all | +0.1201 | +0.0584 | +0.1114 |
| string_full700 | ot_nolit | +0.1267 | +0.0631 | +0.1221 |
| string_full700 | intogen_temporal_new | +0.0604 | +0.0410 | +0.0459 |
| string_full700 | clingen_label | +nan | +0.0618 | +0.2158 |
| reactome | intogen2024 | +0.0429 | +0.0419 | +0.0286 |
| reactome | ot_all | +0.1103 | +0.0666 | +0.0687 |
| reactome | ot_nolit | +0.0996 | +0.0568 | +0.0725 |
| reactome | intogen_temporal_new | -0.0125 | -0.0122 | -0.0174 |
| reactome | clingen_label | +nan | +0.1290 | +0.2002 |

## How to read it

* The **original label** column is the pre-registered endpoint.
* The **adjusted label** column is Control 7 proper: the publication-heavy half of the positives is trimmed, so the effect that survives is not a pure popularity effect.
* The **matched background** column is the hardest test: the background is matched to the positives one-to-one on log10(n_pubmed+1), so even the background cannot be distinguished by literature volume.

**What the numbers say.** On the IntOGen label the effect survives both controls (median over the six networks +0.0826 on the original label, +0.0804 on the publication-adjusted label, +0.0587 with the literature-matched background), so the prior's contribution there is not a literature effect.

On the **temporal** label the picture is different. Trimming the publication-heavy positives *raises* the effect (+0.0162 -> +0.0308; the trimmed positives were the hard ones), and the literature-matched background then removes it (+0.0022). The matched comparison on this label is small -- 76 positives against 76 matched negatives -- so the right reading is that the advantage on newly discovered drivers **is no longer clearly positive once literature volume is controlled on both sides**, not that it is exactly zero. The genes that are both new and well studied are the ones the prior ranks well; that is a property of the label, not of the method, and it belongs in the paper next to the curator-lag result, because the two explanations of the temporal label are not the same and they partly compete.

## Caveats

1. The matching is done per fold on the log10 publication score only; a richer matching (on several covariates) would need more degrees of freedom than a fold of 17–120 positives allows.
2. The temporal label's adjusted version keeps 76 positives, so its folds hold 15–16 positives: the matched null is noisy there by construction.
3. This is an **exploratory** control (amendment 4): it does not enter the H1 decision.
