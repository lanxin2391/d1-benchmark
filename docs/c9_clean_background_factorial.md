# C9 -- the factorial re-run with a cleaned negative set

**Date:** 2026-11-12  
**Script:** `scripts/run_c9_clean_background_factorial.py`  
**Table:** `results/tables/c9_clean_background_comparison.tsv`

## What was done

G8 measured the contamination of our negative set post hoc, on the evaluation set of the existing runs. This is the stronger version: **the benchmark is re-run from the splits with the curated drivers removed from the universe.** The NCG genes (all three tiers) that are negatives of a label leave the split files; the positives are untouched, so the comparison isolates the background. The learner, the folds, the repeats and the alphas are unchanged.

- negative rows removed: **5,072,800** (2280 distinct genes)

## The margin, arm by arm

| network | label | standard | cleaned | change |
|---|---|---|---|---|
| string_phys700 | ot_nolit | +0.1711 | +0.1741 | +0.0030 |
| string_phys700 | ot_all | +0.1582 | +0.1624 | +0.0042 |
| string_full700 | ot_nolit | +0.1267 | +0.1269 | +0.0002 |
| string_full700 | ot_all | +0.1201 | +0.1214 | +0.0013 |
| string_phys700 | intogen2024 | +0.1189 | +0.1289 | +0.0100 |
| string_full700 | intogen2024 | +0.1133 | +0.1183 | +0.0050 |
| reactome | ot_all | +0.1103 | +0.1119 | +0.0016 |
| funmap | intogen2024 | +0.1001 | +0.1046 | +0.0045 |
| reactome | ot_nolit | +0.0996 | +0.1028 | +0.0032 |
| funmap | ot_all | +0.0853 | +0.0889 | +0.0036 |
| intact | ot_all | +0.0842 | +0.0853 | +0.0011 |
| rna_coexp | ot_nolit | +0.0834 | +0.0832 | -0.0002 |
| intact | ot_nolit | +0.0784 | +0.0793 | +0.0009 |
| funmap | ot_nolit | +0.0780 | +0.0813 | +0.0033 |
| rna_coexp | ot_all | +0.0767 | +0.0768 | +0.0001 |
| intact | intogen2024 | +0.0651 | +0.0671 | +0.0020 |
| string_full700 | intogen_temporal_new | +0.0604 | +0.0681 | +0.0077 |
| rna_coexp | intogen2024 | +0.0521 | +0.0533 | +0.0012 |
| reactome | intogen2024 | +0.0429 | +0.0466 | +0.0037 |
| string_phys700 | intogen_temporal_new | +0.0378 | +0.0490 | +0.0112 |
| funmap | intogen_temporal_new | +0.0196 | +0.0239 | +0.0043 |
| rna_coexp | intogen_temporal_new | +0.0127 | +0.0152 | +0.0025 |
| intact | intogen_temporal_new | +0.0071 | +0.0089 | +0.0018 |
| reactome | intogen_temporal_new | -0.0125 | -0.0099 | +0.0026 |

## The verdict

- positive arms: **23/24** before, **23/24** after; sign flips: **0**
- the median margin: **+0.0809** before, **+0.0822** after

| interval | arms excluding 0 (standard) | arms excluding 0 (cleaned) |
|---|---|---|

## Reading

This is the sensitivity analysis the negative-set audit demands: if the headline survives with the curated drivers removed from the background, then the contamination of the negative set is not what carries the result.

## Caveats

1. The exploration family (amendment 4): the pre-registered endpoint is the 24-arm standard design, untouched.
2. Removing the NCG genes from the universe changes the walk's normalisation and the fold sizes slightly, so the cleaned margins are not the standard margins with a few negatives deleted; they are a separate, equally valid run.
3. The cleaning uses one catalogue (NCG). A driver curated only in another source (for example the Reactome-curated positives of the open-targets labels) is not removed.
