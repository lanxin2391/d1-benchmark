# Control 26 — the label-permutation audit (the ultimate leakage check)

**Date:** 2026-11-25  
**Script:** `scripts/run_c26_label_permutation.py`  
**Tables:** `c26_label_permutation.tsv`, `c26_label_permutation_summary.tsv`

## The question

If the positive labels carried no information about the genes, would the pipeline still report a positive margin? The label vector is permuted inside every (repeat, fold) group -- the positive count and the gene set stay exactly as they are -- and the **whole pipeline is re-run**: the permuted labels decide the seeds, the seeds decide the scores, and the scores are evaluated against the permuted test labels. The comparison statistic is the arm-level mean delta-AUROC, the same statistic the headline uses.

This is the one leakage test that was still missing from the list, and it exercises the real code path end to end. A leak anywhere in that chain would show up as a margin above zero.

## The result

| network | label | the observed mean | the null mean | the null SD | the null 95 % interval | the empirical p |
|---|---|---|---|---|---|---|
| funmap | intogen2024 | +0.1001 | -0.0002 | 0.0096 | [-0.0174, +0.0166] | 0.000 |
| funmap | intogen_temporal_new | +0.0196 | +0.0053 | 0.0232 | [-0.0275, +0.0500] | 0.200 |
| funmap | ot_all | +0.0853 | +0.0032 | 0.0130 | [-0.0149, +0.0288] | 0.000 |
| funmap | ot_nolit | +0.0780 | +0.0026 | 0.0148 | [-0.0177, +0.0326] | 0.000 |
| intact | intogen2024 | +0.0651 | -0.0003 | 0.0043 | [-0.0077, +0.0069] | 0.000 |
| intact | intogen_temporal_new | +0.0071 | +0.0037 | 0.0120 | [-0.0144, +0.0278] | 0.400 |
| intact | ot_all | +0.0842 | +0.0037 | 0.0083 | [-0.0085, +0.0202] | 0.000 |
| intact | ot_nolit | +0.0784 | -0.0010 | 0.0060 | [-0.0113, +0.0083] | 0.000 |
| reactome | intogen2024 | +0.0429 | -0.0028 | 0.0114 | [-0.0239, +0.0155] | 0.000 |
| reactome | intogen_temporal_new | -0.0125 | +0.0106 | 0.0348 | [-0.0407, +0.0797] | 0.750 |
| reactome | ot_all | +0.1103 | +0.0005 | 0.0142 | [-0.0199, +0.0256] | 0.000 |
| reactome | ot_nolit | +0.0996 | +0.0053 | 0.0172 | [-0.0245, +0.0365] | 0.000 |
| rna_coexp | intogen2024 | +0.0521 | +0.0007 | 0.0162 | [-0.0200, +0.0336] | 0.000 |
| rna_coexp | intogen_temporal_new | +0.0127 | -0.0020 | 0.0231 | [-0.0430, +0.0343] | 0.250 |
| rna_coexp | ot_all | +0.0767 | -0.0016 | 0.0117 | [-0.0237, +0.0165] | 0.000 |
| rna_coexp | ot_nolit | +0.0834 | -0.0022 | 0.0092 | [-0.0160, +0.0143] | 0.000 |
| string_full700 | intogen2024 | +0.1133 | +0.0001 | 0.0096 | [-0.0189, +0.0164] | 0.000 |
| string_full700 | intogen_temporal_new | +0.0604 | +0.0009 | 0.0149 | [-0.0277, +0.0250] | 0.000 |
| string_full700 | ot_all | +0.1201 | +0.0030 | 0.0090 | [-0.0107, +0.0151] | 0.000 |
| string_full700 | ot_nolit | +0.1267 | +0.0021 | 0.0091 | [-0.0111, +0.0169] | 0.000 |
| string_phys700 | intogen2024 | +0.1189 | +0.0012 | 0.0112 | [-0.0184, +0.0171] | 0.000 |
| string_phys700 | intogen_temporal_new | +0.0378 | +0.0007 | 0.0274 | [-0.0415, +0.0507] | 0.100 |
| string_phys700 | ot_all | +0.1582 | -0.0013 | 0.0083 | [-0.0162, +0.0128] | 0.000 |
| string_phys700 | ot_nolit | +0.1711 | +0.0007 | 0.0133 | [-0.0200, +0.0210] | 0.000 |

## Reading

**The null is centred on zero in every one of the 24 arms** (the null means run from -0.0028 to +0.0106), so the pipeline does not manufacture a margin from random labels: the positive result is a property of the labels, not of the machinery.

The observed margin is beyond the null's 95 % interval in **19 of the 24 arms** (in those arms none of the 20 permutations reached it, so p < 1/20). In the other **5 arms** the observed margin lies inside the null's range, and **all of them are the newest-driver (temporal) label** -- the arm family where the margin is smallest and negative in one case. Read with F.25 that is the expected shape: the newest-driver labels are where the network and the catalogue are least in step, and the permutation null cannot separate a real but small margin from a random one.

The table also calibrates what a margin is worth: the null's SD runs from 0.0043 to 0.0348, so **a margin of +0.01 on its own is not evidence**. The smallest observed margin that clears its null is +0.0429 (reactome x intogen2024), about 4.0 null SDs above its null mean.

Two things this does **not** say. It does not say that the biological claim is free of bias -- the bait-matched null (F.12) and the publication-count control (F.14) do that work. And it does not license the margin's size: the permutation null says the margin is not an artefact of the pipeline, not that it is large.

## Caveats

1. 20 permutations per arm, so the empirical p is a multiple of 1/20; a p of 0 means that none of the 20 permutations reached the observed value, that is p < 1/20.
2. The permutation is at the gene level: the label of a gene is permuted among the genes of the split file and applied to every row of that gene, so the positive count is unchanged and the seeds of a fold are permuted with the test labels. The null is 'the label is random', not 'the evaluation is random'.
3. The statistic is the arm-level mean delta-AUROC over the 200 folds, the same statistic as the headline.
4. **Exploratory** (amendment 4): no pre-registered endpoint is touched.
