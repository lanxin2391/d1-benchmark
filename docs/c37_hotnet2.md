# Control 37 — HotNet2 head-to-head

**Date:** 2026-11-26  
**Script:** `scripts/run_c37_hotnet2.py`  
**Tables:** `c37_hotnet2.tsv`

## The design

The mutation matrix is the **MC3 TCGA MAF** (the pan-cancer somatic mutation aggregation) and the network is one of our six, so the mutation evidence is independent of the labels: this arm is a genuine mutation-first strategy, not the label re-read. The evaluation is the same endpoint as everything else: the AUROC of the HotNet2 score for the label's positives against the degree baseline on the same test genes.

The released HotNet2 is Python 2 code; the pipeline runs a 2to3-ported copy.

## The result

| network | label | the HotNet2 AUROC | the degree AUROC | the margin |
|---|---|---|---|---|

## Reading

## Caveats

1. HotNet2's score depends on the mutation matrix and the network, not on the label, so one run per network serves all four labels; the test genes are the label's positives against the network's other genes.
2. The permutation count for the significance test is reduced from the paper's 100 (the runtime), which makes the test less powerful and the scores noisier than the published setting.
3. The C and Fortran extensions are not compiled, so the similarity matrix falls back to the NumPy path (the same result, slower).
4. **Exploratory** (amendment 4): no pre-registered endpoint is touched.
