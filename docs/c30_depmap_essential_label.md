# Control 30 — the DepMap essential-gene catalogue as a fifth label

**Date:** 2026-11-25  
**Script:** `scripts/run_c30_depmap_essential_label.py`  
**Tables:** `c30_depmap_essential_label.tsv`, `c30_depmap_essential_summary.tsv`

## The question

The teacher asked whether the network prior helps on the functional side. Control 24 answered that as a ranking question; this control answers it in the benchmark's own terms: the DepMap essential genes become a **label** -- a fifth catalogue at four intervals -- and the whole protocol runs on it (the same networks, the same folds, the same degree baseline, the same endpoint).

The label: the genes essential (CERES < -0.5 and expressed) in at least X % of the 1,361 cancer cell lines, for X in {1, 5, 10, 20}.

## The result

| network | 1 % of the cell lines | 5 % of the cell lines | 10 % of the cell lines | 20 % of the cell lines |
|---|---|---|---|---|
| funmap | +0.1254 | +0.1372 | +0.1395 | +0.1362 |
| rna_coexp | +0.2653 | +0.3357 | +0.3481 | +0.3526 |
| string_phys700 | +0.0893 | +0.1154 | +0.1187 | +0.1218 |
| intact | +0.0204 | +0.0482 | +0.0557 | +0.0630 |
| string_full700 | +0.0878 | +0.1103 | +0.1131 | +0.1110 |
| reactome | +0.1402 | +0.1889 | +0.2007 | +0.2144 |

Across the 24 arms the mean margin against degree is **+0.1236 at the median**, positive in **24 of 24 arms**, with the range +0.0204 to +0.3526.

## Reading

**The H1 verdict is not met on a functional answer key**, which is the answer the teacher asked for: on the curated driver catalogues the prior beats degree in 23 of 24 arms, and on the functional catalogue it does not. This is the same conclusion Control 24 reached from the ranking side (the AUROC of the prior inside the candidate pool is 0.501, i.e. chance), and it is now stated in the protocol's own vocabulary. The two controls agree, which is worth saying explicitly: the prior's contribution is a property of the curated driver labels, not of the network's ability to point at functionally important genes.

The caveats are part of the result: the label is derived from the same expression matrix that defines the essentiality filter, so the RNA co-expression network has a mild dependency on it, and the essential genes are a housekeeping-like set that no curation decision reproduces. The label is therefore a harder and a different kind of answer key than the driver catalogues, and a null result on it is a statement about that key, not about the network.

## The mechanism: is the margin a module-propagation effect?

If the effect is a propagation of the seeds' neighbourhood, the margin should track how module-like the positive set is in that network. The measure is the mean degree of the positives against the background (the ratio):

| network | 1 % of the cell lines | 5 % of the cell lines | 10 % of the cell lines | 20 % of the cell lines |
|---|---|---|---|---|
| funmap | 2.12 | 2.76 | 2.86 | 2.92 |
| rna_coexp | 1.28 | 1.21 | 1.18 | 1.11 |
| string_phys700 | 2.89 | 3.25 | 3.35 | 3.41 |
| intact | 2.12 | 2.22 | 2.21 | 2.21 |
| string_full700 | 2.57 | 3.07 | 3.22 | 3.34 |
| reactome | 1.60 | 1.50 | 1.50 | 1.45 |

Across the 24 arms the Spearman correlation between the margin and the positives' degree ratio is **rho = -0.59**.

## Caveats

1. 10 repeats x 5 folds per arm (the released design uses 40 repeats), because each arm is a full protocol run on a new label.
2. The splits are built fresh for the new label with the protocol's base seed and written to a new directory (`data/processed/splits_depmap/`); the released panels are not touched.
3. **Exploratory** (amendment 4): no pre-registered endpoint is touched.
