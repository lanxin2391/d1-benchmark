# Control 20 -- the label-difficulty regression (#8)

**Date:** 2026-11-12  
**Script:** `scripts/run_label_difficulty_regression.py`  
**Tables:** `label_difficulty_models.tsv`, `label_difficulty_coefs.tsv`, `label_difficulty_labels.tsv`

## The question

Control 6 attributes 34.78 % of the between-arm spread of the margin to a LABEL main effect. Is that five curation philosophies, or one shared difficulty axis?

## The labels

| label | positives | degree AUROC | RWR AUROC | closeness AUROC | mean margin |
|---|---|---|---|---|---|
| ot_nolit | 585 | 0.6434 | 0.7496 | 0.6692 | +0.1062 |
| ot_all | 589 | 0.6643 | 0.7701 | 0.6939 | +0.1058 |
| ncg_known | 693 | 0.6443 | 0.7319 | 0.6782 | +0.0876 |
| intogen2024 | 623 | 0.6479 | 0.7299 | 0.6880 | +0.0821 |
| intogen_temporal_new | 150 | 0.5827 | 0.6036 | 0.6075 | +0.0209 |

## Models (primary 24 arms)

| model | R2 | adj R2 |
|---|---|---|
| M0 network only (deg_auroc) | 0.3132 | 0.1225 |
| M1 network + difficulty (deg_auroc) | 0.8473 | 0.7934 |
| M2 network + label factor (deg_auroc) | 0.9047 | 0.8539 |

## Difficulty metrics (primary 24 arms)

| difficulty metric | beta | p (asymptotic) | R2 gain over network only | share of the label factor | F-test: label factor beyond this metric, p |
|---|---|---|---|---|---|
| deg_auroc | +1.0685 | 0.000 | +0.5341 | 90.3% | 0.029 |
| rwr_auroc | +0.5300 | 0.000 | +0.5784 | 97.8% | 0.382 |
| close_auroc | +0.9180 | 0.000 | +0.4818 | 81.5% | 0.003 |
| n_pos | +0.0002 | 0.000 | +0.5180 | 87.6% | 0.014 |

## Reading

- the same models on all 30 arms (adding the NCG label): deg_auroc beta +1.0656 (label-factor share 90.4%, F-test p 0.039); rwr_auroc beta +0.5292 (label-factor share 97.8%, F-test p 0.526); close_auroc beta +0.9101 (label-factor share 81.5%, F-test p 0.003); n_pos beta +0.0001 (label-factor share 79.0%, F-test p 0.002)
- the recomputed RWR margins match the published Table 1 on all 24 primary arms (max |diff| 4.4e-05)

## Caveats

1. With four or five labels, one difficulty number and the label factor are close to collinear by construction; the permutation test (relabelling the metric across labels) is the honest inference, not the asymptotic p-value.
2. `deg_auroc` is mechanically related to the margin (a label that the degree baseline already ranks well leaves little headroom), so its slope is partly a ceiling effect; `rwr_auroc` and `close_auroc` are reported as difficulty measures that do not subtract the baseline.
3. Exploratory family (amendment 4): the pre-registered endpoint is untouched.
