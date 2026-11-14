# Control 14 -- the seed-perturbation experiment

**Date:** 2026-11-12  
**Script:** `scripts/run_seed_perturbation.py`  
**Tables:** `seed_perturbation.tsv`, `seed_perturbation_summary.tsv`, `seed_perturbation_predictions.tsv`

## The question

The stability analysis attributes the prior's per-gene instability to the gene's own score variance, and the mechanism to the seed-count law: RWR weights the seeds linearly, so a gene's score is a linear function of how many seeds sit near it. That law predicts that **dropping seeds shrinks the scores together and leaves the ranking nearly intact** -- so the margin over degree should be far more robust to losing seeds than the per-gene score CV would suggest, and the arms with the noisiest scores should lose the most.

## What was run

Every arm, with 100 %, 90 %, 70 % and 50 % of its training-fold seeds (a fresh random subsample per fold, a fixed seed, the same locked learner and the same test folds).

## The result

| seeds kept | median margin | median relative loss | arms losing > 10 % | rho(score CV, loss) | p |
|---|---|---|---|---|---|
| 100 % | +0.0809 | 0 % | 0 | - | - |
| 90 % | +0.0781 | 1.9 % | 1 | -0.037 | 0.862 |
| 70 % | +0.0726 | 5.9 % | 5 | +0.245 | 0.248 |
| 50 % | +0.0673 | 11.7 % | 18 | +0.263 | 0.214 |

## Reading

**The margin is far more robust to losing seeds than the per-gene score CV would suggest.** Dropping 10 % of the seeds costs a median of about 2 % of the margin, and only one of the 24 arms loses more than 10 % of it; halving the seeds costs a median of about 12 %. That is the seed-count law in practice: the scores shrink together, so the ranking -- and with it the margin -- is stable. The answer to the practical question (how many seeds does an arm need?) is therefore not the gene-level noise but the *number* of seeds: the arm with the fewest seeds (Reactome x the newest drivers, about 50) loses 29 % of its margin when half the seeds are dropped, while an arm with about 470 seeds (IntAct x IntOGen) loses 8 %.

**A prediction that did not survive.** We expected the arms whose per-gene score CV is highest to lose the most margin; the Spearman correlation between the score CV and the loss is small and not significant (-0.04 / +0.25 / +0.26 at 90 %, 70 % and 50 % of the seeds). So the gene-level noise does not predict the loss -- the seed count does. We report it as a negative result: the seed-count law explains the *mechanism* of the instability, but the practical cost of losing seeds is set by how many seeds the label has.

## Caveats

1. The subsample is drawn independently per fold, so this is the seed-noise of the estimator, not a change of the training set's size across repeats.
2. The seeds are the label's training positives; the test positives are untouched.
3. Exploratory family (amendment 4): the pre-registered endpoint is untouched.
