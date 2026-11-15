# F3b -- the DepMap validation three ways: prior, degree, random

**Date:** 2026-11-12  
**Script:** `scripts/run_depmap_degree_baseline.py`  
**Tables:** `results/tables/depmap_three_way.tsv`, `depmap_prior_vs_degree.tsv`

## Why this control was missing

The released validation (F3) asked whether the genes the prior puts at the top of the ranking are more likely to be required in cancer cell lines. It reported no enrichment beyond chance for the top 30 candidates and a clear one for the top 100 -- but against a nominal 5 % level, and without the two comparisons that matter: **the degree baseline** (the strongest simple baseline, and the one our own stability analysis says is strong) and **random candidates** (the actual base rate, which is not 5 % when a cell line has a few hundred essential genes among ~18,000).

All three rankings are computed on the same test folds, the same universe (the prior's candidate pool, so the prior's column reproduces the released numbers), the same K and the same cell lines.

### K = 0.01

- the prior: **2.8%** of the cell lines enriched (median p 0.60)
- the degree baseline: **65.1%** (median p 0.05)
- random candidates: **2.5%** (median p 0.68)
- the prior beats random in **13/24** arms; the degree beats the prior in **21/24**

### K = 30

- the prior: **5.8%** of the cell lines enriched (median p 0.52)
- the degree baseline: **69.2%** (median p 0.04)
- random candidates: **3.1%** (median p 0.65)
- the prior beats random in **18/24** arms; the degree beats the prior in **20/24**

### K = 100

- the prior: **47.6%** of the cell lines enriched (median p 0.06)
- the degree baseline: **97.8%** (median p 0.00)
- random candidates: **6.1%** (median p 0.53)
- the prior beats random in **23/24** arms; the degree beats the prior in **21/24**

## Reading

This is the control the earlier statement was missing, and it cuts both ways. The prior is better than random at the larger candidate sets, but **the degree baseline is clearly better than the prior at every scale**. Hubness predicts functional essentiality well -- that is a real biological signal, not an artefact of our benchmark -- and a diffusion prior adds nothing to it on this task.

So the scope of our claim narrows further and it is worth stating plainly: on the curated driver labels the prior adds to degree (23 of 24 arms); on a cell-line functional task degree is the better predictor.

## Caveats

1. The candidates are the held-out negatives of each arm, so the test set is identical in the three columns.
2. The essentiality definition is a proxy for dependency, and it correlates with expression and study intensity.
3. Exploratory family (amendment 4); the pre-registered endpoint is untouched.
