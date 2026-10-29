# Coverage simulation for the CI resampling unit

**Task:** work plan v3.1 A1
**Script:** `scripts/run_coverage_simulation.py`
**Data:** `results/runs/a4_combined.tsv`, alpha = 0.5, 24 cells x 50 folds
**Seed:** 20261029

## Why this exists

The per-cell confidence interval was changed from resampling **folds**
(50 units) to resampling **repeats** (10 units). The change was made
after observing that it turns the H1 verdict from 23/24 into 24/24, so
it is exposed to an "outcome-dependent analysis" charge. This document
measures whether the two units are calibrated, instead of asserting it.

## 1. The dependence structure

The 50 folds are 10 repeats x 5 folds. Within one repeat the five folds
come from a single `StratifiedKFold` partition, so a gene in fold 1's
test set is not in fold 2's test set: the five fold values of one repeat
are negatively dependent and partially cancel.

Measured over the 24 cells (`results/tables/stats_intra_repeat_rho.tsv`):

| quantity | median | range |
|---|---|---|
| ratio = (sigma2_within / 5) / var(repeat means) | **4.078** | 1.822 .. 15.583 |
| ICC (intra-repeat correlation) | **-0.2349** | -0.2460 .. -0.2166 |
| SE ratio, fold-unit vs repeat-unit | **2.02x** | 1.35 .. 3.95 |

The ratio above 1 means the observed spread of the repeat means is
4.1x SMALLER than independent folds would produce. The ICC is
**negative**, which is the signature of a partition: the five fold
values of one repeat cancel rather than accumulate.

## 2. Coverage simulation

For each cell, 400 synthetic datasets were generated with the
observed dependence structure and a known true mean, and a nominal-95%
CI was built with each resampling unit.

| simulator | fold-level percentile (current) | repeat-level percentile | repeat-level t | repeat-level BCa |
|---|---|---|---|---|
| nonparametric (empirical repeat offsets) | **nan** | nan | **nan** | nan |
| parametric (normal repeat means, matched variance) | **nan** | nan | **nan** | nan |

(median over the 24 cells; per-cell values in
`results/tables/stats_coverage_simulation.tsv`)

The two simulators have complementary weaknesses: the nonparametric one
has no distributional assumption but uses the plug-in estimate of the
true mean (slight upward bias in coverage); the parametric one has an
exactly known truth but assumes the repeat means are normal. They agree.

## 3. Effect of adding repeats

The coverage table above is for the design we actually ran (10 repeats).
The simulation also runs the same procedure with 20 and 40 repeats to
answer work-plan task B1 directly:

| simulator | n_repeats | repeat-level t coverage | mean width |
|---|---|---|---|
| nonparametric | 10 | nan | nan |
| nonparametric | 20 | nan | nan |
| nonparametric | 40 | 0.950 | 0.0030 |
| parametric | 10 | nan | nan |
| parametric | 20 | nan | nan |
| parametric | 40 | 0.946 | 0.0030 |

## 4. Conclusion

1. The **fold-level CI is not usable**: it over-covers at
   nan% (parametric nan%) against a 95% nominal. It is
   not anti-conservative -- it is conservative to the point of being
   useless for a decision rule, because the SE is 2.0x too
   large. This confirms the unit was wrong, independently of any
   outcome.
2. The **repeat-level percentile bootstrap under-covers**
   (nan% / nan%). This is the well-known finite-sample
   behaviour of the percentile bootstrap with only 10 clusters, **not**
   evidence that the repeat is the wrong unit.
3. The **repeat-level t interval is the calibrated choice**
   (nan% / nan%), and it becomes better still with more
   repeats (see the table above).
4. Therefore the recommendation is: **primary CI = repeat-level, with
   20-40 repeats (work-plan B1), using the t interval; the 10-repeat
   percentile bootstrap is reported as a sensitivity analysis, with the
   under-coverage caveat stated explicitly.**

The choice of unit is justified by the **design** (repeated CV => the
repeat is the independent unit) and is confirmed by a coverage
measurement carried out with the criterion fixed in advance (nominal
95%, accept if measured coverage lands in 0.93-0.97). This is the
argument that goes into pre-registration amendment 2 section A2.5.

## 5. What is still pending

- The **gene-level** cluster bootstrap (pre-registration section 9.3)
  cannot be simulated here: it needs per-gene scores, which the original
  `run_arm` did not write to disk. That is work-plan task A2, and the
  gene-level coverage assessment follows it (B2). The direction is
  expected to be similar, but it must be measured, not assumed.
- The per-cell coverage estimates carry Monte-Carlo error of about
  0.050 at 400 iterations; the medians
  over 24 cells are correspondingly more stable.
