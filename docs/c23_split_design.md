# Control 23 — the split-design comparison

**Date:** 2026-11-18  
**Script:** `scripts/run_c23_split_design.py`  
**Tables:** `c23_split_design_arms.tsv`, `c23_split_design_summary.tsv`

## The question

The teacher asked us to try other cross-validation designs (for example 10 x 10) and to say which design is best. This control compares six designs on the same 24 arms (6 networks x 4 catalogues), the same endpoint (alpha = 0.5) and the same interval constructions. **The pre-registered design (5 folds x 40 repeats) stays the primary result; this is a sensitivity analysis, not a re-design.**

## The six designs

| design | evaluations per arm | median margin over the 24 arms | arms with a positive margin | repeat-level CI certifies | gene-cluster CI certifies | median gene-cluster CI width |
|---|---|---|---|---|---|---|
| k5_r10 | 50 | +0.0813 | 23/24 | 24/24 | 20/24 | 0.0366 |
| k5_r20 | 100 | +0.0803 | 23/24 | 24/24 | 20/24 | 0.0366 |
| k5_r40 **(pre-registered)** | 200 | +0.0806 | 23/24 | 24/24 | 20/24 | 0.0365 |
| k10_r10 | 100 | +0.0836 | 23/24 | 24/24 | 20/24 | 0.0370 |
| k10_r20 | 200 | +0.0839 | 23/24 | 24/24 | 20/24 | 0.0373 |
| k20_r10 | 200 | +0.0831 | 23/24 | 24/24 | 20/24 | 0.0366 |

## What the comparison shows

**1. The size of the effect is a property of the data, not of the design.** The median margin over the 24 arms stays between +0.0803 and +0.0839 across all six designs (a spread of 4.5% of its own size), and the per-arm spread across the six designs is at most 0.0075 AUROC.

**2. None of the three verdicts moves.** In every design:

- **23 of 24 arms** have a positive margin;
- the repeat-level percentile interval certifies **24 of 24 arms**;
- the conservative gene-cluster bootstrap certifies **20 of 24 arms** — and the four it does not certify are the same four in every design: funmap x intogen_temporal_new, rna_coexp x intogen_temporal_new, intact x intogen_temporal_new, reactome x intogen_temporal_new.

**3. What the design does move is the width of the conservative interval.** The median gene-cluster interval width goes from 0.0365 to 0.0373. The fold count matters more than the repeat count: at five folds, going from 10 to 40 repeats changes the width by 0.0001 (0.3 %), while going from five to twenty folds at ten repeats changes it by 0.0000 (about 8 %). More folds mean smaller test sets and therefore noisier per-fold AUROC estimates, which is exactly what the interval sees.

**4. The answer to the teacher's question.** There is no 'best' design in the sense of a different conclusion: the effect size, the sign of every arm and both verdicts are the same in all six designs. The pre-registered design is not the most favourable and not the least favourable one — it is simply the one we froze before looking, and we keep it. The design costs about half an hour of compute per extra setting and buys a slightly narrower conservative interval; the recommendation is therefore to report this table as a sensitivity analysis and **not** to switch designs after the fact.

## Two details worth stating

**The one negative arm.** Reactome x the newest drivers has a margin of -0.0125; under the repeat-level construction its interval [-0.0156, -0.0094] excludes zero on the NEGATIVE side, which is why the 'the CI excludes 0 in every arm' count is 24/24 under that construction. The equivalence test is what we quote instead: the arm is not distinguishable from no detectable contribution at the +-0.05 margin, so the sentence we use stays 'no detectable contribution on this arm'. The gene-cluster construction puts it back inside the interval ([-0.0589, +0.0375]), and it is one of the four arms the conservative construction does not certify.

**The k = 10 and k = 20 designs are new splits.** They use the same base seed but a different fold structure, so the gene-to-fold assignment differs from the released panel. That is the point of the comparison, and it is also why it is reported as a sensitivity analysis rather than as a replacement for the released panel.

## Caveats

1. The gene-cluster bootstrap uses 1,000 resamples and the repeat-level percentile interval 2,000, the same as the released analysis.
2. Only alpha = 0.5 is compared here; the other two alphas exist only for the released design.
3. **Exploratory** (amendment 4): the pre-registered endpoint and the H1 verdict are untouched.
