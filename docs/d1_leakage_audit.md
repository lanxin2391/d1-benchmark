# D1 -- how much complex leakage do the current splits contain?

**Date:** 2026-11-05  
**Script:** `scripts/run_d1_leakage_audit.py`  
**Task:** Phase D, task D1 (quantify first, then fix)

## What was measured

For every (arm, repeat, fold): the number of held-out (test) positives that have at least one **complex-linked** gene in the training fold. Two genes are complex-linked when they appear in the same curated complex (ComplexPortal identifiers via UniProtKB; `data/processed/complexes/complex_map.tsv`). The same quantity is computed for the held-out negatives as a control, and the whole thing is normalised by the number of test positives that are in a curated complex at all.

## Result 1 -- the channel is wide open

| quantity | median over the 24 arms | arm range |
|---|---|---|
| test positives with a complex-linked gene in the train fold | **0.345** | 0.167 - 0.427 |
| test NEGATIVES with a complex-linked gene in the train fold | 0.218 | 0.179 - 0.313 |
| **among the test positives that are in a complex**, the fraction with a linked gene in the train fold | **0.957** | 0.000 - 1.000 |
| the same among the test negatives | 0.939 | |
| excess (positives - negatives) | **+0.094** | -0.027 - +0.195 |

**Reading.** For every gene that *can* be leaked, the standard CV leaks it: among the held-out positives that are in a curated complex, **96%** have a complex-linked gene in the training fold. That is not a property of the positives -- the held-out **negatives** are linked at the same rate (93.9%), because whether two genes share a complex does not depend on whether they are drivers. The leakage channel is therefore **wide open but symmetric**: it is opened for the positive and the negative classes alike.

The second quantity is the one that matters for the benchmark's claim. The *enrichment* of positives over negatives is +0.094 on average (-0.027 to +0.195 across arms) and is essentially zero for the four `intogen_temporal_new` arms (their positives are few and mostly not in a curated complex). A symmetric channel with a modest enrichment is the profile of a leak that inflates the *absolute* performance of a propagating learner but is unlikely to inflate its *advantage over the degree baseline* much -- but whether it does is exactly what D4 has to measure, not what D1 can assume.

## Result 2 -- the same holds for the control split families

| split family | test positives linked | negatives linked |
|---|---|---|
| native (the 24 primary arms) | 0.345 | 0.218 |
| shared (D-12) | 0.352 | 0.286 |
| capped (D-11) | 0.306 | 0.221 |

## Decision for Phase D

The leakage exists and is measurable, but it is **symmetric** between the classes. Under the plan's Q3 rule the consequence is decided by D4 (the change in delta-AUROC when the splits become complex-aware), not by D1. D1 does establish that the criticism of Picart-Armada et al. 2019 is *mechanically* applicable to this panel, so D2-D4 are worth doing -- and the prediction to test is that the change in delta-AUROC will be **smaller than the raw leakage numbers suggest**, because the channel favours neither class.

## Caveats

1. The complex map covers 18.3 % of the panel's genes and 33-37 % of the positives -- curated complex data is sparse, and a gene in no curated complex cannot be leaked through this channel. This is a property of the resource, not of the benchmark, and it means the audit measures a **lower bound** on the leakage that a complete complex resource would show.
2. The audit measures leakage in the **v1** splits. The complex-aware splits (D3) must show a near-zero excess on the same metric; that check is part of D3's acceptance criteria.
3. 'Complex-linked' is based on curated complex membership, not on interactions: two proteins that interact strongly without being curated as a complex are not counted here. The reviewer's channel is therefore partly, not fully, covered.

## Files

* `results/tables/d1_leakage_audit.tsv` -- one row per (family, arm, repeat, fold)
* `results/tables/d1_leakage_summary.tsv` -- one row per arm
* `data/processed/complexes/complex_map.tsv` and `data/processed/complexes/provenance.md` -- the complex map and where it came from
