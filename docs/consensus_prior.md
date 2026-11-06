# Direction 1 -- a consensus prior across networks

**Date:** 2026-11-05
**Script:** `scripts/run_consensus_prior.py`
**Status:** exploratory (work plan Phase S family; see the note in §6)

---

## 1. Why this model

Our own results say the benchmark's one negative cell is a **coverage
failure**, not a noise artefact: Reactome's curated pathway annotations
lag behind the 2020-2024 driver discoveries, and the failure survives
the three controls that equalise the positive set (D-09, D-11, D-12).
A single stale network cannot be repaired by a better learner -- it can
only be out-voted.

So the model is a **consensus prior over the six networks**, built on the
shared universe (D-12's 1,605 genes, so every network is scored on the
same genes and the same 40-repeat x 5-fold split):

1. within one (repeat, fold) evaluation, each network's RWR score is
   converted to a percentile rank over the shared genes;
2. the six ranks are averaged (rank aggregation is scale-free, so a
   network with a different score scale cannot dominate);
3. the consensus is scored against a degree baseline.

## 2. The comparison that matters

Two degree baselines are possible and both are reported:

* **each network's own degree** -- this is the D-12 convention, and it
  is what the released `d12_table.tsv` contains. The members' numbers
  computed this way reproduce that table exactly (e.g. `string_full700
  x intogen2024` +0.1112, `rna_coexp` +0.0518, `reactome` +0.0581), so
  the pipeline is verified against a released table.
* **one common baseline** -- the mean degree of the six networks on the
  shared universe. Comparing the consensus against one network's degree
  while the members are compared against their own degree would inflate
  the consensus by the difference between the baselines, so the
  **common baseline is used for every fair comparison below**.

The common baseline is the stronger of the two (averaging six degree
sequences is a better predictor than any one of them), which is why
the members' numbers against it are smaller, and sometimes negative.

## 3. Q1 -- is the consensus better than its members?

Repeat-level mean delta-AUROC against the common baseline, 95 % t
interval over the 40 repeats:

| label | members median | members best | members worst | **consensus (6)** |
|---|---|---|---|---|
| intogen2024 | +0.0486 | +0.1071 | −0.0838 | **+0.1128** [0.1115, 0.1142] |
| ot_all | +0.0826 | +0.1574 | −0.0950 | **+0.1421** [0.1409, 0.1434] |
| ot_nolit | +0.0845 | +0.1690 | −0.1093 | **+0.1384** [0.1367, 0.1401] |
| **intogen_temporal_new** | **−0.0110** | +0.0306 | −0.0790 | **+0.0358** [0.0301, 0.0416] |

**The consensus beats the median member on all four labels, and the
best member on two** (intogen2024 and the temporal label). On the two
Open Targets labels it lands between the median and the best member:
the consensus is a robustness device, not a universal winner.

## 4. Q3 -- does the consensus repair the Reactome failure?

This is the question the model was built for, and the answer is yes.

On `intogen_temporal_new` the **median single network is negative**
(−0.0110 against the common baseline), three of six are negative, and
Reactome is the worst of all six (−0.0790). The consensus is
**+0.0358 [0.0301, 0.0416]** -- positive and clearly separated from
zero.

The leave-one-network-out consensus variants sharpen it:

| variant | delta-AUROC (temporal) |
|---|---|
| consensus over all six | +0.0358 |
| **consensus without Reactome** | **+0.0514** [0.0457, 0.0570] |
| consensus without intact | +0.0227 |
| consensus without rna_coexp | +0.0262 |
| consensus without string_full700 | +0.0276 |
| consensus without funmap | +0.0331 |
| consensus without string_phys700 | +0.0340 |

**Reactome's inclusion costs the consensus 0.0156 delta-AUROC on the
temporal label** -- it is not merely failing on its own, it is
actively degrading the aggregate. That is the first *quantitative*
statement of the price of curation lag at the level of an ensemble,
and it is exactly what the curator-lag reading predicts.

(The same leave-one-out table is informative on the other labels:
dropping funmap raises the consensus on intogen2024, ot_all and
ot_nolit, so funmap is the weakest member; the consensus is a cheap way
to rank the networks' reliability.)

## 5. Q2 -- is the consensus more stable? (refuted)

**No.** On the shared universe, the consensus's per-gene rank SD is
*higher* than every member's on three of the four labels:

| label | consensus rank SD | member rank SD (min-max) |
|---|---|---|
| intogen2024 | 0.0409 | 0.0210 – 0.0387 |
| ot_all | 0.0364 | 0.0230 – 0.0338 |
| ot_nolit | 0.0379 | 0.0247 – 0.0339 |
| intogen_temporal_new | 0.0528 | 0.0257 – 0.0439 |

So the consensus trades **per-gene variance for accuracy**: averaging
six networks' ranks does not cancel their noise, it adds the
disagreement between networks as variance, while the disagreement is
also where the signal is. This is a bias-variance trade-off that
**AUROC alone cannot see** -- it is the first concrete payoff of the
stability axis built in S1/S2, and it is the honest headline for this
direction: *the consensus is better and less stable*, and a user who
cares about reproducible per-gene rankings should not use it.

## 6. Caveats

1. **The numbers are not comparable with the released 24-cell table.**
   Everything here is on the shared universe (1,605 genes) with a common
   degree baseline; the released table is on each network's own LCC with
   its own degree baseline. The member numbers against their own degree
   reproduce `d12_table.tsv` exactly, which is the verification.
2. **The consensus is a new quantity**, so it cannot be reported under
   the pre-registered endpoint (the pre-registration fixes RWR as the
   learner). It belongs to the second paper, or to the first paper as an
   exploratory arm registered in the exploration family of amendment 4.
3. **The shared universe is small** (1,605 genes, 5–15 % positives per
   label). The temporal label has 5 positives in the whole universe --
   the consensus's +0.0358 is driven by a handful of positives and its
   interval, although tight, does not reflect that.
4. The consensus weights all six networks equally. Weighting by measured
   reliability (the stability atlas of S1/S3) is the obvious next
   variant and is left to the second paper.
5. This is a **single aggregation rule** (mean of percentile ranks).
   Rank fusion (RRF), mean of z-scored ranks and median of ranks were
   not compared; the sensitivity of the conclusion to the rule is
   untested.

## 7. What this suggests for the next steps

* The consensus is the **first model we can justify from our own
  findings**, and it is a good candidate for the second paper's method.
* It needs one robustness check before it goes anywhere: the
  **weighting by reliability** variant, and the aggregation-rule
  sensitivity (caveat 5).
* The leave-one-out consensus is also a useful *diagnostic* for the
  first paper: it ranks the networks by how much they contribute to the
  aggregate, which is a different question from how good each network is
  on its own.
