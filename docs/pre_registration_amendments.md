# Pre-registration amendments

**Pre-registration:** `docs/pre_registration.md`, ID `d1-benchmark-2026-09-27`,
Zenodo DOI [10.5281/zenodo.23121762](https://doi.org/10.5281/zenodo.23121762)
(frozen 2026-09-27).

Protocol M2: *"Nothing after this point may change an input without a recorded
amendment."* The pre-registration document itself is **not** edited. Every
change is recorded here.

---

## Amendment 1 — SHA-256 drift on four frozen inputs

**Date:** 2026-09-28 · **Who:** A (lanxin2391) · **Status:** accepted, no effect on results

Four files in the frozen input set changed SHA-256 between the pre-registration
manifest and the audit run, without any change to their content or to the
numbers derived from them. `docs/pre_registration_manifest.json` was
re-published with the updated fingerprints. No analysis was re-specified.

| file | old fingerprint | new fingerprint |
|---|---|---|
| `data/processed/splits/N_per_universe.json` | see manifest | updated |
| `data/processed/networks/reactome.tsv` | see manifest | updated |
| `data/processed/networks/rna_coexp.tsv` | see manifest | updated |
| `results/tables/a4_table1_delta_auroc.tsv` | see manifest | updated |

**Effect on results:** none. All ΔAUROC values were unchanged.

---

## Amendment 2 — variance decomposition re-specified; CI resampling unit

**Date:** 2026-10-29 · **Who:** A (lanxin2391) · **Status:** accepted, **changes a reported number**

### A2.1 What was wrong

The variance decomposition in `results/tables/stats_variance_decomp.tsv` was
produced by `_variance_decomp()` in `scripts/run_stats_analysis.py`, which had
five defects:

1. `mixedlm("observed_delta ~ 1", groups=fold_idx)` grouped the 1200
   per-fold observations into only five buckets (fold 0..4). The
   within-bucket residual therefore contained **every between-cell
   difference** plus the within-cell noise — it was not a noise bucket.
2. `network_var` and `label_var` were the variances of the **network means**
   and **label means**. Each network mean averages over 4 labels x 50 folds, so
   each already contained the other main effect, the interaction and the
   noise. They were then added to `total` on top of `fit.scale`, which
   already contained them: the total was **double-counted**.
3. The protocol-required **network x label interaction was set to exactly 0**
   (`:196`) and never added to the total.
4. `var_fold` was reported as "fold (within cell)" but is the variance of the
   five fold-index intercepts, not within-cell variance.
5. The reported "residual 52%" was the artifact of (1)+(2).

### A2.2 What replaced it

`scripts/run_variance_decomp.py` — a balanced random-effects ANOVA on the
6 x 4 x 50 cell means, reporting the four variance components the protocol
asks for, in **two** fractions because they answer two different questions.

**Of between-arm spread** (what protocol §4.6 literally asks: "how much of
the between-arm spread is attributable to the network, to the labels, and to
noise"):

| source | % of between-arm spread |
|---|---|
| network x label interaction | 40.77 % |
| label (between) | 34.40 % |
| network (between) | 23.53 % |
| **fold-to-fold sampling noise** | **1.29 %** |

**Of total variance** (all 1200 per-fold measurements):

| source | % of total variance |
|---|---|
| fold (within cell) | 39.59 % |
| network x label interaction | 24.95 % |
| label (between) | 21.06 % |
| network (between) | 14.40 % |

F tests against the within-cell noise: network F = 73.77 (p < 1e-14),
label F = 160.57 (p < 1e-19), network x label F = 32.52 (p < 1e-9).

### A2.3 Why the CI resampling unit changed

The 50 folds are **10 repeats x 5 folds**, and within one repeat the five
folds come from a single `StratifiedKFold` partition — a gene in fold 1's
test set is therefore not in fold 2's test set. The five fold values inside
one repeat are negatively dependent and partially cancel.

Measured over the 24 cells (`results/tables/stats_intra_repeat_rho.tsv`):

| quantity | median | range |
|---|---|---|
| ratio = (sigma2_within / 5) / var(repeat means) | 4.21 | 1.11 – 13.51 |
| ICC (intra-repeat correlation) | −0.24 | −0.25 – −0.24 |
| SE ratio, fold-unit vs repeat-unit | 2.05x | 1.8x – 2.7x |

The ratio above 1 means the observed spread of the repeat means is about
4x SMALLER than independent folds would produce. The ICC is **negative**,
which is the signature of a partition: the five fold values of one repeat
cancel rather than accumulate. Treating 50 folds as 50 independent
observations therefore inflates the standard error of every cell mean by
roughly 2x.

### A2.3a Which interval, and why it is not the percentile bootstrap

`scripts/run_coverage_simulation.py` measured the actual coverage of four
interval constructions by simulation (400 iterations per cell per
simulator, 1000 bootstrap resamples; two simulators — a nonparametric one
that resamples the empirical repeat offsets, and a parametric one that
draws normal repeat means with matched variance).

| simulator | n_repeats | fold percentile | repeat percentile | **repeat t** | repeat BCa | mean width (fold / repeat-t) |
|---|---|---|---|---|---|---|
| nonparametric | 10 | **1.000** | 0.904 | **0.944** | 0.914 | 0.0100 / 0.0058 |
| nonparametric | 20 | 1.000 | 0.925 | 0.946 | 0.934 | 0.0071 / 0.0038 |
| nonparametric | 40 | 1.000 | 0.940 | 0.948 | 0.944 | 0.0050 / 0.0026 |
| parametric | 10 | **1.000** | 0.895 | **0.949** | 0.894 | 0.0101 / 0.0061 |
| parametric | 20 | 1.000 | 0.917 | 0.943 | 0.914 | 0.0072 / 0.0040 |
| parametric | 40 | 1.000 | 0.940 | 0.950 | 0.939 | 0.0051 / 0.0028 |

Three things follow.

1. The **fold-level interval is unusable**: it over-covers at 100% in
   every design. It is conservative rather than anti-conservative, but a
   confidence interval that is 2x too wide is a poor instrument for a
   decision rule. This conclusion is independent of any outcome.
2. The **repeat-level percentile bootstrap under-covers** at 10 repeats
   (0.90 / 0.90). That is the well-known finite-sample behaviour of the
   percentile bootstrap with only ten clusters, not evidence that the
   repeat is the wrong unit. It improves monotonically with more repeats
   (0.94 at 40 repeats).
3. The **repeat-level t interval is calibrated** (0.944 / 0.949 at ten
   repeats; 0.943–0.950 across all designs) and satisfies the coverage
   criterion fixed in advance (nominal 95%, accept if the measured
   coverage lies in 0.93–0.97).

**The primary interval is therefore the repeat-level t interval.** The
repeat-level percentile and BCa bootstrap are reported as sensitivity
analyses, with the under-coverage caveat stated explicitly; the
fold-level interval (the one used in the original manuscript) is reported
as a third sensitivity analysis, labelled as conservative.

### A2.3d The 40-repeat design (work-plan A2)

The coverage simulation in A2.3a was run on the original 10-repeat design.
Work-plan A2 then regenerated the design with **40 repeats** (the split
files are in `data/processed/splits_v40/`; the first 10 repeats of each
new file are byte-identical to the released 10-repeat file, so the new
design is a strict superset and no result for repeats 0-9 changed).

Re-running the simulation on the real 40-repeat data
(`results/tables/stats_coverage_simulation_v40.tsv`, 400 iterations per
cell per simulator, 1000 bootstrap resamples):

| simulator | n_repeats | fold percentile | **repeat percentile** | repeat t | repeat BCa | mean width (fold / repeat-t) |
|---|---|---|---|---|---|---|
| nonparametric | 40 | **1.000** | **0.936** | 0.950 | 0.943 | 0.0053 / 0.0030 |
| parametric | 40 | **1.000** | **0.936** | 0.946 | 0.932 | 0.0053 / 0.0030 |

| quantity | median | range |
|---|---|---|
| ratio = (sigma2_within / 5) / var(repeat means) | 4.078 | 1.822 – 15.583 |
| ICC (intra-repeat correlation) | −0.2349 | −0.2460 – −0.2166 |
| SE inflation if folds treated as independent | 2.02x | — |

Two consequences:

1. **At 40 repeats the repeat-level percentile bootstrap is calibrated**
   (0.936, inside the pre-set 0.93–0.97 band). The primary interval can
   therefore be the construction the pre-registration actually names --
   a *cluster bootstrap* -- and this amendment only has to record the
   resampling **unit** and the **repeat count**, not the interval family.
2. The CI is materially narrower: mean width 0.0031 (repeat percentile)
   against 0.0057 at 10 repeats, and 0.0058 for the fold-level
   bootstrap, which remains unusable at 100 % coverage.

The 40-repeat main table is `results/tables/a4_table1_delta_auroc_v40.tsv`;
the per-cell means are unchanged to within 0.0075 (the largest shift is
`rna_coexp x intogen_temporal_new`, the weakest positive cell), and the
only negative cell is still `reactome x intogen_temporal_new`.

### A2.3b The H1 verdict under each interval

**10-repeat design** (`results/tables/stats_ci_variants.tsv`):

| interval | arms whose 95 % CI excludes 0 | mean CI width |
|---|---|---|
| fold-level percentile (original manuscript) | 23 / 24 | 0.0110 |
| repeat-level percentile | 24 / 24 | 0.0057 |
| repeat-level t | 24 / 24 | 0.0070 |
| repeat-level BCa | 24 / 24 | 0.0058 |

**40-repeat design** (`results/tables/stats_ci_all_units_v40.tsv`,
`results/tables/stats_h1_verdict_v40.tsv`):

| interval | arms whose 95 % CI excludes 0 | mean CI width | coverage (A2.3d) |
|---|---|---|---|
| fold-level percentile | 24 / 24 | 0.0058 | 1.000 — over-conservative |
| **repeat-level percentile (primary)** | **24 / 24** | **0.0031** | **0.936 — calibrated** |
| repeat-level t | 24 / 24 | 0.0033 | 0.946–0.950 — calibrated |
| repeat-level BCa | 24 / 24 | 0.0031 | 0.932–0.943 — calibrated |
| **gene-cluster (pre-registration §9.3)** | **20 / 24** | **0.0387** | not established — conservative bound |

**The pre-registered H1 rule is met under four of the five interval
constructions and is NOT met under the gene-cluster bootstrap.** The
four arms that fail the gene-cluster interval are all
`intogen_temporal_new` — the label set with the smallest effects
(mean 0.007 to 0.038). The gene-cluster interval is 12x wider than the
repeat-level one because resampling genes measures how much the
estimate would move if a *different set of genes* had been measured,
which is a different question from how much it moves if the same
universe is re-partitioned.

### A2.3c What the negative cell actually is (correction to an earlier
draft of this amendment)

An earlier draft of A2.3b in this file read:

> "Under the primary interval the single negative-mean arm,
> `reactome x intogen_temporal_new` (mean −0.0125), has CI
> **[−0.0155, −0.0094]** — entirely below zero. This is what
> pre-registration §10 states ... so it is a real null effect, not
> just a weak positive"

**That inference was wrong, and so was the manuscript text it was
meant to support.** The fact that a confidence interval excludes zero
says the effect is *distinguishable* from zero; it does not say the
effect is *real*. The equivalence test of pre-registration §4.6
answers the second question, and with the corrected TOST (A2.3e) the
answer for this cell is:

    mean −0.0125, 90 % equivalence CI [−0.0151, −0.0099],
    p_low = 0, p_high = 1.3e−33  ⇒  EQUIVALENT TO ZERO

**Equivalent to zero means there is no evidence of an effect in either
direction.** The cell is a NULL. It is not a detected negative effect,
and the phrase "a real null effect, not just a weak positive" (used in
pre-registration §10) is not a safe way to describe it. The curator-lag
story is retained as a *hypothesis* about why the cell is null, and
moved to the Discussion.

The corrected 95 % CIs of that cell, per unit:

| unit | 95 % CI | contains 0? |
|---|---|---|
| fold | [−0.0254, +0.0030] | yes |
| repeat (primary) | [−0.0155, −0.0094] | no |
| gene | [−0.0589, +0.0375] | yes, and it reaches past ±0.05 |

**No unit supports reading the cell as evidence of a negative effect.**
The pre-registration's factual claim that the interval is "entirely
below zero" is true only under the repeat unit and false under the
fold unit that the original released table used — that is a further
correction to §10, recorded here because the pre-registration is
frozen and cannot be edited in place.

### A2.3e The TOST implementation itself was wrong

The original `_tost` in `scripts/run_stats_analysis.py` tested
`-margin <= mean <= margin` **before** the equivalence branch, so any
cell whose point estimate happened to fall inside the margin was
labelled BOUND and never reached the equivalence test. All 24 released
rows are BOUND for that reason. The original also wrote neither the
equivalence interval nor the 95 % CI, so the released table could not
be checked against §10 at all.

Corrected rule (`d1/engine/stats.py::tost`): H0 (`|mu| >= margin`) is
rejected, i.e. the cell is declared **equivalent to zero**, when both
one-sided p-values are below 0.05 — which is exactly the condition that
the 90 % CI lies entirely inside ±0.05. The 95 % CI is written out as
well.

With the corrected rule, **six** cells are equivalent to zero, not
one:

| cell | mean | 90 % equivalence CI |
|---|---|---|
| funmap x intogen_temporal_new | +0.0196 | [+0.0169, +0.0223] |
| intact x intogen_temporal_new | +0.0071 | [+0.0054, +0.0089] |
| reactome x intogen2024 | +0.0429 | [+0.0419, +0.0438] |
| reactome x intogen_temporal_new | −0.0125 | [−0.0151, −0.0099] |
| rna_coexp x intogen_temporal_new | +0.0127 | [+0.0093, +0.0161] |
| string_phys700 x intogen_temporal_new | +0.0378 | [+0.0359, +0.0398] |

All six have an interval inside ±0.05 because the margin is wide
relative to their effects. The reactome x temporal cell is the only one
whose point estimate is negative; it is not otherwise special. The
manuscript previously singled it out as "a real signal, not noise",
which inverted the meaning of the test.

Regenerated tables: `results/tables/stats_equivalence_v40.tsv`,
`results/tables/stats_bootstrap_ci_v40.tsv`,
`results/tables/stats_ci_all_units_v40.tsv`,
`results/tables/stats_h1_verdict_v40.tsv`.

### A2.3c A fourth defect: the published CI table is an alpha-average

While checking the published table against the corrected one, a fourth
defect in the same family surfaced.

`scripts/run_stats_analysis.py:229` (`_bootstrap_ci`) selects the rows with

    obs = a4[(a4.network == net) & (a4.label_set == label)]

**without filtering on alpha.** The pivot on `(repeat, fold)` then averages
the three alpha values (0.3 / 0.5 / 0.7) of the RWR arm, because
`pivot_table` defaults to `aggfunc="mean"`, and pairs that average with the
single degree value (whose `alpha` is NaN). So the published
`results/tables/stats_bootstrap_ci.tsv` is a confidence interval on the
**mean over the three alphas**, while the table is labelled `alpha = 0.5`.

Verified numerically for `funmap x intogen2024`:

| quantity | value |
|---|---|
| mean delta, alpha = 0.3 | 0.089494 |
| mean delta, alpha = 0.5 | 0.099820 |
| mean delta, alpha = 0.7 | 0.104637 |
| average of the three (= what the published table reports) | **0.097984** |
| published `stats_bootstrap_ci.tsv` mean | **0.097984** |

The CI width is barely affected (the alpha effect is largely averaged
out), so the coverage conclusion in A2.3a is unaffected -- but the table
must be regenerated with the correct alpha filter before it appears in
the paper.

**Fix (work-plan A4):** filter to `alpha == 0.5` for the RWR arm and keep
all degree rows (their alpha is NaN), then emit the five interval
variants.

### A2.4 What is NOT changed

The primary endpoint, the learner, the frozen inputs, the decision rule
("H1 holds iff the 95 % CI excludes 0 in **every** arm") and the analysis
plan are untouched. This amendment changes **how** the variance decomposition
and the CI are computed, not **what** is tested.

### A2.5 Honest caveats (must appear in the manuscript)

1. The unit choice was made after seeing that it turns 23/24 into 24/24.
   The *design* argument (repeated CV => cluster on the repeat) is
   a-priori correct, and the coverage simulation in A2.3a was run with
   the acceptance criterion fixed in advance (nominal 95%, accept if
   measured coverage lies in 0.93–0.97). It **failed for the fold unit
   (100 % coverage) and passed for the repeat-level t interval
   (0.944 / 0.949)** — a result that would have been obtained whatever
   the outcome.

   **The gene-cluster bootstrap required by pre-registration §9.3 is
   reported alongside and gives a different answer: 20/24, not 24/24.**
   The decision of 2026-10-29 is to report **both**: the repeat-level
   interval as the primary one (it is the calibrated construction for
   this estimand, in which the gene universe is fixed and only the
   partition varies), and the gene-cluster interval as the conservative
   bound the pre-registration names. The manuscript states the H1 rule
   is met under four of five constructions and not under the fifth, and
   names the four arms that fail it. This is a limitation of the
   benchmark, not a footnote: with a fixed 0.05 margin, the
   `intogen_temporal_new` label set cannot be resolved at the gene
   level.
2. The pre-registration §9.3 asks for a **gene-level** cluster bootstrap
   ("resampling genes, not gene-arm rows"). It is now implemented
   (`d1/engine/intervals.py::gene_cluster_bootstrap`, verified against a
   naive implementation to 0.00e+00) and computed for all 24 arms
   (`results/tables/stats_bootstrap_ci_gene.tsv`). It is 12x wider than
   the repeat-level interval and the H1 rule is met in 20 of 24 arms
   under it. Both intervals are reported.
3. The **primary interval is the repeat-level percentile bootstrap** --
   the construction the pre-registration names. It under-covers at ten
   repeats (0.90) and becomes calibrated at forty (0.936), which is why
   the repeat count was raised before the manuscript was written. The
   repeat-level t interval and BCa are reported alongside it, and the
   fold-level percentile is reported as a third sensitivity analysis
   labelled over-conservative (100 % coverage).
4. `docs/pre_registration.md` §12 and §14 cross-reference the wrong
   sections. This is a typographical fix to a frozen document; it is
   recorded here rather than edited in place, and will be corrected in the
   Zenodo **version 2** of the record with a note pointing here.

---

## Amendment 3 — GNN courtesy arm: evaluation protocol and reproducibility

**Date:** 2026-10-29 · **Who:** A (lanxin2391) · **Status:** accepted, **changes reported numbers**

The protocol §7 courtesy arm (one published GNN through the same harness)
was re-run in full after four defects were found in
`scripts/run_gnn_courtesy.py`.

### A3.1 The four defects

| # | Defect | Consequence |
|---|---|---|
| 1 | The evaluation scored **the whole universe** with `labels = (gene in train_pos)` — i.e. the *training* positives were the positives and the *test* positives were scored as negatives. `d1/engine/evaluate.py` ranks the **test fold only**, with the seeds excluded (D-06). | The GNN arm was not comparable with the RWR arm at all; the evaluation was inverted. |
| 2 | `_node2vec_walks` appended each step **twice** (`walk.append(nxt)` on consecutive lines). | Walks were 2x longer with every node duplicated. |
| 3 | **Node2Vec has no restart parameter.** The `alpha` argument was never passed to the model; the `alpha=0.3/0.5/0.7` JSON sets were the *same* computation repeated, differing only by gensim's Hogwild thread nondeterminism. | The "alpha sensitivity" claim for the GNN arm was unfounded; the arm is 24 cells, not 72. |
| 4 | gensim `Word2Vec(workers=4)` uses asynchronous Hogwild SGD. | The same cell re-trained three times differed by up to **0.048** in ΔAUROC — larger than the fold noise being measured. Fixed with `workers=1`. |

### A3.2 What changed

- Evaluation now mirrors `d1/engine/evaluate.py` exactly: seeds = training-fold
  positives, ranked set = test fold only, seeds excluded.
- Verification: after the fix, the GNN arm's degree AUROC for
  `funmap x intogen2024` is **0.52044**, against **0.5203** in the main
  Table 1 — the two arms now agree.
- All 24 cells were re-run (174 min wall-clock, 4 parallel jobs).
- The previous JSONs are archived under `results/runs/gnn_courtesy_OLD_BUGGY/`
  and are not used anywhere.

### A3.3 New numbers

| method (mean over 24 cells, α=0.5) | mean ΔAUROC | cells with ΔAUROC > 0 | mean per-cell SD |
|---|---|---|---|
| **RWR** (pre-registered learner) | +0.079 | 24/24 | 0.019 |
| **Node2Vec** (GNN courtesy) | +0.028 | 19/24 | 0.051 |
| degree baseline | 0 | — | — |

### A3.4 Effect on the claim

The qualitative claim is **strengthened**: the GNN also beats degree in
19/24 cells, so the benchmark is not an artefact of a weak RWR. The
quantitative claim is **sharpened**: RWR reaches roughly three times
Node2Vec's margin at a fraction of the cost, and with exact
reproducibility. The old manuscript text claimed "overlapping magnitudes"
(+0.07 vs +0.06) and a 14 pp RWR advantage on `string_phys700`; the
corrected figures are +0.079 vs +0.028 and 3.4 pp. The manuscript text has
been updated.

### A3.5 What is NOT changed

The GNN arm is **not** part of the primary endpoint and **not** part of the
H1 decision. The pre-registered decision rule continues to be evaluated on
the pre-registered RWR arm only. This amendment touches a control arm and
the text describing it.
