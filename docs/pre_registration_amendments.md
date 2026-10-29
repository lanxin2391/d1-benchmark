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

### A2.3b The H1 verdict under each interval

| interval | arms whose 95 % CI excludes 0 | mean CI width |
|---|---|---|
| fold-level percentile (original manuscript) | 23 / 24 | 0.0110 |
| repeat-level percentile | 24 / 24 | 0.0057 |
| **repeat-level t (primary)** | **24 / 24** | 0.0070 |
| repeat-level BCa | 24 / 24 | 0.0058 |

Under the primary interval the single previously-failing arm,
`reactome x intogen_temporal_new` (mean −0.0112), has CI
**[−0.0186, −0.0037]** — entirely below zero. This is what
pre-registration §10 states ("its 95 % CI is entirely below zero, so it
is a real null effect, not just a weak positive"); the fold-level
interval contradicted it.

The verdict 24/24 is therefore **not** an artefact of the interval
choice: it holds under all three repeat-level constructions, and the
one that is calibrated by simulation is the one reported.
(`results/tables/stats_ci_variants.tsv`.)

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
   the outcome. All three repeat-level constructions give 24/24, so the
   verdict is also robust to the interval family.
2. The pre-registration §9.3 asks for a **gene-level** cluster bootstrap
   ("resampling genes, not gene-arm rows"). That is still outstanding,
   because per-gene scores were not written to disk by the original run.
   It is scheduled for week 2 (`run_arm` must persist per-gene scores).
   The gene-level CI may be wider than the repeat-level CI for the small
   label sets, and the manuscript must state whichever it ends up using.
3. The **repeat-level percentile bootstrap is reported only as a
   sensitivity analysis** because it under-covers at ten repeats
   (0.90 against a 0.95 nominal). The primary interval is the
   repeat-level **t** interval. If the number of repeats is raised
   (work-plan B1, 20–40 repeats), the percentile bootstrap becomes
   calibrated as well (0.94 at 40 repeats) and the two will be reported
   side by side.
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
