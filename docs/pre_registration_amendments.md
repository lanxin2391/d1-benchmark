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

Measured: the ratio `(sigma2_within / 5) / var(repeat means)` has
**median 4.21** (range 1.11 – 13.51) across the 24 cells. Treating the 50
folds as 50 independent observations therefore inflates the standard error
of every cell mean by roughly **2x**.

Recomputing the per-cell bootstrap CI with the repeat as the resampling unit:

| unit | mean CI width | arms whose CI excludes 0 |
|---|---|---|
| fold-level (old) | 0.0111 | **23 / 24** |
| repeat-clustered (new) | 0.0057 | **24 / 24** |

The single arm that previously failed — `reactome x intogen_temporal_new`,
mean −0.0112 — gets CI **[−0.0174, −0.0056]** under repeat clustering: the
interval lies entirely **below** zero. This is what the pre-registration
§10 states ("its 95 % CI is entirely below zero, so it is a real null effect,
not just a weak positive"); the fold-level CI contradicted it.

### A2.4 What is NOT changed

The primary endpoint, the learner, the frozen inputs, the decision rule
("H1 holds iff the 95 % CI excludes 0 in **every** arm") and the analysis
plan are untouched. This amendment changes **how** the variance decomposition
and the CI are computed, not **what** is tested.

### A2.5 Honest caveats (must appear in the manuscript)

1. The unit choice was made after seeing that it turns 23/24 into 24/24.
   The *design* argument (repeated CV => cluster on the repeat) is
   a-priori correct, but a coverage simulation must be run to confirm the
   repeat-clustered CI has ~95 % coverage, and **both** CIs must be
   reported. This is scheduled as week-0 work.
2. The pre-registration §9.3 asks for a **gene-level** cluster bootstrap
   ("resampling genes, not gene-arm rows"). That is still outstanding,
   because per-gene scores were not written to disk by the original run.
   It is scheduled for week 2 (`run_arm` must persist per-gene scores).
   The gene-level CI may be wider than the repeat-level CI for the small
   label sets, and the manuscript must state whichever it ends up using.
3. `docs/pre_registration.md` §12 and §14 cross-reference the wrong
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
