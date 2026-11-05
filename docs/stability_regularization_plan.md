# Stability-Regularised Scoring (SRS) — pre-registered analysis plan

**Status:** pre-registered (written and committed **before** any result of this analysis exists).
**Branch:** `feature/stability-regularized-scoring`
**Date:** 2026-11-05 · **Author:** analysis agent, on behalf of A
**Family:** **EXPLORATORY** (see §0.3). This analysis does **not** enter the H1 decision.

---

## 0. Why this analysis, and what it is not

### 0.1 The question

S2 (stability attribution, upstream repo) reported that the instability of a network prior is
driven by *crowding*: genes whose prior scores are packed into a narrow band swap ranks whenever
the seed set changes. If that mechanism is right, the instability should be **removable**: shrink a
gene's score toward the median of its own score neighbourhood, and its rank should stop moving.

This analysis therefore asks a single, falsifiable question:

> **Can the instability of a network prior be reduced by a stability-regularised scorer, without
> paying for it in discrimination?**

The deliverable is not "a better model". It is a **trade-off curve**: how much ranking stability a
network prior can be bought with, and what that costs in ΔAUROC. That curve is, as far as we can
tell, not in the published record — Guala & Sonnhammer (2017) report *which* method wins, Hill et
al. (2019) report *degeneracy* under random input, Wright et al. (2025) report *robust Z* against a
degree-matched null, and Picart-Armada et al. (2019) report *CV-scheme* effects, but none of them
reports stability of a prior as a quantity that can be traded against discrimination.

### 0.2 Relation to the locked learner (protocol §4.3)

Protocol §4.3 fixes the learner (RWR, α = 0.5 primary) and states that it "is deliberately simple
and deliberately never tuned. The study measures what the *prior* contributes; any learner
improvement would confound that."

SRS is a **learner-side** operation, so:

* it is **not** applied to the confirmatory factorial. The 24 primary arms remain untouched.
* it is reported as an **exploratory arm** in its own section, with its own ΔAUROC numbers that are
  never mixed with the primary table.
* it is registered in the pre-registration amendment as exploratory (§0.3).

This is the same handling as Phase E's hub-capped RWR (work_plan v3.1, E3) and Phase H's own
models (H6): they are controls/second-paper material, not confirmatory results.

### 0.3 Family assignment (work_plan v3.1 E4 / R19)

| item | family | enters H1? |
|---|---|---|
| SRS scorer, all λ/k settings | exploratory | **no** |
| Stability metrics (rank SD, score CV, crowding) | exploratory | **no** |
| Random-smoothing and global-median controls | exploratory | **no** |
| 40-repeat split regeneration (`splits_v40/`) | infrastructure; reproduces the design already adopted upstream (work_plan A2) | primary numbers re-derived, not re-decided |
| Existing 24-arm ΔAUROC | **confirmatory** | yes — unchanged in substance |

---

## 1. Pre-registered predictions

Stated before any SRS number exists. Each is falsifiable and each has a stated decision rule.

| # | Prediction | Decision rule (fixed now) |
|---|---|---|
| **P1** | Median per-gene **rank SD** decreases as λ increases. | The λ-monotone decrease must hold in **≥ 24/24 arms** at the primary setting (k = 20, confidence-weighted, transductive neighbourhood). If it holds in fewer than 24/24, P1 is **not** supported and is reported as such. |
| **P2** | **ΔAUROC does not fall.** | Median across arms of \|Δ(ΔAUROC)\| between λ = 0 and the λ chosen by the stability objective **≤ 0.010 AUROC** (≈ 1 AUROC point). Also reported as the full distribution per arm. |
| **P3** | A **monotone trade-off curve** exists: stability improves with λ while ΔAUROC degrades slowly. | Quantified as the *price of stability*: Δ(ΔAUROC) per 10 % reduction of median rank SD, per arm; summarised as median and IQR across 24 arms. The curve must be reported whether or not it is monotone; a non-monotone curve is a negative result and will be reported as such. |
| **P4** | **Specificity**: the stability gain of score-neighbourhood smoothing exceeds a *matched-strength random* smoothing control. | For each arm, the observed reduction in median rank SD must exceed the **95th percentile** of the matched random-smoothing null (30 replicates, same k, same λ) in **≥ 24/24 arms**. |
| **P5** | The gain is **not** an artefact of flattening the scores. | (a) Rank SD and AUROC are invariant to any monotone transformation of the score vector — verified numerically to machine precision on real folds; therefore any observed change is a *ranking* change, not a scale change. (b) The AUROC is preserved (P2). (c) The score-space dispersion (SD, IQR) after smoothing is reported, and the fraction of genes whose rank changes at all is reported. |
| **P6** | The effect is **crowding-specific**, not uniform. | The per-gene reduction in rank SD must be larger for crowded genes than for isolated genes (paired test over genes within each arm; reported as Spearman ρ between crowding and Δrank-SD, and as a mixed-model coefficient with arm random effects). |
| **P7** | The gain **survives leaving the fold** (is not an artefact of transduction). | The *inductive* variant (neighbourhood built from training-fold genes only) reproduces the direction of P1/P4 with an effect size not smaller than half of the transductive one. |

**Declared circularity risk (the reason P4/P5/P7 exist).** SRS optimises a stability metric with a
stability metric. Of course smoothing lowers the metric it optimises. The protections are:

1. **P4** — a matched-strength null that has the same operator and the same strength but no
   crowding information. If score-neighbourhood smoothing were merely "shrinking the score vector",
   the random control would shrink it just as much.
2. **P5** — the invariance argument plus the AUROC check and the dispersion report.
3. **P7** — an inductive variant whose neighbourhood contains no test-fold gene.
4. **Out-of-sample confidence** — the per-gene confidence weights used by the operator are estimated
   **leave-one-repeat-out** (see §3.3), so no test repeat contributes to its own weights.
5. **The stability metric is never the only metric.** Every stability number is reported next to
   ΔAUROC, ΔAUPRC, pAUC@FPR≤0.02 and MedRR.

---

## 2. Data and design

### 2.1 Arms

The full primary factorial: **6 networks × 4 label sets = 24 arms**, native universe, α = 0.5
(primary) with α ∈ {0.3, 0.7} as sensitivity. Networks: `funmap`, `rna_coexp`, `string_phys700`,
`intact`, `string_full700`, `reactome`. Labels: `intogen2024`, `ot_all`, `ot_nolit`,
`intogen_temporal_new`.

### 2.2 Splits

**40 repeats × 5 folds** (work_plan v3.1 A2), generated into a **new** directory
`data/processed/splits_v40/`; the released v1 splits are never modified (work_plan risk R3).

**Mandatory verification before use:** for every arm, rows belonging to repeats 0–9 of the v40 file
must be byte-identical to the corresponding rows of the released v1 file (`make_splits` uses
`random_state = base_seed + r`, so repeats 0–9 are shared by construction; this is checked, not
assumed).

### 2.3 Per-gene score dump

For every arm and every α, the RWR score of every gene in every test fold is written to
`results/scores_v40/<network>__<label>__alpha<a>.npz` (float64), together with the gene index,
the repeat/fold index and the label vector. The directory is derived data (regenerable in minutes)
and is git-ignored; a SHA-256 manifest of every file is committed.

**Mandatory verification:** the AUROC recomputed from the dumped scores must equal the AUROC
computed by the harness to ≤ 1e-12 (A3 in work_plan v3.1 established this check; it is repeated
here).

### 2.4 Why 40 repeats

Per-gene stability is estimated from the across-repeat variability of a gene's rank. With 10
repeats each gene contributes 10 rank observations; with 40 it contributes 40. The stability
metrics are the *dependent variables* of this analysis, so their own precision matters. 40 repeats
is the design already adopted upstream (work_plan A2), so results are directly comparable.

---

## 3. Method

### 3.1 The scorer

Notation for one evaluated fold. Let `G` be the fold's test genes (|G| = n), `s_g` their RWR
scores, `y_g` their labels.

**Step 1 — score-space neighbourhood.** For each gene `g ∈ G`, `N(g)` is the set of the `k` genes
in `G \ {g}` with the smallest `|s_j − s_g|`. Ties (equal scores, or equal distance) are broken by
ascending gene symbol, matching the D-07 determinism convention. `k` is a hyper-parameter
(default 20; grid {5, 20, 50}).

**Step 2 — confidence-weighted local median.** Let `c_j ∈ (0, 1]` be the *confidence* of gene `j`
(§3.3). The local median `m_g` of `g` is the minimiser of `Σ_{j ∈ N(g) ∪ {g}} c_j · |s_j − v|`
over `v` (a weighted median; for uniform weights it is the ordinary median).

**Step 3 — shrinkage.** `s'_g = s_g + λ_g · (m_g − s_g)`.

* **Variant A (primary, global λ):** `λ_g = λ`, `λ ∈ {0, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0}`.
* **Variant B (crowding-adaptive λ):** `λ_g = λ · d_g / max_{j∈G} d_j`, where `d_g` is the local
  score density at `s_g` (Gaussian kernel, Scott's bandwidth within the fold). Crowded genes are
  shrunk fully; isolated genes are barely touched. This is the "denoise only where it is noisy"
  version.

**Step 4 — evaluation.** `s'` is ranked with the D-07 tie-break (score descending, symbol ascending)
and scored with the *unchanged* metrics: AUROC (primary), AUPRC, pAUC@FPR≤0.02 (Guala &
Sonnhammer's top-of-list endpoint), MedRR (median rank ratio of true positives), p@k ∈ {50, 100,
500}. The degree-only baseline is **not** smoothed.

The operator reads only the score vector of the fold. It never reads `y`.

### 3.2 Neighbourhood source: transductive (primary) vs inductive

* **Transductive (primary).** `N(g)` is formed inside the evaluated fold. This is the literal
  reading of "shrink toward the median of its score neighbourhood". It is transductive: the
  neighbourhood contains other test genes.
* **Inductive.** `N(g)` is formed among the **training-fold genes only** (excluding the training
  positives, which are the seeds). No test gene enters any neighbourhood.

Both are reported. P7 states the decision rule for the inductive variant.

### 3.3 Confidence weights (leave-one-repeat-out)

For gene `g` in arm `a`, let `CV_g` be the coefficient of variation of `g`'s RWR score across the
repeats in which `g` was in the test set, estimated **on the repeats other than the repeat being
evaluated**:

```
CV_g^(−r) = sd_{r' ≠ r}( s_g^{r'} ) / max( mean_{r' ≠ r}( s_g^{r'} ), eps )
c_g = 1 / (1 + CV_g^2)
```

`c_g ∈ (0, 1]`; a gene whose score is stable relative to its level gets weight 1, a gene whose
score fluctuates gets a low weight. Genes with no measurement in the other repeats get the arm's
median confidence.

Two weighting schemes are run:

* **W1 (primary):** `c_g` as above.
* **W0 (ablation):** uniform weights (`c_g ≡ 1`) — the plain local median.

### 3.4 Stability metrics

For gene `g` and repeat `r` in which `g` was in the test fold, let `ρ_{g,r}` be the percentile rank
of `g` inside that fold (0 = best, 1 = worst). Then

* **primary:** `SDrank_g = sd_r( ρ_{g,r} )` — per-gene rank SD;
* **secondary:** `CVscore_g = sd_r(s_g^r) / mean_r(s_g^r)`;
* **arm level:** median and mean of `SDrank_g` over genes; fraction of genes with `SDrank_g` above
  the arm's 75th percentile.

### 3.5 Crowding

For gene `g` in a fold, `d_g` = Gaussian kernel density of the fold's score distribution at `s_g`
(Scott's bandwidth). The arm-level crowding of `g` is the median of `d_g` over the folds in which
`g` was tested. Crowding is computed on the **unsmoothed** scores (it is a property of the prior,
not of the operator).

### 3.6 Controls

| control | what it is | what it rules out |
|---|---|---|
| **C0** | RWR baseline, λ = 0 | reference |
| **C1** | **random-neighbourhood smoothing**: `N(g)` = k uniformly drawn genes of the fold (same k, same λ, same weights, 30 replicates) | "any shrinkage of the score vector lowers the metric" |
| **C2** | **global-median shrinkage**: `m_g` = the median of *all* fold scores, same λ | "any shrinkage toward a typical value lowers the metric" |
| **C3** | **monotone flattening**: apply a monotone transform to the score vector (e.g. `tanh`) and recompute | "the metric moved because the scores were compressed, not because the ranking changed" |
| **C4** | **degree baseline left unsmoothed** (all arms) | keeps the ΔAUROC margin meaningful |

C1 is compared to the observed arm with a **robust Z** statistic based on the median absolute
deviation of the C1 null (the estimator used by Wright et al. 2025 for exactly this purpose).

### 3.7 Statistics

* Per-arm stability change: paired Wilcoxon signed-rank test over genes (baseline vs smoothed).
* Specificity vs C1: robust Z (MAD scale) of the observed change against the 30-replicate null.
* Multiplicity across the 24 arms: Benjamini–Hochberg within the exploratory family (the
  implementation in `scripts/run_stats_analysis.py`).
* Effect sizes are reported with the raw numbers; the BH q-values are secondary (this is an
  exploratory family).

### 3.8 The λ chosen by "our own stability metric as the training objective"

Because the brief asks for the stability metric to be the *training objective*, λ is selected by
**leave-one-repeat-out cross-fitting**:

1. For held-out repeat `r`: estimate the confidence weights from repeats `≠ r`.
2. On repeats `≠ r`, compute the arm-level stability objective `J(λ) = median_g SDrank_g(λ)` and
   the discrimination `D(λ) = mean over folds of ΔAUROC(λ)`.
3. Choose `λ_r = argmin_λ J(λ)` subject to `D(λ) ≥ D(0) − 0.010` (the P2 tolerance).
4. Report `λ_r`, and the out-of-sample stability and ΔAUROC of the selected scorer on repeat `r`.

This is the honest version of "train on stability": the reported stability gain of the *selected*
scorer is out-of-sample.

The full λ-grid curve is reported in addition (a curve is a description, not a selected model).

---

## 4. Outputs

| artefact | content |
|---|---|
| `d1/engine/stability.py` | the operators and metrics (crowding, confidence, CLMS, controls) |
| `scripts/run_stability_dump.py` | per-gene score dump + SHA-256 manifest + AUROC cross-check |
| `scripts/run_stability_regularized.py` | the λ × k grid with C1/C2/C3 over all 24 arms |
| `scripts/make_figure_stability_tradeoff.py` | the trade-off figure |
| `results/tables/stability_per_gene.tsv` | per-gene rank SD, score CV, crowding |
| `results/tables/stability_by_arm.tsv` | per-arm stability summary |
| `results/tables/stability_crowding_attribution.tsv` | crowding → instability attribution |
| `results/tables/stability_regularization_grid.tsv` | every (arm, α, variant, k, λ, control) cell |
| `results/tables/stability_tradeoff_curve.tsv` | the aggregated trade-off curve |
| `results/tables/stability_control_comparison.tsv` | C1/C2/C3 comparisons, robust Z, BH q |
| `results/tables/stability_selected_scorer.tsv` | the cross-fitted λ selection and its out-of-sample performance |
| `figures/figure7_stability_tradeoff.png` | the trade-off figure |
| `docs/stability_regularization.md` | method + results + limitations |
| `docs/manuscript_a_sections.md` | new exploratory section |
| `tests/test_stability.py` | unit tests (invariance, operator behaviour, no-label-leakage) |

## 5. Known limitations (stated now, not after the result)

1. **The operator is post-hoc.** It is not part of the pre-registered benchmark. It is reported as
   exploratory, in its own section, and cannot be used to support H1/H2/H3.
2. **SRS is a transductive ranker in its primary form.** It uses the score distribution of the
   evaluated fold. The inductive variant is reported alongside, but the transductive variant is
   the one that matches the brief's wording.
3. **AUROC is a rank metric** and therefore insensitive to monotone score transforms; this is
   exactly why P5 needs the additional dispersion and rank-change diagnostics.
4. **Crowding is measured in score space**, not in network space. A network-space neighbourhood
   variant is possible but is out of scope here.
5. **The stability metric is a proxy.** A gene's rank SD across seed-set resamples is not the same
   as its biological correctness; lowering it is a property of the scorer, not evidence of
   biology.
6. **The primary factorial is not re-decided.** Whatever SRS does to ΔAUROC, the confirmatory
   numbers remain the RWR numbers.
