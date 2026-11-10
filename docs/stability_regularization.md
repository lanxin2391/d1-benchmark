# Stability of a network prior: what it costs, what buys it, and what cannot be smoothed away

**Status:** complete (2026-11-05). **Family:** EXPLORATORY — it does not enter the H1 decision.
**Pre-registration:** `docs/stability_regularization_plan.md` (committed before any result existed;
Amendment A1 at the end of that file records the one design change that was made after the first
result was seen, and why).
**Artefacts:** `results/tables/stability_*.tsv`, `figures/figure7_stability_tradeoff.png`,
`scripts/run_stability_*.py`, `d1/engine/stability.py`.

---

## 0. One-paragraph summary

We pre-registered a **stability-regularised scorer** (SRS): shrink each gene's prior score toward
the confidence-weighted median of its *score neighbourhood*, and use our own stability metric
(per-gene rank SD across cross-validation repeats) as the training objective. The predictions were
that the per-gene rank SD would fall, that ΔAUROC would not, and that this would trace out a
stability-versus-discrimination trade-off curve. **The first prediction failed, and the reason is
structural rather than incidental.** The operator is *order-preserving*: the k nearest genes in
score space have, by construction, nearly identical scores, so their weighted median is the gene's
own score, the shrinkage compresses the scores without touching the ranking, and the rank SD
cannot move (max change 5.3 × 10⁻⁵ across 24 arms and eight shrinkage strengths). Both matched
controls confirm the direction: every operator that *does* move the ranking makes the prior
**less** stable (random-neighbourhood smoothing raises the median rank SD from 0.028 to 0.109 at
λ = 0.4 and collapses ΔAUROC to −0.11 at λ = 1). What we did find is a quantitative law for the
instability itself: the per-gene rank SD of a network prior falls as **n_seeds^−0.49** (24/24 arms
monotone), i.e. it is the sampling noise of a mean over seeds, not a property of a crowded score
profile; and because RWR is *linear* in the seed weights, bagging over seed resamples is provably a
no-op for it. The trade-off curve therefore exists, but it does not run along the axis the brief
proposed: it runs along the **number of seeds**, and on that axis stability and discrimination
improve *together*, so the unsmoothed prior is the Pareto point.

---

## 1. Design

| element | value |
|---|---|
| arms | 6 networks × 4 label sets = **24** (native universe, uncapped label positives, as in the released primary panel) |
| networks | `funmap`, `rna_coexp`, `string_phys700`, `intact`, `string_full700`, `reactome` |
| labels | `intogen2024`, `ot_all`, `ot_nolit`, `intogen_temporal_new` |
| learner | RWR (locked, protocol §4.3), α = 0.5 primary; α ∈ {0.3, 0.7} as sensitivity |
| CV | **40 repeats × 5 folds**, new panel `data/processed/splits_v40/` |
| primary endpoint | ΔAUROC against the degree baseline (unchanged) |
| primary stability endpoint | **median per-gene rank SD** (percentile rank inside the fold, SD across the 40 repeats) |

Two mandatory verifications were run before any result was interpreted:

1. **Split identity.** Repeats 0–9 of the new 40-repeat panel are byte-identical to the released
   10-repeat panel for all 24 arms (`scripts/build_splits_v40.py` fails loudly otherwise). The
   released panel is never modified.
2. **Score-dump identity.** The AUROC recomputed from the dumped per-gene score matrices equals the
   AUROC returned by the harness `d1.engine.evaluate.run_arm` with **max |Δ| = 0.0** for all
   24 arms × 3 α (`scripts/run_stability_dump.py`; manifest with SHA-256 of all 72 dumps in
   `results/tables/stability_per_gene_scores_manifest.tsv`).

## 2. The operator

For one evaluated fold with test-gene scores `s_g`:

1. **Score-space neighbourhood.** `N(g)` = the `k` genes with the smallest `|s_j − s_g|` (self
   excluded), ties broken by ascending gene symbol. `k ∈ {5, 20, 50}`; primary `k = 20`.
2. **Confidence-weighted local median.** `m_g` = the minimiser of `Σ_{j ∈ N(g) ∪ {g}} c_j|s_j − v|`,
   where `c_j = 1/(1 + CV_j²)` and `CV_j` is gene `j`'s score coefficient of variation across the
   repeats **other than the repeat being evaluated** (leave-one-repeat-out, so no evaluated repeat
   contributes to its own weights).
3. **Shrinkage.** `s'_g = s_g + λ_g (m_g − s_g)`. Variant A: `λ_g = λ` for all genes (primary).
   Variant B (crowding-adaptive): `λ_g = λ · crowd_g / max_g crowd_g`.
4. **Evaluation.** `s'` is ranked with the D-07 tie-break and scored with the unchanged metrics
   (AUROC, AUPRC, pAUC@FPR ≤ 0.02, MedRR, p@k). The degree baseline is **never** smoothed. The
   operator never reads a label vector.

Controls: **C1** random neighbourhood of the same size (20 replicates, matched strength),
**C2** shrinkage toward the fold's global median, **C3** a monotone flattening of the score vector.

## 3. Results

### 3.1 The stability atlas (the reproducible part)

Per-gene statistics for all 24 arms are in `results/tables/stability_per_gene.tsv`
(rank SD, score CV, score mean/SD, crowding) and `results/tables/stability_by_arm.tsv`.

| arm family | median per-gene rank SD (percentile units) |
|---|---|
| all 24 arms | 0.0187 (intact × temporal) – 0.0437 (reactome × temporal) |
| median over arms | **0.0277** |

A gene's rank moves, across seed-set resamples, by about **2.8 % of the list** (median arm).

### 3.2 The pre-registered operator does not change the ranking

`results/tables/stability_prediction_scorecard.tsv` evaluates every pre-registered prediction
against the rule fixed in advance.

| prediction | rule | observed | verdict |
|---|---|---|---|
| **P1** rank SD decreases | monotone decrease in ≥ 24/24 arms | max \|Δ rank SD\| = **5.3 × 10⁻⁵** over all arms and λ; no arm decreases | **NOT SUPPORTED** |
| **P2** ΔAUROC does not fall | median \|Δ(ΔAUROC)\| ≤ 0.010 | max \|Δ\| = **4.1 × 10⁻⁵** | supported, vacuously |
| **P3** a trade-off curve | report the price of stability | the SRS curve is a **single point**; both controls move to *worse* stability *and* *worse* discrimination | **no trade-off in the smoothing direction** |
| **P4** specificity vs C1 | SRS below the 95th percentile of the random null in ≥ 24/24 arms | **72/72** (arm, λ) cells; the random control raises the median rank SD from 0.0277 to **0.163** at λ = 0.4 | supported, vacuously |
| **P5** not a flattening artefact | monotone transforms leave rank statistics unchanged | a monotone power transform changes the AUROC by ≤ **3.0 × 10⁻⁵** (float ties only) and the rank of no gene in the checked fold; the SRS operator compresses the score SD by up to 40 % while leaving every rank metric unchanged | supported |
| **P6** crowding-specific | crowded genes gain more | score-space crowding: ρ = **−0.10** (negative in 21/24 arms); rank-space crowding: ρ = **+0.23** (positive in 24/24); the gene's own score CV: ρ = **+0.82** (positive in 24/24) | **NOT SUPPORTED** for score-space crowding, partial for rank-space crowding |
| **P7** survives leaving the fold | inductive variant keeps ≥ half the effect | the inductive variant gives the *identical* rank SD (0.02771) and ΔAUROC | not applicable — there is no effect |

### 3.3 Why the operator cannot work (the mechanism)

Three facts, all measured:

1. **The neighbourhood median is the gene's own score.** The `k` nearest genes *in score space*
   have, by construction, almost the same score, so the confidence-weighted local median `m_g ≈ s_g`
   and `s'_g ≈ s_g`. At λ = 1 the score dispersion falls to 49–63 % of the original, yet
   0.03–0.5 % of the ranks change and the rank SD moves by ≤ 5.3 × 10⁻⁵.
2. **Shrinkage cannot separate signal from noise in a gap.** For two genes `g, g'` in the same
   crowded block, `m_g ≈ m_g'`, so the gap becomes `(1 − λ)(s_g − s_g') + λ(m_g − m_g') ≈ (1 − λ)(s_g − s_g')`.
   The *signal* and the *noise* in the gap are scaled by the same factor, so the probability that
   the pair swaps order is unchanged. Shrinking toward a neighbourhood therefore cannot reduce
   rank volatility: it can only leave it unchanged (λ below the flip threshold) or make it worse.
3. **Every ranking-changing operator makes it worse.** C1 (random neighbourhood, same strength)
   moves the median rank SD from 0.0277 to 0.109 (λ = 0.4), 0.168 (λ = 0.7) and 0.287 (λ = 1), and
   drives ΔAUROC from +0.081 to −0.108. C2 (global median) is flat until λ = 1, where every score
   collapses onto one value, the D-07 symbol tie-break takes over, and the rank SD jumps to 0.111
   while ΔAUROC falls to −0.147.

**Conclusion.** The per-gene rank SD of a network prior is already at its minimum for the score
vector it is computed from: the ranking the prior produces is the most informed ranking available,
and any departure from it injects independent noise. Post-hoc smoothing of the score vector
therefore cannot buy stability.

### 3.4 What does predict the instability — and what buys it

`results/tables/stability_crowding_attribution.tsv` (per-arm Spearman ρ against the per-gene rank
SD, 24 arms):

| predictor | median ρ | sign consistency |
|---|---|---|
| the gene's own score CV | **+0.824** (0.742 – 0.861) | positive in 24/24 |
| rank-space crowding (±1 % of ranks) | +0.229 (0.054 – 0.288) | positive in 24/24 |
| score-space crowding (± Scott bandwidth) | **−0.101** (−0.331 – +0.093) | negative in 21/24 |
| node degree | +0.023 | none |

So the upstream premise "crowding causes rank instability" is **not** supported in its score-space
form (the sign is if anything reversed), and only weakly supported in rank space; what predicts
instability is the gene's **own score noise**.

That points at the seed set, and the seed-set experiment
(`results/tables/stability_seed_size_scaling.tsv`) closes the argument. Truncating the training
seed set to a fraction `f` of the available positives and re-running the whole evaluation:

| f | median seeds | median per-gene rank SD | median ΔAUROC |
|---|---|---|---|
| 0.05 | 15 | 0.1339 | 0.0388 |
| 0.10 | 30 | 0.1206 | 0.0442 |
| 0.25 | 76 | 0.0917 | 0.0541 |
| 0.50 | 152 | 0.0632 | 0.0674 |
| 1.00 | 304 | **0.0277** | **0.0809** |

* the rank SD decreases monotonically in **24/24 arms**;
* the ΔAUROC increases monotonically in **20/24 arms**;
* the fitted exponent of the rank SD against the seed count is **−0.488** (median; IQR −0.51 to
  −0.46; `results/tables/stability_seed_size_exponents.tsv`) — i.e. the textbook 1/√n law for the
  sampling error of a mean.

**Therefore the instability of a network prior is the sampling noise of a mean over seeds.** Two
consequences follow, and one of them is a theorem:

* *Consequence 1.* RWR is linear in the seed weights, `P = α(I − (1−α)Wᵀ)⁻¹P₀`. Bagging over
  bootstrap resamples of the seed set therefore returns exactly `P₀` — **no amount of seed
  resampling can stabilise an RWR prior by averaging**. The only stabiliser is *more distinct
  seeds*, i.e. a larger training set.
* *Consequence 2.* The stability/discrimination trade-off the brief asked for does not exist on
  this axis either: more seeds improve **both**. What the seed axis does buy is a lever on the
  benchmark's own design — the training fraction of the CV (5 folds ⇒ 80 % of the positives are
  available as seeds). Halving it doubles the rank SD (0.028 → 0.063) and costs only ≈ 10 % of
  the ΔAUROC (0.100 → 0.090 for `funmap × intogen2024`).

### 3.5 The circularity guard, and what it caught

The brief asked for the stability metric to be the *training objective*, so λ was selected by
leave-one-repeat-out cross-fitting over four groups of repeats: the objective (median per-gene rank
SD) is minimised on three groups subject to `ΔAUROC ≥ ΔAUROC(0) − 0.010`, and the selected scorer
is then scored on the held-out group.

The procedure selected **λ = 1 in 44/96 and λ = 0 in 49/96** held-out groups. Since the objective
is flat in λ to within 10⁻⁵, the selection is decided by noise — and the held-out groups confirm
it: the out-of-sample stability and discrimination of the selected scorer are identical to the
unsmoothed prior. This is the circularity guard doing its job: *optimising the stability metric on
itself produces a model that optimises nothing.*

## 4. What this changes for the benchmark

1. **Do not smooth the prior.** There is no stability to buy and discrimination to pay; the
   unsmoothed prior is the Pareto point. Any "stability-regularised" prior built by post-processing
   the score vector is either a no-op or a degradation.
2. **Report the instability instead of hiding it.** The per-gene rank SD (0.019 – 0.044 percentile
   units across the 24 arms) is a cheap, reproducible property of a prior and belongs in the
   benchmark's output: a third party who wants a stable prior can read it off.
3. **The seed count is a design parameter.** The 1/√n law means the benchmark's headline ΔAUROC is
   partly a function of the training-set size, and the ranking stability is *mostly* a function of
   it. Both should be reported next to the margin.
4. **The "crowding causes instability" mechanism should not be repeated** in its score-space form.
   The measured driver is the gene's own score noise; crowding contributes weakly and only in rank
   space.

## 5. Limitations

1. **Exploratory.** This analysis changes the learner's output and therefore never enters H1/H2/H3.
   The confirmatory factorial is untouched; the RWR numbers are the reference.
2. **The stability metric is a proxy.** A gene's rank SD across seed-set resamples is not the same
   as its biological correctness. Lowering it would be a property of the scorer, not evidence of
   biology — which is exactly why the two controls were mandatory, and why they are decisive here.
3. **The operator was specified in score space.** A network-space neighbourhood (the gene's
   1-hop neighbours) is a different operator and was only probed, not pre-registered; it behaves
   like the random control (it changes the ranking and raises the rank SD).
4. **One label panel, one universe.** The arms use the native universe and the uncapped label
   positives, exactly as the released primary panel; the shared-universe (D-12) and capped (D-11)
   arms were not re-run.
5. **α sensitivity is at the headline setting only.** The α = 0.3 and 0.7 runs reproduce the
   qualitative result (order-preserving operator, flat curve) but do not carry the full control
   stack.
6. **The seed-count law is measured, not modelled.** The fitted exponent (−0.49) matches the 1/√n
   prediction; a formal test of the linearity of RWR in the seed weights was not run (it is a
   one-line algebraic check on the closed form).
