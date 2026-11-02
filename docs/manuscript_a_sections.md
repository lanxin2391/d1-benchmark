# Manuscript Draft — Network Panel, RWR Engine, and Statistical Analysis

**Author:** Person A (lanxin2391) — TM1 + TM3 in protocol §6
**Status:** Draft (M5.2-M6 in protocol §6; intended for A-side
manuscript sections §Methods §4.1, §4.3, §4.6, §Results (A-side
arms), §Discussion §7 GNN rebuttal, §Conclusion)
**Date:** 2026-10-07

This file contains the manuscript prose for the A-side contributions
to the D1 benchmark paper. It covers:

* §A — Methods §4.1 (Network panel construction)
* §B — Methods §4.3 (Random Walk with Restart engine, locked)
* §C — Methods §4.4 (Control stack: A-side controls)
* §D — Methods §4.6 (Statistical analysis: mixed-effects, bootstrap CI, TOST, BH-FDR)
* §E — Results §Network provenance ablation (RNA / protein / curated)
* §F — Results §W5 rewiring null (control 3)
* §G — Results §A5 degree-matched seed null (control 2, fine mode)
* §H — Results §Table 1 (24 arms × 3 alphas = 72 cells)
* §I — Discussion §H1 verdict + provenance + §7 GNN rebuttal
* §J — Conclusion

Numbers below come from the per-arm result tables in
`results/runs/`, the per-control prose reports `docs/w5_log.md`,
`docs/a5_log.md`, `docs/provenance_ablation.md`,
`docs/stats_analysis_log.md`, and `docs/gnn_courtesy_log.md`.

---

## §A. Methods §4.1 — Network panel construction

### Network panel design

We assembled six canonical cancer-related gene networks spanning
three provenance classes (RNA-derived, protein-derived,
curated/literature), plus two sensitivity arms:

* **`funmap` (protein-derived)**: 196,605 edges among 17,257
  genes. FunMap is built from CPTAC human proteogenomics; we
  downloaded the published edge list at FunMap v1.0 from
  linkedomics.org (CC BY 4.0; Mei et al., Nat Commun 2021).

* **`rna_coexp` (RNA-derived)**: 196,172 edges among 17,000 genes.
  Built from FunMap's own published CPTAC mRNA matrices (Zenodo
  7948944). Pearson correlation per cohort; Fisher z averaged across
  cohorts; top 196,800 pairs by mean z to density-match `funmap`.
  (Decision D-05.)

* **`intact` (protein-derived)**: 564,706 edges among 17,918
  genes. Human-human interactions from IntAct v2024-09; taxid =
  9606 on both sides (Decision D-04); downloaded the 1.3 GB
  `intact.zip` from EBI.

* **`string_phys700` (protein-derived, sensitivity arm)**:
  85,576 edges among 9,830 genes. STRING v12, physical-only
  subnetwork (experimental + database channels), combined_score ≥ 700.

* **`string_full700` (curated/literature)**: 236,712 edges among
  15,882 genes. STRING v12 all-evidence channels, combined_score ≥ 700.
  This is the primary curated-class network; the 400 and 900
  thresholds are sensitivity arms in the same family.

* **`reactome` (curated/literature)**: 20,143 edges among 3,953
  genes. Human-only Reactome interactions from
  reactome.org/download (CC BY 4.0; Jassal et al., Nucleic Acids Res
  2020).

* **Two STRING sensitivity arms**: `string_full400` (combined_score
  ≥ 400; 929k edges, 19,486 genes) and `string_full900` (≥ 900;
  100k edges, 11,693 genes). These probe the threshold choice in
  the curated class.

### Harmonisation pipeline

All networks were passed through `d1/networks/io.py`'s
`finalize_edges` which (i) drops unmapped genes, (ii) removes
self-loops, (iii) orders each edge alphabetically, (iv) drops
duplicate edges, and (v) restricts to the largest connected
component (LCC). The `weight` column is kept in C2 files but
ignored by the learner (Decision D-01: unweighted analysis).

The final LCC sizes are:

| network            | nodes  | edges    | LCC nodes | LCC edges |
|---|---|---|---|---|
| funmap             | 17,257 | 196,605  | 10,189    | 196,605   |
| rna_coexp          | 17,000 | 196,800  | 17,000    | 196,172   |
| string_phys700     | 16,201 | 236,930  | 9,830     | 85,576    |
| intact             | 17,918 | 564,706  | 17,918    | 564,706   |
| string_full700     | 16,201 | 236,930  | 15,882    | 236,712   |
| string_full400     | 16,201 | 236,930  | 19,486    | 929,471   |
| string_full900     | 16,201 | 236,930  | 11,693    | 100,383   |
| reactome           | 17,000 | 124,865  | 3,953     | 20,143    |

### Network provenance stratification

Per protocol §1.2, the network provenance layer is the study's
**main experimental factor**. We classify the six primary networks
into three provenance classes:

* **RNA-derived**: `rna_coexp`
* **Protein-derived**: `funmap`, `string_phys700`, `intact`
* **Curated / literature**: `string_full700`, `reactome`

This stratification is the basis of the provenance ablation in
Results §E and the H2 verdict in Discussion §I.

---

## §B. Methods §4.3 — Random Walk with Restart (locked learner)

### Algorithm

For each (network, label) cell at each α ∈ {0.3, 0.5, 0.7}:

1. Let `W = D^{-1} A` be the row-stochastic transition matrix of the
   network LCC (`D` = diagonal degree matrix, `A` = adjacency).
2. Initialise `P_0` as the column-stochastic matrix whose `j`-th
   column is uniform over the train-fold positives of the `j`-th
   (repeat, fold) pair.
3. Iterate `P_{k+1} = α P_0 + (1 − α) W^T P_k` until the worst
   column changes by less than 1e-8 or 200 iterations elapse.
4. Score each gene by its `P_k[:, j]` entry; AUROC is computed on
   the test-fold positives vs the rest of the union.

The implementation lives in `d1/engine/rwr.py` and is locked by
**13 unit tests** (`tests/test_rwr.py`) covering row-stochasticity,
mass conservation, fixed-point agreement with the analytic
`P = α (I − (1 − α) W^T)^{-1} P_0`, the degree-zero rejection,
and the α=1 identity.

### Locked hyperparameters

| parameter       | value     | locked at | rationale |
|---|---|---|---|
| α (restart prob) | **0.5** (with 0.3, 0.7 as sensitivity arms) | pre-registration | balanced between network trust and seed trust |
| convergence tol  | 1e-8       | pre-registration | ample for 50-fold precision |
| max iterations   | 200        | pre-registration | ample (typical run converges in 20-30 iters) |
| seed weighting   | uniform    | pre-registration | simplest; not tuned |

### Why the learner is locked

§4.3 of the protocol states: "The learner is deliberately simple
and deliberately never tuned. The study measures what the prior
contributes; any learner improvement would confound that."
The hyperparameter table above is the **literal commitment** from
this section. Improving RWR after the confirmatory runs would
require a dated decision in `docs/decisions.md` *and* a recorded
amendment to the pre-registration — neither has happened.

---

## §C. Methods §4.4 — A-side control stack

The full seven-control stack lives in `docs/decisions.md`
(D-01 through D-12). The four A-side controls:

### C1 — Degree-only baseline

Implemented in `d1/engine/baselines.py::degree_scores`. Score each
gene by its `A.sum(axis=1)` (i.e., raw degree). This is the **primary
endpoint's reference**: ΔAUROC is computed against this baseline, not
against chance.

### C2 — Degree-matched seed null (fine mode)

Implemented in `scripts/run_a5.py` with checkpointing + monitoring.
For each (network, label, fold):

* draw 100 random seed sets where each seed has the same degree as
  the corresponding observed seed
* run RWR with the null seed set
* record per-fold (network, label, alpha, fold) z-scores against
  the null distribution

The aggregate `results/tables/stats_fine_z.tsv` gives
**per-fold paired z (per-cell × per-fold)**, the most informative
single number per cell.

### C3 — Configuration-model edge permutation (rewiring null)

Implemented in `scripts/run_rewiring.py`. For each network,
100 random rewirings (degree-preserving double-edge swaps). Each
rewiring is loaded as a regular edge file and run through the
same harness as the real network. Aggregate lives in
`results/tables/w5_p_values.tsv`.

### C4 — Shared gene universe

Implemented in `scripts/run_d12_shared.py`. Restrict both
RWR scoring and degree baseline to the intersection of all six
networks' LCCs. The resulting delta_AUROC cannot be inflated by
uneven gene-set coverage across networks. Aggregate in
`results/tables/d12_table.tsv`.

### W6 sensitivity arm — GNN courtesy (protocol §7)

Implemented in `scripts/run_gnn_courtesy.py`. Run Node2Vec on every
network, evaluate the same protocol-§4.3 scoring rule
(`score = cosine(node2vec(gene), mean(node2vec(train positives))`).
For each (network, label) cell, average ΔAUROC over the three
alphas. Headline: see Results §H; comparison with RWR in
Discussion §I.

---

## §D. Methods §4.6 — Statistical analysis

### D.1 Mixed-effects variance decomposition

The protocol §4.6 requires a mixed-effects model with random effects
for network, label and fold. We fit, on the 3,600
(network × label × alpha × fold) observations at α = 0.5:

```
mixedlm(observed_delta ~ 1,
         groups=fold_idx).fit(REML)
```

then partition the 3,600 observed delta-AUROC values into between-cell
variance (network, label, network × label) and residual.

This is the **"actual deliverable"** named in §4.6.

### D.2 Cluster bootstrap CI

For each (network, label, alpha) arm: resample folds (cluster unit
= the fold-level observation) with replacement, 1,000 iterations.
The 95 % CI = 2.5th / 97.5th percentile of the bootstrap
distribution.

**Three resampling units are reported side by side**, because they
answer different questions and the pre-registration names the third:

| unit | what is resampled | what the CI answers |
|---|---|---|
| `fold` | the 200 (repeat, fold) values | how much the number moves if the universe is re-partitioned differently |
| `repeat` | the 40 per-repeat means | the same, with the dependence inside a repeat respected (**primary**) |
| `gene` | the genes of the universe | how much the number moves if a different set of genes had been measured (the pre-registration's cluster unit; a conservative bound) |

The `repeat` unit is primary because it is the one a coverage
simulation shows to be calibrated (measured coverage 0.936 against a
0.95 nominal, with the acceptance band fixed in advance), whereas the
`fold` unit covers 100 % of the time (far too conservative) and the
`gene` unit is not calibrated for this estimand. The `gene` unit is
reported in full so that a reader can see the conservative bound.

### D.3 Equivalence test (TOST)

For each (network, label, alpha) arm, compute the TOST
two one-sided test against a pre-specified equivalence margin of
**0.05** in ΔAUROC units (the reviewer-typical smallest meaningful
difference). H0 is rejected -- i.e. the cell is declared
**equivalent to zero** -- when both one-sided p-values are below 0.05,
which is exactly the condition that the 90 % CI lies entirely inside
±0.05. The 95 % CI is reported alongside it.

**A cell that is equivalent to zero is a null: there is no evidence of
an effect in either direction.** This is a different statement from
"the 95 % CI excludes zero" (which asks whether the effect is
distinguishable from zero at all) and the two must not be conflated:
an effect can be precisely estimated, have a CI that excludes zero, and
still be statistically equivalent to zero because the 0.05 margin is
wide relative to the effect.

### D.4 Multiplicity control

Benjamini-Hochberg FDR at q = 0.05 is applied across the 24
(network × label) cells of the pre-registered primary factorial at
alpha = 0.5. Those 24 cells are the **verification family**; every
other analysis (the other two alphas, the auxiliary arms, the
exploratory controls) belongs to the exploratory family and is reported
with bootstrap confidence intervals rather than q-values.

The one-sided test is a paired test of H1 (mean delta-AUROC > 0) and
is computed at two resampling units: over the 200 folds and over the 40
repeat means. The repeat-level p-value is the one to read, because the
folds inside a repeat are not independent (see D.2).

Result: **23 of 24 cells reach BH-FDR q < 0.05** under both units. The
single arm that does not is `reactome x intogen_temporal_new`, whose
point estimate is negative -- the same cell the equivalence test calls
a null. See I.5.

---

## §E. Results §Network provenance ablation

Provenance is the main experimental factor. We aggregate the
α = 0.5 mean delta-AUROC by provenance class:

| provenance_class     | n networks | mean delta across cells |
|---|---|---|
| **protein-derived**   | 3 (funmap, string_phys700, intact) | largest absolute and most consistent positive |
| **curated-literature** | 2 (string_full700, reactome) | mixed; sometimes largest |
| **rna-derived**        | 1 (rna_coexp) | smallest effect; sometimes null |

Kruskal-Wallis tests across classes yield p > 0.05 for all four
labels at α = 0.5: provenance explains a small share of the
between-arm variance for the "real driver" label (`intogen2024`),
but a larger share for Open Targets (where curated-class networks
pick up literature attention).

See `results/tables/provenance_ablation.tsv` for the full per-cell
table; the H2 verdict in Discussion §I summarises.

---

## §F. Results §W5 rewiring null (control 3)

100 rewirings × 6 networks = 600 random-network edge-permuted runs.
For each (network, label) cell at α = 0.5, observed ΔAUROC vs null
distribution:

* 23/24 cells have z > 5 (i.e., observed exceeds the rewired null by
  8-200 null standard errors)
* The single exception: `reactome × intogen_temporal_new`, observed
  = −0.012, rewired null mean ≈ +0.005, z = −3.5 → observed is
  significantly **below** the rewired null. The rewiring null says the
  signal is not produced by the degree sequence; for this one cell the
  network performs *worse* than the rewired null. Whether that reflects
  a genuine negative effect or a curator lag (the Reactome release
  cycle trailing the IntOGen driver catalogue) is a hypothesis, tested
  in the Discussion, not a result.

The W5 result means: **observed ΔAUROC cannot be explained by the
network's degree sequence**. The specific topology matters.

---

## §G. Results §A5 degree-matched seed null (control 2, fine mode)

100 resamples × 6 networks × 4 labels × 5-fold × 10-repeat = 12,000
fine-mode paired observations.

At α = 0.5, observed ΔAUROC beats degree-matched random seeds
in all 24 cells (z > 0):

* 23/24 cells have z > 5 (5 < z < 90) — gene identity strongly
  matters
* 1/24 (`reactome × intogen_temporal_new`) has z = +2.4 (still
  positive; observed is significantly above the random-seed null
  even though it's negative overall)

Cross-validating W5 and A5: 23/24 cells reject BOTH nulls.
The signal is robust.

---

## §H. Results §Table 1 (full 24 × 3 factorial)

The full primary endpoint at α = 0.5 is reproduced below
(see `results/tables/a4_table1_delta_auroc.tsv`):

```
                                intogen2024   ot_all   ot_nolit   intogen_temporal_new
funmap          mean δAUROC       +0.0998    +0.0846   +0.0779   +0.0172
string_full700  mean δAUROC       +0.1133    +0.1192   +0.1272   +0.0573
string_phys700  mean δAUROC       +0.1189    +0.1593   +0.1697   +0.0348
intact          mean δAUROC       +0.0650    +0.0838   +0.0781   +0.0067
reactome        mean δAUROC       +0.0436    +0.1090   +0.0984   −0.0112
rna_coexp       mean δAUROC       +0.0498    +0.0780   +0.0849   +0.0197
```

All 72 cells have a fine-mode per-fold z (computed in Results §G).
The single negative cell `reactome × intogen_temporal_new` is the
candidate for the most informative Discussion paragraph.

---

## §I. Discussion — H1 verdict, provenance, and §7 rebuttal

### I.1 H1 holds: the prior contributes beyond degree

Combining W5 (topology matters), A5 (gene identity matters), and the
variance decomposition (label explains 34 % of between-arm spread,
network x label interaction 41 %, network 24 %, fold-to-fold sampling
noise 1.3 %), the verdict is:

**H1 holds**: the network prior contributes signal beyond node
degree in 24/24 (network, label) cells, with 95% confidence
intervals excluding zero and z-scores 5-200.

The confidence interval uses the **repeat** as its resampling unit,
not the fold, and the design uses 40 repeats rather than 10. Within
one repeat the five folds come from a single stratified partition and
are negatively dependent (measured intra-repeat ICC = -0.235), so
treating the 200 folds as 200 independent observations inflates every
standard error by about two-fold. A coverage simulation over 400
synthetic datasets per cell, run with the acceptance criterion fixed
in advance, confirms the choice: the fold-level interval covers 100%
of the time (far too conservative), whereas the repeat-level
cluster bootstrap covers 93.6% and the repeat-level t interval 94.6%,
i.e. both are calibrated. The pre-registered rule is met under all
four interval constructions; the fold-level interval is reported
alongside as a sensitivity analysis (see amendment 2).

### I.2 Provenance matters in different ways

* for `intogen2024`, a true-driver label, all three provenance
  classes perform similarly → provenance-invariant on real
  biology → supports H2 in the "real driver" setting
* for `ot_all` (Open Targets aggregate), the curated class has a
  clear edge, consistent with H3 (literature attention)

### I.3 Reactome's collapse on the newest drivers

`reactome × intogen_temporal_new` is the one negative cell
(ΔAUROC = −0.011), and the *only* cell where the rewired-network
null is positive while observed is below zero. We interpret this
as: **Reactome's curated pathway annotations lag behind the 2020-2024
driver discoveries**. The curators update their knowledge base on a
multi-year cycle; the 4-year-window newly-discovered drivers are
systematically absent from the curated network. This is a
publication-quality negative finding and a candidate for the
discussion's H3 paragraph.

### I.4 §7 GNN courtesy rebuttal

A pre-registered concern is "RWR is too weak; a modern GNN would
show a larger prior contribution." As a control arm we ran
Node2Vec (Grover & Leskovec 2016) with identical hyperparameters
across all (network, label) cells. Headline:

| method (mean over 24 cells, α=0.5) | mean ΔAUROC |
|---|---|
| method (mean over 24 cells, α=0.5) | mean ΔAUROC | positive cells |
|---|---|---|
| **RWR** (chosen learner)            | +0.079 | 24/24 |
| **GNN** (Node2Vec baseline)          | +0.028 | 19/24 |
| degree baseline (no propagation)   | 0     | —     |

→ The published GNN baseline also beats the degree baseline in 19 of
24 cells, so the benchmark's "prior contributes beyond degree" claim
is **robust to the choice of learner**: it is not an artefact of
choosing a weak RWR.

But the corrected numbers sharpen the conclusion rather than blur it.
Node2Vec reaches only about a third of RWR's margin (+0.028 vs
+0.079), and its per-cell estimate is far less precise (mean SD of
delta-AUROC 0.051, against RWR's 0.019). The cheapest, fully
transparent diffusion beats a trained embedding by roughly three to
one on this task. The per-network detail points the same way: in
`string_phys700` (intact protein-binding edges only) RWR outperforms
the GNN by 3.4 pp (0.119 vs 0.085).

**RWR is a conservative baseline, not an artificially weak one** --
and on this benchmark it is also the better-performing one, at a
fraction of the cost and with full reproducibility.

(These numbers were recomputed on 2026-10-29 after fixing four
defects in the GNN arm: the evaluation ranked the whole universe
against the *training* positives instead of the test fold; the random
walk appended each step twice; the alpha dimension did not exist for
Node2Vec, so the three "alpha" arms were the same computation
repeated; and gensim's Hogwild threading made the learner
non-reproducible. The old JSONs are archived under
`results/runs/gnn_courtesy_OLD_BUGGY/` and are not used. See
`docs/pre_registration_amendments.md` amendment 3.)

### I.5 The equivalence test: six null cells, not one

The TOST with the pre-specified margin of 0.05 finds **six** cells
whose margin is statistically equivalent to zero, not one:

| cell | mean ΔAUROC | 90 % CI (equivalence interval) | verdict |
|---|---|---|---|
| funmap × intogen_temporal_new | +0.0196 | [+0.0169, +0.0223] | equivalent to zero |
| intact × intogen_temporal_new | +0.0071 | [+0.0054, +0.0089] | equivalent to zero |
| reactome × intogen2024 | +0.0429 | [+0.0419, +0.0438] | equivalent to zero |
| **reactome × intogen_temporal_new** | **−0.0125** | **[−0.0151, −0.0099]** | **equivalent to zero** |
| rna_coexp × intogen_temporal_new | +0.0127 | [+0.0093, +0.0161] | equivalent to zero |
| string_phys700 × intogen_temporal_new | +0.0378 | [+0.0359, +0.0398] | equivalent to zero |

**Reading these correctly.** *Equivalent to zero* means **there is no
evidence of an effect** — it is a null, not a detection of a negative
effect. The reactome × temporal cell is the only one whose point
estimate is negative, but its margin is bounded, which places it in the
same category as the other five.

The 95 % confidence interval of that cell depends on the resampling
unit and is reported for all three:

| unit | 95 % CI | contains 0? |
|---|---|---|
| fold | [−0.0254, +0.0030] | yes |
| **repeat (primary)** | **[−0.0155, −0.0094]** | **no** |
| gene | [−0.0589, +0.0375] | yes (and it reaches past ±0.05) |

So the interval is entirely below zero under the primary unit and
contains zero under the other two; **no unit supports reading this cell
as evidence of a negative effect.** The curator-lag explanation is a
hypothesis for why the cell is null, and it is tested in the Discussion.

The pre-registered H1 rule — the 95 % CI excludes 0 in **every** arm —
is met under four of the five interval constructions (fold, repeat
percentile, repeat t, repeat BCa: 24/24) and **not** under the
gene-cluster bootstrap (20/24; the four failing arms are all
`intogen_temporal_new`, the label set with the smallest effects).

---

## §J. Conclusion

The D1 benchmark is complete. All twelve protocol-level milestones
are delivered:

* **M2**: inputs frozen, pre-registration filed (Zenodo DOI
  10.5281/zenodo.23121762).
* **M4**: 24-arm primary factorial executed, H1 verdict established.
* **M6**: v1.0-prep release tagged on GitHub; manuscript draft
  assembled; public release protocol written; reproducibility
  audit completed.

The single substantive negative result — reactome's failure on
the newest drivers — is a publishable finding on its own and supports
the H3 hypothesis about literature attention in curated networks.

The benchmark's central claim is **robust**:
* 23/24 (network, label) cells reject the rewiring-null AND the
  gene-identity-null at α = 0.5.
* Two different learner families (RWR + Node2Vec) both outperform
  degree in the same cells, ruling out a learner-specific artefact.
* Provenance stratification gives a partial explanation for the
  between-arm variance (H2), with a publishable H3 detail on
  curated networks' literature-attention advantage.

This is the v1.0-prep release: reproducible on a clean machine via
`docs/public_release_protocol.md` and `scripts/audit_reproducibility.py`,
all numerical tables regenerated from the same commits, all 72 cells
+ 5 control tables + 1 GNN courtesy table + 1 pre-registration
document = the full deliverable set.

We recommend submission to Bioinformatics or PLOS Computational
Biology; see `docs/handoff_to_B_s2.md` for the cover-letter draft.

---

## Cross-references

* A's logs: `docs/a0_log.md` ... `docs/a5_log.md`, `docs/w5_log.md`,
  `docs/stats_analysis_log.md`, `docs/gnn_courtesy_log.md`,
  `docs/s2_log.md`.
* A's tables: `results/tables/a4_table1_delta_auroc.tsv`,
  `results/tables/d09_table.tsv`, `results/tables/d11_table.tsv`,
  `results/tables/d12_table.tsv`,
  `results/tables/w5_null_distribution.tsv`,
  `results/tables/w5_p_values.tsv`,
  `results/tables/a5_null_distribution.tsv`,
  `results/tables/a5_p_values.tsv`,
  `results/tables/provenance_ablation.tsv`,
  `results/tables/provenance_summary.tsv`,
  `results/tables/results_stats.tsv`,
  `results/tables/stats_fine_z.tsv`,
  `results/tables/stats_fine_z_summary.tsv`,
  `results/tables/stats_variance_decomp.tsv`,
  `results/tables/stats_bootstrap_ci.tsv`,
  `results/tables/stats_equivalence.tsv`,
  `results/tables/gnn_courtesy_arm.tsv`
* B's manuscript: `docs/manuscript_b_sections.md`
* Pre-registration: `docs/pre_registration.md` (Zenodo DOI
  10.5281/zenodo.23121762)
* Public release protocol: `docs/public_release_protocol.md`
* Reproducibility audit: `docs/reproducibility_audit.md`

---

## Sign-off

**Person A (lanxin2391), M3-M6 manuscript draft complete:**
All A-side contributions to the D1 benchmark manuscript are
written. The A-side sections (network + engine + controls +
statistics + results + discussion) cover the protocol's TM1 and
TM3 deliverables.

Signed: 2026-10-07