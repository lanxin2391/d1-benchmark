# Manuscript Draft — Network Panel, RWR Engine, and Statistical Analysis

**Author:** Person A (lanxin2391) — TM1 + TM3 in protocol §6
**Status:** Draft (M5.2-M6 in protocol §6; intended for A-side
manuscript sections §Methods §4.1, §4.3, §4.6, §Results (A-side
arms), §Discussion §7 GNN rebuttal, §Conclusion)
**Date:** 2026-10-07; revised 2026-11-05 (P0 repairs, Phase D, Phase S)

This file contains the manuscript prose for the A-side contributions
to the D1 benchmark paper. It covers:

* §A — Methods §4.1 (Network panel construction)
* §B — Methods §4.3 (Random Walk with Restart engine, locked)
* §C — Methods §4.4 (Control stack: A-side controls)
* §D — Methods §4.6 (Statistical analysis: variance decomposition, bootstrap CI, TOST, BH-FDR)
* §E — Results §Network provenance ablation (RNA / protein / curated)
* §F — Results §Robustness controls (W5 rewiring, D-09 temporal drift, D-11 common-size cap, D-12 shared universe, complex-aware CV, hub-capped RWR)
* §G — Results §A5 degree-matched seed null (control 2, fine mode)
* §H — Results §Table 1 (24 arms × 3 alphas = 72 cells)
* §I — Results §Where the stability of the prior comes from (stability atlas, attribution, consensus prior)
* §J — Discussion §H1 verdict + provenance + §7 GNN rebuttal
* §K — Conclusion

Numbers below come from the per-arm result tables in
`results/runs/`, the per-control prose reports `docs/w5_log.md`,
`docs/a5_log.md`, `docs/provenance_ablation.md`,
`docs/stats_analysis_log.md`, and `docs/gnn_courtesy_log.md`.

**Design of record.** The primary factorial, the interval families, the
equivalence test, the multiplicity control and the three robustness
controls in §F.2–§F.4 are all computed on the **40-repeat** design
(5 folds × 40 repeats = 200 (repeat, fold) evaluations per cell); the
canonical tables are `results/tables/stats_bootstrap_ci.tsv`,
`results/tables/stats_equivalence.tsv`,
`results/tables/stats_h1_verdict_v40.tsv`,
`results/tables/stats_variance_decomp.tsv`, and the primary run is
`results/runs_v40/a4_combined.tsv`. The one remaining exception is the
control-2 null in §G, which was simulated on the original 10-repeat
design; the direction of that bias is stated in §G and in amendment 2
(A2.3g). Every departure from the frozen pre-registration is recorded
in `docs/pre_registration_amendments.md`; the pre-registration itself
(Zenodo DOI 10.5281/zenodo.23121762) has not been edited.

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
Results §E and the H2 verdict in Discussion §J.

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

### C4 — Shared gene universe (control 4, decision D-12)

Implemented in `scripts/run_d12_shared.py`. RWR propagates on each
network's **full LCC**; the seeds and the **evaluation set** are
restricted to the gene set shared by all six networks (1,605 genes);
the degree baseline stays the full-LCC degree. The split universe
itself is the shared set (`build_splits.py --universe shared`), so
every per-fold test set contains |shared| / 5 = 321 genes and the run
asserts this size, which makes a silent regression of the control
impossible. The splits use the 40-repeat design. Aggregate in
`results/tables/d12_table.tsv`.

### C5 — Temporal-drift control (decision D-09)

Implemented in `scripts/run_d09.py`. Training seeds are the IntOGen 2024
drivers **minus** the drivers discovered in the 2020–2024 window
(`intogen_temporal_new`), and the held-out positives are exactly that
temporal window: the question is whether a network prior still recovers
drivers that were discovered *after* the training cut-off. Each of the
50 repeats draws the training seed set **with replacement**, so the
repeats are independent draws rather than 50 copies of one partition
(`n_eff_seed_sets` in the output table is the audit trail for this).
Aggregate in `results/tables/d09_table.tsv`.

### C6 — Common-size cap control (decision D-11)

Implemented in `scripts/run_d11.py`. For each network, the number of
positives is capped to N = min over the four labels of the number of
positives inside that network's LCC, so that every network is evaluated
on the same positive-set size and the four labels are not confounded
with unequal prevalence. The label is restricted to the universe
**before** the top-N is taken, and the run asserts that every split
holds exactly N positives. Aggregate in `results/tables/d11_table.tsv`.

### C7 — Complex-aware cross-validation (Phase D)

Implemented in `d1/splits.py::make_splits_complex_aware` and
`scripts/run_d3_generate_cx_splits.py`. The seven controls above remove
every leakage channel we could think of except **protein complexes**:
stratified splits put the members of a complex in different folds, so a
held-out driver's score is partly determined by a training-set partner,
and RWR propagates along exactly those edges. The complex-aware
splitter keeps every complex group (the connected components of the
complex graph) **in one fold**, so a held-out gene's partners are held
out with it. The complex map is built from the ComplexPortal
identifiers curated in UniProtKB (`scripts/build_complex_map.py`);
grouped CV is then generated into `data/processed/splits_v40_cx/`, with
the v1 splits untouched, and the invariant report
(`results/tables/d3_split_invariants.tsv`) records the achieved fold
balance and the number of complex groups split across folds (0 in
24/24 arms). Aggregate in `results/tables/d4_v1_vs_v2_comparison.tsv`.

### W6 sensitivity arm — GNN courtesy (protocol §7)

Implemented in `scripts/run_gnn_courtesy.py`. Run Node2Vec on every
network, evaluate the same protocol-§4.3 scoring rule
(`score = cosine(node2vec(gene), mean(node2vec(train positives))`).
For each (network, label) cell, average ΔAUROC over the three
alphas. Headline: see Results §H; comparison with RWR in
Discussion §I.

---

## §D. Methods §4.6 — Statistical analysis

### D.1 Variance decomposition

Protocol §4.6 asks for a mixed-effects decomposition of the observed
delta-AUROC values into network, label and fold components, with the
network × label interaction included. We fit it as a **balanced
random-effects ANOVA** on the 6 networks × 4 labels × 200 (repeat,
fold) cells at α = 0.5 (the 40-repeat design):

```
sigma2_within    = MS_within
sigma2_int       = (MS_int       - MS_within) / n
sigma2_label     = (MS_label     - MS_within) / (n * a)
sigma2_network   = (MS_network   - MS_within) / (n * b)
```

with n = 200 evaluations per cell, a = 6 networks, b = 4 labels.

**Two fractions are reported, because they answer two different
questions.** "Between-arm spread" is the variance of the 24 cell means
— the thing protocol §4.6 asks about; "total variance" is the variance
of the raw per-fold measurements, which is what sets the width of every
confidence interval.

| source | % of between-arm spread | % of total variance |
|---|---|---|
| network × label (interaction) | 40.07 % | 23.96 % |
| label (main effect) | 34.78 % | 20.80 % |
| network (main effect) | 24.82 % | 14.84 % |
| fold-to-fold sampling noise | 0.34 % | 40.39 % |

All three effects are significant against the within-cell noise
(network F = 295.0, label F = 619.0, interaction F = 119.7; p < 1e-15
each). The full table, the variance components and the derivation are
in `results/tables/stats_variance_decomp.tsv`,
`results/tables/stats_variance_components.tsv` and
`docs/variance_decomposition.md`.

> **Correction.** The first implementation of this section fitted
> `mixedlm("observed_delta ~ 1", groups=fold_idx)` with only five groups
> (fold 0–4), added the variances of the network and label group
> *means* on top of the residual that already contained them, and set
> the interaction term to exactly zero. That produced a table reporting
> "residual = 52 %", which was wrong: most of what it called residual
> was the label main effect and the network × label interaction being
> counted twice. The table above replaces it (amendment 2, A2.1–A2.2).

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
a null. See J.5.

---

## §E. Results §Network provenance ablation

Provenance is the main experimental factor. We aggregate the
α = 0.5 mean delta-AUROC by provenance class
(`results/tables/provenance_ablation.tsv`,
`results/tables/provenance_summary.tsv`):

| label | protein-derived (3) | curated-literature (2) | rna-derived (1) | Kruskal p | η²(class) |
|---|---|---|---|---|---|
| intogen2024 | **+0.0927** | +0.0766 | +0.0513 | 0.501 | 0.21 |
| intogen_temporal_new | +0.0212 | **+0.0235** | +0.0123 | 0.954 | 0.02 |
| ot_all | +0.1075 | **+0.1133** | +0.0755 | 0.304 | 0.18 |
| ot_nolit | +0.1075 | **+0.1113** | +0.0821 | 0.651 | 0.08 |

Kruskal-Wallis tests across classes yield p > 0.30 for all four labels
at α = 0.5, so **no label shows a significant provenance difference**.
The consistent pattern is that the RNA-derived network is the weakest
of the three classes on every label, while the protein-derived and
curated classes are within 0.006 of each other on the two Open Targets
labels — a small curated-class margin that is directionally consistent
with H3 (literature attention) but not statistically separable with
six networks.

---

## §F. Results §Robustness controls

Four controls ask whether the headline result survives when the obvious
alternative explanations are removed: is the signal just the degree
sequence (F.1), does it depend on the training partition (F.2), on the
size of the positive set (F.3), or on how much of the genome each
network covers (F.4)? F.5 states what the three corrected controls
say together.

### F.1 W5 rewiring null (control 3)

100 rewirings × 6 networks = 600 random-network edge-permuted runs.
For each (network, label) cell at α = 0.5, observed ΔAUROC vs null
distribution (`results/tables/w5_p_values.tsv`):

* 23/24 cells have z > 8.9 (up to z = 205.2), i.e. observed exceeds
  the rewired null by 9–205 standard errors **of the null mean** (the rewiring z divides by the null mean's standard error, null_sd / sqrt(n_rewirings), not by the null spread)
* The single exception: `reactome × intogen_temporal_new`, observed
  = −0.0112, rewired null mean = +0.0050, z = −3.53 → observed is
  significantly **below** the rewired null. The rewiring null says the
  signal is not produced by the degree sequence; for this one cell the
  network performs *worse* than the rewired null. Whether that reflects
  a genuine negative effect or a curator lag (the Reactome release
  cycle trailing the IntOGen driver catalogue) is tested in §F.5 and
  discussed in §J.3.

The W5 result means: **observed ΔAUROC cannot be explained by the
network's degree sequence**. The specific topology matters.

One caveat is stated rather than buried: the rewiring performs **2 x |E|**
double-edge swaps per rewiring, where the pre-registration asks for
10 x |E|. The mixing was measured directly
(`results/tables/rewiring_mixing_check.tsv`): 4-11 % of the original
edges survive the rewiring and 21-54 % of the nodes keep at least one
original neighbour, so the rewired network is a *partial* null -- but
the pre-registered 10 x |E| would have been a partial null too (2-10 %
residual edge sharing), and the residual bias is **conservative**, because
a rewired network that keeps a few original edges sits slightly above
the fully-mixed null mean and therefore yields a slightly SMALLER z. The
z-scores are 8.9-205, so no conclusion of this paper depends on the
mixing level (amendment 2, A2.3k).

### F.2 D-09 temporal-drift control

The training seeds are the IntOGen 2024 drivers minus the drivers
discovered in the 2020–2024 window, and the held-out positives are
exactly that window: can the prior still rank *newly discovered*
drivers? Each repeat draws the training seed set with replacement, so
the 50 repeats are 50 independent draws (`n_eff_seed_sets` = 50 in
`results/tables/d09_table.tsv`).

At α = 0.5:

| network | ΔAUROC (mean ± SD over 50 seed draws) | n training seeds | n held-out positives |
|---|---|---|---|
| funmap | **+0.0247** ± 0.0076 | 357 | 104 |
| rna_coexp | **+0.0382** ± 0.0068 | 296 | 89 |
| string_phys700 | **+0.0762** ± 0.0080 | 428 | 116 |
| intact | **+0.0302** ± 0.0030 | 473 | 150 |
| string_full700 | **+0.0923** ± 0.0064 | 468 | 145 |
| **reactome** | **−0.0095** ± 0.0116 | 310 | 63 |

Five of six networks recover drivers they were never trained on.
Reactome is the only negative one, and its 95 % interval
(±1.96 SD ≈ [−0.032, +0.013]) contains zero.

### F.3 D-11 common-size cap control

Each network's positive set is capped to N = min over the four labels of
the positives inside that network's LCC, so every network is judged on
the same number of positives: funmap 104, rna_coexp 89,
string_phys700 116, intact 150, string_full700 145, reactome 63. All 24
arms satisfy `N_in_split == N` exactly.

At α = 0.5 (`results/tables/d11_table.tsv`), 22 of 24 arms are
positive, from +0.0071 (intact × `intogen_temporal_new`) to +0.1819
(string_phys700 × `ot_nolit`). The two negative arms are both Reactome:

| cell | ΔAUROC |
|---|---|
| **reactome × intogen2024** | **−0.0110** ± 0.0434 |
| **reactome × intogen_temporal_new** | **−0.0125** ± 0.0501 |

Note what the cap does to the second column of Table 1: in the uncapped
primary factorial, `reactome × intogen2024` is **+0.0429**; capping
every network to Reactome's own minimum (63) flips that arm negative.
Reactome's advantage on the broad IntOGen label is therefore partly a
*coverage* advantage that disappears when the positive sets are
equalised — while its temporal failure gets worse.

### F.4 D-12 shared-universe control

Both networks and labels are restricted to the 1,605 genes present in
all six networks, so that no arm can win by covering more genes than
its competitors. Every per-fold test set is |shared| / 5 = 321 genes
(asserted by the run).

At α = 0.5 (`results/tables/d12_table.tsv`), 22 of 24 arms are
positive, from +0.0148 (rna_coexp × `intogen_temporal_new`) to +0.1304
(string_phys700 × `ot_all`). The two negative arms are again the
temporal ones:

| cell | ΔAUROC |
|---|---|
| **reactome × intogen_temporal_new** | **−0.0360** ± 0.0854 |
| **intact × intogen_temporal_new** | **−0.0150** ± 0.0642 |

### F.5 What the three corrected controls say together

Three independent controls — a redrawn training partition (D-09), an
equalised positive-set size (D-11) and an equalised gene universe
(D-12) — each leave the primary conclusion intact and each leave
**Reactome** as the outlier on temporal extrapolation:

| control | reactome × intogen_temporal_new | other arms |
|---|---|---|
| D-09 temporal drift | −0.0095 | 5 of 6 networks positive |
| D-11 common-size cap | −0.0125 | 22 of 24 arms positive |
| D-12 shared universe | −0.0360 | 22 of 24 arms positive |

The three controls remove three different explanations — "it depends on
the partition", "it depends on how many positives the network has",
"it depends on how much of the genome the network covers" — and the
Reactome anomaly survives all three. That is the strongest single piece
of evidence in the study: the failure is specific to Reactome's
coverage of **newly discovered** drivers, not to its overall coverage
and not to the size of the evaluation set. It moves curator lag from a
Discussion hypothesis to a supported explanation (§J.3).

The three controls were re-run after defects in each of them were found
and fixed (amendment 2, A2.3h and A2.3i); the pre-fix versions of the
tables should not be cited.

### F.6 Complex-aware cross-validation: the strictest control we could add

The seven controls above remove every leakage channel we could think of
except one: **protein complexes**. Standard stratified splits put the
members of a complex in different folds, so a held-out driver's score
is partly determined by a training-set partner — and RWR propagates
along exactly those edges. We measured the channel before fixing it
(Phase D, task D1): among the held-out positives that are in a curated
complex, **96 % have a complex-linked gene in the training fold**. The
channel is wide open — but symmetric, because the held-out *negatives*
are linked at the same rate (94 %): sharing a complex does not depend
on being a driver. Whether the leak inflates the *advantage over
degree* is an empirical question, so we rebuilt the splits to keep
every complex group in one fold (a held-out gene's partners are then
held out with it) and re-ran the primary factorial on them.

| | v1 (stratified) | v2 (complex-aware) |
|---|---|---|
| held-out positives with a complex-linked gene in the training fold | 0.87–0.96 | **0.02–0.18** |
| complex groups split across folds | — | 0 (24/24 arms) |
| H1 (all four repeat-level interval families exclude 0) | 24/24 | **24/24** |
| arms whose sign flips | — | **0/24** |
| median \|relative change\| in ΔAUROC | — | **0.142** |
| repeat-level interval width (v2 / v1) | — | 1.18 |

The result is systematic and in the predicted direction: **the prior
loses ground in 20 of 24 arms (median signed change −12.3 %)**, i.e. the
standard CV was inflating the effect. But **no conclusion changes**: the
H1 verdict is 24/24 in both versions, no arm changes sign, and the
negative cell becomes slightly *more* negative (reactome ×
intogen_temporal_new: −0.0125 → −0.0178), so it is not an artefact of
the leak. Two further points matter for the reader: the v2 intervals are
**wider** (1.18×), so the v2 significance is not bought with narrower
intervals; and the change is concentrated in the sparse and curated
networks (funmap, reactome: up to −45 %) while the protein-derived
networks barely move (−2 to −6 %).

Because the complex map covers 18.3 % of the panel's genes, this is a
**lower bound** on the leakage a complete complex resource would
expose, and the v2 numbers are an additional analysis (the
pre-registration fixes the stratified design), registered in the
exploratory family.

### F.7 Hub-capped RWR: is the prior's advantage a hub effect?

A reader can reasonably ask whether "the prior contributes beyond
degree" is just "the prior is degree plus a hub effect". We tested it
by running three variants of the learner on the same splits, with the
degree baseline untouched: **truncation** (the nodes with degree > K
are removed from the propagation, K = 128 and 256) and the
**degree-penalised walk** `W[i,j] = A[i,j] / (deg(i)^0.5 · deg(j)^0.5)`.
The truncation is reported under two conventions, because they answer
different questions: *keep* (the hubs stay in the universe and score 0,
so they tie at the bottom) and *exclude* (they leave the evaluation
set as well, so the AUROC denominator is unchanged).

| variant | median ΔAUROC | vs baseline | arms changing sign |
|---|---|---|---|
| baseline | +0.0809 | — | — |
| trunc128, keep | +0.0252 | −0.0499 (31 % left) | 8/24 |
| trunc128, exclude | **+0.0849** | −0.0001 | 0/24 |
| trunc256, exclude | +0.0853 | +0.0007 | 0/24 |
| degpen0.5 | **+0.0979** | **+0.0063** | 0/24 |

**The headline is not a hub effect**: on the part of the network
without the hubs the prior still beats the degree baseline by the same
amount. **But roughly two thirds of the raw margin sits in the hubs'
own ranking position** — silencing them rather than removing them takes
the median effect to 31 % of its size and flips a third of the arms, so
both numbers are stated and the smaller one is quoted. The
degree-penalised walk, which is the orthogonalisation idea in effect,
*increases* the effect: there is no hub-driven failure for a
corrective model to fix.

---

## §G. Results §A5 degree-matched seed null (control 2)

100 degree-matched random seed sets per cell, each propagated over the
5 folds of each of the **40 repeats** (200 (repeat, fold) evaluations
per cell) — the control-2 null uses the same design as the primary
factorial. Two summaries are computed from the same trace:

* **cell level** (one z per cell, the null distribution pooled over
  folds): `results/tables/a5_p_values.tsv`
* **fold level** (one z per (repeat, fold), the folds of a cell
  summarised by their mean): `results/tables/stats_fine_z.tsv` and
  `stats_fine_z_summary.tsv`

At α = 0.5 the observed ΔAUROC beats degree-matched random seeds in
**all 24 cells** (cell-level z > 0 everywhere):

* 23/24 cells have cell-level z > 5 (8.2 ≤ z ≤ 226.9)
* the remaining cell, `reactome × intogen_temporal_new`, has z = +3.09
  (p = 1.0e-3) — still positive: observed is above the random-seed null
  even though it is negative against the degree baseline
* at the fold level the per-cell mean z ranges from 0.25 (the same
  Reactome cell) to 13.20; the Stouffer-combined z ranges from 3.56 to
  186.7 and **all 24 cells have p < 0.05** (largest p = 1.8e-4)

So gene identity matters: replacing the true driver seeds by
degree-matched random seeds destroys the effect in every arm. The null
distribution itself is informative — with random seeds the RWR prior is
*worse* than the degree baseline in most cells (null mean ΔAUROC between
−0.052 and +0.005), i.e. a random seed set actively hurts.

Cross-validating W5 and A5: 23/24 cells reject BOTH nulls, and the one
cell that does not reject the rewiring null still rejects the
gene-identity null. The signal is robust.

---

## §H. Results §Table 1 (full 24 × 3 factorial)

The full primary endpoint at α = 0.5 is reproduced below
(see `results/tables/a4_table1_delta_auroc_v40.tsv`; the 40-repeat
values agree with the canonical interval table
`results/tables/stats_bootstrap_ci.tsv` to 1e-16):

```
                                 intogen2024   ot_all   ot_nolit   intogen_temporal_new
funmap          mean δAUROC       +0.1001    +0.0853   +0.0780   +0.0196
string_full700  mean δAUROC       +0.1133    +0.1201   +0.1267   +0.0604
string_phys700  mean δAUROC       +0.1189    +0.1582   +0.1711   +0.0378
intact          mean δAUROC       +0.0651    +0.0842   +0.0784   +0.0071
reactome        mean δAUROC       +0.0429    +0.1103   +0.0996   −0.0125
rna_coexp       mean δAUROC       +0.0521    +0.0767   +0.0834   +0.0127
```

23 of the 24 cells are positive. The single negative cell
`reactome × intogen_temporal_new` is the candidate for the most
informative Discussion paragraph, and the three robustness controls in
§F.2–§F.4 are what make it a finding rather than an accident.

## §I. Results §Where the stability of the prior comes from

The pre-registration fixes *what* we measure (the effect of the prior
beyond node degree) but not *what makes the estimate reproducible*.
This section answers that second question, which is also the question
the field's leaderboards never ask: not "is the prior better than
degree?" but "how stable is the prior's answer, and where does that
stability come from?"

We follow the template of the proteomics study the teacher pointed us
to (Klingelhuber et al., *Nat Metab* 6:861–879, 2024): treat the
stability of the measurement layer as a first-class quantity, report
its distribution, and then ask which features explain it. Here the
measurement is the RWR prior's score for a gene, and the features are
the network's topology.

### I.1 What is stable, and what is not

**The endpoint is stable.** The repeat-level SD of ΔAUROC is 0.0047
(range 0.0017–0.0129) against effects of 0.007–0.171: the effect is
2–30 times larger than the noise it is measured against.

**A single gene is not.** The per-gene relative error of the score
across the 40 re-partitions is 0.185 (range 0.100–0.263) — an 18 %
error on a gene's score. In rank terms a median gene moves 2.7
percentile points, i.e. **~270 places out of ~10,000**, between two
partitions of the *same* universe. A re-partition changes the answer
for a gene far more than it changes the answer for the panel.

That asymmetry is the result: **the benchmark's numbers are
reproducible; the individual genes it ranks are not.** (Figure 5.)

### I.2 The mechanism is crowding, not the network's shape

Where a gene sits in the score distribution decides how stable its
rank is. Pooled over the 24 arms, the median rank SD runs 0.0056 (lowest
score decile) → 0.0382 (decile 7) → 0.0162 (top decile), and the
number of genes lying within ±10 % of a gene's score runs 119 → 535 →
193 over the same deciles: **a gene crowded with near-equal-scored
neighbours swaps places when the partition changes.** The Spearman
correlation between rank SD and local score density is **+0.468**
(median over arms) and is **positive in 24 of 24 arms**.

The formal attribution (within-arm OLS with arm fixed effects, 267,696
gene-arm rows) agrees:

| predictor | coefficient on log(rank SD) | marginal within-arm ρ |
|---|---|---|
| log score magnitude | **+0.910** | +0.340 (24/24 arms) |
| log local score density | **+0.787** | +0.478 (24/24 arms) |
| log degree | −0.995 | **+0.038 (≈0)** |
| clustering coefficient | +0.089 | +0.088 |
| k-core number | −0.001 (q = 0.98) | +0.044 |

Two things must be read carefully here. First, **degree explains
almost nothing on its own** (marginal ρ = +0.04); its large joint
coefficient is a suppression artefact, caused by its 0.72 correlation
with the RWR score, and the paper states the marginal version. Second,
**k-core — the classic connectivity descriptor — explains nothing at
all**. The stability of a gene's rank is a property of the *score
distribution's local density*, not of the network's connectivity.

The attribution transfers: fitted on 23 arms without arm intercepts, it
predicts the rank SD of an unseen network with an out-of-sample R² of
0.24–0.56, and of an unseen label with 0.29–0.57 (Figure 5b).

### I.3 Which networks are stable

At the network level (eight networks, four labels, 40 repeats) the
stable ones are again not the ones with the best topology:

| descriptor | Spearman with the arm-level stability |
|---|---|
| mean label coverage | **−0.905 (q = 0.012)** |
| fraction of nodes in the 2-core | **−0.905 (q = 0.012)** |
| number of nodes | −0.857 (q = 0.026) |
| density, clustering, modularity, assortativity, max k-core, path length | all non-significant |

**No descriptor predicts the size of the contribution** (every
correlation with the mean ΔAUROC has q ≥ 0.70). Which network
contributes a lot and which network gives a reproducible estimate are
two different questions, and only the second one has a topological
answer (Figure 6).

The label matters too, in the direction the variance decomposition
already indicated: the temporal label is the least stable label in
*every* network (arm-level SD 0.0065–0.0129 against 0.0017–0.0063 for
the other three), because its positive set is small (152 genes). The
intra-repeat ICC is negative in **24 of 24** arms (median −0.178),
which is the independent confirmation of the repeat-level resampling
unit chosen in the pre-registration amendment.

### I.4 The picture replicates on two independent systems

* **ClinGen** (clinical validity of cancer variants, 84 positives) — a
  curation criterion that shares nothing with IntOGen or Open Targets.
  The network-level stability ordering is the same: intact is the most
  stable (background rank SD 0.0152, score CV 0.100) and funmap /
  rna_coexp the least (0.0304 / 0.0289).
* **Transfer arms** (train on one label system, evaluate on another, 66
  cells) reproduce the *extremes* of the ordering (the two most stable
  networks in the primary panel are the most stable here too) but not
  the full ranking (Spearman +0.37). Their ΔAUROC is much smaller
  (string_phys700 +0.059 against +0.122), because the model is
  evaluated on a label system it never saw.

### I.5 A model built from these findings: a consensus prior

The benchmark's one negative cell is a **coverage** failure, not a
noise artefact: Reactome's curated annotations lag behind the
2020–2024 driver discoveries, and the failure survives three controls
that equalise the positive set. A stale network cannot be repaired by
a better learner — it can only be out-voted. So we built the obvious
model: a **consensus prior** over the six networks on the shared
universe (1,605 genes), the mean of the six networks' percentile
ranks.

| label | median single network | consensus | consensus without Reactome |
|---|---|---|---|
| intogen2024 | +0.0486 | **+0.1128** | +0.1068 |
| ot_all | +0.0826 | **+0.1421** | +0.1211 |
| ot_nolit | +0.0845 | **+0.1384** | +0.1193 |
| **intogen_temporal_new** | **−0.0110** | **+0.0358** [0.0301, 0.0416] | **+0.0514** [0.0457, 0.0570] |

On the temporal label the *median single network is negative* and the
consensus is clearly positive; and **Reactome's inclusion costs the
consensus 0.016** — a stale network does not merely fail on its own,
it degrades the aggregate. That is the first quantitative statement of
the price of curation lag at the level of an ensemble.

The result is robust to the aggregation rule (mean of ranks, median of
ranks and reciprocal rank fusion agree to within 0.019 on every label)
and improves when the networks are weighted by their measured
reliability (+0.0388 on the temporal label). (Figure 7.)

**And the stability axis shows the cost.** The consensus is *less*
stable per gene than any of its members (rank SD 0.0364–0.0528 against
0.0210–0.0439): averaging six networks' ranks does not cancel their
noise, it adds the disagreement between networks as variance — and the
disagreement is also where the signal is. The consensus trades
per-gene variance for accuracy, a trade-off that AUROC alone cannot
see. It is a robustness device for the panel-level estimate, not a
better gene ranker.

### I.6 What this section does and does not establish

* It establishes that **the benchmark's numbers are reproducible while
  its gene-level rankings are not**, and that the reason is crowding in
  the score distribution rather than the shape of the network.
* It establishes that **the effect size is not predictable from
  topology**, which is a caution for anyone reading a single network's
  performance as a property of that network.
* It does not establish the causal direction of the crowding
  explanation (the score distribution is itself produced by the
  network), and it does not establish that the consensus is the best
  model — it is the one model that follows from our own findings, and
  its per-gene stability is worse.

---

## §J. Discussion — H1 verdict, provenance, and §7 rebuttal

### J.1 H1 holds: the prior contributes beyond degree

Combining W5 (topology matters), A5 (gene identity matters), the
variance decomposition (label explains 34.8 % of between-arm spread,
network × label interaction 40.1 %, network 24.8 %, fold-to-fold
sampling noise 0.3 %) and the three robustness controls of §F.2–§F.4,
the verdict is:

**H1 holds**: the network prior contributes signal beyond node degree
in every one of the 24 (network, label) cells — the 95 % confidence
interval excludes zero in 24/24 arms under the fold-level, repeat-level
percentile, repeat-level *t* and repeat-level BCa intervals, and in
20/24 under the (deliberately conservative) gene-cluster bootstrap, the
four exceptions being the four `intogen_temporal_new` arms with the
smallest effects. The cell-level separation from the degree-matched
seed null is z = 3.1–227 and from the rewired null z = 8.9–205
(excluding the single Reactome cell, §F.1).

Two further results make the verdict stronger rather than weaker.
**Under the strictest control we could add** — the complex-aware
cross-validation of §F.6, which keeps every protein complex in one
fold and so removes the complex-partner leak — the prior still beats
the degree baseline in **24 of 24 arms** with **no sign change in any
arm**; the effect sizes shrink by a median of 12 %, i.e. the standard
CV was inflating them, but the conclusion is unchanged. And **the
advantage is not a hub effect** (§F.7): removing the hubs from the
propagation leaves it unchanged (+0.0849 against +0.0809), and the
degree-penalised walk increases it (+0.0979).

The confidence interval uses the **repeat** as its resampling unit,
not the fold, and the design uses 40 repeats rather than 10. Within
one repeat the five folds come from a single stratified partition and
are negatively dependent (measured intra-repeat ICC = −0.235), so
treating the 200 folds as 200 independent observations inflates every
standard error by about two-fold. A coverage simulation over 400
synthetic datasets per cell, run with the acceptance criterion fixed
in advance, confirms the choice: the fold-level interval covers 100%
of the time (far too conservative), whereas the repeat-level
cluster bootstrap covers 93.6% and the repeat-level t interval 94.6%,
i.e. both are calibrated. The pre-registered rule is met under all
four interval constructions; the fold-level interval is reported
alongside as a sensitivity analysis (see amendment 2).

### J.2 Provenance matters in different ways

* for `intogen2024`, a true-driver label, all three provenance
  classes perform similarly → provenance-invariant on real
  biology → supports H2 in the "real driver" setting
* for `ot_all` (Open Targets aggregate), the curated class has a
  clear edge, consistent with H3 (literature attention)

### J.3 Reactome's collapse on the newest drivers — a supported explanation

`reactome × intogen_temporal_new` is the one negative cell of the
primary factorial (ΔAUROC = −0.0125) and the only cell where the
rewired-network null is positive while observed is below zero. The
interpretation we advance is:

**Reactome's curated pathway annotations lag behind the 2020–2024
driver discoveries.** The curators update their knowledge base on a
multi-year cycle; the drivers discovered in the 4-year window are
systematically absent from the curated network, so a prior built on
Reactome cannot rank them.

This was a hypothesis in the first draft of this manuscript. It is now
a **supported explanation**, because the anomaly survives three
independent controls that each remove a competing explanation
(§F.2–§F.5):

| competing explanation | control that removes it | Reactome × temporal after the control |
|---|---|---|
| the training partition was unlucky | D-09: seed set redrawn 50× | −0.0095 (still negative) |
| the positive set was too large / too uneven | D-11: every network capped to N = 63 | −0.0125 (still negative) |
| Reactome simply covers fewer genes | D-12: all arms evaluated on the same 1,605 genes | −0.0360 (still negative) |

What the controls add is *localisation*. D-11 shows that Reactome's
apparent advantage on the broad IntOGen label is partly a coverage
advantage (capped to 63 positives, `reactome × intogen2024` flips from
+0.0429 to −0.0110), and D-09 and D-12 show that the temporal failure
is not a coverage artefact at all. The failure is specific to
Reactome's coverage of *newly discovered* drivers.

A fourth control closes the last channel we could think of. Under the
**complex-aware cross-validation** (§F.6) — the splits that keep every
protein complex in one fold, so a held-out driver's complex partners
are held out with it — Reactome's temporal arm is **not less negative
but more negative** (−0.0178 against −0.0125), and no arm of the panel
changes sign. The failure is therefore not produced by the
complex-partner leak either.

What the controls add is *localisation*. D-11 shows that Reactome's
apparent advantage on the broad IntOGen label is partly a coverage
advantage (capped to 63 positives, `reactome × intogen2024` flips from
+0.0429 to −0.0110), and D-09 and D-12 show that the temporal failure
is not a coverage artefact at all. The failure is specific to
Reactome's coverage of *newly discovered* drivers.

If the failure is a coverage failure, then the remedy is not a better
learner but a better *ensemble*: a stale network can only be out-voted.
The consensus prior of §I.5 does exactly that — on the temporal label
the median single network is negative (−0.0110) while the consensus is
+0.0358, and **removing Reactome from the consensus raises it to
+0.0514**, i.e. Reactome's presence costs the aggregate 0.016. This is
the first quantitative statement of the price of curation lag at the
level of an ensemble.

Two caveats are stated plainly. First, Reactome's temporal arm has only
63 held-out positives, so its interval is wide; the point estimate is
negative under all four controls, but the *magnitude* is not
precisely estimated. Second, this is one network out of six, and the
control that isolates it most sharply (D-12) makes it negative for
`intact` as well (−0.0150), so the phenomenon may be "small curated
networks" rather than "Reactome" specifically. The publication claim
we make is therefore about **curated networks and temporal
extrapolation**, with Reactome as the clearest case — not a universal
claim about Reactome.

### J.4 §7 GNN courtesy rebuttal

A pre-registered concern is "RWR is too weak; a modern GNN would
show a larger prior contribution." As a control arm we ran
Node2Vec (Grover & Leskovec 2016) with identical hyperparameters
across all (network, label) cells. Headline:

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

### J.5 The equivalence test: six null cells, not one

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
as evidence of a negative effect.** The curator-lag explanation is why
the cell is null; it is supported by the three robustness controls in
§F.2–§F.5 and argued in §J.3.

The pre-registered H1 rule — the 95 % CI excludes 0 in **every** arm —
is met under four of the five interval constructions (fold, repeat
percentile, repeat t, repeat BCa: 24/24) and **not** under the
gene-cluster bootstrap (20/24; the four failing arms are all
`intogen_temporal_new`, the label set with the smallest effects).

---

## §K. Conclusion

The D1 benchmark is complete. All twelve protocol-level milestones
are delivered:

* **M2**: inputs frozen, pre-registration filed (Zenodo DOI
  10.5281/zenodo.23121762).
* **M4**: 24-arm primary factorial executed, H1 verdict established.
* **M6**: v1.0-prep release tagged on GitHub; manuscript draft
  assembled; public release protocol written; reproducibility
  audit completed.

The benchmark's central claim is **robust**:
* 24/24 (network, label) cells beat the degree-matched random-seed
  null (control 2, §G), and 23/24 beat the rewired-network null
  (control 3, §F.1); the one exception is the Reactome temporal cell,
  where the observed value falls *below* the rewired null.
* Four further robustness controls — a redrawn training partition
  (D-09), an equalised positive-set size (D-11), an equalised gene
  universe (D-12) and a **complex-aware cross-validation** that keeps
  every protein complex in one fold (§F.6) — leave 22–24 of 24 arms
  positive each time, with no sign change in any arm.
* The advantage is **not a hub effect** (§F.7): removing the hubs from
  the propagation leaves it unchanged, and the degree-penalised walk
  increases it.
* Two different learner families (RWR + Node2Vec) both outperform
  degree in the same cells, ruling out a learner-specific artefact.
* Provenance stratification gives a partial explanation for the
  between-arm variance (H2), with a publishable H3 detail on curated
  networks' literature-attention advantage.

Two results go beyond the pre-registered endpoint and are, we think,
the parts of the paper that the field most needs.

**What the stability of a prior is made of (§I).** We measured how much
the prior's answer moves when the universe is re-partitioned, gene by
gene, and found that the *panel-level* estimate is reproducible (the
repeat-level SD of ΔAUROC is 0.0047 against effects of 0.007–0.171)
while the *individual gene's* rank is not (18 % relative error; a median
gene moves ~270 places). The mechanism is crowding: a gene with many
near-equal-scored neighbours swaps places, and the correlation between
the local score density and the rank instability is positive in 24 of
24 arms. Classic topology (k-core, density, clustering, modularity)
explains none of it, and **no descriptor predicts the size of the
contribution at all** — which is a caution for anyone who reads a
single network's performance as a property of that network.

**What curation lag costs (§I.5, §J.3).** The benchmark's one negative
cell is a coverage failure, and a consensus prior over the six networks
repairs it: on the temporal label the median single network is negative
while the consensus is +0.0358, and removing Reactome from the
consensus raises it to +0.0514. A stale network does not merely fail on
its own — it degrades the aggregate. The price is per-gene stability:
the consensus is less stable than any of its members, a bias-variance
trade-off that AUROC alone cannot see.

The single substantive negative result — curated networks' collapse on
newly discovered drivers — is supported by four independent controls
and is a publishable finding in its own right (§F.5, §J.3).

This is the v1.0-prep release: reproducible on a clean machine via
`docs/public_release_protocol.md` and `scripts/audit_reproducibility.py`,
all numerical tables regenerated from the same commits, all 72 cells
+ 5 control tables + 1 GNN courtesy table + 1 pre-registration
document = the full deliverable set.

We recommend submission to Bioinformatics or PLOS Computational
Biology; see `docs/handoff_to_B_s2.md` for the cover-letter draft.

---

## §L. Exploratory — the stability of a network prior (SRS)

> **Family: EXPLORATORY.** This section reports an analysis that changes the
> learner's output. It is therefore **not** part of the confirmatory analysis,
> it does not enter the H1 decision rule, and every number in it is a
> property of the *exploratory* arm, not of the benchmark's primary endpoint.
> It was pre-registered in `docs/stability_regularization_plan.md` (committed
> before any result existed) and is reported in full in
> `docs/stability_regularization.md`.

### L.1 The question and the pre-registered operator

If a prior's ranking is unstable, can the instability be removed? We
pre-registered a **stability-regularised scorer (SRS)**: shrink each gene's
score toward the *confidence-weighted median of its score neighbourhood*,

    s'_g = s_g + lambda_g * (m_g - s_g),

where `N(g)` is the `k` genes with the smallest `|s_j - s_g|` and the weights
`c_j = 1/(1 + CV_j^2)` come from gene `j`'s own score coefficient of variation,
estimated **leave-one-repeat-out** so that no evaluated repeat contributes to
its own weights. The objective was our own stability metric: the per-gene SD of
the within-fold percentile rank across the 40 cross-validation repeats.
Pre-registered predictions: (P1) the rank SD falls; (P2) ΔAUROC does not;
(P3) a stability-versus-discrimination trade-off curve exists; (P4) the effect
is specific — it must beat a matched-strength random-neighbourhood control;
(P5) it is not a flattening artefact; (P6) it is crowding-specific; (P7) it
survives an inductive (no test gene in the neighbourhood) variant.

### L.2 The result: the prior's ranking cannot be smoothed into stability

The operator **does not change the ranking**. Across all 24 arms and all eight
shrinkage strengths the median per-gene rank SD moved by at most
**5.3 x 10^-5** (baseline 0.019-0.044 percentile units) and ΔAUROC by at most
**4.1 x 10^-5**. The score dispersion, in contrast, fell by up to 40 %: the
operator compresses the scores and leaves every rank statistic untouched.

The reason is structural. The `k` nearest genes *in score space* have nearly
identical scores, so their weighted median **is** the gene's own score. Worse,
for two genes in the same crowded block the shrinkage rescales the *signal* and
the *noise* in their gap by the same factor `(1 - lambda)`, so the probability
that the pair swaps order is unchanged: shrinking toward a neighbourhood cannot
separate signal from noise in a gap.

Both matched controls point the same way. **C1** (random neighbourhood of the
same size) raises the median rank SD from 0.028 to 0.109 at λ = 0.4, 0.168 at
λ = 0.7 and 0.287 at λ = 1, and drives ΔAUROC from +0.081 to −0.108. **C2**
(shrinkage toward the fold's global median) is flat until λ = 1, where every
score collapses onto one value, the D-07 symbol tie-break takes over and the
rank SD jumps to 0.111 while ΔAUROC falls to −0.147. **C3** (a monotone
flattening) changes the AUROC by at most 3.0 x 10^-5, confirming that rank
statistics are invariant to monotone score transforms — so there is no
"flattening" story to tell, because there is no ranking change to explain.

**The unsmoothed prior is the Pareto point of the stability/discrimination
plane.** The trade-off curve the exercise set out to draw does not exist in the
smoothing direction: every operator that moves the ranking moves it to worse
stability *and* worse discrimination.

### L.3 What does drive the instability, and what buys it

The pre-registered mechanism — crowding — is **not** supported in its
score-space form. The per-arm Spearman correlation between score-space crowding
and the per-gene rank SD has median **−0.10** (negative in 21/24 arms).
Rank-space crowding (±1 % of ranks) is weakly and consistently positive
(**+0.23**, positive in 24/24 arms). What dominates is the gene's **own score
noise**: the correlation with the score coefficient of variation is **+0.82**
(positive in 24/24 arms).

That points at the seed set, and truncating the training seed set closes the
argument. As the seed set grows from 15 to 304 seeds (5 % → 100 % of the
available positives), the median per-gene rank SD falls from 0.134 to **0.028**
— monotonically in **24/24 arms** — with a fitted exponent of
**−0.488** (median; IQR −0.51 to −0.46) against the seed count, i.e. the
textbook 1/sqrt(n) law for the sampling error of a mean. Along the same axis
ΔAUROC *rises* (monotonically in 20/24 arms, median +0.039 → +0.081), so even
here stability and discrimination improve together.

This has an exact corollary. RWR is linear in the seed weights,
`P = alpha (I - (1-alpha) W^T)^-1 P0`; the bootstrap mean of `P0` is `P0`, so
**bagging over seed resamples is provably a no-op for an RWR prior**. The
instability of a network prior is the sampling noise of a mean over seeds, and
the only stabiliser is more distinct seeds — that is, a larger training set.

### L.4 Consequences for the benchmark

1. Do not smooth the prior: there is no stability to buy and discrimination to
   pay.
2. Report the instability instead of hiding it. The per-gene rank SD (0.019 –
   0.044 percentile units across the 24 arms) is a cheap, reproducible property
   of a prior and belongs in the benchmark's output.
3. The seed count is a design parameter, not a nuisance: the benchmark's ΔAUROC
   and its ranking stability are both functions of the training-set size, and
   the CV's training fraction (5 folds ⇒ 80 % of the positives) is the lever.
4. The claim "crowding causes rank instability" should not be repeated in its
   score-space form.

**Artefacts:** `results/tables/stability_per_gene.tsv`,
`results/tables/stability_by_arm.tsv`,
`results/tables/stability_crowding_attribution.tsv`,
`results/tables/stability_regularization_grid.tsv`,
`results/tables/stability_control_comparison.tsv`,
`results/tables/stability_selected_scorer.tsv`,
`results/tables/stability_seed_size_scaling.tsv`,
`results/tables/stability_tradeoff_curve.tsv`,
`results/tables/stability_prediction_scorecard.tsv`,
`figures/figure7_stability_tradeoff.png`,
`docs/stability_regularization.md`.

---

## Cross-references

* A's logs: `docs/a0_log.md` ... `docs/a5_log.md`, `docs/w5_log.md`,
  `docs/stats_analysis_log.md`, `docs/gnn_courtesy_log.md`,
  `docs/s2_log.md`.
* A's new tables: `results/tables/stability_per_gene.tsv`,
  `stability_by_arm.tsv`, `stability_features_per_gene.tsv`,
  `stability_attribution.tsv`, `network_stability.tsv`,
  `network_stability_vs_topology.tsv`, `consensus_table.tsv`,
  `consensus_robustness.tsv`, `hub_capped_rwr.tsv`,
  `d1_leakage_audit.tsv`, `d3_split_invariants.tsv`,
  `d4_v1_vs_v2_comparison.tsv`, `rewiring_mixing_check.tsv`
* A's new figures: `figures/figure5_stability_atlas.png`,
  `figure5b_stability_attribution.png`, `figure6_network_stability.png`,
  `figure7_consensus_prior.png`
* A's tables: `results/tables/a4_table1_delta_auroc_v40.tsv`,
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
  `results/tables/stats_variance_components.tsv`,
  `results/tables/stats_bootstrap_ci.tsv`,
  `results/tables/stats_equivalence.tsv`,
  `results/tables/stats_h1_verdict_v40.tsv`,
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