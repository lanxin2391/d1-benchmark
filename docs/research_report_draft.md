# What a network prior contributes to cancer-driver gene prioritisation: a leakage-controlled benchmark

**Draft research report (SCI-standard manuscript, v1.0-rc3)**

**Authors:** lanxin2391 (Networks & Engine), kakamiku (Labels & Data)
**Pre-registration:** Zenodo DOI [10.5281/zenodo.23121762](https://doi.org/10.5281/zenodo.23121762) (frozen 2026-09-27; amendments 1–5 in `docs/pre_registration_amendments.md`)
**Code and data:** <https://github.com/lanxin2391/d1-benchmark>, tag `v1.0-rc3`
**Reproducibility:** `python scripts/audit_reproducibility.py` → 50 PASS / 0 FAIL; `python scripts/verify_claims.py` → 57 PASS / 0 FAIL
**Date:** 2026-11-11

---

## Abstract

**Background.** Cancer-driver gene prioritisation routinely uses a
protein–protein interaction or co-expression network as a prior, on the
guilt-by-association premise that two genes whose proteins interact are
more likely to share driver status. Random walk with restart (RWR) is the
canonical formalisation. Whether the *specific topology* of a network
contributes to driver-gene prediction **beyond what node degree alone
explains** is a question the field has asked for a decade and answered
without a leakage-controlled design.

**Methods.** We assembled a six-network panel stratified by molecular
provenance (one RNA-derived, three protein-derived, two curated or
literature-derived) and a four-label panel spanning mutation-based and
clinical-evidence curation philosophies. We held the learner fixed —
RWR at α = 0.5, no hyperparameter tuning — so that the quantity of
interest is the network's contribution and not the learner's quality.
The primary endpoint was ΔAUROC = AUROC(RWR) − AUROC(degree-only) over
5-fold × 40-repeat cross-validation (200 folds per cell, 24 primary
cells). Three hypotheses were pre-registered (H1: the prior contributes
beyond degree; H2: the contribution is provenance-invariant; H3:
between-network differences reflect literature attention) together with
seven protocol controls and three added after a systematic read of the
reference corpus. A systematic read of a 41-paper reference corpus then
motivated three further controls — seed-degree stratification, network
size and coverage adjustment, and random-input degeneracy — which were
run and reported without being allowed to alter the pre-registered
decision rule.

**Results.** Across the 24 (network × label) cells at α = 0.5, ΔAUROC
ranged from −0.0125 (Reactome × newly-discovered IntOGen drivers) to
+0.1711 (STRING physical × Open Targets no-literature); 23 of 24 cells
were positive. Both null models were rejected overwhelmingly: 23 of 24
cells beat the degree-matched random-seed null and 23 of 24 beat the
rewired-network null at z > 5. On the two endpoints the published
evidence actually uses — partial AUC below 2 % FPR (Guala & Sonnhammer
2017) and precision in the top 10 % — the margin is smaller but still
positive in 20 of 24 arms (median +0.031 and +0.041 against +0.081 for
the full AUROC), and the ordering of the arms is unchanged
(Spearman ρ = 0.80 and 0.90): the prior's advantage is a
global-ranking advantage, not a top-of-list artefact. The 95 % confidence interval on ΔAUROC
excluded zero in 24 of 24 arms under four of the five interval
constructions and in 20 of 24 under the pre-registered gene-cluster
bootstrap. The margin was flat across the degree tercile of the seeds
(mean ΔAUROC +0.053, +0.061 and +0.054 for low, mid and high seed
degree; mean high-minus-low +0.0010), so it is not a hub artefact; it
survived adjustment for network size and label coverage (23 of 24 cells
positive after adjustment; pooled Spearman ρ between ΔAUROC and log₁₀
network size −0.021, permutation p = 0.93); and it survived a
literature-matched null (mean difference +0.0010). The learner's
ranking collapsed onto an almost unrelated gene list under a random
input (top-100 Jaccard 0.019), confirming that it uses its input rather
than returning a fixed hub set. H2, however, is **not testable** with
this panel: provenance is confounded with network size and coverage, and
the three provenance classes hold three, one and two networks.

**Conclusions.** A network prior contributes to cancer-driver
prioritisation beyond node degree, and that contribution is not
explained by the degree of the seed genes, by network size, or by the
degree distribution of the labels. The single negative cell is a
*coverage failure* — Reactome's network contains only 41 % of the
newly-discovered drivers — not evidence of a negative effect; the
equivalence test places it in the same category as five cells with
bounded positive margins. Whether the contribution depends on the
network's molecular provenance remains open and requires a panel that
separates provenance from size and coverage.

**Keywords:** cancer driver genes; network prior; random walk with
restart; degree baseline; leakage control; benchmark; pre-registration;
gene prioritisation.

---

## 1. Introduction

Cancer is a disease of dysregulated genes, and the first step of any
molecular oncology study is to identify the small set of genes that
drive it. The cancer genomics community has curated several
driver-gene catalogues from independent curation philosophies: mutation
recurrence (IntOGen; Bailey et al. 2018), clinical gene–disease
validity (ClinGen; Rehm et al. 2015) and aggregated target–disease
association (Open Targets; Ochoa et al. 2023). These catalogues
disagree in instructive ways. The drivers discovered between 2020 and
2024, for instance, are largely absent from the Open Targets
target–disease associations, because clinical and literature evidence
takes time to accumulate.

A parallel literature has argued that **network structure** can help
identify driver genes beyond what the catalogues alone provide. The
intuition is the guilt-by-association principle: a gene close to many
known drivers in a protein–protein interaction (PPI) network is itself
more likely to be a driver. This has been formalised as random walk with
restart (RWR; Kohler et al. 2008; Vanunu et al. 2010), as label
propagation (Zhang et al. 2018), and as a family of graph neural network
methods (Boiarsky et al. 2021, 2023).

The field has converged on three empirical observations. First, that RWR
and its relatives beat a degree-only ranking. Second, that different
networks give different absolute discriminative performance — FunMap
outperforms STRING, STRING outperforms Reactome — and that this
0.25-AUROC gap is the basis of the "best network" question. Third, and
critically, that **no published study has measured whether the network
prior contributes anything beyond degree in a leakage-controlled
setting**. Existing benchmarks either fix a network without comparing
against a degree baseline (Picart-Armada et al. 2019; Hill et al.
2019), compare networks without reporting the margin over degree
(Zhang et al. 2021), or use overlapping train/test labels that confound
the degree signal with label leakage (Chu et al. 2023; DriverNet,
Bashashati et al. 2012; DawnRank, Hou & Ma 2014).

That gap is not merely a matter of a missing comparison. The published
evidence on the confounders is unusually specific, and each of the
following has been shown to threaten exactly the quantity we want to
measure:

* **Degree confounding is real and severe.** Driver genes are hubs:
  the median out-degree of Cancer Gene Census genes is 12 against 3 for
  non-CGC genes (Hou & Ma 2014), and the enrichment is largely an
  artefact of research intensity — bait usage correlates with degree at
  r = 0.52, and correcting for it removes the cancer-versus-control
  degree difference (Schaefer et al. 2015). Zhang et al. (2021) showed
  directly that RWR's apparent performance tracks the **degree
  centrality of its seed nodes**, warning that the advantage "may look
  better simply because its seed nodes have higher degree centrality".
* **Network size and label coverage confound the comparison.** Over 45
  human interactomes, performance tracks network size, and the effect
  largely collapses once gene-set coverage is controlled (Wright et al.
  2025).
* **A learner that ignores its input returns the hubs.** Over 10,000
  random gene lists, algorithms differ enormously in how much of their
  output is input-independent (Hill et al. 2019).
* **Complex leakage inflates the margin.** Holding protein complexes
  together in the cross-validation fold drops a strong network classifier
  from roughly 12 to fewer than 4.5 hits in the top 20 (Picart-Armada et
  al. 2019).
* **The choice of endpoint changes the ranking of methods.** AUROC and
  top-k give different orderings (Picart-Armada et al. 2019), and a
  degree-controlled local method can beat RWR on a full ranking while
  RWR wins only in the top 1–10 % (Guala & Sonnhammer 2017).

The D1 benchmark is designed to answer the question these papers leave
open. Its design principles are: **a fixed learner** (so the finding is
about the network, not the method); **a panel stratified by molecular
provenance** (so the network's origin is an experimental factor rather
than a nuisance); **a pre-registered endpoint and decision rule**
(including the resampling unit, the equivalence margin and the
multiplicity family, so that the analysis cannot be re-tuned after the
fact); **released split files** rather than regenerated seeds (Boutros
et al. 2014); and **a control stack** that removes each alternative
explanation in turn, including the three that the literature identifies
and that a standard design does not cover.

---

## 2. Materials and methods

### 2.1 Network panel

Six primary networks spanning three provenance classes, plus two
STRING threshold sensitivity arms. The harmonisation pipeline drops
unmapped genes, removes self-loops, orders each edge alphabetically,
drops duplicates and restricts each network to its largest connected
component (LCC). The primary analysis is unweighted (decision D-01).

| net_id | provenance class | source | nodes (pre-LCC) | LCC nodes | LCC edges |
|---|---|---|---:|---:|---:|
| `funmap` | protein-derived | FunMap v1.0, CPTAC pan-cancer proteogenomics (Shi et al. 2025) | 10,525 | 10,189 | 196,605 |
| `intact` | protein-derived | IntAct, human–human, taxid 9606 (del Toro et al. 2022) | 17,918 | 17,918 | 564,706 |
| `string_phys700` | protein-derived | STRING v12, physical subnetwork, combined_score ≥ 700 | 10,746 | 9,830 | 85,576 |
| `rna_coexp` | RNA-derived | Built from FunMap's published CPTAC mRNA matrices, density-matched | 9,928 | 9,152 | 196,172 |
| `string_full700` | curated / literature | STRING v12, all-evidence channels, combined_score ≥ 700 | 16,201 | 15,882 | 236,712 |
| `reactome` | curated / literature | Reactome v86, human only (Milacic et al. 2024) | 4,993 | 3,953 | 20,143 |
| `string_full400` | sensitivity | STRING v12, combined_score ≥ 400 | 19,488 | 19,486 | 929,471 |
| `string_full900` | sensitivity | STRING v12, combined_score ≥ 900 | 12,323 | 11,693 | 100,383 |

**Provenance overlaps, stated because they bound the interpretation.**
FunMap's own gold standard is Reactome-derived: a gene pair is labelled
functionally associated when the two genes share a "detailed" Reactome
pathway of ≤ 200 proteins (Shi et al. 2025). Reactome is also an input
to STRING's `database` channel. The "curated/literature" class is
therefore one upstream knowledge base reaching the panel by two routes,
and the "protein-derived" class contains a network built on that same
knowledge base. This is one reason the provenance hypothesis cannot be
settled by this panel (§2.6, §4.3).

### 2.2 Label panel

| label_id | catalogue | positives |
|---|---|---:|
| `intogen2024` | IntOGen 2024.09.20 (Martínez-Jiménez et al. 2020; the compendium of Bailey et al. 2018) | 633 |
| `ot_all` | Open Targets Platform 26.06, all-evidence | 600 |
| `ot_nolit` | Open Targets Platform 26.06, one literature extractor removed | 600 |
| `intogen_temporal_new` | IntOGen 2024.09.20 \ IntOGen 2020.02.01 | 152 |
| `clingen_label` (optional) | ClinGen gene–disease validity (Rehm et al. 2015) | 84 |

All gene identifiers are HGNC-approved symbols; previous symbols and
aliases are used only when they resolve to exactly one approved symbol,
and ambiguous identifiers are dropped and counted (decision D-02). Label
files contain positives only, ranked by descending evidence with ties
broken by symbol.

Two properties of the label panel matter for the interpretation of
§4.4. First, the two Open Targets label sets share **430 of 600**
positives (72 %), so `ot_all → ot_nolit` is close to a re-analysis of
the same label and should not be read as a catalogue-transfer result.
Second, the `intogen_temporal_new` positives are almost disjoint from
the Open Targets positives (24 of 600 and 18 of 600): the newly
discovered drivers are precisely the genes the target–disease
aggregation has not yet picked up.

**`ot_nolit` is not a literature-free label.** It is Open Targets 26.06
with the score of a single evidence extractor removed — the `europepmc`
text-mining extractor (Ochoa et al. 2023). The platform's molecular
interaction widget re-packages STRING, IntAct, Reactome and SIGNOR, and
its other evidence channels (genetic association, somatic mutation,
known drug, expression, animal model) are unaffected. The difference
between the two label sets therefore measures the sensitivity of the
Open Targets score to **one** text-mining extractor, not the
contribution of literature as a whole. This paper's H3 argument rests on
Control 7, which controls literature volume on both sides of the
comparison, and not on `ot_all` versus `ot_nolit`.

### 2.3 The learner (locked)

For each (network, label) cell and each α ∈ {0.3, 0.5, 0.7}:

1. `W = D⁻¹A` is the row-stochastic transition matrix of the network LCC.
2. `P₀` is the column-stochastic matrix whose *j*-th column is uniform
   over the train-fold positives of the *j*-th (repeat, fold) pair.
3. `P_{k+1} = α P₀ + (1 − α) Wᵀ P_k` until the worst column changes by
   less than 1e-8, or 200 iterations.
4. Genes are scored by `P_k[:, j]`; AUROC is computed on the test-fold
   positives against the rest of the universe, with the seeds excluded
   from the test set (decision D-06).

The learner is never tuned. α = 0.5 is the primary value; 0.3 and 0.7
are sensitivity arms. The implementation is covered by 13 unit tests
including fixed-point agreement with the analytic solution
`P = α(I − (1 − α)Wᵀ)⁻¹P₀`.

### 2.4 Cross-validation

5-fold × 40 repeats = 200 (repeat, fold) pairs per cell, with
`StratifiedKFold(shuffle=True, random_state=20260926 + repeat)`. The
40-repeat design (rather than the originally pre-registered 10) was
adopted after a coverage simulation showed the repeat-level interval to
be the calibrated construction (amendment 2); split files are released
as artefacts, not regenerated from seeds.

### 2.5 Controls

| # | control | manuscript label | what it removes | implementation |
|---|---|---|---|---|
| 1 | degree-only baseline | Control 1 | the degree sequence itself | `d1/engine/baselines.py` |
| 2 | degree-matched random-seed null | Control 2 | "the seeds are just hubs" | `scripts/run_a5.py` |
| 3 | configuration-model rewiring null | Control 3 | "any edge set of this degree sequence would do" | `scripts/run_rewiring.py` |
| 4 | shared gene universe | Control 4 | "the networks differ because they cover different genes" | `scripts/run_d12_shared.py` |
| 5 | cross-catalogue label transfer | Control 5 | "the catalogue is the answer, not the network" | `scripts/run_transfer_all.py` |
| 6 | temporal split | Control 6 | "the prior only knows the drivers of 2020" | `scripts/run_d09.py` |
| 7 | publication-count adjustment with a literature-matched null | Control 7 | "popular genes are popular drivers" | `scripts/build_pubcount_adjusted.py`, `scripts/run_control7_matched_null.py` |
| 8 | complex-aware cross-validation | Phase D | "held-out drivers inherit their complex partners' signal" | `scripts/run_d3_generate_cx_splits.py` |
| 9 | hub-capped and degree-penalised RWR | E3 | "the prior is degree plus a hub effect" | `scripts/run_hub_capped_rwr.py` |
| 10 | seed-degree stratification (Zhang et al. 2021) | **Control 8** | "the real seeds are hubs" | `scripts/run_control8_seed_degree.py` |
| 11 | network size and coverage adjustment (Wright et al. 2025) | **Control 9** | "bigger networks and better coverage win" | `scripts/run_control9_size_adjustment.py` |
| 12 | random-input degeneracy (Hill et al. 2019) | **Control 10** | "the learner ignores its input and returns the hubs" | `scripts/run_control10_random_input.py` |

(The two numbering systems are kept apart deliberately: the table's first
column is the *position* in the control stack, and the "manuscript label"
column is the name the manuscript, the abstract and the amendment log use.
The three literature-mandated controls are **Controls 8, 9 and 10**
everywhere except this table's position column.)

Controls 8–10 (the three shaded in bold) were specified **after** the
confirmatory runs existed, in response to a systematic read of the
41-paper reference corpus. They are registered in the exploratory family
(amendment 5) and cannot alter the pre-registered verdict, but they are
reported whatever they show.

**Control 8 (seed-degree stratification)** splits the training-fold
positives of every arm into three equal degree terciles and recomputes
the primary endpoint with the seeds drawn from one tercile only. The
test set, the degree baseline and the fold structure are identical
across the three strata, so the three numbers are directly comparable.

**Control 9 (size and coverage adjustment)** regresses the per-cell
endpoint on the label set (a fixed effect), log₁₀(LCC edges) and the
label coverage, and reports the fitted value as the size- and
coverage-adjusted endpoint; the size correlation is additionally tested
against a 10,000-permutation null.

**Control 10 (random-input degeneracy)** runs the learner twice on the
same fold universe — once with the real training-fold positives, once
with a random seed set of the same size drawn from the network's genes
**excluding** the fold's own test genes — and reports Hill et al.'s
degeneracy statistic, the hub concentration of the top-K under both
regimes, and the input sensitivity (Spearman and top-K Jaccard between
the two rankings of the same fold).

### 2.5b Co-primary endpoints

The published evidence does not use the full AUROC. Guala & Sonnhammer
(2017) compare AUCs only up to an FPR of 0.02 and report that the
diffusion methods win only in the **top 10 %** of the ranking; Picart-
Armada et al. (2019) find that AUROC and top-k order the methods
differently. We therefore compute two further endpoints on the same
folds, the same seeds and the same degree baseline — the partial AUC
below 2 % FPR, normalised so that a random scorer gives 0 and a perfect
one 1, and the precision in the top 10 % — and report them as
**co-primary**. The pre-registered endpoint remains the full AUROC and
the H1 verdict is unchanged; the repeat is the unit of record
throughout (`scripts/run_pauc_top10.py`,
`results/tables/pauc_top10_summary.tsv`).

### 2.6 Statistics

The primary endpoint is ΔAUROC at α = 0.5 over the 200 (repeat, fold)
pairs of each cell. The **pre-registered decision rule** is: H1 holds
if and only if the 95 % bootstrap confidence interval on ΔAUROC excludes
zero in **every** (network, label) arm.

Five interval constructions are reported side by side: the fold-level
percentile bootstrap, and the repeat-level percentile, *t* and BCa
intervals, plus the gene-cluster bootstrap of pre-registration §9.3, in
which genes are the resampling unit (the gene draw is shared across
folds, so a gene moves as a block). The repeat-level percentile is
designated primary because a coverage simulation showed it to be the
calibrated construction; the gene-cluster interval is the
pre-registered construction and is reported as the conservative bound.

Multiplicity across the 24 cells is controlled with Benjamini–Hochberg
FDR at q = 0.05, applied to the cell-level one-sided p-values.
Equivalence is tested with two one-sided tests (TOST) at a margin of
±0.05 in ΔAUROC units, per Schuirmann (1987): the verdict "equivalent
to zero" requires both one-sided p-values below 0.05, which is
equivalent to the 90 % interval lying entirely inside [−0.05, +0.05].
Between-arm variance is partitioned by a balanced random-effects
analysis into network, label, network × label and fold components.

### 2.7 Reproducibility

`scripts/audit_reproducibility.py` verifies the schema and referential
integrity of every artefact, checks the SHA-256 fingerprints against the
release manifest and runs an end-to-end smoke arm (50 PASS / 0 FAIL at
the time of writing). `scripts/recompute_ground_truth.py` recomputes
every headline number from the released tables into
`review/ground_truth.tsv`, and `scripts/verify_claims.py` checks the
manuscript against that table (57 PASS / 0 FAIL). The second mechanism
was added after a review pass found eleven places where the manuscript
prose had drifted from the released tables; it is a mechanical guard
against that class of error recurring.

---

## 3. Results

### 3.1 The primary endpoint

Across the 24 (network × label) cells at α = 0.5, ΔAUROC ranged from
**−0.0125** (Reactome × `intogen_temporal_new`) to **+0.1711** (STRING
physical × `ot_nolit`); **23 of 24 cells were positive** (Fig. 2,
Table 1). The median cell was +0.0809 and the mean +0.0787.

Per label, the mean ΔAUROC over the six networks was +0.0821
(`intogen2024`), +0.1058 (`ot_all`), +0.1062 (`ot_nolit`) and +0.0208
(`intogen_temporal_new`); the optional ClinGen label gave +0.1149.
The temporal label is the weakest by a factor of four, and the rest of
the paper is largely an account of why.

### 3.2 The two null models

**Degree-matched random-seed null (Control 2).** 100 degree-matched
random seed sets per cell, each propagated over the 40 repeats. The
observed ΔAUROC beat the null in **all 24 cells**: 23 cells at cell-level
z > 5 (8.2 ≤ z ≤ 226.9) and the remaining cell at z = +3.09
(p = 1.0 × 10⁻³). The null mean was itself ≤ 0 in 23 of 24 cells: with
degree-matched random seeds, RWR is *worse* than the degree baseline.
The observed margin is therefore not a degree artefact.

**Rewired-network null (Control 3).** 598 degree-preserving edge
permutations (100 per network, 98 for Reactome; the released table
records the count per network). The observed ΔAUROC exceeded the
rewired null in **23 of 24 cells** at z > 5 (up to z = 205.2). The single
exception was the Reactome temporal cell, where the observed value fell
*below* the rewired null (z = −3.53).

Both nulls say the same thing: the margin cannot be produced by the
degree sequence, nor by any edge set with that degree sequence.

### 3.2b The margin on the endpoints the literature uses

The full AUROC is the metric on which propagation methods look worst, so
the same 24 arms were re-analysed on the partial AUC below 2 % FPR and
on the precision in the top 10 % (`results/tables/pauc_top10_summary.tsv`):

| endpoint | median margin over degree | arms positive |
|---|---:|---|
| full AUROC (pre-registered primary) | **+0.0809** | 23/24 |
| partial AUC, FPR ≤ 0.02 | **+0.0308** | 20/24 |
| precision in the top 10 % | **+0.0412** | 20/24 |

The prior's advantage is a **global-ranking advantage**: it is 38 % (partial
AUC) to 51 % (top-10 % precision) smaller in those regions than in the
full AUROC, and 4 of the 24 arms turn negative on each. The *ordering* of
the arms is nevertheless stable — Spearman ρ = 0.80 (AUROC vs partial AUC)
and 0.90 (AUROC vs top-10 % precision) — so the endpoint choice changes
how large the margin is, not which networks look best.

The sharpest single statement is `rna_coexp × intogen2024`: ΔAUROC
+0.0521 [+0.0504, +0.0538] but Δpartial-AUC −0.0072 [−0.0088, −0.0055],
i.e. the prior is *worse* than degree below 2 % FPR while being clearly
better over the whole ranking. The negative cell (Reactome ×
newly-discovered drivers) is negative under all three endpoints, so that
finding does not depend on the metric.

### 3.3 The margin is not a hub effect, and not a seed-degree effect

**Hub-capped RWR (E3).** Removing the hubs from the propagation
while keeping the degree baseline untouched left the median margin
unchanged (+0.0849 against a baseline of +0.0809); the degree-penalised
walk *increased* it (+0.0979). But silencing the hubs rather than
removing them took the median effect to 31 % of its size, so roughly
two thirds of the raw margin sits in the hubs' own ranking position.
Both numbers are reported and the smaller one is quoted.

**Seed-degree stratification (Control 8).** Recomputing the endpoint
with the seeds drawn from the low, mid or high degree tercile of the
training positives gave mean ΔAUROC of **+0.0527, +0.0609 and +0.0537**
respectively — a mean high-minus-low gap of **+0.0010**, one fiftieth of
the margin itself, with the sign of the gap inconsistent across the
panel (negative in 11 of 24 arms). The margin does not depend on the
degree of the seeds. This is the direct answer to Zhang et al.'s
warning, and it is the opposite of what a degree-driven artefact would
look like.

**Random-input degeneracy (Control 10).** Under a random input the
learner's ranking collapsed onto an almost unrelated gene list: the
top-100 sets of the real-input and random-input rankings of the same
fold overlapped by **1.9 %**, and only **0.6 genes** persisted in the
top 100 across independent random runs, against 123.6 under the real
input. For scale, Hill et al. reported 239 for ToppNet-HITS, 200 for
SigNet, 4 for Random Walk and 0 for Network Propagation. The learner
uses its input.

### 3.4 The margin is not a network-size effect — but coverage matters

**Size and coverage adjustment (Control 9).** A bigger network does
not give a bigger margin: the pooled Spearman correlation between ΔAUROC
and log₁₀(LCC edges) was **−0.021** (10,000-permutation p = 0.93), and
within a label set the six networks gave correlations of +0.26, +0.26,
−0.31 and −0.43 (all p > 0.39). After regressing on the label set,
log₁₀ edges and the label coverage, the adjusted endpoint was positive
in **23 of 24 cells** — the same count as the raw endpoint — with the
sign preserved in 22 of 24.

Two of the model terms deserve attention because they cut in opposite
directions. `log₁₀(LCC edges)` was significant (coefficient −0.0393,
p = 0.013, partial R² = 0.30) but had the *wrong* sign for a size
artefact: larger networks did not do better. **Coverage** was the
strongest single covariate (coefficient +0.1449, p = 8.6 × 10⁻⁴, partial
R² = 0.47): an arm is easier when more of its positives are present in
the network. This is the mechanism behind the Reactome temporal cell
(§3.7) and it is why the provenance hypothesis cannot be settled by this
panel (§3.8).

### 3.5 The controls that equalise the design

| control | what it equalises | result |
|---|---|---|
| Shared gene universe (Control 4) | the evaluation universe, across all six networks | 1,605 shared genes, 321 test genes per fold; ΔAUROC −0.036 to +0.130, 22/24 positive |
| Cross-catalogue transfer (Control 5) | the catalogue the prior was trained on | mean ΔAUROC +0.0313, **39.8 %** of the within-label magnitude; 59 of 66 cells positive |
| Temporal split (Control 6) | the recency of the drivers | ΔAUROC −0.0095 to +0.0923; 5 of 6 networks positive, on 50 independent training-seed draws per arm |
| Complex-aware CV (Phase D) | the protein-complex leakage channel | 0 of 24 arms change sign; median absolute relative change 0.142; the prior *loses* ground (median signed change −12.3 %) |
| Common-size cap (D-11) | the number of positives per arm | 22 of 24 arms positive after capping every network to its per-network minimum |
| Publication-count adjustment (Control 7) | the literature volume of the positives and of the background | the literature-matched null leaves the margin essentially unchanged (mean difference **+0.0010**); the matched endpoint is positive in 28 of 30 arms |

The cross-catalogue retention of 39.8 % is the honest headline for
Control 5: the prior carries **partial** catalogue-independent signal.
The floor is set by the `intogen_temporal_new` target, whose positives
are almost disjoint from the Open Targets catalogues.

The complex-aware control is the one that argues *against* the
headline, and it is reported as such: the standard stratified CV inflates
the effect by a median of 12 %, because a held-out driver's score is
partly determined by a complex partner that sits in the training fold.
The verdict survives without the leak, and the shrinkage is itself the
finding.

### 3.6 The variance decomposition

A balanced random-effects partition of the between-arm variance gives
**40.1 % to the network × label interaction**, **34.8 % to the label
set**, **24.8 % to the network** and **0.3 % to fold noise**. The
interaction is the largest component: which network helps most depends
on which catalogue you are trying to predict. That is the quantitative
form of the "best network" question, and it is why the panel cannot
answer a provenance question.

### 3.7 The one negative cell is a coverage failure

`reactome × intogen_temporal_new` (ΔAUROC = −0.0125) is the only
negative cell and the only cell that falls below the rewired null. It
is **not** a demonstration of a negative effect:

* The equivalence test (TOST, margin ±0.05) places it in the same
  category as five cells with *positive* but bounded margins. "Equivalent
  to zero" means there is no evidence of an effect — in either
  direction. Six of the 24 cells are in that category, not one.
* Its 95 % interval excludes zero under four of the five interval
  constructions and **contains zero under the pre-registered
  gene-cluster bootstrap** ([−0.0589, +0.0375]). No construction
  supports reading the cell as evidence of a negative effect.
* Reactome's LCC contains only **63 of the 152** newly-discovered
  drivers (41 %), the lowest coverage in the panel. The arm is measured
  on a universe that largely lacks the genes it is asked to find.
* The coverage term is the strongest covariate in the size model
  (§3.4), which is the quantitative version of the same point.
* The cell survives every control that removes a competing explanation:
  the temporal split (−0.0095), the common-size cap (−0.0125), the
  shared universe (−0.0360, on 5 test positives) and complex-aware CV
  (−0.0178, *more* negative). It is null in every seed-degree stratum
  as well.

The remedy, if the failure is a coverage failure, is not a better
learner but a better ensemble: a stale network can only be out-voted.
The consensus prior of §I.5 confirms this — on the temporal label the
median single network is negative (−0.0110) while the consensus is
+0.0358, and removing Reactome from the consensus raises it to +0.0514.
That is the first quantitative statement of the price of curation lag at
the level of an ensemble.

### 3.8 Provenance: an untested hypothesis

**H2 is not tested by this design, and the paper does not claim it.** The
provenance ablation has three classes holding three, one and two
networks. A Kruskal–Wallis test on such a design cannot return a
p-value below 0.1173, the smallest value attainable over all 3! = 6
assignments of the six networks to the three classes. The published
p = 0.50 is the absence of evidence, and the absence is *structural*
rather than empirical: with this panel the hypothesis could not have
been rejected even if the classes were perfectly separated.

Control 11 quantifies why. Coverage is the strongest single covariate
in the panel (partial R² = 0.47) and provenance is correlated with it
— Reactome is both the smallest network and the least curated. Testing
H2 properly needs a panel in which the provenance classes are separated
from network size and coverage: more networks per class than three, one
and two. That is a follow-up, not a result of this paper.

The descriptive pattern, which we report as description and not as
evidence, is that the RNA-derived network is the weakest class on every
label.

---

## 4. Discussion

### 4.1 What the benchmark establishes

The headline is a negative-space claim with a positive-space
consequence: **a network prior contributes to cancer-driver
prioritisation beyond node degree**, and the contribution survives every
control we ran.

What makes the claim stronger than the published alternatives is not the
effect size but the number of alternative explanations that have been
removed:

* It is not the **degree sequence** — the rewiring null holds the
  degree sequence fixed and the margin survives (23/24 cells, z up to
  205).
* It is not **degree-matched seeds** — the degree-matched random-seed
  null has a mean ΔAUROC ≤ 0 in 23 of 24 cells, while the observed
  margin is positive in all 24.
* It is not **the hubs** — removing them from the propagation leaves the
  margin unchanged, and stratifying the seeds by degree leaves it flat
  (mean high-minus-low +0.0010).
* It is not **the size of the network** — the correlation with log₁₀
  edges is −0.021 (permutation p = 0.93) and the size-adjusted endpoint
  is positive in 23/24 cells.
* It is not **literature volume** — the literature-matched null leaves
  the margin unchanged (mean difference +0.0010).
* It is not **the label catalogue** — cross-catalogue transfer retains
  39.8 % of the within-label magnitude, and the temporal split, the
  common-size cap, the shared universe and the complex-aware split each
  leave 22–24 of 24 arms positive.
* It is not **a degenerate learner** — the ranking under a random input
  overlaps the real-input ranking by 1.9 % in the top 100, so the learner
  is genuinely reading its input.

### 4.2 What it does not establish

**Provenance invariance (H2) is untested.** §3.8. This is the most
important limitation and it is a design limitation, not a data
limitation: a three-networks-per-class panel with a 4.5-fold size range
cannot separate provenance from size and coverage.

**The publication-bias story is label-specific.** Control 7 shows the
margin surviving a literature-matched null on the IntOGen and Open
Targets labels (intogen2024 keeps 71 % of its original margin), but
*not* on the temporal label, where trimming the publication-heavy
positives raises the effect and the matched background then removes it.
The matched comparison there is small — 76 positives against 76 matched
negatives — so the right reading is "no longer clearly positive once
literature volume is controlled on both sides", not "exactly zero".
This is a property of the *label*, and it competes with the coverage
reading of the Reactome cell; the two explanations are not the same, and
we report both.

**The advantage is global, and it is smaller where the literature
looks.** On the endpoints Guala & Sonnhammer (2017) actually use — the
partial AUC below 2 % FPR — and on the top-10 % precision, the margin
is 38–51 % smaller than the full-AUROC margin, and 4 of the 24 arms
turn negative (§3.2b). This is the most important qualification on the
headline number, and it is quantified rather than asserted. Two things
follow. The claim survives: the median margin is still +0.031 (partial
AUC) and +0.041 (top 10 %) against 20 of 24 arms positive, and the
ordering of the arms is unchanged (ρ = 0.80 and 0.90), so a reader who
cares about the top of the list still gets the same answer about which
network to use. And the claim is narrower than "the prior helps most
where it matters": it helps *globally*, and it helps least in the
extreme tail of the ranking, where degree is already a strong baseline
and a handful of hubs carry most of the signal. A reader who needs the
top of the list should know that before quoting +0.08.

**One learner family.** The design fixes RWR and measures the network's
contribution under it. The GNN courtesy arm (§J.4) shows that an
unsupervised embedding also beats degree in 19 of 24 cells, but reaches
only about a third of RWR's margin (+0.028 against +0.079). That arm
uses **p = q = 1.0**, which per the appendix of Hamilton et al. (2017)
is equivalent to DeepWalk (Perozzi et al. 2014); we therefore name it a
DeepWalk-style unsupervised embedding rather than claiming a node2vec
implementation we did not make, and we note explicitly that the
feature-using graph neural networks (GCN, GAT, GraphSAGE) are **not**
represented by it and are untested here on featureless graphs. The
node2vec original is also absent from the reference corpus this study
was built from.

**Circularity cannot be eliminated.** The gold standards are literature
products, and so are several of the networks. Bailey et al. (2018) is
the most independent source (its core is derived from mutation
recurrence), and IntOGen is defined by positive-selection statistics
rather than literature co-occurrence; NCG-type catalogues would be the
least independent. The shared upstream source between FunMap and
Reactome (§2.1) is a specific instance. We report the overlap rather
than claim it away, and Control 7 is the control that addresses it.

### 4.3 Implications

**For method developers.** A network prior should not be evaluated in
absolute terms. The margin over a degree baseline is the only
meaningful measure, and it must be reported against the *pre-registered*
resampling unit rather than a post-hoc choice of one — here the choice
of unit moved the H1 count between 20/24 and 24/24, and we report all
five constructions rather than the flattering one.

**For users choosing a network.** The choice among FunMap, STRING
physical, STRING full, IntAct, Reactome and an RNA co-expression network
is largely a wash for mutation-based catalogues; the published 0.25-AUROC
gap between the best and the worst network is mostly a degree effect.
What *does* matter is coverage: a network that lacks the genes you want
to find cannot rank them, and the network × label interaction (40 % of
the between-arm variance) says the best network is catalogue-dependent.

**For the field's benchmarks.** The complex-aware control shows that
standard stratified CV inflates a network prior's measured advantage by
a median of 12 % when protein complexes are split across folds, and
that the inflation is concentrated in sparse and curated networks
(Reactome and FunMap, up to −45 %) while protein-derived networks barely
move. Benchmark suites that use ordinary stratified k-fold on
gene-level labels should expect their network-based methods to look
better than they are.

### 4.4 Limitations

1. **H2 is untested** (structural; §3.8).
2. **Provenance is confounded with coverage**, and coverage is the
   strongest covariate in the panel (§3.4, §3.8).
3. **The `intogen_temporal_new` label is partly a cohort-size slice**
   rather than a pure time slice: driver detection is confounded with
   cohort size (Martínez-Jiménez et al. 2020; Repana et al. 2019), and
   the label's coverage of the newly discovered drivers is the lowest in
   the panel.
4. **The negative cell's D-12 estimate rests on 5 test positives**; it
   localises the failure but does not measure it.
5. **The graph-neural-network comparison is partial** (§4.2).
6. **The cross-catalogue transfer grid has 66 of 72 cells**; the six
   missing cells are reported as missing.
7. **The primary endpoint is ΔAUROC**, and the literature warns that
   this is the metric on which the effect is smallest (§4.2). The
   co-primary low-FPR and top-10 % endpoints quantify the gap (§3.2b):
   the margin is 38–51 % smaller there and 4 of 24 arms turn negative.
8. **No external prospective validation.** Boutros et al. (2014) argue
   that the gold standard must be held out from the developers of the
   methods being assessed. We pre-registered the analysis and released
   the splits, which addresses the reproducibility half of that
   requirement; the independence half would need an external catalogue
   that none of our four labels has been allowed to see.

---

## 5. Conclusion

A network prior contributes to cancer-driver gene prioritisation beyond
what node degree alone explains. The contribution is not a degree
artefact, not a hub artefact, not a network-size artefact, not a
literature-volume artefact, not a catalogue artefact, and not a
degenerate-learner artefact: each of these was removed by a
pre-specified control and the margin survived. The single negative cell
is a coverage failure rather than evidence of a negative effect, and
the provenance question remains open — not because the data are
inadequate but because the panel is too small to separate provenance
from size and coverage.

---

## Data and code availability

* **Pre-registration.** Zenodo DOI [10.5281/zenodo.23121762](https://doi.org/10.5281/zenodo.23121762), frozen at commit `d9c68e5`. Amendments 1–5 in `docs/pre_registration_amendments.md`.
* **Code, split files, result tables and figures.** <https://github.com/lanxin2391/d1-benchmark>, tag `v1.0-rc3`. Release manifest with SHA-256 of 33 tables and 8 figures: `docs/release_manifest_v1.0-rc3.json`; verify with `python scripts/build_release_manifest.py --release v1.0-rc3 --check`.
* **Reproduction.** `conda env create -f environment.yml && conda activate d1 && pip install -e . && pytest` (131 passed, 3 skipped); `python scripts/audit_reproducibility.py` (50 PASS / 0 FAIL); `python scripts/verify_claims.py` (57 PASS / 0 FAIL).
* **Co-primary endpoints.** `results/tables/pauc_top10.tsv`, `pauc_top10_summary.tsv`, `docs/pauc_top10.md`.
* **Per-control tables.** `results/tables/control8_seed_degree_summary.tsv`, `control9_size_adjustment.tsv`, `control9_size_model.tsv`, `control10_random_input_summary.tsv`; the pre-existing controls in `a5_p_values.tsv`, `w5_p_values.tsv`, `d09_table.tsv`, `d11_table.tsv`, `d12_table.tsv`, `d3_split_invariants.tsv`, `d4_v1_vs_v2_comparison.tsv`, `control7_matched_null.tsv`.
* **Ground truth and manuscript verification.** `review/ground_truth.tsv` (242 recomputed quantities) and `review/verification_report.txt`.

## Author contributions

**lanxin2391** — network panel construction, RWR engine, evaluation
harness, controls 1–3 and 8–12, statistical analysis, manuscript
sections A–C and F–J. **kakamiku** — label panel construction, ID
resolution, cross-catalogue and publication-bias controls, ClinGen
integration, manuscript sections D–E and the front matter.

## Ethics statement

Not applicable. The study analyses published molecular interaction data
and published gene–disease associations. No human subjects, no
personal data, and no new wet-laboratory work were involved.

## Competing interests

The authors declare none.

---

### Appendix A. Reference corpus

The manuscript's claims about the published literature are drawn from a
41-paper corpus read in full (methods, thresholds and quantitative
findings extracted per paper into a structured compendium). The
load-bearing citations are:

* **Benchmarking of network methods.** Picart-Armada et al. (2019)
  *PLoS Comput Biol* 15:e1007276; Hill et al. (2019) *PLoS Comput Biol*
  15:e1007403; Guala & Sonnhammer (2017) *Sci Rep* 7:46598; Zhang et
  al. (2021) *Brief Bioinform* 22:bbab006; Wright et al. (2025) *Mol
  Syst Biol* 21:e10077; Bridges et al. (2025) *BMC Bioinformatics*
  26:87; Choobdar et al. (2019) *Nat Methods* 16:843.
* **Cancer-driver methods.** Ghiassian et al. (2015) *PLoS Comput Biol*
  11:e1004120 (DIAMOnD); Guo et al. (2019) *PLoS Comput Biol* 15:e1007520
  (PNC); Erten et al. (2021) *BMC Bioinformatics* 22:62 (BetweenNet);
  Bashashati et al. (2012) *Genome Biol* 13:R124 (DriverNet); Hou & Ma
  (2014) *Genome Med* 6:56 (DawnRank); Chu et al. (2023) *BMC Genomics*
  24:426; Wu et al. (2022) *Interdiscip Sci*; Nourbakhsh et al. (2024),
  review.
* **Statistics and machine learning.** Schuirmann (1987) *J
  Pharmacokinet Biopharm* 15:657 (TOST); Benjamini & Hochberg (1995)
  *J R Stat Soc B* 57:289 (FDR); Perozzi et al. (2014) *KDD* (DeepWalk);
  Grover & Leskovec (2016) *KDD* (node2vec — original absent from the
  corpus); Veličković et al. (2018) *ICLR* (GAT); Kipf & Welling (2017)
  *ICLR* (GCN); Hamilton et al. (2017) *NIPS* (GraphSAGE); Chen &
  Guestrin (2016) *KDD* (XGBoost).
* **Data resources.** Shi et al. (2025) *Nat Cancer* 6:205 (FunMap);
  Szklarczyk et al. (2023, 2025) *NAR* 51:D638, 53:D730 (STRING);
  del Toro et al. (2022) *NAR* 50:D648 (IntAct); Milacic et al. (2024)
  *NAR* 52:D672 (Reactome); Martínez-Jiménez et al. (2020) *Nat Rev
  Cancer* 20:555 (IntOGen); Ochoa et al. (2023) *NAR* 51:D1353 (Open
  Targets); Rehm et al. (2015) *N Engl J Med* 372:2235 (ClinGen); Seal
  et al. (2023, 2026) *NAR* (HGNC); Maglott et al. (2016) *NAR* 44:D733
  (RefSeq); Koutrouli et al. (2024), FAVA.
* **Gold standards and benchmark design.** Bailey et al. (2018) *Cell*
  173:371 (MC3); Repana et al. (2019) *Genome Biol* 20:1 (NCG 6.0);
  Boutros et al. (2014) *Genome Biol* 15:462; Schaefer et al. (2015)
  *Front Genet* 6:260.
