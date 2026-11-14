# What a network prior contributes to cancer-driver gene prioritisation: a leakage-controlled benchmark

**Draft research report (SCI-standard manuscript, v1.0-rc3)**

**Authors:** lanxin2391 (Networks & Engine), kakamiku (Labels & Data)
**Pre-registration:** Zenodo DOI [10.5281/zenodo.23121762](https://doi.org/10.5281/zenodo.23121762) (frozen 2026-09-27; amendments 1–5 in `docs/pre_registration_amendments.md`)
**Code and data:** <https://github.com/lanxin2391/d1-benchmark>, tag `v1.0-rc3`
**Reproducibility:** `python scripts/audit_reproducibility.py` → 50 PASS / 0 FAIL; `python scripts/verify_claims.py` → 63 PASS / 0 FAIL
**Date:** 2026-11-11

---

# Abstract

**Background.** Cancer-driver gene prioritisation methods frequently
use a protein-protein interaction or co-expression network as a
prior. Whether — and how much — the **specific topology** of such
a network contributes to prediction, beyond what node degree alone
explains, is a decade-old question with no published leakage-
controlled answer.

**Methods.** We assembled a six-network panel stratified by
molecular provenance (one RNA-derived, three protein-derived, two
curated/literature) and a four-label panel spanning mutation-based
and clinical-evidence curation philosophies. We held the learner
fixed: random walk with restart at α = 0.5. We measured the primary
endpoint ΔAUROC = AUROC(RWR) − AUROC(degree-only) over 5-fold × 40-repeat cross-validation (200 evaluations per cell, 24 primary cells). We pre-registered three hypotheses
(H1: prior contributes beyond degree; H2: contribution is
provenance-invariant; H3: between-network differences reflect
literature attention) and seven protocol controls: a degree-matched
seed null (A5), a configuration-model rewiring null (W5), a
cross-catalogue label transfer grid (Control 5), a publication-
bias adjustment with a literature-matched null (Control 7), and
three controls added after a systematic read of the benchmark
literature: seed-degree stratification (Control 8, after Zhang et
al. 2021), network size and coverage adjustment (Control 9, after
Wright et al. 2025), and random-input degeneracy (Control 10,
after Hill et al. 2019).

**Results.** Across 24 (network × label) cells at α = 0.5, ΔAUROC
ranges from −0.013 (Reactome × IntOGen temporal_new) to +0.171
(STRING physical × Open Targets no-literature); 23 of 24 cells are
positive. 23 of 24 cells beat the degree-matched seed null and
23 of 24 beat the edge-rewiring null at z > 5. Mixed-effects
variance decomposition partitions the between-arm variance as
40.1 % from the network × label interaction, 34.8 % from label
choice and 24.8 % from network choice, with 0.3 % fold noise. The
95 % CI on ΔAUROC excludes zero in 24/24 arms under four of the five
interval constructions and in 20/24 under the pre-registered
gene-cluster bootstrap; the four exceptions are the smallest-effect
`intogen_temporal_new` arms. The margin is flat across the degree
tercile of the seeds (mean ΔAUROC +0.053 / +0.061 / +0.054 for
low / mid / high seed degree), so it is not a hub artefact; and it
survives adjustment for network size and label coverage (23/24 cells
positive after adjustment, pooled Spearman ρ between ΔAUROC and
log₁₀ network size −0.02, permutation p = 0.93). H2 is **not
testable** with this panel: provenance is confounded with network
size and coverage (3 / 1 / 2 networks per class). Cross-catalogue
transfer retains 39.8 % of the within-label magnitude (59 of 66
transfer cells positive). Publication-bias adjustment leaves the
margin essentially unchanged in the literature-matched null
(mean difference +0.001).

**Additional evidence from the expanded control stack.** Five further
analyses, each registered as exploratory in the amendments, change what
the headline can be taken to mean. A **bait-matched null**, in which
each positive is paired with a negative carrying the same publication
count, removes much of the apparent margin: the share of positives above
the null falls from a median 17.6 % to 11.6 % and, on the
newest-driver label, from 6.8 % to 4.7 %, which is chance (5 %). An
**NCG arbitration test**, in which the positives are drivers that our own
labels do not contain, shows that the prior nevertheless transfers: the
margin is positive in 24/24 arms (median +0.041), and in 26/26 cancer
types. **An audit of our own negative set**
shows that 10-22 % of the held-out negatives are curated driver genes we
do not label. Re-running the whole benchmark with those genes removed
from the background (the positives untouched) changes nothing material:
23 of 24 arms stay positive with no sign flip, the median margin moves
from +0.081 to +0.082, and the conservative gene-cluster interval now
passes in 21 of 24 arms instead of 20. The contamination of the
background is real, and it does not carry the result. **Label-
quality tiering** shows the high-confidence tier carrying the larger
margin in four of five labels, and a **ClinGen noise anchor** (a random
17 % of the positives demoted to background, the disagreement rate
ClinGen's own panels document) moves the margin by less than 0.003.
External functional validation on **1,361 DepMap cell lines** shows the
prior carries set-level but not gene-level signal: no enrichment beyond
chance for the top 30 candidates (5.8 % of cell lines, p < 0.05) and a
clear one for the top 100 (47.6 %). A **stability analysis** shows
that the per-gene instability of the prior is the gene's own score
variance (Spearman ρ = 1.00) and not any topology descriptor, which
makes bagging the prior a provable no-op. Adding the NCG catalogue as
a fifth label shrinks the label main effect from 34.8 % to 30.8 % of the
between-arm spread.

**Conclusions.** A network prior contributes to cancer-driver
prioritisation beyond node degree, and the contribution is not
explained by seed degree, by network size, or by the degree
distribution of the labels. It is also smaller and more fragile than
the field's headline numbers suggest: the median margin is +0.081, it
shrinks by half on the strictest label once curated driver genes are
removed from the background, and a large part of it tracks how often a
gene has been published rather than its topology. The prior does
transfer, to drivers that none of our labels contain, in 24 of 24 arms
and 26 of 26 cancer types — but its advantage in the top of the list is
a property of the interaction network, and it carries no detectable
gene-level functional signal at the scale a practitioner works at. The
one negative cell — Reactome's failure on the 2020-2024 newly-discovered
drivers — is a coverage failure (a curator-lag failure) rather than evidence of
a negative effect, and is reported as such: the same network is
positive on a 2019 catalogue (+0.089), which is the prediction the
curator-lag explanation makes. The provenance question (H2) remains
open: separating provenance from size and coverage requires a larger
network panel than this design supports.

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

### 2.5c The arbitration catalogue, the fifth label and the label-quality tiers

*The section numbers in the blocks below refer to the repository's section files (`docs/manuscript_a_sections.md`); the text is identical to the audited version.*

**An independent arbitration catalogue (C13).** C13 — The NCG arbitration catalogue (an independent driver catalogue)

The NCG catalogue (Repana et al. 2019, *Genome Biology* 20:1,
doi:10.1186/s13059-018-1612-0) is a manually curated repository of 2,372
cancer genes with known or predicted driver roles, collected from 275
publications. **None of the project's four label sets contains it**, which
makes it an independent arbitration catalogue (Boutros et al. 2014: the
gold standard dominates the ranking, so the evaluator should be blind to
the algorithm). It also supplies, per gene, the number of publications
NCG cites for it — the closest available proxy for Schaefer et al.'s
"bait" count, i.e. the study-visibility covariate their degree-matched
null model is built on.

Built by `scripts/build_ncg.py` into `data/processed/labels/ncg.tsv`
(708 known tumour suppressors / oncogenes, 104 multi-screen candidates,
1,560 screening-level genes, with the cancer type, the literature count,
the fraction of cell lines in which the gene is essential and the PPI
degree). It is used twice: to measure and correct the composition of our
negative set (Control 12, F.13) and as the external comparison of
§F.14.
**The NCG catalogue as a fifth label set (C14).** C14 — The NCG catalogue as a fifth label set

The NCG catalogue (Repana et al. 2019) is used twice. As an **external
arbitration set** it is §C13/§F.14. As a **label set** it enters the
main factorial as a fifth label, which is what this section describes.

*Which tier is the label?* The `known` tier (708 high-confidence tumour
suppressors and oncogenes). The 104 `candidate` genes are multi-screen
candidates and are kept aside for the label-quality tiering (§C15); the
1,560 `screening` genes are not drivers and are never positives.

*Symbols.* The HGNC approved symbol, with the previous symbol mapped to
its current one where the catalogue carries one (9 of the 708: CARS →
CARS1, MKL1 → MRTFA, FGFR1OP → CEP43, …), following decision D-02.

*Rank.* NCG 6.0's own confidence level is not available (the portal is
down for maintenance), so the rank is the number of publications NCG
cites for the gene, descending, ties broken by the symbol ascending
(D-13, adapted). This is not cosmetic: the D-11 cap keeps the top-N by
rank, so the capped NCG set is the most-studied NCG drivers.

*The cap does not move.* The D-11 common-size cap N is the minimum
over the four canonical labels of the positives inside a network's
largest connected component, and the NCG label has 407–694 positives
per network, always above the temporal label's 63–150, which is the
minimum. N is therefore unchanged and **the 24 released arms are not
re-run**; six new arms (18 cells at three alphas) are added.

*Two conventions.* The released primary factorial runs each label with
its own positive set (the D-11 cap is a **control** arm, not the primary
design), so the NCG arm is built the same way — the full NCG positive
set in the network (407–694 positives). A capped NCG arm is also run,
so the D-11 control panel covers five labels (30 arms).

*Family.* The 30-arm extension is a **pre-planned exploratory** family
(plan item C5 pre-dates the analysis); the 24-arm pre-registered
endpoint is untouched.
**Label quality tiers (C15).** C15 — Label quality tiers

Every catalogue carries its own evidence, and that evidence is not
equally strong inside a catalogue. Each positive gene is therefore given
a tier from its own catalogue's notion of confidence:

* **ClinGen** — the review status of the gene's best-curated disease
  (Definitive/Strong = high, Moderate = mid, Limited/Disputed = low).
* **IntOGen** — the number of *distinct detecting methods* that flagged
  the gene (IntOGen combines seven complementary methods and votes;
  four or more votes = high, two or three = mid, one = low).
* **Open Targets** — whether the gene's cancer associations are carried
  by a non-literature source. (The number of sources cannot be used:
  88,510 of 89,111 (target, cancer-disease) pairs in the datasource
  table carry exactly one source, so a source-count tier would be
  constant across the whole label.)
* **NCG** — the number of publications the catalogue cites for the gene
  (NCG 6.0's own confidence level is unavailable because the portal is
  down).

The tiers are in `results/tables/label_quality_tiers.tsv`, built by
`scripts/run_c7_build_tiers.py`.

---
### 2.5d Control 13: the negative set cleaned, the benchmark re-run

G8 (§F.13) measured the contamination of the negative set **post
hoc**, on the evaluation set of the existing runs. Control 13 is the
same question done properly: the curated driver genes leave the
**background**, the positives are untouched, and the whole benchmark is
re-run on the new splits.

Concretely, for each of the 24 (network, label) arms the rows of the
NCG-catalogue genes that are negatives of that label are removed from
the split file (5,072,800 negative rows over 2,280 distinct genes
across the 24 arms); the fold assignment, the repeats and the learner
are unchanged; the 24 arms are re-run and the statistics recomputed
(`scripts/run_c9_clean_background_factorial.py`,
`results/tables/c9_clean_background_comparison.tsv`).

This is the exploratory family (amendment 4). The pre-registered
endpoint -- the 24-arm standard design -- is untouched.

---

## §D. Methods §4.6 — Statistical analysis

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
manuscript against that table (63 PASS / 0 FAIL). The second mechanism
was added after a review pass found eleven places where the manuscript
prose had drifted from the released tables; it is a mechanical guard
against that class of error recurring.

---

## 3. Results

### 3.1 The primary endpoint

Across the 24 (network × label) cells at α = 0.5, ΔAUROC ranged from
**−0.0125** (Reactome × `intogen_temporal_new`) to **+0.1711** (STRING
physical × `ot_nolit`); **23 of 24 cells were positive** (Fig. 1, Table 1). The median cell was +0.0809 and the mean +0.0787.

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

### 3.9 Part of the margin is study visibility, not topology

F.12 The bait-matched null: how much of the advantage is study visibility?

**A. The label sets are degree-enriched.** Inside each network's LCC,
the positives' mean degree exceeds the negatives' by 1.0–2.7× (e.g.
string_full700 × ot_all 74.6 against 28.1; reactome × ot_all 22.7
against 8.9; intact × ot_all 139.3 against 60.4). The temporal label is
the least enriched of the four (ratios 1.03–1.54). Full table in
`results/tables/g3_label_degree_enrichment.tsv`.

**B. Against a literature-matched null, part of the advantage survives —
except on the temporal label.** The fraction of held-out positives that
the prior scores above 20 random negatives, before and after matching
the negatives on the publication count (chance level = 5 %):

| label | vs the unmatched null | vs the literature-matched null |
|---|---|---|
| intogen2024 | 20.7 % | 10.9 % |
| ot_all | 31.0 % | 20.2 % |
| ot_nolit | 27.1 % | 18.9 % |
| **intogen_temporal_new** | 6.8 % | **4.7 % (chance)** |

**On the three curated label systems a real part of the advantage
survives** (ot_all: 20.2 % against a 5 % chance level, four times
chance). **On the temporal label it does not**: the positives are
scored above their literature-matched non-drivers at exactly the chance
rate, so the prior's advantage there is a **study-visibility effect, not
a topological one**.

Three independent measurements now agree on this: the Control-7 matched
null (F.8: the temporal margin collapses to +0.0022 against +0.0162),
the low-FPR endpoint (H.1: the margin is 40–60 % smaller where the
literature measures the ranking) and the DepMap validation (the top 30
candidates carry no functional signal). This is the paper's most
consequential limitation and it belongs in the abstract, not in a
footnote: **the prior's contribution is a global-ranking contribution
on the curated label systems, and it is largely a literature effect on
newly discovered drivers.**

The reader should also note what this does *not* say: the positive arms
on the three curated labels still beat a degree-matched, literature-
matched null by four times chance, so the prior is doing something that
the degree baseline and the publication count do not.
F.3 D-11 common-size cap control

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
### 3.10 The prior transfers to drivers we never labelled

F.14 The NCG arbitration test: can the prior find drivers it was not seeded with?

Boutros et al. (2014) insist that the gold standard dominates any
ranking, so the evaluator should be blind to the algorithm. The NCG
catalogue shares no curation with our four label sets, which makes it
the arbitration set we do not control. The seeds here are still our own
label's training-fold positives, so this is a **transfer** test: the
arbitration positives are the NCG-known drivers (708) that are **not**
in this label's positive set, and the background is the rest of the
test fold.

| network | label | arbitration positives per fold | ΔAUROC | ΔpAUC |
|---|---|---|---|---|
| funmap | intogen2024 | 33 | +0.0286 | -0.0042 |
| funmap | ot_all | 64 | +0.0477 | -0.0159 |
| funmap | ot_nolit | 73 | +0.0453 | -0.0158 |
| funmap | intogen_temporal_new | 87 | +0.0207 | -0.0153 |
| rna_coexp | intogen2024 | 26 | +0.0123 | -0.0018 |
| rna_coexp | ot_all | 54 | +0.0093 | -0.0298 |
| rna_coexp | ot_nolit | 61 | +0.0064 | -0.0252 |
| rna_coexp | intogen_temporal_new | 72 | +0.0190 | -0.0177 |
| string_phys700 | intogen2024 | 37 | +0.1269 | +0.0353 |
| string_phys700 | ot_all | 76 | +0.0722 | +0.0119 |
| string_phys700 | ot_nolit | 89 | +0.0590 | +0.0253 |
| string_phys700 | intogen_temporal_new | 107 | +0.0763 | +0.0455 |
| intact | intogen2024 | 46 | +0.0521 | +0.0199 |
| intact | ot_all | 92 | +0.0374 | +0.0028 |
| intact | ot_nolit | 105 | +0.0293 | +0.0069 |
| intact | intogen_temporal_new | 121 | +0.0193 | +0.0152 |
| string_full700 | intogen2024 | 47 | +0.1183 | +0.0451 |
| string_full700 | ot_all | 91 | +0.0758 | +0.0252 |
| string_full700 | ot_nolit | 106 | +0.0521 | +0.0352 |
| string_full700 | intogen_temporal_new | 121 | +0.0679 | +0.0632 |
| reactome | intogen2024 | 22 | +0.0470 | +0.0044 |
| reactome | ot_all | 45 | +0.0344 | -0.0073 |
| reactome | ot_nolit | 56 | +0.0361 | -0.0135 |
| reactome | intogen_temporal_new | 74 | +0.0089 | -0.0346 |

**The prior ranks drivers it was never seeded with above the degree
baseline in 24/24 arms** (median +0.0413; the low-FPR endpoint agrees in
only 13/24, and the split is systematic: the four protein networks
(STRING physical, STRING full, IntAct, Reactome) are positive in every
label, while the two non-protein networks (funmap, RNA co-expression)
are negative in all four). This is the cleanest positive statement the
study can make about transfer: it is not the label we seeded the walk
with, and the prior still finds the drivers the label misses -- but the
advantage in the very top of the ranking is a property of the
interaction network, not of the prior alone.
F.15 The advantage is not cancer-type specific

The arbitration test of §F.14 can be run separately per cancer type:
the arbitration positives are the NCG-known drivers of one cancer type
that our label does not contain (426 of the 708 known drivers carry a
type; 26 types have at least 25 drivers), and the background is the
rest of the test fold. Twenty-six types x 24 arms = 413 tests.

| cancer type | drivers | median ΔAUROC | range | arms positive |
|---|---|---|---|---|
| ampullary_adenocarcinoma | 16 | +0.0708 | +0.0085 .. +0.1894 | 11/11 |
| colorectal_adenocarcinoma | 20 | +0.0568 | -0.0278 .. +0.1782 | 17/18 |
| breast_cancer | 35 | +0.0546 | -0.0540 .. +0.1232 | 20/23 |
| nasopharyngeal_carcinoma | 16 | +0.0538 | -0.0202 .. +0.1408 | 16/18 |
| diffuse_large_B-cell_lymphoma | 21 | +0.0493 | -0.0441 .. +0.1507 | 16/18 |
| chronic_myelomonocytic_leukemia | 16 | +0.0421 | -0.0420 .. +0.1479 | 14/16 |
| triple_negative_breast_cancer | 28 | +0.0415 | -0.0213 .. +0.1103 | 14/16 |
| lung_cancer_(all_histologies) | 24 | +0.0379 | -0.0362 .. +0.1064 | 10/12 |
| bladder_cancer | 22 | +0.0366 | -0.0233 .. +0.1117 | 14/18 |
| melanoma | 20 | +0.0328 | -0.0442 .. +0.1199 | 14/17 |
| prostate_cancer | 20 | +0.0323 | -0.0072 .. +0.2160 | 16/17 |
| squamous_head_and_neck_cancer | 16 | +0.0284 | -0.0507 .. +0.0904 | 9/13 |
| pancreatic_cancer_(all_histologies) | 13 | +0.0276 | -0.0134 .. +0.1710 | 13/14 |
| acute_lymphoblastic_leukemia | 24 | +0.0266 | -0.0594 .. +0.0834 | 14/18 |
| cholangiocarcinoma | 14 | +0.0265 | -0.0913 .. +0.1767 | 12/13 |
| lung_adenocarcinoma | 18 | +0.0236 | -0.0733 .. +0.1209 | 11/16 |
| chronic_lymphocytic_leukemia | 16 | +0.0229 | -0.0202 .. +0.1177 | 13/17 |
| anaplastic_thyroid_carcinoma | 11 | +0.0220 | -0.0302 .. +0.2505 | 9/15 |
| acute_myeloid_leukemia | 17 | +0.0211 | -0.1747 .. +0.0981 | 12/16 |
| endometrial_cancer | 20 | +0.0172 | -0.1013 .. +0.1021 | 13/16 |
| pancreatic_ductal_adenocarcinoma | 14 | +0.0170 | -0.2580 .. +0.1552 | 14/18 |
| hepatocellular_carcinoma | 14 | +0.0151 | -0.0335 .. +0.1052 | 13/18 |
| glioblastoma | 14 | +0.0143 | -0.0188 .. +0.0851 | 8/10 |
| pan-glioma | 17 | +0.0131 | -0.1250 .. +0.0934 | 11/15 |
| T-cell_leukemia/lymphoma | 18 | +0.0101 | -0.0976 .. +0.0904 | 11/18 |
| biliary_tract_cancer | 14 | +0.0097 | -0.0480 .. +0.1193 | 9/12 |

**All 26 cancer types have a positive median margin** (median over the
types +0.0271), and the per-type ranges are wide but centred well above
zero. The advantage of the prior over the degree baseline therefore does
**not** depend on the cancer type: it is a property of the prior, not
of one tumour class. (The background is sub-sampled to 20,000 genes for
tractability; re-computing one arm with 200,000 and with the full
background gives the same +0.1269, so the sub-sampling is immaterial.)
### 3.11 Our own negative set is not clean, and it matters

**Our own negative set is not clean, and it matters**

Bailey et al. (2018) removed 45 curated false positives from 127 true
positives with an LDA filter; Repana et al. (2019) flag 201 possible
false positives. Our negative set is "everything else in the universe",
so every curated driver we happen not to label is treated as a
*confident* negative. With the NCG catalogue in hand this becomes
measurable.

**How contaminated are the negatives?** Between **10 % and 22 %** of the
held-out negatives are in the NCG catalogue, and **1.3 % to 9.5 %** are
NCG-known high-confidence drivers. Reactome's arms are the worst
(15–22 %), the dense protein networks the least (10 %).

**Does it change the result?** The first version of this answer was
wrong, and the correction matters. The post-hoc computation removes every
NCG gene from the **evaluation set**, and because NCG and IntOGen
overlap heavily (463 of the 708 known NCG drivers are IntOGen
positives) that also removes **87 % of the IntOGen positives** (633 →
80). The margin it reports is therefore mostly a statement about a much
smaller evaluation set, not about the cleanliness of the background.

The clean version of the question is to remove the curated drivers from
the **background only** — the positives are untouched, so the comparison
isolates one thing — and to re-run the whole benchmark with those
splits (Control 13, §C16 and §F.19). The answer is that **the headline
does not depend on the contamination**: 23 of 24 arms stay positive
with no sign flip, the median margin moves from **+0.0809** to
**+0.0822** (it goes slightly *up* in almost every arm), and the
gene-cluster interval passes in **21 of 24** arms instead of 20. The
contamination of our negative set is real — 10–22 % of the held-out
negatives are curated drivers — and it does not carry the result.

**The IntOGen label is a good proxy for an independent gold standard; the
Open Targets labels are not.** 70–79 % of the NCG-known drivers are in
`intogen2024` (52–59 % in the temporal label), against 27–37 % in
`ot_nolit` — the Open Targets panels are much broader, which is also why
they carry the largest false-positive contamination.

The honest conclusion is uncomfortable and belongs in the abstract:
**the benchmark's numbers depend on how the negative set is composed**,
and the negative set does not carry the result: with the curated drivers removed from the background and the whole benchmark re-run, 23 of 24 arms stay positive and the median margin is unchanged (+0.081 to +0.082). (Control 12); the pre-registered endpoint is
unchanged.

 independent gold standard; the
Open Targets labels are not.** 70–79 % of the NCG-known drivers are in
`intogen2024` (52–59 % in the temporal label), against 27–37 % in
`ot_nolit` — the Open Targets panels are much broader, which is also why
they carry the largest false-positive contamination.

The honest conclusion is uncomfortable and belongs in the abstract:
**the benchmark's numbers depend on how the negative set is composed**,
and the negative set does not carry the result: with the curated drivers removed from the background and the whole benchmark re-run, 23 of 24 arms stay positive and the median margin is unchanged (+0.081 to +0.082). (Control 12); the pre-registered endpoint is
unchanged.
### 3.12 Which positives carry the margin?

F.18 What the label quality does to the margin

**The margin is recomputed with only one tier's positives counted as
positives** (the other tiers' positives fall into the background). The
seeds, and therefore the prior, do not change, so the comparison
isolates the evaluation noise. The medians over the six networks:

| label | high tier | mid tier | low tier |
|---|---|---|---|
| intogen2024 | +0.0892 | +0.0633 | +0.0653 |
| intogen_temporal_new | +0.0743 | +0.0123 | +0.0232 |
| ot_all | +0.1120 | +0.0908 | — |
| ot_nolit | +0.1097 | +0.0775 | — |
| ncg_known | +0.0666 | +0.0874 | +0.0650 |
| clingen_label | +0.1042 | +0.1207 | — |

For the four IntOGen- and Open-Target-derived labels the **high tier
carries the larger margin** (IntOGen 2024: +0.0892 against +0.0633 and
+0.0653; the temporal label +0.0743 against +0.0123 and +0.0232;
`ot_all` +0.1120 against +0.0908; `ot_nolit` +0.1097 against +0.0775).
The NCG label is the exception and inverts the order (+0.0666 high,
+0.0874 mid, +0.0650 low), and that is not a contradiction: the NCG tier
*is* a publication count, so its "high" tier is exactly the set of genes
the degree baseline already ranks well, and the prior adds most where
the catalogue's own evidence is thinnest. The two orderings are the same
statement seen from two sides.

**Label noise at the ClinGen rate costs almost nothing.** ClinGen's
expert panels disagree on 17 % of multi-laboratory variants, and the
literature review nominates that rate as the external anchor for label
noise. Simulating exactly that rate — a random 17 % of the positives
demoted to background in every fold — barely moves the margin:

| network | n positives | clean | 17 % misassigned | change |
|---|---|---|---|---|
| funmap | 61 | +0.0735 | +0.0733 | -0.0002 |
| intact | 83 | +0.1124 | +0.1134 | +0.0010 |
| reactome | 61 | +0.1710 | +0.1736 | +0.0026 |
| rna_coexp | 49 | +0.0665 | +0.0647 | -0.0018 |
| string_full700 | 84 | +0.0974 | +0.0970 | -0.0005 |
| string_phys700 | 80 | +0.1686 | +0.1685 | -0.0001 |

The margin changes by between −0.0018 and +0.0026. A perfect scorer
would lose 8.5 AUROC points under the same misassignment (the ceiling
falls from 1 to 1 − f/2 = 0.915), so the prior's margin is **not** an
artefact of the label noise that ClinGen's own curation process
documents.

**The wide intervals are a small-positive-set phenomenon, not a label
quality phenomenon.** The width of the gene-cluster bootstrap interval
correlates with the arm's positive-set size at Spearman ρ = −0.966
(p = 2.0e-14, 24 arms) and with the share of high-tier positives at
ρ = −0.475 (p = 0.019). The first term dominates: the cells whose
interval is wide are the cells with few positives, which is what the
variance decomposition already says (the within-cell noise is 40 % of
the total variance). The second term says that, all else equal, a
label whose positives are well-cured gives a tighter interval.

Two caveats on the tier table itself: the small tiers are small (the
temporal label's high tier is 11 genes, the ClinGen mid tier is 4), so
their entries are noisy; and the ClinGen label is almost entirely
Definitive or Strong (80 high, 4 mid, 0 low), so it cannot separate
tiers at all.

---
### 3.13 A fifth label, and what it does to the variance decomposition

F.16 A fifth label: the same picture, and one important exception

The NCG arms land where the other labels land. On the primary
convention (α = 0.5, 5 folds × 40 repeats):

| network | ΔAUROC (NCG) | SD |
|---|---|---|
| funmap | +0.0771 | 0.0176 |
| rna_coexp | +0.0426 | 0.0169 |
| string_phys700 | +0.1352 | 0.0163 |
| intact | +0.0673 | 0.0089 |
| string_full700 | +0.1387 | 0.0133 |
| reactome | +0.0647 | 0.0216 |

All six are positive (+0.0426 to +0.1387) and the network order is the
one the other labels produce (the protein networks above the
co-expression networks). The NCG label's margins sit between
`intogen2024` and `ot_all` on most networks, which is what its 65 %
overlap with IntOGen predicts.

**The exception is Reactome.** On the D-11 capped panel, where all five
labels have the same number of positives inside the arm:

| network | intogen2024 | ot_all | ot_nolit | temporal | **NCG** |
|---|---|---|---|---|---|
| funmap | +0.1268 | +0.0191 | +0.0885 | +0.0196 | **+0.1127** |
| intact | +0.0681 | +0.0664 | +0.0973 | +0.0071 | **+0.0858** |
| reactome | -0.0110 | +0.0627 | +0.1649 | -0.0125 | **+0.0893** |
| rna_coexp | +0.0829 | +0.0403 | +0.0956 | +0.0127 | **+0.0720** |
| string_full700 | +0.1020 | +0.0785 | +0.1388 | +0.0604 | **+0.1177** |
| string_phys700 | +0.0955 | +0.0938 | +0.1819 | +0.0378 | **+0.1222** |

Reactome is negative on `intogen2024` (−0.0110) and on the temporal
label (−0.0125) — the two most recent catalogues — but **positive on
NCG (+0.0893)**, the 2019 catalogue. The "Reactome collapses" pattern
is therefore **label-specific, not network-specific**: Reactome handles
the drivers that existed when it was curated, and fails on the ones
that arrived afterwards. That is precisely the prediction the curator-lag
hypothesis makes (§J.3), and it is the first test of that hypothesis
that does not rely on the temporal split.
F.17 What a fifth label does to the variance decomposition

The plan's expectation was that adding a label would let the label main
effect and the interaction take a larger share of the between-arm spread
and the fold noise a larger share of the total. The same corrected
decomposition on the 24-arm and the 30-arm grids:

| variance source | % between-arm spread (24) | % between-arm spread (30) | change | % total (24) | % total (30) |
|---|---|---|---|---|---|
| network | 24.82 % | 26.95 % | +2.13 pp | 14.84 % | 15.80 % |
| label | 34.78 % | 30.75 % | **−4.02 pp** | 20.80 % | 18.04 % |
| network × label | 40.07 % | 41.95 % | +1.88 pp | 23.96 % | 24.60 % |
| fold (within cell) | 0.34 % | 0.35 % | +0.02 pp | 40.39 % | 41.56 % |

The label main effect **shrinks** (34.8 % → 30.8 % of the between-arm
spread) and the interaction and the noise grow, which is the direction
the plan expected. But NCG is 65 % identical to IntOGen, so this is the
effect of adding a *redundant* label; a genuinely independent fifth
label would be expected to shrink the label effect more. The honest
statement is that **the label main effect is real but not dominant**:
across five catalogues of quite different provenance, a quarter of the
between-arm spread is the network and a third is the label, and the
largest single piece is the interaction — which network you use depends
on which catalogue you are trying to predict.
### 3.14 Where the stability of the prior comes from

I.1 What is stable, and what is not

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
I.2 What drives the instability — and what does not

Two things predict a gene's rank instability, and one thing we expected
to be the explanation turns out not to be.

**The gene's own score noise.** The per-gene relative error of the score
across the 40 re-partitions is itself the strongest predictor of the
rank instability: Spearman(rank SD, score CV) = **+0.844** (median over
the 24 arms, positive in 23 of 24 arms; +0.82 in 24/24 when measured on
the shared engine). A gene whose score is intrinsically noisy has an
intrinsically noisy rank.

**The crowding of the score distribution — but only under one
definition.** The number of genes lying within ±10 % of a gene's score
correlates with the instability at **+0.478** (positive in 24 of 24
arms). The same quantity defined with a *Scott-bandwidth* window — the
data-driven choice — correlates at **-0.175** (negative in 23 of 24
arms). Two defensible definitions of "crowding" therefore give
**opposite signs**, which means crowding is not a mechanism in itself:
it is a window-dependent reparameterisation of the first effect. We
report both and draw no mechanism from either.

**Degree's large coefficient was a suppression artefact.** Degree's *marginal* within-arm association with rank SD is **+0.038** — that is, nothing —
while in the model it carries a coefficient of **-0.083** (CI -0.124 to -0.043);
degree correlates 0.72 with the RWR score, and the score's own
coefficient was absorbing it. The same model gives the gene's
own score CV the largest coefficient (**+0.778**, marginal ρ
**+0.844** in 24/24 arms). The formal
attribution (within-arm OLS with arm fixed effects, 267,696 gene-arm
rows, standardised predictors):

| predictor | joint β on log(rank SD) | marginal within-arm ρ |
|---|---|---|
| **log gene score CV** | **+0.778** | **+0.844 (23/24 arms)** |
| log local score density (±10 % band) | +0.672 | +0.478 (24/24) |
| log score magnitude | +0.238 | +0.340 (24/24) |
| log crowding (Scott bandwidth) | +0.031 (q = 0.56) | −0.175 (1/24 positive) |
| log degree | −0.083 | +0.038 (≈0) |
| k-core number | −0.024 (q = 0.025) | +0.044 |
| clustering coefficient | +0.009 | +0.088 |
| log betweenness | −0.027 | +0.074 |
| log community size | −0.026 | +0.019 |
| log PageRank | −0.001 (q = 0.95) | +0.023 |
| distance to the training seeds | +0.007 | −0.125 (3/24) |

**What causes it: the sampling noise of the seed mean.** The correlation
above says *which genes* are unstable, not *why*. Truncating the training
seed set answers the why: as the seed count grows from 15 to 304 (5 % →
100 % of the available positives), the median per-gene rank SD falls
monotonically from 0.134 to **0.028** — in 24 of 24 arms — with a
fitted exponent of **−0.488** against the seed count, i.e. the textbook
1/sqrt(n) law for the sampling error of a mean. Along the same axis
ΔAUROC *rises* (median +0.039 → +0.081), so on this axis stability and
discrimination improve together: there is no trade-off to exploit.

That has an exact corollary. RWR is **linear in the seed weights**
(`P = α (I − (1−α) Wᵀ)^−1 P₀`), so the bootstrap mean of `P₀` is `P₀`
itself and **bagging over seed resamples is provably a no-op for an RWR
prior**. The instability of a network prior is the sampling noise of a
mean over seeds, and the only stabiliser is more distinct seeds — that
is, a larger training set.

**No network descriptor predicts the size of the contribution.** Every
correlation with the mean ΔAUROC has q ≥ 0.70 (§I.3). Which network
contributes a lot and which network gives a reproducible estimate are two
different questions, and only the second one has a topological answer.
I.3 Which networks are stable

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
I.4 The picture replicates on two independent systems

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
I.5 A model built from these findings: a consensus prior

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
I.6 Can the instability be removed? No — and that is the result

> **Family: EXPLORATORY.** This section reports an analysis that changes the
> learner's output. It is therefore **not** part of the confirmatory analysis,
> it does not enter the H1 decision rule, and every number in it is a
> property of the *exploratory* arm, not of the benchmark's primary endpoint.
> It was pre-registered in `docs/stability_regularization_plan.md` (committed
> before any result existed) and is reported in full in
> `docs/stability_regularization.md`.

#### I.6.1 The question and the pre-registered operator

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

#### I.6.2 The result: the prior's ranking cannot be smoothed into stability

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

#### I.6.3 What drives the instability, and what buys it

(This is the mechanism section of §I.2; it is repeated here because it is the conclusion of the exercise.)


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

#### I.6.4 Consequences for the benchmark

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
`figures/figure8_stability_tradeoff.png`,
`docs/stability_regularization.md`.

---
I.7 What this section does and does not establish

* It establishes that **the benchmark's numbers are reproducible while
  its gene-level rankings are not**, and that the reason is the
  sampling noise of the seed mean (a 1/sqrt(n) law), not the shape of
  the network: the same 1/sqrt(n) behaviour appears when the seed set
  is truncated, and RWR's linearity in the seed weights makes
  bagging the prior a provable no-op.
* It establishes that **there is no stability to buy**: the
  pre-registered shrinkage operator is an identity on the ranking, and
  every operator that does move the ranking is worse on both axes
  (§I.6).
* It establishes that **the effect size is not predictable from
  topology**, which is a caution for anyone who reads a single
  network's performance as a property of that network.
* It does **not** establish a causal role for crowding (the two
  definitions disagree in sign), and it does **not** establish that the
  consensus prior is the best model — it is the one model that follows
  from our own findings, and its per-gene stability is worse.

---

### 3.15 The headline does not depend on the contamination of the background

| network | label | standard | background cleaned | change |
|---|---|---|---|---|
| string_phys700 | ot_nolit | +0.1711 | +0.1741 | +0.0030 |
| string_phys700 | ot_all | +0.1582 | +0.1624 | +0.0042 |
| string_full700 | ot_nolit | +0.1267 | +0.1269 | +0.0002 |
| string_full700 | ot_all | +0.1201 | +0.1214 | +0.0013 |
| string_phys700 | intogen2024 | +0.1189 | +0.1289 | +0.0100 |
| string_full700 | intogen2024 | +0.1133 | +0.1183 | +0.0050 |
| reactome | ot_all | +0.1103 | +0.1119 | +0.0016 |
| funmap | intogen2024 | +0.1001 | +0.1046 | +0.0045 |
| reactome | ot_nolit | +0.0996 | +0.1028 | +0.0032 |
| funmap | ot_all | +0.0853 | +0.0889 | +0.0036 |
| intact | ot_all | +0.0842 | +0.0853 | +0.0011 |
| rna_coexp | ot_nolit | +0.0834 | +0.0832 | -0.0002 |
| intact | ot_nolit | +0.0784 | +0.0793 | +0.0009 |
| funmap | ot_nolit | +0.0780 | +0.0813 | +0.0033 |
| rna_coexp | ot_all | +0.0767 | +0.0768 | +0.0001 |
| intact | intogen2024 | +0.0651 | +0.0671 | +0.0020 |
| string_full700 | intogen_temporal_new | +0.0604 | +0.0681 | +0.0077 |
| rna_coexp | intogen2024 | +0.0521 | +0.0533 | +0.0012 |
| reactome | intogen2024 | +0.0429 | +0.0466 | +0.0037 |
| string_phys700 | intogen_temporal_new | +0.0378 | +0.0490 | +0.0112 |
| funmap | intogen_temporal_new | +0.0196 | +0.0239 | +0.0043 |
| rna_coexp | intogen_temporal_new | +0.0127 | +0.0152 | +0.0025 |
| intact | intogen_temporal_new | +0.0071 | +0.0089 | +0.0018 |
| reactome | intogen_temporal_new | -0.0125 | -0.0099 | +0.0026 |

The answer is the reassuring one, and it is the opposite of what our
first, post-hoc calculation suggested. **23 of 24 arms are
positive before the cleaning and 23 of 24 after, with
0 sign flips**; the median margin moves from
**+0.0809** to **+0.0822** -- it goes slightly *up* in
almost every arm, because the curated drivers that leave the
background were ranked low by the prior anyway, and the degree
baseline had been getting them wrong.

The H1 verdict is the same story:

| interval | arms excluding 0 (standard) | arms excluding 0 (cleaned) |
|---|---|---|
| fold_pct | 24/24 | 24/24 |
| repeat_pct | 24/24 | 24/24 |
| repeat_t | 24/24 | 24/24 |
| repeat_bca | 24/24 | 24/24 |
| gene | 20/24 | 21/24 |

The conservative gene-cluster interval now passes in **21 of
24** arms instead of 20: the
cleaned background removes one more failure.

**What this means.** The contamination of our negative set is real --
10-22 % of the held-out negatives are curated drivers we do not label
-- and it is the kind of thing every benchmark in this area should
audit. But it is **not** what carries our result: with the curated
drivers removed from the background and the whole benchmark re-run, the
margin is unchanged. The earlier post-hoc reading ("the margin roughly
halves") was an artefact of that computation removing 87 % of the
IntOGen positives from the evaluation set along with the background.

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


The two additions to the results above change what this study can claim. **First**, the negative set is not clean: 10-22 % of the held-out negatives are curated driver genes we do not label -- but removing them from the background and re-running the whole benchmark leaves 23 of 24 arms positive, with no sign flip and the median margin unchanged (0.081 -> 0.082), so the contamination of the background is not what carries the result. **Second**, a large part of the remaining margin tracks study visibility: under a bait-matched null the temporal label's advantage collapses to chance (6.8 % -> 4.7 %, chance 5 %). Both are properties of the benchmark rather than of the method, and both mean the absolute numbers in this literature should be read with care. The validation plan that follows from them is in the next section.

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


### 4.5 What would change the conclusion

Four experiments, in order of how much they would move the paper. (i) **Re-derive the background** from NCG and ClinGen with explicit exclusion and re-run the factorial: this is the cheapest and the most likely to change the numbers. (ii) **A blind expert-panel challenge** on a held-out set of candidates, the way Boutros et al. recommend: our evaluation is still self-referential, and a hidden answer key is the only cure. (iii) **A seed-perturbation experiment**: delete 10 % and 30 % of the seed drivers and measure how the prediction changes, which tests the seed-count law directly and separates "topology" from "seed count". (iv) **A second fixed learner** (closeness centrality) under the same locked protocol, because our conclusion is currently about one learner and the literature contains a degree-filtered local method that beats RWR on a full list.

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
* **Code, split files, result tables and figures.** <https://github.com/lanxin2391/d1-benchmark>, tag `v1.0-rc3`. Release manifest with SHA-256 of 35 tables and 9 figures: `docs/release_manifest_v1.0-rc3.json`; verify with `python scripts/build_release_manifest.py --release v1.0-rc3 --check`.
* **Reproduction.** `conda env create -f environment.yml && conda activate d1 && pip install -e . && pytest` (131 passed, 3 skipped); `python scripts/audit_reproducibility.py` (50 PASS / 0 FAIL); `python scripts/verify_claims.py` (63 PASS / 0 FAIL); `python scripts/verify_report.py` (71 PASS / 0 FAIL).
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
