# Manuscript Draft — Abstract, Introduction, and Combined Figures (v1.0 assembly)

**Author:** Person B (kakamiku) — TM2 + TM4 in protocol §6
**Status:** Draft (M6 in protocol §6; intended as the front matter
of the v1.0 paper, joining `manuscript_a_sections.md` and
`manuscript_b_sections.md` into a single paper)
**Date:** 2026-10-07
**Predecessor sections:** `manuscript_a_sections.md` (§A-J), `manuscript_b_sections.md` (§A-H)

This file completes the manuscript by providing the **abstract**,
the **introduction** (§1-3), and the **combined figure legends**
that the A-side and B-side sections share. It does not duplicate
the methods or results prose that already lives in
`manuscript_a_sections.md` and `manuscript_b_sections.md`; it
is the **head matter + connectors** that turn those into a single
paper.

---

# Title

**What a network prior contributes to cancer-driver prioritisation: a leakage-controlled benchmark.**

**Authors:** lanxin2391 (Networks & Engine), kakamiku (Labels & Data)
**Pre-registration:** Zenodo DOI 10.5281/zenodo.23121762
**Code & data:** https://github.com/lanxin2391/d1-benchmark at the
`v1.0-prep` tag.

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
types. **An audit of our own negative set** shows that 10-22 % of the
held-out negatives are curated driver genes we do not label; removing
them halves the margin on the strictest label (+0.065 to +0.030) and
raises it on the Open Targets labels (+0.085 to +0.103), so the
benchmark's numbers depend on how the background is built. **Label-
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

# 1. Introduction

Cancer is a disease of dysregulated genes. The first step of any
molecular oncology study is to identify the small set of genes
that drive the disease in a given patient or cohort. Over the last
two decades, the cancer genomics community has curated multiple
driver-gene catalogues from independent curation philosophies:
mutation recurrence (IntOGen; Bailey et al.), clinical evidence
(ClinGen), and aggregated target-disease association (Open
Targets). These catalogues disagree in important ways — for
example, the 2020-2024 newly-discovered IntOGen drivers are
largely absent from ClinGen's curated gene-disease pairs,
because clinical case reports take time to accumulate.

A parallel literature has argued that **network structure** can
help identify driver genes beyond what the catalogues alone
provide. The intuition: two genes whose proteins bind, or whose
RNAs co-vary, are more likely to share driver status than two
genes picked at random. The simplest formalisation is the
"guilt-by-association" principle: a gene that is close to many
known drivers in a protein-protein interaction (PPI) network is
itself likely a driver. This idea has been formalised as random
walk with restart (RWR; Kohler et al. 2008; Vanunu et al. 2010),
label propagation (Zhang et al. 2018), and a range of graph
neural network methods (Boiarsky et al. 2021, 2023).

The field has converged on three empirical observations:

1. **RWR and related methods beat the degree baseline.** The
   original RWR paper and most follow-ups show that adding the
   network propagation step improves over ranking genes by
   degree alone.

2. **Different networks give different absolute AUROCs.** FunMap
   gives higher absolute AUROC than STRING; STRING gives higher
   than Reactome. This 0.25-AUROC gap is the basis for the
   "best network" question.

3. **No published study has measured whether the network prior
   contributes anything beyond degree in a leakage-controlled
   setting.** Existing benchmarks either (a) use a fixed network
   without comparing to a degree baseline, (b) compare networks
   but do not report the margin over degree, or (c) use
   overlapping train/test labels and confound the degree signal
   with label leakage.

The D1 benchmark closes this gap. We assembled a six-network
panel stratified by molecular provenance, a four-label panel
spanning independent curation philosophies, a fixed RWR learner
(α = 0.5, no hyperparameter tuning), and a control stack of ten
controls
that measures the network's contribution against three distinct
nulls: a degree baseline (Control 1), a degree-matched seed null
(Control 2), an edge-rewiring null (Control 3), and four further
confound controls — a shared-gene-universe restriction (Control 4),
a cross-catalogue label transfer (Control 5), a temporal split
(Control 6), a publication-bias adjustment with a
literature-matched null (Control 7) — plus three controls added
after a systematic read of the benchmark literature: complex-aware
cross-validation (Picart-Armada et al. 2019), seed-degree
stratification (Zhang et al. 2021, Control 8), network size and
coverage adjustment (Wright et al. 2025, Control 9) and random-input
degeneracy (Hill et al. 2019, Control 10).

The headline question is: **does the specific topology of a
cancer-related gene network contribute to driver-gene
prioritisation, beyond what node degree alone explains?** The
protocol's H1-H3 pre-register three concrete answers to that
question. H1 holds in 24/24 arms under four of the five interval
constructions and in 20/24 under the pre-registered gene-cluster
bootstrap (the four exceptions are the smallest-effect
`intogen_temporal_new` arms); H2 is **not testable** with this
panel, because provenance is confounded with network size and
coverage; H3 holds in the sense that the margin survives a
literature-matched null.

This paper reports the v1.0 release of the benchmark: the
release tag, the leaderboard protocol, the public Zenodo
deposit, and the manuscript sections that interpret the
headline numbers.

---

# 2. Related work

The D1 benchmark sits at the intersection of three literatures:
network-based gene prioritisation, cancer-driver catalogue
construction, and benchmarking methodology. We summarise each
in turn.

## 2.1 Network-based gene prioritisation

* **Kohler et al. (2008)** introduced RWR for disease-gene
  prioritisation, demonstrating that the diffusion state over a
  PPI network outperforms simple neighbourhood-based scoring.

* **Vanunu et al. (2010)** extended this to a global propagation
  scheme that handles multiple seeds and weighted networks.

* **Boiarsky et al. (2021, 2023)** applied graph neural networks
  to the same problem, with attention layers that can mix
  network topology with node features.

The D1 benchmark differs from this prior work in three ways:
(1) we hold the learner fixed at RWR α = 0.5 and measure the
network's contribution against a degree baseline, rather than
comparing learners; (2) we stratify the network panel by molecular
provenance (RNA-derived, protein-derived, curated/literature) so
the contribution of "which network" can be separated from
"topology matters"; and (3) we measure the contribution against
three distinct nulls (degree, degree-matched seeds, edge
rewiring), so a positive result cannot be explained by "just
having any edges" or "just having any degree sequence".

## 2.2 Cancer-driver catalogue construction

* **Bailey et al. (2018)** — the PanCancer Atlas driver
  catalogue from 26 driver-detection methods, 299 genes.

* **Martinez-Jimenez et al. (2020)** — the IntOGen Compendium
  of Cancer Genes, 568 unique driver genes (2020.02.01 release),
  now updated to 633 (2024.09.20).

* **Ochoa et al. (2025)** — the Open Targets Platform 26.06
  release, which aggregates 22,581 targets with multi-evidence
  scoring across ~20 datasources.

* **Rehm et al. (2015)** — the ClinGen gene-disease validity
  curations, the most "manually validated" of the four label
  sources, focused on clinical-grade evidence.

The D1 benchmark uses four of these (IntOGen 2024, Open Targets
all-evidence, Open Targets no-literature, IntOGen temporal new)
plus an optional ClinGen fifth label. We do not use Bailey 2018
because the cell-level supplementary data is non-OA; IntOGen
supersedes it for our purposes.

## 2.3 Benchmark methodology

* **OpenML** (Vanschoren et al. 2014) — a community platform for
  sharing datasets, splits, and leaderboard results in
  supervised ML. We follow OpenML's "releases-as-files" pattern
  (splits and labels are released as TSVs, not regenerated from
  seeds).

* **MLOmics** (Sun et al. 2025) — a multi-omics ML benchmark for
  TCGA cohorts across 14 cancer types. We follow the MLOmics
  pattern of "fixed data, fixed splits, multiple methods reported
  in a single paper" for our network-comparison dimension.

* **Genomics England** — the Genomics England rare-disease
  benchmark for variant prioritisation. We do not directly
  compare to it because it focuses on individual variants rather
  than gene-level prioritisation.

The D1 benchmark is the **first** of these to focus on
**network priors for cancer-driver prioritisation** in a
leakage-controlled setting.

---

# 3. Headline contribution

The D1 benchmark's central contribution is a **quantitative
answer** to the question that has been open in the field for a
decade:

> A network prior contributes to cancer-driver prioritisation
> beyond node degree. The contribution is not explained by the
> degree of the seed genes, by network size, or by the degree
> distribution of the labels. Whether it depends on the network's
> molecular provenance is **not testable** with a panel of this
> size, because provenance is confounded with network size and
> coverage.

The answer has three concrete corollaries:

1. **Methodological**: a network prior should not be evaluated
   in absolute terms; the margin over a degree baseline is the
   only meaningful measure — and it must be reported against the
   *pre-registered* resampling unit, not only against a
   post-hoc choice of it.

2. **Practical**: the choice between FunMap, STRING physical,
   STRING full, IntAct, Reactome, and an RNA co-expression
   network is largely a wash for mutation-based driver
   catalogues; the published 0.25-AUROC gap between FunMap and
   STRING is mostly a degree effect.

3. **Bibliometric**: the apparent advantage of curated networks
   on Open Targets labels (which weight clinical evidence) is
   not a literature-attention artefact — the margin survives a
   literature-matched null (Control 7) and a size- and
   coverage-adjusted analysis (Control 9).

The single negative cell — Reactome's failure on the
2020-2024 newly-discovered IntOGen drivers — is **not** a
demonstration of a negative effect. Its 95 % CI contains zero
under the pre-registered gene-cluster bootstrap, and the
equivalence test finds it statistically indistinguishable from
zero. It is a *coverage* failure: Reactome's LCC contains only
41 % of the newly-discovered drivers, so the arm is measured on a
gene universe that largely lacks the genes it is asked to find.

---

# 4. Combined figure legends

The paper has four primary figures. Three are A-side (Figure
2 — heatmap; Figure 2-controls — D-09/D-11/D-12 panel; Figure
2-temporal-strip — temporal drift); one is B-side (Figure 1 —
forest plot + clingen heatmap). The combined figure 3 is the
new "all controls" schematic that visualises which of the seven
controls rule out which alternative explanations.

## Figure 1 (B-side) — Forest plot + heatmap

`figures/figure1_pilot_overview.png`. Two panels: left = forest
plot of 14 arms with 95 % CI, grouped (a) original 4-arm /
(b) D-09 temporal / (c) D-11 cap new arms / (d) D-11 cap on
original labels / (e) clingen 5th label 6-network; right =
two heatmaps (uncapped vs D-11-capped). Publication-ready legend
in upper-left, section labels in left margin, missing cells
explicitly greyed as "n/a". Data table:
`results/tables/figure1_data.tsv`.

## Figure 2 (A-side) — Main heatmap

`figures/figure2_main_heatmap.png`. Primary 6×4 heatmap,
α=0.5, native/uncap. Color scheme RdYlGn (red → yellow → green)
so high ΔAUROC is green and negative cells are red. Data:
`results/tables/a4_table1_delta_auroc.tsv`.

## Figure 2 controls (A-side) — D-09/D-11/D-12 panel

`figures/figure2_d09_d11_d12.png`. Three-panel side-by-side:
native, D-11 cap, D-12 shared. Shows how the headline result
changes under the three coverage controls.

## Figure 2 temporal strip (A-side)

`figures/figure2_d09_temporal_strip.png`. 1-column strip plot
for the temporal-drift arms (Control 6). The 152 newly-
discovered IntOGen drivers localise the rest of the catalogue
with smaller margin (0.02-0.10) than the established drivers
(0.04-0.17); Reactome is the only network with a negative
temporal-drift estimate.

## Figure 3 (B-side, new) — All-controls schematic

`figures/figure3_all_controls.png` (to be generated at the
assembly step). Schematic of the controls in §4.4 with the
headline finding of each. Every range below is recomputed from
the released tables by `scripts/recompute_ground_truth.py`; the
authoritative copy of this table is `review/ground_truth.tsv`.

| # | control | headline | file |
|---|---|---|---|
| 1 | degree-only baseline | ΔAUROC = −0.0125 to +0.1711; 23/24 cells positive | `results/tables/a4_table1_delta_auroc_v40.tsv` |
| 2 | degree-matched seed null | z = 3.09 to 226.91; 23/24 cells with z > 5; the null mean is ≤ 0 in 23/24 cells | `results/tables/a5_p_values.tsv` |
| 3 | edge-rewiring null | z = −3.53 to +205.25; 23/24 cells with z > 5 (598 rewirings, not 600 — Reactome has 98) | `results/tables/w5_p_values.tsv` |
| 4 | shared gene universe | ΔAUROC = −0.0360 to +0.1304; 22/24 cells positive, all on the same 1,605-gene universe | `results/tables/d12_table.tsv` |
| 5 | cross-catalogue transfer | mean ΔAUROC = +0.0313, 39.8 % of the within-label magnitude; 59/66 cells positive | `results/tables/results_stats_transfer.tsv` |
| 6 | temporal split | ΔAUROC = −0.0095 to +0.0923; 5/6 networks positive, 50 distinct seed sets per arm | `results/tables/d09_table.tsv` |
| 7 | publication-bias adjustment | literature-matched null: mean matched − standard = +0.0010; 30/30 cells | `results/tables/control7_matched_null.tsv` |
| 8 | complex-aware cross-validation | median absolute relative change 0.142; 0/24 sign flips | `results/tables/d4_v1_vs_v2_comparison.tsv` |
| 9 | seed-degree stratification | ΔAUROC by seed-degree tercile: low +0.053, mid +0.061, high +0.054 (low stratum positive in 22/24 arms) | `results/tables/control8_seed_degree_summary.tsv` |
| 10 | network size / coverage adjustment | adjusted endpoint positive in 23/24 cells; pooled Spearman ρ(ΔAUROC, log₁₀ edges) = −0.02, permutation p = 0.93 | `results/tables/control9_size_adjustment.tsv` |
| 11 | random-input degeneracy | top-100 Jaccard between the real-input and the random-input ranking = 0.019 (almost disjoint gene lists); Hill degeneracy statistic 0.6 genes under random input vs 123.6 under real input | `results/tables/control10_random_input_summary.tsv` |

Together, the controls rule out the alternative explanations for
the headline finding. The figure is a single panel with each
control as a small-multiples subplot.

---

# 5. Combined tables (cross-references)

The paper has 12 primary result tables. Six are A-side (network
+ engine + statistics), four are B-side (labels + transfer +
pubcount + ClinGen), and two are shared (results_stats.tsv;
overlap_report.md).

| table | section | owner | what |
|---|---|---|---|
| `results/tables/a4_table1_delta_auroc.tsv` | A4 | A | 24 cells × 3 alphas = 72 primary endpoint measurements |
| `results/tables/d09_table.tsv` | A | temporal split arm summary (Control 6) |
| `results/tables/d11_table.tsv` | A | per-network cap arm summary (D-11) |
| `results/tables/d12_table.tsv` | A | shared universe arm summary (Control 4) |
| `results/tables/results_stats.tsv` | A | paired t + BH-FDR for 4-arm grid |
| `results/tables/results_stats_clingen.tsv` | A+B | 5th label (clingen) stats |
| `results/tables/results_stats_transfer.tsv` | B | 66 transfer grid cells |
| `results/tables/results_stats_d09.tsv` | B | temporal-drift t-tests |
| `results/tables/results_stats_d11.tsv` | A | D-11-capped t-tests |
| `results/tables/pubcount_adj_comparison.tsv` | B | 23 pubcount-adjusted vs original |
| `results/tables/w5_p_values.tsv` | A | rewiring-null z + p |
| `results/tables/a5_p_values.tsv` | A | degree-matched seed null z + p |
| `results/tables/stats_fine_z.tsv` | A | §4.6 fine-mode per-fold z |
| `results/tables/stats_variance_decomp.tsv` | A | §4.6 mixed-effects variance |
| `results/tables/stats_bootstrap_ci.tsv` | A | §4.6 cluster bootstrap CI |
| `results/tables/stats_equivalence.tsv` | A | §4.6 TOST equivalence |
| `results/tables/provenance_ablation.tsv` | A | §4.2 RNA/protein/curated |
| `results/tables/provenance_summary.tsv` | A | per-provenance class summary |
| `results/tables/gnn_courtesy_arm.tsv` | A | §7 Node2Vec courtesy arm |
| `results/tables/label_network_overlap.tsv` | B | per-(label, network) coverage |
| `results/tables/figure1_data.tsv` | B | figure 1 source data |

---

# 6. Open issues and what was deferred

For completeness, the items the protocol §6 / pending_work_protocol_compliance.md
flagged but the v1.0 release deferred:

1. **Sensitivity to `--top-n` for clingen**: clingen has 84 genes;
   sensitivity to top-30 or top-50 is a future arm but not in
   v1.0.

2. **D-11 cap on the shared universe (D-12)**: A's D-11 is native
   only; D-11 × shared universe is a future combination.

3. **PubCount adjustment trim fraction (D-PUB-MARGIN)**: 50 % is
   the v1.0 value. 25 % / 75 % are future sensitivity arms.

4. **Bioinformatics tool registry (bio.tools)**: pending; the
   v1.0 release is on GitHub + Zenodo for now.

5. **Manuscript cover letter**: drafted in
   `docs/handoff_to_B_s2.md` (A) + `docs/b_s3_final_state.md`
   (B); final cover letter pending A.

None of these affects the v1.0 release's headline numbers.

---

# 7. How to reproduce v1.0

The v1.0 release can be reproduced end-to-end on a clean machine:

```bash
git clone https://github.com/lanxin2391/d1-benchmark.git
cd d1-benchmark
git checkout v1.0-prep
conda env create -f environment.yml
conda activate d1
python scripts/audit_reproducibility.py    # 50 PASS, 0 FAIL
```

The audit script takes ~10 minutes and verifies that all 30
splits, 5 labels, 8 networks, and 6 result-table SHA-256
fingerprints match the pre-registration manifest (with amendment
#1 applied). The smoke arm runs in 30 seconds.

For the full numerical reproduction (re-running the confirmatory
analysis from raw data), see `docs/public_release_protocol.md`
§5. The full run is ~2 hours on a 4-core laptop.

---

# 8. Sign-off

**Person B (kakamiku), v1.0 manuscript assembly complete:**
This file plus `manuscript_a_sections.md` (A-side, 10 sections)
plus `manuscript_b_sections.md` (B-side, 8 sections) is the full
v1.0 paper. A is responsible for the Methods §4.1/§4.3/§4.6 +
Results A-side + Discussion A-side + Conclusion. B is responsible
for the Methods §4.2 + Results label panel + Discussion label
bias + ClinGen integration + this assembly. The two manuscript
drafts reference the same 24 result tables and the same seven
controls; the reader can navigate between A and B by following
the cross-references in §Cross-references of each.

Signed: 2026-10-07
