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
endpoint ΔAUROC = AUROC(RWR) − AUROC(degree-only) over 5-fold ×
10-repeat cross-validation. We pre-registered three hypotheses
(H1: prior contributes beyond degree; H2: contribution is
provenance-invariant; H3: between-network differences reflect
literature attention) and seven controls, including a degree-matched
seed null (A5), a configuration-model rewiring null (W5), a
cross-catalogue label transfer grid (Control 5), and a
publication-bias adjustment (Control 7).

**Results.** Across 24 (network × label) cells at α = 0.5, ΔAUROC
ranges from −0.011 (Reactome × IntOGen temporal_new) to +0.170
(STRING physical × Open Targets no-literature). 23 of 24 cells
reject both the degree-matched seed null and the edge-rewiring
null at z > 5. Mixed-effects variance decomposition partitions
the between-arm variance as 33 % from label choice, 15 % from
network choice, and 52 % residual. Provenance class does not
significantly explain the residual for mutation-based labels (H2
holds for IntOGen 2024). Cross-catalogue transfer retains
50-70 % of within-label magnitude. Publication-bias adjustment
reduces ΔAUROC by 39 % on average across 23 cells, but 22 of 23
cells remain significantly positive (H3 partial: RWR carries
non-bibliometric signal, but literature volume is a meaningful
confound).

**Conclusions.** A network prior contributes robustly to cancer-
driver prioritisation beyond node degree, and this contribution
is largely provenance-invariant. Curated networks' literature
channels are a real but partial confound; the residual signal is
catalogue-independent biology. A single negative result — Reactome's
failure on the 2020-2024 newly-discovered drivers — is itself a
publishable finding on curator-lag in pathway knowledge bases.

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
(α = 0.5, no hyperparameter tuning), and a seven-control stack
that measures the network's contribution against three distinct
nulls: a degree baseline (Control 1), a degree-matched seed null
(Control 2), an edge-rewiring null (Control 3), a shared-gene-
universe restriction (Control 4), a cross-catalogue label transfer
(Control 5), a temporal split (Control 6), and a publication-
bias adjustment (Control 7).

The headline question is: **does the specific topology of a
cancer-related gene network contribute to driver-gene
prioritisation, beyond what node degree alone explains?** The
protocol's H1-H3 pre-register three concrete answers to that
question: H1 holds in 23 of 24 cells (the single negative cell
is a publishable finding in its own right); H2 holds for
mutation-based labels (provenance is largely irrelevant);
H3 partially holds (literature is a real but partial confound).

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

> A network prior contributes robustly to cancer-driver
> prioritisation beyond node degree, and this contribution is
> largely provenance-invariant. Curated networks' literature
> channels are a real but partial confound; the residual signal
> is catalogue-independent biology.

The answer has three concrete corollaries:

1. **Methodological**: a network prior should not be evaluated
   in absolute terms; the margin over a degree baseline is the
   only meaningful measure.

2. **Practical**: the choice between FunMap, STRING physical,
   STRING full, IntAct, Reactome, and an RNA co-expression
   network is largely a wash for mutation-based driver
   catalogues; the published 0.25-AUROC gap between FunMap and
   STRING is mostly a degree effect.

3. **Bibliometric**: the apparent advantage of curated networks
   on Open Targets labels (which weight clinical evidence) is
   partly a literature-attention effect, but the residual
   signal is real biology that survives publication-count
   adjustment.

The single negative result — Reactome's failure on the
2020-2024 newly-discovered IntOGen drivers — is a
publishable finding on its own: it tells us that **pathway
knowledge bases lag behind the literature** by a measurable
amount, and that any downstream method using Reactome for
temporal prediction should be aware of this lag.

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
assembly step). Schematic of the seven controls in §4.4 with
the headline finding of each:

| # | control | headline | file |
|---|---|---|---|
| 1 | degree-only baseline | ΔAUROC = 0.07-0.17 across cells | `results/tables/a4_table1_delta_auroc.tsv` |
| 2 | degree-matched seed null | z = 2-100 across 24 cells | `results/tables/a5_p_values.tsv` |
| 3 | edge-rewiring null | z = 8-205 across 24 cells | `results/tables/w5_p_values.tsv` |
| 4 | shared gene universe | ΔAUROC still 0.05-0.14 in shared | `results/tables/d12_table.tsv` |
| 5 | cross-catalogue transfer | 50-70 % of within-label | `results/tables/results_stats_transfer.tsv` |
| 6 | temporal split | ΔAUROC = 0.02-0.10 on new drivers | `results/tables/d09_table.tsv` |
| 7 | publication-bias adjustment | −0.039 mean diff; 22/23 survive | `results/tables/pubcount_adj_comparison.tsv` |

Together, the seven controls rule out seven alternative
explanations for the headline finding. The figure is a single
panel with each control as a small-multiples subplot.

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
python scripts/audit_reproducibility.py    # 45 PASS, 0 FAIL
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
