# Stability atlas of the network prior (work plan S1)

**Date:** 2026-11-05
**Script:** `scripts/run_stability_atlas.py`
**Inputs:** `results/scores_v40/*.npz` (24 arms, per-gene RWR scores for
3 alphas x 40 repeats x 5 folds), `results/runs_v40/a4_combined.tsv`
**Outputs:** `results/tables/stability_per_gene.tsv` (267,696 rows),
`results/tables/stability_by_arm.tsv` (24 rows),
`figures/figure5_stability_atlas.png`
**Status:** exploratory (work-plan Phase S is registered in the
exploratory family, see plan R19)

---

## 1. Why this analysis

The teacher (2026-11-04) pointed at a proteomics paper — Klingelhuber et
al., *A spatiotemporal proteomic map of human adipogenesis*, Nature
Metabolism 6:861–879 (2024), DOI 10.1038/s42255-024-01025-8 — and asked
for the **quantitative source of our network's stability**. That paper's
move is to treat the stability of the measurement layer as a first-class
quantity: it reports the per-protein coefficient of variation, the
replicate correlation, and the distribution of both, and then asks which
biological features explain the stable core.

Our analogue: the "measurement" is the RWR prior's score for a gene.
For every arm we have 40 re-partitions of the same universe, so we can
measure how much a gene's score and its place in the ranking move, and
then ask what predicts that movement.

## 2. Definitions

For each (network, label) arm at alpha = 0.5:

| quantity | definition |
|---|---|
| **evaluated score** | for each gene and each of the 40 repeats, the score of the **one fold in which that gene was in the test set** — the score that actually entered the AUROC |
| **all-fold score** | the mean over the 5 folds of a repeat (the repeat-level estimate used by the primary CI); reported as a sensitivity variant |
| **score CV** | SD / mean of the evaluated score across the 40 repeats — the direct analogue of the paper's per-protein CoV |
| **percentile rank** | the gene's rank among all genes of the universe within one repeat, scaled to [0, 1] |
| **rank SD** | SD of the percentile rank across the 40 repeats — "does the gene keep its place in the ranking?" |
| **local score density** | the number of genes whose score lies within +-10 % of this gene's score (measured on the repeat-mean scores) — the crowding of the gene's neighbourhood |
| **top-k hits** | in how many of the 40 repeats the gene is inside the top 1 % / 5 % of the universe |

Arm-level stability (the endpoint's own movement) is added: the SD of
the 40 repeat means of delta-AUROC, and the intra-repeat ICC of the
fold-level deltas.

## 3. Results

### 3.1 The endpoint is stable; the per-gene scores are not

| quantity | median per arm | range |
|---|---|---|
| repeat-level SD of delta-AUROC | **0.0047** | 0.0017 – 0.0129 |
| gene-level score CV (background) | **0.185** | 0.100 – 0.263 |
| gene-level rank SD (background) | **0.0267** | 0.0174 – 0.0412 |

The arm-level endpoint moves by about +-0.005 when the universe is
re-partitioned, which is 2–30x smaller than the effects themselves
(0.007 – 0.171): **the effect is much larger than the re-partition
noise**. At the level of a single gene, however, the prior is far from
deterministic: the relative error of a gene's score is ~18 %, and a
median gene moves ~2.7 percentile points, i.e. **~270 places out of
~10,000**, between two re-partitions of the same universe.

### 3.2 The source of the gene-level instability is crowding, not the network

| correlation (background genes) | median over 24 arms | sign |
|---|---|---|
| rank SD vs **local score density** | **+0.468** | positive in **24/24** arms |
| rank SD vs score magnitude | +0.393 | positive in 24/24 |
| rank SD vs degree | +0.074 | ~0 |
| score CV vs degree | **-0.137** | negative in 22/24 |
| score CV vs local score density | +0.083 | positive in 19/24 |

**A gene is unstable in the ranking exactly when it has many
near-equal-scored neighbours.** The density of the score distribution
explains the instability in every arm, while degree explains almost
nothing of it. Higher-degree genes are, if anything, *more* stable in
relative terms (CV vs degree is negative in 22 of 24 arms) — the
opposite of the intuition that hubs are noisy because they have many
ties.

The rank-SD profile over score deciles is an inverted U (pooled over the
24 arms): 0.0056 at the bottom decile, 0.0382 at decile 7, 0.0162 at the
top. Two effects contribute: (i) crowding (median local density rises
from 119 genes in the bottom decile to 535 in decile 6 and falls back to
193 in the top decile), and (ii) a boundary effect — a gene pinned near
the bottom or the top of the ranking cannot move much in percentile
terms, so the two ends of the inverted U are partly artefactual. The CV profile
is much flatter (0.134 -> 0.179 -> 0.155), which is why **score CV is
the primary metric** and rank SD the secondary one.

### 3.3 Driver genes are not intrinsically more stable

| comparison | result |
|---|---|
| drivers have a lower rank SD than background | 19/24 arms (median advantage 0.0033) |
| drivers have a lower score CV than background | 16/24 arms |
| within score deciles 1–7 | drivers are **slightly less** stable |
| within score deciles 8–10 | drivers are more stable (0.0083 vs 0.0162 in the top decile) |

The apparent driver advantage is a **position effect**: drivers sit
higher in the score distribution, and the top of the distribution is
stable. Stratified by score decile the advantage largely disappears.
So the honest statement is: *the prior ranks driver genes about as
unstablely as it ranks background genes; what is stable is the
endpoint, not the individual gene.*

Related observation worth reporting: the **top 1 % of the ranking is
not where drivers are**. Mean top-1 % hits over the 40 repeats are 0.073
for drivers versus 0.413 for background genes (0.18x), and only 0.6 % of
drivers are ever in the top 1 %. The RWR prior improves the *global*
ranking (which is what AUROC measures); it does not concentrate drivers
at the extreme top. This is the same warning as protocol §9 ("candidate
genes look like hubs, not drivers").

### 3.4 Which networks and labels are stable

| network | median background score CV | median background rank SD | median arm SD of delta-AUROC |
|---|---|---|---|
| intact | 0.103 | 0.0196 | 0.0026 |
| string_full700 | 0.138 | 0.0209 | 0.0023 |
| rna_coexp | 0.192 | 0.0251 | 0.0061 |
| funmap | 0.183 | 0.0269 | 0.0057 |
| string_phys700 | 0.183 | 0.0278 | 0.0040 |
| reactome | 0.206 | 0.0323 | 0.0050 |

| label | median background score CV | median arm SD of delta-AUROC |
|---|---|---|
| intogen2024 | 0.184 | 0.0044 |
| ot_all | 0.180 | 0.0046 |
| ot_nolit | 0.180 | 0.0047 |
| **intogen_temporal_new** | **0.232** | **0.0099** |

Two clear patterns:

* **intact and string_full700 are the most stable networks** at both the
  gene level (lowest CV) and the arm level (smallest repeat-level SD);
  **reactome is the least stable** (CV 0.206, rank SD 0.0323).
* **The temporal label is the least stable label in every network**
  (arm-level SD 0.0065–0.0129 versus 0.0017–0.0063 for the other
  three). The 152 newly discovered drivers are a small positive set
  (5–15 % of each universe), so both the effect and its noise are
  larger — the same heteroscedasticity that motivated the pre-registered
  gene-level interval (work plan B3).

The intra-repeat ICC is negative in **24/24** arms (median -0.178,
range -0.230 to -0.099), which is the independent confirmation of the
decision to use the repeat as the resampling unit (work plan A1).

### 3.5 The Reactome anomaly, in stability terms

`reactome x intogen_temporal_new` is simultaneously the **least stable**
arm (background rank SD 0.0412, arm-level SD 0.0099, drivers' rank SD
0.0479 — *less* stable than its own background) and the only negative
arm. The three robustness controls (D-09/D-11/D-12) already showed that
the negative delta-AUROC is not a coverage artefact; the stability atlas
adds that it is not a precision artefact either: the arm is unstable
for the same reason it is negative — the temporal drivers are exactly
the genes the network cannot rank.

## 4. Caveats

1. The evaluated score (one fold per repeat) and the all-fold mean are
   two defensible definitions of "the score of that repeat"; they agree
   only moderately (per-arm Spearman of the two stability measures:
   median 0.28). The all-fold variant is reported as a sensitivity
   check and should be added to the released table.
2. The boundary effect in rank SD (section 3.2) means the U-shape must
   not be read as "extreme genes are stable"; score CV is the metric to
   quote.
3. The local-density explanation is a **description**, not a causal
   test: density is itself a function of the score distribution. The
   formal version (regressing the stability of a gene on topological
   features) is S2.
4. Everything here is **exploratory** (work plan Phase S, R19): it does
   not touch the pre-registered endpoint or the H1 verdict.

## 5. What comes next

* **S2** — regress the per-gene stability on topological features
  (degree, k-core, clustering, betweenness, module membership, path
  count to the seed set) with network and label random effects, and
  check whether the crowding explanation survives as a property of the
  *neighbourhood* rather than of the score value.
* **S3** — relate the arm-level stability to network-level descriptors
  (needs the two STRING sensitivity arms for n = 8 networks).
* **S4** — ML attribution with stability selection (the PCP + SVM
  analogue), reusing the complex-aware splitter from Phase D.
* **S5** — repeat the stability analysis on the ClinGen label and the
  transfer arms, as independent systems.
