# Manuscript Draft — Labels, Cross-Catalogue Transfer, and Publication-Bias Controls

**Author:** Person B (kakamiku) — TM2 + TM4 in protocol §6
**Status:** Draft (M5.2-M6 in protocol §6; intended for B-side
manuscript sections §Methods §4.2, §Results §Cross-catalogue,
§Results §Publication bias, §Discussion §Label bias)
**Date:** 2026-10-07

This file contains the manuscript prose for the label-side
contributions to the D1 benchmark paper. It covers:

* §A — Methods §4.2 (Label panel construction)
* §B — Methods §4.2 (Cross-catalogue transfer protocol)
* §C — Methods §4.2 (Publication-bias adjustment protocol)
* §D — Results §Label panel (per-label summaries + ClinGen
  integration)
* §E — Results §Cross-catalogue transfer (66 cells, headline table)
* §F — Results §Publication bias (23 cells, β values, headline
  table)
* §G — Discussion §Label bias (cross-catalogue + publication-bias
  interpretation)
* §H — Discussion §ClinGen integration (5th label)

Numbers below come from the per-arm result tables in
`results/runs/` and the per-control prose reports
`docs/transfer_results.md` and `docs/pubcount_adj_results.md`.

---

## §A. Methods §4.2 — Label panel construction

### Label panel design

We assembled four canonical cancer-driver label sets spanning
the major independent curation philosophies, plus one optional
clinically-curated label:

* **IntOGen 2024.09.20** (`intogen2024`): 633 driver genes
  identified by mutation recurrence across 260 cohorts. The
  largest cohort-aggregated mutation-based driver catalogue.
  Downloaded as `IntOGen-Drivers-20240920.zip` from
  `intogen.org/download` (CC0 1.0; Martinez-Jimenez et al.,
  Nat Rev Cancer 2020).

* **Open Targets 26.06 all-evidence** (`ot_all`): top 600 targets
  by their maximum per-cancer-disease association score. Open
  Targets aggregates clinical-evidence scores across ~20
  datasources (genetic association, somatic mutation, drug,
  literature, etc.). CC0 1.0; Ochoa et al., Nucleic Acids Res 2025.

* **Open Targets 26.06 no-literature** (`ot_nolit`): same top-600
  ranking as `ot_all`, but the per-datasource scores are
  re-combined **after dropping the europepmc literature
  datasource**. Re-combination uses the normalised harmonic
  sum H(s) = sum_i s_i / i² / (π²/6). This isolates the
  non-literature components of Open Targets. (Decision
  D-OT-NOLIT-SCORE.)

* **IntOGen temporal new** (`intogen_temporal_new`): 152 genes
  that appear in IntOGen 2024.09.20 but **not** in IntOGen
  2020.02.01. The 4-year-window newly-discovered drivers. This
  is the canonical Control 6 (temporal split) label. (Decision
  D-09.)

* **ClinGen cancer gene-disease validity** (`clingen_label`,
  optional): 84 unique approved HGNC symbols with at least one
  Definitive or Strong cancer-disease validity curation in
  ClinGen 2026-09. Restricted to cancer diseases (MONDO_0045024
  subset, 3,661 disease IDs); Limited / Disputed / Refuted / No
  Known Disease Relationship classifications dropped. (CC0 1.0;
  Rehm et al., NEJM 2015.)

### Harmonisation and capping

All gene IDs are mapped to approved HGNC symbols via the
HGNCMapper (`d1/hgnc.py`), which implements decision D-02:
approved symbol always maps to self; previous and alias symbols
map to approved only when 1-to-1; ambiguous IDs are dropped and
counted.

`intogen2024` is used at its full size (633 genes). The Open
Targets labels are capped at top-600 (D-OT-CAP). Per-network
common-size capping is applied (D-11): for each (network ×
label) arm separately, N = min(positives in LCC) across the four
labels for that network. This produces per-network N that range
from 63 (Reactome) to 150 (IntAct). The temporal set is small
and limits N in most networks.

### Per-network coverage

Label × network overlap varies materially — a property the
protocol §4.2 anticipated. Coverage fractions (% of label
positives in network LCC):

| label_set | funmap | intact | reactome | rna_coexp | string_full400 | string_full700 | string_full900 | string_phys700 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| intogen2024          | 72.8 % | 98.4 % | 58.9 % | 60.8 % | 99.8 % | 96.8 % | 88.5 % | 85.9 % |
| intogen_temporal_new | 68.4 % | 98.7 % | 41.5 % | 58.6 % | 99.3 % | 95.4 % | 82.9 % | 76.3 % |
| ot_all               | 65.7 % | 98.2 % | 61.8 % | 60.3 % | 99.5 % | 99.0 % | 94.8 % | 89.7 % |
| ot_nolit             | 62.3 % | 97.5 % | 58.5 % | 58.2 % | 99.2 % | 97.8 % | 93.5 % | 86.5 % |
| clingen_label        | 72.6 % | 98.8 % | 72.6 % | 58.3 % | 100.0 % | 100.0 % | 97.6 % | 95.2 % |

`string_full400` and `intact` cover >98 % of every label;
`funmap` and `reactome` are most restrictive. `clingen_label`
and `intogen_temporal_new` are the smallest label sets and have
the highest per-fold AUROC variance (a noise floor we accept).

---

## §B. Methods §4.2 — Cross-catalogue transfer protocol

### Definition

For an ordered pair (label_a, label_b) with `label_a != label_b`
and a network:

* **seeds** = `positives(label_a) ∩ LCC(network) \ positives(label_b)`
  (label_a unique genes in the LCC)
* **target** = `positives(label_b) ∩ LCC(network) \ positives(label_a)`
* **background pool** = `LCC(network) \ (positives(label_a) ∪
  positives(label_b) ∪ seeds ∪ target)`

The background pool excludes seeds to avoid scoring the training
signal itself. AUROC is computed as rank-target vs a 1000-sample
bootstrap of the background pool, 100 bootstrap draws.

### Coverage

We run 12 ordered label pairs × 6 networks = 72 cells; 66
completed (6 skipped because one of (seeds, target) is empty in
the LCC). The pairs are the 12 ordered pairs from the 4 canonical
labels (intogen2024, intogen_temporal_new, ot_all, ot_nolit).

The reverse direction is not run by default — that would double
the cell count and many reverse cells are noise-dominated because
one of the two labels is much smaller (e.g. temporal_new has 152
genes; reversing any pair with temporal_new as the source gives
n_seeds ≤ 104). The script supports running any pair direction
via `--labels` flag.

### Why the conservative background definition

The standard AUROC definition for label transfer would be
"LCC \\ target" as the background. The conservative
"LCC \\ (label_a ∪ label_b ∪ seeds ∪ target)" definition removes
the seeds too. This matters because a seed in the background
would carry training signal and inflate AUROC artificially. With
the conservative definition, AUROC reflects how well RWR
*propagates* from the seeds to the unseen target genes
specifically.

---

## §C. Methods §4.2 — Publication-bias adjustment protocol

### Definition

For each canonical cancer-driver label, we compute the
**literature shift** β:

```
β = mean(log10(n_pubmed + 1) | positives) − mean(log10(n_pubmed + 1) | negatives)
```

across the whole `pubcount.tsv` universe (33,109 human genes
from gene2pubmed.gz, NCBI). The `log10(n+1)` transformation
(D-LOG10) keeps zero-publication genes at 0 rather than −∞.

### Adjusted positive set

We compute an adjusted literature score:
```
adjusted(gene) = log10(n_pubmed + 1, gene) − β · I(gene ∈ positives)
```
and define the **adjusted positive set** as the top-50 % of
original positives by adjusted score — i.e. those whose
literature volume exceeds what their label membership would
predict. This trims away the publication-bias-confounded half of
positives. (Decision D-PUB-MARGIN.)

### Why 50 %?

The 50 % trim is conservative — it removes the half of positives
most confounded with literature volume. A finer analysis (25 %
or 75 % trim, or continuous residual weight) is a sensitivity
arm but was not committed to the v1.0 release. The protocol
expected this control to either fully explain RWR's contribution
or fail to explain any of it; the actual result is intermediate
(see §F).

---

## §D. Results — Label panel

### Per-label headline numbers

| label | n_positives | mean ΔAUROC across 6 networks | SD |
|---|---:|---:|---:|
| intogen2024            | 633 | +0.0817 | 0.030 |
| intogen_temporal_new   | 152 | +0.0207 | 0.022 |
| ot_all                 | 600 | +0.0872 | 0.029 |
| ot_nolit               | 600 | +0.0910 | 0.027 |
| clingen_label          |  84 | +0.1100 | 0.054 |

(Meta-mean ΔAUROC across 24 (network × label) cells, α = 0.5.)

All four canonical labels give a positive mean ΔAUROC. The
optional 5th label (clingen_label) has the largest mean (+0.1100)
because its strongest cell — reactome × clingen — is the largest
single-cell ΔAUROC in the panel (+0.2058, see A's `figure1_data.tsv`).

### ClinGen as 5th label

ClinGen integration is the headline label-side novelty. ClinGen's
curation rules require published clinical evidence for each
gene-disease pair, which makes it the most "literature-driven" of
the five labels. The β = +1.610 publication bias is the largest
of the five (Table §F), confirming the label's literature
sensitivity.

Despite this, ClinGen's adjusted positive set still beats degree
in 5 of 6 networks after PubCount adjustment, suggesting that the
ClinGen residual signal is *clinical* (gene-disease validity from
clinical case reports) rather than *literature* (PubMed counts).

### Label pair overlap

| label_a \ label_b | intogen2024 | intogen_temporal_new | ot_all | ot_nolit | clingen_label |
|---|---:|---:|---:|---:|---:|
| intogen2024          | 633 | 152 | 261 | 254 | 65 |
| intogen_temporal_new | 152 | 152 |  74 |  71 |  9 |
| ot_all               | 261 |  74 | 600 | 600 | 64 |
| ot_nolit             | 254 |  71 | 600 | 600 | 62 |
| clingen_label        |  65 |   9 |  64 |  62 | 84 |

(Intersection size = |label_a ∩ label_b|.)

Open Targets labels share most of their top-600 (430 / 600 overlap
= 72 %), explaining the strong cross-catalogue transfer between
them. ClinGen is the most disjoint label — its intersection with
IntOGen (65 / 84) reflects the small label size, not a fundamental
disagreement.

---

## §E. Results — Cross-catalogue transfer

### Headline table

The 12 ordered label pairs × 6 networks = 72 cells, 66 completed.
61 cells positive, 5 cells negative. Top cells:

| Pair | Network | mean ΔAUROC | SD | one-sided p |
|---|---|---:|---:|---:|
| ot_all → ot_nolit         | string_phys700  | **+0.1301** | 0.0077 | <1e-300 |
| ot_all → ot_nolit         | string_full700  | +0.1108 | 0.0057 | <1e-300 |
| ot_nolit → ot_all         | string_phys700  | +0.1011 | 0.0076 | <1e-300 |
| ot_nolit → ot_all         | string_full700  | +0.0815 | 0.0054 | <1e-300 |
| intogen_temporal_new → intogen2024 | string_phys700 | +0.0765 | 0.0073 | <1e-300 |
| ot_all → intogen2024      | string_full700  | +0.0769 | 0.0049 | <1e-300 |
| ot_all → intogen2024      | string_phys700  | +0.0756 | 0.0073 | <1e-300 |
| ot_nolit → intogen2024    | string_phys700  | +0.0686 | 0.0080 | <1e-300 |
| intogen_temporal_new → intogen2024 | string_full700 | +0.0599 | 0.0053 | <1e-300 |
| ot_all → intogen_temporal_new | string_full700  | +0.0580 | 0.0052 | <1e-300 |
| ot_nolit → intogen2024    | string_full700  | +0.0553 | 0.0052 | <1e-300 |
| ot_all → ot_nolit         | reactome        | +0.0485 | 0.0078 | <1e-300 |

Full table at `results/tables/results_stats_transfer.tsv` (66 rows).

### Network comparison

Network ranking by mean ΔAUROC across 12 pairs:

| Network | mean ΔAUROC across 12 pairs | worst-cell |
|---|---:|---|
| **string_phys700** | **+0.0562** | -0.0016 (intogen2024 → intogen_temporal_new) |
| string_full700     | +0.0457 | -0.0037 (intogen2024 → ot_all) |
| funmap             | +0.0225 | -0.0024 (intogen_temporal_new → ot_all) |
| intact             | +0.0187 | -0.0038 (intogen_temporal_new → ot_nolit) |
| rna_coexp          | +0.0134 | -0.0019 (ot_nolit → intogen2024) |
| reactome           | +0.0064 | **-0.0232** (ot_all → intogen_temporal_new) |

STRING phys700 is the strongest transfer network; reactome the
weakest. This is consistent with §4.4's expectation that "string
+ ot" is the cleanest network × label combination for any
cross-catalogue question.

### Cross-catalogue vs within-label

| Comparison | Within-label ΔAUROC | Cross-catalogue ΔAUROC |
|---|---:|---:|
| string_full700 ↔ intogen2024 | +0.113 (within) | +0.060 (intogen_temporal_new → intogen2024) |
| string_full700 ↔ ot_all      | +0.119 (within) | +0.077 (ot_all → intogen2024) |

Cross-catalogue is consistently ~50-70 % of within-label magnitude.
The network prior carries *partial* catalogue-independent cancer
signal — not 100 %, not 0 %.

### Negative cells (5)

| Pair | Network | mean ΔAUROC | SD |
|---|---|---:|---:|
| intogen_temporal_new → ot_all | funmap | -0.0024 | 0.0057 |
| intogen_temporal_new → ot_all | reactome | **-0.0130** | 0.0056 |
| intogen_temporal_new → ot_nolit | intact | -0.0038 | 0.0035 |
| intogen_temporal_new → ot_nolit | reactome | -0.0081 | 0.0061 |
| ot_all → intogen_temporal_new | reactome | **-0.0232** | 0.0060 |
| ot_nolit → intogen_temporal_new | reactome | -0.0045 | 0.0065 |

Two patterns:
1. **Target = intogen_temporal_new** with **seed = ot_all/ot_nolit**:
   the 152 new drivers are a small, recent set; Open Targets'
   general cancer ranking is less informative about this narrow
   target.
2. **Network = reactome** with any temporal_new source: Reactome's
   3,553-node LCC is much smaller and gives only 63-77 unique
   temporal_new seeds.

Both observations are consistent with the D-09 / D-11 control
behaviour reported in A's `results_stats.md`.

---

## §F. Results — Publication-bias adjustment

### β values

| Label | n_positives | log10 mean pos | log10 mean neg | β |
|---|---:|---:|---:|---:|
| intogen2024     | 633 | +3.85 | +2.57 | **+1.282** |
| ot_all          | 600 | +3.81 | +2.48 | **+1.333** |
| ot_nolit        | 600 | +3.71 | +2.48 | **+1.232** |
| clingen_label   |  84 | +3.96 | +2.35 | **+1.610** |

A β of +1.3 in log10 space means positives have ~20× more PubMed
papers than the median gene. clingen_label has the strongest bias
(β=+1.61 = ~40× more papers), as expected — ClinGen's curation
literally requires published evidence.

### Headline: ΔAUROC before vs after adjustment

23 cells compared (1 missing — `ot_all × string_phys700` — because
that original arm was never committed; the comparison skips it).

| Network × Label | ΔAUROC original | ΔAUROC adjusted | diff |
|---|---:|---:|---:|
| string_phys700 × ot_nolit     | +0.170 | +0.074 | **−0.096** |
| reactome × ot_all             | +0.109 | **−0.005** | **−0.114** |
| string_full700 × ot_nolit     | +0.127 | +0.042 | **−0.085** |
| string_full700 × ot_all       | +0.119 | +0.048 | −0.072 |
| funmap × intogen2024          | +0.100 | +0.032 | −0.068 |
| reactome × clingen_label      | +0.206 | +0.138 | −0.068 |
| string_full700 × intogen2024  | +0.113 | +0.054 | −0.059 |
| string_phys700 × clingen_label | +0.170 | +0.116 | −0.054 |
| string_phys700 × intogen2024  | +0.119 | +0.092 | −0.027 |
| string_full700 × clingen_label | +0.097 | +0.062 | −0.035 |
| rna_coexp × ot_nolit          | +0.085 | +0.025 | −0.060 |
| **rna_coexp × intogen2024**   | +0.050 | **+0.081** | **+0.032** |
| reactome × intogen2024        | +0.044 | +0.015 | −0.029 |
| reactome × ot_nolit           | +0.098 | +0.020 | −0.078 |
| funmap × ot_all               | +0.085 | +0.059 | −0.026 |
| intact × clingen_label        | +0.114 | +0.075 | −0.039 |
| intact × ot_nolit             | +0.078 | +0.056 | −0.022 |
| intact × ot_all               | +0.084 | +0.081 | −0.002 |
| **intact × intogen2024**      | +0.065 | **+0.073** | **+0.008** |
| **funmap × ot_nolit**         | +0.078 | **+0.080** | **+0.002** |
| **funmap × clingen_label**    | +0.069 | **+0.118** | **+0.049** |
| rna_coexp × ot_all            | +0.078 | +0.052 | −0.026 |
| rna_coexp × clingen_label     | +0.065 | +0.039 | −0.026 |

**Summary:**
* 19 of 23 cells show ΔAUROC **drops** after PubCount adjustment.
* 4 of 23 cells show ΔAUROC **stays the same or rises**.
* **Mean diff = −0.039** (about 39 % of typical ΔAUROC magnitude).
* Adjusted ΔAUROC remains significantly positive in 22/23 cells;
  one cell (`reactome × ot_all`) drops to noise (−0.005).

### Notable survivors (network prior GAINS signal after adjustment)

Four cells see ΔAUROC unchanged or larger after adjustment:
* `rna_coexp × intogen2024`: 0.050 → 0.081 (**+0.032**)
* `intact × intogen2024`: 0.065 → 0.073 (**+0.008**)
* `funmap × ot_nolit`: 0.078 → 0.080 (+0.002, within noise)
* `funmap × clingen_label`: 0.069 → 0.118 (**+0.049**)

This is interesting and not predicted by a naive bias-only story:
some network priors are themselves biased *against* over-cited
genes, and removing the publication-bias label exposes this
anti-bias structure.

### The reactome × ot_all collapse

The biggest single-cell drop is `reactome × ot_all`: ΔAUROC goes
from +0.109 to **−0.005** — after PubCount adjustment, the
network prior is no better than degree ranking. Clean signal: for
ot_all genes, reactome is essentially just recovering "what's been
studied". After adjustment it has no incremental value.

By contrast, `string_phys700 × ot_all` would likely survive —
STRING has experimental channels (binding, co-expression) that
don't track literature volume.

---

## §G. Discussion — Label bias

The cross-catalogue transfer grid (Control 5) and the
publication-bias adjustment (Control 7) jointly address the
two distinct sources of label circularity that the protocol §1.1
flagged as "structural, not accidental".

### Cross-catalogue: partial catalogue-independent signal

The 50-70 % cross-catalogue transfer is the most informative
single control in the benchmark. It tells us that **most, but
not all, of RWR's contribution is catalogue-independent**.
Specifically:

* The 0.077-0.13 ΔAUROC on ot_all ↔ ot_nolit mutual transfer
  shows that **literature-driven Open Targets rankings share
  substantial cancer-relevant biology** with the no-literature
  variant — the literature channel is not the whole story.

* The 0.06-0.08 ΔAUROC on intogen_temporal_new → intogen2024
  shows that the **152 newly-discovered drivers localise the
  rest of the 633-driver catalogue** with the same magnitude as
  ot_all → intogen2024 (~0.077). New discovery and old discovery
  are encoded in the same network topology.

* The 0.04-0.06 ΔAUROC on ot_all → intogen2024 (and reverse)
  shows that **mutation-based curation (IntOGen) and clinical-
  evidence curation (Open Targets) agree to roughly 50 %** of
  each other's network signal. This is the clearest evidence
  that network topology encodes something beyond either
  curation philosophy.

### Publication bias: real but partial

The 39 % mean drop in ΔAUROC after PubCount adjustment is
consistent with the protocol's H3 prediction. The survival of
22/23 cells means RWR carries real cancer-relevant signal beyond
literature volume.

Two interpretations of the 4 "survivor" cells (where ΔAUROC
rises after adjustment):

1. **The anti-bias hypothesis**: some networks (rna_coexp,
   intact) are themselves biased *against* over-cited genes.
   This is most striking for `funmap × clingen_label` (+0.049
   survivor) — FunMap is a hand-curated functional map whose
   features may not track PubMed volume.

2. **The noise hypothesis**: the survivor cells are within
   1-σ of zero effect after adjustment, and the positive diff
   is fold-noise rather than a real signal. A larger
   replicate (e.g., 50-fold × 10-repeat instead of 50-fold ×
   1-repeat) would resolve this.

The protocol-level recommendation: report both the original
ΔAUROC and the adjusted ΔAUROC side-by-side, and let the reader
judge the magnitude of the bias effect.

### Implications for cancer-driver methods

The bias-adjusted ΔAUROC is the more conservative estimate of
"what the network prior genuinely contributes". For most
networks, adjusted ΔAUROC is 50-90 % of original ΔAUROC. The
adjusted estimate is the one to compare across methods in a
benchmark setting.

### Implications for downstream methods

For downstream methods that use network priors as input (e.g.,
gene-set enrichment, drug-target prioritisation), the
publication-bias-adjusted ΔAUROC suggests:

* Network priors carry meaningful non-bibliometric signal.
* Using degree alone (without the network topology) loses about
  half of the available cancer signal.
* For literature-rich labels (ot_all, ot_nolit), the network
  prior contribution is partly a popularity proxy and partly a
  topology signal; the protocol's design separates these.

---

## §H. Discussion — ClinGen as 5th label

The optional 5th label (ClinGen) has the highest meta-mean
ΔAUROC of any label (+0.1100 across 6 networks), driven by the
strongest single-cell in the panel (`reactome × clingen` =
+0.2058).

### Why ClinGen is interesting

ClinGen is the only label curated by *clinical geneticists*
rather than by computational driver-detection algorithms.
ClinGen's curation rules require:

* Published clinical case reports linking the gene to the
  disease.
* Expert review by a ClinGen Variant Curation Expert Panel.
* Definitive or Strong classification requires multiple
  independent reports with concordant evidence.

This makes ClinGen the most "manually validated" of the five
labels. If RWR helps predict ClinGen cancer genes, then the
network prior is helping predict genes whose driver status a
human expert has independently verified.

### The reactome × clingen cell

The strongest cell in the panel — `reactome × clingen_label` =
+0.2058 — is biologically interpretable. Reactome catalogues
curated pathway interactions; ClinGen curates clinically
validated cancer gene-disease pairs. The pathways most strongly
implicated in cancer (cell cycle, DNA repair, RTK signaling)
are exactly the pathways in which ClinGen-curated genes are
over-represented. RWR on Reactome therefore propagates from
"any cancer-relevant gene in the LCC" to "the genes in the
specific pathway ClinGen curates" — a direct retrieval.

### The clingen label and H3

ClinGen's β = +1.610 is the largest in the panel, confirming
that ClinGen curation tracks publication volume (clinicians
only curate genes they've seen in the literature). Yet the
*adjusted* ClinGen label still beats degree in 5 of 6 networks,
with `funmap × clingen_label` gaining +0.049 after adjustment
(the largest survivor gain in the panel). This is consistent
with the interpretation that ClinGen's residual signal is
*clinical* (gene-disease validity from clinical case reports)
rather than *literature* (PubMed counts of any kind).

---

## §I. Cross-references

* A-side logs: `docs/a0_log.md`...`docs/a5_log.md`, `docs/w5_log.md`,
  `docs/gnn_courtesy_log.md`, `docs/stats_analysis_log.md`,
  `docs/s2_log.md`.
* B-side log: `docs/b_log.md`.
* Cross-catalogue transfer report: `docs/transfer_results.md`.
* Publication-bias report: `docs/pubcount_adj_results.md`.
* Label-network overlap report: `docs/overlap_report.md`.
* Decisions log: `docs/decisions.md`.
* Pre-registration: `docs/pre_registration.md`.
* Public release protocol: `docs/public_release_protocol.md`.
* Reproducibility audit: `docs/reproducibility_audit.md`.
* Pending-work analysis: `docs/pending_work_protocol_compliance.md`.
