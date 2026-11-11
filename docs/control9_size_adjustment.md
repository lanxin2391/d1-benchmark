# Control 9 — Network size and coverage adjustment (Wright et al. 2025)

**Purpose.** Test whether the primary endpoint (ΔAUROC = AUROC(RWR) −
AUROC(degree)) is a **network-size** or **label-coverage** artefact.

**Why it is needed.** Wright et al. (2025, *Molecular Systems Biology*
21:e10077) evaluated 45 human interactomes and found that the leading
predictor of a network's disease-gene prioritisation performance is its
**size**, and that the effect largely collapses once the **coverage** of
the gene set in the network is controlled for:

> "performance correlates with network size, but the magnitude of the
> effect is significantly reduced after controlling for gene set
> coverage" (p = 8.8e-34 for the literature-derived gene sets,
> 4.1e-6 for the genetic ones)

They adjust by regressing on log₁₀ of the interaction count and taking
the residual ("size-adjusted performance").

In D1 this is the most dangerous confound in the design: the six
networks span 3,953 to 17,918 LCC nodes (4.5×) and 41 % to 99 % label
coverage (2.4×), and provenance is correlated with both.

**Design.**

1. The per-cell primary endpoint at α = 0.5 comes from
   `results/tables/a4_table1_delta_auroc_v40.tsv` (24 cells).
2. The covariates come from `data/processed/networks/*.meta.json`
   (LCC edges) and `results/tables/label_network_overlap.tsv` (the
   fraction of the label's positives present in the network's LCC).
3. Ordinary least squares:

   ```
   delta_auroc ~ 1 + label_set + log10(LCC edges) + coverage
   ```

   The label set is a **fixed effect**, because the four label sets are
   not comparable across rows (one of them has 152 positives and
   another 633).
4. The **fitted value** is the size- and coverage-adjusted endpoint;
   the **residual** is how far a cell sits from what its size and
   coverage predict.
5. The size correlation is additionally tested against a
   10,000-permutation null: the Spearman correlation between ΔAUROC and
   log₁₀(LCC edges), pooled over the 24 cells and over the 4 labels
   separately.

**Results** (`results/tables/control9_size_model.tsv`,
`results/tables/control9_size_adjustment.tsv`):

| term | coefficient | SE | p | partial R² |
|---|---|---|---|---|
| log₁₀(LCC edges) | **−0.0393** | 0.0142 | 0.0129 | 0.297 |
| label coverage | **+0.1449** | 0.0363 | 0.00086 | 0.469 |
| label = intogen_temporal_new | −0.0528 | 0.0142 | 0.0016 | 0.435 |
| label = ot_all | +0.0235 | 0.0140 | 0.111 | 0.135 |
| label = ot_nolit | +0.0273 | 0.0140 | 0.068 | 0.173 |

| quantity | value |
|---|---|
| pooled Spearman ρ(ΔAUROC, log₁₀ edges) | **−0.021** |
| permutation p (10,000 permutations) | **0.93** |
| within-label Spearman (6 networks each) | +0.257, +0.257, −0.314, −0.429 (all p > 0.39) |
| raw endpoint positive | 23/24 cells |
| size/coverage-adjusted endpoint positive | **23/24 cells** |
| sign preserved by the adjustment | 22/24 cells |
| mean \|residual\| | 0.0168 (raw SD 0.0462) |

**Verdict — two findings, pointing in opposite directions, both
reported.**

1. **The headline result is not a size artefact.** Wright et al.'s
   concern does not materialise here. A bigger network does *not* give a
   bigger margin (pooled ρ = −0.02, permutation p = 0.93), and the
   adjusted endpoint stays positive in 23/24 cells — the same count as
   the raw endpoint.
2. **But size and coverage are individually significant, and coverage
   is the strongest single covariate** (partial R² = 0.47, p = 8.6e-4).
   Coverage, not size, is what moves the endpoint: an arm is easier
   when more of its positives are present in the network. This is the
   mechanism behind the Reactome temporal cell (Results §F.10, §J.3) —
   Reactome's LCC contains only 41 % of the newly-discovered drivers,
   the lowest in the panel — and it is also the reason the provenance
   hypothesis (H2) cannot be settled by this panel (§J.2): a panel in
   which provenance is separated from size and coverage needs more
   networks per class than three, one and two.

**Cost.** Seconds (the endpoint is already computed; only the regression
runs).

**Reproduce.**

```bash
python scripts/run_control9_size_adjustment.py
```

**Status.** Exploratory family (amendment 5). It does not enter the H1
decision rule. It does bear on H2, and the manuscript uses it to state
*why* H2 is not testable rather than to claim a provenance result.
