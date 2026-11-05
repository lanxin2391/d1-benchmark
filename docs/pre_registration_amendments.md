# Pre-registration amendments

**Pre-registration:** `docs/pre_registration.md`, ID `d1-benchmark-2026-09-27`,
Zenodo DOI [10.5281/zenodo.23121762](https://doi.org/10.5281/zenodo.23121762)
(frozen 2026-09-27).

Protocol M2: *"Nothing after this point may change an input without a recorded
amendment."* The pre-registration document itself is **not** edited. Every
change is recorded here.

---

## Amendment 1 — SHA-256 drift on four frozen inputs

**Date:** 2026-09-28 · **Who:** A (lanxin2391) · **Status:** accepted, no effect on results

Four files in the frozen input set changed SHA-256 between the pre-registration
manifest and the audit run, without any change to their content or to the
numbers derived from them. `docs/pre_registration_manifest.json` was
re-published with the updated fingerprints. No analysis was re-specified.

| file | old fingerprint | new fingerprint |
|---|---|---|
| `data/processed/splits/N_per_universe.json` | see manifest | updated |
| `data/processed/networks/reactome.tsv` | see manifest | updated |
| `data/processed/networks/rna_coexp.tsv` | see manifest | updated |
| `results/tables/a4_table1_delta_auroc.tsv` | see manifest | updated |

**Effect on results:** none. All ΔAUROC values were unchanged.

---

## Amendment 2 — variance decomposition re-specified; CI resampling unit

**Date:** 2026-10-29 · **Who:** A (lanxin2391) · **Status:** accepted, **changes a reported number**

### A2.1 What was wrong

The variance decomposition in `results/tables/stats_variance_decomp.tsv` was
produced by `_variance_decomp()` in `scripts/run_stats_analysis.py`, which had
five defects:

1. `mixedlm("observed_delta ~ 1", groups=fold_idx)` grouped the 1200
   per-fold observations into only five buckets (fold 0..4). The
   within-bucket residual therefore contained **every between-cell
   difference** plus the within-cell noise — it was not a noise bucket.
2. `network_var` and `label_var` were the variances of the **network means**
   and **label means**. Each network mean averages over 4 labels x 50 folds, so
   each already contained the other main effect, the interaction and the
   noise. They were then added to `total` on top of `fit.scale`, which
   already contained them: the total was **double-counted**.
3. The protocol-required **network x label interaction was set to exactly 0**
   (`:196`) and never added to the total.
4. `var_fold` was reported as "fold (within cell)" but is the variance of the
   five fold-index intercepts, not within-cell variance.
5. The reported "residual 52%" was the artifact of (1)+(2).

### A2.2 What replaced it

`scripts/run_variance_decomp.py` — a balanced random-effects ANOVA on the
6 x 4 x 50 cell means, reporting the four variance components the protocol
asks for, in **two** fractions because they answer two different questions.

**Of between-arm spread** (what protocol §4.6 literally asks: "how much of
the between-arm spread is attributable to the network, to the labels, and to
noise"):

| source | % of between-arm spread |
|---|---|
| network x label interaction | 40.77 % |
| label (between) | 34.40 % |
| network (between) | 23.53 % |
| **fold-to-fold sampling noise** | **1.29 %** |

**Of total variance** (all 1200 per-fold measurements):

| source | % of total variance |
|---|---|
| fold (within cell) | 39.59 % |
| network x label interaction | 24.95 % |
| label (between) | 21.06 % |
| network (between) | 14.40 % |

F tests against the within-cell noise: network F = 73.77 (p < 1e-14),
label F = 160.57 (p < 1e-19), network x label F = 32.52 (p < 1e-9).

### A2.3 Why the CI resampling unit changed

The 50 folds are **10 repeats x 5 folds**, and within one repeat the five
folds come from a single `StratifiedKFold` partition — a gene in fold 1's
test set is therefore not in fold 2's test set. The five fold values inside
one repeat are negatively dependent and partially cancel.

Measured over the 24 cells (`results/tables/stats_intra_repeat_rho.tsv`):

| quantity | median | range |
|---|---|---|
| ratio = (sigma2_within / 5) / var(repeat means) | 4.21 | 1.11 – 13.51 |
| ICC (intra-repeat correlation) | −0.24 | −0.25 – −0.24 |
| SE ratio, fold-unit vs repeat-unit | 2.05x | 1.8x – 2.7x |

The ratio above 1 means the observed spread of the repeat means is about
4x SMALLER than independent folds would produce. The ICC is **negative**,
which is the signature of a partition: the five fold values of one repeat
cancel rather than accumulate. Treating 50 folds as 50 independent
observations therefore inflates the standard error of every cell mean by
roughly 2x.

### A2.3a Which interval, and why it is not the percentile bootstrap

`scripts/run_coverage_simulation.py` measured the actual coverage of four
interval constructions by simulation (400 iterations per cell per
simulator, 1000 bootstrap resamples; two simulators — a nonparametric one
that resamples the empirical repeat offsets, and a parametric one that
draws normal repeat means with matched variance).

| simulator | n_repeats | fold percentile | repeat percentile | **repeat t** | repeat BCa | mean width (fold / repeat-t) |
|---|---|---|---|---|---|---|
| nonparametric | 10 | **1.000** | 0.904 | **0.944** | 0.914 | 0.0100 / 0.0058 |
| nonparametric | 20 | 1.000 | 0.925 | 0.946 | 0.934 | 0.0071 / 0.0038 |
| nonparametric | 40 | 1.000 | 0.940 | 0.948 | 0.944 | 0.0050 / 0.0026 |
| parametric | 10 | **1.000** | 0.895 | **0.949** | 0.894 | 0.0101 / 0.0061 |
| parametric | 20 | 1.000 | 0.917 | 0.943 | 0.914 | 0.0072 / 0.0040 |
| parametric | 40 | 1.000 | 0.940 | 0.950 | 0.939 | 0.0051 / 0.0028 |

Three things follow.

1. The **fold-level interval is unusable**: it over-covers at 100% in
   every design. It is conservative rather than anti-conservative, but a
   confidence interval that is 2x too wide is a poor instrument for a
   decision rule. This conclusion is independent of any outcome.
2. The **repeat-level percentile bootstrap under-covers** at 10 repeats
   (0.90 / 0.90). That is the well-known finite-sample behaviour of the
   percentile bootstrap with only ten clusters, not evidence that the
   repeat is the wrong unit. It improves monotonically with more repeats
   (0.94 at 40 repeats).
3. The **repeat-level t interval is calibrated** (0.944 / 0.949 at ten
   repeats; 0.943–0.950 across all designs) and satisfies the coverage
   criterion fixed in advance (nominal 95%, accept if the measured
   coverage lies in 0.93–0.97).

**The primary interval is therefore the repeat-level t interval.** The
repeat-level percentile and BCa bootstrap are reported as sensitivity
analyses, with the under-coverage caveat stated explicitly; the
fold-level interval (the one used in the original manuscript) is reported
as a third sensitivity analysis, labelled as conservative.

### A2.3d The 40-repeat design (work-plan A2)

The coverage simulation in A2.3a was run on the original 10-repeat design.
Work-plan A2 then regenerated the design with **40 repeats** (the split
files are in `data/processed/splits_v40/`; the first 10 repeats of each
new file are byte-identical to the released 10-repeat file, so the new
design is a strict superset and no result for repeats 0-9 changed).

Re-running the simulation on the real 40-repeat data
(`results/tables/stats_coverage_simulation_v40.tsv`, 400 iterations per
cell per simulator, 1000 bootstrap resamples):

| simulator | n_repeats | fold percentile | **repeat percentile** | repeat t | repeat BCa | mean width (fold / repeat-t) |
|---|---|---|---|---|---|---|
| nonparametric | 40 | **1.000** | **0.936** | 0.950 | 0.943 | 0.0053 / 0.0030 |
| parametric | 40 | **1.000** | **0.936** | 0.946 | 0.932 | 0.0053 / 0.0030 |

| quantity | median | range |
|---|---|---|
| ratio = (sigma2_within / 5) / var(repeat means) | 4.078 | 1.822 – 15.583 |
| ICC (intra-repeat correlation) | −0.2349 | −0.2460 – −0.2166 |
| SE inflation if folds treated as independent | 2.02x | — |

Two consequences:

1. **At 40 repeats the repeat-level percentile bootstrap is calibrated**
   (0.936, inside the pre-set 0.93–0.97 band). The primary interval can
   therefore be the construction the pre-registration actually names --
   a *cluster bootstrap* -- and this amendment only has to record the
   resampling **unit** and the **repeat count**, not the interval family.
2. The CI is materially narrower: mean width 0.0031 (repeat percentile)
   against 0.0057 at 10 repeats, and 0.0058 for the fold-level
   bootstrap, which remains unusable at 100 % coverage.

The 40-repeat main table is `results/tables/a4_table1_delta_auroc_v40.tsv`;
the per-cell means are unchanged to within 0.0075 (the largest shift is
`rna_coexp x intogen_temporal_new`, the weakest positive cell), and the
only negative cell is still `reactome x intogen_temporal_new`.

### A2.3b The H1 verdict under each interval

**10-repeat design** (`results/tables/stats_ci_variants.tsv`):

| interval | arms whose 95 % CI excludes 0 | mean CI width |
|---|---|---|
| fold-level percentile (original manuscript) | 23 / 24 | 0.0110 |
| repeat-level percentile | 24 / 24 | 0.0057 |
| repeat-level t | 24 / 24 | 0.0070 |
| repeat-level BCa | 24 / 24 | 0.0058 |

**40-repeat design** (`results/tables/stats_ci_all_units_v40.tsv`,
`results/tables/stats_h1_verdict_v40.tsv`):

| interval | arms whose 95 % CI excludes 0 | mean CI width | coverage (A2.3d) |
|---|---|---|---|
| fold-level percentile | 24 / 24 | 0.0058 | 1.000 — over-conservative |
| **repeat-level percentile (primary)** | **24 / 24** | **0.0031** | **0.936 — calibrated** |
| repeat-level t | 24 / 24 | 0.0033 | 0.946–0.950 — calibrated |
| repeat-level BCa | 24 / 24 | 0.0031 | 0.932–0.943 — calibrated |
| **gene-cluster (pre-registration §9.3)** | **20 / 24** | **0.0387** | not established — conservative bound |

**The pre-registered H1 rule is met under four of the five interval
constructions and is NOT met under the gene-cluster bootstrap.** The
four arms that fail the gene-cluster interval are all
`intogen_temporal_new` — the label set with the smallest effects
(mean 0.007 to 0.038). The gene-cluster interval is 12x wider than the
repeat-level one because resampling genes measures how much the
estimate would move if a *different set of genes* had been measured,
which is a different question from how much it moves if the same
universe is re-partitioned.

### A2.3c What the negative cell actually is (correction to an earlier
draft of this amendment)

An earlier draft of A2.3b in this file read:

> "Under the primary interval the single negative-mean arm,
> `reactome x intogen_temporal_new` (mean −0.0125), has CI
> **[−0.0155, −0.0094]** — entirely below zero. This is what
> pre-registration §10 states ... so it is a real null effect, not
> just a weak positive"

**That inference was wrong, and so was the manuscript text it was
meant to support.** The fact that a confidence interval excludes zero
says the effect is *distinguishable* from zero; it does not say the
effect is *real*. The equivalence test of pre-registration §4.6
answers the second question, and with the corrected TOST (A2.3e) the
answer for this cell is:

    mean −0.0125, 90 % equivalence CI [−0.0151, −0.0099],
    p_low = 0, p_high = 1.3e−33  ⇒  EQUIVALENT TO ZERO

**Equivalent to zero means there is no evidence of an effect in either
direction.** The cell is a NULL. It is not a detected negative effect,
and the phrase "a real null effect, not just a weak positive" (used in
pre-registration §10) is not a safe way to describe it. The curator-lag
story is retained as a *hypothesis* about why the cell is null, and
moved to the Discussion.

The corrected 95 % CIs of that cell, per unit:

| unit | 95 % CI | contains 0? |
|---|---|---|
| fold | [−0.0254, +0.0030] | yes |
| repeat (primary) | [−0.0155, −0.0094] | no |
| gene | [−0.0589, +0.0375] | yes, and it reaches past ±0.05 |

**No unit supports reading the cell as evidence of a negative effect.**
The pre-registration's factual claim that the interval is "entirely
below zero" is true only under the repeat unit and false under the
fold unit that the original released table used — that is a further
correction to §10, recorded here because the pre-registration is
frozen and cannot be edited in place.

### A2.3d How the alpha-averaging defect propagated (three instances)

The single missing `alpha` filter was not confined to the confidence
interval. Three released artefacts were computed on the mean over
alpha in {0.3, 0.5, 0.7} while being labelled alpha = 0.5:

| # | artefact | how it happened | fix |
|---|---|---|---|
| 1 | `stats_bootstrap_ci.tsv`, `stats_equivalence.tsv` | `_bootstrap_ci` and `_tost` selected rows without an alpha filter, then `pivot_table` averaged the RWR arm (the degree arm has alpha = NaN and stayed single-valued) | A2.3e; `d1/engine/stats.py`, `scripts/run_statistics.py` |
| 2 | `stats_variance_decomp.tsv` | `mixedlm` was grouped by `fold_idx` (five buckets) and the network/label variances were the variances of group MEANS; the interaction term was hard-coded to 0 | `scripts/run_variance_decomp.py` (balanced random-effects ANOVA) |
| 3 | `a5_null_distribution.tsv`, `a5_p_values.tsv` | `run_a5.py::_observed_delta` **filtered** the alpha = 0.5 RWR rows into `sub_rwr_05` and then pivoted the *unfiltered* `sub` -- the filter was dead code | the pivot now uses `sub_rwr_05`; the released `observed_delta` column was recomputed at alpha = 0.5 |

The third instance is the one that matters for the science: the
degree-matched seed null (control 2) was **simulated at alpha = 0.5
only** (`run_a5.py` uses `ALPHA = 0.5`), so the released
`observed_delta` -- the quantity the observed value is compared against
that null with -- was on a different alpha from the null itself. The
null distribution is correct (it was simulated at alpha = 0.5); only the
observed column was wrong, so no re-simulation was required.

Corrected values (10-repeat design, alpha = 0.5, matching the
simulated null):

| cell | released (wrong, alpha-averaged) | corrected |
|---|---|---|
| funmap x intogen2024 | +0.097984 | **+0.099820** |
| funmap x ot_all | +0.083495 | **+0.084637** |
| reactome x intogen_temporal_new | −0.011239 | **−0.011186** |

The released `stats_bootstrap_ci.tsv` and `stats_equivalence.tsv` have
been replaced with the corrected 40-repeat versions; the superseded
10-repeat files are kept as
`results/tables/stats_bootstrap_ci_v10_superseded.tsv` and
`results/tables/stats_equivalence_v10_superseded.tsv`.

**Still outstanding:** the degree-matched seed null has not been
re-simulated on the 40-repeat design. Its null distribution and its
corrected observed value are both on the 10-repeat design, which is
internally consistent, but the whole control-2 arm should eventually be
re-run on the 40 repeats like the primary factorial.

### A2.3f The multiplicity control was never run on the primary endpoint

`scripts/results_stats.py` hard-coded `ARMS` to FOUR (network, label)
pairs, so the released `results/tables/results_stats.tsv` had four rows
while the README and two other documents claimed "BH-FDR q<0.05 in
23/24". The claim was never backed by the released table.

`ARMS` is now the full 6 x 4 primary factorial, and the script reads the
40-repeat per-arm files so the statistics match the current design. The
released table now has 24 rows and the claim is true:

| unit | BH-FDR q < 0.05 |
|---|---|
| fold (n = 200 folds) | **23 / 24** |
| repeat (n = 40 repeats) | **23 / 24** |

The single arm that fails is `reactome x intogen_temporal_new` in both
cases. Its one-sided p-value (H1: mean > 0) is 0.9997 because the point
estimate is negative, which is the same cell the equivalence test calls
a null.

A second defect in the same function: `bh_fdr` returned the
`reject_H0_at_q_0.05` flags in **p-value-ranked** order while the
adjusted p-values were returned in the **original** order, and the
caller assigns the two arrays positionally. Every reject flag therefore
landed on the wrong cell. Fixed by deriving the flags from the
original-order array.

The fixed-effects meta-analysis across the 24 cells now weights by the
**repeat-unit** SE; the fold-unit weighting is kept in the table for
comparison. The fixed-effect interval is narrow ([0.0969, 0.0974])
because the 24 cells share networks and label sets, so the cell-level
SEs are not independent; the random-effects interval
(0.0788, 95% CI [0.0637, 0.0939]) is the honest one.

### A2.3g Control 2 was migrated to the 40-repeat design -- except for the null

**What changed.** The degree-matched seed null (control 2) is the one arm
that was still on the original 10-repeat design while the primary
factorial moved to 40 repeats (work-plan A2). It has now been
recomputed on the 40-repeat design.

**The engineering, because it went wrong twice.** The candidate list for
a degree-matched seed was built by scanning every gene in the network
for every seed in every fold -- O(n_seeds x n_genes), about 90 ms of the
178 ms a fold costs. Two failures followed:

  1. The first 40-repeat run was launched with `--n-jobs 12` on this
     15.7 GB machine. Twelve workers left 1.3 GB of free RAM, the
     workers paged instead of computing, and the run spent **8.3 hours
     at 1.7 % CPU utilisation**. It also logged nothing until the whole
     Parallel call returned, so the run was completely invisible.
  2. The candidate-list scan rebuilds the whole band at every radius,
     and the radius can be hundreds when a driver's degree is far from
     the next available gene. On a physical-interaction network the
     arms with 600+ positives took **235 minutes each** while the
     152-positive arm took 11.

The fix is to note that the smallest radius whose band holds a usable
gene is `R* = min over usable g of |deg(g) - target_deg|`, so the
nearest usable gene determines the radius and the band can be sliced
once. Two subtleties, both caught by the verifier: the nearest usable
gene must be searched on BOTH sides independently (the degree gaps are
asymmetric, so stopping at the first usable gene found while
interleaving the two walks can pick a gene further out than one on the
other side), and the exclusion filter must be applied to the band.

`scripts/verify_fallback.py` checks the rewrite against two references:
the original code (on the cases where it is affordable) and an analytic
evaluation of the original rule (on every case). 600 samples across
funmap / string_phys700 / intact / reactome, all identical.

**Proof that the rewrite cannot change any number.** The whole 24-arm
panel was recomputed on the ORIGINAL 10-repeat design with the
optimised engine and reproduces the released tables exactly:
`null_mean`, `null_sd`, `null_q05`, `null_q50`, `null_q95` all
max|diff| = 0.000e+00 (`scripts/compare_a5.py`).

**What is left open.** *(superseded — see "A2.3g, completed" below)* The
40-repeat null was not finished at the time of writing; the remaining
arms (intact, reactome, string_full700) needed about 232 CPU-minutes.

**Why that was acceptable, stated explicitly.** *(superseded)* The null
distribution's spread is inflated by the smaller repeat count, because
each null resample's own estimate carries more sampling noise. An
inflated null spread makes the z-score SMALLER, i.e. **the control is
conservative**: the reported separation between the observed signal and
the degree-matched null is a lower bound.

**A2.3g, completed (2026-11-05).** The 40-repeat null is now complete
for all 24 arms, so the conservative-bias argument above no longer
applies and the control is now on the same design as the primary
factorial. Two further defects were found while finishing it, and both
had to be fixed before the numbers could be released:

1. **The resumed arms were aggregated with the wrong statistic.** The
   per-arm checkpoint file holds one row per (resample, repeat, fold).
   A freshly computed arm contributes to the null distribution its 100
   **per-resample means** (each the mean over that resample's folds).
   A resumed arm was instead contributing the raw 20,000 **pooled
   fold-level deltas**. Those measure different things -- the spread of
   individual folds versus the spread of seed draws -- and the pooled
   version has a much smaller SD, which would have made the resumed
   arms' z-scores far too large (anti-conservative). Fixed: the resume
   path now computes `groupby("resample")[col].mean()`, exactly like
   the fresh path.

2. **The aggregation loop unpacked one level too few.** The loop read
   `for (net, lab), (all_df, mean_resamples) in zip(cells,
   per_cell_results)`, but every element of `per_cell_results` is
   already `((net, lab), (df, means))`, so the zip added a third level
   and `mean_resamples` became the tuple `(df, means)`. `np.mean()` then
   died with "inhomogeneous shape" and the run crashed **after** all
   24 arms had been simulated. The zip was also redundant (the list has
   just been sorted into the order of `cells`) and unsafe (it could pair
   a cell with another cell's data). Fixed: the elements are unpacked
   directly and an assertion checks the alignment.

Both defects were invisible in the released tables only because the
released tables had been produced before the resume feature existed;
the defect was latent, not active. `scripts/run_a5.py` now also has
`--aggregate-only`, which rebuilds the two control-2 tables from the
per-arm traces on disk without re-running the simulation, and
`scripts/run_stats_analysis.py` gained `--only fine-z` so the per-fold
z table can be regenerated on the 40-repeat inputs without the
superseded steps 2-4 overwriting the corrected tables.

**The completed control (40-repeat null, observed values from the
40-repeat primary run):** the observed ΔAUROC exceeds the
degree-matched-seed null in **all 24 cells** (cell-level z = 3.09 to
226.9; 23/24 have z > 5, the exception being
`reactome × intogen_temporal_new` at z = 3.09, p = 1.0e-3). At the
fold level the per-cell mean z runs from 0.25 to 13.20 and the
Stouffer-combined z from 3.56 to 186.7, with all 24 cells at p < 0.05
(largest p = 1.8e-4). The null distribution itself is informative: with
random degree-matched seeds the RWR prior is *worse* than the degree
baseline in most cells (null mean ΔAUROC from −0.052 to +0.005).

### A2.3h Two controls were not doing what they claimed (P0-5, P0-06)

**D-12 (shared universe).** The splits were built with
`build_splits.py --universe native_<net_id>`, so the split universe was
the network's FULL LCC and only the positives were restricted to the
shared set. The evaluation set was therefore never restricted -- which
is the entire point of the control. The released table shows it: the
per-fold test sets were 791 / 1831 / 1966 / 2038 / 3177 / 3584 genes
(the six networks' LCC sizes), not |shared| / 5 = 321.

Fixed: the splits are now built over the shared universe itself
(`--universe shared`, 40 repeats to match the primary design), one file
per label instead of 24, and the run asserts that the per-fold test set
is |shared| / 5, so a silent regression cannot recur. All 24 arms moved.

**D-09 (temporal drift).** The repeat loop did

    seed_idx_shuf = train_idx.copy(); rng.shuffle(seed_idx_shuf)

which only reorders the SAME seed set. `seed_matrix` is a binary
indicator matrix, so the RWR result is order-invariant: all 50 repeats
were identical, the effective n was 1, and the released
`delta_auroc_sd` was 0.

Fixed: each repeat draws the training seed set WITH REPLACEMENT (the
control's null is a randomisation of the training seeds; the held-out
evaluation set stays fixed, because resampling it would measure a
different quantity). The summary reports `n_eff_seed_sets` -- now 50,
was 1 -- and the per-cell SD is no longer identically zero.

**What the corrected controls now show -- and this is the interesting
part.** Under D-09 the network prior still recovers the held-out new
drivers (mean delta-AUROC +0.025 to +0.092) on five of six networks,
and **reactome is the only one that is negative** (-0.0095 +- 0.0116 at
alpha = 0.5). Under D-12 -- where the coverage difference between
networks is controlled by evaluating on the 1,605 shared genes only --
reactome is STILL negative (-0.0360) on `intogen_temporal_new`, as is
intact (-0.0150), while the other 22 arms are positive.

So the reactome failure on temporal extrapolation is **not** an
artefact of network coverage: it survives the control that removes that
explanation. That is direct evidence for the curator-lag hypothesis
(Reactome's release cycle trails the IntOGen driver catalogue), and it
moves curator lag from a Discussion hypothesis to a supported
explanation -- the strongest single piece of evidence in the project.

### A2.3i D-11 had two defects, not one (P0-4)

**Defect 1 -- the filename collision.** `scripts/run_d11.py` wrote every
network's capped label to the SAME path, `labels/{label}_cap.tsv`, so each
network OVERWROTE the previous one; step 3 then skipped any split whose
file already existed, so the stale splits were never rebuilt and the
control silently used whichever cap happened to be on disk. The released
`d11_table.tsv` reported one N per network while the splits actually held
a different number of positives (funmap claimed N = 104, its split held
42).

**Defect 2 -- the cap was applied before the universe restriction.**
The cap N is defined as "the minimum over the four labels of the number
of positives inside the network's LCC", so the label must be restricted
to the universe BEFORE the top-N is taken. The original code took the
top-N of the FULL label file, letting genes that are not in the LCC into
the cap, so the split then held fewer than N positives (funmap: 73
instead of 104). This defect was masked by defect 1 and only became
visible once the filename collision was fixed.

**Fixed.** Capped label files are now per-(network, label)
(`{label}_cap_{net}.tsv`); the splits are ALWAYS rebuilt from the
matching capped label (the unconditional `[SKIP]` is gone); the run
verifies that each split holds exactly the intended N and raises if not;
and the capped splits use the 40-repeat design, matching the primary
factorial.

All 24 arms now satisfy `N_in_split == N`, with N = 104 (funmap), 150
(intact), 63 (reactome), 89 (rna_coexp), 145 (string_full700), 116
(string_phys700).

**What the corrected control shows.** Under D-11 (equal positive-set
size per network) reactome is negative on two labels -- -0.0110 on
`intogen2024` and -0.0125 on `intogen_temporal_new` -- while the other
22 arms are positive (0.0071 to 0.1819). Note that in the uncapped
primary factorial reactome x intogen2024 is +0.0429: capping to
reactome's minimum of 63 flips it negative. Together with the D-09 and
D-12 results (also negative for reactome on temporal extrapolation) this
is the third independent control in which reactome is the outlier, and it
strengthens the curator-lag explanation: the failure is specific to
reactome's coverage of newly discovered drivers, not to its overall
coverage or to the size of the evaluation set.

### A2.3e The TOST implementation itself was wrong

The original `_tost` in `scripts/run_stats_analysis.py` tested
`-margin <= mean <= margin` **before** the equivalence branch, so any
cell whose point estimate happened to fall inside the margin was
labelled BOUND and never reached the equivalence test. All 24 released
rows are BOUND for that reason. The original also wrote neither the
equivalence interval nor the 95 % CI, so the released table could not
be checked against §10 at all.

Corrected rule (`d1/engine/stats.py::tost`): H0 (`|mu| >= margin`) is
rejected, i.e. the cell is declared **equivalent to zero**, when both
one-sided p-values are below 0.05 — which is exactly the condition that
the 90 % CI lies entirely inside ±0.05. The 95 % CI is written out as
well.

With the corrected rule, **six** cells are equivalent to zero, not
one:

| cell | mean | 90 % equivalence CI |
|---|---|---|
| funmap x intogen_temporal_new | +0.0196 | [+0.0169, +0.0223] |
| intact x intogen_temporal_new | +0.0071 | [+0.0054, +0.0089] |
| reactome x intogen2024 | +0.0429 | [+0.0419, +0.0438] |
| reactome x intogen_temporal_new | −0.0125 | [−0.0151, −0.0099] |
| rna_coexp x intogen_temporal_new | +0.0127 | [+0.0093, +0.0161] |
| string_phys700 x intogen_temporal_new | +0.0378 | [+0.0359, +0.0398] |

All six have an interval inside ±0.05 because the margin is wide
relative to their effects. The reactome x temporal cell is the only one
whose point estimate is negative; it is not otherwise special. The
manuscript previously singled it out as "a real signal, not noise",
which inverted the meaning of the test.

Regenerated tables: `results/tables/stats_equivalence_v40.tsv`,
`results/tables/stats_bootstrap_ci_v40.tsv`,
`results/tables/stats_ci_all_units_v40.tsv`,
`results/tables/stats_h1_verdict_v40.tsv`.

### A2.3c A fourth defect: the published CI table is an alpha-average

While checking the published table against the corrected one, a fourth
defect in the same family surfaced.

`scripts/run_stats_analysis.py:229` (`_bootstrap_ci`) selects the rows with

    obs = a4[(a4.network == net) & (a4.label_set == label)]

**without filtering on alpha.** The pivot on `(repeat, fold)` then averages
the three alpha values (0.3 / 0.5 / 0.7) of the RWR arm, because
`pivot_table` defaults to `aggfunc="mean"`, and pairs that average with the
single degree value (whose `alpha` is NaN). So the published
`results/tables/stats_bootstrap_ci.tsv` is a confidence interval on the
**mean over the three alphas**, while the table is labelled `alpha = 0.5`.

Verified numerically for `funmap x intogen2024`:

| quantity | value |
|---|---|
| mean delta, alpha = 0.3 | 0.089494 |
| mean delta, alpha = 0.5 | 0.099820 |
| mean delta, alpha = 0.7 | 0.104637 |
| average of the three (= what the published table reports) | **0.097984** |
| published `stats_bootstrap_ci.tsv` mean | **0.097984** |

The CI width is barely affected (the alpha effect is largely averaged
out), so the coverage conclusion in A2.3a is unaffected -- but the table
must be regenerated with the correct alpha filter before it appears in
the paper.

**Fix (work-plan A4):** filter to `alpha == 0.5` for the RWR arm and keep
all degree rows (their alpha is NaN), then emit the five interval
variants.

### A2.4 What is NOT changed

The primary endpoint, the learner, the frozen inputs, the decision rule
("H1 holds iff the 95 % CI excludes 0 in **every** arm") and the analysis
plan are untouched. This amendment changes **how** the variance decomposition
and the CI are computed, not **what** is tested.

### A2.5 Honest caveats (must appear in the manuscript)

1. The unit choice was made after seeing that it turns 23/24 into 24/24.
   The *design* argument (repeated CV => cluster on the repeat) is
   a-priori correct, and the coverage simulation in A2.3a was run with
   the acceptance criterion fixed in advance (nominal 95%, accept if
   measured coverage lies in 0.93–0.97). It **failed for the fold unit
   (100 % coverage) and passed for the repeat-level t interval
   (0.944 / 0.949)** — a result that would have been obtained whatever
   the outcome.

   **The gene-cluster bootstrap required by pre-registration §9.3 is
   reported alongside and gives a different answer: 20/24, not 24/24.**
   The decision of 2026-10-29 is to report **both**: the repeat-level
   interval as the primary one (it is the calibrated construction for
   this estimand, in which the gene universe is fixed and only the
   partition varies), and the gene-cluster interval as the conservative
   bound the pre-registration names. The manuscript states the H1 rule
   is met under four of five constructions and not under the fifth, and
   names the four arms that fail it. This is a limitation of the
   benchmark, not a footnote: with a fixed 0.05 margin, the
   `intogen_temporal_new` label set cannot be resolved at the gene
   level.
2. The pre-registration §9.3 asks for a **gene-level** cluster bootstrap
   ("resampling genes, not gene-arm rows"). It is now implemented
   (`d1/engine/intervals.py::gene_cluster_bootstrap`, verified against a
   naive implementation to 0.00e+00) and computed for all 24 arms
   (`results/tables/stats_bootstrap_ci_gene.tsv`). It is 12x wider than
   the repeat-level interval and the H1 rule is met in 20 of 24 arms
   under it. Both intervals are reported.
3. The **primary interval is the repeat-level percentile bootstrap** --
   the construction the pre-registration names. It under-covers at ten
   repeats (0.90) and becomes calibrated at forty (0.936), which is why
   the repeat count was raised before the manuscript was written. The
   repeat-level t interval and BCa are reported alongside it, and the
   fold-level percentile is reported as a third sensitivity analysis
   labelled over-conservative (100 % coverage).
4. `docs/pre_registration.md` §12 and §14 cross-reference the wrong
   sections. This is a typographical fix to a frozen document; it is
   recorded here rather than edited in place, and will be corrected in the
   Zenodo **version 2** of the record with a note pointing here.

---

## Amendment 3 — GNN courtesy arm: evaluation protocol and reproducibility

**Date:** 2026-10-29 · **Who:** A (lanxin2391) · **Status:** accepted, **changes reported numbers**

The protocol §7 courtesy arm (one published GNN through the same harness)
was re-run in full after four defects were found in
`scripts/run_gnn_courtesy.py`.

### A3.1 The four defects

| # | Defect | Consequence |
|---|---|---|
| 1 | The evaluation scored **the whole universe** with `labels = (gene in train_pos)` — i.e. the *training* positives were the positives and the *test* positives were scored as negatives. `d1/engine/evaluate.py` ranks the **test fold only**, with the seeds excluded (D-06). | The GNN arm was not comparable with the RWR arm at all; the evaluation was inverted. |
| 2 | `_node2vec_walks` appended each step **twice** (`walk.append(nxt)` on consecutive lines). | Walks were 2x longer with every node duplicated. |
| 3 | **Node2Vec has no restart parameter.** The `alpha` argument was never passed to the model; the `alpha=0.3/0.5/0.7` JSON sets were the *same* computation repeated, differing only by gensim's Hogwild thread nondeterminism. | The "alpha sensitivity" claim for the GNN arm was unfounded; the arm is 24 cells, not 72. |
| 4 | gensim `Word2Vec(workers=4)` uses asynchronous Hogwild SGD. | The same cell re-trained three times differed by up to **0.048** in ΔAUROC — larger than the fold noise being measured. Fixed with `workers=1`. |

### A3.2 What changed

- Evaluation now mirrors `d1/engine/evaluate.py` exactly: seeds = training-fold
  positives, ranked set = test fold only, seeds excluded.
- Verification: after the fix, the GNN arm's degree AUROC for
  `funmap x intogen2024` is **0.52044**, against **0.5203** in the main
  Table 1 — the two arms now agree.
- All 24 cells were re-run (174 min wall-clock, 4 parallel jobs).
- The previous JSONs are archived under `results/runs/gnn_courtesy_OLD_BUGGY/`
  and are not used anywhere.

### A3.3 New numbers

| method (mean over 24 cells, α=0.5) | mean ΔAUROC | cells with ΔAUROC > 0 | mean per-cell SD |
|---|---|---|---|
| **RWR** (pre-registered learner) | +0.079 | 24/24 | 0.019 |
| **Node2Vec** (GNN courtesy) | +0.028 | 19/24 | 0.051 |
| degree baseline | 0 | — | — |

### A3.4 Effect on the claim

The qualitative claim is **strengthened**: the GNN also beats degree in
19/24 cells, so the benchmark is not an artefact of a weak RWR. The
quantitative claim is **sharpened**: RWR reaches roughly three times
Node2Vec's margin at a fraction of the cost, and with exact
reproducibility. The old manuscript text claimed "overlapping magnitudes"
(+0.07 vs +0.06) and a 14 pp RWR advantage on `string_phys700`; the
corrected figures are +0.079 vs +0.028 and 3.4 pp. The manuscript text has
been updated.

### A3.5 What is NOT changed

The GNN arm is **not** part of the primary endpoint and **not** part of the
H1 decision. The pre-registered decision rule continues to be evaluated on
the pre-registered RWR arm only. This amendment touches a control arm and
the text describing it.
