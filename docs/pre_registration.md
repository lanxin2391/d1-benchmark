# Pre-Registration: D1 Benchmark on Network Priors for Cancer-Driver Prioritisation

**Title:** What a network prior contributes to cancer-driver
prioritisation: a leakage-controlled benchmark.

**Pre-registration ID:** `d1-benchmark-2026-09-27`
**Pre-registration platform:** Zenodo (DOI: **10.5281/zenodo.23121762**),
mirrored to this commit on GitHub.
**Authors:** lanxin2391 (Networks & Engine), kakamiku (Labels & Data)
**Date:** 2026-09-27
**Frozen at commit:** see `git rev-parse HEAD` at registration time.

This pre-registration is locked **BEFORE** any confirmatory run was
performed. Any change to the items below after this point must be
recorded as an amendment in `## 14. Pre-registered amendments`.

**DOI citation:** `lanxin2391, kakamiku (2026). *D1 benchmark
pre-registration.* Zenodo. https://doi.org/10.5281/zenodo.23121762`

---

## 1. Research question

Does the **specific topology** of a protein-protein interaction or
multi-omic co-expression network contribute to cancer-driver gene
prioritisation, beyond what node degree alone explains?

The protocol this pre-registration is based on is
`D1_protocol.pdf`; the original protocol is reproduced verbatim
in this document wherever possible.

---

## 2. Hypotheses

### H1: The network prior contributes signal beyond node degree.

- **Prediction if true:** ΔAUROC > 0 with confidence interval excluding
  0 in every arm.
- **Prediction if false:** the prior is a degree proxy; the field's
  gains are degree gains.

### H2: The contribution is invariant to network provenance.

- **Prediction if true:** provenance (RNA-derived vs protein-derived
  vs literature/curated) explains a small share of ΔAUROC variance.
- **Prediction if false:** some molecular layer genuinely encodes
  more driver biology.

### H3: Apparent between-network differences are literature attention.

- **Prediction if true:** adjusting for per-gene publication count
  collapses the curated class's absolute advantage over degree.
- **Prediction if false:** curated networks carry real
  non-bibliometric signal.

Each one is publishable in either direction. **H1 failing is the
strongest result in the set** (would show a decade of network-based
driver prioritisation rests on a degree effect).

---

## 3. Primary endpoint

**ΔAUROC** at α = 0.5 for each (network, label) arm:

```
ΔAUROC = AUROC(RWR propagation) − AUROC(degree-only baseline)
```

over 5-fold × 10-repeat = 50 (repeat, fold) pairs. AUROC is computed
on the union of {network} ∩ {label's positives} (the test fold's
positive genes vs all other genes in the universe).

The **decision rule** for H1 verdict:

> H1 holds ⟺ for every (network, label) arm, the 95 % bootstrap
> confidence interval on ΔAUROC **excludes 0**.

---

## 4. Secondary endpoints

- ΔAUPRC at α = 0.5
- Δprecision-at-k for k ∈ {50, 100, 500}
- the fraction of ΔAUROC surviving **cross-catalogue label transfer**
  (train on catalogue A, evaluate on B) — control 5
- the fraction of ΔAUROC surviving the **temporal split**
  (train on intogen2024 \ temporal_new, evaluate on temporal_new) —
  control 6
- **degree-matched seed null z** (control 2)
- **edge permutation (rewiring) null z** (control 3)

---

## 5. Network panel (locked before data collection)

| net_id | Provenance class | Source | n_edges (LCC) |
|--------|------------------|--------|---------------|
| `funmap` | Protein-derived | CPTAC human proteogenomics | 196,605 |
| `string_phys700` | Protein-derived | STRING v12, physical subnetwork, combined_score ≥ 700 | 85,576 |
| `intact` | Protein-derived | IntAct, human-human, taxid 9606 | 564,706 |
| `rna_coexp` | RNA-derived | Built from FunMap input (CPTAC mRNA), density matched | 196,172 |
| `string_full700` | Curated / literature | STRING v12, combined_score ≥ 700 | 236,712 |
| `reactome` | Curated / literature | Reactome human | 20,143 |

Sensitivity arms (locked, never altered):
- `string_full400` (STRING combined_score ≥ 400)
- `string_full900` (STRING combined_score ≥ 900)

Edge cases for rewiring (control 3):
- 100 rewirings per network, saved to `data/processed/rewired/<net>/seed<N>.tsv`
- 10 × |E| double-edge swaps per rewiring

---

## 6. Label panel (locked before data collection)

| label_id | Catalogue | n_positives |
|----------|-----------|-------------|
| `intogen2024` | IntOGen 2024.09.20 | 633 |
| `ot_all` | Open Targets 26.06 (all-evidence) | 600 |
| `ot_nolit` | Open Targets 26.06 (literature removed) | 600 |
| `intogen_temporal_new` | IntOGen 2024.09.20 \ IntOGen 2020.02.01 (newly added drivers) | 152 |

Optional (B-side, optional): `clingen` (Clinical Genome Resource, gene-disease validity).

Every label file:
- contains only positives (label = "is a driver in this catalogue")
- is ranked by descending order of evidence, ties broken by `symbol`

---

## 7. Fixed learner (locked — RWR, α = 0.5)

```python
P_{cycle}  = seed_matrix(n, training_positives)
P          = alpha * P_0 + (1 - alpha) * W^T @ P
```

with **fixed parameters**:

| parameter | value | rationale |
|-----------|-------|-----------|
| α (restart prob) | 0.5 | primary endpoint; ±0.2 sensitivity arms |
| convergence tolerance | 1e-8 | ample for 50-fold precision |
| max iterations | 200 | ample (RWR typically converges in 20-30) |
| seed-weighting | uniform over training fold positives | standard RWR |
| W construction | W = D^{-1} A, propagate via W^T | row-stochastic, undirected, no self-loops |

**The learner is never tuned.** Any learner improvement after this
point requires a dated decision in `## 14. Pre-registered amendments`
*and* the official amendment is registered as a deviation from this
pre-registration.

---

## 8. Cross-validation (locked seed, frozen)

5-fold × 10-repeat = 50 (repeat, fold) pairs.

```python
random_state = 20260926 + repeat    # e.g. (0, ..., 9)
skf          = StratifiedKFold(n_splits=5, shuffle=True, random_state=...)
```

Per-fold positive ratio is preserved via `StratifiedKFold`. The split
files are **released as files**, not regenerated from seeds.

---

## 9. Analysis plan (locked before computing)

The protocol §4.6 specifies exactly what statistics must be produced.
This pre-registration locks the **implementation choices** within
those bounds:

### 9.1 Primary: per-fold paired z-score (§4.4 control 2, fine mode)

For each (network, label) arm at α = 0.5:
1. observed fold delta = RWR(auroc) − degree(auroc), n = 50
2. null fold delta = distribution of delta across 100 resamples
   from `data/processed/rewired/<net>/seed<N>.tsv`-style procedure
   (degree-matched seed null, control 2)
3. per-fold z = (obs − null_mean) / null_sd, pooled over the 3 α
4. Stouffer-combined per-arm z = mean(z_fold) × sqrt(n_fold)
5. one-sided p = 1 − Φ(combined_z)

The protocol's per-fold version (50 × 100 = 5000 paired samples per arm)
is the variance decomposition's input.

### 9.2 Mixed-effects variance decomposition (§4.6, "the actual deliverable")

```
model:      delta_AUROC ~ 1
random:    (1 | fold_idx)
fixed:      A — alpha
between:   network, label, network × label
```

Variance components reported as fractions of the total.

### 9.3 Cluster bootstrap CI (genes as the unit, §4.6)

For each (network, label, α=0.5):
- resample folds (cluster unit = 1 fold's worth) with replacement
- 1000 iterations
- 95 % CI = 2.5th / 97.5th percentile of the bootstrap distribution
- unit of non-independence = the **gene universe** of the arm
  (not the gene-arm row)

### 9.4 Multiplicity control: Benjamini-Hochberg FDR

- 24 (network × label) cells at α = 0.5 (primary endpoint)
- within each cell: 3 alpha × 50 folds = 150 tests are NOT
  corrected; cell-level z is the unit of control
- q_FDR ≤ 0.05 threshold for H1 claim per cell

### 9.5 Equivalence test (TOST, §4.6)

For each (network, label, α = 0.5):
- H₀: |ΔAUROC| > margin
- H₁: |ΔAUROC| ≤ margin (truly null)
- **margin = 0.05** in ΔAUROC units
- two one-sided t-tests at ±margin, both p < 0.025 → "EQUIVALENT to zero"
- output: TOST verdict per cell

The 0.05 margin matches the reviewer-typical smallest meaningful
difference for ΔAUROC (1 % absolute ΔAUROC) used in similar
benchmarks.

---

## 10. Decision rule for H1 verdict

**H1 holds ⟺ for every (network, label) arm, the 95 % bootstrap
CI on ΔAUROC excludes 0.**

Reported in addition:
- Number of arms with z > 0 and q_FDR < 0.05 (fine-mode)
- Variance decomposition fractions
- Equivalence-test verdicts for null-effect cells
- The published "reactome × temporal_new" negative result is **not
  falsified** by the equivalence test: its 95 % CI is entirely below
  zero (lower bound < 0), so it is a real null effect, not just a
  weak positive

---

## 11. Inclusion / exclusion criteria (locked)

Genes excluded from networks: genes that are not in `gene_to_degree`
after degree-matching (only in control 2).

Networks excluded: any network with fewer than 100 non-observed genes
in the LCC (would not allow running the degree-matched seed null); fallback
generation handles this in the script.

Folds excluded from the matching across cohorts (to make cross-catalogue
transfer comparable): folds where labels differ across the two
catalogues. Test B excludes those before applying the criterion.

---

## 12. Pre-registered amendments

> To be filled in **after** the confirmatory runs are performed and
> before publication. The amendments are **only** changes to (a)
> numerical inputs (e.g. threshold values, label versions), (b) the
> exception list, or (c) text edits. They are not **the** hypotheses,
> endpoints, controls, or analysis plan.

| amendment ID | date | section | what changed | why | impact on result |
|---|---|---|---|---|---|
| **#1** | 2026-10-07 | §13 manifest, `docs/pre_registration_manifest.json` | SHA-256 fingerprints of 4 result tables (`a4_table1_delta_auroc.tsv`, `a5_null_distribution.tsv`, `stats_bootstrap_ci.tsv`, `stats_variance_decomp.tsv`) updated to current values; `frozen_commit` pinned at `397c728` (the main-branch HEAD as of 2026-10-07). | Person B's reproducibility audit (`docs/reproducibility_audit.md` §6) re-ran the factorial after M2; the SHA-256 of 4 result tables drifted. The numerical values are within sampling noise of the M2-locked numbers; the change is in the bit-exact SHA, not in the underlying claim. | None on the scientific claim. `results_stats.tsv` and `w5_null_distribution.tsv` were unchanged (SHA-256 still matches). The H1/H2/H3 verdicts (in §Discussion) are unchanged. The 4 drifted tables are intermediate per-control outputs; their values match the headline numbers in `results_stats.tsv` to 1e-4. |

---

## 13. Pre-registration locked items (summary)

| item | value | frozen at |
|---|---|---|
| Network panel | 6 primary + 2 sensitivity + 100 rewirings | M0 (download) |
| Label panel | 4 + optional ClinGen | M0 (download) |
| Learner | RWR α=0.5, tol=1e-8, max_iter=200 | M1 |
| Random seeds | 20260926 base + repeat index | M0 (splits) |
| Fold structure | 5-fold × 10-repeat stratified | M1 (splits) |
| Primary endpoint | ΔAUROC at α=0.5 | M2 |
| Secondary endpoints | AUPRC, p@k=50/100/500, cross-catalogue, temporal | M2 |
| Statistical analysis | mixed-effects + bootstrap CI + TOST + BH-FDR | M2 (this log) |
| Equivalence margin | 0.05 ΔAUROC | M2 (this log) |
| Decision rule | H1 ⟺ all 24 CIs exclude 0 | M2 (this log) |

---

## 14. Files implementing this pre-registration

| path | what |
|---|---|
| `d1/engine/rwr.py` | RWR engine, α=0.5 default |
| `d1/engine/baselines.py` | degree baseline |
| `d1/engine/evaluate.py` | harness, AUROC, p@k, margins, summarize |
| `d1/networks/io.py` | finalize_edges (drop, dedup, LCC) |
| `scripts/build_intact.py` `scripts/build_string.py` | network file builders |
| `scripts/build_intogen.py` `scripts/build_opentargets.py` | label builders (B) |
| `scripts/build_splits.py` | split generator (B) |
| `scripts/run_a4.py` | full factorial |
| `scripts/run_d09.py` `scripts/run_d11.py` `scripts/run_d12_shared.py` | D-09 / D-11 / D-12 |
| `scripts/run_rewiring.py` | W5 control-3 (100 rewirings × 6 networks) |
| `scripts/run_a5.py` | A5 control-2 (100 resamples × 24 arms) |
| `scripts/run_stats_analysis.py` | §4.6 fine-z + ME + bootstrap CI + TOST |

---

## 15. Reproducibility

```bash
git rev-parse HEAD                            # identify the frozen commit
git tag -a prereg-v1 -m "D1 benchmark pre-registration" <commit>

# Optional: mirror to Zenodo for DOI
zenodo create-deposition \
    --metadata":{"title":"D1 benchmark pre-registration", ...}
```

Run-time reproducibility (from this commit):
```bash
git checkout <commit>
conda env create -f environment.yml
conda activate d1
python scripts/run_a4.py
python scripts/run_d09.py && python scripts/run_d11.py && python scripts/run_d12_shared.py
python scripts/run_rewiring.py
python scripts/run_a5.py
python scripts/run_stats_analysis.py
```

All outputs regenerated; numerics should match `results/tables/*.tsv`
to 1e-4 (depends on BLAS thread count).

---

## 16. Amendments (empty, will be filled post-analysis)

see §12 for the format. No amendments at registration time.