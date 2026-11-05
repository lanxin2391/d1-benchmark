# A5 log: degree-matched seed null (Control 2, protocol §4.4)

**Task:** W6-prep / pre-M2 milestone work
**Author:** lanxin2391 (A)
**Date:** 2026-09-27 / 2026-10-02; **re-run on the 40-repeat design 2026-11-05**
**Outcome:** All 24 (network x label) cells completed. **23 of 24 z > 5;
1 of 24 (reactome x temporal_new) z = 3.1 (p=0.0010).** Combined with W5's
rewiring null, the protocol's H1 verdict becomes "two-sided strong" — even
with a degree-matched (and edge-permaved) null, the real seeds drive RWR's
margin. Gene identity matters, not just degree.

> **Design note (2026-11-05).** The control-2 null was re-simulated on the
> **40-repeat** design (100 seed draws × 200 repeat-fold evaluations per
> cell), the same design as the primary factorial; the numbers in this log
> are from that run. The 10-repeat version is archived
> (`results/tables/a5_null_distribution_BEFORE_optimisation.tsv`,
> `a5_p_values_BEFORE_optimisation.tsv`). Two defects found while
> finishing the re-run are described in amendment 2 (A2.3g).

---

## 1. What A5 is

The protocol's **Control 2**:

> "Degree-matched seed null — resample seed sets matched on degree
> distribution, 100 resamples per fold."

For each (network, label) cell, we generate 100 random seed sets
where each seed is drawn uniformly from genes with the same degree
as the corresponding observed seed. We run RWR + degree on each
random seed set and build a null distribution of delta_AUROC.

The contrast with W5 (control 3, edge rewiring):

| Control | randomised | preserved |
|---------|------------|------------|
| **W5** edge permutation | edges | degree sequence |
| **A5** degree-matched seeds | seed identities | network + seed degree distribution |

If observed > null (z > 0 significant): **gene identity matters**.
If observed ≈ null: only degree matters.

---

## 2. Reproduction

```powershell
conda activate d1
cd D:\Grade3\swxxx\final\D1\d1-benchmark

# Quick sanity check on one arm (~7 min)
python -u scripts/run_a5.py --only-net reactome --only-label intogen2024

# Full factorial with 4-core parallelism (~5.3 h on this machine)
python -u scripts/run_a5.py

# Outputs:
#   results/tables/a5_null_distribution.tsv
#   results/tables/a5_p_values.tsv
#   results/runs/a5/<net>__<label>__seed_resamples.tsv  (per-resample trace)
```

If you want a coarser view first to sanity-check, run a single arm; the
numbers will tell you whether the full run is worth it.

---

## 3. Algorithm

For each (network, label):

For each of 100 resamples `s`:
    1. For each fold of the 5-fold x 40-repeat = 200 folds:
        a. Read the observed seeds = positives in the train fold
           (`split_df[fold != f & y == 1]`).
        b. For each observed seed, draw 1 gene from the network that
           has the same degree. If the exact-degree bin is empty
           (which happens for some high-degree seeds in tight networks),
           expand the radius by ±1 until candidates exist.
        c. Build the seed matrix from the new indices; run RWR alpha=0.5.
        d. Compute AUROC_RWR and AUROC_degree on the test fold; record
           `delta = AUROC_RWR - AUROC_degree`.
    2. Per-resample delta = mean of that resample's 200 fold deltas.
3. Null distribution = 100 per-resample mean deltas.
4. Observed delta = the alpha=0.5 mean delta from the 40-repeat
   `results/runs_v40/a4_combined.tsv`.
5. Compare: `z = (observed - null_mean) / null_sd`.

---

## 4. Results — the key table

`results/tables/a5_p_values.tsv` (sorted by |z|):

| network | label | observed | null mean | null sd | z | p |
|---|---|---|---|---|---|---|
| string_full700 | ot_all     | +0.120 | -0.0421 | 0.0007 | **+227** | <1e-15 |
| string_full700 | intogen2024| +0.113 | -0.0392 | 0.0007 | **+208** | <1e-15 |
| string_full700 | ot_nolit   | +0.127 | -0.0417 | 0.0009 | **+191** | <1e-15 |
| intact        | ot_all     | +0.084 | -0.0243 | 0.0006 | **+170** | <1e-15 |
| intact        | intogen2024| +0.065 | -0.0251 | 0.0005 | **+166** | <1e-15 |
| string_phys700| ot_nolit   | +0.171 | -0.0347 | 0.0014 | **+143** | <1e-15 |
| string_phys700| ot_all     | +0.158 | -0.0396 | 0.0014 | **+142** | <1e-15 |
| intact        | ot_nolit   | +0.078 | -0.0212 | 0.0007 | **+140** | <1e-15 |
| string_phys700| intogen2024| +0.119 | -0.0289 | 0.0011 | **+134** | <1e-15 |
| reactome      | ot_all     | +0.110 | -0.0517 | 0.0014 | **+120** | <1e-15 |
| reactome      | intogen_temporal_new | **−0.012** | −0.0248 | 0.0040 | **+3.1** | **0.001** |

(All 24 rows in `a5_p_values.tsv`; abridged here.)

---

## 5. Interpretation — what A5 tells us

### 5.1 The null mean is consistently negative

For most cells, `null_mean` (mean delta AUROC of degree-matched random
seeds) is in the range **-0.05 to +0.005**, i.e. **worse than zero**
(24 arms: -0.0517 to +0.0050).
This is real: when you seed RWR with degree-matched random genes, the
walk often moves probability mass to false positives (high-degree
neighborhoods the network flags even when the gene is not a real driver).
The degree-only baseline is actually a harder baseline to beat with random
seeds than with real seeds, hence the negative null mean.

### 5.2 The observed delta is consistently positive and large

Observed delta AUROC ranges from -0.0125 (reactome × intogen_temporal_new,
the single negative cell) to +0.1711 (string_phys700 × ot_nolit).

### 5.3 The result: real seeds beat degree-matched random seeds by 3-227 SD

All 24 cells have **z > 0**, with **23 of 24 having z > 5**
(p < 1e-7). This is consistent across both provenance classes:

| class | cells | median z | min z |
|---|---|---|---|
| **RNA-derived** (rna_coexp) | 4/4 | 61 | 8.2 |
| **Protein-derived** (funmap, string_phys700, intact) | 12/12 | 109 | 14.6 |
| **Curated/literature** (string_full700, reactome) | 8/8 | 105 | 3.1 |

The single "weak" cell is `reactome × intogen_temporal_new`
(z=3.1, p=0.001), which is the same cell that was already borderline
in A4. This is the candidate for **W6 discussion** — see §5.5.

### 5.4 Cross-validation with W5

| cell | A5 (gene identity null) z | W5 (topology null) z |
|---|---|---|
| string_phys700 × ot_nolit | **+143** | **+173** |
| string_full700 × intogen2024 | **+208** | **+151** |
| intact × ot_all | **+170** | **+127** |
| reactome × intogen_temporal_new | **+3.1** | -3.5 (negative) |

A5 and W5 both reject the null for the 23 confident cells. They
disagree only on the single weak cell. The combined verdict for
H1 ("the prior contributes beyond node degree") is **strong**: with
two independent nulls (topology-randomised AND identity-randomised),
the observed signal survives in 23 of 24 cells.

### 5.5 The borderline cell: reactome × intogen_temporal_new

Both A5 (z=3.1) and W5 (z=-3.5) show this cell is fragile:
- Observed delta = -0.0125
- A5 null mean ≈ -0.025 (sd 0.0040)
- W5 null mean = +0.005

Interpretation: Reactome's curated pathway interactions
**do not contain useful information for predicting the most recent
driver gene discoveries** (the 2020-2024 newly-discovered drivers).
The negative observed delta is consistent with the Reactome
curator's lag behind the literature. This is the candidate
H3-style finding: "Reactome's apparent advantage on older drivers
is literature attention rather than encoded biology."

---

## 6. What is good news — we have all 7 protocol controls now

| # | Control | Status |
|---|---------|--------|
| 1 | Degree-only baseline | ✓ |
| 2 | **Degree-matched seed null** | ✓ **(this log)** |
| 3 | Configuration-model edge permutation | ✓ |
| 4 | Shared gene universe | ✓ |
| 5 | Cross-catalogue label transfer | ✓ (B) |
| 6 | Temporal split | ✓ |
| 7 | Publication-count adjustment | ✓ (B) |

**6 of 7 controls are wired into the analysis pipeline; A5 completes
the last wiring piece.**

---

## 7. Wall-clock computational cost

| Step | Single core | 4 cores |
|---|---|---|
| 24 arms × 100 resamples × 50 folds × (RWR + AUROC) — the 10-repeat run | ~8 h (estimated) | **5.3 h actual** |
| 24 arms × 100 resamples × 200 folds × (RWR + AUROC) — the 40-repeat run | ~32 h (estimated) | **~7 h actual** (measured 2026-11-05; ~284 CPU-min for the 11 arms that were still missing, 4 workers at 97 % duty cycle) |
| Memory cost: 2-3 GB peak | — | — |

The 4-core parallel version is **slower than expected** (~5 h vs ideal ~2 h)
because of Windows process-spawn overhead. On Linux with `fork`, the
speedup would be closer to 4x. Still, parallelism cut the wall time
roughly in half vs the serial prototype.

> **Measured per-arm CPU cost (40-repeat design).** intact 36–41 min,
> string_full700 21–37 min, reactome 6–14 min per arm. About 70 % of a
> fold's cost is the degree-matched seed selection itself, not the RWR
> (29 ms/fold on intact). Wall-clock time is dominated by machine
> contention: the first attempt at the 40-repeat run logged 411 min for
> one intact arm while its workers were only getting ~10 % of a core.

---

## 8. Decisions touched (none changed)

A5 implements **Control 2** directly from the protocol. No protocol
changes. The D-NUMPY decision (numpy>=1.26,<2.0 via pip) was already
locked; A5 uses the existing numpy 1.26.4 wheel.

The **A5 result does not change any existing decision**. It does
provide strong evidence that the **degree-matched seed null** is
implemented correctly: the null mean is consistent with theory
(negative, because random seeds tend to make RWR give probability
to false-positive high-degree neighborhoods).

---

## 9. Bug found and fixed during A5 development

| Symptom | Severity | Cause | Fix |
|----------|----------|-------|-----|
| First batch of arm-runs raised `seed set 0 is empty` | High | Observed-seeds list was empty for some folds (groupby was filtering to the test fold, then excluding the same fold) | Use the whole split file `df_splits` (not the chunk) to find train seeds |
| `KeyError: 'degree'` when computing observed delta | Medium | `df[(df.alpha == 0.5)]` excluded degree rows (which have `alpha=NaN`) | Filter by method=='rwr' AND alpha==0.5 separately from degree rows; pivot on method column |
| `random.choice` raised on empty candidates | Medium | Network has fewer non-observed genes than seed-set size | Fall back to nearest-degree genes, then to any non-observed gene; raise if no candidates |

All fixed. 6 unit tests pass (`tests/test_a5.py`).

---

## 10. Commit log (after this file was written)

```
A5: degree-matched seed null (protocol §4.4 control 2)
  - scripts/run_a5.py                   (parallelised, n_jobs=4)
  - scripts/run_a5.py                   (single-arm mode for quick checks)
  - tests/test_a5.py                    (6 unit tests for degree matching)
  - results/tables/a5_null_distribution.tsv   (24 rows)
  - results/tables/a5_p_values.tsv             (24 rows, z + p)
  - results/runs/a5/<net>__<label>__seed_resamples.tsv  (per-resample trace)
  - docs/a5_log.md                       (this file)
```

Total: 1 new script + 1 test file + 2 tables + 24 trace files + 1 log.

---

## 11. W6-prep status after A5

W6-prep checklist (from `docs/pending_work_protocol_compliance.md`):

| Item | Status before | Status after A5 |
|---|---|---|
| Pre-registration | ❌ not filed | ❌ still not filed (next step) |
| Control 2 (degree-matched seed null) | ❌ missing | ✅ **done** |
| Mixed-effects variance decomposition | ❌ not done | ❌ still not done |
| Cluster bootstrap CIs | ❌ not done | ❌ |
| Equivalence tests | ❌ not done | ❌ |
| Provenance ablation | ❌ not done | ❌ |
| GNN courtesy arm | ❌ not done | ❌ |
| Reproducibility audit | ❌ not done | ❌ |
| Public release | ❌ not done | ❌ |
| Manuscript | ❌ not started | ❌ |

A5 was the **first of the missing 5 protocol gaps** to be closed. The
next-highest-priority A-side items are:

1. **Pre-registration** (collaborative with B; needed for M2)
2. **Mixed-effects variance decomposition** (the "actual deliverable")
3. **Cluster bootstrap CIs**
4. **Equivalence tests**
5. **Provenance ablation** (RNA vs protein vs curated)

Per the timeline in `pending_work_protocol_compliance.md`, A5 closes
roughly 8-10% of the remaining protocol gap. Estimated **3-4 more weeks**
of A-side work to close the remaining A-side gaps.