# A1 log: understand the evaluation harness

**Task:** Day 3 (protocol section 6). Person A (Networks & Engine).
**Author:** lanxin2391
**Date:** 2026-09-26
**Outcome:** Toy report at `docs/toy_results.md`. 27 tests still pass.

This log is the reproducer for A1. It walks through what I read, what numbers
came out of the toy benchmark, what they mean, and where the artefacts live.

---

## 1. What A1 is

A0 locked the *learner*. A1 verifies the *evaluation harness*: does
`run_arm()` correctly run RWR + degree on every fold, compute AUROC, and
report the margin? The harness is what produces every number in the W5
paper, so it must be trusted before any real-data run.

The toy benchmark is the **calibration**: three synthetic networks with
known answers. If `run_arm` cannot tell them apart, the harness is broken.
If it can, we trust it for real data.

## 2. Step-by-step reproduction

```powershell
conda activate d1
cd D:\Grade3\swxxx\final\D1\d1-benchmark

# 1. Read the harness (the comments are most of the explanation)
code d1\engine\evaluate.py
code d1\networks\io.py
code d1\toy.py
code d1\engine\baselines.py

# 2. Run the toy benchmark (matches Appendix B of D1_W0_protocol_PersonA.docx)
python scripts/run_toy_benchmark.py
# Expected last line: Done in ~3-9 s. Results in results/toy/

# 3. Verify numbers
Get-Content results\toy\toy_summary.tsv

# 4. Verify the harness unit tests still pass (these test the same logic)
pytest tests\test_evaluate.py -v
# Expected: 11 passed

# 5. Verify the engine and IO tests too (nothing should regress)
pytest
# Expected: 27 passed
```

## 3. What I read (key functions)

| File | Function | Role |
|------|----------|------|
| `d1/engine/evaluate.py` | `precision_at_k(scores, y, genes, k)` | Tie-break by gene symbol (D-07): when two scores tie, alphabetic order resolves |
| `d1/engine/evaluate.py` | `evaluate(scores, y, genes)` | One fold: AUROC + AUPRC + precision-at-k; returns NaN if only one class is in y |
| `d1/engine/evaluate.py` | `check_splits(df)` | Enforces contract C4: columns, gene uniqueness within repeat, fold has positives |
| `d1/engine/evaluate.py` | `run_arm(network, splits, ...)` | **The harness**. Loads network, propagates all 50 seed vectors in one RWR call per alpha, returns C5-format rows |
| `d1/engine/evaluate.py` | `margins(results)` | Pairs RWR rows with their degree rows, computes ΔAUROC and ΔAUPRC per fold |
| `d1/engine/evaluate.py` | `summarize(results)` | Groups by (network, label_set, universe, alpha) and reports mean and SD |
| `d1/networks/io.py`   | `finalize_edges(df, net_id, ...)` | Clean: drop unmapped → drop self-loops → order alphabetically → dedupe → LCC only → write TSV + meta.json |
| `d1/networks/io.py`   | `load_network(path)` | Read C2 file → (genes, A) tuple. Weights ignored (D-01) |
| `d1/networks/io.py`   | `check_network_file(path)` | Assert C2: header, no NaN, no self-loops, alphabetical, no duplicates, single component |
| `d1/toy.py`           | `toy_network(scenario, ...)` | Build one of the three synthetic networks; return edges + positive set |
| `d1/toy.py`           | `_plant_module(G, pos, n_swaps, rng)` | Degree-preserving swaps to create a planted module around the positives |
| `d1/engine/baselines.py` | `degree_scores(A)` | Number of distinct neighbours per gene (control 1) |

## 4. Toy benchmark numbers and where they came from

```
3 networks × 10 repeats × 5 folds × (3 RWR alphas + 1 degree) = 600 rows in results/toy/toy_runs.tsv
```

Cross-check of the three scenarios against the unit tests in
`tests/test_evaluate.py`:

| Scenario | Toy numbers (ΔAUROC at α=0.5) | Unit test that locks the same logic |
|---|---|---|
| `toy_random`  | −0.012 (no signal) | `test_random_labels_no_signal` |
| `toy_hubs`    | −0.120 (degree wins) | `test_hub_labels_degree_wins` |
| `toy_module`  | +0.292 (RWR wins) | `test_planted_module_rwr_beats_degree` |

All three pass. See `docs/toy_results.md` for the full interpretation.

## 5. Reproducer artefacts

| Path | What it is |
|---|---|
| `scratch/toy_bench.log` | Verbatim output of `python scripts/run_toy_benchmark.py` |
| `results/toy/toy_runs.tsv` | 600 per-fold rows, full C5 schema |
| `results/toy/toy_summary.tsv` | 9 rows: 3 networks × 3 alphas, with mean AUROC + SD + delta |
| `results/toy/networks/{toy_random,toy_hubs,toy_module}.tsv` | Clean C2 edge lists |
| `results/toy/networks/{...}.meta.json` | C2b metadata |
| `results/toy/splits/toy__native_{...}.tsv` | Toy C4 split tables (Person B will write the real ones) |

## 6. Where this fits in the larger project

The harness has three moving parts:

```
                +-------------------+
                | d1.engine.rwr     |  (locked at engine-v1, A0)
                +---------+---------+
                          |
                          | uses
                +---------v---------+        +---------------------+
                | d1.networks.io   |<------>| data/raw + scripts/ |
                |  finalize_edges  |       | (A2-W3 builds the   |
                |  load_network    |       |  six networks)      |
                +---------+--------+       +---------------------+
                          |
                          | produces (genes, A)
                +---------v---------+
                | d1.engine.evaluate|  (this task: A1)
                |  run_arm          |
                |  margins / summary|
                +---------+---------+
                          |
                          | produces results/runs/{run_id}.tsv
                +---------v---------+
                | d1 (figures,      |  (W5-W6: real-data arms)
                |  paper text)      |
                +-------------------+
```

A1 verifies the bottom-left box. A2-W3 will populate the data side;
W4-W5 will run real arms using the same `run_arm`. Everything from W4
on is plugging numbers into `run_arm`; the harness will not change.

## 7. Day-3 protocol checklist

| Item | Status |
|---|---|
| 6.1 Read `evaluate.py` and `io.py` | ✓ |
| 6.2 Toy benchmark numbers match Appendix B | ✓ |
| 6.3 Open `results/toy/toy_runs.tsv` in Excel / pandas | ✓ (600 rows: 150 degree + 450 RWR; 30 positives per test fold) |
| 6.4 Write `docs/toy_results.md` with table + 3 sentences + 1 sentence on margin | ✓ |
| 6.5 Optional exercises (n_swaps, shared universe test) | skipped — the W0 protocol calls them optional |
| Commit + merge + push to a-networks | pending — see commit log below |
| Update `docs/sync_notes.md` | pending — to be written when S0 happens |

## 8. Commit log (after this file was written)

```
A1: toy benchmark verified and interpreted
  - docs/toy_results.md         (new)
  - docs/a1_log.md              (this file)
```

Will be merged into `main` and pushed once written. The engine-v1 tag is
unchanged (A1 does not touch the engine).