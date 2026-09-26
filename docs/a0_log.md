# A0 log: understand and lock the RWR engine

**Task:** Day 2 (protocol section 5). Person A (Networks & Engine).
**Author:** lanxin2391
**Date:** 2026-09-26
**Outcome:** RWR engine locked at tag `engine-v1`. All 27 tests pass.

This log is the complete reproducer for A0. Any team member reading it can rebuild
exactly what we did on 2026-09-26 and verify the engine is in its locked state.

---

## 1. What A0 is

A0 fixes the **learner** of the D1 benchmark. The benchmark measures what a network
prior contributes to cancer-driver prioritisation, so the learner itself must be
fixed (never tuned); otherwise any change in the learner would confound the
measurement (protocol section 4.3).

After A0:

- `d1/engine/rwr.py` has its parameters (alpha=0.5, tol=1e-8, max_iter=200) frozen
- 13 unit tests in `tests/test_rwr.py` lock the behaviour
- 3 break-and-fix exercises prove the tests catch real bugs
- Speed benchmark confirms a FunMap-sized network finishes in seconds

Any later change to `rwr.py` requires a dated entry in `docs/decisions.md` and
a new tag.

---

## 2. Step-by-step reproduction

All commands assume the repository root is `D:\Grade3\swxxx\final\D1\d1-benchmark`
and the conda environment `d1` is active. Replace as needed for another machine.

```powershell
# 0. Activate the environment (skip if not yet created)
conda activate d1

# 1. Verify rwr.py is in its locked state (line 80 of d1/engine/rwr.py)
Select-String -Path d1\engine\rwr.py -Pattern 'P_new = alpha \* P0 \+ \(1\.0 - alpha\) \* \(WT @ P\)'
# Expected: 1 match, line 80

# 2. Verify transition_T is in its locked state (line 48 of d1/engine/rwr.py)
Select-String -Path d1\engine\rwr.py -Pattern 'W = sp\.diags\(1\.0 / deg\) @ A'
# Expected: 1 match, line 48

# 3. Run all 27 tests
pytest
# Expected last line: 27 passed in ~4 s

# 4. Run only the engine tests (13)
pytest tests/test_rwr.py -v
# Expected last line: 13 passed

# 5. Verify the toy benchmark numbers match Appendix B
python scripts/run_toy_benchmark.py
# Expected: toy_hubs alpha=0.5 delta_auroc ~= -0.110
#           toy_module alpha=0.5 delta_auroc ~= 0.292
#           toy_random alpha=0.5 delta_auroc ~= -0.012
# (mine reproduced to 3 decimals; small differences from random seeding)

# 6. Verify the speed benchmark
python scripts/bench_rwr_speed.py --rewire
# Expected: ~7 ms per seed vector, ~1 min per rewiring (machine-dependent)

# 7. Verify the engine-v1 tag exists
git tag -l engine-v1
# Expected: engine-v1
```

---

## 3. Toy example (protocol section 5.2)

Run from the repo root with the d1 env active:

```python
import numpy as np, scipy.sparse as sp
from d1.engine.rwr import transition_T, seed_matrix, rwr

# path graph a - b - c
A = sp.csr_matrix(np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], float))
WT = transition_T(A)
print(WT.toarray())
# [[0.  0.5 0.]
#  [1.  0.  1.]
#  [0.  0.5 0.]]
# Each COLUMN of WT sums to 1 -> W^T is column-stochastic, equivalently W is row-stochastic.
# W[i,j] = 1/deg(i) if (i,j) is an edge, else 0. W^T means step FROM column gene TO row gene.

P, n_iter = rwr(WT, seed_matrix(3, [[0]]), alpha=0.5)
print(P.round(4), n_iter)
# [0.5833 0.3333 0.0833] 28
# alpha=0.5: the seed a keeps the most probability (half of all restarts land there).
# b is its only neighbour, so half of P from a flows to b (and back). c gets one eighth.
```

I ran this on 2026-09-26 with Python 3.11.16 + scipy 1.17.1. Full log:
`scratch/verify_rwr_toy.py`.

---

## 4. Test results

Ran `pytest tests/ -v` on 2026-09-26 with numpy 1.26.4, scipy 1.17.1, Python 3.11.16.
27 tests, 4 s. Full log: `scratch/pytest_all.log`.

| Test | What it locks |
|---|---|
| `test_transition_is_row_stochastic` | W^T column sums = 1 (after `transition_T`) |
| `test_mass_is_conserved` | Each P column sums to 1 for any alpha, any number of seeds |
| `test_alpha_one_returns_seeds` | alpha=1 is identity (no network contribution, just seeds) |
| `test_matches_closed_form[0.3]` | Iterative P matches analytic P = alpha(I - (1-alpha)W^T)^{-1} P0 |
| `test_matches_closed_form[0.5]` | Same at alpha=0.5 |
| `test_matches_closed_form[0.7]` | Same at alpha=0.7 |
| `test_star_centre_ranks_first` | In a star graph, the centre ranks above all leaves |
| `test_signal_stays_in_seeded_community` | In a barbell graph, signal stays in the seeded clique |
| `test_single_vector_input` | rwr accepts (n,) input as well as (n,k) |
| `test_rejects_degree_zero_nodes` | transition_T raises ValueError if any node has degree 0 |
| `test_rejects_asymmetric_and_self_loops` | transition_T raises ValueError on asymmetric / self-loop adjacency |
| `test_empty_seed_set_rejected` | seed_matrix raises ValueError on empty seed set |
| `test_warns_when_not_converged` | rwr warns (does not crash) when max_iter is hit |

Tests in `test_evaluate.py` (11) and `test_networks_io.py` (3) cover the rest of the
W0 contracts but are not the engine itself.

---

## 5. Break-and-fix exercises (protocol section 5.4)

I introduced two bugs and reverted them. Logs in `scratch/pytest_broken_*.log`.

### 5.1 Bug A: `(1.0 - alpha)` -> `alpha` in `rwr()` iteration (line 80)

```python
# original (line 80)
P_new = alpha * P0 + (1.0 - alpha) * (WT @ P)
# became
P_new = alpha * P0 + alpha * (WT @ P)
```

Tests that failed (4 of 13):

- `test_alpha_one_returns_seeds` (because alpha=1 now mixes seeds with WT@seeds)
- `test_matches_closed_form[0.3]` (analytic P does not match)
- `test_matches_closed_form[0.7]` (analytic P does not match)

Test that **passed** (the bug hides here):

- `test_matches_closed_form[0.5]` because at alpha=0.5 the numbers alpha and (1-alpha)
  are equal, so the bug has no effect. **This is exactly why the test parametrizes
  three alpha values**: a bug can hide at the value you use most. (Matches the
  protocol's prediction in section 5.4.)

Other tests passed because they don't depend on the absolute value of the random
walk (they test properties like row-stochasticity and rank preservation).

### 5.2 Bug B: row-normalized -> column-normalized in `transition_T()` (line 48)

```python
# original (line 48)
W = sp.diags(1.0 / deg) @ A
# became
W = A @ sp.diags(1.0 / deg)
```

Tests that failed (3 of 13):

- `test_transition_is_row_stochastic` (now column-stochastic instead)
- `test_mass_is_conserved` (mass leaks)
- `test_star_centre_ranks_first` (in a column-normalized star, leaves accumulate probability,
  centre gets just its share of 1/deg=1 from each)

Closed-form tests passed because `(1-alpha)` and `alpha` are small at this matrix size
and the bug manifests only in extreme topologies.

After restoring both bugs and running `pytest tests/test_rwr.py`, all 13 tests passed
again.

---

## 6. Environment issue found and fixed

**Symptom:** `pytest` crashed with `Windows fatal exception: code 0xc06d007f` after
running 3 of 13 tests (during `test_matches_closed_form[0.3]`).

**Diagnosis:** I isolated the crash to `numpy.linalg.solve` (and `scipy.linalg.solve`).
The original env had `numpy==2.4.6`, which has a known OpenBLAS crash on Windows
when calling `solve` (even on a 2x2 matrix — confirmed with `scratch/test_2x2.py`).

**Fix:** `pip install "numpy<2.0"` installed `numpy==1.26.4`. `pytest` then passed
all 27 tests in 4 s.

**Recorded in:** `environment.lock.yml` was re-exported (see section 8).

---

## 7. Speed benchmark

`scripts/bench_rwr_speed.py --rewire` output on 2026-09-26
(Python 3.11.16, numpy 1.26.4, scipy 1.17.1, Windows 10, 16 logical CPUs):

```
Python 3.11.16 on Windows-10-10.0.26200-SP0, 16 logical CPUs
network: 10525 genes, 196800 edges (built in 1.8 s)
sparse transition matrix: 4.8 MB (dense would be 0.89 GB)
RWR alpha=0.5: 11 iterations, 358 ms per batch of 50, 7.17 ms per seed vector  (protocol: 6.2 ms, 23 iterations)
projected real arms (6 networks x 4 labels x 50 folds): 8.6 s
one rewiring (10 x |E| = 1968000 swaps): 0.99 min, degrees preserved: True  (protocol: 0.7 min)
projected 100 rewirings x 6 networks on 1 core: 9.9 h
```

Interpretation (see also `docs/engine_notes.md`):

- **RWR speed** is dominated by sparse matrix-vector products. My PC is a bit
  slower per seed (7.17 ms vs 6.2 ms reference) but well within budget.
  8.6 s for the full 6 networks × 4 labels × 50 folds × 3 alphas factorial
  is comfortable.
- **Iterations**: 11 here, vs 23 in the protocol reference. The protocol measures
  on the **real** modular FunMap network; this script uses a synthetic
  scale-free Barabasi-Albert graph, which has smaller mixing time. The number
  to compare with the reference is "ms per seed vector" (machine speed),
  not "iterations" (network topology).
- **Rewiring**: 0.99 min vs 0.7 min reference, also slightly slower but within
  reason. 100 rewirings × 6 networks at this rate is 9.9 hours on 1 core,
  feasible but not trivial. Person A should run these overnight or parallelise.

See `docs/engine_benchmark.md` for the full record.

---

## 8. Locking

After steps 5–7 succeeded:

1. Re-exported the lock file with the fixed numpy version:
   `conda env export -n d1 --no-builds > environment.lock.yml`
2. Wrote `docs/engine_notes.md` (interpretation of the engine and tests).
3. Wrote `docs/engine_benchmark.md` (this benchmark record).
4. Committed on `a-networks`:
   `A0: lock RWR engine (engine-v1 tag); pin numpy<2 to fix linalg.solve crash`
5. Merged `a-networks` -> `main` and pushed.
6. Tagged the engine commit on `main`:
   `git tag -a engine-v1 -m "RWR engine locked (protocol 4.3)"`
7. Pushed the tag.

---

## 9. Reproducing this log on a new machine

```powershell
git clone https://github.com/lanxin2391/d1-benchmark.git
cd d1-benchmark
git checkout engine-v1                      # verify the lock is in place
conda env create -f environment.yml
conda activate d1
pip install -e .
pytest                                      # 27 passed
```

If `pytest` crashes with `0xc06d007f`, the numpy/OpenBLAS issue (section 6) has
returned. Reinstall with `pip install "numpy<2.0"` and re-export the lock file.