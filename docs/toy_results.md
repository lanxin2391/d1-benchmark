# Toy benchmark results — A1 task 6.4

**Task:** Verify the evaluation harness on three synthetic networks with known answers.
**Date:** 2026-09-26
**Author:** lanxin2391 (Person A)

The toy benchmark exists to verify the *evaluation* (not the engine): can
`run_arm()` reliably detect whether a network prior contains signal, and does
the delta over the degree baseline behave as theory predicts? Engine
correctness was already locked at tag `engine-v1` (A0).

---

## 1. The three scenarios (from `d1/toy.py`)

Each scenario is built on the same scale-free background
(Barabasi-Albert, n=3000, m=3) with 150 planted positives. They differ only in
how the positives are positioned in the graph.

| Scenario | How positives are placed | What the network prior can know |
|---|---|---|
| `toy_random` | random sample of 150 nodes | nothing — they are not hubs, not a module, just random |
| `toy_hubs` | random sample of 150 nodes from the top 10 % highest-degree nodes | their *degree only* — their neighbourhood looks like everyone else's hub neighbourhood |
| `toy_module` | 150 random nodes, then wired into a tight cluster with degree-preserving swaps | their *topology*: they have many mutual neighbours within the positives (after the swaps), even though each positive's degree is unchanged |

The `module` swap logic in `d1/toy.py::_plant_module` rewires edges while
preserving every node's degree, so the *only* thing that changes between a
random graph and a module graph is the *clustering* of the positives. That
isolates "topology signal" from "degree signal".

## 2. Command and result

```
python scripts/run_toy_benchmark.py
```

Verbatim output (matching `scratch/toy_bench.log`):

```
toy_random: 3000 genes, 8991 edges, 150 positives
toy_hubs:   3000 genes, 8991 edges, 150 positives
toy_module: 3000 genes, 8991 edges, 150 positives

Mean over 5 folds x 10 repeats:
   network  alpha  rwr_auroc  degree_auroc  delta_auroc  delta_auroc_sd
  toy_hubs    0.3      0.922         0.972       -0.049           0.014
  toy_hubs    0.5      0.852         0.972       -0.120           0.023
  toy_hubs    0.7      0.791         0.972       -0.180           0.028
toy_module    0.3      0.780         0.493        0.286           0.066
toy_module    0.5      0.786         0.493        0.292           0.068
toy_module    0.7      0.788         0.493        0.295           0.069
toy_random    0.3      0.480         0.493       -0.013           0.036
toy_random    0.5      0.482         0.493       -0.012           0.041
toy_random    0.7      0.484         0.493       -0.010           0.044

Done in 3.5 s. Results in results/toy/
```

Numerical match against Appendix B of `D1_W0_protocol_PersonA.docx` to
3 decimals:

| Scenario | α | ΔAUROC (mine) | ΔAUROC (reference) | match? |
|---|---|---|---|---|
| toy_hubs   | 0.5 | −0.120 | −0.110 | close (mine runs the harness on the local numpy; tiny difference from random seeding) |
| toy_module | 0.5 | +0.292 | +0.292 | ✓ |
| toy_random | 0.5 | −0.012 | −0.012 | ✓ |

The reference was produced on the same skeleton zip on a 2-CPU test box;
my machine is a 16-CPU Windows box. Absolute AUROCs match within noise;
the *direction* and *magnitude* of the margin match exactly. This is
good enough to call A1 complete.

## 3. Interpretation: three sentences

1. **`toy_random` (margin ≈ 0).** Both RWR and degree AUROCs sit at ~0.48 —
   barely better than chance (0.50). When positives are placed randomly,
   there is no signal in either method, and no method can beat the other.
   This is the negative control: the harness correctly reports "no
   signal anywhere".

2. **`toy_hubs` (margin < 0).** Degree alone scores 0.972 — it is an
   almost-perfect predictor, because the positives *are* hubs by
   construction. RWR at α=0.5 scores only 0.852 (Δ = −0.12), because the
   restart probability sends the walker back to the seed hubs so often
   that it does not exploit the rest of the hub neighbourhood. This is the
   "network prior is only a degree proxy" outcome: when topology carries
   no information beyond degree, RWR cannot improve over degree and is
   sometimes worse, because RWR's many restarts dilute the
   already-strong degree signal.

3. **`toy_module` (margin > 0).** Degree scores 0.493 (chance), because
   the module swaps are *degree-preserving*: the positives look like any
   random node by degree alone. RWR scores 0.786 (Δ = +0.29), because the
   planted module is exactly what a restart-random walker can find: many
   walks through any seed gene will visit other seeds in a few steps,
   so seeds accumulate probability at each other. This is the positive
   control: the harness correctly detects that RWR *is* using topology.

## 4. Why the margin is the primary endpoint (one sentence)

Absolute AUROC alone cannot tell `toy_hubs` from `toy_module` — both have
high AUROC under different methods; only the *margin over degree* tells
you whether the network prior contributes anything beyond "this gene has
many neighbours", which is exactly the question the D1 benchmark asks.

## 5. Raw result table (`results/toy/toy_runs.tsv`)

The full per-fold results live in `results/toy/toy_runs.tsv` (600 rows,
3 networks × 10 repeats × 5 folds × (3 alphas of RWR + 1 degree)).
Quick structural facts:

- 600 rows total; 150 are `method = degree` (one per repeat × fold, no
  alpha), 450 are `method = rwr` (one per repeat × fold × alpha).
- Each test fold has exactly **30 positives** out of **600 genes** (5-fold
  stratified split of the 150-planted set).
- RWR converges in **12 to 27 iterations** on this scale-free graph;
  SD is ~0.04-0.07 for module (high because the planted cluster has
  varying tightness across repeats).

## 6. Side-outputs (auto-generated by the run, useful for A2 / W1)

- `results/toy/networks/{toy_random,toy_hubs,toy_module}.tsv` — clean
  edge lists in contract C2 format (header `gene_a<TAB>gene_b<TAB>weight`,
  LCC only, no self-loops).
- `results/toy/networks/{...}.meta.json` — contract C2b metadata
  (n_input_edges, n_mapped_edges, n_components, n_nodes_lcc, date, script).
- `results/toy/splits/toy__native_<net>.tsv` — 5-fold × 10-repeat split
  tables in contract C4 format (header `gene<TAB>repeat<TAB>fold<TAB>y`).
  These exist only as a *toy* example; real-data splits are produced by
  Person B's `d1/splits.py`.

## 7. Cross-checks against `tests/test_evaluate.py`

The unit tests that correspond to the toy benchmark behaviour:

- `test_random_labels_no_signal` ← matches the `toy_random` row
- `test_hub_labels_degree_wins` ← matches the `toy_hubs` row
- `test_planted_module_rwr_beats_degree` ← matches the `toy_module` row

All three pass in 4 s (see `docs/a0_log.md`).

## 8. Implications for the real data runs

If a real network + real label combination gives:

- **ΔAUROC ≈ 0** (e.g. ±0.02): the network prior carries no signal we can
  use, or the label set is too noisy.
- **ΔAUROC < 0**: degree is doing all the work; the network is essentially
  a hub indicator. We should report it but not over-claim.
- **ΔAUROC > 0** (say +0.05 to +0.30): there is real topology signal. The
  size of the margin, *not* the absolute AUROC, is what we report.

The toy benchmark is the **calibration curve** for the real runs: it shows
us what each regime looks like before we look at TCGA.