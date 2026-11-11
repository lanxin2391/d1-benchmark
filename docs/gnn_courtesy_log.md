# GNN courtesy arm log (protocol §7)

**Task:** §7 control; rebuttal of the "RWR is too weak" reviewer objection.
**Author:** lanxin2391 (A)
**Date:** 2026-10-03 / 2026-10-04; **recomputed 2026-10-29**;
**this log rewritten 2026-11-11** against the released table.
**Status:** ✅ 24/24 (network, label) cells complete, on the
40-repeat design.

> **Why this log was rewritten (2026-11-11).** The previous version of
> this file reported the PRE-FIX numbers of the arm. Amendment 3 fixed
> four defects (the evaluation ranked the whole universe against the
> *training* positives instead of the test fold; the random walk
> appended each step twice; the alpha dimension did not exist for
> node2vec, so the three "alpha" arms were the same computation
> repeated; gensim's Hogwild threading made the learner
> non-reproducible). The numbers below are recomputed from
> `results/tables/gnn_courtesy_arm.tsv`, which is the authoritative
> released table. The old JSONs are archived under
> `results/runs/gnn_courtesy_OLD_BUGGY/` and are not used.

---

## 1. Method, and what this arm actually is

A **node2vec-style unsupervised embedding** (Grover & Leskovec 2016),
implemented in pure Python + gensim Word2Vec:

| hyperparameter | value |
|---|---|
| dimensions      | 64 |
| walk_length      | 40 |
| num_walks        | 80 per node |
| window           | 10 |
| p, q (return, in-out) | 1.0, 1.0 |
| min_count       | 1 |
| epochs           | 1 |
| seed             = 0 |
| alpha            | **not applicable** (`alpha_applies: false` in the released table) |

Per fold (5-fold × 40-repeat = 200 per cell):
1. take the train-fold positives;
2. compute the mean of their embeddings;
3. score every gene by its cosine similarity to that mean;
4. compute AUROC on the test fold's positives vs all other genes;
5. delta = embedding_AUROC − degree_AUROC.

This mirrors the existing A4 / D-09 / W5 / A5 pipeline, so the same
**ΔAUROC vs degree** axis is comparable.

### The naming caveat (stated because the name would otherwise mislead)

The walk uses **p = q = 1.0**, i.e. the uniform first-order walk.
Per the appendix of Hamilton et al. (2017, *NIPS*, Appendix C), node2vec
with p = q = 1 is **equivalent to DeepWalk** (Perozzi et al. 2014,
*KDD*): the second-order bias is inert, and the only difference from
DeepWalk is the embedding dimension and the corpus sampler. This arm is
therefore a **DeepWalk-style unsupervised embedding**, and the paper
names it that way rather than claiming a node2vec implementation it did
not make.

Two consequences are stated plainly:

1. This is a **no-feature** embedding. The graph neural networks that
   *do* use node features (GCN, GAT, GraphSAGE) are **not** implemented
   here and are **not** represented by this arm. Their published
   behaviour on featureless graphs is untested by this study.
2. The node2vec original (Grover & Leskovec 2016) is **not in the local
   reference library** this study was built from; the equivalence above
   rests on the GraphSAGE appendix and on the descriptions of node2vec
   in the other papers of that corpus, not on the original text. This is
   a limitation of the citation base, not of the measurement.

---

## 2. Results — 24 (network, label) cells at α = 0.5

The script caches each (network, label) cell at
`results/runs/gnn_courtesy/{net}__{label}__alpha0.5.json`, so re-runs
are idempotent. The aggregate
`results/tables/gnn_courtesy_arm.tsv` reports the α = 0.5 results, one
row per (network, label).

**24/24 cells complete.** (The "72 cells" of the pre-fix log counted
the three alpha arms, but node2vec has no alpha parameter: the three
"alpha" arms were the same computation repeated. That defect was fixed
in amendment 3, and the arm is a 24-cell design.)

---

## 3. Comparison with A4 (RWR) — the §7 rebuttal table

The single table that goes into the discussion section (mean over the
same 24 cells, α = 0.5; source `results/tables/gnn_courtesy_arm.tsv`
and `results/tables/a4_table1_delta_auroc_v40.tsv`):

| method (at α=0.5) | mean ΔAUROC | positive cells |
|---|---|---|
| **RWR** (the benchmark's chosen learner) | **+0.079** | 24/24 |
| **embedding baseline** (node2vec p=q=1) | **+0.028** | 19/24 |
| **degree baseline** (no propagation) | 0 | — |

→ Both non-trivial learners beat the degree baseline in most cells, so
the benchmark's claim "the prior contributes beyond degree" holds across
learner families. The **more expressive** learner does *not* win: the
embedding reaches roughly a third of RWR's margin, and its per-cell
estimate is far less precise (mean SD of ΔAUROC 0.051 against RWR's
0.019).

### 3.1 Per-network comparison (α = 0.5, mean over the 4 labels)

| network            | embedding ΔAUROC | RWR ΔAUROC |
|---|---|---|
| string_phys700     | +0.019 | **+0.122** |
| string_full700     | +0.057 | **+0.105** |
| intact             | +0.048 | +0.059 |
| funmap             | +0.065 | +0.071 |
| rna_coexp          | +0.017 | +0.056 |
| reactome           | −0.035 | +0.060 |

**Reading the table**

- **The embedding beats degree in 5 of 6 networks** (+0.017 to +0.065).
  The exception is `reactome`, where the embedding is negative
  (−0.035) while RWR is positive (+0.060) — the same Reactome weakness
  on the newest drivers that the W5 / A5 controls localise.
- **RWR leads the embedding on every network**, and by 10.3 pp on
  `string_phys700` (0.122 vs 0.019). When the network carries intact
  protein-binding edges only, degree-aware propagation extracts
  substantially more signal than a generic walk-based embedding.
- **The embedding's estimate is noisier**: mean SD 0.051 against RWR's
  0.019, and 5 of its 24 cells are negative (against 0 for RWR).

The conclusion for the paper is therefore *stronger* than the
pre-registered §7 wording anticipated: RWR is not merely a
conservative baseline that is hard to beat, it is the better-performing
learner on this task, at a fraction of the cost and with full
reproducibility.

---

## 4. Why the single test arm matters

The test arm (funmap × intogen2024) shows the embedding at
**+0.148** against RWR's **+0.100** in that cell — the embedding is
ahead in this one cell. If a reviewer argues that RWR's +0.100 is
artificial, the arm shows that an *unsupervised embedding* also
produces a positive margin, i.e. the benchmark's conclusion does not
depend on RWR. The panel-level comparison then shows that the margin is
*larger* with RWR, so the benchmark's choice of learner is conservative
in both directions.

---

## 5. Reproducibility

```bash
conda activate d1
cd <repo root>

# full re-run (uses checkpointing; will skip already-saved cells)
python -u scripts/run_gnn_courtesy.py --n-jobs 4

# aggregate
python scripts/aggregate_gnn_courtesy.py
```

Outputs:
- `results/runs/gnn_courtesy/{net}__{label}__alpha0.5.json` — per-cell
  per-fold AUROCs (idempotent)
- `results/tables/gnn_courtesy_arm.tsv` — aggregate α = 0.5 table
  (24 rows) — **the authoritative table**
- `results/runs/gnn_courtesy/status.json` — progress JSON

---

## 6. Outputs committed

```
scripts/run_gnn_courtesy.py             (the learner and the evaluation)
scripts/aggregate_gnn_courtesy.py        (the aggregate table)
results/runs/gnn_courtesy/*.json         (24 cells, per-fold)
results/runs/gnn_courtesy/status.json    (progress JSON)
results/tables/gnn_courtesy_arm.tsv      (24 rows, α=0.5 aggregate)
docs/gnn_courtesy_log.md                 (this file)
```

---

## 7. Cross-references

* Manuscript §J.4 (discussion of the §7 objection) — the numbers in
  §3 above.
* `docs/pre_registration_amendments.md` amendment 3 — the four defects
  and their fixes.
* `results/tables/gnn_courtesy_arm.tsv` — the `params` column records
  `p: 1.0, q: 1.0, alpha_applies: false` for every row, so the
  equivalence above is auditable from the released table itself.
