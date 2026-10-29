# Per-gene score dumps (A3)

**Task:** work-plan v3.1 A3
**Script:** `scripts/run_a4.py --scores-dir results/scores_v40`
**Code:** `d1/engine/evaluate.py::dump_per_gene_scores`
**Size:** 1133 MB across 24 files (float64)

## Why

Pre-registration section 9.3 asks for a **cluster bootstrap resampling
*genes*, not gene-arm rows**, because gene identity is the unit of
non-independence. The C5 result table only carries **per-fold metrics**
(one AUROC per fold), so the gene-cluster bootstrap could not be
computed from it: to resample genes one needs the score of *every gene*
in *every fold*.

`run_arm` now writes that out. The cost is one extra array per alpha and
no extra computation -- RWR already produces the full
`(n_genes, n_folds)` score matrix.

## How to regenerate

```powershell
uv run python scripts/run_a4.py \
    --splits-dir data/processed/splits_v40 \
    --runs-dir   results/runs_v40 \
    --table-suffix _v40 \
    --n-rep 40 \
    --scores-dir results/scores_v40
```

Wall-clock on this machine: **5.1 minutes** for the whole 24-arm
factorial plus the 24 dumps. The RWR propagation itself is ~75 seconds;
the rest is writing and compressing 1.1 GB.

## File format

One `.npz` per (network, label) arm, named `<net>__<label>.npz`:

| array | shape | dtype | meaning |
|---|---|---|---|
| `genes` | (N,) | object | universe gene names -- the row order of every matrix below |
| `y` | (N,) | int8 | label of each gene (constant across repeats by contract C4) |
| `degree` | (N,) | f64 | degree score of each gene |
| `rwr` | (A, N, F) | f64 | RWR score of each gene in each fold; A = number of alphas, F = number of (repeat, fold) columns |
| `in_test` | (F, N) | bool | True where the gene is in that fold's test set |
| `alphas` | (A,) | f32 | the alpha of each `rwr` slice |
| `repeat_idx` | (F,) | i16 | repeat of each column |
| `fold_idx` | (F,) | i8 | fold of each column |

`in_test` has exactly one True per gene within each repeat (the five
folds of a repeat partition the universe) -- verified for all 24 files.

## How the gene-cluster bootstrap uses this

1. Draw N gene indices with replacement from the universe.
2. For each fold j, take the drawn genes that fall in fold j's test
   set, keeping their multiplicities.
3. Recompute AUROC for RWR and for the degree baseline on that resampled
   set, take the difference.
4. Average over folds; repeat 1000 times; take the 2.5 / 97.5
   percentiles.

The gene is the cluster unit, so a gene's contribution to the estimate
moves together across folds -- which is what the pre-registration
asks for and what the row-level bootstrap gets wrong.

## Validation

`scratch/validate_scores.py` recomputes the AUROC of every
(network, label, alpha, fold) from the dumps and compares it with
`results/runs_v40/a4_combined.tsv`:

```
files checked        : 24
(fold, alpha) checks : 14,400
max |AUROC mismatch| : 0.000e+00
max |delta mismatch|: 0.000e+00
PASS
```

**Exact**, not approximate. (The scores are stored in float64
deliberately: a float32 dump reproduced the table only to 1.1e-4, which
is 3 % of the 0.0031-wide confidence interval at 40 repeats.)

## Version control

`results/scores_v40/` is **git-ignored** -- 1.1 GB of derived data does
not belong in the repository, and it is fully regenerable in 5 minutes
from the released splits plus the released code. The SHA-256 of every
file is recorded in `results/tables/stats_per_gene_scores_manifest.tsv`
(which *is* tracked), so a dump can still be checked for tampering.
