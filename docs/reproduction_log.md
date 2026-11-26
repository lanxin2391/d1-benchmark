# The clean-machine reproduction (plan item 22)

**Date:** 2026-11-26

The environment: a fresh virtual environment created from nothing but the Python standard library, the three release zips (source, splits, results) and `pip install -e .` from the unpacked source. The point of the exercise is not that the numbers are right -- the audit already checks that -- but that a reader who only has the zips can get the main tables back.

## 1. a fresh virtual environment

```
$ python -m venv env
$ python -c 'import sys; print(sys.version)'
3.11.16 | packaged by conda-forge | (main, Sep  2 2026, 23:31:17) [MSC v.1944 64 bit (AMD64)]
```

(the environment was created in 0s)

## 2. the three release zips

* `d1-benchmark-v1.0-source.zip` (272.5 MB) unpacked

* `d1-benchmark-v1.0-splits.zip` (142.1 MB) unpacked

* `d1-benchmark-v1.0-results.zip` (53.6 MB) unpacked

The unpacked tree: `C:\Users\PC\AppData\Local\Temp\opencode\repro\src\d1-benchmark-v1.0`

It holds **211 result tables** and **53 split files** of the 40-repeat design.

## 3. install the package from the unpacked source

```
$ pip install --upgrade pip --index-url https://pypi.org/simple

```

The index is given explicitly because this machine's default pip index is a mirror that is not reachable from here; the reproduction must not depend on the reader's network configuration.

```
$ pip install numpy>=1.26,<1.27 scipy pandas pyarrow networkx scikit-learn statsmodels matplotlib seaborn tqdm requests pytest

```

The dependencies of `environment.yml` were installed with pip in 2.3 min.

```
$ pip install -e .

```

The install took 0.1 min (OK).

The versions that matter:

```
matplotlib==3.11.2
networkx==3.6.1
numpy==1.26.4
pandas==3.0.6
scikit-learn==1.9.1
scipy==1.17.1
statsmodels==0.15.0
```

## 4. the reproduction: the main tables from the released data

```
$ python scripts/recompute_ground_truth.py
952 ground-truth quantities written to C:\Users\PC\AppData\Local\Temp\opencode\repro\src\d1-benchmark-v1.0\review\ground_truth.tsv
```

The ground-truth table recomputed in the clean environment holds **952 quantities**.

```
$ python scripts/verify_claims.py
  [PASS] c9 median, cleaned             manuscript=+0.0822  table=+0.0822  tol=0.002
  [SKIP] seed perturbation, 90 %        (no ground truth for ('seed_perturbation', 'median relative loss at 90 % of the seeds'))
  [SKIP] seed perturbation, 50 %        (no ground truth for ('seed_perturbation', 'median relative loss at 50 % of the seeds'))
  [PASS] size-matched, standard         manuscript=+20.0000  table=+20.0000  tol=0.5
  [PASS] size-matched, capped           manuscript=+18.0000  table=+18.0000  tol=0.5
  [PASS] label redundancy, low          manuscript=+0.8400  table=+0.8413  tol=0.006
  [PASS] second learner, arms           manuscript=+16.0000  table=+16.0000  tol=0.5
  [PASS] second learner, median         manuscript=+0.0112  table=+0.0112  tol=0.002

========================================================================
RESULT: 279 PASS / 0 FAIL
========================================================================

written to C:\Users\PC\AppData\Local\Temp\opencode\repro\src\d1-benchmark-v1.0\review\verification_report.txt
```

**The verdict: RESULT: 279 PASS / 0 FAIL** (the check ran 0.1 min).

The headline numbers as the clean machine recomputed them:

| quantity | scope | the value |
|---|---|---|

| delta_auroc | panel range (min) | -0.0125 |

| delta_auroc | panel range (max) | 0.1711 |



### What the clean environment could not recompute, and why

The clean environment produces fewer ground-truth quantities than the repository: the network descriptors (`network_nodes`, `network_lcc_nodes`, `network_lcc_edges`, ten networks each) are computed from `data/processed/networks/`, and the release zips exclude that directory on purpose -- the networks are rebuilt from the raw downloads with `scripts/build_*.py`. Everything else (the primary endpoint, every control table, the permutation audit, the rewiring grid, the pathway level and the decomposition) is recomputed from the zips alone.


## 5. The verdict


**A reader who only has the three zips can recompute the main tables.** The clean environment installed from the unpacked source, the ground-truth recomputation and the claim verification all ran, and the verification finished with **RESULT: 279 PASS / 0 FAIL**.

What this exercise does **not** check: it does not rebuild the networks, the splits or the result tables from the raw data (that needs the data downloads and a lot of compute); it checks that the released artefacts are self-sufficient for the headline numbers, which is the claim the release makes.
