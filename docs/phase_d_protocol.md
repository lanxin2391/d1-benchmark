# Phase D protocol — complex-aware cross-validation

**Owner:** Person B
**Status:** ready to execute (written by Person A, 2026-11-05)
**Working plan:** `D:\Grade3\swxxx\final\D1\work_plan_v3.1.md`, Phase D
(D1–D4) + C1
**Context:** `docs/pre_registration_amendments.md` (amendment 2), plan risks
R5 / R20, decision point Q3

---

## 0. Why this exists (one paragraph, for the record)

The benchmark's seven controls cover degree leakage, label leakage,
seed identity, edge rewiring, network coverage, positive-set size and
temporal drift. They do **not** cover **protein-complex leakage**.
Standard stratified k-fold splits put the members of a protein complex
in the same fold, and any method that propagates along the network
(RWR, and every diffusion method) then reads the held-out complex's
members' scores off their training-side partners. This makes the
standard CV optimistic for topology-based methods. Picart-Armada et al.
2019 (*PLoS Comput Biol* 15(12):e1007259) made this criticism of the
field explicitly; our first paper must either answer it or declare it
as a limitation. Phase D answers it.

**Note on the estimation contract:** complex-aware splits change the
*estimator* (what the held-out set means), so they are a **new
analysis**, not a fix. The pre-registration (DOI 10.5281/zenodo.23121762)
is frozen; the complex-aware result is reported as an additional,
exploratory analysis with its own amendment entry, and the pre-registered
result stays as it is.

## 1. Scope and non-goals

**In scope**

* Acquire a gene → protein-complex mapping and map it onto HGNC symbols.
* Measure how much complex leakage the current v1 splits contain.
* Implement a complex-aware splitter, with tests.
* Generate a v2 split set (new directory, v1 untouched).
* Re-run the primary 24-cell factorial on v2 and compare it with v1
  under the **same** statistical machinery.
* Deliver the comparison tables, the figures and the amendment text.

**Out of scope (do not do these)**

* Do not change the learner, the α values, the labels or the networks.
* Do not re-run the control arms (W5 rewiring, A5 seed null, D-09/D-11/D-12)
  unless §6.3 says so.
* Do not touch `data/processed/splits/` (v1), `data/processed/splits_v40/`,
  `results/runs/`, `results/runs_v40/` or any released table.
* Do not edit `docs/pre_registration.md`.

## 2. Deliverables checklist

| # | artefact | path | acceptance |
|---|---|---|---|
| 1 | complex map | `data/processed/complexes/complex_map.tsv` | one row per (complex, gene), columns `complex_id, source, gene, hgnc_id, n_members` |
| 2 | provenance note | `data/processed/complexes/provenance.md` | source URL, release/date, licence, download date, SHA-256 of the raw file, number of complexes and genes |
| 3 | leakage audit | `results/tables/d1_leakage_audit.tsv` | per (network, label, repeat, fold): test positives, of which N have ≥1 complex partner in the training fold; plus per-arm summaries |
| 4 | leakage summary | `results/tables/d1_leakage_summary.tsv` | per-arm: mean and max fraction of test positives with a training-side complex partner; the same for the *negatives* |
| 5 | splitter | `d1/splits.py::make_splits_complex_aware` | unit tests, see §5.3 |
| 6 | v2 splits | `data/processed/splits_v40_cx/` | 24 files (6 networks × 4 labels), 40 repeats × 5 folds, plus `build_log.tsv` |
| 7 | split invariant report | `results/tables/d3_split_invariants.tsv` | per file: n_genes, n_positives, positives per fold, **complex members split across folds** |
| 8 | v2 run | `results/runs_v40_cx/` + `results/runs_v40_cx/a4_combined.tsv` | 24 arms, same protocol as v1 |
| 9 | comparison | `results/tables/d4_v1_vs_v2_comparison.tsv` | per arm: ΔAUROC v1, ΔAUROC v2, difference, both 95 % CIs (repeat-level), both H1 verdicts |
| 10 | statistics on v2 | `results/tables/stats_bootstrap_ci_v2_cx.tsv`, `stats_equivalence_v2_cx.tsv`, `stats_h1_verdict_v2_cx.tsv`, `results_stats_v2_cx.tsv` | same machinery as the v1 tables, same resampling units |
| 11 | figure | `figures/figure3_complex_aware_cv.png` | v1 vs v2 per arm (ΔAUROC with CIs) + the leakage audit panel |
| 12 | write-up | `docs/complex_aware_cv.md` | methods, the leakage numbers, the comparison, the decision |
| 13 | amendment | `docs/pre_registration_amendments.md`, new section **A2.3j** | what was added, that the pre-registration is unchanged, the decision under Q3 |

## 3. Task C1 — the complex data (Person B, ~1–1.5 days)

**Sources** (either is enough; use both if reachable)

* **CORUM** (Comprehensive Resource of Mammalian protein complexes),
  current release, `http://mips.helmholtz-muenchen.de/corum/` —
  organism = *Homo sapiens*, "all complexes" (not only manually curated).
* **ComplexPortal** (EMBL-EBI), `https://www.ebi.ac.uk/complexportal/` —
  human complexes, "As shown" (default) set.

If a source cannot be downloaded from this machine, say so in
`provenance.md` and hand the task back to A — do not substitute a
guess.

**Steps**

1. Download the raw file(s) into `data/raw/complexes/`, record the URL,
   the release date and the SHA-256.
2. Map every complex member to an HGNC symbol:
   * CORUM gives UniProt or NCBI gene identifiers; ComplexPortal gives
     UniProt or Ensembl.
   * `data/raw/reference/hgnc_complete_set.txt` is already in the repo and
     contains `uniprot_ids`, `entrez_id`, `ensembl_gene_id`, `symbol`,
     `alias_symbol`, `prev_symbol`. Build the identifier → symbol map from
     it (one row per identifier, aliases included).
   * Drop members that cannot be mapped; record how many were dropped and
     why. **Do not guess symbols.**
3. Merge the two sources; keep the `source` column so that a
   source-specific analysis stays possible.
4. Write `data/processed/complexes/complex_map.tsv` with one row per
   (complex, gene) and the number of members of that complex.

**Definition to use downstream:** a **complex group** is a set of genes
that are members of the same complex. Two genes are *complex-linked* if
they belong to at least one common complex. Overlapping complexes are
allowed; a gene may belong to several.

**Acceptance:** ≥ 100 human complexes and ≥ 1,000 distinct genes mapped
into the six-network universe; the provenance note explains any shortfall.

## 4. Task D1 — quantify the leakage before fixing it (Person B, ~0.5–1 day)

The plan's rule is **"quantify first, then fix"**. Measure the leakage
that the *current* v1 splits contain.

**Metric.** For each (network, label, repeat, fold), take the held-out
(test) positives and count how many of them have **at least one
complex-linked gene in the training fold**. Report

* `n_test_pos`, `n_test_pos_with_train_partner`, and the fraction,
* the same for the test **negatives** (a control: if the fraction is
  equally high for negatives, the effect is not about the positives),
* the size of the largest complex in the panel (a complex with more than
  k = 5 members cannot be kept whole in 5 folds; that is the hard case
  for D2 and must be reported).

**Also measure** the same quantity for the *shared universe* splits
(D-12) and for the D-11 capped splits, because those are two of the
three split families in the panel.

**Deliverables:** `d1_leakage_audit.tsv`, `d1_leakage_summary.tsv`.

**Acceptance:** the audit covers all 24 v1 arms (6 networks × 4 labels,
40 repeats × 5 folds), and the summary states, per arm, the median and
the maximum fraction of held-out positives with a training-side complex
partner. **If the median fraction is below 5 %, say so explicitly** — that
is itself an important result: it would mean the standard CV is not
materially optimistic for this panel, and Phase D's conclusion would be
"the criticism does not apply here, and here is the measurement".

## 5. Task D2 — the complex-aware splitter (Person B, ~2–3 days)

### 5.1 What it must guarantee

For every (repeat, fold) and every complex group, **the members of the
group must not all be in the same fold**. The hard case is a complex
with more members than folds (a 12-member complex in 5 folds): it *must*
be broken up. The rule below breaks it deterministically.

### 5.2 The algorithm (proposed; B may refine it, but must document any change)

> **CORRECTION (2026-11-05, before implementation).** The first version
> of this section said the opposite of what is needed: it proposed to
> **spread** a complex's members across the folds ("no complex group
> entirely inside one fold"). That is the *worst* possible choice for
> this control, and it was wrong because the leakage has a direction.
> The leak is: *a held-out positive has a complex-linked gene in the
> training fold*, and RWR then propagates that partner's signal into
> the held-out gene. Spreading a complex across folds **maximises**
> that leak. The correct control keeps each complex **together in one
> fold**: then a held-out gene's complex partners are held out with it
> and cannot leak into its score. This is the standard group-preserving
> (grouped) CV, and it is what Picart-Armada et al. 2019 ask for. The
> D1 audit (`docs/d1_leakage_audit.md`) is the check: after D3 the same
> metric must be ~0 for the complex-bearing held-out positives.

1. Build the **complex graph** on the universe: an edge between two
   genes if they share a complex. (If the graph is too dense, use a
   *sparse* graph: keep, for each gene, at most the 10 complex-links
   of its largest complexes — document the choice.)
2. Partition the genes into **connected components** of the complex
   graph. A component is a "group" that must stay together. Overlapping
   complexes merge into one group, which is the conservative choice.
3. Genes **not** in any complex are singleton groups.
4. **Assign whole groups to folds**, balancing the number of positives
   per fold: process the groups in descending order of positive count
   (ties broken by size, then randomly per repeat) and give each group
   to the fold that currently holds the fewest positives (ties: the
   fold with the fewest genes, then the smallest index). The tie
   breaking is randomised with `base_seed + repeat` so that the 40
   repeats are genuinely different partitions.
5. **Balance is a goal, not a hard constraint.** Group-preserving CV
   cannot balance the positives exactly when a complex is large
   relative to the positive set (the temporal label has 12–30 positives
   per fold). The greedy assignment keeps the deviation as small as it
   can; the achieved per-fold positive counts and their deviations are
   **recorded in the invariant report**, and the acceptance criterion is
   "recorded and reported", not "zero deviation".
6. Genes **not** in any complex are assigned by the same greedy
   procedure, so the fold sizes stay comparable with v1.

This keeps the v1 behaviour for everything that is not complex-linked,
which makes the v1-versus-v2 comparison interpretable.

### 5.3 The interface contract

```python
# d1/splits.py
def make_splits_complex_aware(universe, positives, *, complex_map,
                              base_seed=20260926, n_rep=40, k=5,
                              max_links_per_gene=10):
    """Same output schema as make_splits(): columns gene, repeat, fold, y.
    Guarantees: (a) every complex group (connected component of the
    complex graph) lies ENTIRELY in one fold -- the leakage channel
    "a held-out gene has a complex partner in the training fold" is
    therefore closed; (b) the positives per fold are as balanced as
    group preservation allows, and the achieved deviation is reported;
    (c) byte-identical output for identical inputs.
    """
```

and a CLI wrapper in `scripts/build_splits.py`:

```
python scripts/build_splits.py --label intogen2024 --universe native_funmap \
       --n-rep 40 --k 5 --out-dir data/processed/splits_v40_cx \
       --complex-aware --complex-map data/processed/complexes/complex_map.tsv
```

### 5.4 Tests (in `tests/test_splits_complex_aware.py`)

1. **No complex group is split across folds** — every member of a
   complex (of a connected component of the complex graph) is in the
   same fold; checked on a synthetic example with a 12-member complex
   and k = 5.
2. **Determinism** — two runs give byte-identical files.
3. **Positives per fold** are as balanced as group preservation allows
   (the deviation is reported; on the real panel it must be no worse
   than the deviation the greedy assignment achieves).
4. **v1 parity** — with an EMPTY complex map, `make_splits_complex_aware`
   reproduces `make_splits` exactly (the strongest test that the new
   path does not silently change anything else).
5. **Every gene of the universe appears exactly once per repeat**, and
   every fold has at least one positive.

## 6. Task D3 + D4 — generate v2, re-run, compare (Person B, ~2–3 days)

### 6.1 D3 — generate the v2 splits

```
python scripts/build_splits.py --label intogen2024 --universe native_funmap \
       --n-rep 40 --k 5 --out-dir data/processed/splits_v40_cx \
       --complex-aware --complex-map data/processed/complexes/complex_map.tsv
# ... for the 4 labels x 6 networks (and, if time allows, the shared
# universe and the D-11 capped arms)
```

Write the invariant report `d3_split_invariants.tsv`: per file, the
number of genes, the number of positives, positives per fold, the number
of complex groups, the number of groups larger than k, and the number of
fold assignments where a group is entirely inside one fold (**must be 0**
except for the documented exceptions).

### 6.2 D4 — re-run the primary factorial

```
python scripts/run_a4.py --n-rep 40 --splits-dir data/processed/splits_v40_cx \
       --runs-dir results/runs_v40_cx --table-suffix _v40_cx
```

**Compute cost:** the 24 arms of the primary factorial take about 15
minutes of wall clock on the 40-repeat design; the statistics below take
a few minutes. Machine: `E:\anaconda\python.exe` (numpy 1.26, pandas
2.2, scipy 1.13). **Use `--n-jobs` ≤ 4** on this 15.7 GB machine; higher
values cause paging (an earlier run with `--n-jobs 12` spent 8 hours at
1.7 % CPU).

### 6.3 The comparison and the statistics

Compute, for the v2 splits and with **exactly the same machinery as v1**
(same resampling units, same interval families, same TOST margin, same
BH-FDR family):

* `stats_bootstrap_ci_v2_cx.tsv` — the five interval families
  (fold percentile, repeat percentile, repeat t, repeat BCa, gene
  cluster) for the 24 arms.
* `stats_equivalence_v2_cx.tsv` — the TOST at margin ±0.05.
* `stats_h1_verdict_v2_cx.tsv` — the pre-registered H1 rule per interval
  family.
* `results_stats_v2_cx.tsv` — the 24-cell BH-FDR table (fold and repeat
  units, as in `results_stats.tsv`).
* `d4_v1_vs_v2_comparison.tsv` — per arm: ΔAUROC under v1 and under v2,
  the difference, both repeat-level 95 % CIs, and both verdicts.

**Re-use the existing scripts** (`scripts/run_statistics.py`,
`scripts/run_ci_variants.py`, `scripts/run_gene_cluster_ci.py`,
`scripts/results_stats.py`); add a `--splits-dir`/`--tag` argument
rather than writing new estimators. Any change to an estimator must be
applied to v1 as well and the v1 numbers re-verified — a v1/v2
difference must never be an artefact of different code paths.

### 6.4 Decision rule (plan Q3)

| ΔAUROC change (v2 vs v1, repeat-level) | consequence |
|---|---|
| < 5 % | **the complex-aware result goes into the first paper** as a positive selling point ("our control is stricter than the field's standard") |
| 5–10 % | both are reported; the paper states the standard CV is optimistic by X % and uses the conservative number for the headline |
| > 10 % | the complex-aware result is **not** in the first paper; it becomes the second paper's main result, and the first paper lists it as a limitation |

In all three cases the decision, the numbers and the date go into
amendment 2 as **A2.3j** and into `docs/complex_aware_cv.md`.

## 7. Do-not-touch list

```
data/processed/splits/            # v1 splits -- frozen
data/processed/splits_v40/         # the 40-repeat design of record -- frozen
results/runs/  results/runs_v40/   # released run artefacts -- frozen
results/tables/*.tsv              # released tables -- do not overwrite;
                                   #   write new files with a _cx suffix
docs/pre_registration.md           # frozen pre-registration
d1/engine/*.py                    # the learner and the interval code
```

## 8. Time plan and critical path

| task | effort | depends on |
|---|---|---|
| C1 complex data + symbol mapping | 1–1.5 d | — |
| D1 leakage audit (quantify) | 0.5–1 d | C1 |
| D2 splitter + tests | 2–3 d | C1 (the map), can start on synthetic data in parallel |
| D3 v2 split generation + invariants | 0.5 d | D2 |
| D4 re-run + statistics + comparison | 1–1.5 d | D3 |
| write-up + amendment A2.3j + figure | 1 d | D4 |
| **total (critical path)** | **6–9 working days** | |
| minimal version (D1 + D2 + primary factorial only, no controls) | 3–4 d | |

**Machine time is not the bottleneck**: the 24-arm factorial is ~15 min
and the statistics a few minutes. The bottleneck is C1 (the data) and D2
(the splitter), both of which need care rather than compute.

If the full version is not affordable, the **minimal version is enough
to answer the reviewer**: D1 (how much leakage was there), D2/D3 (the
new splits), and the primary factorial re-run. The control arms only
need re-running if the decision rule says the conclusion changes.

## 9. What to hand back to Person A

1. `docs/complex_aware_cv.md` (methods + numbers + the decision).
2. The tables in §2 (rows 3, 4, 7, 9, 10) and the figure.
3. The split invariant report, **with the count of groups that could not
   be split** — if that count is large, the result must be read with
   that caveat.
4. One line for the manuscript: the measured leakage and the change in
   ΔAUROC.
5. The amendment text for A2.3j (do not edit the pre-registration).
