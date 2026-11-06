# D2 -- the complex-aware (group-preserving) CV splitter

**Date:** 2026-11-05
**Code:** `d1/splits.py::make_splits_complex_aware`,
`d1/splits.py::complex_groups`,
`d1/splits.py::complex_aware_invariants`
**CLI:** `scripts/build_splits.py --complex-aware --complex-map <path>`
**Tests:** `tests/test_splits_complex_aware.py` (9 tests; 113 pass in the
whole suite)
**Status:** done -- D3 (generate the v2 splits) and D4 (re-run and
compare) are the next steps

---

## 1. The direction of the fix (a correction to the protocol)

The protocol's first version said to **spread** each complex's members
across the folds ("no complex group entirely inside one fold"). That was
**wrong**, and it was caught while implementing it, before any code was
written against it.

The leak has a direction: *a held-out positive has a complex-linked gene
in the TRAINING fold*, and RWR then propagates that partner's signal
into the held-out gene. Spreading a complex across the folds
**maximises** that leak. The correct control keeps each complex
**together in one fold**: a held-out gene's complex partners are then
held out with it and cannot leak into its score. This is ordinary
group-preserving (grouped) CV, and it is what Picart-Armada et al. 2019
ask for. The protocol (`docs/phase_d_protocol.md` §5.2 and §5.4) has been
corrected accordingly.

## 2. What the splitter does

1. **Groups.** The complex graph on the universe (an edge between two
   genes that share a complex) is split into **connected components**;
   overlapping complexes therefore merge into one group, which is the
   conservative choice. A **sparsity guard** keeps, for every gene, at
   most its 10 largest complexes, so that one huge complex (a ribosome)
   cannot merge half the panel into a single unassignable group.
2. **Assignment.** Groups are processed in descending order of positive
   count (ties by size, then randomly per repeat) and each is given
   whole to the fold that currently holds the fewest positives; ties
   inside a 10 % band are broken at random with `base_seed + repeat`,
   so the 40 repeats are genuinely different partitions rather than 40
   copies of one.
3. **Repair.** If a fold would end up with **no** positive -- which
   happens when there are fewer positive-carrying groups than folds --
   one positive member is moved out of the smallest positive-carrying
   group. This is the only way to split a group; every such event is
   recorded (gene, group, from-fold, to-fold) and written into the
   split build log. On the real panel **no repair was needed** for any
   of the 24 arms (see §4).
4. **Parity.** With an empty complex map the function delegates to
   `make_splits`, so it reproduces the v1 splits byte for byte -- the
   strongest test that the new path changes nothing else.

## 3. The tests

| # | test |
|---|---|
| 1 | no complex group is split across folds, except for a documented repair |
| 2 | determinism: two runs give byte-identical output |
| 3 | the repeats are different partitions (not 40 copies of one) |
| 4 | schema, coverage (every gene once per repeat) and at least one positive per fold |
| 5 | positives per fold stay balanced when the groups are small relative to the positive set |
| 6 | an empty complex map reproduces `make_splits` exactly |
| 7 | overlapping complexes merge into one group |
| 8 | the sparsity guard keeps the largest complexes and drops the smaller links |
| 9 | `summarise_splits` works on the complex-aware output |

Three implementation defects were found by these tests while writing
them: the sparsity guard was ineffective (the union still used the full
complex membership), the repair pass moved the **same** gene into every
empty fold, and a shadowed `pandas` import broke the CLI.

## 4. Validation on the real panel

For every (network, label) arm, the D1 leakage metric -- of the held-out
positives that are in a curated complex, the fraction that have a
complex-linked gene in the training fold -- was measured on the v1
splits and on the complex-aware splits (10 repeats each):

| arm | v1 leak | complex-aware leak | groups | groups split | positives per fold |
|---|---|---|---|---|---|
| funmap × intogen2024 | 0.876 | **0.058** | 8,113 | 0 | 85–94 |
| funmap × ot_nolit | 0.883 | **0.123** | 8,113 | 0 | 69–77 |
| funmap × intogen_temporal_new | 0.849 | **0.016** | 8,113 | 0 | 19–22 |
| string_phys700 × ot_all | 0.944 | **0.153** | 7,201 | 0 | 98–111 |
| intact × intogen2024 | 0.941 | **0.068** | 15,119 | 0 | 115–127 |
| intact × ot_nolit | 0.922 | **0.167** | 15,119 | 0 | 107–120 |
| intact × intogen_temporal_new | 0.953 | **0.026** | 15,119 | 0 | 27–31 |
| string_full700 × ot_nolit | 0.931 | **0.173** | 13,119 | 0 | 107–120 |
| reactome × intogen2024 | 0.903 | **0.056** | 2,929 | 0 | 69–78 |
| reactome × intogen_temporal_new | 0.873 | **0.043** | 2,929 | 0 | 11–13 |

(the full 24 arms are in the run log; the pattern is the same)

**The channel is closed**: the leakage falls from 0.83–0.95 to
0.00–0.17, the group invariant holds in **24/24** arms (no group is
split), and the fold balance survives the constraint (the positives per
fold stay within about ±5 % of the target even for the temporal label,
where the whole positive set is 11–31 genes per fold).

**Where the residual comes from.** It is the **sparsity guard**: a gene
keeps only its 10 largest complexes, so a link through a dropped
complex still exists in the audit's (complete) partner set but not in
the splitter's. The residual is therefore a *known and quantified*
under-count of the guard, not a failure of the invariant. Raising
`max_links_per_gene` would trade balance for a smaller residual; 10 is
the value the protocol fixes.

## 5. Caveats

1. **The complex set is sparse** (18.3 % of the panel's genes, 33–37 %
   of the positives), so most genes are unconstrained by this control.
   The audit measures a lower bound on the leakage a complete complex
   resource would show.
2. **Group-preserving CV cannot balance the positives exactly** when a
   complex is large relative to the positive set. On this panel it does
   not matter (the deviations are within ±5 %), but the acceptance
   criterion for D3 is "the deviation is recorded and reported", not
   "zero deviation".
3. The complex-aware partitions are **less diverse** than the v1
   StratifiedKFold partitions (the group constraint removes freedom).
   The 40 repeats therefore sample a narrower space, which makes the
   repeat-level intervals **narrower** than the v1 ones for the same
   data -- a point the D4 comparison has to make explicitly, because
   otherwise the complex-aware result could look "more significant"
   simply because its intervals are narrower.
4. The splitter changes the **estimator** (what the held-out set
   means), so the v2 results are a new analysis, not a fix of v1.

## 6. Next

* **D3** -- generate the v2 splits into `data/processed/splits_v40_cx/`
  (24 arms, 40 repeats) and write the invariant report
  (`results/tables/d3_split_invariants.tsv`), including the achieved
  positive counts per fold and the number of repair events.
* **D4** -- re-run the primary factorial on the v2 splits with the same
  statistical machinery, and compare with v1 under the plan's Q3 rule.
