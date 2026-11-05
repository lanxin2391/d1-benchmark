# S3 -- network-level stability and what it depends on

**Date:** 2026-11-05
**Script:** `scripts/run_network_stability.py`
**Inputs:** `results/runs_v40/a4_combined.tsv` (6 primary networks),
`results/runs_v40_sens8/a4_combined.tsv` (the two STRING sensitivity
arms, run for this phase), `data/processed/networks/*.tsv`,
`results/tables/label_network_overlap.tsv`,
`results/tables/stability_attribution_loo.csv`
**Outputs:** `results/tables/network_descriptors.tsv` (8),
`results/tables/network_stability.tsv` (8),
`results/tables/network_stability_by_arm.tsv` (32),
`results/tables/network_stability_vs_topology.tsv`,
`results/tables/network_stability_loo.tsv`,
`figures/figure6_network_stability.png`
**Status:** exploratory (work plan Phase S, R19; risk R18 applies: eight
networks is a descriptive sample, not an inferential one)

---

## 1. Why the two extra networks

S1 and S2 worked inside an arm. S3 asks the question the teacher
actually asked at the network level: **which networks are stable, and
what makes them so?** With six networks a regression is hopeless, so
the two STRING sensitivity arms (`string_full400`, 19,486 nodes /
929,471 edges, and `string_full900`, 11,693 nodes / 100,383 edges) were
run on the 40-repeat design for this phase, giving n = 8 networks and
32 arms. They are in the curated/literature class, so they also probe
the class label itself.

The new arms were produced with

```
python scripts/build_splits.py --label <label> --universe native_<net> \
       --n-rep 40 --k 5 --out-dir data/processed/splits_v40_sens8
python scripts/run_a4.py --only-net string_full400 --only-net string_full900 \
       --n-rep 40 --splits-dir data/processed/splits_v40_sens8 \
       --runs-dir results/runs_v40_sens8 --table-suffix _v40_sens8
```

(the split directory is regenerable and therefore not committed; the
six-network released tables are untouched).

## 2. The eight networks

| network | nodes | edges | mean ΔAUROC | arm SD (stability) | cross-label spread | mean coverage | 10-core fraction | S2 transfer R² |
|---|---|---|---|---|---|---|---|---|
| string_full400 | 19,486 | 929,471 | +0.0732 | **0.0021** | **0.0129** | 0.996 | 0.906 | – |
| string_full700 | 15,882 | 236,712 | +0.1051 | 0.0034 | 0.0303 | 0.978 | 0.539 | 0.270 |
| intact | 17,918 | 564,706 | +0.0587 | 0.0035 | 0.0353 | 0.983 | 0.762 | 0.129 |
| string_full900 | 11,693 | 100,383 | **+0.1307** | 0.0043 | 0.0398 | 0.915 | 0.353 | – |
| string_phys700 | 9,830 | 85,576 | +0.1215 | 0.0044 | **0.0600** | 0.867 | 0.317 | 0.036 |
| **reactome** | **3,953** | 20,143 | +0.0601 | 0.0054 | 0.0567 | **0.587** | **0.231** | **−1.177** |
| funmap | 10,189 | 196,605 | +0.0708 | 0.0065 | 0.0353 | 0.684 | 0.431 | 0.213 |
| rna_coexp | 9,152 | 196,172 | +0.0562 | **0.0077** | 0.0320 | 0.592 | 0.473 | 0.172 |

## 3. What predicts the stability (Spearman, n = 8, BH-FDR within each target)

| descriptor | arm-level SD (stability) | cross-label spread | mean ΔAUROC (effect size) |
|---|---|---|---|
| **mean label coverage** | **−0.905 (q = 0.012)** | −0.571 (q = 0.167) | +0.310 (q = 0.760) |
| **fraction in the 2-core** | **−0.905 (q = 0.012)** | −0.548 (q = 0.175) | +0.286 (q = 0.760) |
| n nodes | −0.857 (q = 0.026) | −0.643 (q = 0.128) | +0.238 (q = 0.760) |
| n edges | −0.619 (q = 0.239) | **−0.857 (q = 0.039)** | −0.167 (q = 0.832) |
| fraction in the 10-core | −0.595 (q = 0.239) | **−0.905 (q = 0.024)** | −0.238 (q = 0.760) |
| mean degree | −0.310 (q = 0.547) | −0.762 (q = 0.084) | −0.429 (q = 0.760) |
| density | +0.071 (q = 0.911) | −0.643 (q = 0.128) | −0.690 (q = 0.696) |
| mean clustering | +0.452 (q = 0.347) | +0.714 (q = 0.112) | +0.286 (q = 0.760) |
| modularity | +0.476 (q = 0.347) | +0.833 (q = 0.041) | +0.524 (q = 0.760) |
| degree assortativity | +0.500 (q = 0.347) | +0.167 (q = 0.693) | +0.333 (q = 0.760) |
| approximate mean path length | +0.714 (q = 0.140) | +0.571 (q = 0.167) | +0.024 (q = 0.955) |
| max k-core | −0.048 (q = 0.911) | −0.667 (q = 0.128) | −0.048 (q = 0.955) |

Three statements, and the third is the one that matters:

1. **Stability is predicted by how much of the label set the network
   contains and by how much of the network is in the 2-core** (both
   ρ = −0.905, q = 0.012): a network that holds most of the positives,
   inside a well-connected core, gives a reproducible estimate. A
   network that holds a third of the positives gives a noisy one.
2. **Classic topology does not predict stability**: density, clustering,
   modularity, assortativity, max k-core and path length are all
   non-significant after FDR. The stability of a network's contribution
   is not a property of its shape.
3. **No descriptor predicts the SIZE of the contribution.** Every
   correlation with `net_delta_mean` is non-significant (q ≥ 0.70).
   Which network contributes a lot and which network contributes a
   stable estimate are two different questions, and only the second one
   has a topological answer.

A leave-one-network-out check (fit on seven networks, predict the
eighth's stability from the single best descriptor) predicts the
arm-level SD to within 0.0000–0.0014 — but with n = 8 this is a
sanity check, not a result.

## 4. Reactome, in structural terms

Reactome is the smallest network (3,953 nodes), has the smallest
well-connected core (23 % of the nodes in the 10-core) and the lowest
label coverage (0.587 mean, 0.415 minimum). By the rule in section 3
it is therefore *expected* to be one of the least stable networks, and
it is (arm SD 0.0054, second worst of eight; its temporal arm is the
least stable arm of all 32 at 0.0099).

That is a different statement from the negative result, and the two must
be kept apart:

* Reactome's **instability** is explained by coverage and core size --
  an evaluation-set property, shared with rna_coexp and funmap.
* Reactome's **negative temporal delta** is not: the D-09 / D-11 / D-12
  controls equalise the positive set and the coverage and the sign
  stays negative. And S2 showed that Reactome is the only network whose
  *stability pattern* does not transfer (R² = −1.18): its temporal
  failure is a different kind of failure, not a noisier version of the
  same thing.

**The confound to state explicitly**: Reactome's temporal arm has only
63 positives, the smallest positive set in the whole panel, so its
interval is wide for that reason alone; the controls in §F.2–§F.4 exist
precisely because "small positive set" is the obvious alternative
explanation for a negative, noisy cell.

## 5. Caveats

1. **n = 8 networks.** Every number in section 3 is a descriptive
   correlation with a wide interval; the two that survive FDR are
   reported because they are the ones the leave-one-out check
   corroborates, not because they are established. Work plan risk R18.
2. The two STRING sensitivity arms are **not independent** of
   `string_full700` (same data source, different thresholds), so the
   effective sample size is smaller than 8. They are marked as
   sensitivity arms in every table and figure.
3. "Coverage" and "size" are themselves consequences of the LCC
   restriction and of the label choice; the correlation is between
   quantities the benchmark controls, which is why it is a useful
   caution (a network's contribution is more stable when it covers more
   of the label) and not a deep topological law.
4. The arm-level SD is a stability measure of the *estimate*, not of the
   *prior*: it mixes the sampling error of the estimator with the
   intrinsic re-partition sensitivity of the score (which S1 quantifies
   separately).
5. The stability arms' per-arm runs were produced for this phase only;
   they are not part of the released six-network panel, and the
   manuscript's Table 1 and intervals remain the six-network ones.

## 6. What comes next

* **S4** -- ML attribution with stability selection (the PCP + SVM
  analogue), reusing the complex-aware splitter from Phase D.
* **S5** -- repeat S1–S3 on the ClinGen label and the transfer arms.
* **S6** -- manuscript section "§K Where the stability of the prior
  comes from", with Figures 5, 5b and 6.
