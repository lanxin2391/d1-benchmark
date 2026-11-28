# Pre-Registration v2 (draft): D1 Benchmark on Network Priors for Cancer-Driver Prioritisation

**Status:** draft, to be uploaded to Zenodo as version 2 of the record whose
version 1 is frozen at DOI 10.5281/zenodo.23121762.
**Prepared:** 2026-11-26.
**Relation to version 1:** version 1 is unchanged and remains the
pre-registration of the confirmatory analysis. This version 2 adds (a) the
complete register of the exploratory analysis family, (b) the register of
the corrections and deviations, and (c) the register of the controls
added after the defence. **It does not change the pre-registered
endpoint, the decision rule, the network panel, the label panel, the
learner or the cross-validation design.**

---

## 1. What version 1 fixed, and what it still says

Version 1 (frozen at DOI 10.5281/zenodo.23121762, commit `d9c68e5`) fixed
the following before any data were collected, and none of it changes
here:

* the research question and the three hypotheses (H1: the network prior
  contributes signal beyond node degree; H2: the contribution is invariant
  to network provenance; H3: apparent between-network differences are
  literature attention);
* the primary endpoint: the per-fold paired z-score of the margin
  delta-AUROC = AUROC(RWR at alpha = 0.5) − AUROC(degree), with the
  5-fold x 40-repeat design, the label-level Bonferroni and the
  gene-cluster bootstrap as the two resampling units;
* the network panel (six networks) and the label panel (four catalogues);
* the learner (RWR, alpha = 0.5, no tuning), the cross-validation
  (StratifiedKFold, the locked base seed) and the inclusion criteria;
* the decision rule for the H1 verdict and the multiplicity control
  (Benjamini-Hochberg over the 24 cells).

**The confirmatory result of version 1 is unchanged by this document.**

## 2. Why a version 2 exists

Two reasons, both about transparency rather than about the endpoint.

1. **The exploratory family was not registered.** Version 1 registered the
   confirmatory analysis; the controls that were run afterwards (the
   leakage controls, the design controls, and the twenty-odd controls
   added after the defence) were run in the exploratory family but never
   written down as a register. Section 4 is that register.
2. **Several corrections and deviations happened after the freeze**, and
   the honest thing is to list them in the record itself rather than only
   in the amendment log (Section 5) and the release notes.

## 3. The pre-registered endpoint (unchanged)

| item | the value |
|---|---|
| the primary endpoint | delta-AUROC = AUROC(RWR, alpha = 0.5) − AUROC(degree), per fold |
| the design | 5-fold x 40-repeat, StratifiedKFold, base seed 20260926 |
| the network panel | funmap, rna_coexp, string_phys700, intact, string_full700, reactome |
| the label panel | IntOGen 2024, Open Targets (all), Open Targets (no literature channel), IntOGen (newest drivers) |
| the resampling units | the label level (Bonferroni over the 24 cells) and the gene cluster (the bootstrap) |
| the multiplicity control | Benjamini-Hochberg over the 24 cells |
| the decision rule for H1 | the interval excludes zero under both units; the equivalence test is reported as a null cell |

## 4. The register of the exploratory analysis family

Everything in this section is **exploratory**: it is reported with its own
intervals, it is not part of the H1 verdict, and it cannot change it.
The families, with the script that produces each one.

### 4.1 The design and the implementation repairs (pre-flight)

| the control | the script | the question |
|---|---|---|
| the complex-aware cross-validation | `scripts/run_d3_generate_cx_splits.py`, `scripts/run_d4_compare_v1_v2.py` | does the leakage from protein complexes change the endpoint? |
| the degree-capped label (D-11) | `scripts/build_splits_v40.py` (the capped panels) | is the margin an artefact of the label's size? |
| the size-matched panel | `scripts/run_control15_16.py` | does the margin survive equal-positive-size labels? |
| the seed perturbation | `scripts/run_control8_seed_degree.py` | how much of the margin survives losing seeds? |
| the size adjustment | `scripts/run_control9_size_adjustment.py` | is the margin explained by the network's size? |
| the random-input control | `scripts/run_control10_random_input.py` | is the result an artefact of the input? |
| the gene-identity null | `scripts/run_a5.py` | is the result an artefact of the gene universe? |

### 4.2 The controls against alternative explanations

| the control | the script | the question |
|---|---|---|
| Control 2, the rewiring null | `scripts/run_rewiring.py` | is the increment carried by the degree sequence? |
| Control 3, the degree-matched seed null | `scripts/run_a5.py` | is the increment carried by the seeds' degree? |
| Control 7, the publication-matched null | `scripts/run_control7_matched_null.py` | is the increment carried by the literature attention? |
| Control 8, the seed perturbation | `scripts/run_control8_seed_degree.py` | how stable is the increment? |
| Control 13, the cleaned background | `scripts/run_c9_clean_background_factorial.py` | is the increment an artefact of the negative set? |
| Control 14, the no-database channel | `scripts/build_string_nodb.py` | is the increment an artefact of the STRING literature channel? |
| Control 17, the second learner | `scripts/run_second_learner.py` | is the increment specific to the diffusion? |
| Control 19, the baseline family | `scripts/run_baseline_family.py` | is degree the strongest classical baseline? |
| Control 20, the label-difficulty regression | `scripts/run_label_difficulty_regression.py` | are the four catalogues four voices or one? |
| Control 21, the visibility decomposition | `scripts/run_c21_visibility_second_learner.py` | how much of the increment is visibility? |
| Control 22, the rewiring null under the second learner | `scripts/run_c22_rewiring_second_learner.py` | does the rewiring conclusion hold under another learner? |
| Control 23, the split-design comparison | `scripts/run_c23_split_design.py` | does the verdict depend on the cross-validation design? |
| Control 24, the functional task read as a ranking | `scripts/run_c24_depmap_within_degree.py` | can the network find functionally important genes? |
| Control 25, mutation-first versus network-first | `scripts/run_c25_mutation_vs_network.py` | does the network add to the mutation evidence? |
| Control 26, the label-permutation audit | `scripts/run_c26_label_permutation.py` | does the pipeline manufacture the margin? |
| Control 27, the rewiring grid under the second learner | `scripts/run_c27_rewiring_second_learner.py` | does the rewiring conclusion hold across the panel? |
| Control 28, the topology prior suite | `scripts/run_c28_topology_suite.py` | is the increment generic to topology? |
| Control 30, the DepMap essential-gene catalogue | `scripts/run_c30_depmap_essential_label.py` | what does the protocol say on a functional answer key? |
| Control 31, the cross-method complementarity | `scripts/run_c31_cross_method.py` | do the methods pick the same genes? |
| Control 32, the pathway level | `scripts/run_c32_pathway_level.py` | do the candidates concentrate on the drivers' pathways? |
| Control 33, the decomposition of the network's contribution | `scripts/run_c33_decomposition.py` | is there a network contribution beyond degree and visibility? |
| Control 38, the H2 rescue | `scripts/run_c38_h2_rescue.py` | can the class effect be separated from the size? |

### 4.3 The register of the corrections

| the correction | the amendment |
|---|---|
| the SHA-256 drift on the four frozen inputs | amendment 1 |
| the variance decomposition re-specified; the CI resampling unit changed | amendment 2 |
| the GNN courtesy arm's four protocol defects | amendment 3 |
| the family registration of every analysis | amendment 4 |
| the three literature-mandated controls | amendment 5 |
| the two adversarial controls on the headline; the H2 verdict restated | amendment 6 |
| the DepMap reading corrected (the saturated within-stratum statistic) | amendment 7 (to be filed with this version) |
| the two z-score conventions separated | amendment 7 (to be filed with this version) |
| the label-permutation audit, the temporal transfer matrix, the topology suite, the pathway level and the decomposition | amendment 7 (to be filed with this version) |

## 5. The deviations, and why they were made

| the deviation | the reason |
|---|---|
| the released split panel was replaced by the 40-repeat design | the gene-cluster bootstrap needs the repeats as its sample size; the amendment records that the first ten repeats are identical to the released panel, row by row |
| the confirmatory endpoint's interval is reported under two units | the folds inside a repeat are not independent; both are reported |
| the "H2 is not testable" verdict | six networks cannot separate provenance from size and coverage; the collection of eleven networks (Control 38) quantifies it rather than resolving it |
| the HotNet2 arm is a one-arm feasibility result | the method's personalised-PageRank matrix is O(n^2) and only fits the smallest network; the successor NetMix2 is registered as the intended comparison |

## 6. The data and the code (unchanged)

* the pre-registration manifest is frozen (`docs/pre_registration_manifest.json`);
* the released tables, the split manifests and the release manifest are
  fingerprinted in `docs/release_manifest_v1.0-rc14.json`;
* the audit scripts: `scripts/verify_claims.py` (the headline numbers
  against the released tables), `scripts/verify_report.py` (the report
  against the tables) and `scripts/audit_reproducibility.py` (the
  fingerprints);
* the clean-machine reproduction log: `docs/reproduction_log.md`.

## 7. The Zenodo metadata for the upload

| the field | the value |
|---|---|
| the title | D1 Benchmark: a leakage-controlled decomposition of the network prior's contribution to cancer-driver prioritisation -- pre-registration version 2 |
| the version | 2 (a new version, with version 1 linked as the related identifier) |
| the creators | the four authors, with their ORCIDs |
| the related identifiers | `10.5281/zenodo.23121762` (version 1) |
| the description | the abstract of the paper, plus the note that version 2 adds the exploratory register and the corrections and does not change the confirmatory endpoint |
| the licence | the same as version 1 |
| the files | `docs/pre_registration_v2.md` and the amendment log |

## 8. The one thing that still needs a human hand

The upload itself: Zenodo's API refuses this machine (HTTP 403 ✓), so the
record has to be created in a browser ✓ and the minted DOI pasted back
into the README, the release note and this document ✓.
