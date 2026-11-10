# D1 Benchmark Public Release Protocol

**Author:** Person B (kakamiku)
**Date:** 2026-10-07
**Status:** Draft (B2 deliverable; M2.8-M4.5 in protocol §6)
**DOI reference:** Pre-registration: 10.5281/zenodo.23121762

This document specifies how a third party runs the D1 benchmark on a
new method and reports a number comparable to ours. It is the public
face of the benchmark — anyone reading the released artefacts should
be able to reproduce the headline claims without contacting us.

---

## 1. What is released

The D1 benchmark release v1.0 consists of:

1. **Pre-registration document** (DOI 10.5281/zenodo.23121762)
   — frozen before any confirmatory run.
2. **Pre-registration submission manifest** with SHA-256 fingerprints
   for every frozen input file.
3. **Six primary networks** (FunMap, IntAct, Reactome, rna_coexp,
   STRING v12.0 full at score ≥ 700, STRING v12.0 physical at
   score ≥ 700) + 2 sensitivity arms (STRING full at scores ≥ 400
   and ≥ 900).
4. **Four primary label sets** (IntOGen 2024.09.20, Open Targets
   26.06 all-evidence, Open Targets 26.06 no-literature,
   intogen_temporal_new) + optional ClinGen 5th label.
5. **Cross-validation splits** (5-fold × 10-repeat = 50 splits per
   (label × universe) pair; **released as files**, not regenerated
   from seeds — the protocol requires this for reproducibility).
6. **The fixed RWR learner** (α = 0.5, tol = 1e-8, max_iter = 200,
   row-stochastic transition matrix, restart probability α).
7. **The evaluation harness** (RWR + degree baseline + AUROC +
   AUPRC + precision-at-k, with D-06 test-set definition and D-07
   precision-at-k tie-break).
8. **All 7 protocol §4.4 controls** wired into the pipeline.
9. **All 4 §4.6 statistical-analysis outputs** (mixed-effects
   variance decomposition, cluster bootstrap CI, equivalence tests,
   BH-FDR).

## 2. What is *not* released (yet)

1. The raw input data (FunMap, STRING, IntAct, Reactome, Open
   Targets parquets, gene2pubmed, IntOGen, HGNC) — these are
   available from the upstream sources listed in
   `resource_and_dataset_table.csv`. The benchmark assumes the user
   has downloaded them following the protocol.
2. The benchmark's trained models — the "models" are the **networks**
   and the **labels**, both of which are released.
3. The 600-cell-per-CPU-hour compute slot — the protocol projects
   the full factorial at ~2 hours on a 4-core laptop; the user
   provides their own compute.

## 3. Release channels

### 3.1 Zenodo

A **single DOI** (`10.5281/zenodo.d1-benchmark-v1.0`) will be minted
on Zenodo for the benchmark release v1.0. This DOI is distinct from
the pre-registration DOI (`10.5281/zenodo.23121762`) and is the one
users should cite in their papers when reporting a number on the
benchmark.

The Zenodo deposit will contain:

* `d1-benchmark-v1.0.zip` — git-tag snapshot of the repo at the
  release commit (full source + data/processed/ artefacts + tests).
* `d1-benchmark-splits-v1.0.zip` — just the cross-validation splits
  (for users who already have their own networks/labels and only
  want the splits).
* `d1-benchmark-results-v1.0.zip` — the 24-cell primary factorial
  results (Table 1 + D-09 + D-11 + D-12 + W5 + A5 + Control 5 +
  Control 7 + provenance ablation + GNN courtesy) for users who
  want to compare without re-running.

Zenodo is preferred over figshare because:

* Zenodo mints a DOI per release; figshare requires manual DOI
  registration.
* Zenodo has a 50 GB limit per record; figshare is 5 GB. Our
  release is ~200 MB compressed, well under both.
* Zenodo is hosted by CERN; figshare is hosted by Digital Science.
  Zenodo is more common for benchmark releases in computational
  biology (e.g., OpenML, UCI ML, Genomics England).

### 3.2 GitHub

The `lanxin2391/d1-benchmark` repository at the `v1.0-prep` tag
is the canonical source. Anyone who clones the tag and follows
the `README.md` instructions will reproduce the headline numbers
to 1e-4. The tag name `v1.0-prep` is intentional: it is the
"pre-publication" form of the v1.0 release; a future
`v1.0-released` tag will be cut at the M6 manuscript-acceptance
milestone and will reference the same commit unless the journal
review process requires a fix.

### 3.3 Bioinformatics resource registry

A submission to the **bio.tools** registry (Elixir tools registry)
will be made under the tool name `d1-benchmark`. This gives the
benchmark a discoverability surface for biologists who would not
naturally find it on Zenodo.

## 4. Leaderboard protocol

The leaderboard is the canonical place where third parties report
their numbers. The protocol is:

### 4.1 What a third-party submission must contain

1. The **method name** (e.g., "GraphSAGE with custom loss").
2. The **hyperparameters** used (any deviation from defaults must
   be reported; defaults are listed below).
3. The **commit hash** of the d1-benchmark release they ran on.
4. The **mean ΔAUROC** for each of the 24 cells in the primary
   factorial (6 networks × 4 labels, α = 0.5, native universe).
5. The **per-fold AUROC values** (50 per cell) so we can compute
   the variance decomposition ourselves.
6. The **runtime** (wall-clock) on their machine.
7. The **environment** (Python version, library versions, OS).

### 4.2 Defaults a third party must use unless they say otherwise

| parameter | value | source |
|-----------|-------|--------|
| α (restart prob) | 0.5 | protocol §4.3 |
| convergence tolerance | 1e-8 | protocol §4.3 |
| max iterations | 200 | protocol §4.3 |
| edge weights | unweighted (binary) | D-01 |
| split seed | 20260926 + repeat | D-BASE |
| split structure | 5-fold × 10-repeat StratifiedKFold | protocol §4.3 |
| test fold scoring | train-positives seeds; test fold is test-positives + test-negatives; seeds excluded from ranked test set | D-06 |
| precision-at-k tie-break | score desc, symbol asc | D-07 |
| HGNC mapping | approved HGNC; previous iff 1-to-1; alias iff unique | D-02 |
| label cap | per-network min across 4 labels | D-11 |
| numpy | >=1.26,<1.27 | D-NUMPY |

### 4.3 What is NOT comparable

Third-party numbers are NOT directly comparable to ours if the
third party:

* Uses a different label version (e.g., IntOGen 2024.06 instead of
  2024.09.20). The pre-reg locks the version.
* Uses a different split seed. The pre-reg locks the seed.
* Uses weighted RWR. The pre-reg locks unweighted.
* Uses a custom universe (e.g., only their protein-coding subset).
  The pre-reg locks the universe per network.
* Re-tunes the hyperparameters on the test fold. The pre-reg
  locks hyperparameters before any test-fold evaluation.

### 4.4 Submission process

Third parties submit to the **d1-benchmark-leaderboard GitHub
issue tracker** at `lanxin2391/d1-benchmark/issues` with the
label `leaderboard-submission`. Each submission triggers:

1. An automatic reproducibility check: the d1-benchmark CI runs
   the submission on a clean machine and verifies the reported
   numbers match.
2. A maintainer review (Person A or B, currently) to check that
   the submission follows §4.2 defaults.
3. Publication to the leaderboard at
   `docs/leaderboard.md` (auto-generated from
   `submissions/*.json`).

### 4.5 Leaderboard format

```markdown
| Method | ΔAUROC funmap × intogen2024 | ... 24 cells ... | Avg | Runtime |
|---|---|---|---|---|
| RWR (locked) | 0.0998 | ... | ... | ... |
| GraphSAGE (XYZ 2024) | TBD | ... | ... | ... |
| Node2Vec (Grover 2016) | TBD | ... | ... | ... |
```

Each cell reports **mean ± SD over 50 folds**. The "Avg" column is
the **fixed-effects inverse-variance meta mean** across the 24
cells (same as A's W2 stats).

## 5. Citation

Anyone reporting a number on the D1 benchmark should cite both:

* The pre-registration:
  `lanxin2391, kakamiku (2026). D1 benchmark pre-registration.
  Zenodo. https://doi.org/10.5281/zenodo.23121762`
* The benchmark release (when published):
  `lanxin2391, kakamiku (2026). D1 benchmark v1.0: a
  leakage-controlled benchmark of network priors for cancer-driver
  prioritisation. Zenodo. https://doi.org/10.5281/zenodo.d1-benchmark-v1.0`

And the relevant method paper for any new method they ran.

## 6. Versioning

The benchmark follows [semantic versioning](https://semver.org/):

* **Major version** (v2.0): a frozen input changes (a new network
  replaces an old one; a label version is upgraded). The pre-reg
  is re-filed.
* **Minor version** (v1.1): new controls are added; existing
  results are unchanged.
* **Patch version** (v1.0.1): bug fixes that do not change
  results (e.g., test fixes, doc fixes).

The v1.0 release is locked at the commit `d9c68e5` (the
post-cleanup merge commit on `main` as of 2026-10-08, tagged
`v1.0-prep`). The README release index points to all 30+
deliverables. The intervening commit `bbf87e1` (the M2-M4
merge) is a strict ancestor of `d9c68e5`; if you need to
reproduce a number from a milestone between M2 and v1.0, use
`bbf87e1`. Any change beyond v1.0.0 requires:

1. A new pre-registration filing (with a new DOI).
2. A new release with a new v1.x or v2.0 tag.
3. The old v1.0 release remains on Zenodo for reproducibility.

## 7. Reproducibility contract

When a third party reports a number on the D1 benchmark v1.0:

* The commit hash on `lanxin2391/d1-benchmark` must be a descendant
  of `d9c68e5` (the `v1.0-prep` tag) or be `d9c68e5` itself.
* The default parameters in §4.2 must be used unless explicitly
  noted as a deviation.
* The pre-registration DOI must be cited.
* The third party must agree to share their per-fold AUROC values
  with us so we can independently verify their variance decomposition.

This is the same contract we ourselves accept: the d1-benchmark
v1.0 paper's headline numbers must be reproducible from the
released artefacts by a third party who was not involved in the
work. The reproducibility audit at `docs/reproducibility_audit.md`
records our end of the contract: 45 PASS / 0 FAIL on the
independent end-to-end replay.

## 8. Maintenance

The benchmark is currently maintained by lanxin2391 and kakamiku.
Long-term, the maintenance plan is:

* **For the duration of the D1 paper submission (M6 + 12 months):**
  the current maintainers respond to issues within 1 week.
* **After acceptance:** the benchmark moves to a community-
  maintained model (e.g., a GitHub org with a steering committee).
  The current maintainers stay on for 6 months as transition.

## 9. Open issues at the time of release v1.0

1. **Zenodo DOI for the benchmark release** is pending; the
   pre-registration DOI is already minted. The release DOI will
   be minted at the M6 milestone (manuscript submission).
2. **The leaderboard page** at `docs/leaderboard.md` is a
   placeholder; the first leaderboard-submission CI hookup is
   pending.
3. **The bioinformatics registry submission** to bio.tools is
   pending; this is a low-effort item.

These are documented in `pending_work_protocol_compliance.md` §13
(remaining B tasks for M4-M6) and will be addressed before the
M6 milestone.

## 10. Cross-references

* Pre-registration: `docs/pre_registration.md` (DOI
  10.5281/zenodo.23121762)
* Submission manifest with SHA-256 fingerprints:
  `docs/pre_registration_manifest.json`
* Resource table: `resource_and_dataset_table.csv`
* Decisions: `docs/decisions.md`
* Pending-work analysis: `docs/pending_work_protocol_compliance.md`
* A-side and B-side logs: `docs/a*_log.md`, `docs/b_log.md`
