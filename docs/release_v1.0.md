# D1 benchmark — v1.0-rc14 release note

**Release commit:** `1e9e381638b66b4cf676c703b7f7393eca563f45`
**Built:** 2026-11-26
**Zenodo DOI:** *pending* — the three zips have to be uploaded from a
browser, and Zenodo mints the DOI at that point. The pre-registration
(DOI 10.5281/zenodo.23121762) is the frozen v1 and is unchanged.

## What changed since v1.0-prep

| area | change |
|---|---|
| P0 repairs | all 11 done: the alpha filter, the TOST rule, the D-11 filename collision (plus a second bug: the cap was applied before the universe restriction), the D-12 shared universe, the D-09 seed resampling, the variance decomposition, the gene-cluster bootstrap, the 24-cell BH-FDR, the H1 rule as re-computed, and the control-2 null moved to the 40-repeat design |
| Phase D | the complex-aware cross-validation: the complex map (ComplexPortal via UniProtKB), the leakage audit of the v1 splits (96 % of the complex-bearing held-out positives have a partner in the training fold), the group-preserving splitter, the v2 splits and the v1-vs-v2 comparison (the effect shrinks by a median of 12 %, the H1 verdict is unchanged) |
| Phase S | the stability atlas: per-gene stability, its attribution to the sampling noise of the seed mean, the network level, and the replication on ClinGen |
| direction 1 | the consensus prior over the six networks, and its robustness checks |
| E3 | the hub-capped RWR: the headline is not a hub effect, but two thirds of the raw margin sits in the hubs' ranking position |
| P1 | the rewiring metadata now records the swaps that were executed, and the mixing was measured |

### The three zips (v1.0-rc14)

| file | the size | the SHA-256 |
|---|---|---|
| `release/d1-benchmark-v1.0-results.zip` | 53.6 MB | `66e4b71244e44885d2da522dc9e5b9014f3973e4807bc9ea7899e53cda798cb6` |
| `release/d1-benchmark-v1.0-source.zip` | 272.5 MB | `20fff70add854ffff73f178e34e59ae1b12d00605b0760b444f7d82b84975a73` |
| `release/d1-benchmark-v1.0-splits.zip` | 142.1 MB | `4b14408c28324a62454f2e2f2b6f9d0c813551f2dc6c86f36d50ae5328921a95` |

They are built from the tag with `python scripts/build_zenodo_release.py --tag v1.0-rc14`; the upload to Zenodo (and the DOI that Zenodo mints) is a manual step.

## What changed since v1.0-rc13

* **Controls 24-33** were added and written into the manuscript (F.26-F.34),
  the report draft, the facts sheet, the question bank and the amendment
  register: the functional task read as a ranking (24), mutation-first
  versus network-first (25), the label-permutation audit (26), the rewiring
  null under the second learner over the full grid (27), the topology prior
  suite (28), the label-permutation audit (26), the DepMap essential-gene
  catalogue as a fifth label (30), the cross-method complementarity (31),
  the pathway level (32) and the decomposition of the network's
  contribution (33).
* **Items 13 and 14** (the executable protocol and the clinical-use
  boundary) were written as `docs/protocol_network_to_experiment.md` and
  as the manuscript's section M.
* The DepMap reading was corrected (the released "the degree baseline
  69.2 %" is a saturated within-stratum statistic, not predictive power),
  and the two z-score conventions were separated.

## What changed since v1.0-rc9

| area | change |
|---|---|
| the title | one sentence, no punctuation inside, written in the shape the source papers use ("Comprehensive Characterization of...", "Mapping...", "Benchmarking..."), and the Chinese and English correspond word for word: **《癌症驱动基因排序中网络先验贡献的定量分解》** / **Quantitative decomposition of the contribution of network priors to cancer-driver gene prioritisation**. The paper, the deck cover, the spoken script and the PDF all carry the same name |

## What changed since v1.0-rc8

| area | change |
|---|---|
| the title | one sentence, no punctuation inside, the way a paper title is written: **"What a network prior contributes to cancer-driver gene prioritisation"** (the Chinese name is 《网络先验到底贡献了多少》). The "leakage-controlled benchmark" part moved into the abstract and into the deck's subtitle, where it already was |
| the 20-minute script | the spoken text no longer carries document-style labels ("一句话结论：" and friends) -- everything inside the spoken paragraphs can be read out as it stands; the opening now starts with the conclusion itself |
| the length | the 20-minute talk was measured and trimmed: **3,793 Chinese characters + 242 numbers**, about 19-20 minutes at 240-260 characters per minute, with the rule "start the closing at 18:30" written at the top |

## What changed since v1.0-rc7

| area | change |
|---|---|
| the paper draft, as a printable PDF | `docs/research_report_draft.pdf` (39 pages, about 17,500 words) is now in the repository and in the report folder -- the TA needs to print it |
| the AI-fusion section | a new section in the manuscript and the integrated draft, **"What an AI-fused design would look like"**, with three designs that follow from the findings: (1) a **visibility-aware residual prior**, in which a gradient-boosted model acts as a null generator and only the out-of-fold residual is scored; (2) **difficulty-calibrated** cross-catalogue comparison, in which every arm is a standardised effect against its own null; (3) an **LLM-assisted evidence layer** with a pre-registered rubric and blinded adjudication, which is the only route to a label axis that is not derived from the catalogues being tested. Each design states what it fixes AND what it does not fix |
| the 20-minute script | `docs/report_1118_speaker_script_20min.md` -- a Chinese word-for-word script for the 18-20 minute slot the teacher asked for, with a table mapping the teacher's ten required elements to the sections and the slides, a slide plan (15 body slides, the evidence pages as backup), 13 number cards, a three-person allocation, a 10-minute compression plan and a rehearsal checklist. The pre-defence sweep now checks it as well |
| the materials guide | `docs/report_1118_materials_guide.md` -> `材料使用说明.docx`: rewritten for the two slots (20 and 40 minutes), the PDF, the new slide list (A1-A15) and the 32 number cards |
| the deck | 31 slides: appendix A15 carries the three AI-fusion designs, and appendix A14 now states the rewiring feasibility check |

## What changed since v1.0-rc6

| area | change |
|---|---|
| Control 22 -- the rewiring null under the second learner | a **feasibility check on one arm** (funmap x intogen2024), reading the 20 rewired networks Control 2 had already produced: observed **+0.0804**, rewired null **-0.0124 +- 0.0072**, **z = +12.9** (the same arm under RWR: z = +94.2). Every rewired network gives a negative margin, and the rewiring preserves the degree sequence **exactly** (the largest degree difference over the 20 rewirings is 0). The full grid (6 networks x 4 labels x ~100 rewirings) is the first item of the validation plan. `scripts/run_c22_rewiring_second_learner.py`, `results/tables/c22_rewiring_second_learner*.tsv` |
| the claim verification | `scripts/recompute_ground_truth.py` records **443** quantities (was 436) and `scripts/verify_claims.py` runs **144 checks** (was 140) |
| the release manifest | it covers **189 tables + 14 figures** (was 187 + 14) |

## What changed since v1.0-rc5

| area | change |
|---|---|
| Control 21 -- the visibility decomposition under the second learner | the bait-matched null recomputed with the proximity learner as the scorer (same splits, folds, publication counts and matching rule; the matcher reproduces the published G3 statistic exactly on the RWR dumps). RWR 17.6 % -> 11.6 %, proximity 7.7 % -> 4.8 %, and **the sign of the drop agrees in 24 of 24 arms**; the newest-driver label falls to chance under both learners. **The size of the topology component is learner-dependent; the visibility component is not.** `scripts/run_c21_visibility_second_learner.py`, `results/tables/c21_visibility_*.tsv` |
| the headline | it is now an interval rather than a point: **+0.081** against degree, **+0.046** against the strongest classical baseline, **+0.011** under a second locked learner |
| a methodological correction inside Control 21 | the proximity learner's score is an integer distance, so almost every test gene ties; G3's strict comparison counts a tie as "not above", which inflates the statistic (91 % against an AUROC of 0.627). The comparison is made on the protocol's ranks (decision D-07). RWR has few ties, so its published numbers are unchanged, and the self-check confirms it |
| the speaker script | rewritten in Chinese (the slides stay English): 7,530 characters plus 326 numbers, about 38-40 minutes; the English draft is kept as `docs/report_1118_speaker_script_en.md` |
| the claim verification | `scripts/recompute_ground_truth.py` records **436** quantities (was 415) and `scripts/verify_claims.py` runs **140 checks** (was 136) |
| the release manifest | it covers **187 tables + 14 figures** (was 185 + 14) |

## What changed since v1.0-rc4

| area | change |
|---|---|
| Control 19 -- the baseline family sweep | three more classical baselines on the same folds and the same test genes (closeness, eigenvector centrality, PageRank) plus the literature's default PageRank damping. **Closeness beats degree in 24/24 arms**, and against the strongest baseline in each arm the prior is positive in **19/24 arms** (median **+0.0461** against +0.0809 against degree). The headline is now stated with that qualification. `scripts/run_baseline_family.py`, `results/tables/baseline_family_*.tsv`, `docs/baseline_family.md` |
| Control 20 -- the label-difficulty regression | the label main effect of the variance decomposition is **one difficulty axis, not five curation philosophies**: a single label-level number recovers **97.8 %** of what the full label factor explains (F-test p = 0.38). `scripts/run_label_difficulty_regression.py`, `results/tables/label_difficulty_*.tsv`, `docs/label_difficulty_regression.md` |
| the H2 disclosure | the manuscript has always said that H2 (provenance invariance) is **not testable** on this panel (3/1/2 networks per class, provenance confounded with size and coverage, smallest attainable Kruskal-Wallis p = 0.1173). The defence material now says the same, in one place each: deck slides 2, 14 and 15, the speaker script's opening and number card 30, Q&A bank C22, and the plain-language report's 3.4 |
| the claim verification | `scripts/recompute_ground_truth.py` now records **415** quantities (was 296) and `scripts/verify_claims.py` runs **136 checks** (was 104). Two of the new checks are FORBIDDEN entries for the superseded "the margin shrinks by half" sentence, which had survived the C9 correction in the abstract and the integrated draft; both are corrected now, and a revert fails the audit |
| the release manifest | it covers **185 tables + 14 figures** (was 144 + 14) |
| the repository | the per-gene score dumps of the cap and no-database arms (682 MB of derived data) are untracked again, and `.gitignore` carries one general rule for every `results/scores_*` directory. The source zip is 216.7 MB again |
| the defence material | the deck has 29 slides (the body is still 16 pages; two evidence pages were added), the speaker script runs to 39.5 minutes with 30 number cards, the Q&A bank has 51 + 3 questions and eleven defensive facts, and the plain-language report is a separate deliverable |

## What changed since v1.0-rc3

| area | change |
|---|---|
| the G/C series | the bait-matched null (G3), the NCG arbitration test (G7), the negative-set audit (G8), the cancer-type stratification (GR-9), the NCG circularity check (C2), the fifth label in the main factorial (C5), the variance decomposition on the 30-arm grid (C6), the label-quality tiers and the ClinGen noise anchor (C7) — all exploratory, all registered in amendment 4 |
| the fifth label | the NCG known tier (708 genes) enters the factorial as a pre-planned exploratory extension (6 arms, 18 cells); the D-11 cap does not move, so the 24 released arms were **not** re-run |
| the stability analysis | the Phase S attribution is now in the manuscript with the released table's numbers; two stale numbers in section I.2 (the crowding definitions) and one untraceable coefficient were corrected, and `docs/stability_attribution.md` now carries the released table |
| the claim verification | `scripts/recompute_ground_truth.py` now records 296 quantities (was 55) and `scripts/verify_claims.py` runs 104 checks (was 63) — the Phase S and G/C numbers are machine-checked against the tables, and three inconsistencies were found and fixed in the process |
| the release manifest | it now covers every tracked table and figure (144 + 14, was 35 + 9), so the tamper-evidence record describes the whole release |
| the defence material | the 1118 deck (26 slides), the 40-minute speaker script, the Q&A bank and the public dataset list are in the repository under `docs/` |

## What is in the three zips

| zip | size | content |
|---|---|---|
| `d1-benchmark-v1.0-source.zip` | 210.7 MB | the full repository at the release tag, minus `data/raw/`, `data/cache/`, `data/processed/networks/` and `data/processed/rewired/` (all regenerable from the build scripts) |
| `d1-benchmark-v1.0-splits.zip` | 142.1 MB | every split family: the v1 10-repeat set, the 40-repeat design of record, the complex-aware v2 and the ClinGen 40-repeat splits |
| `d1-benchmark-v1.0-results.zip` | 47.1 MB | `results/tables/` and `figures/` |

## Verifying the release

`release_manifest.json (inside the source zip) and docs/release_manifest_v1.0-rc1.json (in the repository)` (inside the source zip) holds the
SHA-256 of the 23 key tables, the 5 figures and the 3 zips:

| zip | SHA-256 |
|---|---|
| `d1-benchmark-v1.0-results.zip` | `2e6eae57617310e7...` |
| `d1-benchmark-v1.0-source.zip` | `b45dc086e345dbee...` |
| `d1-benchmark-v1.0-splits.zip` | `10ed09a3ea608841...` |

```bash
python -c "import hashlib,sys;print(hashlib.sha256(open(sys.argv[1],'rb').read()).hexdigest())" d1-benchmark-v1.0-source.zip
```

## Reproducing the headline numbers

```bash
pytest                                        # 113 tests
python scripts/run_a4.py --n-rep 40 --splits-dir data/processed/splits_v40 \
       --runs-dir results/runs_v40 --table-suffix _v40
python scripts/run_statistics.py --combined results/runs_v40/a4_combined.tsv \
       --tables results/tables --suffix _v40
python scripts/audit_reproducibility.py        # 39/44 (the 5 are the
                                               # expected fingerprint
                                               # mismatches, section 6)
```
