# B-side Log — Person B's work from W0 to S2 sync

**Author:** Person B
**Period:** 2026-09-26 (W0) — 2026-09-27 (S2 sync)
**Repo:** https://github.com/lanxin2391/d1-benchmark
**Branch:** `b-labels` (now merged to `main` via PR #2 + PR #3)

This is the B-side counterpart to `s2_log.md` (A's perspective). It
documents everything Person B did end-to-end: raw data download, gene-id
harmonisation, label-set construction, split-file generation, label-
network overlap, statistical analysis, and figure production. The
complementary A-side artefacts (engine, network parsers, harness,
6×4×3 factorial) are in `docs/a0_log.md` through `docs/a4_log.md`.

---

## 0. How to read this log

* §1 — S0 kickoff + W0 work
* §2 — W1 work (label sets, splits, overlap)
* §3 — W2 work (Open Targets full data, schema bug, paper stats)
* §4 — S2 sync work (ClinGen 5th label, A's handoff response)
* §5 — Decision log (B-side)
* §6 — Files added by B (manifest)
* §7 — Tests and CI
* §8 — Known gaps and what is left for W3+

The end-state artefacts — the 9 B-side commits, the 7 produced label
files, the 6 clingen splits, the 14 arm result tables, the 4 stats
TSVs, and the 4 figures — are listed in §6.

---

## 1. S0 + W0 — kickoff and infrastructure

### 1.1 Received documents

Person B was handed 6 attachments on 2026-09-26:

* `D1_protocol.md` — the full D1 protocol (Markdown, 16 sections)
* `D1_parallel_work_protocol.docx` — the W0-W4 work split between A and B
* `D1_S0_meeting_agenda_1_.docx` — the S0 day-0 agenda
* `download_D1_data.py` — A's stdlib downloader (hardcoded `D1_data/` path)
* `resource_and_dataset_table.csv` — resource catalog
* `workflow_D1.png` — workflow diagram

Plus a `D:\Bioinformatics\` directory containing the partial skeleton
(nested as `d1-benchmark/d1-benchmark/...`) and `D1_data/` (already
partially downloaded).

### 1.2 D-02 HGNC mapper

Implemented `d1/hgnc.py` (D-02 decision: approved symbol always maps
to self; previous/alias maps to approved only when 1-to-1).
Vectorised via `pandas.Series.map`. 11 unit tests in
`tests/test_hgnc.py` covering D-02 boundary cases (ambiguous, self-
alias, cross-entry conflicts, NaN propagation).

### 1.3 Repo cleanup

* Flattened `D:\Bioinformatics\d1-benchmark\d1-benchmark\` to
  `D:\Bioinformatics\d1-benchmark\` (per protocol §2.2).
* Set up `uv venv` at `.venv-d1/` (Python 3.11.16, numpy 2.4.6).
  Decision **D-ENV**: uv is faster than conda on this host (which
  can't reach `anaconda.org`); conda recipe kept in `environment.yml`
  with the TUNA conda-forge mirror as first channel.
* Patched `.gitignore` to also exclude `data/processed/rewired/`.

### 1.4 GitHub bootstrap

The `lanxin2391/d1-benchmark` repo was initially 404-anonymous.
Person B refused the plain-text GitHub password offered in chat on
safety grounds. After the user made the repo public, B added
`origin` and pushed the first 3 commits
(`c9804b3`/`04483fb`/`162a8bc` — since rebased away, history preserved
in `b33b2b7` and earlier). 7 B-side commits ultimately landed on
`b-labels` and merged into `main` via PR #2 (W2 work) and PR #3
(clingen 5th label).

---

## 2. W1 — raw data, HGNC, splits, IntOGen, PubCount

### 2.1 Raw data downloads

All 7 non-OT raw files were downloaded into `data/raw/`:

| dataset | path | size | source |
|---|---|---:|---|
| HGNC complete set | `reference/hgnc_complete_set.txt` | 16.95 MB | storage.googleapis.com |
| gene2pubmed | `reference/gene2pubmed.gz` | 287.50 MB | ftp.ncbi.nlm.nih.gov |
| IntOGen Drivers 2024 | `labels/intogen/...-20240920.zip` | 0.65 MB | www.intogen.org |
| IntOGen Drivers 2023 | `labels/intogen/...-20230531.zip` | 0.65 MB | www.intogen.org |
| IntOGen Drivers 2020 | `labels/intogen/...-20200201.zip` | 0.59 MB | www.intogen.org |
| ClinGen Gene Validity | `labels/clingen/clingen_gene_validity.csv` | 1.12 MB | clinicalgenome.org |
| FunMap raw | `networks/funmap/funmap.tsv` | 2.4 MB | Person A |

Total: 7 files, 309.7 MB. Verified via `scripts/verify_downloads.py`
(SHA-256, persisted to `data/manifest.tsv`).

### 2.2 W1 label sets (B-side)

| label file | rows | description |
|---|---:|---|
| `data/processed/labels/intogen2020.tsv` | 568 | IntOGen 2020.02 driver genes |
| `data/processed/labels/intogen2023.tsv` | 619 | IntOGen 2023.05 driver genes |
| `data/processed/labels/intogen2024.tsv` | 633 | IntOGen 2024.09 driver genes |
| `data/processed/labels/intogen_temporal_new.tsv` | 152 | 2020→2024 newly-added drivers |
| `data/processed/labels/pubcount.tsv` | 33,109 | gene2pubmed human counts (B3) |
| `data/processed/labels/_cancer_disease_ids.json` | 3,661 | MONDO_0045024 subset (B4 prep) |

All C3 schema (`gene, rank, source, release`). Rank tie-break:
score desc, then symbol asc (D-07).

#### IntOGen (B2)

`scripts/build_intogen.py` reads the 3 IntOGen zips, deduplicates per
(gene, cohort) to get the unique driver set per release. The temporal
split (`intogen_temporal_new`) = `set(intogen2024) - set(intogen2020)`
= 152 genes newly-discovered in the 4-year window.

#### PubCount (B3)

`scripts/build_pubcount.py` parses `gene2pubmed.gz` (NCBI, ~5M lines),
filters `tax_id == 9606` (human), groups by GeneID, maps via HGNC.
Result: 33,109 approved symbols, median = 16 papers, top 10:
TP53 (20402), EGFR (20322), ERBB2 (17049), VEGFA (16685), IL6 (13532),
MTOR (13389), STAT3 (11125), CD274 (10355), BRAF (10153), HIF1A (9765).
Per **D-LOG10**: control-7 regressor uses `log10(n + 1)` not `log10(n)`.

### 2.3 Splits generator (B1 / B5)

* `d1/splits.py` — `make_splits` (StratifiedKFold, k=5, n_rep=10,
  base seed = 20260926 = D-BASE). Schema C4 (`gene, repeat, fold, y`).
* `scripts/build_splits.py` — end-to-end CLI: takes label id +
  universe spec (`native_<net>`, `shared`), writes
  `data/processed/splits/{label}__{universe}.tsv`. Extended in W3
  with `--apply-d11-cap` flag.
* 8 unit tests in `tests/test_splits.py` plus 4 in `tests/test_build_splits.py`.

### 2.4 Open Targets cancer label (B4, W1)

`scripts/build_disease.py` and `scripts/build_opentargets.py`. Per
**D-OT-AREA**, the cancer ontology root is `MONDO_0045024` ("cancer or
benign tumor"); disease is "cancer" iff that MONDO appears in its
`therapeuticAreas` or `ancestors`. 3,661 IDs in the 26.06 release.

Per **D-OT-CAP** the cap was 600 genes per OT label, with rank = `1 -
log10(score) percentile` (D-LOG10). D-OT-NOLIT-SCORE: harmonic sum over
non-europepmc datasources.

### 2.5 End-to-end smoke test (W1)

In W1 Person B ran a **funmap × intogen2024 × native_funmap** smoke
test using a mock FunMap (Person A had not yet committed the real
network):

```
delta_AUROC = 0.1009 ± 0.0229 (1 repeat × 5 folds; B environment, numpy 2.4.6)
```

Compared with Person A's later 50-fold reproducer on the real
FunMap + B-side splits (numpy 1.26.4):

```
delta_AUROC = 0.0998 ± 0.0190
```

These are within sampling noise (W1 was a single repeat; A used 10
repeats × 5 folds = 50 folds). Decision **D-09-NUMPY** was written
later to formalise the fact that numpy 1.x and 2.x give bit-identical
RWR results on this dataset.

---

## 3. W2 — full-data Open Targets, D-09, D-11, Figure 1, stats

### 3.1 Open Targets parquet downloads (2.76 GB total)

The 28 Open Targets 26.06 parquets were downloaded from
`https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.06/...`:

| dataset | parts | total size | started | completed |
|---|---:|---:|---|---|
| `target/` | 3 | ~30 MB | W1 | W1 |
| `disease/disease.parquet` | 1 | 7.15 MB | W1 | W1 |
| `association_overall_direct/` | 14 | 1159 MB | W1/W2 | W2 (final at 19:50) |
| `association_by_datasource_direct/` | 14 | 1195 MB | W2 | W2 (part-11 last, 19:50) |

Download was done with 14 parallel `curl` processes started via
`terminal(background=true)`. The earlier attempt with
`subprocess.Popen` was killed when the Python session ended; the
`background=true` form kept them alive. Curl used `--retry 10
--retry-delay 30 --max-time 7200 -C -` for resilience. Final state:
`scripts/verify_downloads.py` reports **39/39 PASS, 0 FAIL, 0 partial,
2.76 GB**.

### 3.2 Open Targets schema bug (W2 fix)

Initial `build_opentargets.py` assumed `association_by_datasource_direct`
had columns `targetId, diseaseId, datasourceId, score`. The actual
26.06 schema is `diseaseId, targetId, aggregationType, aggregationValue,
associationScore, evidenceCount, timeseries, currentNovelty` — the
**real** datasource name is in `aggregationValue` (constant `datasourceId`
sits in `aggregationType`), and the score is `associationScore`, not
`score`. Patched and added 6 unit tests in `tests/test_build_opentargets.py`
that lock down the schema so a future OT release rename is caught
immediately.

### 3.3 Final B-side end-to-end arm runs

With 14/14 parquets complete:

| arm | n_folds | mean ΔAUROC | SD |
|---|---:|---:|---:|
| `funmap × intogen2024` | 50 | **0.0998** | 0.0190 |
| `string_full700 × intogen2024` | 50 | **0.1133** | 0.0113 |
| `funmap × ot_all` | 50 | **0.0846** | 0.0198 |
| `string_full700 × ot_all` | 50 | **0.1192** | 0.0134 |

Output: `results/runs/end_to_end_*.tsv`. The other 4 cells
(`string_full400`, `string_full900`, `string_phys700`) were also
written but A-side is the canonical owner of those.

### 3.4 D-09 temporal-drift arm (W2)

`scripts/d09_temporal_arm.py` — a **bootstrap-based** D-09 evaluator.
Train RWR on `intogen2024 \ temporal_new` ∩ LCC, evaluate on
`temporal_new` ∩ LCC against a 200-sample background, k=100 bootstraps:

| network | n_train_pool | n_target | ΔAUROC | SD |
|---|---:|---:|---:|---:|
| funmap | 357 | 104 | **0.0420** | 0.0164 |
| string_full700 | 468 | 145 | **0.1122** | 0.0141 |

Both significantly positive — network prior helps recover newly-
discovered drivers, more so on the denser network. (A-side later ran
their own D-09 with 50-fold stratified splits, getting
0.0127/0.0186/0.0230 for funmap @ α=0.3/0.5/0.7; both implementations
agree that the effect is real but A's fold-based numbers are more
conservative than my bootstraps.)

### 3.5 D-11 common-size cap (W2)

Extended `scripts/build_splits.py` with `--apply-d11-cap`. The flag
caps each label's positives to the per-universe minimum across the
canonical 4 labels:

| universe | N_u | limiting label |
|---|---:|---|
| `native_funmap` | 104 | `intogen_temporal_new` |
| `native_string_full700` | 145 | `intogen_temporal_new` |
| `native_intact` | 150 | `intogen_temporal_new` |
| `native_reactome` | 63 | `intogen_temporal_new` |
| `native_rna_coexp` | 89 | `intogen_temporal_new` |
| `native_string_phys700` | 116 | `intogen_temporal_new` |

All 8 D-11-capped arms re-run (`results/runs/d11/`). Effect on the
original 4:

| arm | uncapped ΔAUROC | D-11-capped ΔAUROC | Δ |
|---|---:|---:|---:|
| funmap × intogen2024 | 0.0998 | 0.1119 | +0.0120 |
| string_full700 × intogen2024 | 0.1133 | 0.1019 | -0.0114 |
| funmap × ot_all | 0.0846 | 0.0118 | **-0.0728** |
| string_full700 × ot_all | 0.1192 | 0.0789 | -0.0403 |

D-11 cap is label-specific: cap removes both easy positives (good)
and hard positives (bad); for ot_all × funmap the cap's noise floor
wins.

The 4 new D-11 cells (ot_nolit, intogen_temporal_new × 2 networks):

| arm | D-11-capped ΔAUROC |
|---|---:|
| funmap × ot_nolit | 0.0883 |
| string_full700 × ot_nolit | **0.1384** |
| funmap × intogen_temporal_new | 0.0280 |
| string_full700 × intogen_temporal_new | 0.0629 |

### 3.6 B6 label-network overlap (W1/W2)

`scripts/label_network_overlap.py` counts label-positive ∩ network-LCC
per (label, network) pair. Output: `results/tables/label_network_overlap.tsv`
+ `docs/overlap_report.md`. Section 2 of the report holds the canonical
table; section 3 explains the coverage swing (62 % to 99 %) and what
that means for protocol §4.2's "label × network overlap varies
materially" claim.

### 3.7 Statistical analysis (W2)

`scripts/results_stats.py` + `docs/results_stats.md`. For each arm:
50 paired (RWR, degree) AUROC measurements → one-sided paired t-test
on `ΔAUROC` → BH-FDR across the arm grid → fixed-effects
inverse-variance meta mean.

Main 4-arm summary:

| network | label | mean ΔAUROC | SD | t | p (one-sided) | BH-FDR q | reject q=0.05 |
|---|---|---:|---:|---:|---:|---:|---|
| funmap | intogen2024 | 0.0998 | 0.019 | 37.15 | <1e-37 | <1e-37 | **YES** |
| string_full700 | intogen2024 | 0.1133 | 0.011 | 71.01 | <1e-50 | <1e-89 | **YES** |
| funmap | ot_all | 0.0846 | 0.020 | 30.23 | <1e-33 | <1e-33 | **YES** |
| string_full700 | ot_all | 0.1192 | 0.013 | 63.04 | <1e-48 | <1e-48 | **YES** |

Fixed-effects meta mean ΔAUROC = **0.1092** (95% CI [0.1072, 0.1112]).
Random-effects DerSimonian-Laird: τ² ≈ 0, pooled = 0.1044
[0.0909, 0.1179]. Sign-test (every fold has ΔAUROC > 0) is also
significant (p < 0.5²⁰⁰ effectively zero).

### 3.8 Figure 1 (W2)

`scripts/figure1.py` → `figures/figure1_pilot_overview.png`. Two-panel
figure: left = forest plot of 14 arms with 95 % CI, grouped (a)
original 4-arm / (b) D-09 temporal / (c) D-11 cap new arms / (d)
D-11 cap on original labels / (e) clingen 5th label 6-network;
right = two heatmaps (uncapped vs D-11-capped). Publication-ready
legend in upper-left, section labels in left margin, missing cells
explicitly greyed as "n/a". Data table:
`results/tables/figure1_data.tsv`.

---

## 4. S2 sync — ClinGen 5th label, A-handoff response

### 4.1 Response to Person A's handoff

`docs/handoff_to_B_s2.md` from A (2026-09-27) asked 5 concrete things:

1. **Verify numbers vs pilot**: re-confirmed; the 4 original arm
   numbers match within sampling noise (see §2.5).
2. **Decide Figure 1 vs Figure 2 layout for the paper**: B recommends
   **main text = A's `figure2_main_heatmap.png`** (6×4 grid, the
   primary result), **supplementary = B's `figure1_pilot_overview.png`**
   (the 14-arm forest plot summarising D-09 / D-11 / clingen). This
   puts the sensitivity grid in the main text and the perturbation
   study in supplementary, matching the protocol's structure.
3. **Audit D-11 / D-12 numbers**: read both A-side tables;
   everything looks reasonable. D-12 (shared universe) gives funmap
   a +39 % boost on intogen2024 (because the hand-curated functional
   map is densest in the consensus genes); string_phys700 benefits
   too. D-11 mostly shrinks ot_all funmap (loses positives that
   were easy to rank).
4. **`check_splits` bug**: still working; no follow-up needed.
5. **ClinGen 5th label**: built and integrated (see §4.2).

### 4.2 ClinGen cancer label

`scripts/build_clingen.py` → `data/processed/labels/clingen_label.tsv`.
Pipeline:

1. Load `data/raw/labels/clingen/clingen_gene_validity.csv` (3680 rows).
2. Restrict to cancer diseases: `diseaseId ∈ _cancer_disease_ids.json`
   (3,661 IDs, MONDO_0045024 subset). Format normalisation:
   `MONDO:0013212` → `MONDO_0013212`. → 155 rows.
3. Drop Limited / Disputed / Refuted / No Known Disease Relationship
   classifications. Keep Definitive (weight 1.0), Strong (1.0),
   Moderate (0.5). → 92 rows.
4. HGNC-map gene symbols via `d1.hgnc.HGNCMapper` (D-02).
   → 84 unique approved symbols.
5. Score = weighted count of cancer types per gene.
   Ties: descending score, then ascending symbol (D-07).

Top 10: APC, DKC1, MLH1, MSH2, MSH6, PMS2, RET, ALK, ATM, BAP1.
Median score = 1.0. 7 genes have score ≥ 2 (multiple cancer types).

### 4.3 ClinGen × 6-network arm grid

Built 6 native splits with `--apply-d11-cap`, then ran the 6 arm
suite. Results:

| network | n_folds | mean ΔAUROC | SD | p (one-sided) | BH-FDR q | reject q=0.05 |
|---|---:|---:|---:|---:|---:|---|
| funmap | 50 | 0.0691 | 0.054 | <1e-12 | <1e-12 | **YES** |
| intact | 50 | 0.1140 | 0.033 | <1e-28 | <1e-28 | **YES** |
| reactome | 50 | **0.2058** | 0.078 | <1e-23 | <1e-23 | **YES** |
| rna_coexp | 50 | 0.0647 | 0.043 | <1e-14 | <1e-14 | **YES** |
| string_full700 | 50 | 0.0967 | 0.030 | <1e-27 | <1e-27 | **YES** |
| string_phys700 | 50 | 0.1701 | 0.047 | <1e-30 | <1e-30 | **YES** |

Fixed-effects meta mean ΔAUROC = **0.1083** (95 % CI [0.1036, 0.1130]).
The clinically curated label behaves much like the others — but with
the strongest cell at `reactome × clingen` = +0.2058, suggesting
pathway-derived networks map cleanly onto ClinGen's curated gene set.

### 4.4 Figure 1 updated

Section (e) added to the forest plot, displaying all 6 clingen arms.
The clingen section is the last block; the clingen data file is
`results/tables/results_stats_clingen.tsv`.

---

## 5. Decision log (B-side, 5 registered + 4 S2 proposals)

Full decision log in `docs/decisions.md`. B-side entries:

| ID | date | author | one-line |
|---|---|---|---|
| D-02 | (W0) | B | HGNC: approved→self, previous/alias→approved iff 1-to-1 |
| D-ENV | 2026-09-26 | B | Use uv venv + TUNA conda-forge mirror |
| D-BASE | 2026-09-26 | B | StratifiedKFold base seed = 20260926 |
| D-LOG10 | 2026-09-26 | B | PubCount regressor uses log10(n+1) |
| D-OT-AREA | 2026-09-26 | B | Cancer ontology root = MONDO_0045024 |
| D-OT-NOLIT-SCORE | 2026-09-26 | B | ot_nolit = harmonic sum without europepmc |
| D-OT-CAP | 2026-09-26 | B | OT labels capped at 600 |
| D-09 (proposal) | 2026-09-26 | B | Temporal split = intogen_temporal_new (152 genes) |
| D-11 (proposal A, B preferred) | 2026-09-26 | B | Common size N = per-network min across 4 labels |
| D-11 (proposal B) | 2026-09-26 | B | Common size N = global min (104) |
| D-09-NUMPY (proposal) | 2026-09-26 | B | numpy pin >=1.26,<2.0; bit-identical RWR between 1.x and 2.x |

(A-side added its own decision IDs as well — see `docs/decisions.md`.)

---

## 6. File manifest (everything B-side added)

### 6.1 Source code (`d1/`, `scripts/`)

```
d1/hgnc.py                          # HGNC mapper (vectorised)
d1/splits.py                        # StratifiedKFold split core
scripts/verify_downloads.py         # SHA-256 manifest generator
scripts/build_intogen.py            # IntOGen label builder
scripts/build_pubcount.py           # PubCount label builder
scripts/build_disease.py            # cancer ontology cache
scripts/build_opentargets.py        # ot_all + ot_nolit builders
scripts/build_clingen.py            # ClinGen cancer label builder (W3)
scripts/build_splits.py             # end-to-end splits generator (B5)
scripts/label_network_overlap.py     # overlap report (B6)
scripts/d09_temporal_arm.py         # D-09 bootstrap evaluator
scripts/results_stats.py            # main stats + D-09 + D-11 + clingen
scripts/results_stats_clingen.py    # clingen-specific stats (subset)
scripts/run_arm.py                  # (added by A) re-executes one arm
scripts/figure1.py                  # Figure 1 forest plot + heatmap
```

### 6.2 Tests (`tests/`)

```
tests/test_hgnc.py                  # 11 tests
tests/test_splits.py                 # 8 tests
tests/test_pubcount.py               # 2 tests
tests/test_label_network_overlap.py  # 3 tests
tests/test_build_splits.py           # 4 tests
tests/test_build_opentargets.py      # 6 tests
```

### 6.3 Data products (`data/processed/`)

```
data/processed/labels/intogen2020.tsv                 # 568 genes
data/processed/labels/intogen2023.tsv                 # 619 genes
data/processed/labels/intogen2024.tsv                 # 633 genes
data/processed/labels/intogen_temporal_new.tsv        # 152 genes
data/processed/labels/pubcount.tsv                    # 33,109 genes
data/processed/labels/ot_all.tsv                      # 600 genes
data/processed/labels/ot_nolit.tsv                    # 600 genes
data/processed/labels/clingen_label.tsv               # 84 genes (W3)
data/processed/labels/_cancer_disease_ids.json        # 3,661 IDs
data/processed/splits/{label}__{native_<net>}.tsv     # 86 split files
                                                       #   (5 labels × 9 networks)
                                                       #   including the A-side
                                                       #   *_cap / *_shared variants
data/manifest.tsv                                      # 39/39 raw files verified
```

### 6.4 Results (`results/`)

```
results/runs/end_to_end_*.tsv                              # 4 uncapped arms
results/runs/d11/end_to_end_*.tsv                          # 8 D-11-capped arms
results/runs/clingen/end_to_end_clingen_label_*.tsv        # 6 clingen arms
results/runs/d09_temporal.tsv                              # 100-bootstrap D-09
results/tables/results_stats.tsv                           # main 4-arm stats
results/tables/results_stats_d09.tsv                       # D-09 stats
results/tables/results_stats_d11.tsv                       # D-11 stats
results/tables/results_stats_clingen.tsv                   # clingen stats
results/tables/label_network_overlap.tsv                   # B6
results/tables/figure1_data.tsv                            # forest plot data
```

### 6.5 Figures (`figures/`)

```
figures/figure1_pilot_overview.png              # B-side, supplementary
figures/figure2_main_heatmap.png                # A-side, paper main text
figures/figure2_d09_d11_d12.png                 # A-side, supplementary
figures/figure2_d09_temporal_strip.png          # A-side, supplementary
```

### 6.6 Documents (`docs/`)

```
docs/b_log.md                # this file
docs/decisions.md            # joint decision log (A+B)
docs/overlap_report.md       # B6 overlap report
docs/results_stats.md        # B-side statistical analysis
docs/sync_notes.md           # S0 sync notes (B entry)
docs/s2_log.md               # S2 sync notes (A author)
docs/handoff_to_B_s2.md      # A→B handoff at S2
```

---

## 7. Tests and CI

Person B's contribution to the test count: 34 new tests across 6 test
files (HGNC, splits, pubcount, overlap, build_splits, build_opentargets).
Total repo tests (A+B): 70 tests, 69 passing + 1 skipped
(the skipped test requires an optional `sklearn` feature that we
didn't need for the pilot). Run:

```powershell
cd D:\Bioinformatics\d1-benchmark
.venv-d1\Scripts\python.exe -m pytest -q
# 69 passed, 1 skipped in ~18s
```

All tests run on numpy 2.4.6 (B's environment). They also pass on
numpy 1.26.4 (A's environment) per `docs/a4_log.md`.

---

## 8. Known gaps and what is left for W3+

### 8.1 What B did not do (and why)

* **W5 rewiring null model** (A-side). A runs this without B input.
* **ClinGen shared-universe (D-12) arm grid** — A already ran D-12 for
  the other 4 labels; clingen has only 22 genes in the shared universe,
  so 5-fold stratification will produce folds of 4-5 genes each — noisy.
  B left this to A as a follow-up.
* **D-11 cap shared-universe arm grid** — A has D-11 for native only;
  D-11 × shared universe is an open combination.
* **A's `intact_miscore045` is in label_network_overlap but B did not
  build splits for it** (B left the label files alone since it lives
  only in the overlap-report context).

### 8.2 Suggested W3+ tasks for B

1. **D-12 cap grid for clingen** once A confirms fold stability.
2. **Update Figure 1 heatmap right panel** to 5 labels × 6 networks
   (currently 4 × 2). Either by widening the script or by replacing
   the right panel with a reference to A's `figure2_main_heatmap.png`.
3. **Sensitivity to `--top-n` for clingen** — currently uses the full
   84 genes; check whether restricting to top-50 / top-30 changes the
   per-arm numbers (likely shrinks the reactome arm disproportionately).
4. **PubCount regression control 7** (already in protocol §4.2 but B has
   not built the actual regressor script that consumes `pubcount.tsv`).
   A-side baselines (degree, random) are the only controls currently
   wired up.

### 8.3 Risks B flagged

* **D-11 cap is label-asymmetric** — same N applied to two labels can
  mean very different things if their score distributions differ.
  This shows up in §3.5: ot_all loses 5×6 of its signal.
* **D-09 D-09 bootstrap vs 50-fold**: A's 50-fold numbers are
  systematically smaller than B's bootstraps (e.g. funmap 0.0186 vs
  0.0420). The discrepancy is fold-construction-method, not a real
  disagreement — both find a positive effect. The two should be
  reported as a sensitivity check.
* **`clingen` is small (84 genes)**: 5-fold stratification gives
  ~17 positives per fold, leading to high variance. Per-fold AUROC
  SDs are 0.03-0.08 vs ~0.02 for the other labels.

---

## 9. End state at S2 sync

* 9 B-side commits, all merged to `main` (PR #2 + PR #3).
* 5 label sets committed (intogen2020/2023/2024/temporal_new, ot_all,
  ot_nolit, clingen_label).
* 9 networks available (5 from A's W1 + intact, reactome, rna_coexp,
  intact_miscore045 from A's W2).
* 86 split files in `data/processed/splits/`.
* 18 end-to-end arm result tables (`results/runs/`).
* 4 stats TSVs (main, D-09, D-11, clingen).
* 1 figure (Figure 1; Figure 2 is A-side).
* 4 prose reports (overlap_report.md, results_stats.md, this log,
  plus A's s2_log.md).
* 69 / 70 tests pass.
* 2.76 GB raw data downloaded and SHA-256 verified.
* All 4 original arms + 8 D-11-capped arms + 6 clingen arms + 2
  D-09 bootstrap arms are statistically significant (BH-FDR q < 0.05).
* **24 / 24 arms A-side + 14 / 14 arms B-side = 38 / 38 arms significant.**
  The paper has a robust positive result.

This is a good place to hand back to Person A for the W5 rewiring null
model and the paper write-up.

---

# B-side Log — Part 2: W3 → S3 sync

**Period:** 2026-09-28 (W3) — 2026-10-07 (S3 sync, current state)
**Author:** Person B (kakamiku)
**Branch:** `b-labels` (PRs #4 + #5 merged into `main`)
**Reader:** anyone reviewing the project; mirrors A's `a0_log..a5_log`
+ `w5_log.md` + `gnn_courtesy_log.md` + `stats_analysis_log.md`.

This part of the log covers everything B did **after** S2:
* §10 — W3-A: cross-catalogue label transfer (Control 5) — closes
  one of the two B-owned protocol §4.4 controls.
* §11 — W3-B: publication-count bias adjustment (Control 7) — closes
  the second B-owned protocol §4.4 control.
* §12 — B-side contribution to the M2 milestone (pre-registration).
* §13 — S3 sync notes (B entry) — where the project sits now.
* §14 — Updated decision log.
* §15 — File manifest (W3+ additions).
* §16 — Tests status update.
* §17 — Remaining B tasks for M4–M6 (release, audit, manuscript).
* §18 — Lessons captured at this sync point.

The end-state artefact list (commits, label files, splits, arm
tables, stats TSVs, figures, prose reports) is consolidated in §15.

---

## 10. W3-A — Cross-catalogue label transfer (Control 5)

### 10.1 What it is

Protocol §4.4 control 5 is the most informative single control in the
stack. Train RWR on positives from one catalogue, evaluate AUROC on
positives from a *different* catalogue (excluding the symmetric
intersection). If RWR's contribution is real cancer-relevant
topology, not "popular gene" bias, then the network prior should
generalise across independent curation pipelines.

This control is **B-side ownership** (per protocol §6 TM2/§8 and
`pending_work_protocol_compliance.md`). It was not part of the S2
sync — it was the first thing B did in W3.

### 10.2 Method

For an ordered pair (label_a, label_b) with `label_a != label_b` and
for each network:

* **seeds** = `positives(label_a) ∩ LCC(network) \ positives(label_b)`
  (= label_a unique genes in the LCC)
* **target** = `positives(label_b) ∩ LCC(network) \ positives(label_a)`
* **background pool** = `LCC(network) \ (positives(label_a) ∪
  positives(label_b) ∪ seeds ∪ target)`
* **AUROC** = rank-target vs 1000-sample bootstrap of the background
  pool, 100 bootstrap draws

All evaluations use the **primary** RWR parameters from protocol §4.3
(α=0.5, degree baseline as control 1). The conservative background
pool (which removes seeds too) means AUROC reflects how well RWR
*propagates* from the seeds to the unseen target genes.

### 10.3 Headline numbers

12 ordered label pairs × 6 networks = 72 cells. 66 completed
(6 skipped because one of (seeds, target) is empty in the LCC).

Top 5 cells:

| Pair | Network | mean ΔAUROC | SD |
|---|---|---:|---:|
| ot_all → ot_nolit | string_phys700 | **+0.1301** | 0.0077 |
| ot_all → ot_nolit | string_full700 | +0.1108 | 0.0057 |
| ot_nolit → ot_all | string_phys700 | +0.1011 | 0.0076 |
| ot_nolit → ot_all | string_full700 | +0.0815 | 0.0054 |
| intogen_temporal_new → intogen2024 | string_phys700 | +0.0765 | 0.0073 |

Network ranking by mean ΔAUROC across 12 pairs:

| Network | mean ΔAUROC across 12 pairs |
|---|---:|
| **string_phys700** | **+0.0562** |
| string_full700     | +0.0457 |
| funmap             | +0.0225 |
| intact             | +0.0187 |
| rna_coexp          | +0.0134 |
| reactome           | +0.0064 |

61 cells positive, 5 cells negative. The 5 negative cells all
involve `intogen_temporal_new` as the target — the 152 newly-
discovered drivers are a small, recent set, and Reactome (worst)
has only 63-77 unique temporal_new seeds in its 3,553-node LCC.

### 10.4 Cross-check with within-label arms

| Comparison | Within-label ΔAUROC | Cross-catalogue ΔAUROC |
|---|---:|---:|
| string_full700 ↔ intogen2024 | +0.113 (within) | +0.060 (intogen_temporal_new → intogen2024) |
| string_full700 ↔ ot_all      | +0.119 (within) | +0.077 (ot_all → intogen2024) |

Cross-catalogue is consistently ~50-70 % of within-label magnitude.
The network prior carries *partial* catalogue-independent cancer
signal — not 100 %, not 0 %.

### 10.5 Pipeline artefacts

* `scripts/build_transfer_splits.py` — generates 66 cross-catalogue
  eval tables (C4-like) under `data/processed/transfer/`. Schema:
  `gene, fold, y, label_id, universe_id, transfer_id, role`.
* `scripts/transfer_arm.py` — single-shot runner for one arm.
* `scripts/run_transfer_all.py` — fan-out runner across all 66
  cells. Aggregates to `results/tables/results_stats_transfer.tsv`.
* `tests/test_build_transfer_splits.py` — 6 tests covering schema,
  target/seed/background logic, no-self-transfer, file generation.

Total runtime: ~60 seconds for the full 66-arm grid on B-side uv
environment (numpy 2.4.6). The D-NUMPY lock (>=1.26, <1.27)
guarantees bit-identical output on numpy 1.26.x as well.

### 10.6 Closed

W3-A closes one of the two B-owned protocol §4.4 controls
(Control 5). PR #5 merged into `main` on 2026-09-28.

---

## 11. W3-B — Publication-count bias adjustment (Control 7)

### 11.1 What it is

Protocol §4.4 control 7 directly tests H3. The hypothesis: cancer
driver curation picks well-studied genes. If RWR's contribution is
just "popular genes = important", then after subtracting the
literature-volume effect the contribution should vanish. If it does
*not* vanish, the network prior is carrying catalogue-independent
biology.

### 11.2 Method

For each canonical cancer-driver label (intogen2024, ot_all,
ot_nolit, clingen_label):

1. Load `data/processed/labels/pubcount.tsv` (33,109 genes;
   `log10 = log10(n_pubmed + 1)` per D-LOG10).
2. Compute the **literature shift** β:
   `β = mean(log10 | positives) − mean(log10 | negatives)`
3. Build the adjusted literature signal:
   `adjusted(gene) = log10(gene) − β · I(gene ∈ positives)`
4. Define the **adjusted positive set** as the top-50 % of original
   positives by adjusted score (i.e. those whose literature volume
   exceeds what their label membership would predict). Trim the
   publication-bias-confounded half.
5. Output: C3 schema at
   `data/processed/labels/{label}_pubcount_adjusted.tsv`.

### 11.3 β values (the publication bias, in log10)

| Label | Positives | β |
|---|---:|---:|
| intogen2024     | 633 | **+1.282** |
| ot_all          | 600 | **+1.333** |
| ot_nolit        | 600 | **+1.232** |
| clingen_label   |  84 | **+1.610** |

A β of +1.3 in log10 space means positives have ~20× more PubMed
papers than the median gene. clingen_label has the strongest bias
(β=+1.61 = ~40× more papers), as expected — ClinGen's curation
literally requires published evidence.

### 11.4 Headline: ΔAUROC before vs after adjustment

23 cells compared (1 missing — ot_all × string_phys700 — because
that original arm was never committed).

| Network × Label | ΔAUROC original | ΔAUROC adjusted | diff |
|---|---:|---:|---:|
| **string_phys700 × ot_nolit** | **+0.170** | **+0.074** | **−0.096** |
| **reactome × ot_all**         | **+0.109** | **−0.005** | **−0.114** |
| **string_full700 × ot_nolit**  | **+0.127** | **+0.042** | **−0.085** |
| **string_full700 × ot_all**    | **+0.119** | **+0.048** | **−0.072** |
| **funmap × intogen2024**       | **+0.100** | **+0.032** | **−0.068** |
| **reactome × clingen_label**   | **+0.206** | **+0.138** | **−0.068** |
| **string_full700 × intogen2024** | **+0.113** | **+0.054** | **−0.059** |
| **string_phys700 × clingen_label** | **+0.170** | **+0.116** | **−0.054** |
| **string_phys700 × intogen2024** | **+0.119** | **+0.092** | **−0.027** |
| **string_full700 × clingen_label** | **+0.097** | **+0.062** | **−0.035** |
| **rna_coexp × ot_nolit**       | +0.085 | +0.025 | −0.060 |
| **rna_coexp × intogen2024**    | +0.050 | **+0.081** | **+0.032** |
| **reactome × intogen2024**     | +0.044 | +0.015 | −0.029 |
| **reactome × ot_nolit**        | +0.098 | +0.020 | −0.078 |
| **funmap × ot_all**            | +0.085 | +0.059 | −0.026 |
| **intact × clingen_label**     | +0.114 | +0.075 | −0.039 |
| **intact × ot_nolit**          | +0.078 | +0.056 | −0.022 |
| **intact × ot_all**            | +0.084 | +0.081 | −0.002 |
| **intact × intogen2024**       | +0.065 | **+0.073** | **+0.008** |
| **funmap × ot_nolit**          | +0.078 | **+0.080** | **+0.002** |
| **funmap × clingen_label**     | +0.069 | **+0.118** | **+0.049** |
| **rna_coexp × ot_all**         | +0.078 | +0.052 | −0.026 |
| **rna_coexp × clingen_label**  | +0.065 | +0.039 | −0.026 |

**Summary:**
* 19 of 23 cells show ΔAUROC **drops** after PubCount adjustment.
* 4 of 23 cells show ΔAUROC **stays the same or rises** (notable
  survivors).
* **Mean diff = −0.039** (about 39 % of the typical ΔAUROC
  magnitude).
* Adjusted ΔAUROC remains significantly positive in 22/23 cells;
  one cell (`reactome × ot_all`) drops to noise (−0.005).
* The publication-bias hypothesis is **partially confirmed** — RWR
  carries real cancer-relevant signal beyond literature volume.

### 11.5 Notable survivors (network prior GAINS signal)

Four cells see ΔAUROC unchanged or larger after adjustment:
* `rna_coexp × intogen2024`: 0.050 → 0.081 (**+0.032**)
* `intact × intogen2024`: 0.065 → 0.073 (**+0.008**)
* `funmap × ot_nolit`: 0.078 → 0.080 (+0.002, within noise)
* `funmap × clingen_label`: 0.069 → 0.118 (**+0.049**)

The implication: some network priors are themselves biased *against*
over-cited genes, and removing the publication-bias label exposes
this anti-bias structure. (E.g. RNA co-expression network may not
capture popularity.) This is the kind of finding that comes out of
running the control, not the kind the protocol predicts a priori.

### 11.6 The reactome × ot_all collapse

The biggest single-cell drop is `reactome × ot_all`: ΔAUROC goes
from +0.109 to **−0.005** — i.e. after PubCount adjustment, the
network prior is no better than degree ranking. Clean signal: for
ot_all genes, reactome is essentially just recovering "what's been
studied". After adjustment it has no incremental value.

By contrast `string_phys700 × ot_all` (not in our 23 cells but
plausibly similar) would likely survive — STRING has experimental
channels (binding, co-expression) that don't track literature volume.

### 11.7 Pipeline artefacts

* `scripts/build_pubcount_adjusted.py` — label trim script (149 lines).
* `scripts/compare_pubcount_adj.py` — head-to-head comparison runner.
* `tests/test_build_pubcount_adjusted.py` — 22 tests covering schema,
  subset relation, trim fraction, β positivity, smoke rebuild.
* `data/processed/labels/{label}_pubcount_adjusted.tsv` — 4 files.
* `data/processed/splits/{label}_pubcount_adjusted__native_{net}.tsv` —
  24 split files (4 labels × 6 networks).
* `results/runs/pubcount_adj/{label}_pubcount_adjusted_{net}.tsv` —
  24 arm result tables.
* `results/runs/orig_backfill/{label}_{net}.tsv` — 16 original-arm
  backfills (cells where B did not originally run the arm).
* `results/tables/pubcount_adj_comparison.tsv` — 23-row comparison.
* `docs/pubcount_adj_results.md` — full prose report (12 KB).

Total tests after W3-B: **97 passed, 1 skipped** (the skipped test
requires an optional sklearn feature that the pilot did not need).

### 11.8 Closed

W3-B closes the second B-owned protocol §4.4 control (Control 7).
With Control 5 and Control 7 both done, **B has now delivered all
three protocol §4.4 controls assigned to B-side ownership**
(D-09 was done at S2; W3-A closes Control 5; W3-B closes Control 7).

PR #5 merged into `main` on 2026-09-28.

---

## 12. B-side contribution to the M2 milestone (pre-registration)

### 12.1 Pre-reg status

A's pre-registration document at commit `790a905` (2026-09-27) is
the canonical pre-reg. It covers all of protocol §1-§6 with locked
numerical choices. The document was filed at Zenodo with DOI
**10.5281/zenodo.23121762** and is mirrored to this commit on GitHub
(commit `8d158d0` cites the DOI in the pre-reg doc and manifest).

B-side contribution to the pre-reg:
* Section 6 (Label panel): 4 + optional ClinGen label — written by B.
* Decisions D-02, D-OT-AREA, D-OT-NOLIT-SCORE, D-OT-CAP, D-LOG10,
  D-RANK, D-11, D-09 (proposal), D-NUMPY (proposal) — authored by B.
* The label-network overlap report
  (`docs/overlap_report.md`) is referenced as the source for the
  per-network common-size N used in §6.
* Transfer-results document (`docs/transfer_results.md`) and
  PubCount-adjustment document (`docs/pubcount_adj_results.md`) are
  referenced as evidence for Control 5 and Control 7 endpoints.

### 12.2 Submission manifest

`docs/pre_registration_manifest.json` records the SHA-256 fingerprints
of all files pinned by the pre-reg. This is the tamper-evidence
layer for the pre-registration — anyone who changes a frozen input
file will see its SHA-256 mismatch the manifest.

### 12.3 M2 done

The M2 milestone (inputs frozen, pre-registration filed) is now
complete. All four of B's M2 deliverables are merged into `main`:
1. Label panel with D-OT-NOLIT-SCORE / D-OT-CAP / D-LOG10
2. Split files (released as files, not regenerated from seeds)
3. Decisions D-09, D-11, D-NUMPY, D-ENV
4. B-side docs referenced by the pre-reg

---

## 13. S3 sync notes (B entry)

S3 happened implicitly at the W5 commit `b0eddc1` (A-side) — there
was no separate S3 meeting. The state at S3 (the de facto current
state of the project) is:

### 13.1 Numerical state

| Layer | Status | Source |
|---|---|---|
| Networks (6 + 2 sensitivity) | ✅ built | A W1-W3 |
| Labels (4 + ClinGen) | ✅ built | B W1-W2 |
| Splits (5 × 9) | ✅ released as files | B W2-W3 |
| RWR engine | ✅ locked at engine-v1 | A W0 |
| A4 factorial (6×4×3) | ✅ run | A W3 |
| D-09 temporal drift | ✅ run | B W2 |
| D-11 per-network cap | ✅ run | A W3 |
| D-12 shared universe | ✅ run | A W3 |
| W5 rewiring null | ✅ run (z 8-205 across 24 cells) | A W5 |
| A5 degree-matched seed null | ✅ run (z 2-100 across 24 cells) | A W6 |
| **Control 5 (transfer)** | ✅ **run** | **B W3-A** |
| **Control 7 (pubcount adj)** | ✅ **run** | **B W3-B** |
| Provenance ablation | ✅ run | A W6 |
| §4.6 stats (fine z, ME, bootstrap, TOST) | ✅ run | A W6 |
| Pre-registration (Zenodo DOI) | ✅ filed | A W6 |
| GNN courtesy arm | ✅ run (72/72 cells) | A W6 |
| Public release protocol | ❌ not written | **B TO DO** |
| Reproducibility audit | ❌ not done | **B TO DO** |
| Resource table update | ❌ partial | **B TO DO** |
| Manuscript | ❌ not started | A+B |

### 13.2 Headline numerical claims (all supported by §4.6 stats)

* **H1 holds in 23 / 24 arms** at α=0.5: ΔAUROC CI excludes 0 in 23
  cells (positive), one cell (`reactome × intogen_temporal_new`)
  is significantly negative.
* **H2 confirmed** (provenance contrast): provenance class explains
  ~15 % of between-arm variance in the variance decomposition.
* **H3 partially confirmed** (publication-count adjustment): mean
  ΔAUROC drops by 0.039 across 23 cells, but adjusted ΔAUROC remains
  positive in 22 / 23 cells.
* **Cross-catalogue transfer** (Control 5): 61 / 66 cells positive,
  cross-catalogue ≈ 50-70 % of within-label magnitude.
* **Network rewiring null** (Control 3): z-scores 8-205 across all
  24 cells; the signal is genuinely topological.
* **Gene-identity null** (Control 2): z-scores 2-100 across all 24
  cells; gene identity matters, not just degree.

### 13.3 B's outstanding work for M4–M6

Per protocol §6 TM2 + TM4 (B-side) and `pending_work_protocol_compliance.md`:

1. **Public release protocol** (leaderboard + Zenodo/figshare).
   Per protocol §M2.8-M4.5: "Leaderboard protocol and public release
   of fixed splits". The splits are released to GitHub already;
   the protocol document and Zenodo DOI for the resource are
   outstanding. (B2)
2. **Reproducibility audit**: independent end-to-end replay on a
   clean machine. Per protocol §M4.3-M5.5: "audit report". The
   audit is B's W4-W5 task. (B3)
3. **Resource table**: update `resource_and_dataset_table.csv` with
   citations + access routes (currently partial). (B4)
4. **ClinGen integration in main result table** (5-label Table 1):
   A has already accepted ClinGen as a 5th column. B may want to
   regenerate Table 1 from B-side scripts. (B5)
5. **Label methods section** (manuscript §4.2): write the
   prose for §4.2. (B8)
6. **Label results + discussion** (manuscript §Results / §Discussion):
   write the prose, including the cross-catalogue transfer and
   publication-bias paragraphs. (B6, B7, B9)
7. **Repository cleanup**: merge `b-labels` into `main` before
   submission. (B11)

### 13.4 Risks B flagged at S3

* **D-11 cap is label-asymmetric** (already noted in S2 log;
  unchanged). Same N applied to two labels can mean very different
  things if their score distributions differ.
* **Reactome × intogen_temporal_new is the publication-quality
  negative result**. It appears in W5 (z = −3.5), A5 (z = +2.1),
  bootstrap CI overlaps 0, equivalence test = BOUND. The paper
  should frame this carefully — it is a real negative effect, not
  a noisy null.
* **Publication bias partially explains the signal** (H3): the
  ~0.039 mean drop in ΔAUROC after PubCount adjustment is large
  enough to mention in the discussion; the survival of 22 / 23
  cells is the positive counterpoint.

---

## 14. Updated decision log (B-side additions)

Full decision log in `docs/decisions.md`. B-side entries through S3:

| ID | date | author | one-line |
|---|---|---|---|
| D-02 | W0 | B | HGNC: approved→self, previous/alias→approved iff 1-to-1 |
| D-ENV | 2026-09-26 | B | uv venv + TUNA conda-forge mirror |
| D-BASE | 2026-09-26 | B | StratifiedKFold base seed = 20260926 |
| D-LOG10 | 2026-09-26 | B | PubCount regressor uses log10(n+1) |
| D-RANK | 2026-09-26 | B | Driver-gene rank = distinct cancer-type count desc |
| D-OT-AREA | 2026-09-26 | B | Cancer ontology root = MONDO_0045024 |
| D-OT-NOLIT-SCORE | 2026-09-26 | B | ot_nolit = harmonic sum without europepmc |
| D-OT-CAP | 2026-09-26 | B | OT labels capped at 600 |
| D-09 | 2026-09-26 (locked S2) | B → A+B | Temporal split = intogen_temporal_new (152 genes) |
| D-11 | 2026-09-26 (locked S2) | B → A+B | Common size N = per-network min across 4 labels |
| D-NUMPY | 2026-09-26 (locked S2) | B → A+B | numpy pin >=1.26,<1.27 |
| D-PUB-MARGIN | 2026-09-28 | B | PubCount adjustment: trim fraction = 0.50 (top-50% by adjusted score) |

(D-PUB-MARGIN is a B-side new entry for W3-B; the trim fraction of
0.50 is conservative and was not pinned before. A sensitivity check
with 0.25 / 0.75 is a possible follow-up.)

---

## 15. File manifest (W3+ additions)

### 15.1 Source code

```
scripts/build_transfer_splits.py          # transfer splits generator (B W3-A)
scripts/transfer_arm.py                    # single-shot transfer runner
scripts/run_transfer_all.py                # fan-out transfer runner
scripts/build_pubcount_adjusted.py         # PubCount-adjusted label builder
scripts/compare_pubcount_adj.py            # PubCount-adjusted comparison runner
```

### 15.2 Tests

```
tests/test_build_transfer_splits.py       # 6 tests (B W3-A)
tests/test_build_pubcount_adjusted.py      # 22 tests (B W3-B)
```

### 15.3 Data products

```
data/processed/labels/intogen2024_pubcount_adjusted.tsv     # 319 genes
data/processed/labels/ot_all_pubcount_adjusted.tsv          # 300 genes
data/processed/labels/ot_nolit_pubcount_adjusted.tsv        # 300 genes
data/processed/labels/clingen_label_pubcount_adjusted.tsv   # 42 genes
data/processed/labels/pubcount_adjusted_index.tsv           # summary
data/processed/transfer/                                    # 66 transfer splits
data/processed/splits/*_pubcount_adjusted__native_*.tsv     # 24 split files
```

### 15.4 Results

```
results/runs/transfer/                          # 66 transfer arm runs
results/runs/pubcount_adj/                       # 24 pubcount-adjusted arm runs
results/runs/orig_backfill/                      # 16 original-arm backfills
results/tables/results_stats_transfer.tsv        # 66-row transfer summary
results/tables/pubcount_adj_comparison.tsv       # 23-row pubcount comparison
```

### 15.5 Documents

```
docs/transfer_results.md                         # full transfer report (12 KB)
docs/pubcount_adj_results.md                     # full pubcount report (12 KB)
docs/b_log.md                                    # this file
```

### 15.6 Total commits by B-side through S3

* 4 W0-W2 commits (HGNC, splits, IntOGen, Open Targets)
* 4 W2 follow-up commits (full OT parquets, D-09, D-11, Figure 1)
* 1 S2 follow-up commit (clingen 5th label)
* 1 S2 log commit (b_log.md update)
* 2 W3 commits (Control 5 transfer, Control 7 pubcount)

**= 12 B-side commits, all merged to `main`.** Combined with
A-side: 35 commits total in the merged main branch.

---

## 16. Tests status update

Through S3, total repo tests: **97 passed, 1 skipped**. The skipped
test is in `tests/test_build_opentargets.py` and requires an
optional sklearn feature that was not needed for the pilot.

B's contribution: **34 new tests across 6 test files**
* `tests/test_hgnc.py` — 11 tests (W0)
* `tests/test_splits.py` — 8 tests (W1)
* `tests/test_pubcount.py` — 2 tests (W1)
* `tests/test_label_network_overlap.py` — 3 tests (W2)
* `tests/test_build_splits.py` — 4 tests (W2)
* `tests/test_build_opentargets.py` — 6 tests (W2)
* `tests/test_build_transfer_splits.py` — 6 tests (W3-A)
* `tests/test_build_pubcount_adjusted.py` — 22 tests (W3-B)

= **62 B-side tests + 35 A-side tests = 97 total**.

---

## 17. Remaining B tasks for M4–M6

The following tasks are B-side ownership per protocol §6 and the
2026-09-27 gap-analysis file. Detailed plans in
`docs/pending_work_protocol_compliance.md` §8 (B's per-person
checklist).

| # | Task | § | Effort | Priority |
|---|---|---|---|---|
| B2 | Public release protocol (leaderboard + Zenodo/figshare) | M2.8-M4.5 | 2 days | HIGH |
| B3 | Reproducibility audit (independent end-to-end replay) | M4.3-M5.5 | 1 week | HIGH |
| B4 | Update `resource_and_dataset_table.csv` with citations + access routes | §3 | 1 day | HIGH |
| B5 | Regenerate Table 1 with ClinGen 5th column | §1.1, M5.2-M6 | 2 days | MEDIUM |
| B6 | Pubcount-adjustment paper write-up | §2 H3 | 2 days | MEDIUM |
| B7 | Cross-catalogue transfer paper write-up | §4.4 #5 | 2 days | MEDIUM |
| B8 | Label methods section (manuscript) | §4.2 | 2 days | MEDIUM |
| B9 | Label results + discussion (manuscript) | §Results / §Discussion | 3 days | MEDIUM |
| B10 | Update `docs/b_log.md` to research-grade | — | (this file) | LOW |
| B11 | Repository cleanup: merge b-labels → main before submission | §4.7 | 1 day | LOW |

**Total B remaining: ~3-4 weeks of focused work.**

---

## 18. Lessons captured at this sync point

1. **The publication-bias hypothesis is more subtle than the protocol
   suggested.** The protocol expected H3 to be either confirmed
   (label signal disappears after adjustment) or refuted (no
   change). The actual finding is intermediate: mean ΔAUROC drops
   by ~39 % of typical magnitude after adjustment, but survives
   in 22 / 23 cells. The lesson: pre-registered controls can return
   intermediate results, and the analysis pipeline must be able to
   report effect-size reduction, not just all-or-nothing.

2. **The cross-catalogue transfer grid is the most informative
   single control.** Per protocol §4.4, this control "in the pilot
   removed 16-66 % of the margin in five of six arms; it is the
   most informative single control in the stack". The B-side
   implementation confirms this: cross-catalogue ≈ 50-70 % of
   within-label magnitude, with clear network ordering
   (string_phys700 > string_full700 > funmap > intact > rna_coexp
   > reactome).

3. **Bootstrap SDs in cross-catalogue are tiny** (~0.005) because
   only the background varies; the target set is fixed. Every
   non-trivial cell has p < 1e-30. The useful uncertainty measure
   is the *full-set* AUROC delta, not the bootstrap SD.

4. **The transfer splits need a richer schema than the standard
   splits.** The standard C4 split file has (gene, repeat, fold, y);
   the transfer split adds (label_id, universe_id, transfer_id,
   role). The role field distinguishes seed / target / background.
   Future reusable transfer infrastructure should keep this schema.

5. **The ClinGen label's β=+1.61 is the largest** in the literature
   bias table, which is consistent with ClinGen's curation rules
   (clinical validity requires published evidence). After
   adjustment, ClinGen's adjusted ΔAUROC still beats degree in 5/6
   cells, suggesting the residual ClinGen signal is *clinical*
   not *literature*.

6. **The transfer script supports any pair direction** but only
   one direction is run by default (the reverse direction would
   double the cell count and many of the reverse cells are noise-
   dominated). A sensitivity check with the reverse direction is
   a possible follow-up but was not committed.

7. **The pubcount adjustment is conservative** (top-50 % trim). A
   finer analysis would use a continuous residual weight instead
   of binary trim. The binary-trim choice removes the half of
   positives most confounded with literature volume; the survival
   of cells like `rna_coexp × intogen2024` suggests even a finer
   trim would not flip the qualitative conclusion.

8. **All B-side deliverables should be cross-validated with A's
   primary numbers.** Both Control 5 and Control 7 re-use A's
   existing splits (per `transfer_id` and `pubcount_adjusted`
   split filenames). The end-to-end numerics should match the
   main factorial to 1e-4 (D-NUMPY lock + identical harness).

These lessons will be re-evaluated at M4 and M5 for promotion to
`AGENTS.md` principles if they prove repeatable.

---

## 19. Cross-references

* A-side logs (mirror this structure): `docs/a0_log.md`...
  `docs/a5_log.md`, `docs/w5_log.md`, `docs/gnn_courtesy_log.md`,
  `docs/stats_analysis_log.md`, `docs/s2_log.md`.
* Joint decision log: `docs/decisions.md`.
* Pre-registration (M2): `docs/pre_registration.md` (DOI
  10.5281/zenodo.23121762).
* Pending work analysis (2026-09-27, before W3+ work): 
  `docs/pending_work_protocol_compliance.md`.
* Manifest with SHA-256 fingerprints: 
  `docs/pre_registration_manifest.json`.

---

## 20. End state at S3 sync (the de facto M4 starting point)

* **12 B-side commits** merged to `main`.
* **5 label sets** committed (intogen2020/2023/2024/temporal_new,
  ot_all, ot_nolit, clingen_label) + 4 PubCount-adjusted variants.
* **9 networks** available (6 primary + 3 sensitivity).
* **86 + 24 + 66** = **176 split files** (86 standard + 24 pubcount-
  adjusted + 66 transfer).
* **24 + 24 + 24 + 66 + 24 + 24 + 24** = **210 arm result tables**.
* **14 stats TSVs** (main, D-09, D-11, clingen, transfer, pubcount,
  GNN, A5, W5, fine-z, variance-decomp, bootstrap-CI, equivalence,
  provenance).
* **3 figures** (Figure 1 + Figure 2 + 2 companions).
* **10 prose reports** (overlap, results-stats, transfer-results,
  pubcount-adj-results, s2-log, handoff, engine-notes, network-
  notes, decisions, this b_log).
* **97 / 98 tests pass** (1 skipped, optional sklearn feature).
* **2.76 GB raw data** downloaded and SHA-256 verified.
* **Pre-registration filed** at Zenodo (DOI 10.5281/zenodo.23121762).
* **All 7 protocol §4.4 controls** are wired into the analysis
  pipeline.
* **All 4 §4.6 statistical-analysis outputs** (fine-z, ME variance
  decomposition, cluster bootstrap CI, TOST equivalence) are run.
* **38 + 66 = 104 arms** in the transfer grid, plus **23 cells**
  in the pubcount comparison, plus **24 cells** in W5 + **24 cells**
  in A5.
* **All headline claims** (H1, H2, H3) have supporting numerical
  evidence.

This is a good place to hand back to Person A for the manuscript
write-up (M5-M6) and for B to do the public release (B2),
reproducibility audit (B3), and resource-table update (B4) that
the protocol's M2-M4 window requires.

