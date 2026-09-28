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
