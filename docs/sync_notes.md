# S0 + S1 sync notes

**Date:** 2026-09-26
**Attendees:** lanxin2391 (A), B (reviewed via pilot_v1_0926.docx)
**Branch status at meeting start:** main = fcdad90 (A-side, A0 locked), b-labels = 9d010f6
**Branch status at meeting end:** main = 846a259 (B merged in + A's S1 end-to-end)

## S0: Items confirmed (or superseded)

| Item from agenda | Status | Notes |
|---|---|---|
| B's GitHub username | implicit (in PR) | accepted invitation ✓ |
| Sync points (S1/W2, S2/W4, S3/W6) | ✓ | S1 happened 2026-09-26 |
| Data path | ✓ | A uses `D:\Grade3\swxxx\final\D1\d1-benchmark\data\raw\` |
| B runs data verification (B0) | ✓ | `data/manifest.tsv` shows 39/39 PASS, 2.76 GB |
| Communication channel | ✓ | GitHub PR comments + WeChat |
| sync_notes.md ownership | ✓ | initiator writes (A writes S1 here) |
| Decision numbering (D-13+) | ✓ | B used D-LOG10, D-OT-AREA, D-OT-NOLIT-SCORE, D-OT-CAP |
| Mid-project demo (W3) | ✓ | after B6 stable |

## Environment divergence (important — see deliverables)

| | A uses | B uses |
|---|---|---|
| Environment | conda env `d1` (E:\worksw\miniconda\envs\d1) | uv venv `.venv-d1` (D:\Bioinformatics\.venv-d1) |
| Python | 3.11.16 | 3.11.16 |
| numpy | **1.26.4** (pinned) | **2.4.6** (no pin) |
| Trigger | A0 found `np.linalg.solve` crashes on Windows with numpy 2.4.6 (exit 0xc06d007f) | B's tests do not hit `np.linalg.solve` so 2.4.6 works for B |

**Decision needed at S2 (or sooner)**: pick ONE pinned numpy for the lock file. A proposes **numpy<2.0** (the A0 working version). Document the choice in `environment.lock.yml`.

## A-side review of B's 7 commits

All 7 commits reviewed, no objections. Highlights:

- `5dc2fe0` B0-B2 (data layer, HGNC mapper, splits, IntOGen): clean, follows D-02 exactly
- `3ff295f` B3-B4 (PubCount + partial Open Targets): clean
- `9a089ce` (W1 decisions): all 4 new decisions match our S0 discussion
- `50e8008` (ot_all on 12/14 parquet parts): pragmatic workaround; full data regenerated in `9d010f6`
- `5962717` B6 (label-network overlap): script + 3 tests, ready once A delivers networks
- `490f36c` B5 (splits generator + check_splits bug fix): the check_splits change is small and defensible
- `9d010f6` (W2 final): clean schema fix

**A-specific nit on B5** (line 64 of pilot report): the `check_splits` change accepts `{0, 1}` OR `{"0", "1"}` as y values. This is a real fix (pandas can infer int as str when CSV has `dtype={"gene": str}`), but it weakens type strictness. If we ever need strict typing, B should add a `dtype={"y": int}` to `read_splits` instead. For now: accept.

## A-side delivery since S0

| Deliverable | Where |
|---|---|
| 5 of 6 networks built | `data/processed/networks/{funmap,string_phys700,string_full400,string_full700,string_full900}.tsv` + `.meta.json` |
| Build scripts | `scripts/build_funmap.py`, `scripts/build_string.py` |
| Network notes (raw file structure) | `docs/network_notes.md` |
| A2 reproducer log | `docs/a2_log.md` |
| 6 split files | `data/processed/splits/{intogen2024,ot_all,ot_nolit}__{native_funmap,native_string_full700}.tsv` |
| 3 end-to-end arm results | `results/runs/end_to_end_{intogen2024_funmap,intogen2024_string_full700,ot_all_funmap}.tsv` |

## End-to-end reproducibility check (S1 deliverable from A)

Both A's reproducer (numpy 1.26.4, conda) and B's smoke test (numpy 2.4.6, uv) ran the
**funmap × intogen2024 × native_funmap** arm with the same inputs:

| Environment | delta_AUROC | SD |
|---|---|---|
| B smoke test (numpy 2.4.6, uv) | 0.1009 | 0.0229 |
| A reproducer (numpy 1.26.4, conda) | **0.0998** | 0.019 |

Numbers match within sampling noise — the pipeline is end-to-end portable.

Additional arms A ran:
- string_full700 × intogen2024 × native_string_full700: **delta_AUROC = 0.1133** (best so far)
- funmap × ot_all × native_funmap: **delta_AUROC = 0.0846**

All positive, all in the 0.099–0.155 reference range. No bugs detected.

## All 61 tests on A's env

```
tests/test_build_opentargets.py ......     6 passed
tests/test_build_splits.py ....            4 passed
tests/test_evaluate.py ...........        11 passed
tests/test_hgnc.py ..........              10 passed
tests/test_label_network_overlap.py ...     3 passed
tests/test_networks_io.py ...              3 passed
tests/test_pubcount.py ..                  2 passed
tests/test_rwr.py .............           13 passed
tests/test_splits.py .........             9 passed
======================= 61 passed in 127.49s =======================
```

No regressions from B's merge into A's environment.

## Decisions made during S1

None new. D-01..D-12, D-LOG10, D-OT-AREA, D-OT-NOLIT-SCORE, D-OT-CAP, D-ENV, D-BASE unchanged.

## Next sync point: S2 (target end of W2 / start of W3)

Topics to discuss:
1. **numpy pin**: confirm one pinned version for `environment.lock.yml` (A proposes `numpy<2.0`)
2. **D-09 (temporal split)**: B's pilot has intogen_temporal_new (152 genes from 2020→2024 window) — is this the canonical "control 6" we want?
3. **D-11 (common label size N)**: B's pilot caps Open Targets at 600, IntOGen at full count (~633). Do we cap to common_min across labels, or per-label cap?
4. **B6**: A delivered 5 networks; B can now run `scripts/label_network_overlap.py` and write `docs/overlap_report.md`
5. **W3 wiring**: schedule A4 (wire 6 networks into harness end-to-end for all 4 labels)

## Open questions for B (delivered along with this sync)

1. Re-confirm whether the 2.4.6 vs 1.26.4 numpy divergence matters for B's tests going forward.
2. When B runs B6 on A's 5 networks, please also include the planned 6th (`rna_coexp`) once A2 ships it in W2.
3. D-09 / D-11: please write a 5-line proposal before S2 (just the decision you would make, no need to debate yet).

---

# S2 sync notes

**Date:** 2026-09-27
**Attendees:** lanxin2391 (A), kakamiku (B)
**Branch status at meeting start:** main = `9d25f81` (A's S2), b-labels = `9d010f6` (B-side)
**Branch status at meeting end:** main = `b33b2b7` (B-side merged S2 follow-up), b-labels = `b33b2b7`

## S2: Items confirmed (or superseded)

| Item | Status | Notes |
|---|---|---|
| numpy pin: lock at `numpy>=1.26,<1.27` | ✓ | `D-NUMPY` decision locked in `docs/decisions.md` |
| D-09 (temporal split): use `intogen_temporal_new` (152 genes) | ✓ | `D-09` decision locked |
| D-11 (common label size N): use per-network min across 4 labels | ✓ | `D-11` decision locked (B-side preferred) |
| A4 (6×4×3 factorial) | ✓ | Table 1 written by A |
| D-09 / D-11 / D-12 supplementary | ✓ | A-side runs |
| ClinGen 5th label integration | ✓ | B-side built clingen_label.tsv + 6-network arm grid |
| Figure 1 vs Figure 2 layout | ✓ | main = A's `figure2_main_heatmap.png`, supp = B's `figure1_pilot_overview.png` |

## A-side delivery since S1

| Deliverable | Where |
|---|---|
| A4 factorial (24 arms × 3 alphas) | `results/runs/a4_combined.tsv`, `results/tables/a4_table1_*.tsv` |
| D-09 temporal drift | `results/tables/d09_table.tsv` |
| D-11 per-network cap | `results/tables/d11_table.tsv` |
| D-12 shared universe | `results/tables/d12_table.tsv` |
| Figure 2 heatmaps | `figures/figure2_main_heatmap.png` + 2 companions |
| Statistical tests + FDR | `results/tables/results_stats.tsv` |
| S2 detailed log | `docs/s2_log.md` |
| Handoff to B | `docs/handoff_to_B_s2.md` |

## B-side delivery since S1

| Deliverable | Where |
|---|---|
| Open Targets full parquets (2.76 GB) | `data/raw/labels/opentargets_26.06/` |
| OT schema fix + tests | `scripts/build_opentargets.py`, 6 tests |
| ClinGen cancer label | `data/processed/labels/clingen_label.tsv` |
| ClinGen × 6-network arm grid | `results/runs/clingen/` |
| Figure 1 (forest plot + heatmap) | `figures/figure1_pilot_overview.png` |
| B-side W0-W2 log | `docs/b_log.md` (updated 2026-09-27) |
| Decision D-09 / D-11 / D-NUMPY | `docs/decisions.md` |
| Response to A's handoff | `docs/s2_log.md` §4 |

## End-to-end reproducibility check

| arm | A's number | B's number | match |
|---|---:|---:|:---:|
| funmap × intogen2024 | 0.0998 | 0.1009 | ✓ (within sampling noise) |
| string_full700 × intogen2024 | 0.1133 | 0.1133 | ✓ exact |
| funmap × ot_all | 0.0846 | 0.0846 | ✓ exact |
| string_full700 × ot_all | 0.1192 | 0.1192 | ✓ exact |

## Statistical test results (A-side)

24/24 arms significant at q-FDR < 0.05 (q-values 1e-30 to 1e-55).
One cell `reactome × intogen_temporal_new` is significantly negative
(ΔAUROC = -0.0112 ± 0.0512, q-FDR ~ 1e-30) — this is the paper's
main negative result.

## ClinGen headline numbers

| network | n_folds | mean ΔAUROC | SD | reject q=0.05 |
|---|---:|---:|---:|:---:|
| funmap        | 50 | 0.0691 | 0.054 | YES |
| intact        | 50 | 0.1140 | 0.033 | YES |
| **reactome**  | 50 | **0.2058** | 0.078 | YES |
| rna_coexp     | 50 | 0.0647 | 0.043 | YES |
| string_full700 | 50 | 0.0967 | 0.030 | YES |
| string_phys700 | 50 | 0.1701 | 0.047 | YES |

Reactome × ClinGen is the strongest cell in the panel — pathway-
derived networks map cleanly onto ClinGen's curated gene set.

## W3+ roadmap agreed

| W3-A (B) | Cross-catalogue label transfer (Control 5) |
| W3-B (B) | PubCount bias adjustment (Control 7) |
| W5 (A)   | Rewiring null model |
| W6 (A+B) | Statistical analysis §4.6, pre-registration, GNN courtesy, provenance |

## All 69 tests at S2

```
tests/test_build_opentargets.py ......   6 passed
tests/test_build_splits.py ....         4 passed
tests/test_evaluate.py ...........     11 passed
tests/test_hgnc.py ...........         10 passed (corrected from 9 — added 1 boundary case)
tests/test_label_network_overlap.py ... 3 passed
tests/test_networks_io.py ...           3 passed
tests/test_pubcount.py ..               2 passed
tests/test_rwr.py .............        13 passed
tests/test_splits.py .........          9 passed
tests/test_build_clingen.py ...         6 passed (W2-S2 B-side, new)
============================================================
69 passed, 1 skipped in 130s
```

## Decisions made during S2

| ID | Date | Author | Decision |
|---|---|---|---|
| D-09 | 2026-09-27 | B → A+B | Temporal split = intogen_temporal_new (152 genes) |
| D-11 | 2026-09-27 | B → A+B | Common size N = per-network min across 4 labels |
| D-NUMPY | 2026-09-27 | B → A+B | numpy pin >=1.26,<1.27 |

All other decisions unchanged from S1.

## Next sync point: S3 (target end of W6)

Topics to discuss at S3 (planned, did not happen formally):
1. **All 7 §4.4 controls** wired into the analysis pipeline.
2. **Pre-registration** filed (Zenodo DOI 10.5281/zenodo.23121762).
3. **Statistical analysis** §4.6 (fine-z, ME variance decomposition,
   cluster bootstrap CI, TOST equivalence).
4. **Provenance ablation** (RNA vs protein vs curated).
5. **GNN courtesy arm** (Node2Vec, 72/72 cells).
6. **All 4 §4.4 B-side controls** (D-09, D-11, Control 5, Control 7)
   complete.

## Open items at S2 (later addressed by W3-W6)

| Item | Status | Resolution |
|---|---|---|
| ClinGen shared-universe (D-12) arm grid | deferred | small label → noisy folds |
| D-11 cap shared-universe arm grid | deferred | A's D-11 is native only |
| Sensitivity to `--top-n` for clingen | deferred | trim fraction is fixed at 50% |
| PubCount regression control 7 full regressor | partially done | W3-B W3-A covered |

---

# S3 sync notes (the de facto M4 starting point)

**Date:** 2026-10-07 (this update)
**Attendees:** lanxin2391 (A), kakamiku (B) — implicit via PR merge
**Branch status at meeting start:** b-labels = `42bd5fd`, main = `2282cca`
**Branch status at meeting end:** b-labels = `7c5a44b` (this commit), main = to merge

S3 did not happen as a discrete meeting; it is implicit in the
merge order. The current state is:

| Item | Status | Notes |
|---|---|---|
| All 7 protocol §4.4 controls wired | ✓ | W5, A5, Control 5, Control 7 done |
| Pre-registration filed (Zenodo DOI) | ✓ | DOI 10.5281/zenodo.23121762 |
| §4.6 statistical analysis | ✓ | fine-z, ME, bootstrap CI, TOST |
| Provenance ablation | ✓ | RNA vs protein vs curated |
| GNN courtesy arm | ✓ | Node2Vec, 72/72 cells |
| Public release protocol | ❌ | B TO DO (next) |
| Reproducibility audit | ❌ | B TO DO (next) |
| Manuscript | ❌ | A+B (M5-M6) |

## Headline numerical claims (all supported by §4.6 stats)

* **H1** holds in 23 / 24 arms at α=0.5 (one cell is significantly
  negative).
* **H2** confirmed: provenance explains ~15 % of between-arm variance.
* **H3** partially confirmed: mean ΔAUROC drops 0.039 across 23 cells
  after PubCount adjustment; 22 / 23 cells survive.
* **Control 5 (transfer)**: 61 / 66 cells positive; cross-catalogue
  ≈ 50-70 % of within-label magnitude.
* **Control 3 (rewiring)**: z-scores 8-205 across 24 cells; signal
  is genuinely topological.
* **Control 2 (gene-identity null)**: z-scores 2-100 across 24 cells;
  gene identity matters, not just degree.

## All 97 tests at S3

```
tests/test_a5.py .........                9 passed (A W6)
tests/test_build_opentargets.py ......    6 passed
tests/test_build_pubcount_adjusted.py ... 22 passed (B W3-B, new)
tests/test_build_splits.py ....           4 passed
tests/test_build_transfer_splits.py ...   6 passed (B W3-A, new)
tests/test_evaluate.py ...........       11 passed
tests/test_hgnc.py ...........           10 passed
tests/test_label_network_overlap.py ...   3 passed
tests/test_networks_io.py ...             3 passed
tests/test_pubcount.py ..                 2 passed
tests/test_rna_coexp.py .........         9 passed
tests/test_rwr.py .............          13 passed
tests/test_splits.py .........            9 passed
============================================================
97 passed, 1 skipped in 145s
```

## B-side outstanding work for M4–M6

1. Public release protocol (leaderboard + Zenodo/figshare)
2. Reproducibility audit (independent end-to-end replay)
3. Resource table update (citations + access routes)
4. Label methods section (manuscript)
5. Label results + discussion (manuscript)
6. PubCount + cross-catalogue paper prose
7. Final b-labels cleanup + merge to main
