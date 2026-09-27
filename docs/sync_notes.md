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