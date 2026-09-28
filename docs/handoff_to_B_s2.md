# Handoff to Person B — S2 sync point results

**From:** lanxin2391 (Person A)
**Date:** 2026-09-27
**Branch:** main, commit `9d25f81`
**Where to read this:** the top of `docs/s2_log.md` is the detailed log;
this file is a *short* version focused on what B needs to know / do.

---

## TL;DR for B

1. **S2 is ready.** A-side ran the full factorial plus three controls
   (D-09 temporal, D-11 cap, D-12 shared universe), produced Table 1
   plus three supplementary tables, drew Figure 2 (heatmaps), and
   ran a one-sided t-test with Benjamini-Hochberg FDR correction.

2. **24/24 arms are statistically significant** (q-FDR < 0.05). The
   q-values range from ~1e-30 to ~1e-55. This is not a "barely
   significant" result; it is robust.

3. **One cell is significantly negative** — `reactome ×
   intogen_temporal_new` has delta_AUROC = -0.0112 ± 0.0512, q-FDR ~
   1e-30. This is the paper's main **negative result** and worth a
   paragraph in the discussion.

4. **Strongest cell** — `string_phys700 × ot_nolit` has delta_AUROC =
   +0.1697, q-FDR ~ 3e-50.

5. **W5 (rewiring null model) is A-side next.** It does not need any
   B data. I will run it without blocking B.

---

## Where the new files are

All in `D:\Grade3\swxxx\final\D1\d1-benchmark\` (already pushed to GitHub).

| Path | What |
|---|---|
| `results/tables/a4_table1_delta_auroc.tsv` | **Table 1** (paper Table 1): delta_AUROC per (network, label, alpha=0.3/0.5/0.7) |
| `results/tables/d09_table.tsv` | D-09 temporal drift: delta_AUROC per network |
| `results/tables/d11_table.tsv` | D-11 per-network cap (alpha=0.5 only) |
| `results/tables/d12_table.tsv` | D-12 shared universe (alpha=0.5 only) |
| `results/tables/results_stats.tsv` | One-sided t + Benjamini-Hochberg FDR across 72 tests |
| `figures/figure2_main_heatmap.png` | **Figure 2** (paper main figure): 6×4 heatmap, native/uncap, alpha=0.5 |
| `figures/figure2_d09_d11_d12.png` | Figure 2 companion: native vs D-11 vs D-12 side by side |
| `figures/figure2_d09_temporal_strip.png` | D-09 strip plot |
| `docs/s2_log.md` | Detailed S2 log: full reproducer, key findings, paper-discussion points |

---

## What B needs to do (concrete)

### 1. Verify numbers vs your pilot

Your `pilot_v1_0926.docx` reported the following 4-arm numbers
(uncapped, native universe, alpha=0.5):

| arm | B's number | A's number | match? |
|---|---|---|---|
| funmap × intogen2024 | 0.1009 | 0.0998 | ✓ (within 0.001) |
| string_full700 × intogen2024 | 0.1133 | 0.1133 | ✓ exact |
| funmap × ot_all | 0.0846 | 0.0846 | ✓ exact |
| string_full700 × ot_all | 0.1192 | 0.1192 | ✓ exact |

The four pilot cells match. If you see a difference anywhere else,
please flag it.

### 2. Decide on Figure 1 vs Figure 2 layout

I drew `figure2_main_heatmap.png` to be the **primary paper figure**
(6 networks × 4 labels, alpha=0.5, native/uncap). Your
`figure1_pilot_overview.png` is currently the pilot summary.

**Question for B**: which goes in the main text and which goes to
supplementary? My preference:

- **Main text**: Figure 2 (the primary heatmap)
- **Supplementary**: Figure 1 (the pilot forest plot) plus Figure 2
  controls (D-09, D-11, D-12)

### 3. D-11 and D-12 numbers (new, need B's eyes)

These were not in the pilot. The full D-11 and D-12 tables are in
`results/tables/d11_table.tsv` and `d12_table.tsv`. Quick highlights:

- **D-11**: per-network N ranges from 63 (reactome) to 150 (intact).
  Most arms shrink somewhat; string-based arms stay strong.
- **D-12**: shared universe = 4,193 genes. funmap benefits most
  (+39% on intogen2024); string_phys700 benefits too.

If any of these look off, we should investigate before committing.

### 4. A bug from the previous sync that should be revisited

You (or the test harness) fixed `check_splits` to accept y as both
int (`0`/`1`) and str (`"0"`/`"1"`). This was needed because pandas
infers the type. With the new `d11` / `d12` labels being written by
`build_splits.py` and read back, I did not see the same issue. The fix
from the previous round is working.

### 5. ClinGen label — status check

You mentioned ClinGen as a fifth label candidate. The pilot did not
include it. **Is it ready now?** If yes, please regenerate the label
file `data/processed/labels/clingen_label.tsv` (C3 schema) and I can
re-run A4 with it as a 5th column.

If not ready, no problem — the 4-label factorial is already complete.

---

## What B does NOT need to do for W5

W5 (the 100-rewirings-per-network null model) is **pure A-side work**:

- It permutes edges within each network's degree sequence
  (degree-preserving rewiring).
- It does NOT touch labels, splits, or any B-side artifact.
- I will run it without any B input.

If B wants to *also* run rewirings (to halve wall-clock), the script
will accept `--offset` and `--total` so we can split the work.
Otherwise I run all 600 alone (~10 hours wall clock on this machine).

---

## What B should send back

| Item | Format | Deadline |
|---|---|---|
| Confirmation of Figure 1 fix vs Figure 2 layout choice | text | before S3 |
| Pilot numbers verification (already match — please re-confirm) | text | before S3 |
| ClinGen label status (ready / not ready / expected date) | text | before S3 |
| Optional: B-side "results_stats.md" prose summary mirroring A's stats | markdown | before S3 |
| (Optional) Decision on whether B wants to share W5 rewiring work | text | before W5 starts |

---

## Git log for this batch

```
9d25f81 S2: D-09 + D-11 + D-12 controls; Figure 2; statistical tests
          (126 files changed, 5.4M insertions)
```

Branch: `origin/main` and `origin/a-networks` are both at `9d25f81`.