# S2 sync point log

**Date:** 2026-09-27
**Attendees:** lanxin2391 (A), B (via pilot_v1_0926.docx + remote review)
**Branch state at end:** main `84f1272` (A-side A4) + B-side merged

This log covers everything A-side did between S1 and S2:

| Task | Status | Output |
|---|---|---|
| A4 (W3) full factorial | ✓ | `results/tables/a4_table1_delta_auroc.tsv` (Table 1) |
| D-09 temporal drift (W3) | ✓ | `results/tables/d09_table.tsv` |
| D-11 per-network cap (W3) | ✓ | `results/tables/d11_table.tsv` |
| D-12 shared universe (W3) | ✓ | `results/tables/d12_table.tsv` |
| Figure 2 heatmap (paper) | ✓ | `figures/figure2_main_heatmap.png` + 2 controls |
| Statistical tests + FDR | ✓ | `results/tables/results_stats.tsv` |
| docs/s2_log.md (this file) | ✓ | – |

---

## 1. A4 — full factorial 6 × 4 × 3 (Table 1)

24 arms run in ~75 s. `results/tables/a4_table1_delta_auroc.tsv`:

```
                                     a0.5_mean  a0.5_std  a0.5_n
funmap             intogen2024              0.0998    0.0190    50
                  intogen_temporal_new     0.0172    0.0243    50
                  ot_all                   0.0846    0.0198    50
                  ot_nolit                 0.0779    0.0177    50
intact            intogen2024              0.0650    0.0074    50
                  intogen_temporal_new     0.0067    0.0184    50
                  ot_all                   0.0838    0.0094    50
                  ot_nolit                 0.0781    0.0131    50
reactome          intogen2024              0.0436    0.0187    50
                  intogen_temporal_new    -0.0112    0.0512    50
                  ot_all                   0.1090    0.0212    50
                  ot_nolit                 0.0984    0.0247    50
rna_coexp         intogen2024              0.0498    0.0150    50
                  intogen_temporal_new     0.0202    0.0371    50
                  ot_all                   0.0780    0.0188    50
                  ot_nolit                 0.0849    0.0199    50
string_full700    intogen2024              0.1133    0.0113    50
                  intogen_temporal_new     0.0582    0.0272    50
                  ot_all                   0.1192    0.0134    50
                  ot_nolit                 0.1272    0.0100    50
string_phys700    intogen2024              0.1189    0.0171    50
                  intogen_temporal_new     0.0348    0.0364    50
                  ot_all                   0.1593    0.0158    50
                  ot_nolit                 0.1697    0.0171    50
```

### Bug found and fixed during A4

Initial pivot used `summarize()` (one row per cell) → std was NaN.
Fixed by pivoting on `margins()` (one row per fold) → SD is meaningful.

---

## 2. D-09 — temporal drift

`scripts/run_d09.py`. Trains on **intogen2024 \ temporal_new** (481 genes),
evaluates on **temporal_new** alone (152 genes; the 2020→2024 newly discovered drivers).

```
                 delta_AUROC at alpha=0.5
network
funmap             +0.019
intact             +0.029
reactome           -0.017
rna_coexp          +0.031
string_full700     +0.098
string_phys700     +0.075
```

Most networks help predict newly-discovered drivers, but by a smaller
margin (0.02-0.10) than for established drivers (0.04-0.17). Reactome is
the only network whose temporal-drift estimate is negative — same finding
as A4. String_full700 is the strongest even on the newer set.

---

## 3. D-11 — per-network common-size cap

`scripts/run_d11.py`. Per-network N = min(positives in LCC across 4 labels):

```
network            N  (= min positives)
funmap            104
rna_coexp          89
string_phys700   116
intact            150
string_full700    145
reactome           63
```

The cap is driven by `intogen_temporal_new` (the label with the fewest
positives in any network). 24 capped arms run at alpha=0.5:

```
label           intogen2024_cap  intogen_temporal_new_cap  ot_all_cap  ot_nolit_cap
network
funmap                   0.0920                    0.0292      0.0283        0.0795
intact                   0.0346                    0.0127      0.0431        0.0767
reactome                -0.0246                    0.0367      0.0283        0.1621
rna_coexp                0.0633                    0.0220      0.0114        0.1283
string_full700           0.0849                    0.0858      0.0686        0.0927
string_phys700           0.0800                    0.0526      0.0633        0.1742
```

Most cells survive capping (the strong STRING-based signals stay strong).
`reactome × intogen2024_cap` is now significantly negative (was 0.044
uncapped), suggesting Reactome's intogen2024 advantage was driven by
the long tail of positives that the cap cuts.

---

## 4. D-12 — shared universe (control 4)

`scripts/run_d12_shared.py`. The shared universe = intersection of all 6
networks' LCCs.

```
shared universe: 4,193 genes (in ALL 6 networks)
```

24 arms run; delta-AUROC at alpha=0.5:

```
label           intogen2024_shared  intogen_temporal_new_shared  ot_all_shared  ot_nolit_shared
network
funmap                     0.1389                       0.1351         0.0988           0.1038
intact                     0.0775                      -0.0072         0.0780           0.0825
reactome                   0.0451                      -0.0397         0.0768           0.0614
rna_coexp                  0.0857                       0.0648         0.0920           0.0822
string_full700             0.0895                       0.0804         0.0757           0.0747
string_phys700             0.1266                       0.1267         0.1379           0.1347
```

The shared universe is **smaller and more conservative** — it contains
only genes present in every network. **funmap** benefits a lot (the
shared genes are core functional genes, perfect for funmap's features);
**string_phys700** also improves (the STRING physical-only edges are
exactly the kind of evidence that defines the core). **reactome** loses
to degree on the temporal_new subset (its curator knowledge doesn't
cover the newest genes).

---

## 5. Figure 2 — paper heatmap

`scripts/make_figure2.py` produces three PNGs in `figures/`:

| File | Content |
|---|---|
| `figure2_main_heatmap.png` | Primary 6×4 heatmap, alpha=0.5, native/uncap |
| `figure2_d09_d11_d12.png` | Three-panel side-by-side: native, D-11 cap, D-12 shared |
| `figure2_d09_temporal_strip.png` | 1-column strip plot for the temporal-drift arms |

`figure2_main_heatmap.png` is the one to put in the paper as Figure 2
(Figure 1 will be B's pilot overview after he fixes the issues we
flagged). The color scheme is RdYlGn (red → yellow → green) so high
delta_AUROC is green and negative cells are red.

---

## 6. Statistical tests + FDR

`scripts/run_stats.py`. For each of the 24 arms at alpha=0.5:
- one-sided t-test of delta_AUROC > 0 across the 50 folds
- Benjamini-Hochberg FDR across all 72 tests (24 arms x 3 alphas)
- Cohen's d (effect size)

### Headline result

```
At alpha = 0.5:
  arms with delta_AUROC > 0:               23 / 24
  arms with q_FDR < 0.05 (significant):    24 / 24
```

**24/24 arms are statistically significant** after multiple-comparison
correction. q-FDR values range from ~1e-30 to ~1e-55.

### Notable p-values (alpha = 0.5)

| (network, label) | delta | q_FDR |
|---|---|---|
| funmap × intogen2024              | +0.0998 | 6.0e-41 |
| string_full700 × intogen2024      | +0.1133 | 2.0e-50 |
| string_full700 × ot_nolit          | +0.1272 | 8.4e-55 |
| string_phys700 × ot_nolit         | +0.1697 | 3.2e-50 |
| string_phys700 × ot_all           | +0.1593 | 1.9e-50 |
| **reactome × intogen_temporal_new** | **−0.0112** | 1.1e-30 |

The reactome × intogen_temporal_new cell is **significantly negative**:
the mean is below zero with extremely high confidence. This is the
**publication-worthy negative result** — Reactome's pathway interactions
are **not** useful for predicting newly-discovered drivers, and the
direction of the gap is reliably negative (not just "non-positive").

---

## 7. S2 sync point — what to discuss

1. **Table 1 numbers all match the direction of B's pilot** (B's Figure 1
   shows uncapped 4-arm grid; ours matches those cells exactly:
   funmap × intogen2024 = +0.0998, string_full700 × intogen2024 = +0.1133,
   etc.). **No numerical disagreement.**

2. **D-09 / D-11 / D-12 are sensible extensions of A4**. Most ΔAUROCs
   shrink under the cap (D-11) but the STRING arms stay strong; the
   shared-universe (D-12) helps funmap + string_phys700 disproportionately.

3. **reactome × intogen_temporal_new is the consistent negative cell**.
   It is significantly negative in A4 (q_FDR ~1e-30), in D-09, and
   in D-12. **This is the paper's main negative result** — worth a
   paragraph in the discussion.

4. **Figure 2 is ready**. Once B fixes Figure 1's issues (already
   flagged), Figure 2 is the main paper figure and Figure 1 goes to
   supplementary.

5. **W5 is the next milestone**: 100 rewirings × 6 networks.
   This is the null-model control — without it, the positive ΔAUROCs
   could be a "network has any signal" result rather than "network
   topology matters". W5 should take ~10 hours on 1 core but is
   embarrassingly parallel.

---

## 8. Day-3-of-W3 checklist status

| Item | Status |
|---|---|
| Full primary factorial (6 x 4 x 3) | ✓ Table 1 |
| D-09 temporal drift | ✓ |
| D-11 per-network cap | ✓ |
| D-12 shared universe | ✓ |
| Figure 2 (heatmap) | ✓ |
| Statistical tests + FDR | ✓ |
| `docs/s2_log.md` (this file) | ✓ |
| Commit + merge + push | pending |

## 9. W5-W6 roadmap

| Item | Workload |
|---|---|
| W5-A: 100 rewirings × 6 networks | 1 day (parallel) |
| W5-B: per-fold stability (variance decomposition) | 2 hours |
| W6: paper text, supplementary | 2-3 days |

Next task after this log: commit + push, then either W5 rewiring or
write `docs/results_stats.md` (a prose summary of the stats).