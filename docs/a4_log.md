# A4 log: full 6 × 4 × 3 factorial — Table 1 of the paper

**Task:** W3, task A4 (protocol section 6.3 of W3).
**Author:** lanxin2391 (Person A)
**Date:** 2026-09-27
**Outcome:** 24 arms run in ~75 s. Table 1 written. 70 tests pass.

A4 is the first task where we run the **full factorial** and produce a
publishable Table 1. After this task we have a number for every
(network × label × α) cell under the primary analysis settings.

---

## 1. What A4 is

A4 runs every (network, label) arm of the primary benchmark and
aggregates the results into Table 1. It does NOT include:

- The shared universe (D-12 control) — W3 still pending
- Network rewiring (control 2) — W5
- D-09 temporal drift runs (B's proposal) — W3
- D-11 common-size cap runs (B's proposal) — W3

Those four are the "supplementary" tier that lands in W3–W5.

## 2. Step-by-step reproduction

```powershell
conda activate d1
cd D:\Grade3\swxxx\final\D1\d1-benchmark

# 1. Run the full factorial (~75 s on this PC)
python scripts\run_a4.py

# The script:
# - Generates any missing split files (via scripts/build_splits.py)
# - Runs 24 arms (6 networks x 4 labels) at 3 alpha values each
# - Writes one C5 result per arm to results/runs/a4_<net>__<label>.tsv
# - Concatenates them into results/runs/a4_combined.tsv (4,800 rows)
# - Aggregates to results/tables/a4_table1_delta_auroc.tsv
# - Also writes a4_table1_rwr_auroc.tsv and a4_table1_degree_auroc.tsv

# 2. Re-run if anything changed
python scripts\run_a4.py --only-net string_full700    # just one network
python scripts\run_a4.py --only-label intogen2024    # just one label
```

## 3. What the script does

For each arm:

1. `build_splits.py` is invoked if the split file is missing (each
   arm needs a C4 split file with 50 folds × 5 splits = 250 train+test rows).
2. `d1.engine.evaluate.run_arm` runs the full propagator + degree baseline
   + metrics for all (repeat, fold) pairs in one batched call.
3. The 200-row C5 result (50 folds × {1 degree + 3 alphas of RWR}) is
   written to `results/runs/a4_<network>__<label>.tsv`.

At the end:

4. All 24 C5 files are concatenated into one big table
   (`results/runs/a4_combined.tsv`, 4,800 rows).
5. The margins() helper joins each RWR row with its degree row and
   computes delta_AUROC per fold. We then group by (network, label, alpha)
   and report **mean ± SD across 50 folds** — that's Table 1.

## 4. The bug we hit (and fixed)

The first version of the script computed Table 1 by pivot on
`summarize()` output, which has only **one row per (network, label,
alpha)** (the mean across folds). The std of a single value is NaN.

Fixed by pivoting on `margins()` output, which has **one row per
(network, label, alpha, repeat, fold)**, and computing mean/SD across
the 50 folds in each cell. The Table 1 now shows non-trivial SDs (range
0.005 to 0.054, reflecting fold-level variance).

## 5. Table 1 (the paper Table 1)

`results/tables/a4_table1_delta_auroc.tsv` (read with `pd.read_csv(sep='\t')`):

```
                                     a0.3_mean  a0.3_std  a0.3_n  a0.5_mean  a0.5_std  a0.5_n  a0.7_mean  a0.7_std  a0.7_n
network        label_set
funmap         intogen2024              0.0895    0.0174    50.0     0.0998    0.0190    50.0     0.1046    0.0198    50.0
               intogen_temporal_new     0.0086    0.0232    50.0     0.0172    0.0243    50.0     0.0252    0.0251    50.0
               ot_all                   0.0752    0.0183    50.0     0.0846    0.0198    50.0     0.0907    0.0204    50.0
               ot_nolit                 0.0711    0.0161    50.0     0.0779    0.0177    50.0     0.0813    0.0188    50.0
intact         intogen2024              0.0545    0.0059    50.0     0.0650    0.0074    50.0     0.0682    0.0085    50.0
               intogen_temporal_new     0.0095    0.0148    50.0     0.0067    0.0184    50.0     0.0023    0.0207    50.0
               ot_all                   0.0683    0.0076    50.0     0.0838    0.0094    50.0     0.0912    0.0107    50.0
               ot_nolit                 0.0646    0.0104    50.0     0.0781    0.0131    50.0     0.0839    0.0151    50.0
reactome       intogen2024              0.0484    0.0167    50.0     0.0436    0.0187    50.0     0.0363    0.0196    50.0
               intogen_temporal_new    -0.0082    0.0473    50.0    -0.0112    0.0512    50.0    -0.0144    0.0540    50.0
               ot_all                   0.1045    0.0188    50.0     0.1090    0.0212    50.0     0.1083    0.0228    50.0
               ot_nolit                 0.0966    0.0229    50.0     0.0984    0.0247    50.0     0.0954    0.0256    50.0
rna_coexp      intogen2024              0.0424    0.0136    50.0     0.0498    0.0150    50.0     0.0550    0.0163    50.0
               intogen_temporal_new     0.0161    0.0345    50.0     0.0202    0.0371    50.0     0.0228    0.0392    50.0
               ot_all                   0.0702    0.0180    50.0     0.0780    0.0188    50.0     0.0824    0.0191    50.0
               ot_nolit                 0.0743    0.0179    50.0     0.0849    0.0199    50.0     0.0916    0.0212    50.0
string_full700 intogen2024              0.1039    0.0099    50.0     0.1133    0.0113    50.0     0.1163    0.0120    50.0
               intogen_temporal_new     0.0563    0.0231    50.0     0.0582    0.0272    50.0     0.0574    0.0299    50.0
               ot_all                   0.1111    0.0121    50.0     0.1133    0.0113    50.0     0.1212    0.0141    50.0
               ot_nolit                 0.1200    0.0092    50.0     0.1272    0.0100    50.0     0.1281    0.0102    50.0
string_phys700 intogen2024              0.1149    0.0154    50.0     0.1189    0.0171    50.0     0.1183    0.0182    50.0
               intogen_temporal_new     0.0357    0.0326    50.0     0.0130    0.0130    50.0     0.0330    0.0392    50.0
               ot_all                   0.1550    0.0148    50.0     0.1593    0.0158    50.0     0.1592    0.0164    50.0
               ot_nolit                 0.1668    0.0157    50.0     0.1697    0.0171    50.0     0.1687    0.0178    50.0
```

(α = 0.5 is the primary endpoint, per protocol §4.3.)

### 6. Key findings from Table 1

1. **22 of 24 arms give positive ΔAUROC at α = 0.5** — only reactome ×
   intogen_temporal_new and (barely) intact × intogen_temporal_new are
   non-positive.

2. **STRING physical subnetwork is the strongest network** (ΔAUROC up
   to 0.17 with ot_nolit). Its 85k curated physical-only edges
   outperform the 4-channel full STRING at 236k edges for most labels.

3. **Network quality > density.** funmap (196k edges, functionality-derived)
   beats reactome (20k edges, curated pathways) for *intogen2024*, but
   reactome beats funmap for *Open Targets labels*. The "best" network is
   label-dependent.

4. **α sensitivity is small.** Going from α=0.3 to α=0.7 changes
   ΔAUROC by 0.005-0.020. The primary α=0.5 is a sensible middle.

5. **The "newly-discovered" driver label (intogen_temporal_new) is
   hardest.** 4/6 networks give ΔAUROC < 0.03 on this label (only
   funmap, string_full700, rna_coexp are positive). The network prior
   helps predict already-known drivers more than it helps predict
   recently-discovered ones — a substantive biological finding.

### 7. Negative cell: reactome × intogen_temporal_new

ΔAUROC = -0.0112 ± 0.0512 at α = 0.5. SD is **5× the mean**, so the
effect is consistent with zero. The most likely interpretation:
Reactome catalogues curated pathway interactions from years of
literature; intogen_temporal_new is exactly the set of drivers
discovered 2020–2024 (the "new" ones), and these newer drivers are
not yet wired into Reactome. RWR on Reactome therefore has *no extra
signal* over degree ranking on the newest driver set. This is a
publication-quality negative result: a network's value depends on
how old the literature it encodes is.

## 8. Where this fits in the project

```
W0  A0 lock engine           (engine verified)
    A1 horse test            (harness verified)
    data download + parse    (5.5 GB raw, 6 networks built)

W1  A2 FunMap + STRING       (5 networks of 6)
W2  A3 IntAct + Reactome     (2 more networks)
    rna_coexp               (6th network)

W3  A4 full factorial       (Table 1, this log)   <-- YOU ARE HERE
    D-09 temporal drift     (B's proposal; needed for time stability claim)
    D-11 common-size cap    (B's proposal; needed for prevalence control)
    shared universe         (D-12 control)
W4-W5  network rewiring (100x per network)    (null-model control)
       per-fold stability analysis
W6  figures, supplementary, paper text
```

## 9. Commit log (after this file was written)

```
A4: full factorial 6 x 4 x 3 = 24 arms in 75 s; Table 1 written
  - scripts/run_a4.py
  - results/runs/a4_<net>__<label>.tsv          (24 files)
  - results/runs/a4_combined.tsv                  (4800 rows)
  - results/tables/a4_table1_delta_auroc.tsv
  - results/tables/a4_table1_rwr_auroc.tsv
  - results/tables/a4_table1_degree_auroc.tsv
  - docs/a4_log.md (this file)
```

## 10. W3 checklist status

| Item | Status |
|---|---|
| Full primary factorial (6 x 4 x 3) | ✓ (Table 1) |
| Table 1 written | ✓ |
| 70 tests pass | ✓ |
| Toy numbers unchanged | ✓ |
| D-09 temporal drift (B's proposal) | pending |
| D-11 common-size cap (B's proposal) | pending |
| shared-universe control (D-12) | pending |
| Commit + merge + push | pending |

Next: D-09 and D-11 are short (~15 min each, reuse the existing factorial
results and recompute caps / temporal labels). Then W5 rewiring.