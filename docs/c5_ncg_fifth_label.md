# C5/C6 — the NCG catalogue as a fifth label, and what it does to the decomposition

**Date:** 2026-11-05
**Scripts:** `scripts/build_ncg_label.py`, `scripts/run_c5_ncg_arms.py`,
`scripts/run_c6_variance_decomp_v2.py`
**Tables:** `a4_table1_delta_auroc_v2.tsv` (30 combos),
`a4_combined_v2.tsv` (24,000 rows, 30 combos),
`d11_table_v2.tsv` (30 rows), `stats_variance_decomp_comparison.tsv`,
`stats_variance_decomp_v2.tsv`, `stats_variance_components_v2.tsv`

## The design decisions

| decision | what | why |
|---|---|---|
| which tier | `known` (708 genes) | the 104 `candidate` genes are multi-screen candidates, the 1,560 `screening` genes are not drivers |
| symbols | HGNC approved symbol, previous symbols mapped (9 of 708) | decision D-02 |
| rank | the number of publications NCG cites, descending, ties by symbol | NCG 6.0's own confidence level is unavailable (the portal is down), so the literature count stands in (D-13, adapted) |
| the cap | unchanged (N = min over the four canonical labels) | NCG has 407–694 positives per network, always above the temporal label's 63–150 → **the 24 released arms are not re-run** |
| the two conventions | the primary arm is uncapped (as the released 24 are), plus a D-11 capped arm | the cap is a *control*, not the primary design (`scripts/run_d11.py`) |
| family | pre-planned exploratory | plan item C5 pre-dates the analysis; the 24-arm pre-registered endpoint is untouched |

## The primary NCG arm (α = 0.5, 5 folds × 40 repeats)

| network | ΔAUROC | SD |
|---|---|---|
| funmap | +0.0771 | 0.0176 |
| rna_coexp | +0.0426 | 0.0169 |
| string_phys700 | +0.1352 | 0.0163 |
| intact | +0.0673 | 0.0089 |
| string_full700 | +0.1387 | 0.0133 |
| reactome | +0.0647 | 0.0216 |

All six positive, the network order is the one the other labels give,
and the margins sit between `intogen2024` and `ot_all` on most
networks — as the 65 % overlap with IntOGen predicts.

## The D-11 capped panel (five labels, equal positive-set size)

| network | intogen2024 | ot_all | ot_nolit | temporal | NCG |
|---|---|---|---|---|---|
| funmap | +0.1268 | +0.0191 | +0.0885 | +0.0196 | **+0.1127** |
| intact | +0.0681 | +0.0664 | +0.0973 | +0.0071 | **+0.0858** |
| reactome | −0.0110 | +0.0627 | +0.1649 | −0.0125 | **+0.0893** |
| rna_coexp | +0.0829 | +0.0403 | +0.0956 | +0.0127 | **+0.0720** |
| string_full700 | +0.1020 | +0.0785 | +0.1388 | +0.0604 | **+0.1177** |
| string_phys700 | +0.0955 | +0.0938 | +0.1819 | +0.0378 | **+0.1222** |

**The Reactome cell is label-specific.** Reactome is negative on
`intogen2024` and on the temporal label — the two most recent
catalogues — but positive on NCG, the 2019 catalogue. The Reactome
collapse is a statement about *when the drivers were curated*, not
about the network: that is the curator-lag hypothesis's prediction,
tested here without the temporal split.

## The variance decomposition (C6)

| variance source | % between-arm spread (24) | % between-arm spread (30) | change |
|---|---|---|---|
| network | 24.82 % | 26.95 % | +2.13 pp |
| label | 34.78 % | 30.75 % | **−4.02 pp** |
| network × label | 40.07 % | 41.95 % | +1.88 pp |
| fold (within cell) | 0.34 % | 0.35 % | +0.02 pp |

The label main effect shrinks and the interaction and the noise grow —
the direction the plan expected. But NCG is 65 % identical to IntOGen,
so this is the effect of adding a **redundant** label; an independent
fifth label would be expected to shrink the label effect more. The
substantive statement is that across five catalogues of different
provenance the label main effect is real but not dominant, and the
largest single piece of the between-arm spread is the interaction.

## Caveats

1. The primary NCG arm has a much larger positive set (407–694) than
   the temporal arm (63–150). That is the released design's convention,
   not a new inconsistency, and the D-11 panel is the control for it.
2. NCG is **not** an independent label (C2: Jaccard 0.522 with
   IntOGen), so the 30-arm grid measures the effect of a redundant
   label.
3. The NCG arms have no control-2 (degree-matched null) trace of their
   own yet; the fine-z table stays on the 24 arms.
