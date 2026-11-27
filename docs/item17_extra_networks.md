# Item 17 — more networks, so that every provenance class has four members

**Date:** 2026-11-26
**Script:** `scripts/build_extra_networks.py` (the collection),
`scripts/run_extra_networks_panel.py` (the protocol on the new networks)

## Why the class sizes matter

The H2 hypothesis ("do the networks of one provenance class behave
alike?") is not testable with two members in a class and three in
another: the paper's own conclusion was that six networks cannot
separate provenance from size and coverage. The panel therefore needs at
least four networks per class.

| class | before | after |
|---|---|---|
| protein-derived | 2 (string_phys700, string_full700) | **4** (+ string_full400, + string_full900) |
| curated-literature | 3 (funmap, reactome, intact) | 3 |
| RNA-derived | 1 (rna_coexp) | **3** (+ rna_coexp_d150k, + rna_coexp_d250k) |

The new members are **matched in size and coverage** to the panel by
construction: the two STRING cuts are the same source at the two
thresholds the builder already produced for the sensitivity arm (400 and
900, 929,471 and 100,383 edges against the released 236,712 at 700), and
the two co-expression variants are the same TCGA expression at a
different edge budget (149,330 and 249,390 edges against the released
196,172).

## What has been done

* the four new networks are built and live in `data/processed/networks/`
  with the same contract C2 format as the released ones;
* the splits for the new networks are built in a **new directory**
  (`data/processed/splits_extra/`, 10 repeats x 5 folds) — the released
  panels are not touched;
* the protocol is running on the four new networks x the four labels
  (alpha = 0.5), with the results in `results/tables/a4_table1_*_extra.tsv`.

## What the run will answer

The per-class consistency of the margin: with four protein-derived and
three RNA-derived networks, the question "is the provenance effect
separable from size and coverage?" becomes testable with the analysis the
paper already uses (the Kruskal-Wallis test across classes, and the
per-network margin spread within a class). The released panel's six
networks stay as they are until this run is read out.

## Caveats

1. The new networks are **not** in the released panel yet; the release
   stays at the six released networks until the extra run is evaluated.
2. The new networks share their sources with the released ones (the same
   STRING download, the same TCGA expression), so they measure *within-
   source* variation, not independent replication.
3. 10 repeats rather than 40, because the arm count is four times
   larger.
