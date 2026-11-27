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

## The read-out (10 repeats, alpha = 0.5, 16 arms)

All 16 arms are positive (the margin runs from +0.0157 to +0.1642). With
the new members, the **within-class** behaviour can finally be compared
with the released panel:

| label | the protein class (the median over its networks) | the RNA class |
|---|---|---|
| intogen2024 | +0.1099 | +0.0629 |
| intogen_temporal_new | +0.0646 | +0.0308 |
| ot_all | +0.1115 | +0.0718 |
| ot_nolit | +0.1212 | +0.0843 |

Two things follow, and they point the same way as the released panel.

1. **The class effect survives the extra networks.** The protein-derived
   networks have a larger margin than the RNA-derived ones in **all four
   labels**, which is the pattern the released panel already showed (the
   RNA network was the weakest in every label). It is now visible with
   four members in a class rather than with one.
2. **But the class effect is still not separable from the within-class
   spread.** The protein class's own spread inside a label is of the same
   size as the gap between the classes (for intogen2024: the protein class
   runs from +0.0871 to +0.1328, a spread of 0.046, against a gap of
   0.047 between the class medians). That is exactly the reason the
   paper's H2 is stated as *not testable* on six networks, and it is now
   visible with more members rather than less.

**So the collection does not rescue H2; it makes the reason visible.**
The next step for H2 is not more networks of the same kind but networks
that are matched on size and coverage *and* differ in provenance --
which is why the curated class still has only three members.

Note the size effect inside the protein class as well: string_full900
(100,383 edges) gives the largest margins (+0.1425 at the median) and
string_full400 (929,471 edges) a smaller one (+0.0744), so network
density moves the margin by more than the provenance class does. The
released panel's conclusion -- that provenance is confounded with size
and coverage -- is therefore confirmed, and now with a within-class
measurement that has four members in the class.

## Caveats

1. The new networks are **not** in the released panel yet; the release
   stays at the six released networks until the extra run is evaluated.
2. The new networks share their sources with the released ones (the same
   STRING download, the same TCGA expression), so they measure *within-
   source* variation, not independent replication.
3. 10 repeats rather than 40, because the arm count is four times
   larger.
