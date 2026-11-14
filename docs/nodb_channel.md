# Control 18 -- STRING without the `database` channel

**Date:** 2026-11-12  
**Script:** `scripts/build_string_nodb.py`  
**Table:** `results/tables/c18_nodb_channel.tsv`

## Why

Reactome is the gold standard of FunMap's edges, Reactome feeds STRING's `database` channel, and Reactome- and CGC-derived evidence is also what the labels are built from. A large part of what makes STRING useful here may therefore be the same literature knowledge the labels encode.

## How (and the approximation)

The edge weight is the **maximum over the non-database channels** (neighbourhood, fusion, co-occurrence, co-expression, experimental, text-mining) and the edges are cut to the same edge count as the primary arm (236714 edges against 236712), so the two networks are density-matched. **This is not STRING's own combined score** (a Bayesian integration of the seven channels) and is labelled as an approximation; the question it answers is qualitative -- does the margin survive when the database channel is gone?

| label | with the database channel | without it | change |
|---|---|---|---|
| intogen2024 | +0.1133 | +0.0884 | -0.0249 |
| intogen_temporal_new | +0.0604 | +0.0661 | +0.0057 |
| ot_all | +0.1201 | +0.0820 | -0.0381 |
| ot_nolit | +0.1267 | +0.1370 | +0.0103 |

## Caveats

1. The approximation: the maximum over channels is not STRING's Bayesian combination, so the absolute weights differ; the comparison is about the *presence* of the database channel, not about STRING's scoring.
2. The node universe changes with the edge set, so a margin difference mixes the channel effect with the universe effect.
3. Exploratory family (amendment 4).
