# Control 35 — a real GraphSAGE control

**Date:** 2026-11-26  
**Script:** `scripts/run_c35_graphsage.py`  
**Tables:** `c35_graphsage.tsv`, `c35_graphsage_summary.tsv`

## What this is

The existing GNN arm is a self-written Node2Vec with p = q = 1, (the DeepWalk limit), and it was found to have four protocol defects (amendment 3). This control implements GraphSAGE properly (Hamilton, Ying & Leskovec, 2017) -- two layers, mean aggregation, a learnable self-feature and a learnable neighbour feature, an identity embedding table (a gene network has no node features) -- and evaluates it with the **same ranking rule** as the other GNN arm: the cosine between a gene's embedding and the centroid of the training positives' embeddings.

The objective is supervised on the training fold only (the embeddings are pulled towards the seeds' centroid and pushed away from the training negatives), so this arm is a *supervised* GNN: the strongest sensible neural baseline for the same task, and a fair comparison with the diffusion prior rather than a structure-only one.

## The result

| network | label | GraphSAGE | degree | the GraphSAGE margin | the RWR prior's margin |
|---|---|---|---|---|---|
| funmap | intogen2024 | 0.5626 | 0.5080 | +0.0546 | +0.1001 |
| funmap | intogen_temporal_new | 0.5474 | 0.5015 | +0.0460 | +0.0196 |
| funmap | ot_all | 0.6001 | 0.5760 | +0.0241 | +0.0853 |
| funmap | ot_nolit | 0.5787 | 0.5712 | +0.0075 | +0.0780 |
| intact | intogen2024 | 0.7080 | 0.7191 | -0.0110 | +0.0651 |
| intact | intogen_temporal_new | 0.6145 | 0.6388 | -0.0243 | +0.0071 |
| intact | ot_all | 0.6912 | 0.7068 | -0.0156 | +0.0842 |
| intact | ot_nolit | 0.6594 | 0.6733 | -0.0138 | +0.0784 |
| reactome | intogen2024 | 0.6748 | 0.6765 | -0.0018 | +0.0429 |
| reactome | intogen_temporal_new | 0.5442 | 0.5667 | -0.0225 | -0.0125 |
| reactome | ot_all | 0.6877 | 0.6969 | -0.0092 | +0.1103 |
| reactome | ot_nolit | 0.6873 | 0.6850 | +0.0022 | +0.0996 |
| rna_coexp | intogen2024 | 0.6207 | 0.5913 | +0.0294 | +0.0521 |
| rna_coexp | intogen_temporal_new | 0.5457 | 0.5285 | +0.0172 | +0.0127 |
| rna_coexp | ot_all | 0.5535 | 0.5339 | +0.0196 | +0.0767 |
| rna_coexp | ot_nolit | 0.5667 | 0.5439 | +0.0228 | +0.0834 |
| string_full700 | intogen2024 | 0.7355 | 0.7138 | +0.0217 | +0.1133 |
| string_full700 | intogen_temporal_new | 0.6251 | 0.6275 | -0.0023 | +0.0604 |
| string_full700 | ot_all | 0.7752 | 0.7903 | -0.0152 | +0.1201 |
| string_full700 | ot_nolit | 0.7540 | 0.7544 | -0.0003 | +0.1267 |
| string_phys700 | intogen2024 | 0.6693 | 0.6536 | +0.0157 | +0.1189 |
| string_phys700 | intogen_temporal_new | 0.5717 | 0.5875 | -0.0158 | +0.0378 |
| string_phys700 | ot_all | 0.7121 | 0.6901 | +0.0220 | +0.1582 |
| string_phys700 | ot_nolit | 0.6569 | 0.6395 | +0.0174 | +0.1711 |

## Reading

## Caveats

1. 10 folds per arm (sampled from the 200 folds of the 40-repeat design, evenly spaced), because each fold is a full GraphSAGE training run; the RWR column is the released 200-fold mean.
2. The GraphSAGE arm needs PyTorch (2.3.1+cpu here); the rest of the benchmark does not, so the arm is self-contained and its results are reproducible with `pip install torch`.
3. The dense adjacency is used (the networks have at most 17,918 nodes), which is why the arm is a CPU job of a few seconds per fold.
4. **Exploratory** (amendment 4): no pre-registered endpoint is touched.
