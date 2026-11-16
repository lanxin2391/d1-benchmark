# Control 19 -- the baseline family sweep (#1)

**Date:** 2026-11-12  
**Script:** `scripts/run_baseline_family.py`  
**Tables:** `baseline_family_folds.tsv`, `baseline_family_arms.tsv`, `baseline_family_summary.tsv`

## Question

Our endpoint is `delta-AUROC = AUROC(RWR) - AUROC(degree)`. Two objections: (A) *degree is a strawman baseline*; (B) *RWR is a lucky special case*. Both are answered here on the same splits, the same folds and the same test genes.

## Scorers

| scorer | kind | definition |
|---|---|---|
| degree | baseline (control 1) | number of distinct neighbours |
| closeness | baseline | harmonic centrality, estimated from 256 random landmarks (the exact quantity needs O(V*E) all-pairs distances) |
| eigen | baseline | eigenvector centrality by power iteration |
| pagerank | baseline | PageRank, uniform teleportation, damping 0.85 |
| rwr05 | learner (primary) | RWR, restart 0.5 (locked) |
| rwr15 | learner (second damping) | RWR with restart 0.15 = the literature's default PageRank damping 0.85, outside our pre-registered sensitivity range of 0.3-0.7 |

## The baselines against degree (primary 24 arms)

| baseline | arms where it beats degree | median margin vs degree |
|---|---|---|
| closeness | 24/24 | +0.0312 |
| eigen | 16/24 | +0.0098 |
| pagerank | 1/24 | -0.1537 |

## The margins that matter

| comparison | positive arms | median | min | max |
|---|---|---|---|---|
| rwr05_vs_degree | 23/24 | +0.0809 | -0.0125 | +0.1711 |
| rwr15_vs_degree | 23/24 | +0.0571 | -0.0051 | +0.1552 |
| rwr05_vs_closeness | 20/24 | +0.0528 | -0.0245 | +0.1229 |
| rwr05_vs_eigen | 21/24 | +0.0613 | -0.0226 | +0.1959 |
| rwr05_vs_pagerank | 24/24 | +0.2040 | +0.0095 | +0.4147 |
| rwr05_vs_best_baseline | 19/24 | +0.0461 | -0.0245 | +0.1229 |

## Reading

- the strongest baseline per arm: closeness 15/24, eigen 9/24
- Spearman correlation between the RWR-vs-degree margin and the RWR-vs-closeness margin across the 24 primary arms: +0.833
- the RWR alpha = 0.5 margin recomputed here differs from Table 1 (`a4_table1_delta_auroc_v40.tsv`) by at most 4.4e-05

## Caveats

1. The landmark estimate of closeness is an approximation; on reactome (the smallest panel) the exact harmonic centrality agrees with it at Spearman +0.9963.
2. `rwr15` is the same family as `rwr05` with a different damping, not an independent learner; the independent-family check is Control 17 (the proximity learner).
3. Exploratory family (amendment 4): the pre-registered endpoint is untouched.
