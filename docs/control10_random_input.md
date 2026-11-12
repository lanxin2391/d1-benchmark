# Control 10 — Random-input degeneracy (Hill et al. 2019)

**Purpose.** Test whether the learner actually **uses** its input, or
whether it returns a fixed structural set (the hubs) regardless.

**Why it is needed.** Hill et al. (2019, *PLoS Computational Biology*
15(12):e1007403) benchmarked 17 network-contextualisation algorithms and
noticed that the algorithms differ sharply in how they behave when the
input is **meaningless**. They fed 10,000 random gene lists of length
200 and counted, for every node, how often it was ranked into the top
200. The result was algorithm-specific and sometimes extreme:

| algorithm | nodes in the top 200 of > 50 % of the random runs |
|---|---|
| ToppNet-HITS | 239 |
| SigNet | 200 |
| Causal Reasoning | 64 |
| InterConnectivity | 44 |
| Random Walk | 4 |
| Network Propagation | 0 |

The interpretation Hill et al. draw is the one a reviewer will draw
about D1: *an algorithm that returns (almost) the same nodes whatever
the input is not using the input at all — it is ranking a fixed
structural set, in practice the hubs.*

D1's control 2 (A5) answers a related question about the **margin**. It
does not measure the degeneracy of the **ranking**, which is the
quantity Hill et al. measured and the quantity a reader can check by
eye.

**Design.**

For every (network, label) arm, and for every one of the 200
(repeat, fold) evaluation folds, the locked learner is run **twice** on
the same fold universe:

* **real input** — the seeds are the training-fold positives;
* **random input** — the seeds are a random set of the **same size**,
  drawn from the network's genes **excluding the fold's own test genes**
  (the same exclusion that applies to the real seeds, so the random
  seeds are never the genes being scored).

Both regimes therefore yield 200 rankings of the same gene set, and
Hill et al.'s statistic is directly comparable between them.

**Results** (`results/tables/control10_random_input.tsv`,
`results/tables/control10_random_input_summary.tsv`):

| statistic | real input | random input |
|---|---|---|
| Hill degeneracy: genes in the top 100 of > 50 % of the 200 runs | **123.6** | **0.6** |
| mean degree percentile of the top 100 | 0.766 | 0.605 |
| mean overlap of the top 100 with the network's degree top 100 | 0.109 | 0.045 |
| mean pairwise rank Spearman between two runs of the same regime | 0.936 | 0.592 |

and the decisive one:

| input sensitivity (real-input ranking vs random-input ranking of the same fold) | value |
|---|---|
| mean Spearman | 0.592 |
| mean top-100 Jaccard overlap | **0.019** |

**Verdict.** **The learner uses its input.**

* Under a random input the RWR ranking collapses onto almost an
  unrelated gene list: the top-100 sets overlap by **1.9 %**, and only
  **0.6** genes persist in the top 100 across independent random runs —
  a degeneracy statistic in the range of the algorithms Hill et al.
  found to be well-behaved (Random Walk 4, Network Propagation 0) and
  three orders of magnitude below the degenerate ones (ToppNet-HITS
  239).
* Under the real input the ranking is highly reproducible (Spearman
  0.94 between folds) and it is *more* hub-concentrated than under a
  random input (degree percentile 0.766 vs 0.605; overlap with the
  network's degree top-100, 10.9 % vs 4.5 %). **The prior does not
  simply return the hubs** — it returns a different, more hub-enriched
  set than a meaningless input would, which is what a genuine,
  input-driven signal looks like.
* Because the degree baseline is the comparator, the margin is the
  part of that signal that degree does not already explain.

**Cost.** ~7 min for the whole panel on one core (5 random replicates
per arm).

**Reproduce.**

```bash
python scripts/run_control10_random_input.py --n-random 5
```

**Status.** Exploratory family (amendment 5). It does not enter the H1
decision rule; it is reported because the literature says the confound
matters.
