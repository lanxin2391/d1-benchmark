# Control 8 — Seed-degree stratification (Zhang et al. 2021)

**Purpose.** Test whether the primary endpoint (ΔAUROC = AUROC(RWR) −
AUROC(degree)) depends on the **degree of the seed genes**.

**Why it is needed.** Zhang et al. (2021, *Briefings in Bioinformatics*
22(5):bbab006) benchmarked seven network-based gene-prioritisation
methods for cerebral small vessel disease and reported the sharpest
published warning against the RWR family:

> "performance of RWR-like algorithms would appear to be better in LOOCV just by having higher degree centralities for the seed nodes, which does not necessarily reflect true ability for the algorithm to predict new candidate genes."

Their random-seed experiment (1,000 trials × 10 random seeds)
correlated a method's median rank with the degree centrality of its
seeds.

D1's control 2 (A5) already contains a *degree-matched random-seed
null*, which asks a weaker question: is the margin larger than what
degree-matched random seeds would give? It does not ask whether the
**real** seeds being hubs is what produces the margin. Control 8 asks
that directly, inside D1's own design.

**Design.**

For each (network, label) arm at α = 0.5, on the 40-repeat split panel
(`data/processed/splits_v40/`):

1. For every (repeat, fold), take the training-fold positives and split
   them into three equal **degree terciles** by their degree in that
   network. Ties are broken by the index, so the strata are always the
   same size even in a scale-free network where many genes share a
   degree.
2. Propagate **three** seed matrices in one batched RWR call — one per
   tercile. RWR is linear in the seed weights, so the three strata cost
   about the same as one.
3. Evaluate each stratum's ranking with the **same** test fold, the
   **same** degree baseline and the same AUROC routine as the primary
   endpoint. The three numbers are therefore directly comparable.

**Results** (`results/tables/control8_seed_degree.tsv`,
`results/tables/control8_seed_degree_summary.tsv`):

| seed-degree stratum | mean ΔAUROC over the 24 arms | arms with a positive margin |
|---|---|---|
| low (lowest third) | **+0.0527** | 22/24 |
| mid | **+0.0609** | 24/24 |
| high (highest third) | **+0.0537** | 23/24 |

Mean(high − low) = **+0.0010**, i.e. one fiftieth of the margin itself.
The sign of the gap is not consistent across the panel (it is negative in
11 of 24 arms).

**Verdict.** The margin **does not depend on the degree of the seeds.**
If the driver catalogues were degree-enriched and RWR were simply
rewarding hubs, the high-degree stratum would carry a much larger
margin; it does not.

The two arms whose *low*-stratum margin is negative are `funmap ×
intogen_temporal_new` (−0.0037) and `reactome × intogen_temporal_new`
(−0.0044), both of which have a small full-catalogue margin to begin
with (0.0196 and −0.0125). The pattern is consistent with "these
labels carry little signal in these networks", not with a systematic
hub dependence.

**Cost.** ~4 min for the whole panel on one core.

**Reproduce.**

```bash
python scripts/run_control8_seed_degree.py
```

**Status.** Exploratory family (amendment 5). It does not enter the H1
decision rule; it is reported because the literature says the confound
matters and because the result is informative whichever way it went.
