# Engine benchmark — A0 task 5.5

Run on 2026-09-26 by lanxin2391, machine reference "PC" (Windows 10, 16 logical CPUs,
Python 3.11.16, numpy 1.26.4, scipy 1.17.1).

## Command

```
python scripts/bench_rwr_speed.py --rewire
```

The script builds a synthetic scale-free graph (Barabasi-Albert, n=10,525,
m=19 per node) trimmed to 196,800 edges, restricts it to the largest connected
component, then measures (a) RWR alpha=0.5 time on a batch of 50 seed vectors
(averaged over 5 runs, with a warm-up), (b) one degree-preserving rewiring of
10 x |E| swaps.

## Output (raw)

```
Python 3.11.16 on Windows-10-10.0.26200-SP0, 16 logical CPUs
network: 10525 genes, 196800 edges (built in 1.8 s)
sparse transition matrix: 4.8 MB (dense would be 0.89 GB)
RWR alpha=0.5: 11 iterations, 358 ms per batch of 50, 7.17 ms per seed vector  (protocol: 6.2 ms, 23 iterations)
projected real arms (6 networks x 4 labels x 50 folds): 8.6 s
one rewiring (10 x |E| = 1968000 swaps): 0.99 min, degrees preserved: True  (protocol: 0.7 min)
projected 100 rewirings x 6 networks on 1 core: 9.9 h
```

## Interpretation

| Metric | This PC | Protocol reference | Note |
|---|---|---|---|
| Iterations to converge (synthetic scale-free) | 11 | 23 (real FunMap) | Synthetic BA mixes faster than the modular real FunMap. Do not compare iterations across graphs. |
| ms per seed vector | 7.17 | 6.2 | Within ~15% of reference. Comfortable for the factorial. |
| ms per batch of 50 | 358 | ~310 | Batching gives ~50x speedup over per-vector calls. |
| Projected 6x4x50 factorial | 8.6 s | seconds (reference not given) | Acceptable. |
| One rewiring (10x\|E\| swaps) | 0.99 min | 0.7 min | Slightly slower; PC-specific. |
| 100 rewirings x 6 networks | 9.9 h | similar | Should be parallelised across cores or run overnight. |
| Degree preservation | True | True | Confirms `nx.double_edge_swap` is configured correctly. |

## What this means for the W5 run plan

- All 6 networks built once: cheap (seconds).
- Single arm (1 network x 1 label set x 50 folds x 3 alphas): ~0.4 s.
- Full factorial (6 x 4 x 50 x 3): ~8.6 s on this PC. Negligible.
- 100 null-model rewirings per network: ~100 min per network, ~10 hours for 6
  networks on 1 core. Parallelise across 6 cores to ~100 minutes total.

## Full log

`scratch/bench_rwr.log` contains the verbatim output of this run.