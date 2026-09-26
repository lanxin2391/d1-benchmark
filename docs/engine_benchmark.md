# Engine benchmark (A0)

Output of `python scripts/bench_rwr_speed.py --rewire` on each team PC. Paste yours below.

## Reference run (cloud test machine, 2 CPUs, 26 Sep 2026)

```
network: 10525 genes, 196800 edges
sparse transition matrix: 4.8 MB (dense would be 0.89 GB)
RWR alpha=0.5: 11 iterations, 114 ms per batch of 50, 2.29 ms per seed vector  (protocol: 6.2 ms, 23 iterations)
projected real arms (6 networks x 4 labels x 50 folds): 2.7 s
one rewiring (10 x |E| = 1968000 swaps): 0.44 min, degrees preserved: True  (protocol: 0.7 min)
projected 100 rewirings x 6 networks on 1 core: 4.4 h
```

Note: the synthetic scale-free graph mixes faster than the real FunMap network, which is modular,
so it needs fewer iterations (11 vs the protocol's 23). Expect about 23 on real FunMap.

## Person A's PC

```
(paste here)
```

## Person B's PC

```
(paste here)
```
