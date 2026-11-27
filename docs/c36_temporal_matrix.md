# Control 36 — the temporal transfer matrix

**Date:** 2026-11-26  
**Script:** `scripts/run_c36_temporal_matrix.py`  
**Tables:** `c36_temporal_matrix.tsv`, `c36_temporal_summary.tsv`

## The design

IntOGen has no 2021 or 2022 release, so the plan's per-year split cannot be built from it; what the catalogue offers is a four-generation chain (2019-11-12, 2020-02-01, 2023-05-31, 2024-09-20). For every (train year, test year) pair the seeds are the drivers known in the train-year release, the test positives are the drivers present in the test-year release but absent from the train-year release, and the endpoint is the delta-AUROC of the prior against the degree baseline on the same test genes (10 bootstrap draws of the seed set per cell).

## The mean delta-AUROC against degree, per network and pair

| network | 2019->2020 | 2019->2023 | 2019->2024 | 2020->2023 | 2020->2024 | 2023->2024 |
|---|---|---|---|---|---|---|
| funmap | +0.1053 | +0.0601 | +0.0484 | +0.0330 | +0.0353 | +0.0517 |
| rna_coexp | +0.0386 | +0.0316 | +0.0359 | +0.0353 | +0.0399 | +0.0305 |
| string_phys700 | +0.1453 | +0.1178 | +0.0971 | +0.1097 | +0.0873 | +0.0748 |
| intact | +0.0953 | +0.0559 | +0.0479 | +0.0428 | +0.0335 | +0.0371 |
| string_full700 | +0.1016 | +0.1067 | +0.1015 | +0.1029 | +0.1008 | +0.0766 |
| reactome | +0.0320 | +0.0503 | +0.0200 | +0.0453 | +0.0121 | -0.0005 |

## The decay curve (the median over the six networks)

| the pair | the gap in years | the median margin | the arms positive |
|---|---|---|---|
| 2019->2020 | 0.22 | +0.0984 | 6/6 |
| 2019->2023 | 3.55 | +0.0580 | 6/6 |
| 2019->2024 | 4.86 | +0.0482 | 6/6 |
| 2020->2023 | 3.33 | +0.0440 | 6/6 |
| 2020->2024 | 4.64 | +0.0376 | 6/6 |
| 2023->2024 | 1.30 | +0.0444 | 5/6 |

## Reading

## Caveats

1. The 2019->2020 pair is 83 days apart, so it is a consistency check rather than a decay step (the two catalogues are almost the same curation).
2. The test positives are the newly curated genes, which are few in each cell (the count is in the table); the AUROC on a few positives is noisier than the main endpoint, which is why each cell is repeated with bootstrap seed draws.
3. The degree baseline is the same baseline as the primary endpoint, so the numbers are comparable in sign and scale with Table 1.
4. **Exploratory** (amendment 4): no pre-registered endpoint is touched.
