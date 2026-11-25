# Control 25 — mutation-first versus network-first

**Date:** 2026-11-18  
**Script:** `scripts/run_c25_mutation_vs_network.py`  
**Tables:** `c25_mutation_vs_network.tsv`, `c25_recurrence.tsv`, `c25_summary.tsv`

## The question

After the defence the question left open was whether the network prior adds anything to the **mutation** evidence. Two experiments, both on data already in the repository.

**25a, the head-to-head.** The answer key is the NCG-known catalogue, curated independently of IntOGen. On the same test genes and the same folds, four scorers are compared: the IntOGen 2024 mutation evidence (the number of cohorts in which the gene is a driver, and the number of samples), the network prior (RWR at alpha = 0.5, seeded with the training positives), degree, and the rank sum of the prior and the mutation score. Because IntOGen and the NCG overlap, the same comparison is run with the IntOGen driver genes removed from the test universe.

**25b, the recurrence comparison.** Among the held-out IntOGen 2024 positives, the prior's score is compared for the recurrently mutated drivers (seen in at least three samples in some cohort) and the single-mutation ones, with degree and the publication count as the comparators.

## 25a: the AUROC on the NCG-known label (the median over the folds)

| network | universe | the prior | degree | the mutation score (the cohorts) | the mutation score (the samples) | the prior + the mutation score |
|---|---|---|---|---|---|---|
| funmap | all | 0.6152 | 0.5370 | 0.8303 | 0.8290 | 0.8084 |
| intact | all | 0.7769 | 0.7112 | 0.8330 | 0.8319 | 0.8911 |
| reactome | all | 0.7185 | 0.6539 | 0.8596 | 0.8586 | 0.8734 |
| rna_coexp | all | 0.6190 | 0.5750 | 0.8446 | 0.8434 | 0.8119 |
| string_full700 | all | 0.8692 | 0.7314 | 0.8287 | 0.8279 | 0.9347 |
| string_phys700 | all | 0.7938 | 0.6584 | 0.8419 | 0.8412 | 0.9045 |
| funmap | minus_intogen | 0.5929 | 0.5418 | 0.5075 | 0.5076 | 0.5984 |
| intact | minus_intogen | 0.7205 | 0.6494 | 0.5156 | 0.5156 | 0.7297 |
| reactome | minus_intogen | 0.6647 | 0.5648 | 0.5114 | 0.5114 | 0.6672 |
| rna_coexp | minus_intogen | 0.5705 | 0.5435 | 0.5115 | 0.5115 | 0.5791 |
| string_full700 | minus_intogen | 0.8386 | 0.6688 | 0.5120 | 0.5120 | 0.8376 |
| string_phys700 | minus_intogen | 0.7498 | 0.5948 | 0.5148 | 0.5149 | 0.7519 |

## 25b: does the network know which drivers are recurrently mutated? (the AUROC for 'is recurrent' among the held-out positives)

| network | the prior | degree | the publication count |
|---|---|---|---|
| funmap | 0.5513 | 0.5031 | 0.5779 |
| intact | 0.5624 | 0.5444 | 0.5694 |
| reactome | 0.4949 | 0.5015 | 0.5771 |
| rna_coexp | 0.5916 | 0.5585 | 0.5417 |
| string_full700 | 0.5548 | 0.5234 | 0.5660 |
| string_phys700 | 0.5197 | 0.4812 | 0.5618 |

## Caveats

1. The mutation score is a static per-gene quantity; the prior is fold-dependent. Both are evaluated on the same test genes in the same folds.
2. The NCG catalogue and IntOGen overlap (the Jaccard is 0.522), so the `minus_intogen` universe is the one to read.
3. The recurrence status is the IntOGen annotation of the gene; a gene being recurrent is also a gene that has been looked at more, which is why the publication count is the comparator to read alongside.
4. **Exploratory** (amendment 4): nothing in the pre-registered endpoint is touched.
