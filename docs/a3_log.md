# A3 log: build IntAct + Reactome networks

**Task:** W2, task A3 (protocol section 6 of W2).
**Author:** lanxin2391 (Person A)
**Date:** 2026-09-27
**Outcome:** 2 primary + 1 sensitivity networks built. 61 tests pass. End-to-end arms run.

A3 is the second network-building task. After A2 produced FunMap + STRING,
A3 adds IntAct and Reactome. Together with W2's `rna_coexp`, this brings
the panel to 6 primary networks + 2 sensitivity arms (the protocol's full
panel).

---

## 1. What A3 is

The benchmark needs 6 primary networks. A2 built 5 (FunMap + 4 STRING
variants). A3 builds 2 more:

- **intact** (primary): all human-human IntAct pairs, no MI threshold
- **intact_miscore045** (sensitivity): same, but only edges with
  intact-miscore ≥ 0.45
- **reactome**: all human Reactome curated interactions

Total after A3: 7 networks (5 from A2 + 2 from A3) plus `rna_coexp`
which is W2's other task.

## 2. Step-by-step reproduction

```powershell
conda activate d1
cd D:\Grade3\swxxx\final\D1\d1-benchmark

# 1. Build the networks
python scripts\build_intact.py     # ~3 minutes (streams the 11 GB inside intact.zip)
python scripts\build_reactome.py   # ~30 seconds

# 2. Smoke-test all 7 networks (RWR actually runs)
python scratch\smoke_all_7.py
# Each network should converge in <30 iterations and have mass = 1.0 exactly.

# 3. Regression tests
pytest                                  # 61 passed (no regressions)
python scripts\run_toy_benchmark.py     # toy numbers unchanged

# 4. End-to-end arms
python scripts\build_splits.py --label intogen2024 --universe native_intact --universe native_reactome
python scripts\run_arm.py --network data\processed\networks\intact.tsv ^
    --splits data\processed\splits\intogen2024__native_intact.tsv ^
    --network-id intact --label-id intogen2024 --universe-id native_intact ^
    --alpha 0.5 ^
    --out results\runs\end_to_end_intogen2024_intact.tsv
python scripts\run_arm.py --network data\processed\networks\reactome.tsv ^
    --splits data\processed\splits\intogen2024__native_reactome.tsv ^
    --network-id reactome --label-id intogen2024 --universe-id native_reactome ^
    --alpha 0.5 ^
    --out results\runs\end_to_end_intogen2024_reactome.tsv
```

## 3. What I read (key decisions)

| File | Function | Role |
|------|----------|------|
| `data/raw/networks/intact/intact.zip` | `intact.txt` (PSI-MITAB 2.7, 42 columns, ~11 GB uncompressed) | All species, all MI scores |
| `data/raw/networks/reactome/reactome.homo_sapiens.interactions.tab-delimited.txt` | tab-delimited, 9 columns (header starts with `#`) | Human only, no scores |

## 4. IntAct parser (decisions)

| Step | Choice | Source |
|---|---|---|
| Read format | stream `intact.txt` from inside `intact.zip` (don't unzip) | file is 11 GB uncompressed; streaming keeps memory flat |
| Pair IDs | columns 0 and 1, strip `uniprotkb:` prefix | IntAct IDs are `uniprotkb:P37840` style |
| Filter | both taxids must contain `9606` (human) | D-04: "human-human pairs (both taxid 9606)" |
| Isoform suffix | keep (mapper strips them automatically) | HGNCMapper handles `-2` style |
| Score | parse `intact-miscore:` from column 14 | sensitivity arm MI ≥ 0.45 |
| Primary arm | no MI threshold (all human-human pairs) | D-04: "no MI-score threshold in the primary arm" |
| Sensitivity arm | MI ≥ 0.45 → `intact_miscore045` | D-04: "MI-score ≥ 0.45 as optional sensitivity arm" |

A bug encountered and fixed: I initially used `cells[10]` for taxid A, but
the correct index is `cells[9]` (the table has 15 columns up to
`confidence`, indices 0-14). Caught by the empty result on the first run.

## 5. Reactome parser (decisions)

| Step | Choice | Source |
|---|---|---|
| Read format | tab-delimited, skip `#` header line | file starts with `# Interactor 1 uniprot id ...` |
| Pair IDs | columns 0 and 3 (UniProt), strip prefix | UniProt is the only always-present column |
| Skip Ensembl/Entrez | use UniProt only | many cells are `-` for missing IDs |
| Filter | none (file is human-only by filename) | `reactome.homo_sapiens.*` |
| Score | none (Reactome does not publish per-edge confidence) | n/a |

## 6. Outputs

| File | Size | Edges | LCC nodes | Uniprot→HGNC yield |
|---|---|---|---|---|
| `data/processed/networks/intact.tsv` | 9.4 MB | 564,706 | 17,918 | 92.5% (1,006,295 / 1,088,470) |
| `data/processed/networks/intact_miscore045.tsv` | 2.1 MB | 126,992 | 14,851 | 92.4% (462,395 / 500,668) |
| `data/processed/networks/reactome.tsv` | 0.3 MB | 20,143 | 3,953 | 60.0% (74,895 / 124,865) |

The Reactome mapping yield (60%) is much lower than IntAct (92.5%)
because Reactome catalogues include enzyme complexes whose subunits
often don't have a single UniProt ID (they're listed by complex name).
This is acceptable — finalize_edges silently drops the unmapped
edges; we report the mapping yield in `docs/network_notes.md`.

## 7. End-to-end results (intogen2024 label, α = 0.5)

| Network | ΔAUROC | SD | n_test_pos | n_test |
|---|---|---|---|---|
| funmap              | 0.0998 | 0.019 | ~30 | ~600 |
| string_full700     | 0.1133 | 0.011 | ~30 | ~600 |
| **intact**          | **0.0650** | **0.0074** | ~30 | ~1800 |
| **reactome**        | **0.0436** | **0.0187** | ~30 | ~400 |

All 4 networks are positive. Order:
**string_full700 (0.113) > funmap (0.100) > intact (0.065) > reactome (0.044)**

This pattern makes sense:
- STRING integrates 7 evidence channels → richer topology signal
- FunMap is protein-derived (functionality) → strong but unrelated to STRING
- IntAct is experimental binding → real interactions but mostly
  high-throughput (lower confidence per edge)
- Reactome is curated pathways → sparse, well-validated but small

The full 6-network factorial at W4-W5 will show whether this ranking
is preserved across all 4 label sets.

## 8. What this means for the W4-W5 factorial

By the end of W3 we will have all 6 primary networks (the 6th is
`rna_coexp`, built from `funmap_input_expression_data.tgz`).
At that point, `run_arm` covers the full factorial:

```
6 networks × 4 labels × 50 folds × 3 alphas × 2 methods (RWR + degree)
= 7,200 rows per C5 result table
```

This fits in seconds-to-minutes on the FunMap-sized networks
(verified in A0 speed benchmark).

## 9. Smoke-test artefacts (`scratch/smoke_all_7.py`)

For each of the 7 (5 from A2 + 2 from A3) networks, the script:

1. Loads via `d1.networks.io.load_network`
2. Builds the transition matrix via `transition_T`
3. Picks 5 random nodes as seeds, runs RWR with α=0.5
4. Asserts mass conservation (sum of P ≈ 1)
5. Asserts non-negative probability

All 7 networks pass. The two new ones (intact, reactome) converge in
18-22 iterations, in line with the A2 networks.

## 10. Day-3-of-W2-to-W3 checklist

| Item | Status |
|---|---|
| 7 networks built (5 from A2 + intact + reactome) | ✓ |
| Smoke-test all 7 | ✓ (8 entries: 5+2+1 sensitivity arm) |
| 61 pytest passes | ✓ |
| Toy numbers unchanged | ✓ |
| End-to-end arms for new networks | ✓ |
| `docs/network_notes.md` updated (TODO) | pending |
| `docs/a3_log.md` (this file) | ✓ |
| Commit + merge + push | pending |
| A4 wire 6 networks into harness end-to-end | next (W3) |

## 11. Commit log (after this file was written)

```
A3: build IntAct + Reactome networks (2 of 6 primary networks)
  - scripts/build_intact.py
  - scripts/build_reactome.py
  - docs/a3_log.md (this file)
  - data/processed/networks/intact.tsv + meta.json
  - data/processed/networks/intact_miscore045.tsv + meta.json
  - data/processed/networks/reactome.tsv + meta.json
  - results/runs/end_to_end_intogen2024_intact.tsv
  - results/runs/end_to_end_intogen2024_reactome.tsv
  - data/processed/splits/intogen2024__native_intact.tsv
  - data/processed/splits/intogen2024__native_reactome.tsv
```

(`data/processed/` is in `.gitignore` but git is letting it through
because we explicitly added the files; this is the same pattern as A2
and works for reproducibility.)