# A2 log: build FunMap + STRING networks

**Task:** W1, task A2 (protocol section 6.2 of W1).
**Author:** lanxin2391 (Person A)
**Date:** 2026-09-26
**Outcome:** 5 networks built, 27 tests still pass, toy benchmark unchanged.

This is the first task where Person A writes code that touches real data.
A0 and A1 were about locking and verifying the engine / harness. A2 is
about feeding the engine real protein-protein interaction networks.

---

## 1. What A2 is

The benchmark needs 6 real networks. A2 builds the first 5 of them (FunMap
protein network + 4 STRING variants). The 6th, `rna_coexp`, is built in
W2 from `funmap_input_expression_data.tgz` once we figure out the cohort
choice.

All 5 networks are produced in contract C2 format
(`gene_a<TAB>gene_b<TAB>weight`, undirected, no self-loops, no duplicates,
LCC only) and have a contract C2b `.meta.json` next to them.

## 2. Step-by-step reproduction

```powershell
conda activate d1
cd D:\Grade3\swxxx\final\D1\d1-benchmark

# 1. Look at raw inputs first
code data\raw\networks\funmap\funmap.tsv       # 196,800 rows, headerless
code data\raw\networks\string\                 # 4 .gz files

# 2. Build the 5 networks
python scripts\build_funmap.py
python scripts\build_string.py
# Outputs: data/processed/networks/{funmap,string_phys700,string_full400,string_full700,string_full900}.tsv + .meta.json

# 3. Verify each one passes contract C2
python scratch\verify_networks.py
# Loads each network via load_network, checks symmetric / LCC / correct edge count.

# 4. Smoke test: RWR actually runs on each
python scratch\smoke_networks.py
# 5 networks in ~5 s total.

# 5. Make sure nothing regressed
pytest
# Expected: 27 passed

python scripts\run_toy_benchmark.py
# Expected: toy numbers unchanged (still +0.292 for toy_module alpha=0.5).
```

## 3. Outputs produced

| File | Size | Edges | LCC nodes |
|---|---|---|---|
| `data/processed/networks/funmap.tsv`         | 3.4 MB  | 196,605 | 10,189 |
| `data/processed/networks/string_phys700.tsv` | 1.5 MB  |  85,576 |  9,830 |
| `data/processed/networks/string_full400.tsv` | 16 MB   | 929,471 | 19,486 |
| `data/processed/networks/string_full700.tsv` | 4.1 MB  | 236,712 | 15,882 |
| `data/processed/networks/string_full900.tsv` | 1.7 MB  | 100,383 | 11,693 |

FunMap is the densest (most edges per node). STRING at 700 is the
protocol's primary arm. STRING at 400 has ~6x more edges; STRING at 900
~2x fewer.

## 4. Why these 5 networks cover the W1 question

The benchmark asks: does network topology help identify cancer drivers,
and how does that depend on the network? We probe this by varying two
axes on STRING:

- **Type** (full vs physical): full integrates text mining and coexpression;
  physical only uses experimental binding. If the two agree, we trust
  topology; if they disagree, the result depends on what STRING counts as
  "interaction".
- **Threshold** (400/700/900): the protocol's main sensitivity arm. We
  want to know if the conclusion survives dropping weak edges.

FunMap is **independent of STRING** (it comes from a different algorithm,
functionality-based), so it acts as a third reference point.

## 5. What the script does that the harness does not

`finalize_edges` in `d1/networks/io.py` already does all the cleaning
(drop unmapped, drop self-loops, order alphabetically, dedupe, LCC). A2's
contribution is the **STRING-specific preprocessing**:

1. Read STRING's raw files (which use `9606.ENSP...` IDs and a
   `combined_score` column).
2. Map ENSP -> HGNC symbol via `9606.protein.info.v12.0.txt.gz`.
3. Apply the score threshold to drop low-confidence edges.
4. Hand the result to `finalize_edges`.

A2 is intentionally minimal: 60 lines per build script, all the clever
logic stays in `finalize_edges`.

## 6. Lessons for A3 (IntAct) and A4 (Reactome)

- Reactome is already a tab-delimited edge list with HGNC symbols on both
  sides, so `build_reactome.py` will be ~20 lines: load + finalize.
- IntAct uses UniProt IDs (`P53` instead of `TP53`) and needs taxid
  filtering (D-04: both partners must be taxid 9606). It is closer to
  the STRING build in structure (ID mapping + threshold cut).

## 7. Smoke test details (`scratch/smoke_networks.py`)

For each of the 5 networks, the script:

1. Loads via `d1.networks.io.load_network`
2. Builds the transition matrix via `transition_T`
3. Picks 5 random nodes as seeds, runs RWR with alpha=0.5
4. Asserts mass conservation (sum of P ≈ 1)
5. Asserts non-negative probability
6. Times the whole pipeline

All 5 networks succeed in 0.2–2.7 s. The 2.7 s outlier is `string_full400`
(929 k edges, 19 k nodes); the rest are well under 1 s. This is the first
real-data check that the RWR engine scales as expected.

## 8. Day-3-to-W1 checklist status

| Item | Status |
|---|---|
| 6.1 Read `evaluate.py` and `io.py` | ✓ (A1) |
| 6.2 Toy numbers match Appendix B | ✓ (A1) |
| 6.3 Open `toy_runs.tsv` in pandas | ✓ (A1) |
| 6.4 Write `docs/toy_results.md` | ✓ (A1) |
| 7.1 Data into `data/raw/` | ✓ (W0 last step) |
| 7.2 First look at FunMap + STRING | ✓ (`docs/network_notes.md`) |
| A2 FunMap + STRING parsers | ✓ (this task) |
| Commit + merge + push | pending |
| A3 IntAct + Reactome | next (W2) |
| A4 wire 6 networks into harness | next (W3) |

## 9. Commit log (after this file was written)

```
A2: build FunMap + STRING networks (5 of 6)
  - scripts/build_funmap.py
  - scripts/build_string.py
  - docs/network_notes.md
  - docs/a2_log.md (this file)
  - data/processed/networks/funmap.tsv + meta.json
  - data/processed/networks/string_{phys700,full400,full700,full900}.tsv + meta.json
```

(The processed networks ARE committed: 25 MB total, fits GitHub
without LFS. This means anyone who clones the repo can run real-data
analyses without re-running the build scripts. To skip the heavy
build step, edit `.gitignore` to add `data/processed/networks/` and
delete the tracked files with `git rm -r --cached data/processed/networks/`.)