# Network notes: what raw FunMap and STRING look like (W0 section 7.2 + A2 input)

This file is the starting point of task A2 in W1. Anyone reading it should be
able to reproduce the parser scripts without re-reading the FunMap / STRING
documentation.

All commands assume the repo root is `D:\Grade3\swxxx\final\D1\d1-benchmark`
and the `d1` conda env is active.

---

## 1. Raw file layout

```
data/raw/networks/
├── funmap/
│   └── funmap.tsv                                 (3.4 MB, ~196,800 edges)
└── string/
    ├── 9606.protein.info.v12.0.txt.gz             (19,699 proteins, ENSP -> symbol)
    ├── 9606.protein.aliases.v12.0.txt.gz         (19,699 x N aliases, NOT used by A2)
    ├── 9606.protein.links.detailed.v12.0.txt.gz   (13.7 M edges, 7 channels + combined)
    └── 9606.protein.physical.links.v12.0.txt.gz  (1.48 M edges, physical only)
```

## 2. FunMap `funmap.tsv`

- Headerless, two columns of HGNC symbols (e.g. `TP53 MDM2`).
- 196,800 rows = 196,800 edges (it is a simple edge list, no header line).
- No weights. No NaN. No self-loops in the file (we verified).
- `pandas.read_csv(..., header=None, names=["gene_a", "gene_b"])` is the
  correct call; letting pandas auto-infer the header treats the first data
  row as column names and silently drops one edge.

```python
>>> edges = pd.read_csv("data/raw/networks/funmap/funmap.tsv", sep="\t",
                        header=None, names=["gene_a", "gene_b"])
>>> edges.shape
(196800, 2)
>>> edges.iloc[:3]
   gene_a    gene_b
0  C3AR1     CLEC7A
1  DHRS7B    HSD17B12
2  DDX17     HNRNPDL
>>> len(set(edges.gene_a) | set(edges.gene_b))
10525
```

FunMap's edge list is **already in HGNC symbols** — no ID mapping is needed
(decision D-02). The build script `scripts/build_funmap.py` just feeds the
edge list straight into `finalize_edges`.

## 3. STRING `9606.protein.info.v12.0.txt.gz`

| Column | Meaning |
|---|---|
| `#string_protein_id` | ENSP id with species prefix, e.g. `9606.ENSP00000000233` |
| `preferred_name` | The HGNC symbol STRING considers canonical for this protein (mostly matches HGNC approved symbol) |
| `protein_size` | Sequence length, not used |
| `annotation` | Free text, not used |

- 19,699 proteins, all unique preferred_name, 0 NA.
- All IDs start with `9606.` (human).
- Mapping is **1-to-1**: every ENSP has exactly one preferred_name (no NA).
- A2 maps `ENSP -> symbol` by stripping the `9606.` prefix.

## 4. STRING `9606.protein.physical.links.v12.0.txt.gz`

| Column | Meaning |
|---|---|
| `protein1` | ENSP id with species prefix |
| `protein2` | ENSP id with species prefix |
| `combined_score` | Aggregate confidence 150–999; **higher = stronger evidence** |

- 1,477,610 edges (physical interactions only).
- Score thresholds: 400 = medium, 700 = high (primary, D-03), 900 = highest.
- A2 builds `string_phys700` from this file with threshold 700.

## 5. STRING `9606.protein.links.detailed.v12.0.txt.gz`

Same schema as the physical file, plus 7 per-channel scores:

| Column | What it measures |
|---|---|
| `neighborhood` | conserved gene neighbourhood across species |
| `fusion` | gene-fusion events across species |
| `cooccurence` | phylogenetic co-occurrence of homologs |
| `coexpression` | correlated expression across experiments |
| `experimental` | biochemical / genetic experimental evidence |
| `database` | curated pathway / complex database annotations |
| `textmining` | text mining of scientific literature |
| `combined_score` | Bayesian integration of all 7 channels |

- 13,715,404 edges.
- A2 builds `string_full400`, `string_full700`, `string_full900` from this
  file with thresholds 400 / 700 / 900. The protocol makes 700 the primary;
  400 and 900 are sensitivity arms.

## 6. Why three thresholds for STRING full

`combined_score` is a number, but it has no absolute meaning ("high" in STRING
is not the same as "high" in IntAct). The protocol pre-specifies three
thresholds so we can show how ΔAUROC changes as we drop low-confidence
edges:

| Threshold | Edges after finalize | LCC nodes |
|---|---|---|
| 400 (medium)  | 929,471 | 19,486 |
| 700 (high, **primary**) | 236,712 | 15,882 |
| 900 (highest) | 100,383 | 11,693 |

If ΔAUROC is stable across all three, the result is robust. If ΔAUROC
crashes at 900, the result depends on weak edges that we cannot trust.

## 7. Why exactly these networks and not e.g. `string_full500`

The protocol §4.1 lists `string_full400`, `string_full700`, `string_full900`
as the three STRING full arms. `string_full500` and `string_full600` are
NOT in the design — adding more thresholds would inflate the multiple-
testing correction without adding insight (they sit between 400 and 900,
neither boundary).

## 8. Decisions touched by A2

A2 implements D-01, D-02, D-03 from `docs/decisions.md`:

- **D-01 (unweighted)**: edges passing the threshold are kept; the score
  column is dropped in the output TSV. The threshold is the only place the
  score is used.
- **D-02 (HGNC symbols)**: STRING ENSPs are mapped via `preferred_name`;
  unmapped proteins silently drop their edges.
- **D-03 (STRING physical >= 700)**: same threshold as full, applied to the
  physical subnetwork.

## 9. Quick sanity checks

After `python scripts/build_funmap.py && python scripts/build_string.py`,
the user can confirm:

```powershell
pytest                                       # 27 passed
python scripts/run_toy_benchmark.py          # toy numbers unchanged
python scratch/smoke_networks.py              # each built network loads + RWR runs
```

`scratch/smoke_networks.py` loads each of the 5 networks, builds W^T, runs
one RWR with 5 random seeds at alpha=0.5, asserts mass conservation. All
5 networks pass in ~5 s total.