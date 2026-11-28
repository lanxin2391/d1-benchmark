# HotNet2: the porting state

**Date:** 2026-11-26
**Scope:** the plan's item 9 (the HotNet2 head-to-head).

## What was achieved

* the upstream repository (`raphael-group/hotnet2`) cloned and ported with
  2to3;
* the **MC3 TCGA MAF** parsed into HotNet2's SNV format (one line per
  tumour sample, the mutated genes in the columns) -- this is the point
  of the arm: the mutation evidence is *independent* of the labels, so
  the HotNet2 score is a genuine mutation-first strategy and not the
  label re-read;
* the **heat file** built (19,687 genes with a mutation heat);
* the **network files** built: the influence matrix (the personalised
  PageRank over the network) plus the permuted edge lists and their PPR
  matrices -- the expensive part of the pipeline -- in about a minute per
  network for Reactome.

The pipeline lives in `scripts/run_c37_hotnet2.py`; the clone is
patched by the script's own `patch_compat` step, so a fresh clone is
fixed the same way.

## The patches that were needed

| # | the problem | the fix |
|---|---|---|
| 1 | the Python 2 print statements | 2to3 |
| 2 | `nx.to_numpy_matrix` / `nx.from_numpy_matrix` (removed in networkx 3) | `to_numpy_array` / `from_numpy_array` |
| 3 | `mp.cpu_count` is a method in Python 3 | `mp.cpu_count()` |
| 4 | h5py 3 cannot write numpy's unicode dtype, and networkx 3 hands over an EdgeView rather than a list | `save_hdf5` writes every non-numeric value as byte strings |
| 5 | h5py 2 read a dataset's value with `.value` | `f[key][()]` |
| 6 | the upstream code joins the output prefix with the output directory | the prefix has to be *relative* |
| 7 | the network name read back from the HDF5 file is a byte string | decode it before it becomes a dictionary key |

## What still blocks the arm

The main driver (`HotNet2.py`) still contains Python 2 idioms that
surface one at a time: the decoded gene names have to be threaded
through the heat filter, the similarity matrix and the output writer.
Each remaining fix costs one pipeline cycle (the inputs are cached, so a
cycle is a few minutes), and the patching was done by hand on a clone.

**So the head-to-head is feasible on this machine -- the expensive part
runs -- but it is not finished and it is not in the paper.** Item 9 stays
open. When it is finished, the port has to be frozen as a patch file so
that the arm is reproducible without re-porting by hand.

## What it would take

1. finish the driver port (the heat filter's gene-name handling, the
   similarity matrix's index mapping, the output writer's column names);
2. freeze the port as `patches/hotnet2-python3.patch` in the repository;
3. run the six networks (about a minute each for the network files, plus
   the scoring) with the paper's 100 permutations, not the 5 used in the
   pilot;
4. evaluate the scores against the four labels on the same endpoint as
   every other arm.

The arm's design is already fixed and independent of the porting: the
mutation matrix is the MC3 MAF, the network is one of our six, and the
endpoint is the delta-AUROC against the degree baseline on the label's
positives.

## NetMix2 (the successor of HotNet2) — 2026-11-27

* the paper: **Ruan et al., *Journal of Computational Biology* 2022, DOI
  10.1089/cmb.2022.0336** (PMID 36525308, PMC9917315). The full text is
  in `data/papers/Ruan_2022_NetMix2.html` (and the plain-text extraction
  next to it); the publisher's PDF is behind a bot check.
* the code: `raphael-group/netmix2` — **Python 3** (no porting needed) and
  **no O(n^2) matrix** (which is why it can run on the panel's
  networks, where HotNet2 could not).
* the pipeline (`scripts/run_c37_netmix2.py`): the MC3 TCGA MAF → the
  per-gene mutation heat → the z-scores NetMix2's local-FDR step wants
  (one score file per network, covering every gene of that network) → the
  upstream driver → the per-gene score → the same AUROC-against-degree
  endpoint as every other arm.
* the state: the inputs are built and the driver runs into the local-FDR
  step, where the upstream's z-score handling produces NaNs for the genes
  it does not find in the score file. That is the last mile, and it is a
  debugging pass rather than a porting project.
