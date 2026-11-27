# Control 37 — HotNet2: a feasibility result with one arm

**Date:** 2026-11-26
**Script:** `scripts/run_c37_hotnet2_port.py`
**Table:** `results/tables/c37_hotnet2.tsv`
**Porting notes:** `docs/hotnet2_porting_notes.md`

## Why this arm exists

HotNet2 (Deuther et al., *Nature Genetics* 2017) is the standard
path-based method in the cancer literature, so "is it only the diffusion
prior that works?" deserves a head-to-head. The arm is deliberately built
so that the mutation evidence is **independent of the labels**: the
mutation matrix is the **MC3 TCGA MAF** (the pan-cancer somatic mutation
aggregation) and the network is one of our six, while the labels are the
literature-curated driver lists. The HotNet2 score therefore cannot be
"the label re-read" — it is a genuine mutation-first strategy.

## What was achieved

* the upstream repository (`raphael-group/hotnet2`) cloned, ported with
  2to3 and patched for networkx 3 / h5py 3 / numpy 2 (the patch list is
  in `docs/hotnet2_porting_notes.md` and in the script, so the port is
  reproducible);
* the MC3 MAF parsed into HotNet2's SNV format (one line per tumour
  sample) and the heat file built (19,687 genes with a mutation heat);
* the network files built for the five networks that fit in memory: the
  influence matrix (the personalised PageRank) and the permuted edge
  lists with their PPR matrices, 20 permutations each;
* the driver run end to end on the **Reactome** network, producing the
  significant subnetworks — the core subnetway contains TP53, TP53BP2,
  CDKN2A, CDKN2A-adjacent genes and the expansion subnetway contains PTEN,
  AXIN2 and TNKS (the expected biology, which is a good sign that the
  port is faithful).

## The result (Reactome, the mutation matrix from the MC3 MAF)

| label | the HotNet2 AUROC | the degree AUROC | the margin |
|---|---|---|---|
| intogen2024 | 0.5343 | 0.6802 | **-0.1460** |
| ot_all | 0.4848 | 0.7032 | **-0.2184** |
| ot_nolit | 0.4849 | 0.6815 | **-0.1965** |
| intogen_temporal_new | 0.5476 | 0.5988 | **-0.0513** |

**Reading.** In this configuration the mutation-first strategy does **not**
beat the degree baseline on the driver labels: the HotNet2 score ranks
746 of the network's 3,953 genes (those in a significant subnetwork) and
the rest at zero, which gives an AUROC near chance, while the degree
baseline is far above it. The direction is the opposite of the diffusion
prior's, and it is the same direction as the mutation-evidence arm of
Control 25 (where the mutation score dominated *because* it was the label's
own curation source; here the mutations are independent, and an
independent mutation burden does not rank literature-curated drivers).

That is a useful result for the paper's framing: **a mutation-based
priority that is independent of the catalogue does not recover the
catalogue's positives, while a network prior does** (the same arm's
network prior is +0.110 on Reactome × ot_all). It is also a statement
about the label, not about the biology: the driver lists are curation
products, and independent somatic mutations are a different thing.

## What is not done, and what the caveats are

* **One network only.** The other four failed on an HDF5 key mismatch in
  the ported writer, and `intact` is out of reach in any case (the
  similarity matrix is 17,918 squared, which the NumPy fallback cannot
  hold; the C extension would be needed).
* **20 permutations, not the paper's 100**, so the significance test is
  under-powered and the selected subnetworks are noisier than the
  published setting.
* **The score is the heat of the best subnetwork a gene belongs to**, not
  the upstream's per-gene significance (the significance file is empty at
  this permutation count). A gene in no significant subnetwork scores zero.
* The C and Fortran extensions are not compiled, so the similarity matrix
  falls back to the NumPy path (the same result, much slower: the Reactome
  scoring run took 20 s, the `intact` one would not have finished).

**So this is a feasibility result with one arm, not a head-to-head over
the panel, and it is not in the paper.** Plan item 9 stays open: the port
has to be frozen as a patch file, the HDF5 key mismatch fixed, and the
arm run over the panel with the paper's permutation count.
