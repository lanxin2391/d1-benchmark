# From the network to the experiment: an executable protocol

**Date:** 2026-11-26
**Status:** the plan's item 13 (the protocol) and item 14 (the clinical-use
boundary).
**Sources:** the released result tables (`results/tables/`) and the
manuscript sections F.12, F.21–F.34.

This document says what a network prior may be used for, how a candidate
list must be produced and reported, how any claim about the network's
contribution must be matched, what would count as experimental evidence
that a candidate matters, what would be needed before an allele-specific
network could be built, and where the clinical boundary lies. Every
threshold quoted here is a number this project measured, not a convention
imported from somewhere else.

---

## 1. What the network prior is for

**The network prior is a navigation aid: it reorders a candidate list. It
is not evidence that a gene matters.**

Three measurements fix that scope.

* On the curated driver catalogues the prior beats a degree baseline in 23
  of 24 arms (the median margin +0.0809) and beats the strongest
  classical baseline in 19 of 24 arms (+0.0461) (F.22). The increment is
  real under the protocol's own endpoint.
* On a functional answer key built from the DepMap assay (the essential
  genes as a fifth catalogue) the same protocol gives a positive margin in
  24 of 24 arms, with a median of +0.1236 (F.30) -- but the mechanism
  there is **propagation**: the essential genes are a module in every
  network, the seeds of a fold are other essential genes, and the margin
  tracks the *weakness of the degree baseline* rather than any discovery
  of function (the correlation between the margin and the degree
  baseline's own AUROC is -0.77).
* With the degree and the publication count regressed out of the score,
  the network's own component is above chance in 20 of 24 arms but modest
  (the median AUROC 0.5783 against the prior's 0.7051 and the publication
  count's 0.8072) (F.34).

So: **the network can be used to order a candidate list; it cannot be
used to assert that a candidate is functionally relevant, and it must
never be used as clinical evidence** (section 6).

---

## 2. The candidate-generation rule

**Input.** A curated positive set (the seeds) drawn from **at least two
independent curation sources**, with the source of each positive
recorded. One source is not enough: the four catalogues in this project
are one curation dimension, not four independent votes (the label factor
explains 97.8 % of the label main effect, F.23), and a candidate list
built from one source inherits that source's bias.

**Step 1 -- the score.** The RWR at the locked alpha = 0.5 on the
pre-registered six-network panel, with the seeds as the seed set. The
network panel and the alpha are fixed in advance; they are not chosen per
candidate set.

**Step 2 -- the aggregation.** Take the **median rank across the six
networks**, not the best network. The split design does not change the
result (the median margin moves by less than 0.004 across six designs,
F.25), so the panel is stable; but individual networks differ a lot in
how much they add (the median margin runs from +0.0469 to +0.1267 across
the six networks), so a "best network" choice would be a free parameter
fitted to the candidate set.

**Step 3 -- the covariates that must be reported, never used to select.**
Every candidate is reported with three numbers: its degree, its
publication count, and the number of catalogues that list it. They are
**not** selection criteria, and the reason is measured: the publication
count alone predicts the driver labels better than the network score does
(the median AUROC 0.8072 against the prior's 0.7051, F.34). A selection
rule that used the publication count would be a literature-curation
filter wearing a network's clothes.

**Step 4 -- the pathway check.** Report the candidate's pathway context,
and report whether that pathway is among the pathways the drivers are
enriched in. The reference numbers: for the IntOGen 2024 catalogue 257
of 4,294 Reactome pathways are enriched at q < 0.05, and for the
newest-driver catalogue **none** is (F.33). The consequence is a rule:
**for a catalogue with no significant pathway enrichment, the pathway
criterion is uninformative and must not be used to promote or demote a
candidate.**

**Step 5 -- the output.** A ranked list, each row carrying: the rank, the
median rank across the six networks, the degree, the publication count,
the number of catalogues listing the gene, the pathway context, and a
stability flag (the rank's movement under the seed perturbation of
section 4). A candidate list without those columns is not a result of
this project.

---

## 3. The degree-and-publication matching rule

Any statement of the form "the network found this" must be accompanied by
a comparison in which the two confounders are matched, per bait, at the
level of the individual gene -- not only at the level of the arm.

* **The publication count is matched per bait.** The machinery is the
  bait-matched null (F.12): each positive is compared with a
  degree-matched and publication-matched random gene drawn from the same
  network. The reference numbers: the median fraction of positives that
  beat the matched null is 11.6 % (the range 3.9 % to 33.1 % across the
  arms), against 17.6 % for the unmatched null. So roughly one positive
  in eight survives the match, and the rate varies threefold across
  arms; a claim that does not quote the matched rate is not reportable.
* **The degree is matched within a stratum**, not globally (the
  degree-matched null, F.14).
* **The effect size is reported against the matched null**, never
  against the unmatched one, and the unmatched value is given alongside so
  the size of the confound is visible.
* **The claim is scoped to the arms in which it holds.** The publication
  count explains the label better than the network does; that fact
  belongs in the sentence that reports the network's contribution, not in
  a footnote.

---

## 4. The perturbation-validation criteria

Two levels: a computational pre-check that costs nothing, and the
experimental criterion for the wet lab.

### 4.1 The computational pre-check (the seed-edit propagation test)

Delete a fixed fraction of the seeds, recompute the ranking, and measure
what happens to the margin (Control 8). The measured cost of the seed
loss on this panel: the median relative loss of the margin is **1.9 %**
when 10 % of the seeds are removed, **5.9 %** at 30 % and **11.7 %** at
50 %; the worst arm loses 11.6 %, 33.3 % and 29.2 % respectively.

**Rule.** A candidate set is *network-driven* only if its margin loses no
more than 12 % under a 10 % seed perturbation and the candidate's own
score CV is at or below the panel median (0.18). A set that fails this
is reported as **seed-fragile** and its candidates are not described as
network-derived.

### 4.2 The experimental criterion

For a candidate *c* with the predicted neighbourhood *N(c)*, the
pre-registered test is **comparative**, not threshold-based:

1. **Pre-registration.** The perturbation, the predicted direction of the
   change in every member of *N(c)*, and the readout are fixed before
   the experiment. A candidate set without a pre-registered prediction is
   an observation, not a test.
2. **The readout is functional** in the relevant cell line (a viability
   or rescue readout; the DepMap CERES score is the reference scale).
3. **The acceptance criterion is comparative.** The predicted set must
   change *more* than a matched control set of the same size and degree
   distribution, with the effect size and the sample size fixed in
   advance. A per-set threshold would be unfalsifiable and is not
   allowed.
4. **The negative control is mandatory.** A degree-matched candidate set
   must not show the same change. If it does, the effect is about the
   confound and not about the network.
5. **The stability requirement carries over.** The prediction must come
   from a seed set that passes section 4.1; a prediction from a
   seed-fragile set cannot be tested, because the experiment cannot
   distinguish it from noise.

---

## 5. The preconditions for an allele-specific network

An allele-specific network (the interactions of a specific mutant
protein, rather than of the gene) would let the method ask about
*variants* instead of genes. Four preconditions have to hold before such
a network is built, and each of them is a way the project has already
been burned:

1. **Independence from the label's curation.** The allele-specific
   interactions must not come from the same literature mining that
   produced the catalogues. If they do, the circularity of section 2
   returns in a new place, and the whole control stack would have to be
   rebuilt for the new network.
2. **Coverage of the same universe.** The allele-specific network must
   contain the same genes as the panel it is compared with (at least 70 %
   of them, and every gene that appears in the candidate list).
   Otherwise the comparison is between two different universes and the
   degree baseline is not comparable.
3. **A matched degree sequence.** The allele-specific network's degree
   distribution must be matched to the gene-level network's, so that the
   comparison is about topology rather than about hubness.
4. **A pre-registered endpoint.** The endpoint is the same
   delta-AUROC against the same degree baseline, on the same splits,
   registered before the run. The rewiring null (F.1, and the full grid
   under a second learner, F.29) and the label-permutation audit (F.28)
   are re-run on the new network, because they are what make a positive
   margin interpretable.

---

## 6. The clinical-use boundary

Three rules, stated so that they can be checked:

**Rule 1 -- the network ranks; it does not diagnose.** A network prior
may reorder a candidate list inside a research pipeline. It may not be
used as a diagnostic test, a risk score, or a basis for a clinical
decision, and no output of this project is a clinical result.

**Rule 2 -- no clinical claim without an outcome dataset.** This project
has no patient outcomes, no survival data and no treatment response, so
no statement about clinical utility can be derived from it, in either
direction. "The network helps" and "the network does not help" are both
unavailable here; what is available is the ranking-level statement of
section 1.

**Rule 3 -- the publication count is reported with every candidate.**
Because the publication count is the strongest single predictor of a
driver label (F.34), a candidate list that does not report it cannot be
interpreted: the reader cannot tell a network-driven candidate from a
well-studied one. This is the operational form of the disclosure the
whole project has been making about its own labels.

**The wording for the paper.** "The network prior is a research
instrument for ordering candidate lists. It is not a clinical marker, and
nothing in this work supports a clinical claim; the strongest single
predictor of a curated driver label in our own data is the number of
publications on the gene, which is exactly why we report it next to every
candidate."

---

## 7. What this protocol cannot do

* It cannot create biological evidence. It states what would count as
  evidence (section 4) and it is written so that a negative result can
  be recognised as one.
* It cannot be run on a catalogue with no independent curation. The
  pathway criterion, the matching rule and the stability check all need
  at least one independent source, and the newest-driver catalogue is the
  case where none of them can be applied (F.33: no pathway is enriched
  for the 152 newest drivers).
* It does not make the prior's margin larger. It makes the margin
  interpretable, and it makes the cases where the margin is not
  interpretable visible.
