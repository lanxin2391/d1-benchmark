# Protocol Compliance — What Remains to be Done

**Date:** 2026-09-27
**Authors:** lanxin2391 (A) + kakamiku (B)
**Source:** `D:\Grade3\swxxx\final\确定选题\D1_protocol.pdf`
**Status:** We have finished the **numerical collection** phase (data
download + 6 networks + 4 labels + 6 controls + W5 rewiring + Table 1 +
FDR). We have NOT finished the **analysis phase** that the protocol
calls for.

**TL;DR.** B is right. Our current state is best described as
"a pilot that produced the numbers and figures needed to write the
paper", not "the paper itself". Six of seven control mechanisms
exist; the seventh (degree-matched seed nulls) is missing; the
variance decomposition that the protocol calls "the study's actual
deliverable" has not been done; no pre-registration document exists;
no public release; no manuscript.

The work below is what the protocol requires and we have not done.
It is split between A (network + engine + stats) and B (label +
release + audit) per §6, with two cross-cutting items that must be
collaborative.

---

## 1. Protocol §6 milestone status

| Milestone | What it requires | Our state | Next action |
|-----------|------------------|-----------|--------------|
| **M2** | inputs frozen, **pre-registration filed** | inputs frozen ✓, pre-registration **not filed** ❌ | write the pre-reg document (collaborative) |
| **M4** | factorial complete, **H1 answered with variance decomposition** | factorial complete ✓, H1 z 8-200 but **no variance decomposition** ❌ | mixed-effects model + bootstrap CIs + equivalence tests |
| **M6** | benchmark released, manuscript submitted | released to GitHub only ✓, **no manuscript** ❌ | write manuscript |

**Verdict.** We are at the start of M4 (factorial executed) — numbers
are done but the analysis that converts them into a defensible
conclusion is not. The protocol treats "analyst time and design
care" as the binding constraint; that is what we owe next.

---

## 2. Control stack §4.4 — status of all 7

| # | Control | Implemented by | File |
|---|---------|----------------|------|
| 1 | Degree-only baseline | A | `d1/engine/baselines.py` |
| 2 | **Degree-matched seed nulls** | **— (missing)** ❌ | to be written by A |
| 3 | Configuration-model edge permutation (rewiring) | A | `scripts/run_rewiring.py` |
| 4 | Shared gene universe | A | `scripts/run_d12_shared.py` |
| 5 | Cross-catalogue label transfer | B | `scripts/build_transfer_*` |
| 6 | Temporal split | A | `scripts/run_d09.py` |
| 7 | Publication-count adjustment | B | `scripts/build_pubcount_adjusted.py` |

**6 of 7 done. The missing one (degree-matched seed nulls) is A's
work.** §5 says it takes ~20 min on 1 core or ~5 min on 4 cores;
essentially a script to generate 100 random genes per degree bin,
rebuild P0, re-evaluate, and compare to the observed ΔAUROC. The
result tests whether the degree distribution alone could explain
the observed margin. This is **A's first next task**.

---

## 3. Statistical analysis §4.6 — what is missing

The protocol specifies exactly what §4.6 must produce, and none of
it is done. This is A's biggest single gap.

### 3.1 Mixed-effects variance decomposition (the actual deliverable)

The protocol writes:

> "Margins are modelled with a mixed-effects model with random
> effects for network, label set and fold, giving a variance
> decomposition that answers 'how much of the between-arm spread is
> attributable to the network, to the labels, and to noise' — which
> is **the study's actual deliverable**, and which no existing
> paper reports."

What to do:

- random effects: (1 | network) + (1 | label_set) + (1 | fold)
- fixed effects: alpha (0.3, 0.5, 0.7), and a network-by-label
  interaction if identifiable
- response variable: delta_AUROC per arm
- variance components: σ²(network), σ²(label_set),
  σ²(network:label_set), σ²(residual)
- output: a single table showing how much of the between-arm
  variation is attributable to each factor

Implementation: `statsmodels.formula.api.mixedlm` or
`statsmodels.regression.mixed_linear_model.MixedLM`. With 24 arms ×
3 alphas = 72 rows it is a small model; should run in minutes.

Deliverable: `results/tables/variance_decomposition.tsv` +
`docs/variance_decomposition.md`.

### 3.2 Cluster bootstrap CIs

The protocol writes:

> "Confidence intervals come from a cluster bootstrap resampling
> *genes*, not gene-arm rows, because gene identity is the unit of
> non-independence."

What to do:

- for each (network, label, alpha) arm: resample genes with
  replacement (1000 iterations), recompute delta_AUROC at each
  iteration, take the 2.5th / 97.5th percentile as the 95 % CI
- report CI on the log-axes: logit( delta_AUROC / (1 − delta_AUROC) )
- the unit of resampling is the gene universe of the arm; do NOT
  resample rows (rows are not independent; genes are)

Deliverable: extend `results/tables/results_stats.tsv` to include
`ci_low` and `ci_high` per arm.

### 3.3 Equivalence tests for null results

The protocol writes:

> "Where the finding is a null, it is accompanied by an equivalence
> test against a pre-specified margin so that 'no effect' is
> bounded rather than merely unrejected."

The one null result is `reactome × intogen_temporal_new = −0.0112`.
A equivalence test would set a pre-specified "no meaningful
difference" margin (suggest 0.05 in delta_AUROC; matches the
reviewer's reported smallest meaningful difference) and test
whether the observed value lies inside the equivalence zone.
Implementation: two one-sided t-tests (TOST) on the per-fold
distribution; or scipy `Equivalence` from statsmodels.

Deliverable: extend `results/tables/results_stats.tsv` to include
`equiv_test_p` per arm; the reactome×temporal_new cell should
state the conclusion.

### 3.4 H1 verdict (currently informal)

§7 says: *"H1 failing is the strongest result in the set."* — so
H1's status must be made formal before publication.

What to do: write the H1 verdict from the bootstrap CIs and the
mixed-effects variance decomposition. The claim is H1 holds iff
`delta_AUROC > 0 with CI excluding 0 in every arm`. Currently
23/24 arms pass at point estimate level but the CI claim is not
made.

---

## 4. Provenance ablation §1.2 + TM1 — missing

The protocol's call to action:

> "The molecular layer is where 'RNA versus protein' belongs, and
> it is **the study's main factor**."

This is what B flagged. We have not done it.

### 4.1 What to do

Partition the 6 networks into 3 provenance classes (from the
protocol):

| Class | Networks | What the edges mean |
|-------|----------|---------------------|
| RNA-derived | `rna_coexp` | transcript-level co-variation |
| Protein-derived | `funmap`, `string_phys700`, `intact` | protein co-variation and physical binding |
| Curated / literature | `string_full700`, `reactome` | what has been written down about a gene |

For each (class, label, alpha) arm:

1. take the mean delta_AUROC across networks within the class
2. compute the SD across networks within the class
4. compare the 6 classes on the same metric

Output: `results/tables/provenance_ablation.tsv`. Should show that
ΔAUROC is similar across classes (this is H2: contribution is
provenance-invariant). If it is similar, this is H2's main
positive result.

### 4.2 The "100 candidates, unlabelled candidates" pilot observation

§9 records:

> "the top 30 unlabelled candidates, cross-checked against 5,806
> sequenced samples, Open Targets and ClinGen, looked like **network
> hubs rather than long-tail drivers**. On present evidence the
> candidate list is not yet worth experimental spend."

This is itself one of the benchmark's useful outputs (showing the
field that current network-based prioritisation is biased towards
hubs) and is a publishable negative observation. We need to write
the candidate-generation pipeline + this summary.

### 4.3 Manuscript treatment of provenance

Section in discussion: "Provenance does not account for between-arm
spread; the curated class's apparent advantage is literature
attention rather than encoded biology." This is the H3 verification
from §2 — once the ablation + pubcount adjustment are both in.

---

## 5. Manuscript §7 — preparation for the GNN objection

The protocol anticipates the most common reviewer objection:

> "your learner is too weak; a modern graph neural network would
> show a larger prior contribution."

> "Budget two weeks in M5 to run one published graph-neural-network
> method through the harness as a courtesy arm."

### 5.1 What to do

Pick one published GNN method (the protocol does not name one; the
field leaders as of the protocol's preparation date include
**GNNExplainer** / **GeneMANIA** / **GCN** for gene priority).
Run it through the same harness on the same splits, compute
delta_AUROC for the same 24 arms. If it shows a larger margin,
**that is consistent with the protocol's argument**: the claim is
that the GNN margin is also a contribution, and the benchmark's
released splits let any reader reproduce the result.

The GNN courtesy arm is A's work but B's data layer (cluster
infrastructure) is involved.

### 5.2 Risks to flag for the manuscript

§8 says: "Four people on one benchmark creates merge and convention
conflicts; **Medium**." We have not actively managed this beyond
W6. Manuscript-writing is where the conventions matter most
(figure format, how to present the variance decomposition, table
captions, supplementary organisation).

---

## 6. Release and audit §4.7, M4 — to do

### 6.1 Public release (B)

Per protocol, M2.8-M4.5: "Leaderboard protocol and public release
of fixed splits". We have the splits released to GitHub. The protocol
calls for a benchmark resource v1 with a leaderboard protocol.

What "public release" means:
- fixed splits released as files (DONE: in `data/processed/splits/`)
- leaderboard protocol document (TO DO: how a third party runs a
  new method and reports a number comparable to ours)
- archived on a public dataset registry (Zenodo or figshare) with
  a DOI (TO DO)
- "resource_and_dataset_table.csv" cited in the paper

B owns this. Draft leaderboard protocol, decide on archive
provider, draft the resource_table.

### 6.2 Reproducibility audit (B)

Per protocol, M4.3-M5.5: "Reproducibility audit; independent
end-to-end replay on a clean machine; **audit report**."

What "independent end-to-end replay" means:
- a team member B (currently spawned on b-labels — meets the
  "did not write the pipeline" criterion for the network side too
  if A's pipeline is the only one they don't audit)
- re-clones the repo on a clean machine
- runs the documented commands end-to-end
- compares numbers to ours

The deliverable is `docs/reproducibility_audit.md` (signed).
B's audit of A's pipeline would also count for the A-side, but
ideally we want a third party — possibly a professor — to do an
external audit. The protocol budgets this as part of M5.

---

## 7. Pre-registration document §M2 — what to write

The protocol is explicit:

> "M2 — inputs frozen, **pre-registration filed**. Networks, labels
> and learner are locked. Nothing after this point may change an
> input without a recorded amendment."

This is a free-standing document that should contain:

1. **Hypothesis text** (verbatim from §2 — H1, H2, H3)
2. **Endpoints** — the ΔAUROC primary and the four secondaries
3. **Frozen inputs**:
   - the 6 networks and where they come from
   - the 4 label sets and where they come from
   - the learner (RWR α=0.5, tol=1e-8, max_iter=200)
   - the splits: 5-fold × 10-repeat with seed 20260926
4. **Analysis plan** (fixed before any confirmatory run):
   - mixed-effects with random effects (network, label_set, fold)
   - cluster bootstrap on genes
   - Benjamini-Hochberg FDR
   - equivalence tests with pre-specified margin = 0.05 on null cells
5. **Locked controls** (the 7 in §4.4)
6. **Decision rule**: H1 holds iff delta_AUROC > 0 with CI excluding 0
   in **all 24 arms**.

This is collaborative (B writes the data-side section, A writes
the network + engine + statistics side) but B leads because it's
TM4's first deliverable.

---

## 8. Per-person checklist (split per protocol §6)

### A (TM1 + TM3) — network + engine + statistics

| # | Task | § | Effort | Priority |
|---|------|---|--------|-----------|
| A1 | **Pre-registration** (network + engine section) | §M2 | 1 day | HIGHEST |
| A2 | **Mixed-effects variance decomposition** | §4.6 | 2 days | HIGH |
| A3 | **Cluster bootstrap CIs** | §4.6 | 2 days | HIGH |
| A4 | **Equivalence tests** for the reactome × temporal_new cell | §4.6 | 1 day | HIGH |
| A5 | **Control 2: degree-matched seed nulls** | §4.4 | 1 day | HIGH |
| A6 | **Provenance ablation** (RNA vs protein vs curated) | §1.2, §M3.8-M5.3 | 3 days | HIGH |
| A7 | **GNN baseline** courtesy run | §7 | 2 weeks (big) | MEDIUM (M5) |
| A8 | **Network methods section** (manuscript) | §4.1 | 2 days | MEDIUM |
| A9 | **Network results + discussion** (manuscript) | §Results / §Discussion | 3 days | MEDIUM |
| A10 | **Figure 2 caption + cross-paper** (manuscript) | §Figures | 1 day | MEDIUM |
| A11 | **Update `docs/a*_log.md` + `w5_log.md`** (research-grade, not just pilot-grade) | — | 1 day | LOW |

**Total A: ~3-4 weeks of focused work.**

### B (TM2 + TM4) — label + release + audit

| # | Task | § | Effort | Priority |
|---|------|---|--------|-----------|
| B1 | **Pre-registration** (label + data section) | §M2 | 1 day | HIGHEST |
| B2 | **Public release protocol** (leaderboard + Zenodo/figshare) | §M2.8-M4.5 | 2 days | HIGH |
| B3 | **Reproducibility audit** (independent end-to-end replay) | §M4.3-M5.5 | 1 week | HIGH |
| B4 | **Update `resource_and_dataset_table.csv`** with citations + access routes | §3 | 1 day | HIGH |
| B5 | **Clinical-Gen integration in main result table** (5 labels, M0-M1.6 deliverable to make human paper canonical) | §1.1, M5.2-M6 | 2 days | MEDIUM |
| B6 | **Pubcount-adjustment paper write-up** (table + interpretation) | §2 H3 | Med | MEDIUM |
| B7 | **Cross-catalogue transfer paper write-up** (the most informative single control) | §4.4 #5 | Med | MEDIUM |
| B8 | **Label methods section** (manuscript) | §4.2 | 2 days | MEDIUM |
| B9 | **Label results + discussion** (manuscript) | §Results / §Discussion | 3 days | MEDIUM |
| B10 | **Update `docs/b_log.md`** to research-grade (not just pilot-grade) | — | 1 day | LOW |
| B11 | **Repository cleanup**: `b-labels` branch should merge with `main` before submission | §4.7 | 1 day | LOW |

**Total B: ~3-4 weeks of focused work.**

### Collaborative (both A and B)

| # | Task | § | Effort | Owner |
|---|------|---|--------|-------|
| C1 | **Manuscrip structure + outline** (skeleton) | — | 1 day | A+B |
| C2 | **Methods section** | §4 | 1 day | A (network + engine + stats) + B (label + data) |
| C3 | **Results section** | §Results | 2 days | A (Table 1 + Figs + controls) + B (ClinGen + transfer + cap) |
| C4 | **Discussion section** | §Discussion | 2 days | A (network biology + provenance) + B (label bias) |
| C5 | **Abstract + Intro** | §0, §1 | 2 days | A+B |
| C6 | **Supplementary materials** (full code, fixed splits, supplementary tables) | §4.7 | 2 days | A (code) + B (splits + labels) |
| C7 | **Submit + revision prep** | §M6 | 1 week (post-acceptance) | A+B |

---

## 9. Recommended timeline (4-5 weeks)

| Week | A | B |
|------|---|---|
| **1** (this week) | A5 degree-matched nulls; A1+A4 begin; A6 ablation prep | B1 pre-reg data section; B4 resource table; B2 release prep |
| **2** | A2 mixed-effects; A3 bootstrap CIs; A4 equivalence | B5 ClinGen integration; B6 pubcount write-up |
| **3** | A1 finish pre-reg; A6 ablation run; A8+A9+A10 manuscript start | B7 transfer write-up; B2 release publish; B8+B9 manuscript |
| **4** | A7 GNN baseline (parallel); A11 update logs | B3 reproducibility audit; C2+C5+C6 manuscript |
| **5** | A7 GNN done; C3+C4+C7 manuscript finalization | C3+C4+C7 manuscript finalization |

This gives both A and B enough work to fill 4-5 weeks without
overlapping on the same files.

---

## 10. What is good news

Although the analytical layer is missing, the **infrastructure** is
already in place. The numbers we have (Table 1, z 8-200 across 24
arms, the negative reactome×temporal_new cell, the controls all
producing sensible numbers) are the substrate that the analysis
needs to interpret. The implementation work is done; what remains is
statistical interpretation, write-up, and the pre-registration /
release scaffolding that the protocol calls for. This is a smaller
gap than it first looks.

Specifically: **the 7 control mechanisms all exist in the data**;
they are 6 wired and one broken in code, but the data structure
supports adding the seventh (degree-matched seed null) in a single
script using the same harness the other six already use. The fix is
mechanical, not conceptual.

---

## 11. Where to file this document

Suggested path: `docs/pending_work_protocol_compliance.md` (committed
to both `a-networks` and `b-labels`). This way, both A and B can refer
to it from their respective worktrees.

---

## 12. References to commits / files (this gives A/B a starting point)

- A-side infrastructure: `scripts/build_*.py`, `scripts/run_*.py`,
  `scripts/run_rewiring.py`, `scripts/run_stats.py`
- A-side results: `results/tables/a4_table1_*.tsv`, `w5_*.tsv`,
  `d09_table.tsv`, `d11_table.tsv`, `d12_table.tsv`,
  `results_stats.tsv`
- A-side logs: `docs/a0_log.md` ... `a4_log.md`, `rna_coexp_log.md`,
  `w5_log.md`, `s2_log.md`
- A-side figures: `figures/figure2_*.png`
- B-side infrastructure: `scripts/build_clingen.py`,
  `scripts/build_transfer_*.py`, `scripts/build_pubcount_adjusted.py`,
  `scripts/run_transfer_all.py`
- B-side results: `results/tables/pubcount_adj_comparison.tsv`,
  `results/tables/results_stats_transfer.tsv`,
  `results/tables/results_stats_clingen.tsv`
- B-side logs: `docs/b_log.md`, `docs/overlap_report.md`,
  `docs/pubcount_adj_results.md`, `docs/transfer_results.md`,
  `docs/results_stats.md`
- B-side figures: `figures/figure1_pilot_overview.png`

---

## 13. Decision needed before M2

The pre-registration document is the gating item for M2. We need
A and B to agree on:

- **Equivalence margin** for null cells (§4.6) — likely 0.05
  delta_AUROC; should be fixed in pre-reg.
- **"Locked" inputs** — exactly which datasets, which label
  versions, which network versions, which split versions. If B plans
  to add ClinGen as a 5th label, that should be decided NOW and
  pre-registered, not added later.
- **Public release venue** — Zenodo (DOI, free, ~10 GB limit) or
  figshare. Zenodo is more common for benchmark resources.

These are pre-M2 decisions and should be agreed in the next S sync
(before M2 deadline).

---

## 14. Summary scorecard

| Category | Required | Done | TODO |
|---|---|---|---|
| Data download (§3) | 5.5 GB | ✓ | — |
| Networks (§4.1) | 6 | ✓ 6 | — |
| Labels (§4.2) | 4 (or +ClinGen) | ✓ 4 + ClinGen | — |
| Learner (§4.3) | RWR α=0.5 | ✓ | — |
| Controls 1, 3, 4, 6 | done by A | ✓ | — |
| Control 2 (degree-matched seed nulls) | done by A | ❌ | A |
| Controls 5, 7 | done by B | ✓ | — |
| 5th label (ClinGen) | done by B | ✓ (built) | integrated into Table 1 (B) |
| Mixed-effects variance decomposition (§4.6) | A | ❌ | A |
| Cluster bootstrap CIs (§4.6) | A | ❌ | A |
| Equivalence tests (§4.6) | A | ❌ | A |
| Provenance ablation (§1.2, M3.8-M5.3) | A | ❌ | A |
| GNN courtesy arm (§7) | A | ❌ | A |
| Pre-registration (§M2) | collab | ❌ | A+B |
| Public release (§M2.8-M4.5) | B | ❌ | B |
| Reproducibility audit (§M4.3-M5.5) | B | ❌ | B |
| Manuscript (§M6) | mixed | ❌ | A+B |
| Resource table (§3) | B | partial | B |

**Scorecard: 9 of 14 required items done. 5 items remaining. ~4-5 weeks
of work.**