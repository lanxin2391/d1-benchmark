# Decision log

Every choice that the protocol (D1_protocol.pdf) does not fix goes here. Never delete or rewrite
an entry: to change a decision, add a new dated entry that supersedes it. At M2 this file becomes
the analysis-details section of the pre-registration.

Format: `ID | date | who | decision | reason`

| ID | Date | Who | Decision | Reason |
|---|---|---|---|---|
| D-01 | 2026-09-26 | A+B | Primary analysis uses **unweighted** edges in every network. Weighted RWR may be a sensitivity arm. | Weights mean different things in different networks (FunMap probability, STRING score, correlation). |
| D-02 | 2026-09-26 | A+B | Gene symbols: approved HGNC symbol; a previous symbol is used only if it points to exactly one approved symbol; an alias only if unique and not itself an approved symbol. Ambiguous IDs are dropped and counted. | Avoid silently merging different genes. |
| D-03 | 2026-09-26 | A+B | STRING physical subnetwork uses the same combined-score threshold as STRING full (≥ 700). | Keep the two STRING arms comparable. |
| D-04 | 2026-09-26 | A+B | IntAct: human–human pairs (both taxid 9606), UniProt on both sides, isoform suffix removed, no MI-score threshold in the primary arm; MI-score ≥ 0.45 as optional sensitivity arm. | Protocol gives no threshold. |
| D-05 | 2026-09-26 | A+B | RNA co-expression: Pearson per cohort, Fisher z averaged across cohorts, top 196,800 pairs by mean z (density matched to FunMap), then LCC. | Protocol §4.1: match density, not threshold. |
| D-06 | 2026-09-26 | A+B | Seeds = training-fold positives; they are not in the ranked test set. Test set = test-fold positives + test-fold negatives. | Protocol §4.3; avoids scoring the seeds themselves. |
| D-07 | 2026-09-26 | A+B | Precision-at-k tie-break: score descending, then gene symbol ascending. | Deterministic results. |
| D-12 | 2026-09-26 | A (proposed, confirm with B at S0) | Shared-universe control (control 4): RWR propagates on each network's **full LCC**; seeds and evaluation genes are restricted to the shared gene set; degree is the degree in the full LCC. | An induced subgraph on the shared genes can fall apart into pieces and changes topology; restricting seeds + evaluation removes the coverage difference, which is what control 4 targets. |

Reserved IDs (planned in the Phase 1 work protocol, to be filled when decided):
D-08 RNA network tumour-only vs pooled samples · D-09 temporal split definition ·
D-10 Open Targets no-literature scoring · D-11 common label size N
