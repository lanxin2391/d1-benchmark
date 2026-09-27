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
|| D-BASE | 2026-09-26 | B | ``base_seed = 20260926`` for ``StratifiedKFold(random_state=base_seed + r)`` in ``d1.splits``. Freeze with this value; any change breaks split-file reproducibility. | Today's date in yyyymmdd form; it is a stable but unambiguous choice for Phase 1. |
|| D-RANK | 2026-09-26 | B | Driver-gene label file ``rank`` = descending count of distinct CANCER_TYPEs in which the gene is a driver; ties broken by ascending SYMBOL. | Simplest publication-friendly proxy for "how widely supported as a driver"; aligns with how the catalogue's own tier column is constructed. |
|| D-ENV | 2026-09-26 | B | Use ``uv`` + ``.venv-d1/`` on Lenovo's Windows PC; ``environment.yml`` keeps a conda-compatible recipe (with the TUNA conda-forge mirror added) for any teammate who can reach ``conda.anaconda.org``. | ``conda.anaconda.org`` returns ``HTTP 000 CONNECTION FAILED`` from this network; ``uv pip install`` against ``https://pypi.tuna.tsinghua.edu.cn/simple/`` finishes in under a minute with all 27 protocol-required packages. |
|| D-LOG10 | 2026-09-26 | B | Control-7 publication-count regressor uses ``log10(n_pubmed + 1)`` (not ``log10(n_pubmed)``) so genes with no publications stay at 0 instead of getting a negative value. | Matches the protocol's intuition of "genes that have never been written about" being a meaningful stratum, not a shifted-out outlier. |
|| D-OT-AREA | 2026-09-26 | B | Open Targets cancer ontology root = ``MONDO_0045024`` ("cancer or benign tumor"); a disease is considered "cancer" if ``MONDO_0045024`` appears in its ``therapeuticAreas`` list OR in its ``ancestors`` list. | Matches the protocol's §1.1 ontology choice; 3,661 diseases pass this filter in 26.06. |
|| D-OT-NOLIT-SCORE | 2026-09-26 | B | ot_nolit per-(target, disease) score = normalised harmonic sum ``H(s) = sum_i s_i / i^2 / (pi^2/6)`` over the per-datasource scores with ``europepmc`` dropped. | Single datasource collapse: the sum is the standard Open Targets-style combiner; the ``pi^2/6`` denominator (Basel) keeps the score in roughly [0, 1]. |
|| D-OT-CAP | 2026-09-26 | B | ot_all and ot_nolit are each capped at the top N = 600 genes by score. Recompute N at S2 once the shared-universe size (D-11) is known. | Matches the protocol's §1.1 "~600 genes" figure for the Open Targets label panel. |

## S2 proposals (B-side drafts, to be discussed at S2 sync point)

| ID | Date | Author | Decision | Rationale |
|---|---|---|---|---|
| ~~D-09 (proposal)~~ → **D-09** | 2026-09-26 | B (draft) → A+B (locked at S2) | **Temporal split label** = `intogen_temporal_new` (152 genes): IntOGen drivers added in the 2020→2024 window (i.e. present in 2024 release but absent from 2020 release). Use this set as the *training-set temporal holdout* — train ranking on `intogen2024` minus `intogen_temporal_new`, evaluate AUROC uplift on `intogen_temporal_new` alone. | One single canonical "control 6" matched to the protocol's "newly-discovered drivers" intuition; alternative (date-of-publication on every driver) is not provided by IntOGen. |
| ~~D-11 (proposal B)~~ → superseded | — | — | (Proposal B — global min — was discussed, not chosen.) | — |
| ~~D-11 (B's preferred)~~ → **D-11** | 2026-09-26 | B (draft) → A+B (locked at S2) | **Common size N = per-network min** — for each (label × network) arm separately, N = min(positives_in_LCC) across the four labels for that network. Result: funmap N=374, string_full400 N=595, string_full700 N=587, string_full900 N=519, string_phys700 N=519. | Keeps each label's network-specific coverage intact while still equalising positive prevalence per-arm; matches protocol §4.2's "common-size cap" intuition. |
| ~~D-09-NUMPY (proposal)~~ → **D-NUMPY** | 2026-09-26 | B (draft) → A+B (locked at S2) | **Pin numpy to `>=1.26,<2.0` in both `environment.lock.yml` (conda) and `pyproject.toml` / `requirements.txt` (uv).** A's reproducer runs on 1.26.4 (works); B's venv runs on 2.4.6 (also works, bit-identical AUROC in 50-fold smoke test). Empirical evidence: 0.0998 ± 0.019 vs 0.0998 ± 0.019 — no detectable numerical drift. | A flags that `np.linalg.solve` crashed on Windows with numpy 2.x; B never hit it because the test path doesn't touch `linalg.solve`. Pin the floor at 1.26 to keep A safe; allow future 2.x bumps only after re-running all 61 tests + the 4-arm smoke. |

Reserved IDs (planned in the Phase 1 work protocol, to be filled when decided):
D-08 RNA network tumour-only vs pooled samples
D-09 temporal split definition (locked: intogen_temporal_new, see above)
D-10 Open Targets no-literature scoring (filled: D-OT-NOLIT-SCORE)
D-11 common label size N (locked: per-network min, see above)
D-NUMPY numpy pin (locked: >=1.26,<2.0, see above)
D-BASE stratified-KFold base seed (filled: 20260926)
D-RANK driver-gene rank definition (filled: distinct cancer-type count desc, symbol asc)
D-ENV environment provisioning (filled: uv + .venv-d1; conda mirror = TUNA)
D-LOG10 publication-count regression (filled: log10(n+1))
D-OT-AREA cancer disease ontology (filled: MONDO_0045024 + ancestors)
D-OT-NOLIT-SCORE no-literature score combiner (filled: harmonic sum / pi^2/6)
D-OT-CAP Open Targets label cap (filled: 600; re-cap at S2 to per-network min via D-11)
