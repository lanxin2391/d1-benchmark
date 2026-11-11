# D1 benchmark

**What a network prior contributes to cancer-driver prioritisation: a leakage-controlled benchmark.**

**Release:** v1.0-rc3 (tag `v1.0-rc3`) · Zenodo DOI pending ·
pre-registration: [Zenodo 10.5281/zenodo.23121762](https://doi.org/10.5281/zenodo.23121762)
**Audit:** 50 PASS / 0 FAIL on `python scripts/audit_reproducibility.py`
**Manuscript verification:** 57 PASS / 0 FAIL on `python scripts/verify_claims.py`

A fixed learner (random walk with restart, RWR), a panel of six gene
networks and four driver-gene label sets. The primary endpoint is the
margin of RWR over a degree-only ranking
(ΔAUROC = AUROC(RWR) − AUROC(degree)) under identical cross-validation folds.
Design: `D1_protocol.pdf`. Work plan: `D1_parallel_work_protocol.docx`.
**Manususcript:** `docs/research_report_draft.md` (SCI-standard draft)
plus the three section files it is assembled from.

| Person | Line | Owns code in |
|---|---|---|
| A | Networks & Engine (TM1 + TM3) | `d1/engine/`, `d1/networks/`, `scripts/build_*network*`, `scripts/run_*` |
| B | Labels & Data (TM2 + TM4) | `d1/hgnc.py`, `d1/splits.py`, `d1/labels/`, `d1/data/` |

Rule: never edit a file the other person owns. Ask them to change the
script and regenerate.

---

## Verifying the release

Three mechanisms, all of which must pass before a release:

```bash
python -m pytest                                  # 131 passed, 3 skipped
python scripts/audit_reproducibility.py            # 50 PASS / 0 FAIL
python scripts/verify_claims.py                    # 57 PASS / 0 FAIL
python scripts/verify_report.py                    # 67 PASS / 0 FAIL
python scripts/build_release_manifest.py --release v1.0-rc3 --check
```

`scripts/verify_claims.py` is the mechanism that keeps the prose honest:
`scripts/recompute_ground_truth.py` recomputes every headline number
from the released tables into `review/ground_truth.tsv`, and the verifier
fails if the manuscript disagrees with it. The report draft has its own
verifier (`scripts/verify_report.py`).

---

## Release v1.0-rc3 — index of deliverables

| File | Owner | What it is |
|---|---|---|
| `docs/pre_registration.md` | A | H1/H2/H3 hypotheses + locked inputs + analysis plan; **Zenodo DOI 10.5281/zenodo.23121762** |
| `docs/pre_registration_manifest.json` | A | SHA-256 fingerprints of pre-reg + key tables for tamper-evidence |
| `docs/manuscript_a_sections.md` | A | Manuscript sections §A-J: network + engine + controls + stats + results + discussion |
| `docs/manuscript_b_sections.md` | B | Manuscript sections §A-H: label panel + transfer + pubcount + ClinGen |
| `docs/manuscript_intro_abstract.md` | B | Abstract, introduction, related work, combined figure legends |
| `docs/research_report_draft.md` | A+B | **The SCI-standard manuscript draft** (title → abstract → methods → results → discussion → data availability → appendix) |
| `docs/control8_seed_degree.md` | A | Control 8 (Zhang et al. 2021): seed-degree stratification |
| `docs/control9_size_adjustment.md` | A | Control 9 (Wright et al. 2025): network size and coverage adjustment |
| `docs/control10_random_input.md` | A | Control 10 (Hill et al. 2019): random-input degeneracy |
| `docs/pre_registration_amendments.md` | A+B | Amendments 1–5: every deviation from the frozen pre-registration |
| `review/ground_truth.tsv` | A | Every headline number, recomputed from the released tables |
| `review/verification_report.txt` | A | The manuscript-vs-tables check (57 PASS / 0 FAIL) |
| `docs/public_release_protocol.md` | B | Release protocol v1.0 — leaderboard, Zenodo DOI, citation, reproducibility contract |
| `docs/reproducibility_audit.md` | B | Independent reproducibility audit (M4.3-M5.5) |
| `results/tables/a4_table1_delta_auroc.tsv` | A | Primary endpoint Table 1: 6 nets × 4 labels × 3 alphas = 72 cells |
| `results/tables/w5_*.tsv` | A | Control 3 (rewiring null): 23/24 cells reject null at z > 5 |
| `results/tables/a5_*.tsv` | A | Control 2 (degree-matched seed null): 23/24 cells positive |
| `results/tables/d09_table.tsv` | A | Control 6 (temporal split) |
| `results/tables/d11_table.tsv` | A | Control 5 (per-network cap) |
| `results/tables/d12_table.tsv` | A | Control 4 (shared universe) |
| `results/tables/results_stats.tsv` | A | Primary endpoint stats, 24 cells: BH-FDR q<0.05 in **23/24** (the one failure is `reactome x intogen_temporal_new`, whose point estimate is negative) |
| `results/tables/stats_fine_z.tsv` | A | Fine-mode per-fold z (A5) |
| `results/tables/stats_variance_decomp.tsv` | A | Mixed-effects variance partition |
| `results/tables/stats_bootstrap_ci.tsv` | A | Cluster bootstrap 95 % CIs |
| `results/tables/stats_equivalence.tsv` | A | TOST equivalence tests (margin=0.05) |
| `results/tables/provenance_ablation.tsv` | A | Provenance ablation (RNA/protein/curated) |
| `results/tables/gnn_courtesy_arm.tsv` | A | §7 control: Node2Vec vs degree |
| `figures/figure1_pilot_overview.png` | B | Pilot overview |
| `figures/figure2_*.png` | A | Main paper figures: heatmap, control panel, temporal strip |
| `docs/w5_log.md` | A | W5 reproducer + interpretation |
| `docs/a5_log.md` | A | A5 reproducer + interpretation |
| `docs/s2_log.md` | A | S2 (4-arm grid + Table 1 + Figure 2 + stats) reproducer |
| `docs/stats_analysis_log.md` | A | §4.6 (fine z + mixed-effects + bootstrap CI + TOST) |
| `docs/gnn_courtesy_log.md` | A | §7 GNN courtesy arm (Node2Vec vs degree, 72/72 cells) |
| `docs/provenance_summary.tsv` | A | Provenance ablation TSV |
| `resource_and_dataset_table.csv` | B | All external resources with citations + access routes |
| `docs/decisions.md` | A+B | Decision log (D-01..D-NUMPY) |
| `docs/sync_notes.md` | A+B | S0..S3 sync notes |
| `docs/b_log.md` | B | B-side work narrative |
| `docs/a0_log.md` ... `a4_log.md` | A | W0-W2 build reproducer per network |
| `docs/rna_coexp_log.md` | A | 6th network reproducer |

**Total: ~127 KB B-side + ~100 KB A-side docs + ~250 KB tables.**

---

## Exploratory arm (2026-11-05) — stability of a network prior

> **Not part of the confirmatory analysis.** This arm changes the learner's
> output, so it never enters the H1 decision rule and its numbers are not mixed
> with the primary table. It was pre-registered in
> `docs/stability_regularization_plan.md` (committed before any result existed).

| File | What it is |
|---|---|
| `docs/stability_regularization_plan.md` | Pre-registration: the stability-regularised scorer, the predictions P1–P7 with their decision rules, the control stack, and Amendment A1 |
| `docs/stability_regularization.md` | Method + results + limitations: what the operator did, why it could not work, and what does control the instability |
| `results/tables/stability_per_gene.tsv` | Per-gene rank SD, score CV, crowding (24 arms, 40 repeats) |
| `results/tables/stability_by_arm.tsv` | Per-arm stability summary + the crowding correlation |
| `results/tables/stability_crowding_attribution.tsv` | Which per-gene property predicts the rank instability (24 arms) |
| `results/tables/stability_regularization_grid.tsv` | The full operator grid (variant × k × λ × control) |
| `results/tables/stability_control_comparison.tsv` | C1 replicate-level random null + robust Z, per arm and λ |
| `results/tables/stability_control_c1_null.tsv` | The 20-replicate random-smoothing null values |
| `results/tables/stability_control_c3_flattening.tsv` | The monotone-flattening invariance check |
| `results/tables/stability_selected_scorer.tsv` | Cross-fitted λ selection and its out-of-sample performance |
| `results/tables/stability_seed_size_scaling.tsv` | Stability and discrimination as a function of the training seed count |
| `results/tables/stability_seed_size_exponents.tsv` | The fitted 1/√n exponent per arm |
| `results/tables/stability_tradeoff_curve.tsv` | The (stability, discrimination) curve of every operator family |
| `results/tables/stability_prediction_scorecard.tsv` | Every pre-registered prediction, its rule, the observed value and the verdict |
| `results/tables/stability_splits_v40_manifest.tsv` | SHA-256 of the 40-repeat split panel (the panel itself is derived data) |
| `results/tables/stability_per_gene_scores_manifest.tsv` | SHA-256 of the 72 per-gene score dumps |
| `figures/figure8_stability_tradeoff.png` | The trade-off figure (smoothing axis, seed axis, mechanism) |
| `d1/engine/stability.py` | The operators and the stability metrics |
| `scripts/build_splits_v40.py` | Builds and *verifies* the 40-repeat split panel |
| `scripts/run_stability_dump.py` | Per-gene score dump + AUROC cross-check against the harness |
| `scripts/run_stability_regularized.py` | The operator grid with the controls, over all 24 arms |
| `scripts/run_stability_controls.py` | The replicate-level random-neighbourhood null |
| `scripts/run_stability_seedsize.py` | The seed-count axis |
| `scripts/run_stability_attribution.py` | Which per-gene property predicts the instability |
| `scripts/summarize_stability.py` | The trade-off curve and the prediction scorecard |

```bash
python scripts/build_splits_v40.py            # 40-repeat panel + byte-identity check
python scripts/run_stability_dump.py --alpha 0.3 0.5 0.7
python scripts/run_stability_regularized.py --alpha 0.5
python scripts/run_stability_controls.py
python scripts/run_stability_seedsize.py
python scripts/run_stability_attribution.py
python scripts/summarize_stability.py
python scripts/make_figure_stability_tradeoff.py
```

**Headline.** The pre-registered stability-regularised scorer does not change the
ranking at all (the median per-gene rank SD moves by ≤ 5.3 × 10⁻⁵ in 24/24 arms),
because shrinking toward the *score* neighbourhood is order-preserving; every
operator that does change the ranking makes the prior *less* stable. The
instability of a network prior is instead the sampling noise of a mean over
seeds — the per-gene rank SD falls as n_seeds^−0.49 — and since RWR is linear in
the seed weights, bagging over seed resamples is provably a no-op for it.

---

## 1. Setup (once per computer)

Use the **Anaconda Prompt** or **Anaconda PowerShell Prompt** (Windows Start menu).

```bash
git clone https://github.com/lanxin2391/d1-benchmark.git
cd d1-benchmark
git checkout v1.0-rc3     # or main
conda env create -f environment.yml
conda activate d1
pip install -e .          # makes `import d1` work everywhere
pytest                    # all tests must pass (currently 131 tests, 3 skipped)
```

After adding a package: `conda env export --no-builds > environment.lock.yml` and commit it
(the this is the locked environment required by protocol §4.7).

**Line endings.** The released tables are LF (`*.tsv text eol=lf` in
`.gitattributes`). `pandas.to_csv` writes `os.linesep`, so on Windows a
regenerated table would otherwise carry CRLF and its SHA-256 would not
match a checkout. Every writer in `scripts/` therefore passes
`lineterminator="\n"`; `scripts/normalise_line_endings.py` repairs any
file that slips through, and the release manifest hashes the
CRLF-normalised bytes.

## 2. Everyday commands

```bash
conda activate d1
git checkout a-networks          # or b-labels (main has the merged v1.0-rc3)
git pull origin main             # get the other person's merged work
pytest                           # before every commit (currently 131 tests)
python scripts/verify_claims.py  # the manuscript still matches the tables
git add -A
git commit -m "A1: evaluation harness"
git push origin a-networks       # then open a Pull Request on GitHub
```

Useful scripts:

```bash
python scripts/run_toy_benchmark.py      # end-to-end check on synthetic networks
python scripts/bench_rwr_speed.py        # RWR speed on a FunMap-sized graph
python scripts/run_arm.py --help         # run one network x label arm (real data)
python scripts/run_a5.py --n-jobs 4     # Control 2 (degree-matched seed null)
python scripts/run_rewiring.py           # Control 3 (edge permutation null)
python scripts/run_stats_analysis.py     # §4.6 (mixed-effects + bootstrap + TOST)
python scripts/run_gnn_courtesy.py --n-jobs 4  # §7 GNN courtesy baseline
```

## 3. Folder layout

```
data/raw/             downloads, never edited, not in git (networks/ labels/ cohorts/ reference/)
data/processed/       generated by scripts: networks/ labels/ splits/ reference/
data/manifest.tsv     file, bytes, md5, download date (B)
d1/                   python package
scripts/              one script per output file
tests/                pytest unit tests (70 tests passing)
results/              run tables, figures
docs/                 release deliverables (see index above)
docs/decisions.md     every choice not fixed by the protocol
docs/sync_notes.md    five-line note after every sync point
```

---

## 4. Interface contracts (do not change without both people agreeing)

**C1 Gene ID.** HGNC approved symbol (e.g. `TP53`), always produced by `d1.hgnc.HGNCMapper`.
IDs that map to more than one approved symbol are dropped and logged.

**C2 Network file.** `data/processed/networks/{net_id}.tsv`, tab-separated, header
`gene_a<TAB>gene_b<TAB>weight`. Undirected, one row per edge, `gene_a < gene_b` (alphabetical),
no self-loops, no duplicates, largest connected component (LCC) only. Written only by
`d1.networks.io.finalize_edges`. `weight` is kept for reference; the primary analysis is
unweighted (D-01).

**C2b Network meta.** `{net_id}.meta.json` with at least: `net_id`, `source_file`, `threshold`,
`n_input_edges`, `n_mapped_edges`, `n_edges_dedup`, `n_components`, `n_nodes_lcc`,
`n_edges_lcc`, `date`, `script`.

**C3 Label file.** `data/processed/labels/{label_id}.tsv`, header `gene rank source release`.
Positives only; `rank` 1 = strongest evidence (used for common-size capping). Negatives are all
other genes of the evaluation universe.

**C4 Split file.** `data/processed/splits/{label_id}__{universe_id}.tsv`, header
`gene repeat fold y`. `repeat` 0–9, `fold` 0–4, `y` 1/0. Each gene exactly once per repeat, same
gene set in every repeat, stratified by `y`. Checked by `d1.engine.evaluate.check_splits`.

**C5 Result table.** `results/runs/{run_id}.tsv`, one row per
network × label_set × universe × alpha × repeat × fold × method, columns
`network label_set universe alpha repeat fold method n_iter auroc auprc p_at_50 p_at_100 p_at_500 n_test_pos n_test`.
Degree rows have `alpha` = NaN (degree does not depend on alpha).

**C6 Mapper API.** `from d1.hgnc import HGNCMapper; m = HGNCMapper(path); m.map(series, id_type)`
returns a Series of symbols (NaN if unmapped/ambiguous); `id_type` ∈ `symbol, entrez, ensembl, uniprot`.

### Fixed IDs

- Networks: `funmap`, `rna_coexp`, `string_phys700`, `intact`, `string_full700`
  (+ `string_full400`, `string_full900` sensitivity), `reactome`
- Labels: `intogen2024`, `ot_all`, `ot_nolit`, `intogen_temporal_new`
- Universes: `native_{net_id}` (that network's LCC) and `shared` (genes in all six networks)

---

## 5. Reproducing the release end-to-end

The v1.0-rc3 release was regenerated on a clean machine via the
B-side reproducibility audit (`docs/reproducibility_audit.md`) and the
line-ending normalisation of the released tables
(`scripts/normalise_line_endings.py`).
The audit script (`scripts/audit_reproducibility.py`) reproduces all
24 primary cells (× 3 alphas) in ~30 minutes on a 16-core machine.

---

## 6. Citing

If you use this benchmark, please cite:

* The protocol: *D1 benchmark: What a network prior contributes to
  cancer-driver prioritisation.*
* The pre-registration: DOI 10.5281/zenodo.23121762.
* The v1.0-rc3 release: this GitHub repository at the
  `v1.0-rc3` tag.
* The manuscript draft: `docs/research_report_draft.md`.