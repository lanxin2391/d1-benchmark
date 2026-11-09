# D1 benchmark — v1.0-rc1 release note

**Release commit:** `093efe1b604e7f97b3f6a57504575dcc1fe0eabe` (tag `v1.0-rc1`)
**Built:** 2026-11-05
**Zenodo DOI:** *pending* — the three zips have to be uploaded from a
browser, and Zenodo mints the DOI at that point. The pre-registration
(DOI 10.5281/zenodo.23121762) is the frozen v1 and is unchanged.

## What changed since v1.0-prep

| area | change |
|---|---|
| P0 repairs | all 11 done: the alpha filter, the TOST rule, the D-11 filename collision (plus a second bug: the cap was applied before the universe restriction), the D-12 shared universe, the D-09 seed resampling, the variance decomposition, the gene-cluster bootstrap, the 24-cell BH-FDR, the H1 rule as re-computed, and the control-2 null moved to the 40-repeat design |
| Phase D | the complex-aware cross-validation: the complex map (ComplexPortal via UniProtKB), the leakage audit of the v1 splits (96 % of the complex-bearing held-out positives have a partner in the training fold), the group-preserving splitter, the v2 splits and the v1-vs-v2 comparison (the effect shrinks by a median of 12 %, the H1 verdict is unchanged) |
| Phase S | the stability atlas: per-gene stability, its attribution to the crowding of the score distribution, the network level, and the replication on ClinGen |
| direction 1 | the consensus prior over the six networks, and its robustness checks |
| E3 | the hub-capped RWR: the headline is not a hub effect, but two thirds of the raw margin sits in the hubs' ranking position |
| P1 | the rewiring metadata now records the swaps that were executed, and the mixing was measured |

## What is in the three zips

| zip | size | content |
|---|---|---|
| `d1-benchmark-v1.0-source.zip` | 204.2 MB | the full repository at the release tag, minus `data/raw/`, `data/cache/`, `data/processed/networks/` and `data/processed/rewired/` (all regenerable from the build scripts) |
| `d1-benchmark-v1.0-splits.zip` | 142.1 MB | every split family: the v1 10-repeat set, the 40-repeat design of record, the complex-aware v2 and the ClinGen 40-repeat splits |
| `d1-benchmark-v1.0-results.zip` | 43.7 MB | `results/tables/` and `figures/` |

## Verifying the release

`release_manifest.json (inside the source zip) and docs/release_manifest_v1.0-rc1.json (in the repository)` (inside the source zip) holds the
SHA-256 of the 23 key tables, the 5 figures and the 3 zips:

| zip | SHA-256 |
|---|---|
| `d1-benchmark-v1.0-results.zip` | `55914fa4186bd69c...` |
| `d1-benchmark-v1.0-source.zip` | `235da8dddd246437...` |
| `d1-benchmark-v1.0-splits.zip` | `ec8e9616cc02f016...` |

```bash
python -c "import hashlib,sys;print(hashlib.sha256(open(sys.argv[1],'rb').read()).hexdigest())" d1-benchmark-v1.0-source.zip
```

## Reproducing the headline numbers

```bash
pytest                                        # 113 tests
python scripts/run_a4.py --n-rep 40 --splits-dir data/processed/splits_v40 \
       --runs-dir results/runs_v40 --table-suffix _v40
python scripts/run_statistics.py --combined results/runs_v40/a4_combined.tsv \
       --tables results/tables --suffix _v40
python scripts/audit_reproducibility.py        # 39/44 (the 5 are the
                                               # expected fingerprint
                                               # mismatches, section 6)
```
