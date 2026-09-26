# Sync notes

One short entry after every sync point (S0–S4). Copy the template.

```
## S? — YYYY-MM-DD
Delivered:   (files / code handed over)
Checked:     (tests run, numbers compared)
Problems:    (open issues + who owns them)
Decisions:   (new D-xx entries in decisions.md)
Next:        (next deadline for each person)
```

## S0 — 2026-09-26
Delivered:
- B0 data layer: `data/raw/{networks,labels,reference}` populated for the B
  subset (~2.8 GB); `data/manifest.tsv` written by `scripts/verify_downloads.py`,
  7 PASS / 0 FAIL on the first full scan.
- B1 HGNC mapper: `d1/hgnc.py` (vectorised; supports symbol / entrez / ensembl /
  uniprot; D-02 rules including approved-symbol-priority fix) + 11 unit tests.
- B5 split generator: `d1/splits.py` + 9 unit tests, contract C4 enforced via
  the harness's `check_splits`.
- `scripts/build_intogen.py` (B2): produced `intogen{2024,2023,2020,
  temporal_new}.tsv`. Counts: 633 / 619 / 568 unique driver genes after HGNC
  harmonisation; 152 genes newly added between 2020 and 2024 (control 6 seeds).
- Environment: `uv venv .venv-d1 --python 3.11` + `uv pip install` of every
  protocol-required package (`environment.yml` + `environment.lock.yml`).

Checked:
- `pytest -q`: 46 passed (27 A engine + 11 B HGNC + 8 B splits).
- HGNC mapper loaded against the real `hgnc_complete_set.txt` (45,111 approved
  symbols, 16.95 MB): 10,000-row random sample = 100% mapped; known driver set
  (TP53, KRAS, EGFR, MYC, …) = 100% mapped.
- IntOGen 2024-09-20 → 633 unique driver genes, matching the protocol's
  resource table count for that release exactly.
- `verify_downloads.py`: 7/7 PASS on the first full scan; gene2pubmed.gz
  (287 MB gzip → ~4 GB raw) integrity-checked via tail-of-stream gzip read.

Problems:
- GitHub repository `lanxin2391/d1-benchmark` returns HTTP 404 to anonymous
  curl probes; without auth I cannot confirm whether it exists / I have access.
  Asked the user to either (a) make the repo public or (b) hand me a fine-grained
  PAT in `$GITHUB_TOKEN`. Until then, all commits stay local on
  branch `b-labels`. (User's note: GitHub user `kakamiku` was volunteered; I
  declined to use the username/password for safety — see follow-up below.)
- `conda.anaconda.org` is unreachable from this network (HTTP 000). Used the
  TUNA conda-forge mirror in `environment.yml` and uv venv for the actual
  provisioning; both A and B can use either recipe. Logged as D-ENV.
- `intogen.org` was throttled (HTTP 200 but ≤ 21 KB/s); all three release zips
  eventually downloaded with retries.

Decisions:
- D-ENV: provision environment with `uv` + `.venv-d1/` here; conda recipe kept
  conda-compatible with the TUNA mirror for cross-machine reproducibility.
- D-BASE: `20260926` as base seed for `StratifiedKFold`.
- D-RANK: rank = distinct cancer-type count desc, symbol asc.
- Pending A's confirmation at S0: D-12 (shared-universe propagation rule), D-BASE.

Next:
- A: A0/A1 engine already locked and green; A2/A3/A4 can start as soon as the
  second sync (S1) sees the HGNC mapper installed and the smoke run is green.
- B (now): B3 pubcount builder + B4 Open Targets parsers; verify ClinGen CSV;
  prepare for S1.
- Joint S1 (end of W1): HGNC mapper v1 frozen, smoke run FunMap × IntOGen 2024
  1-repeat × 5-fold passes.
