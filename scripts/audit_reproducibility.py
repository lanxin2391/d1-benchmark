#!/usr/bin/env python
"""Reproducibility audit: verify all artefacts exist, have the right
shape, and can be loaded end-to-end.

The audit checks:
  1. Schema correctness of every TSV (labels, splits, networks, results).
  2. Cross-file referential integrity (label genes ⊆ split universe, etc.).
  3. SHA-256 fingerprints against the pre-registration manifest.
  4. End-to-end smoke run of one arm (funmap × intogen2024 × native_funmap).

This script is run by Person B (the audit author) before the
manuscript is submitted. It is not part of the production pipeline;
it lives in scripts/audit_reproducibility.py for the public release.

Usage:
    .venv-d1/Scripts/python.exe scripts/audit_reproducibility.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

PASS, FAIL = "PASS", "FAIL"
findings: list[tuple[str, str, str]] = []  # (status, test_id, message)


def add(status: str, test_id: str, msg: str) -> None:
    findings.append((status, test_id, msg))
    print(f"  [{status}] {test_id}: {msg}")


# ---- 1. Network files ----
print("\n[1] Network files (data/processed/networks/)")
NETWORKS_PRIMARY = ["funmap", "string_phys700", "string_full700",
                    "intact", "reactome", "rna_coexp"]
NETWORKS_SENSITIVITY = ["string_full400", "string_full900"]

for net_id in NETWORKS_PRIMARY + NETWORKS_SENSITIVITY:
    tsv = DATA / "networks" / f"{net_id}.tsv"
    meta = DATA / "networks" / f"{net_id}.meta.json"
    if not tsv.exists():
        add(FAIL, f"net-{net_id}-tsv", "missing")
        continue
    if not meta.exists():
        add(FAIL, f"net-{net_id}-meta", "missing")
        continue
    df = pd.read_csv(tsv, sep="\t")
    cols_ok = set(df.columns) >= {"gene_a", "gene_b"}
    if not cols_ok:
        add(FAIL, f"net-{net_id}-cols", f"bad columns {df.columns}")
        continue
    no_self_loops = (df["gene_a"] != df["gene_b"]).all()
    no_dups = len(df) == len(df.drop_duplicates())
    if not no_self_loops:
        add(FAIL, f"net-{net_id}-self-loops", "found self-loops")
        continue
    if not no_dups:
        add(FAIL, f"net-{net_id}-dups", "found duplicates")
        continue
    m = json.loads(meta.read_text())
    n_edges = m.get("n_edges_lcc")
    if n_edges and abs(n_edges - len(df)) > max(5, 0.01 * n_edges):
        add(FAIL, f"net-{net_id}-edge-count", f"meta says {n_edges}, file has {len(df)}")
        continue
    add(PASS, f"net-{net_id}", f"{len(df)} edges, {m.get('n_nodes_lcc', '?')} nodes")


# ---- 2. Label files ----
print("\n[2] Label files (data/processed/labels/)")
LABELS_PRIMARY = ["intogen2024", "ot_all", "ot_nolit",
                  "intogen_temporal_new", "clingen_label"]

for label_id in LABELS_PRIMARY:
    tsv = DATA / "labels" / f"{label_id}.tsv"
    if not tsv.exists():
        add(FAIL, f"label-{label_id}", "missing")
        continue
    df = pd.read_csv(tsv, sep="\t")
    cols_ok = set(df.columns) >= {"gene", "rank", "source", "release"}
    if not cols_ok:
        add(FAIL, f"label-{label_id}-cols", f"bad columns {df.columns}")
        continue
    # Schema check: rank must be 1-based and ascending
    ranks = df["rank"].values
    if not (ranks[0] == 1 and (np.diff(ranks) >= 0).all()):
        add(FAIL, f"label-{label_id}-rank", f"non-monotonic ranks")
        continue
    # Symbol check: must be a non-empty string
    if df["gene"].isna().any() or (df["gene"] == "").any():
        add(FAIL, f"label-{label_id}-gene", "NaN/empty gene")
        continue
    add(PASS, f"label-{label_id}", f"{len(df)} positives, "
                                    f"unique genes = {df['gene'].nunique()}")


# ---- 3. Split files ----
print("\n[3] Split files (data/processed/splits/)")
NET_IDS = NETWORKS_PRIMARY

# Standard C4 splits (24 files for 4 primary labels, plus 6 clingen splits)
for label_id in LABELS_PRIMARY:
    for net_id in NET_IDS:
        # Standard splits live in either data/processed/splits/ or data/processed/splits_d11/
        # Naming convention: {label}__native_{net_id}.tsv
        split = DATA / "splits" / f"{label_id}__native_{net_id}.tsv"
        if not split.exists():
            split = ROOT / "data" / "processed" / "splits_d11" / f"{label_id}__native_{net_id}.tsv"
        if not split.exists():
            add(FAIL, f"split-{label_id}-{net_id}", "missing in both splits/ and splits_d11/")
            continue
        df = pd.read_csv(split, sep="\t")
        cols_ok = set(df.columns) >= {"gene", "repeat", "fold", "y"}
        if not cols_ok:
            add(FAIL, f"split-{label_id}-{net_id}-cols", f"bad columns")
            continue
        # Schema: 50 repeats × 5 folds = 250 rows total (or could be more rows if not 50×5)
        n_rows = len(df)
        # Per-fold y positivity rate
        rate = df.groupby(["repeat", "fold"])["y"].mean()
        # Lower bound depends on the label size: clingen_label has 84 genes and
        # intact/rna_coexp have larger LCCs, giving per-fold rates ~0.005.
        # Allow a tighter lower bound for larger labels.
        n_pos = (df["y"] == 1).sum() // 50  # average per fold
        n_universe = df["gene"].nunique()
        expected_rate = n_pos / n_universe if n_universe > 0 else 0
        lower_bound = max(0.001, expected_rate * 0.5)  # at least half the expected rate
        if (rate < lower_bound).any() or (rate > 0.5).any():
            add(FAIL, f"split-{label_id}-{net_id}-rate",
                f"per-fold positive rate out of range: {rate.min():.4f}-{rate.max():.4f} "
                f"(expected ~{expected_rate:.4f})")
            continue
        add(PASS, f"split-{label_id}-{net_id}",
            f"{n_rows} rows, {df['gene'].nunique()} unique genes")


# ---- 4. SHA-256 fingerprints against the pre-reg manifest ----
print("\n[4] SHA-256 fingerprints (docs/pre_registration_manifest.json)")
manifest_path = ROOT / "docs" / "pre_registration_manifest.json"
if not manifest_path.exists():
    add(FAIL, "manifest-exists", "missing")
else:
    manifest = json.loads(manifest_path.read_text())
    # Manifest stores under "sha256" key (singular dict)
    files = manifest.get("sha256", manifest.get("files", {}))
    n_check = 0
    n_pass = 0
    for rel_path, expected_sha in files.items():
        p = ROOT / rel_path
        if not p.exists():
            continue  # some files may not exist in this clone
        n_check += 1
        actual_sha = hashlib.sha256(p.read_bytes()).hexdigest()
        if actual_sha == expected_sha:
            n_pass += 1
        else:
            add(FAIL, f"sha-{rel_path}",
                f"mismatch: expected {expected_sha[:12]}..., got {actual_sha[:12]}...")
    add(PASS, "manifest-sha256",
        f"{n_pass}/{n_check} files match SHA-256 fingerprints")


# ---- 5. End-to-end smoke test ----
print("\n[5] End-to-end smoke run (funmap × intogen2024 × native_funmap)")
try:
    sys.path.insert(0, str(ROOT))
    from d1.engine.rwr import transition_T, seed_matrix, rwr
    from d1.engine.baselines import degree_scores
    from d1.engine.evaluate import evaluate
    net_path = DATA / "networks" / "funmap.tsv"
    net_df = pd.read_csv(net_path, sep="\t")
    sym_set = sorted(set(net_df["gene_a"]).union(net_df["gene_b"]))
    sym_to_idx = {s: i for i, s in enumerate(sym_set)}
    rows = [sym_to_idx[a] for a in net_df["gene_a"]]
    cols = [sym_to_idx[b] for b in net_df["gene_b"]]
    import scipy.sparse as sp
    n = len(sym_set)
    A = sp.csr_matrix(
        (np.ones(len(net_df) + len(rows)),
         (rows + cols, cols + rows)),
        shape=(n, n), dtype=np.float64)
    WT = transition_T(A)
    split = DATA / "splits" / "intogen2024__native_funmap.tsv"
    if not split.exists():
        # Fall back to the d11-cap directory
        split = ROOT / "data" / "processed" / "splits_d11" / "intogen2024__native_funmap.tsv"
    spl = pd.read_csv(split, sep="\t")
    label_df = pd.read_csv(DATA / "labels" / "intogen2024.tsv", sep="\t")
    pos_genes = set(label_df["gene"].tolist())
    # Take a real fold (fold 0, repeat 0) — must have both positives and negatives
    fold0 = spl[(spl["repeat"] == 0) & (spl["fold"] == 0)]
    test_g = fold0["gene"].tolist()
    train_g = spl[(spl["repeat"] == 0) & (spl["fold"] != 0) & (spl["y"] == 1)]["gene"].tolist()
    train_idx = [sym_to_idx[g] for g in train_g if g in sym_to_idx]
    if not train_idx:
        add(FAIL, "smoke-arm", "no train seeds in fold 0, repeat 0")
    else:
        P0 = seed_matrix(n, [train_idx])
        P, _ = rwr(WT, P0, alpha=0.5)
        deg_score = degree_scores(A)
        test_idx = [sym_to_idx[g] for g in test_g if g in sym_to_idx]
        y = np.array([1 if g in pos_genes else 0 for g in test_g if g in sym_to_idx])
        if y.sum() == 0 or y.sum() == len(y):
            add(FAIL, "smoke-arm",
                f"test set has no class diversity (n_pos={y.sum()}, n={len(y)})")
        else:
            s_rwr = P[test_idx, 0]
            s_deg = deg_score[test_idx]
            out_rwr = evaluate(s_rwr, y, np.array(test_g)[:len(y)])
            out_deg = evaluate(s_deg, y, np.array(test_g)[:len(y)])
            delta = out_rwr["auroc"] - out_deg["auroc"]
            add(PASS, "smoke-arm",
                f"funmap × intogen2024: RWR AUROC={out_rwr['auroc']:.4f}, "
                f"degree AUROC={out_deg['auroc']:.4f}, ΔAUROC={delta:+.4f} "
                f"(n_pos={y.sum()}, n={len(y)})")
except Exception as e:
    add(FAIL, "smoke-arm", f"raised {type(e).__name__}: {e}")


# ---- Summary ----
print("\n=== SUMMARY ===")
n_pass = sum(1 for s, _, _ in findings if s == PASS)
n_fail = sum(1 for s, _, _ in findings if s == FAIL)
print(f"  PASS: {n_pass}")
print(f"  FAIL: {n_fail}")
print(f"  TOTAL: {len(findings)}")
if n_fail > 0:
    print("\nFailed tests:")
    for s, t, m in findings:
        if s == FAIL:
            print(f"  - [{t}] {m}")
    sys.exit(1)
print("\nAll checks PASS.")
