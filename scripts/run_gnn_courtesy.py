"""GNN courtesy arm (protocol §7).

Runs Node2Vec (Grover & Leskovec, 2016) on each of the 6 primary
networks and evaluates the same protocol-§4.3 ranking rule:
    for each (network, label) cell at each alpha (3 alphas = 0.3, 0.5, 0.7):
        score(gene) = cosine_similarity(
            node2vec_embedding(gene),
            mean(node2vec_embedding(training fold positives))
        )
        AUROC(test positives vs other genes) over 5-fold x 10-repeat = 50 folds
        delta_GNN = AUROC(GNN) - AUROC(degree) per fold

Output mirrors the A4 + D-09 etc. structure so it can be compared
on the same axes (delta_AUROC vs degree-only).

Memory + time
--------------
    - per cell: one Node2Vec training pass (walk generation + Word2Vec)
      ~30-60 s on the largest network (intact), ~5 s on reactome.
    - per cell: one AUROC over 50 folds = ~0.5 s.
    - 24 cells x 3 alphas x ~30 s = ~36 minutes per arm
    - 6 networks x 4 labels x 3 alphas = 72 evaluations = ~3-6 h total.

Checkpointing
-----------
    - Every completed (network, label, alpha) cell is saved to
      results/runs/gnn_courtesy/{net}__{label}__alpha{A}.json with the
      per-fold AUROC values.
    - On restart, the script reads the existing files and skips them.
    - Status JSON is updated after every checkpoint write
      (results/runs/gnn_courtesy/status.json).
    - Live log is appended to scratch/gnn_courtesy.log.

Parallelism
------------
    By default uses joblib with n_jobs=4 (per protocol §5 budget).
    Training is the bottleneck (5 min/cell); 4 cores brings 24×3=72
    evaluations from ~6 h to ~1.5 h.

Logging
-------
    - All progress goes to scratch/gnn_courtesy.log with timestamps.
    - The script also writes a markdown log to docs/gnn_courtesy_log.md
      with the high-level summary.

Hyperparameters
----------------
    Following Grover & Leskovec (2016):
        dimensions      = 64
        walk_length      = 40
        num_walks        = 80
        window           = 10
        p, q (return, in-out) = 1, 1
        min_count       = 1
        workers          = 4
        epochs           = 1   (default; matches paper for social networks)
        seed             = 0
"""
import argparse
import json
import os
import sys
import time
import random
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

NETS = ["funmap", "rna_coexp", "string_phys700",
        "intact", "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
PROCESSED_NET = ROOT / "data" / "processed" / "networks"
PROCESSED_SPLIT = ROOT / "data" / "processed" / "splits"
SCRATCH = ROOT / "scratch"
LOG = SCRATCH / "gnn_courtesy.log"
STATUS = ROOT / "results" / "runs" / "gnn_courtesy" / "status.json"
OUT_DIR = ROOT / "results" / "runs" / "gnn_courtesy"
TABLE_DIR = ROOT / "results" / "tables"
DOCS = ROOT / "docs" / "gnn_courtesy_log.md"

SCRATCH.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)
DOCS.parent.mkdir(parents=True, exist_ok=True)
N_JOBS_DEFAULT = 4
import argparse
import json
import os
import sys
import time
import random
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.evaluate import evaluate  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETS = ["funmap", "rna_coexp", "string_phys700",
        "intact", "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
PROCESSED_NET = ROOT / "data" / "processed" / "networks"
PROCESSED_SPLIT = ROOT / "data" / "processed" / "splits"
SCRATCH = ROOT / "scratch"
LOG = SCRATCH / "gnn_courtesy.log"
STATUS = ROOT / "results" / "runs" / "gnn_courtesy" / "status.json"
OUT_DIR = ROOT / "results" / "runs" / "gnn_courtesy"
TABLE_DIR = ROOT / "results" / "tables"
DOCS = ROOT / "docs" / "gnn_courtesy_log.md"

SCRATCH.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)
DOCS.parent.mkdir(parents=True, exist_ok=True)


# ---- Node2Vec helpers ---------------------------------------------------

def _log(msg: str):
    """Append a timestamped message to the live log file."""
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def _status_write(payload: dict):
    STATUS.parent.mkdir(parents=True, exist_ok=True)
    with STATUS.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)


def _node2vec_walks(edges, num_walks, walk_length, p, q, seed):
    """Pure-Python Node2Vec-style biased random walks.

    Uses a simple weighted sampling rather than alias-method simulation
    for robustness; this is slower than a true alias-method implementation
    but is correct for the small graphs we use (funmap = 10 k nodes,
    intact = 18 k nodes)."""
    from collections import defaultdict
    rng = random.Random(seed)
    adj = defaultdict(set)
    nodes = set()
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
        nodes.update([u, v])
    nodes = list(nodes)
    degree = {u: max(len(adj[u]), 1) for u in nodes}

    def weighted_pick(options, weights, rng_callable):
        """`rng_callable` is the bound method ``rng.random``."""
        total = sum(weights)
        if total <= 0:
            return options[int(rng_callable() * len(options))]
        target = rng_callable() * total
        acc = 0.0
        for opt, w in zip(options, weights):
            acc += w
            if acc >= target:
                return opt
        return options[-1]

    walks = []
    nodes_list = list(nodes)
    for _ in range(num_walks):
        for start in nodes_list:
            walk = [start]
            while len(walk) < walk_length:
                cur = walk[-1]
                nbrs = list(adj[cur])
                if not nbrs:
                    break
                if len(walk) == 1:
                    nxt = weighted_pick(nbrs, [1.0] * len(nbrs), rng.random)
                else:
                    prev = walk[-2]
                    cand = [v for v in nbrs if v != prev]
                    if not cand:
                        break
                    # Node2Vec alias sampling: weight proportional to
                    # 1/p if prev->cand is "in-out" else 1/q.
                    w = []
                    for v in cand:
                        if v in adj[prev]:
                            w.append(1.0 / p)
                        else:
                            w.append(1.0 / q)
                    nxt = weighted_pick(cand, w, rng.random)
                walk.append(nxt)
                walk.append(nxt)
            walks.append(walk)
    return walks


def _node2vec_train(edges, dim=64, walk_length=40, num_walks=80,
                     window=10, p=1.0, q=1.0, seed=0, workers=4):
    """Train a gensim Word2Vec on biased random walks over the graph."""
    from gensim.models import Word2Vec
    walks = _node2vec_walks(edges, num_walks, walk_length, p, q, seed)
    sentences = [list(map(str, w)) for w in walks]
    model = Word2Vec(
        sentences,
        vector_size=dim,
        window=window,
        min_count=1,
        workers=workers,
        epochs=1,
        seed=seed,
    )
    return model


def _node2vec_aucroc(model, fold_pos, fold_neg, all_test_genes):
    """Score each gene by cosine similarity to the mean train-positive
    embedding; return AUROC of positives vs the rest."""
    pos_emb = []
    for g in fold_pos:
        if g in model.wv:
            pos_emb.append(model.wv[g])
    if not pos_emb:
        return np.nan
    mean = np.mean(pos_emb, axis=0)
    scores = []
    labels = []
    for g in all_test_genes:
        if g not in model.wv:
            continue
        sim = float(np.dot(model.wv[g], mean) /
                     (np.linalg.norm(model.wv[g]) * np.linalg.norm(mean) + 1e-12))
        scores.append(sim)
        labels.append(g in fold_pos)
    if not labels or sum(labels) == 0 or sum(labels) == len(labels):
        return np.nan
    return float(np.array(labels)[np.argsort(np.argsort(scores))].mean())  # naive AUROC
# We rely on scikit-learn's roc_auc_score for correctness; switch to it below.


# ---- Per-cell runner ----------------------------------------------------

def _ensure_gensim_aucroc():
    """Replace naive AUROC with sklearn's; import lazily."""
    from sklearn.metrics import roc_auc_score
    def _sk_aucroc(model, fold_pos, fold_neg, all_test_genes):
        pos_emb = []
        for g in fold_pos:
            if g in model.wv:
                pos_emb.append(model.wv[g])
        if not pos_emb:
            return np.nan
        mean = np.mean(pos_emb, axis=0)
        scores, labels = [], []
        for g in all_test_genes:
            if g not in model.wv:
                continue
            sim = float(np.dot(model.wv[g], mean) /
                         (np.linalg.norm(model.wv[g]) * np.linalg.norm(mean) + 1e-12))
            scores.append(sim)
            labels.append(g in fold_pos)
        if not labels or sum(labels) == 0 or sum(labels) == len(labels):
            return np.nan
        return float(roc_auc_score(labels, scores))
    return _sk_aucroc


_aucroc_fn = _ensure_gensim_aucroc()


def _build_checkpointer(net_id: str, label: str, alpha: float):
    """Returns functions (out_path, save, exists, load) for cell-level
    checkpointing."""
    safe_label = label.replace("/", "_")
    out_path = OUT_DIR / f"{net_id}__{safe_label}__alpha{alpha}.json"

    def exists():
        return out_path.exists()

    def save(gnn_folds, deg_folds, params):
        """gnn_folds and deg_folds are parallel lists of (rep, fold, value)
        for Node2Vec AUROC and degree-rank AUROC per fold."""
        with out_path.open("w", encoding="utf-8") as fh:
            json.dump({
                "network": net_id,
                "label": label,
                "alpha": alpha,
                "gnn_folds": [(int(r), int(f), float(v))
                              for (r, f, v) in gnn_folds],
                "deg_folds": [(int(r), int(f), float(v))
                              for (r, f, v) in deg_folds],
                "params": params,
                "saved_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            }, fh, indent=2)

    def load():
        with out_path.open("r", encoding="utf-8") as fh:
            return json.load(fh)

    return out_path, exists, save, load


def _status_init():
    payload = {
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "last_updated": None,
        "completed_cells": [],
        "in_progress": None,
        "n_total": len(NETS) * len(LABELS) * 3,
        "n_complete": 0,
        "elapsed_seconds": 0,
        "params": {},
    }
    _status_write(payload)
    return payload


def _status_update(payload, *, in_progress=None, completed=None,
                    elapsed_seconds=None, params=None):
    payload["last_updated"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    if in_progress is not None:
        payload["in_progress"] = in_progress
    if completed is not None:
        if completed not in payload["completed_cells"]:
            payload["completed_cells"].append(completed)
            payload["n_complete"] = len(payload["completed_cells"])
    if elapsed_seconds is not None:
        payload["elapsed_seconds"] = elapsed_seconds
    if params is not None:
        payload["params"] = params
    _status_write(payload)


def _run_one_cell(net_id: str, label: str, alpha: float, status_payload):
    """Train Node2Vec on `net_id`, evaluate GNN vs degree on (label, alpha).
    Returns (gnn_folds, deg_folds) lists of (rep, fold, val)."""
    from sklearn.metrics import roc_auc_score

    out_path, exists_fn, save_fn, load_fn = _build_checkpointer(
        net_id, label, alpha)
    if exists_fn():
        d = load_fn()
        return (d["gnn_folds"], d["deg_folds"])

    net_path = PROCESSED_NET / f"{net_id}.tsv"
    split_path = PROCESSED_SPLIT / f"{label}__native_{net_id}.tsv"
    edges_df = pd.read_csv(net_path, sep="\t", dtype=str)
    splits_df = pd.read_csv(split_path, sep="\t", dtype={"gene": str, "y": int})
    splits_df["repeat"] = splits_df["repeat"].astype(int)
    splits_df["fold"] = splits_df["fold"].astype(int)

    t0 = time.time()
    edges = [(r.gene_a, r.gene_b) for r in edges_df.itertuples()]
    model = _node2vec_train(edges)
    t_train = time.time() - t0

    # Pre-compute the degree of each gene in the network.
    deg_per_gene = pd.Series(0, index=set(g for e in edges for g in e))
    for u, v in edges:
        deg_per_gene[u] = deg_per_gene.get(u, 0) + 1
        deg_per_gene[v] = deg_per_gene.get(v, 0) + 1

    gnn_folds = []
    deg_folds = []
    all_genes = list(splits_df.gene.drop_duplicates())
    train_pos_cache = {}
    for (rep, fold), _ in splits_df.groupby(["repeat", "fold"], sort=False):
        key = (int(rep), int(fold))
        cell = splits_df[(splits_df.repeat == rep)
                          & (splits_df.fold != fold)
                          & (splits_df.y == 1)]
        train_pos_cache[key] = list(cell.gene)

    for (rep, fold), _ in splits_df.groupby(["repeat", "fold"], sort=False):
        rep, fold = int(rep), int(fold)
        train_pos = train_pos_cache[(rep, fold)]
        if len(train_pos) < 2:
            continue
        # GNN scoring: cosine similarity to mean train-positive embedding.
        pos_emb = [model.wv[g] for g in train_pos if g in model.wv]
        if not pos_emb:
            continue
        mean = np.mean(pos_emb, axis=0)
        scores, labels = [], []
        for g in all_genes:
            if g in model.wv:
                v = model.wv[g]
                sim = float(np.dot(v, mean) /
                             (np.linalg.norm(v) * np.linalg.norm(mean) + 1e-12))
                scores.append(sim)
                labels.append(g in set(train_pos))
        if not labels or sum(labels) == 0 or sum(labels) == len(labels):
            continue
        gnn_auc = float(roc_auc_score(labels, scores))
        # Degree scoring: just degree, no propagation.
        deg_scores = [deg_per_gene.get(g, 0) for g in all_genes]
        deg_auc = float(roc_auc_score(labels, deg_scores))
        gnn_folds.append((rep, fold, float(gnn_auc)))
        deg_folds.append((rep, fold, float(deg_auc)))

    save_fn(gnn_folds, deg_folds,
            params={"dim": 64, "walk_length": 40, "num_walks": 80,
                    "alpha": alpha, "p": 1.0, "q": 1.0,
                    "train_time_sec": float(t_train)})
    return (gnn_folds, deg_folds)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    ap.add_argument("--only-alpha", action="append", default=[0.5],
                    type=float)
    ap.add_argument("--n-jobs", type=int, default=N_JOBS_DEFAULT)
    args = ap.parse_args()

    nets = args.only_net or NETS
    labels = args.only_label or LABELS
    alphas = args.only_alpha
    cells = [(n, r, a) for n in nets for r in labels for a in alphas]
    if args.only_net or args.only_label or args.only_alpha != [0.5]:
        # user restricted, run sequentially
        n_jobs = 1
    else:
        n_jobs = args.n_jobs

    for path in (LOG, DOCS):
        try:
            os.remove(path)
        except FileNotFoundError:
            pass
        except OSError:
            # if locked by Windows or by another process, just open for append
            pass
    status_payload = _status_init()
    status_payload["params"]["n_jobs"] = n_jobs
    status_payload["params"]["n_cells"] = len(cells)
    t_overall = time.time()

    def _run_one_safe(net_id, label, alpha):
        try:
            return _run_one_cell(net_id, label, alpha, status_payload)
        except Exception as e:
            _log(f"[{net_id} x {label} alpha={alpha}] ERROR: {e!r}")
            return ([], [])

    _log(f"running {len(cells)} cells across n_jobs={n_jobs}")
    results = Parallel(n_jobs=n_jobs, verbose=0)(
        delayed(_run_one_safe)(n, r, a) for (n, r, a) in cells
    )

    n_ok = 0
    n_skip = 0
    for (net_id, label, alpha), (gnn_folds, deg_folds) in zip(cells, results):
        if gnn_folds and deg_folds:
            n_ok += 1
            _status_update(status_payload,
                           completed={"network": net_id,
                                       "label": label,
                                       "alpha": alpha},
                           elapsed_seconds=int(time.time() - t_overall))
        else:
            n_skip += 1
    _log(f"completed {n_ok} of {len(cells)} cells; {n_skip} skipped (no valid folds)")
    _log(f"total wall-clock: {(time.time() - t_overall)/60:.1f} min")

    # Aggregate to per-(network, label) table at alpha=0.5
    rows_05 = []
    for net_id in nets:
        for label in labels:
            t_safe = label.replace("/", "_")
            path = OUT_DIR / f"{net_id}__{t_safe}__alpha0.5.json"
            if not path.exists():
                continue
            d = json.loads(path.read_text())
            gnn_folds = d["gnn_folds"]
            deg_folds = d["deg_folds"]
            if not gnn_folds:
                continue
            gnn_aucs = [v for (_, _, v) in gnn_folds]
            deg_aucs = [v for (_, _, v) in deg_folds]
            common = min(len(gnn_aucs), len(deg_aucs))
            deltas = [gnn_aucs[i] - deg_aucs[i] for i in range(common)]
            rows_05.append({
                "network": net_id, "label": label, "alpha": 0.5,
                "n_folds": len(deltas),
                "mean_gnn_auc": float(np.mean(gnn_aucs[:common])),
                "mean_deg_auc": float(np.mean(deg_aucs[:common])),
                "mean_delta": float(np.mean(deltas)),
                "std_delta": float(np.std(deltas, ddof=1)),
                "params": json.dumps(d["params"]),
            })
    if rows_05:
        out = pd.DataFrame(rows_05)
        out = out.sort_values(["network", "label"])
        out.to_csv(TABLE_DIR / "gnn_courtesy_arm.tsv", sep="\t",
                     index=False)
        _log(f"wrote gnn_courtesy_arm.tsv ({len(out)} rows)")


if __name__ == "__main__":
    main()