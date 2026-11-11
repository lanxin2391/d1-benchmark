"""Dump the per-gene RWR scores of the 40-repeat panel (results/scores_v40/).

Why
---
The stability analysis needs, for every gene, the prior score it received in
every cross-validation repeat in which it was in the test set. The harness
(:func:`d1.engine.evaluate.run_arm`) only returns per-fold *metrics*, so this
script re-runs the propagation and stores the score matrix itself.

What is stored
--------------
For every (network, label) arm and every alpha:

* ``P``       -- (n_genes, n_folds) RWR score matrix, float64. Column j is the
                 fold (repeat r, fold f) whose seed set is the training-fold
                 positives of that fold.
* ``genes``   -- the network gene symbols, in the row order of ``P``.
* ``degree``  -- the degree baseline score of every gene.
* ``y``       -- the label of every gene (0/1), constant across folds by
                 contract C4.
* ``fold_repeat`` / ``fold_id`` -- the repeat and fold index of every column.
* ``auroc_rwr`` / ``auroc_degree`` -- the per-fold AUROC of the harness.

The test membership of every fold is **not** stored: it is read back from the
released split file, which stays the single authority for the CV design.

Verification (mandatory)
------------------------
The AUROC recomputed from the dumped matrix must equal the AUROC returned by
the harness :func:`d1.engine.evaluate.run_arm` to 1e-12. The script runs the
harness arm and compares fold by fold; a mismatch aborts the run.

Usage
-----
    python scripts/run_stability_dump.py                 # all 24 arms, alpha 0.5
    python scripts/run_stability_dump.py --alpha 0.3 0.7 # sensitivity alphas
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from d1.engine.baselines import degree_scores          # noqa: E402
from d1.engine.evaluate import evaluate, read_splits, run_arm  # noqa: E402
from d1.engine.rwr import rwr, seed_matrix, transition_T  # noqa: E402
from d1.networks.io import load_network                # noqa: E402

NETWORKS = ("funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome")
LABELS = ("intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new")
NET_DIR = REPO_ROOT / "data" / "processed" / "networks"
SPLIT_DIR = REPO_ROOT / "data" / "processed" / "splits_v40"
OUT_DIR = REPO_ROOT / "results" / "scores_v40"
TABLE_DIR = REPO_ROOT / "results" / "tables"
ALPHA_PRIMARY = 0.5


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_arm(net_id: str, label_id: str, alphas, out_dir: Path):
    """Run the propagation for one arm and write one npz per alpha."""
    net_path = NET_DIR / f"{net_id}.tsv"
    split_path = SPLIT_DIR / f"{label_id}__native_{net_id}.tsv"

    genes, A = load_network(str(net_path))
    sp = read_splits(str(split_path))
    index = pd.Series(np.arange(len(genes)), index=genes)
    sp = sp.assign(idx=index.loc[sp.gene].to_numpy())

    combos, seeds, tests = [], [], []
    for (r, f) in sorted(sp.groupby(["repeat", "fold"]).groups):
        d = sp[sp.repeat == r]
        seeds.append(d.loc[(d.fold != f) & (d.y == 1), "idx"].to_numpy())
        tests.append(d[d.fold == f])
        combos.append((r, f))

    WT = transition_T(A)
    deg = degree_scores(A)
    P0 = seed_matrix(len(genes), seeds)
    fold_repeat = np.array([c[0] for c in combos], dtype=np.int32)
    fold_id = np.array([c[1] for c in combos], dtype=np.int32)

    y = (sp.drop_duplicates("gene").set_index("gene").y.reindex(genes).fillna(0)
         .to_numpy().astype(np.int8))

    # harness reference (independent re-run of the same arm) for the AUROC check
    ref = run_arm(network=(genes, A), splits=sp, network_id=net_id,
                  label_id=label_id, universe_id=f"native_{net_id}", alphas=alphas)

    def _ref_auroc(a: float) -> np.ndarray:
        r = ref[(ref.method == "rwr") & (ref.alpha == a)]
        return r.sort_values(["repeat", "fold"]).auroc.to_numpy()

    out_paths = []
    for a in alphas:
        t0 = time.time()
        P, n_iter = rwr(WT, P0, alpha=a)
        auroc = np.empty(len(combos))
        for j, t in enumerate(tests):
            ti = t.idx.to_numpy()
            m = evaluate(P[ti, j], t.y.to_numpy(), t.gene.to_numpy())
            auroc[j] = m["auroc"]
        deg_auroc = np.empty(len(combos))
        for j, t in enumerate(tests):
            ti = t.idx.to_numpy()
            deg_auroc[j] = evaluate(deg[ti], t.y.to_numpy(), t.gene.to_numpy())["auroc"]

        # mandatory cross-check against the harness
        sel = _ref_auroc(a)
        if sel.size != auroc.size:
            raise AssertionError(f"{net_id} x {label_id} alpha={a}: "
                                 f"{sel.size} harness rows vs {auroc.size} dumped")
        max_gap = float(np.nanmax(np.abs(auroc - sel)))
        if not (max_gap <= 1e-12):
            raise AssertionError(
                f"{net_id} x {label_id} alpha={a}: dumped AUROC differs from the "
                f"harness by {max_gap:.3e} (tolerance 1e-12)")

        path = out_dir / f"{net_id}__{label_id}__alpha{a}.npz"
        np.savez_compressed(
            path,
            P=P.astype(np.float64),
            genes=np.asarray(genes, dtype=object),
            degree=deg.astype(np.float64),
            y=y,
            fold_repeat=fold_repeat,
            fold_id=fold_id,
            auroc_rwr=auroc,
            auroc_degree=deg_auroc,
            n_iter=np.int64(n_iter),
            meta=np.array(json.dumps({
                "network": net_id, "label_set": label_id, "alpha": float(a),
                "n_genes": int(len(genes)), "n_folds": int(len(combos)),
                "n_repeats": int(fold_repeat.max() + 1),
                "splits": str(split_path.relative_to(REPO_ROOT)),
                "splits_sha256": _sha256(split_path),
                "network_sha256": _sha256(net_path),
                "max_auroc_gap_vs_harness": max_gap,
            })))
        out_paths.append(path)
        print(f"  [{net_id} x {label_id} alpha={a}] {P.shape} "
              f"max|AUROC-harness|={max_gap:.1e} ({time.time() - t0:.1f}s)")
    return out_paths


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alpha", type=float, nargs="+", default=[ALPHA_PRIMARY])
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    args = ap.parse_args(argv)

    nets = args.only_net or list(NETWORKS)
    labels = args.only_label or list(LABELS)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    manifest = []
    summary = []
    for label in labels:
        for net in nets:
            t0 = time.time()
            paths = build_arm(net, label, args.alpha, OUT_DIR)
            for p in paths:
                meta = json.loads(str(np.load(p, allow_pickle=True)["meta"]))
                manifest.append({"file": p.name, "bytes": p.stat().st_size,
                                 "sha256": _sha256(p), **meta})
            summary.append({"network": net, "label_set": label,
                            "seconds": round(time.time() - t0, 1)})
            print(f"[{net} x {label}] done in {time.time() - t0:.1f}s")

    man = pd.DataFrame(manifest)
    man_path = TABLE_DIR / "stability_per_gene_scores_manifest.tsv"
    man.to_csv(man_path, sep="\t", index=False, lineterminator="\n")
    print(f"\nWrote {len(man)} score dumps; manifest -> {man_path.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
