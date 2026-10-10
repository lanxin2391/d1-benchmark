"""Aggregate all 72 gnn_courtesy JSONs into a single summary table.

Reads results/runs/gnn_courtesy/<net>__<label>__alpha{0.3,0.5,0.7}.json
(72 files = 6 nets x 4 labels x 3 alphas) and writes
results/tables/gnn_courtesy_arm.tsv with one row per (network, label, alpha)
containing mean Node2Vec AUROC, mean degree AUROC, mean delta, and SD of delta.

This fills in the alpha=0.3 and alpha=0.7 cells that the original
run_gnn_courtesy.py only produced for alpha=0.5.
"""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
JSON_DIR = REPO / "results" / "runs" / "gnn_courtesy"
OUT = REPO / "results" / "tables" / "gnn_courtesy_arm.tsv"

ALPHAS = (0.3, 0.5, 0.7)
NETS = (
    "funmap",
    "intact",
    "reactome",
    "rna_coexp",
    "string_full700",
    "string_phys700",
)
LABELS = (
    "intogen2024",
    "intogen_temporal_new",
    "ot_all",
    "ot_nolit",
)


def _load(p: Path) -> dict | None:
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception as e:
        print(f"  WARN: {p.name} failed to parse: {e!r}")
        return None


def main() -> int:
    rows: list[dict] = []
    missing: list[str] = []
    for net in NETS:
        for label in LABELS:
            for alpha in ALPHAS:
                p = JSON_DIR / f"{net}__{label}__alpha{alpha}.json"
                d = _load(p)
                if d is None:
                    missing.append(p.name)
                    continue
                gnn_folds = d.get("gnn_folds", [])
                deg_folds = d.get("deg_folds", [])
                if not gnn_folds or not deg_folds:
                    missing.append(p.name)
                    continue
                gnn_aucs = [float(v) for (_, _, v) in gnn_folds]
                deg_aucs = [float(v) for (_, _, v) in deg_folds]
                common = min(len(gnn_aucs), len(deg_aucs))
                deltas = [gnn_aucs[i] - deg_aucs[i] for i in range(common)]
                rows.append(
                    {
                        "network": net,
                        "label": label,
                        "alpha": alpha,
                        "n_folds": len(deltas),
                        "mean_gnn_auc": float(np.mean(gnn_aucs[:common])),
                        "mean_deg_auc": float(np.mean(deg_aucs[:common])),
                        "mean_delta": float(np.mean(deltas)),
                        "std_delta": float(np.std(deltas, ddof=1)),
                        "params": json.dumps(d.get("params", {})),
                    }
                )

    if not rows:
        print("No cells found. Aborting.")
        return 1

    out = pd.DataFrame(rows).sort_values(["alpha", "network", "label"]).reset_index(drop=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)

    print(f"Wrote {OUT} ({len(out)} rows, {len(missing)} missing)")
    if missing:
        print("Missing cells:")
        for m in missing:
            print(f"  {m}")
    print()
    print("Headline cells (alpha=0.5):")
    for r in rows:
        if r["alpha"] == 0.5:
            print(
                f"  {r['network']:<16} {r['label']:<22} "
                f"gnn={r['mean_gnn_auc']:.4f} deg={r['mean_deg_auc']:.4f} "
                f"delta={r['mean_delta']:+.4f}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
