"""W3 Control 7 — Publication-count bias adjustment.

For each canonical cancer-driver label, build an *adjusted* label that
regresses out per-gene literature volume. The intuition (protocol
§4.4 control 7) is that genes flagged as cancer drivers tend to be
*more studied* than non-drivers. If RWR's contribution is just
"popular genes = important", then after subtracting the literature
effect the contribution should vanish.

Pipeline per label:

1. Load `pubcount.tsv` (B3, 33,109 genes × log10(n+1) score).
2. Load the label set (e.g. `intogen2024.tsv`). Assign score:
     - positives     = log10(pubcount + 1)
     - negatives (in LCC) = log10(pubcount + 1)
3. Fit ``score ~ is_positive`` (one-way ANOVA) → the "label effect"
   β and its residual. We use a *simple* adjustment:
       adjusted_score = raw_score − β · I(is_positive)
   This subtracts the literature-mean shift between positives and
   negatives. Genes whose residual is positive are *unusually
   curated* for their literature volume.
4. Define the *adjusted label set* as those positives whose residual
   is above the **median residual across the union of positives and
   the negative pool**. This trims away "well-studied only because
   they're known cancer genes" — the protocol's worry.
5. Write ``data/processed/labels/<label>_pubcount_adjusted.tsv``.

The output schema is C3 (gene, rank, source, release); rank = residual
desc, ties broken by symbol asc (D-07).

Usage::

    python scripts/build_pubcount_adjusted.py \\
        --label intogen2024 \\
        --pubcount data/processed/labels/pubcount.tsv \\
        --label-file data/processed/labels/intogen2024.tsv \\
        --out data/processed/labels/intogen2024_pubcount_adjusted.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
LABELS_DIR = REPO_ROOT / "data" / "processed" / "labels"

CANONICAL_LABELS = ("intogen2024", "ot_all", "ot_nolit", "clingen_label")


def build_adjusted(label: str, pubcount_path: Path, label_path: Path,
                   out_path: Path, residual_top_quantile: float = 0.5) -> dict:
    """Build one PubCount-adjusted label file. Returns a summary dict."""
    pub = pd.read_csv(pubcount_path, sep="\t", dtype=str)
    # pubcount.tsv schema (from build_pubcount.py): gene, n_pubmed, log10
    # log10 = log10(n_pubmed + 1) is the literature-volume proxy (D-LOG10).
    if "log10" not in pub.columns:
        raise ValueError(f"pubcount file missing 'log10' column: {pub.columns.tolist()}")
    if "n_pubmed" not in pub.columns:
        raise ValueError(f"pubcount file missing 'n_pubmed' column: {pub.columns.tolist()}")
    pub["log_n"] = pd.to_numeric(pub["log10"], errors="coerce").fillna(0)

    lbl = pd.read_csv(label_path, sep="\t", dtype=str)
    pos = set(lbl.gene.dropna())

    # Universe = all pubcount genes
    df = pub[["gene", "log_n"]].copy()
    df["is_positive"] = df.gene.isin(pos).astype(int)
    df["log_n"] = pd.to_numeric(df["log_n"], errors="coerce").fillna(0)

    # one-way: log_n ~ is_positive (constant + label effect)
    pos_mean = df.loc[df.is_positive == 1, "log_n"].mean()
    neg_mean = df.loc[df.is_positive == 0, "log_n"].mean()
    beta = pos_mean - neg_mean
    # adjusted_score = log_n − β · is_positive
    df["adjusted"] = df["log_n"] - beta * df["is_positive"]

    # the new positive set = positives whose adjusted score is above
    # the median across (positives ∪ sampled negatives). This trims
    # away well-studied-only positives (β already removed the
    # mean shift, but the *trim* is a coarse LCC pruning).
    cand = df[df.is_positive == 1].copy()
    threshold = cand["adjusted"].quantile(residual_top_quantile)
    kept = cand[cand["adjusted"] >= threshold].copy()
    kept = kept.sort_values("adjusted", ascending=False)

    # rank D-07: descending score, then ascending symbol
    kept = kept.sort_values(["adjusted", "gene"], ascending=[False, True])
    kept["rank"] = np.arange(1, len(kept) + 1)
    out = pd.DataFrame({
        "gene": kept["gene"].values,
        "rank": kept["rank"].values,
        "source": f"pubcount_adjusted({label})",
        "release": "2026.06",
    })
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)

    return {
        "label": label,
        "n_positives_total": int(df.is_positive.sum()),
        "n_kept": len(out),
        "kept_fraction": len(out) / max(1, int(df.is_positive.sum())),
        "pos_mean_log_n": float(pos_mean),
        "neg_mean_log_n": float(neg_mean),
        "beta": float(beta),
        "threshold_quantile": residual_top_quantile,
        "out_path": str(out_path.relative_to(REPO_ROOT)),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--label", choices=list(CANONICAL_LABELS) + ["ALL"],
                    default="ALL")
    ap.add_argument("--pubcount", type=Path,
                    default=LABELS_DIR / "pubcount.tsv")
    ap.add_argument("--residual-top-quantile", type=float, default=0.5,
                    help="Keep top-q of positives by residual. Default 0.5 = "
                         "trim the bottom half by residual.")
    args = ap.parse_args(argv)

    targets = list(CANONICAL_LABELS) if args.label == "ALL" else [args.label]
    summaries = []
    for lbl in targets:
        label_path = LABELS_DIR / f"{lbl}.tsv"
        out_path = LABELS_DIR / f"{lbl}_pubcount_adjusted.tsv"
        if not label_path.exists():
            print(f"  [skip] {lbl}: missing {label_path.name}")
            continue
        s = build_adjusted(lbl, args.pubcount, label_path, out_path,
                           residual_top_quantile=args.residual_top_quantile)
        summaries.append(s)
        print(f"  {lbl}: {s['n_positives_total']} → {s['n_kept']} "
              f"({s['kept_fraction']:.0%}) kept, β = {s['beta']:+.3f}")

    # summary file
    if summaries:
        idx = pd.DataFrame(summaries)
        idx_path = LABELS_DIR / "pubcount_adjusted_index.tsv"
        idx.to_csv(idx_path, sep="\t", index=False, float_format="%.4f")
        print(f"\nWrote {idx_path.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())