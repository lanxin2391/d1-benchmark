"""W3 Control 5 — Cross-catalogue label transfer.

For each ordered pair (label_a, label_b) with a != b, this script
materialises the gene-set inputs needed to evaluate
"train RWR on label A's positives, evaluate on label B's positives":

* seeds (label_a positives)  : used as RWR training signal
* target (label_b positives) : the held-out positives we want to recover
* background                 : everything else in the network LCC,
                               EXCLUDING the union of the two labels
                               (so the background carries no information
                               about either label — prevents label
                               leakage from the universe itself)

Output: ``data/processed/transfer/{label_a}__to__{label_b}__native_<net>.tsv``
with the C4-like schema:
    gene, fold, y, label_id, universe_id, transfer_id, role
where:
    y = 1 iff gene in label_b (target) and gene in LCC
    y = 0 iff gene in LCC \ (label_a ∪ label_b)
    role = "seed_or_bg" for the y=0 rows (eligible as seeds in the
            protocol's sense but we mark them with role for clarity)
    role = "target"        for the y=1 rows

We also drop pairs where label_a == label_b (trivial self-transfer)
and pairs where one of the labels is "clingen_label" — clingen is too
small (84 genes) to support 5-fold stratification, so its transfer
cells are noise-dominated. (We can revisit if needed.)

NOTE on folds: cross-label evaluation does NOT benefit from per-fold
splitting in the same way as within-label evaluation — there is exactly
one target set per (label_a, label_b) pair, not one per repeat. We
therefore emit k "folds" by sampling k different background draws
(without replacement) so the run_arm-style evaluator has something to
score. The evaluator (`scripts/transfer_arm.py`) takes this further by
optionally doing bootstrap background draws.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
LABELS_DIR = REPO_ROOT / "data" / "processed" / "labels"
NETWORKS_DIR = REPO_ROOT / "data" / "processed" / "networks"
OUT_DIR = REPO_ROOT / "data" / "processed" / "transfer"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# label set used as the transfer-pair pool. We exclude clingen because
# it is too small (84 genes) for stable 5-fold evaluation and would
# make most transfer cells noise-dominated.
DEFAULT_LABEL_POOL = ("intogen2024", "intogen_temporal_new", "ot_all",
                     "ot_nolit")


def _read_label_genes(label_id: str) -> set[str]:
    p = LABELS_DIR / f"{label_id}.tsv"
    if not p.exists():
        raise FileNotFoundError(p)
    return set(pd.read_csv(p, sep="\t", dtype=str).gene.dropna().astype(str))


def _read_native_universe(net_id: str) -> set[str]:
    p = NETWORKS_DIR / f"{net_id}.tsv"
    if not p.exists():
        raise FileNotFoundError(p)
    df = pd.read_csv(p, sep="\t", dtype=str)
    return set(df.gene_a) | set(df.gene_b)


def _make_transfer_table(label_a: str, label_b: str, network: str,
                         n_folds: int = 5, seed: int = 20260926,
                         background_sample_size: int = 1000
                         ) -> tuple[pd.DataFrame, dict]:
    """Return the transfer eval table for one (a, b, network) cell."""
    if label_a == label_b:
        raise ValueError("self-transfer not allowed (label_a == label_b)")

    genes_a = _read_label_genes(label_a)
    genes_b = _read_label_genes(label_b)
    # accept both "funmap" and "native_funmap" as network names
    net_id = network[len("native_"):] if network.startswith("native_") else network
    lcc = _read_native_universe(net_id)
    # explicit split of the universe
    target = (genes_b - genes_a) & lcc      # b-positives not also in a
    seed_pool = genes_a & lcc                # a-positives in LCC
    background_pool = lcc - (genes_a | genes_b)

    if not target:
        raise ValueError(f"empty target: label_b \\ label_a ∩ LCC empty "
                         f"for {label_a}->{label_b} on {network}")
    if not seed_pool:
        raise ValueError(f"empty seed pool for {label_a} on {network}")

    rng = __import__("numpy").random.default_rng(seed)
    rows = []
    fold = 0
    # one fold = one (background sample, target). Stratify so that the
    # target always appears in every fold.
    bg_total = background_sample_size
    if len(background_pool) < bg_total:
        bg_total = len(background_pool)
    if bg_total == 0:
        raise ValueError(f"no background left after exclusions for "
                         f"{label_a}->{label_b} on {network}")
    # sample without replacement; all folds share the same bg pool then
    # assign indices to folds round-robin for stability.
    bg_indices = rng.choice(len(background_pool), size=bg_total, replace=False)
    bg_genes = [sorted(background_pool)[i] for i in bg_indices]
    per_fold = max(1, bg_total // n_folds)
    for f in range(n_folds):
        # target genes all in fold f (they're the small set; the split
        # only really applies to background).
        bg_slice = bg_genes[f * per_fold:(f + 1) * per_fold]
        for g in target:
            rows.append((g, f, 1, label_b, f"native_{network}",
                         f"{label_a}__to__{label_b}", "target"))
        for g in bg_slice:
            rows.append((g, f, 0, label_b, f"native_{network}",
                         f"{label_a}__to__{label_b}", "background"))

    df = pd.DataFrame(rows, columns=["gene", "fold", "y", "label_id",
                                     "universe_id", "transfer_id", "role"])
    summary = {
        "label_a": label_a,
        "label_b": label_b,
        "network": network,
        "n_target": len(target),
        "n_seed_pool": len(seed_pool),
        "n_background_in_lcc": len(background_pool),
        "n_rows": len(df),
        "n_folds": n_folds,
    }
    return df, summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--labels", nargs="+", default=list(DEFAULT_LABEL_POOL))
    ap.add_argument("--networks", nargs="+",
                    default=["funmap", "string_full700", "string_phys700",
                             "intact", "reactome", "rna_coexp"])
    ap.add_argument("--n-folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20260926)
    ap.add_argument("--background-sample", type=int, default=1000)
    args = ap.parse_args(argv)

    written = []
    summaries = []
    pairs = list(itertools.permutations(args.labels, 2))
    print(f"Building {len(pairs)} ordered label pairs × {len(args.networks)} "
          f"networks = {len(pairs) * len(args.networks)} transfer tables\n")
    for la, lb in pairs:
        for net in args.networks:
            try:
                df, summ = _make_transfer_table(
                    la, lb, net,
                    n_folds=args.n_folds,
                    seed=args.seed,
                    background_sample_size=args.background_sample,
                )
            except (FileNotFoundError, ValueError) as e:
                print(f"  [skip] {la}->{lb} on {net}: {e}")
                continue
            out_path = OUT_DIR / f"{la}__to__{lb}__{net}.tsv"
            df.to_csv(out_path, sep="\t", index=False)
            summaries.append(summ)
            written.append(out_path)
            print(f"  {la:>22} -> {lb:<22} on {net:<16}: "
                  f"target={summ['n_target']:>4}  seed_pool={summ['n_seed_pool']:>5}  "
                  f"bg_in_lcc={summ['n_background_in_lcc']:>5}  rows={summ['n_rows']}")

    print(f"\nWrote {len(written)} transfer tables to "
          f"{OUT_DIR.relative_to(REPO_ROOT)}")

    # index file
    idx = pd.DataFrame(summaries)
    idx.to_csv(OUT_DIR / "transfer_index.tsv", sep="\t", index=False)
    print(f"Wrote {OUT_DIR.relative_to(REPO_ROOT) / 'transfer_index.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
