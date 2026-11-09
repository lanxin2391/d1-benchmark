"""Stratified CV split generator (B5, contract C4).

Produces one ``data/processed/splits/{label_id}__{universe_id}.tsv`` file per
(label, universe) pair. Every split file is fully determined by:

* the universe (the list of genes to score)
* the positive set (capped to the universe; everything else is negative)
* a base random seed, written once into decisions.md (D-BASE)
* the number of repeats (default 10) and folds (default 5)

Re-running this module with identical inputs produces byte-identical files. The
split files are released as artefacts rather than regenerated from seeds, per
protocol section 4.7.

Decision D-06 says: ``seeds = training-fold positives; they are not in the
ranked test set.`` The harness enforces that on read. Here we only generate the
gene/repeat/fold/y table.

Decision D-BASE (to be confirmed at S0 with Person A) is the integer seed used
for ``StratifiedKFold(random_state=BASE + r)``; we keep it here as a module
constant for now and write the chosen value into docs/decisions.md.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

SPLIT_COLUMNS = ["gene", "repeat", "fold", "y"]
DEFAULT_BASE_SEED = 20260926   # B's candidate; mark D-BASE in decisions.md
DEFAULT_N_REP = 10
DEFAULT_K = 5


def make_splits(universe, positives, *, base_seed: int = DEFAULT_BASE_SEED,
                n_rep: int = DEFAULT_N_REP, k: int = DEFAULT_K) -> pd.DataFrame:
    """Build a C4 split table.

    Parameters
    ----------
    universe : iterable of str
        Every gene that may appear in any test fold. Genes absent here are
        excluded from the CV entirely.
    positives : iterable of str
        Set of positive genes. A gene is positive iff it appears in *both*
        ``universe`` and ``positives``. Everything else in ``universe`` is a
        negative.
    base_seed, n_rep, k : see DEFAULT_*.

    Returns
    -------
    DataFrame with columns ``gene repeat fold y`` and a stable row order.
    """
    universe = np.array(sorted(set(universe)), dtype=object)
    pos = set(positives) & set(universe)
    y = np.array([g in pos for g in universe], dtype=int)

    rows = []
    for r in range(n_rep):
        skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=base_seed + r)
        for f, (_, test_idx) in enumerate(skf.split(universe, y)):
            test_idx = np.asarray(test_idx)
            for i in test_idx:
                rows.append((str(universe[i]), r, f, int(y[i])))
    df = pd.DataFrame(rows, columns=SPLIT_COLUMNS)
    df = df.sort_values(["repeat", "fold", "gene"], kind="mergesort").reset_index(drop=True)
    return df


def write_splits(df: pd.DataFrame, label_id: str, universe_id: str,
                 out_dir: str | Path = "data/processed/splits") -> Path:
    """Write ``df`` as a C4 file under ``out_dir``.

    Returns the written path.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{label_id}__{universe_id}.tsv"
    df[SPLIT_COLUMNS].to_csv(path, sep="\t", index=False)
    return path


def load_splits(path: str | Path) -> pd.DataFrame:
    """Read a C4 file; returns the DataFrame unchanged."""
    return pd.read_csv(path, sep="\t", dtype=str)


def summarise_splits(df: pd.DataFrame) -> dict:
    """Small per-arm diagnostic summary, useful at S2 to compare the panels."""
    df = df.copy()
    df["y"] = df["y"].astype(int)
    n_total = df.drop_duplicates("gene").shape[0]
    n_pos = int(df.loc[df.y == 1, "gene"].drop_duplicates().shape[0])
    by_fold = df.groupby(["repeat", "fold"]).y.agg(["sum", "size"])
    pos_per_fold = by_fold["sum"]
    return {
        "n_genes": int(n_total),
        "n_positives": n_pos,
        "n_repeats": int(df.repeat.nunique()),
        "n_folds": int(df.fold.nunique()),
        "pos_per_fold_min": int(pos_per_fold.min()),
        "pos_per_fold_max": int(pos_per_fold.max()),
        "pos_per_fold_mean": float(pos_per_fold.mean()),
    }


# ---------------------------------------------------------------------------
# COMPLEX-AWARE SPLITS (Phase D, task D2)
# ---------------------------------------------------------------------------
#
# The leakage channel this closes (Phase D, task D1):
#
#   a held-out positive has a complex-linked gene in the TRAINING fold,
#   and RWR propagates that partner's signal into the held-out gene.
#
# The fix is to keep every complex TOGETHER IN ONE FOLD: then a
# held-out gene's complex partners are held out with it and cannot
# leak into its score. Note the direction -- spreading a complex
# across folds would MAXIMISE the leak, not remove it.
#
# Groups are the connected components of the complex graph (overlapping
# complexes merge, which is the conservative choice), and genes in no
# complex are singleton groups.


def complex_groups(universe, complex_map, max_links_per_gene: int = 10,
                   max_group_size: int = 1000):
    """Connected components of the complex graph restricted to `universe`.

    Parameters
    ----------
    universe : iterable of str
        The genes the split is built over.
    complex_map : mapping gene -> iterable of complex ids
        The gene -> complex membership (see scripts/build_complex_map.py
        for how the real map is built).
    max_links_per_gene : int
        Sparsity guard: a gene keeps at most this many of its
        complexes, the largest ones first. Without it a single huge
        complex (e.g. a ribosome) would merge half the panel into one
        group and make a balanced split impossible.

    Returns
    -------
    list of lists of str
        The groups, each sorted, the list sorted by its first member.
    """
    universe = sorted(set(universe))
    in_universe = set(universe)

    # complex id -> the universe genes it contains
    cpx_members: dict[str, set[str]] = {}
    for gene, cpx_ids in complex_map.items():
        if gene not in in_universe:
            continue
        for c in cpx_ids:
            cpx_members.setdefault(c, set()).add(gene)

    # sparsity: a gene keeps its largest complexes only.
    #
    # The gene -> complexes relation is read from complex_map DIRECTLY.
    # Inverting cpx_members instead (scanning every complex for every
    # gene) is O(n_universe x n_complexes) -- 46 million set lookups on
    # the intact panel -- and made the D3 generation unusably slow.
    gene_cpx: dict[str, list[str]] = {}
    for gene, cids in complex_map.items():
        if gene in in_universe and cids:
            gene_cpx[gene] = list(cids)
    kept: dict[str, list[str]] = {}
    for gene, cids in gene_cpx.items():
        cids.sort(key=lambda c: (-len(cpx_members.get(c, ())), c))
        kept[gene] = cids[:max_links_per_gene]

    # the union must use the SPARSIFIED membership: a complex links only
    # the genes that kept it, otherwise the guard has no effect (a gene
    # that dropped a complex would still be tied to it through the other
    # members)
    cpx_kept: dict[str, list[str]] = {}
    for gene, cids in kept.items():
        for c in cids:
            cpx_kept.setdefault(c, []).append(gene)

    # union-find over the complex links
    parent: dict[str, str] = {g: g for g in universe}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

    # SIZE CAP ON THE GROUPS. The curated complex graph is almost
    # connected: overlapping complexes chain into components of
    # thousands of genes (8,998 of 9,830 on string_phys700). A whole
    # component must go into ONE fold, so without a cap the fold sizes
    # explode (145 vs 8,998 genes) and the degree baseline collapses on
    # the tiny folds -- which made delta-AUROC larger, not smaller, for
    # a reason that has nothing to do with the prior.
    #
    # The fix: process the complexes from the largest down and KEEP a
    # complex only if merging it into the components it touches leaves
    # every component at most max_group_size genes. A complex that would
    # blow the cap is skipped; the link it carried is then not enforced,
    # and the residual leakage that leaves is MEASURED by the D1 audit
    # rather than assumed away.
    comp_size: dict[str, int] = {g: 1 for g in universe}
    for c, members in sorted(cpx_kept.items(),
                             key=lambda kv: (-len(kv[1]), kv[0])):
        members = sorted(members)
        if not members:
            continue
        roots = {find(m) for m in members}
        if sum(comp_size.get(r, 1) for r in roots) > max_group_size:
            continue                      # skip this complex
        base = find(members[0])
        total = comp_size.get(base, 1)
        for other in members[1:]:
            r2 = find(other)
            if r2 == base:
                continue
            total += comp_size.get(r2, 1)
            if base < r2:
                parent[r2] = base
            else:
                parent[base] = r2
                base = r2
        comp_size[base] = total

    comps: dict[str, list[str]] = {}
    for g in universe:
        comps.setdefault(find(g), []).append(g)
    groups = [sorted(v) for v in comps.values()]
    groups.sort(key=lambda g: g[0])
    return groups


def _assign_groups_to_folds(groups, y, k: int, rng):
    """Greedy group assignment that balances the positives per fold.

    Groups are processed in descending order of positive count (ties by
    size, then by the random order) and each is given to the fold that
    currently holds the fewest positives. The tie-break inside the
    near-minimal band is randomised, so the repeats are genuinely
    different partitions instead of 40 copies of one.

    Repair pass: if a fold ends up with NO positive -- which happens
    when there are fewer positive-carrying groups than folds, e.g. the
    temporal label with 12-30 positives spread over fewer than five
    complexes -- one positive member is moved out of the smallest
    positive-carrying group. That is the only way to split a group, it
    is recorded in the returned exceptions, and the leak it re-opens is
    one gene of one complex.
    """
    order = sorted(
        range(len(groups)),
        key=lambda i: (-sum(y[g] for g in groups[i]), -len(groups[i]), i))
    rng.shuffle(order)
    order.sort(key=lambda i: (-sum(y[g] for g in groups[i]),
                              -len(groups[i])))

    fold_pos = [0] * k
    fold_size = [0] * k
    fold_of: dict[str, int] = {}
    group_of: dict[str, int] = {}

    # The assignment must balance BOTH the positives and the fold SIZE.
    # Balancing only the positives (with a band around the current
    # minimum) piles every zero-positive group into the fold that
    # happens to hold the fewest positives -- one fold ended up with
    # 9,030 of 9,830 genes, which destroyed the degree baseline on the
    # tiny folds and inflated delta-AUROC for reasons that had nothing
    # to do with the prior.
    #
    # The cost below is the sum of the two relative imbalances; among
    # the folds within a small tolerance of the best cost the choice is
    # random, so the repeats are genuinely different partitions.
    target_pos = max(sum(y.values()) / k, 1.0)
    target_size = max(sum(len(g) for g in groups) / k, 1.0)
    tol = 0.02

    for gi in order:
        members = groups[gi]
        npos = sum(y[g] for g in members)
        cost = [(fold_pos[f] / target_pos) + (fold_size[f] / target_size)
                for f in range(k)]
        best = min(cost)
        cand = [f for f in range(k) if cost[f] <= best * (1 + tol) + tol]
        f = int(rng.choice(cand))
        for g in members:
            fold_of[g] = f
            group_of[g] = gi
        fold_pos[f] += npos
        fold_size[f] += len(members)

    # ---- repair: every fold must hold at least one positive ----------
    exceptions: list[dict] = []
    moved_genes: set[str] = set()
    for f in range(k):
        if fold_pos[f] > 0:
            continue
        # the positive-carrying group with the fewest positives, and at
        # least two of them so that the donor fold does not become empty
        cands = [(sum(y[g] for g in groups[gi]), len(groups[gi]), gi)
                 for gi in set(group_of.values())
                 if sum(y[g] for g in groups[gi]) >= 2]
        if not cands:
            break
        _, _, gi = min(cands)
        # a gene that was already moved for an earlier empty fold must
        # not be moved again, otherwise the same gene would migrate
        # through every empty fold
        members = [g for g in groups[gi]
                   if y[g] == 1 and g not in moved_genes]
        if not members:
            break
        moved = members[0]
        old = fold_of[moved]
        fold_of[moved] = f
        moved_genes.add(moved)
        fold_pos[old] -= 1
        fold_pos[f] += 1
        fold_size[old] -= 1
        fold_size[f] += 1
        exceptions.append({"gene": moved, "group": gi,
                           "from_fold": old, "to_fold": f,
                           "group_size": len(groups[gi])})
    return fold_of, exceptions


def complex_aware_invariants(df: pd.DataFrame, complex_map,
                              max_links_per_gene: int = 10,
                              groups=None, max_group_size: int = 1000) -> dict:
    """Invariants of a complex-aware split file (for the D3 report).

    Returns the number of groups, the number of groups larger than the
    number of folds, the number of groups that are NOT entirely inside
    one fold (must be 0 except for the documented repairs), and the
    per-fold positive counts.
    """
    universe = sorted(df["gene"].unique())
    if groups is None:
        groups = complex_groups(universe, complex_map, max_links_per_gene,
                             max_group_size)
    per_repeat = []
    n_split_groups = 0
    for r, g in df.groupby("repeat"):
        folds = g.groupby("gene")["fold"].first()
        for members in groups:
            if len(set(folds.loc[members])) > 1:
                n_split_groups += 1
        per_repeat.append(
            g.groupby("fold").y.sum().reindex(range(folds.max() + 1),
                                              fill_value=0).tolist())
    pos = pd.DataFrame(per_repeat)
    return {
        "n_genes": int(len(universe)),
        "n_groups": int(len(groups)),
        "n_groups_larger_than_k": int(sum(
            1 for mem in groups if len(mem) > df["fold"].nunique())),
        "n_groups_split_across_folds": int(n_split_groups),
        "pos_per_fold_min": int(pos.to_numpy().min()),
        "pos_per_fold_max": int(pos.to_numpy().max()),
        "pos_per_fold_mean": float(pos.to_numpy().mean()),
    }


def make_splits_complex_aware(universe, positives, *, complex_map,
                              base_seed: int = DEFAULT_BASE_SEED,
                              n_rep: int = DEFAULT_N_REP,
                              k: int = DEFAULT_K,
                              max_links_per_gene: int = 10,
                              max_group_size: int = 1000,
                              exceptions_out: list | None = None):
    """Group-preserving (complex-aware) CV splits.

    Same output schema as :func:`make_splits` -- columns
    ``gene, repeat, fold, y`` -- with one extra guarantee: **every
    complex group lies entirely in one fold**, which closes the
    "held-out gene has a complex partner in the training fold" leakage
    channel quantified in Phase D task D1.

    With an EMPTY complex map the function delegates to
    :func:`make_splits`, so it reproduces the v1 splits byte for byte.

    The fold sizes are not equal by construction (a complex cannot be
    split), so the per-fold positive counts are reported by
    ``summarise_splits`` rather than assumed.
    """
    universe = sorted(set(universe))
    has_complex = any(complex_map.get(g) for g in universe)
    if not has_complex:
        return make_splits(universe, positives, base_seed=base_seed,
                           n_rep=n_rep, k=k)

    pos = set(positives) & set(universe)
    y = {g: int(g in pos) for g in universe}
    groups = complex_groups(universe, complex_map, max_links_per_gene,
                             max_group_size)

    rows = []
    for r in range(n_rep):
        rng = np.random.default_rng(base_seed + r)
        fold_of, exc = _assign_groups_to_folds(groups, y, k, rng)
        for e in exc:
            e = dict(e)
            e["repeat"] = r
            if exceptions_out is not None:
                exceptions_out.append(e)
        for g in universe:
            rows.append((g, r, int(fold_of[g]), int(y[g])))
    df = pd.DataFrame(rows, columns=SPLIT_COLUMNS)
    df = df.sort_values(["repeat", "fold", "gene"],
                        kind="mergesort").reset_index(drop=True)
    return df
