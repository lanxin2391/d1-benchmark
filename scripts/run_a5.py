"""A5: Degree-matched seed null (Control 2 from protocol §4.4).

For each (network, label) cell:
    1. The observed seeds = the original positives in the train fold (what
       run_arm uses to seed RWR).
    2. For 100 resamples:
            a. Build a NEW seed set of the same size as the observed seeds,
               where each seed is drawn uniformly at random from all genes
               that have the same degree as the corresponding observed seed
               (degree-matched by quantile, not by exact degree, to avoid
               large jumps when the exact-degree bin is empty).
            b. Run RWR (alpha=0.5) with these new seeds.
            c. Compute AUROC_RWR vs AUROC_degree on the test fold.
            d. Take mean across the 50 folds (5-fold x 10-repeat) as the
               one null delta value for this resample.
    3. Null distribution = 100 mean-delta values.
    4. Compare observed delta to the null mean by a one-sided z-test.

This is the "coarse" version of the protocol's "100 resamples per fold"
control: 100 resamples per cell, each averaged over folds, rather than
100 per fold. Both answers are reported in the same table; the paper
should use the per-fold version (a follow-up run) once the coarse
result confirms gene identity matters.

Inputs (all already on disk):
    data/processed/networks/<net>.tsv                  C2 network file
    data/processed/splits/<label>__native_<net>.tsv    C4 splits
    data/processed/labels/<label>.tsv                  C3 labels (rank only)
    d1/engine/rwr.py                                   rwr, seed_matrix, transition_T
    d1/engine/evaluate.py                              evaluate (AUROC, AUPRC, p@k)

Outputs:
    results/tables/a5_null_distribution.tsv       (network, label, null_mean,
                                                   null_sd, null_q05, q50, q95,
                                                   observed_delta, n_resamples)
    results/tables/a5_p_values.tsv                  (network, label, observed_delta,
                                                   null_mean, null_sd, z_score,
                                                   p_one_sided_observed_gt_null)
    results/runs/a5/<net>_<label>_seed_resamples.tsv per-resample trace (large)

Parallelism:
    24 (network x label) cells are run with joblib's Parallel(n_jobs=4).
    Each worker handles 100 resamples for ONE cell. Within a cell, the
    100 resamples are sequential (so the loaded network + WT matrix can
    be reused). Total time should be ~1/4 of the serial run.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.rwr import transition_T, seed_matrix, rwr  # noqa: E402
from d1.networks.io import load_network  # noqa: E402
from d1.engine.evaluate import evaluate  # noqa: E402

NETWORKS = ["funmap", "rna_coexp", "string_phys700",
            "intact", "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
PROCESSED_NET = ROOT / "data" / "processed" / "networks"
PROCESSED_SPLIT = ROOT / "data" / "processed" / "splits"
RUNS_DIR = ROOT / "results" / "runs" / "a5"
TABLE_DIR = ROOT / "results" / "tables"
RUNS_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

N_RESAMPLES = 100
ALPHA = 0.5
N_JOBS = 4                        # per protocol §5 budget (4 cores)


def _log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _load_split(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"gene": str, "y": int})
    return df


def _build_degree_index(gene_to_degree, in_network):
    """Pre-compute the degree -> candidate-gene lookup for one arm.

    Returns (by_degree, sorted_genes, degrees_sorted, orig_pos) where

      by_degree[d]     : genes of degree d, in the SAME order as
                         gene_to_degree.items() would list them
      sorted_genes     : the usable genes sorted by degree (stable, so
                         genes of equal degree keep their original order)
      degrees_sorted   : the degrees of sorted_genes, as an array
      orig_pos[g]      : the position of g in gene_to_degree, used to
                         restore the original order in the fallback path

    Built once per (network, label) instead of once per observed seed.
    The excluded (observed) seeds are filtered out when a candidate list
    is selected, because the excluded set differs per fold.
    """
    usable = [g for g in gene_to_degree if g in in_network]
    orig_pos = {g: i for i, (g, _) in enumerate(gene_to_degree.items())}
    by_degree: dict[int, list[str]] = {}
    for g in usable:
        by_degree.setdefault(gene_to_degree[g], []).append(g)
    sorted_genes = sorted(usable, key=lambda g: (gene_to_degree[g], orig_pos[g]))
    degrees_sorted = np.array([gene_to_degree[g] for g in sorted_genes],
                              dtype=np.int64)
    return by_degree, sorted_genes, degrees_sorted, orig_pos


def _degree_matched_seeds(observed_seeds, gene_to_degree,
                             seed_indices_in_network, rng,
                             index=None):
    """Sample a seed set of the same size as observed, where each new seed
    has the same degree as the corresponding observed seed. Falls back to
    the nearest available degree if the exact-degree bin is empty.
    Excludes already-observed seeds to avoid trivial identity match.

    Raises if the network has fewer non-observed genes than the
    seed-set size (no valid null can be constructed).

    OPTIMISED: the candidate list comes from a pre-computed degree index
    instead of scanning every gene in the network for every seed. The
    candidate LISTS and the RNG call order are unchanged, so the output
    is bit-identical to the original implementation (verified by
    scripts/verify_a5_equivalence.py).
    """
    all_excluded = set(observed_seeds)
    available_in_net = set(gene_to_degree) - all_excluded
    available_in_net &= seed_indices_in_network
    if len(available_in_net) < len(observed_seeds):
        raise ValueError(
            f"network has only {len(available_in_net)} non-observed genes "
            f"but need {len(observed_seeds)} to build a degree-matched null"
        )
    if index is None:
        index = _build_degree_index(gene_to_degree, seed_indices_in_network)
    by_degree, sorted_genes, degrees_sorted, orig_pos = index

    out = []
    for s in observed_seeds:
        target_deg = gene_to_degree[s]
        candidates = _select_candidates(
            index, gene_to_degree, target_deg, all_excluded, seed_indices_in_network)
        if not candidates:
            continue
        out.append(rng.choice(candidates))
    return out


def _select_candidates(index, gene_to_degree, target_deg, all_excluded,
                       seed_indices_in_network):
    """The candidate list for one seed, exactly as the ORIGINAL
    implementation produced it.

    Rule (verbatim from the original): the list is the genes whose degree
    equals target_deg and which are neither an observed seed nor outside
    the network; if that is empty, the smallest radius r >= 1 whose band
    [target_deg - r, target_deg + r] contains at least one usable gene,
    and the list is every usable gene in that band, in the order
    gene_to_degree lists them.

    The implementation avoids the original's O(r^2) rebuild-the-band-every
    -radius scan: the smallest such radius is
        R* = min over usable g of |deg(g) - target_deg|
    so the nearest usable gene determines it, and the band is then sliced
    ONCE. The nearest usable gene must be searched on BOTH sides at
    every step, otherwise a gene further out on the right can be picked
    while a closer one sits on the left.

    This function is exercised directly by
    scripts/verify_fallback.py, so it is the single source of truth.
    """
    by_degree, sorted_genes, degrees_sorted, orig_pos = index
    all_excluded = set(all_excluded)
    # exact-degree bin, minus the observed seeds
    candidates = [g for g in by_degree.get(target_deg, [])
                  if g not in all_excluded and g in seed_indices_in_network]
    if not candidates:
            # Nearest available degree.
            #
            # The original scanned radius = 1, 2, 3, ... and rebuilt the
            # whole band at every radius. For a physical-interaction
            # network the degree of a driver gene can be hundreds away
            # from the next available gene, so the scan ran hundreds of
            # times over a band that grows with the radius -- O(R^2).
            # That is why the arms with 600+ positives took 235 minutes
            # while the 152-positive arm took 11.
            #
            # The smallest radius whose band contains a usable gene is
            #   R* = |deg(nearest usable gene) - target_deg|
            # and the band at R* is every usable gene with degree in
            # [target - R*, target + R*]. So: find the nearest usable
            # gene (a short walk outwards from the insertion point),
            # then slice the band ONCE.
            # The smallest radius whose band holds a usable gene is
            #   R* = min over usable g of |deg(g) - target_deg|.
            # The nearest usable gene on each side must be found
            # INDEPENDENTLY: the degree gaps are asymmetric (one step
            # right can jump six degrees while one step left moves one),
            # so stopping at the first usable gene found while
            # interleaving the two walks can pick a gene further out
            # than one that sits on the other side -- and then the band
            # is too wide and the candidate list differs from the
            # original implementation.
            n = len(sorted_genes)
            pos = int(np.searchsorted(degrees_sorted, target_deg, "left"))
            d_left = None
            lo_i = pos - 1
            while lo_i >= 0:
                g = sorted_genes[lo_i]
                if g not in all_excluded and g in seed_indices_in_network:
                    d_left = abs(gene_to_degree[g] - target_deg)
                    break
                lo_i -= 1
            d_right = None
            hi_i = pos
            while hi_i < n:
                g = sorted_genes[hi_i]
                if g not in all_excluded and g in seed_indices_in_network:
                    d_right = abs(gene_to_degree[g] - target_deg)
                    break
                hi_i += 1
            dists = [d for d in (d_left, d_right) if d is not None]
            if not dists:
                # no usable gene at any degree: cannot build a null here
                return []
            radius = min(dists)
            lo = int(np.searchsorted(degrees_sorted, target_deg - radius, "left"))
            hi = int(np.searchsorted(degrees_sorted, target_deg + radius, "right"))
            band = [g for g in sorted_genes[lo:hi]
                    if g not in all_excluded and g in seed_indices_in_network]
            candidates = sorted(band, key=lambda g: orig_pos[g])
    return candidates


def _run_one_resample(net_id, label, seed):
    """One (network, label, resample_seed) -> list of (repeat, fold, delta)."""
    rng = np.random.default_rng(seed)
    net_path = PROCESSED_NET / f"{net_id}.tsv"
    df_splits = _load_split(PROCESSED_SPLIT / f"{label}__native_{net_id}.tsv")

    genes, A = load_network(str(net_path))
    WT = transition_T(A)
    gene_to_idx = pd.Series(np.arange(len(genes)), index=genes)
    gene_to_degree = dict(zip(genes, np.asarray(A.sum(axis=1)).ravel().astype(int)))
    net_gene_set = set(genes)

    # OPTIMISED: the degree -> candidates index is built ONCE per arm.
    # It holds every network gene per degree bin; the per-fold exclusion
    # filter is applied when the candidate list is selected. Building it
    # per fold meant the cache never hit.
    degree_index = _build_degree_index(gene_to_degree, net_gene_set)

    fold_deltas = []
    for (repeat, fold), split_df_test_only in df_splits.groupby(
            ["repeat", "fold"], sort=False):
        cell = df_splits[(df_splits.repeat == repeat) & (df_splits.fold != fold)]
        obs_seeds_df = cell[cell.y == 1]
        obs_seeds_gene = list(obs_seeds_df.gene)
        null_seeds_gene = _degree_matched_seeds(
            obs_seeds_gene, gene_to_degree, net_gene_set, rng,
            index=degree_index)
        null_seeds_idx = [int(gene_to_idx[g]) for g in null_seeds_gene
                         if g in gene_to_idx.index]
        if len(null_seeds_idx) == 0:
            continue
        P0 = seed_matrix(len(genes), [null_seeds_idx])
        P, _ = rwr(WT, P0, alpha=ALPHA)
        deg = np.asarray(A.sum(axis=1)).ravel().astype(float)
        test_mask = (split_df_test_only.fold == fold)
        if test_mask.sum() == 0:
            continue
        test_idx = gene_to_idx.reindex(
            split_df_test_only.loc[test_mask, "gene"].tolist()).dropna()
        if len(test_idx) == 0:
            continue
        test_y = split_df_test_only.loc[test_mask].set_index("gene").loc[
            test_idx.index].y.to_numpy()
        test_g = test_idx.index.tolist()
        idx_pos = test_idx.to_numpy().astype(int)
        rwr_m = evaluate(P[idx_pos, 0], test_y, test_g)
        deg_m = evaluate(deg[idx_pos], test_y, test_g)
        if not np.isnan(rwr_m["auroc"]) and not np.isnan(deg_m["auroc"]):
            fold_deltas.append((repeat, fold,
                                rwr_m["auroc"] - deg_m["auroc"]))
    return fold_deltas


def _run_one_arm(net_id, label, n_resamples):
    """Returns (DataFrame of fold-level rows, array of per-resample means)."""
    all_rows = []
    means = []
    for s in range(n_resamples):
        rs_deltas = _run_one_resample(net_id, label, s)
        if not rs_deltas:
            continue
        all_rows.extend([(s, r, f, d) for (r, f, d) in rs_deltas])
        means.append(np.mean([d for (_, _, d) in rs_deltas]))
    all_df = pd.DataFrame(all_rows, columns=["resample", "repeat", "fold",
                                            "delta_auroc"])
    return all_df, np.array(means)


def _observed_delta(net_id, label):
    """Mean observed delta across folds for alpha=0.5 (RWR - degree)."""
    p = ROOT / "results" / "runs" / "a4_combined.tsv"
    if not p.exists():
        return np.nan
    df = pd.read_csv(p, sep="\t")
    sub = df[(df.network == net_id) & (df.label_set == label)].copy()
    if sub.empty:
        return np.nan
    # FIX (pre-registration amendment 2, P0-1): the alpha filter below was
    # DEAD CODE -- `pivot` was built from the unfiltered `sub`, so
    # pivot_table averaged the RWR arm over alpha in {0.3, 0.5, 0.7} while
    # the degree arm (alpha = NaN) stayed single-valued. The released
    # observed_delta was therefore the mean over the three alphas, while
    # the null distribution this function is compared against was
    # generated at alpha = 0.5 only (run_a5.py uses ALPHA = 0.5).
    # Null and observed were built on different alphas.
    sub_rwr_05 = sub[(sub.alpha == 0.5) & (sub.method == "rwr")]
    sub_deg = sub[sub.method == "degree"]
    if sub_rwr_05.empty or sub_deg.empty:
        return np.nan
    pivot = sub_rwr_05.pivot_table(index=["repeat", "fold"], values="auroc",
                                    columns="method")
    pivot["delta"] = pivot["rwr"] - sub_deg.set_index(["repeat", "fold"])["auroc"]
    return float(pivot["delta"].mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-resamples", type=int, default=N_RESAMPLES)
    ap.add_argument("--n-jobs", type=int, default=N_JOBS)
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    ap.add_argument("--splits-dir", default=None,
                    help="directory holding the split files (default "
                         "data/processed/splits; use "
                         "data/processed/splits_v40 for the 40-repeat "
                         "design, work plan A2)")
    ap.add_argument("--force", action="store_true",
                    help="recompute arms even if a checkpoint exists")
    args = ap.parse_args()

    nets = args.only_net or NETWORKS
    labels = args.only_label or LABELS
    cells = [(n, r) for n in nets for r in labels]
    _log(f"A5 degree-matched seed null. {len(nets)} nets x {len(labels)} labels "
         f"x {args.n_resamples} resamples = "
         f"{len(cells) * args.n_resamples} total seed-null runs, "
         f"n_jobs={args.n_jobs}")

    global PROCESSED_SPLIT
    if args.splits_dir:
        PROCESSED_SPLIT = ROOT / args.splits_dir

    t_overall = time.time()
    # PER-ARM CHECKPOINTING: the original wrote nothing until the whole
    # Parallel call returned, so a multi-hour run gave no visibility into
    # progress and a crash lost everything. Each arm now writes its
    # per-resample trace as soon as it finishes, and the driver logs it.
    trace_dir = RUNS_DIR / "a5"
    trace_dir.mkdir(parents=True, exist_ok=True)
    done_files = {}
    pending = []
    for net, lab in cells:
        safe = lab.replace("/", "_")
        fp = trace_dir / f"{net}__{safe}__seed_resamples.tsv"
        if fp.exists() and not args.force:
            done_files[(net, lab)] = fp
        else:
            pending.append((net, lab))

    def _run_and_checkpoint(net, lab):
        t0 = time.time()
        all_df, mean_resamples = _run_one_arm(net, lab, args.n_resamples)
        safe = lab.replace("/", "_")
        fp = trace_dir / f"{net}__{safe}__seed_resamples.tsv"
        all_df.to_csv(fp, sep="\t", index=False)
        _log(f"  [{net} x {lab}] done in {(time.time()-t0)/60:.1f} min "
             f"-> {fp.name}")
        return (net, lab), (all_df, mean_resamples)

    per_cell_results = []
    if pending:
        results = Parallel(n_jobs=args.n_jobs, verbose=0)(
            delayed(_run_and_checkpoint)(net, lab) for net, lab in pending
        )
        per_cell_results.extend(results)
    for net, lab in cells:
        if (net, lab) not in done_files:
            continue
        fp = done_files[(net, lab)]
        _log(f"  [{net} x {lab}] resumed from {fp.name}")
        df = pd.read_csv(fp, sep="\t")
        # the per-resample trace column is `delta_auroc` (written by
        # _run_one_arm), not `delta`
        col = "delta_auroc" if "delta_auroc" in df.columns else "delta"
        per_cell_results.append(((net, lab), (df, df[col].to_numpy())))
    order = {c: i for i, c in enumerate(cells)}
    per_cell_results.sort(key=lambda r: order[r[0]])

    null_summary = []
    p_summary = []
    from scipy.stats import norm as _norm
    for (net, lab), (all_df, mean_resamples) in zip(cells, per_cell_results):
        if len(mean_resamples) == 0:
            _log(f"  [{net} x {lab}] no usable resamples")
            continue
        obs_delta = _observed_delta(net, lab)
        null_mean = float(np.mean(mean_resamples))
        null_sd = float(np.std(mean_resamples, ddof=1))
        null_q05 = float(np.quantile(mean_resamples, 0.05))
        null_q50 = float(np.quantile(mean_resamples, 0.50))
        null_q95 = float(np.quantile(mean_resamples, 0.95))
        z = (obs_delta - null_mean) / null_sd if null_sd > 0 else np.nan
        p_one = 1.0 - _norm.cdf(z) if not np.isnan(z) else np.nan
        null_summary.append({
            "network": net, "label": lab, "n_resamples": len(mean_resamples),
            "null_mean": null_mean, "null_sd": null_sd,
            "null_q05": null_q05, "null_q50": null_q50, "null_q95": null_q95,
            "observed_delta": obs_delta,
        })
        p_summary.append({
            "network": net, "label": lab,
            "observed_delta": obs_delta,
            "null_mean": null_mean, "null_sd": null_sd,
            "n_resamples": len(mean_resamples),
            "z_score": z, "p_one_sided_observed_gt_null": p_one,
        })
        trace_path = RUNS_DIR / f"{net}__{lab}__seed_resamples.tsv"
        all_df.to_csv(trace_path, sep="\t", index=False)
        _log(f"  [{net} x {lab}] n_resamples={len(mean_resamples)}  "
             f"null={null_mean:+.4f}  obs={obs_delta:+.4f}  z={z:+.1f}  p={p_one:.4f}")

    null_df = pd.DataFrame(null_summary)
    p_df = pd.DataFrame(p_summary)
    null_df.to_csv(TABLE_DIR / "a5_null_distribution.tsv", sep="\t",
                   index=False)
    p_df.to_csv(TABLE_DIR / "a5_p_values.tsv", sep="\t", index=False)

    _log(f"all {len(cells)} cells done in {(time.time()-t_overall)/60:.1f} min")
    print("\nNull distribution summary:")
    print(null_df.round(4).to_string(index=False))
    print("\np-values summary:")
    print(p_df.round(4).to_string(index=False))


if __name__ == "__main__":
    main()