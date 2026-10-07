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


def _degree_matched_seeds(observed_seeds, gene_to_degree,
                             seed_indices_in_network, rng):
    """Sample a seed set of the same size as observed, where each new seed
    has the same degree as the corresponding observed seed. Falls back to
    the nearest available degree if the exact-degree bin is empty.
    Excludes already-observed seeds to avoid trivial identity match.

    Raises if the network has fewer non-observed genes than the
    seed-set size (no valid null can be constructed)."""
    all_excluded = set(observed_seeds)
    available_in_net = set(gene_to_degree) - all_excluded
    available_in_net &= seed_indices_in_network
    if len(available_in_net) < len(observed_seeds):
        raise ValueError(
            f"network has only {len(available_in_net)} non-observed genes "
            f"but need {len(observed_seeds)} to build a degree-matched null"
        )
    out = []
    for s in observed_seeds:
        target_deg = gene_to_degree[s]
        candidates = [g for g, d in gene_to_degree.items()
                       if d == target_deg
                       and g not in all_excluded
                       and g in seed_indices_in_network]
        if not candidates:
            for radius in range(1, 10000):
                candidates = [g for g, d in gene_to_degree.items()
                               if abs(d - target_deg) <= radius
                               and g not in all_excluded
                               and g in seed_indices_in_network]
                if candidates:
                    break
        out.append(rng.choice(candidates))
    return out


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

    fold_deltas = []
    for (repeat, fold), split_df_test_only in df_splits.groupby(
            ["repeat", "fold"], sort=False):
        cell = df_splits[(df_splits.repeat == repeat) & (df_splits.fold != fold)]
        obs_seeds_df = cell[cell.y == 1]
        obs_seeds_gene = list(obs_seeds_df.gene)
        null_seeds_gene = _degree_matched_seeds(
            obs_seeds_gene, gene_to_degree, net_gene_set, rng)
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
    sub_rwr_05 = sub[(sub.alpha == 0.5) & (sub.method == "rwr")]
    sub_deg = sub[sub.method == "degree"]
    if sub_rwr_05.empty or sub_deg.empty:
        return np.nan
    pivot = sub.pivot_table(index=["repeat", "fold"], values="auroc",
                             columns="method")
    pivot["delta"] = pivot["rwr"] - pivot["degree"]
    return float(pivot["delta"].mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-resamples", type=int, default=N_RESAMPLES)
    ap.add_argument("--n-jobs", type=int, default=N_JOBS)
    ap.add_argument("--only-net", action="append", default=None)
    ap.add_argument("--only-label", action="append", default=None)
    args = ap.parse_args()

    nets = args.only_net or NETWORKS
    labels = args.only_label or LABELS
    cells = [(n, r) for n in nets for r in labels]
    _log(f"A5 degree-matched seed null. {len(nets)} nets x {len(labels)} labels "
         f"x {args.n_resamples} resamples = "
         f"{len(cells) * args.n_resamples} total seed-null runs, "
         f"n_jobs={args.n_jobs}")

    t_overall = time.time()
    per_cell_results = Parallel(n_jobs=args.n_jobs, verbose=0)(
        delayed(_run_one_arm)(net, lab, args.n_resamples) for net, lab in cells
    )

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