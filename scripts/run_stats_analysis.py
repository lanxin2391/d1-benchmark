"""Statistical analysis per protocol §4.6.

Inputs:
    results/runs/a5/<net>__<label>__seed_resamples.tsv
        the per-resample, per-fold null distribution saved by run_a5.py
        (5000 rows per cell: 100 resamples x 50 folds).
    results/runs/a4_combined.tsv
        observed delta_AUROC per (network, label, alpha, repeat, fold).

Outputs:
    results/tables/stats_fine_z.tsv
        per-fold paired z (fine mode) for each (network, label).
    results/tables/stats_variance_decomp.tsv
        mixed-effects variance components (the protocol's actual deliverable).
    results/tables/stats_bootstrap_ci.tsv
        95 % cluster-bootstrap CIs per (network, label) at alpha=0.5.
    results/tables/stats_equivalence.tsv
        TOST equivalence test for the borderline cell
        (reactome x intogen_temporal_new).
"""
import pathlib
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

RUNS_DIR_A5 = ROOT / "results" / "runs" / "a5"
RUN_A4 = ROOT / "results" / "runs" / "a4_combined.tsv"
TABLE_DIR = ROOT / "results" / "tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

NETS = ["funmap", "rna_coexp", "string_phys700",
        "intact", "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]


# ---------------------------------------------------------------------------
# 1. FINE z-score (per-fold paired test, protocol §4.4)
# ---------------------------------------------------------------------------

def _fine_z_one_cell(net: str, label: str, alpha: float = 0.5,
                     a4_path=None, trace_dir=None) -> pd.DataFrame:
    """For each (network, label), compute per-fold paired z using the
    saved null distribution. The fine z for the cell is the mean of
    the per-fold z's (or the z of the paired means -- here we use the
    mean of per-fold z's, which is robust to scale differences).

    The observed fold delta is rwr(auroc) - degree(auroc) for each
    (repeat, fold) at the requested alpha. Degree rows have alpha=NaN
    (degree does not depend on alpha), so we load all rows for the
    (network, label) cell then split by method.

    The null per (repeat, fold) is the 100-resample distribution from
    run_a5.py's saved trace. We pool over the 3 alphas in the trace to
    get 300 null points per (repeat, fold).
    """
    a4 = pd.read_csv(a4_path or RUN_A4, sep="\t")
    obs = a4[(a4.network == net) & (a4.label_set == label)].copy()
    if obs.empty:
        return None
    # delta = rwr_auroc - degree_auroc per (repeat, fold)
    obs["repeat"] = obs["repeat"].astype(int)
    obs["fold"] = obs["fold"].astype(int)
    obs_sub = obs[obs.alpha == alpha]
    if obs_sub.empty:
        return None
    rwr = obs_sub[obs_sub.method == "rwr"].set_index(
        ["repeat", "fold"])["auroc"]
    deg = obs[obs.method == "degree"].set_index(
        ["repeat", "fold"])["auroc"]
    obs_per_fold = (rwr - deg).dropna()   # 50 values per cell
    if len(obs_per_fold) == 0:
        return None

    trace_path = (trace_dir or RUNS_DIR_A5) / f"{net}__{label}__seed_resamples.tsv"
    if not trace_path.exists():
        return None
    null_df = pd.read_csv(trace_path, sep="\t")
    null_df["repeat"] = null_df["repeat"].astype(int)
    null_df["fold"] = null_df["fold"].astype(int)
    # Average over alpha in null (the saved trace has 3 alphas x 50 folds x 100
    # resamples = 15000 rows per cell; for fine z per-fold we use all alphas
    # pooled).
    null_per_fold = (null_df.groupby(["repeat", "fold"]).delta_auroc
                              .agg(["mean", "std"])
                              .reset_index())
    out = []
    for (rep, fold), obs_val in obs_per_fold.items():
        rec_mask = ((null_per_fold.repeat == rep)
                    & (null_per_fold.fold == fold))
        if not rec_mask.any():
            continue
        mu = float(null_per_fold.loc[rec_mask, "mean"].iloc[0])
        sd = float(null_per_fold.loc[rec_mask, "std"].iloc[0])
        z = (obs_val - mu) / sd if sd > 0 else np.nan
        out.append({"network": net, "label": label, "alpha": alpha,
                     "fold": int(fold), "repeat": int(rep),
                     "observed_delta": float(obs_val),
                     "null_mean": mu, "null_sd": sd,
                     "z_fold": z})
    if not out:
        return None
    return pd.DataFrame(out)


def _fine_z_all(a4_path=None, trace_dir=None) -> pd.DataFrame:
    rows = []
    for net in NETS:
        for lab in LABELS:
            for alpha in (0.3, 0.5, 0.7):
                df = _fine_z_one_cell(net, lab, alpha,
                                      a4_path=a4_path, trace_dir=trace_dir)
                if df is not None and len(df):
                    rows.append(df)
    return pd.concat(rows, ignore_index=True)


def _summarize_fine_z(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate per-fold z to per-cell (network, label) z-mean."""
    g = df.groupby(["network", "label", "alpha"]).z_fold.agg(
        ["mean", "std", "count"]
    ).rename(columns={"mean": "z_fine_mean", "std": "z_fine_std",
                       "count": "n_folds"})
    # p-value via the Stouffer combination of independent z's
    from scipy.stats import norm as _norm
    z_combined = g["z_fine_mean"] * np.sqrt(g["n_folds"])
    g["z_combined"] = z_combined
    g["p_one_sided_observed_gt_null"] = 1.0 - _norm.cdf(g["z_combined"])
    return g.reset_index()


# ---------------------------------------------------------------------------
# 2. MIXED-EFFECTS VARIANCE DECOMPOSITION (§4.6, "the actual deliverable")
# ---------------------------------------------------------------------------

def _load_fine_z_or_a4_coarse() -> pd.DataFrame:
    """For the variance decomposition we want fold-level delta_auroc
    values across (network, label, fold). We use the fine_z table when
    present, else fall back to a4_combined."""
    fz_path = TABLE_DIR / "stats_fine_z.tsv"
    if fz_path.exists():
        fz = pd.read_csv(fz_path, sep="\t")
        # observed delta_AUROC per (network, label, repeat, fold)
        out = fz[["network", "label", "alpha", "repeat", "fold", "observed_delta"]]
    else:
        a4 = pd.read_csv(RUN_A4, sep="\t")
        out = a4.copy()
        out["delta"] = out[out.method == "rwr"].set_index(
            ["network", "label_set", "fold"]).auroc - \
            out[out.method == "degree"].set_index(
                ["network", "label_set", "fold"]).auroc
        out = out.reset_index().rename(columns={"label_set": "label",
                                                "delta": "observed_delta"})
        out = out[["network", "label", "repeat", "fold", "observed_delta"]]
        out["alpha"] = 0.5
    # fold (a) as a single "model fold" (5-fold -> one fold index per repeat)
    out = out.copy()
    out["fold_idx"] = out["fold"]   # 0..4
    out["network"] = out["network"].astype("category")
    out["label"] = out["label"].astype("category")
    out["fold_idx"] = out["fold_idx"].astype("category")
    return out


def _variance_decomp() -> pd.DataFrame:
    obs = _load_fine_z_or_a4_coarse()
    obs = obs.dropna(subset=["observed_delta"])
    print(f"variance decomposition on {len(obs)} fold-level observations")

    # MixedLM with random effects for network, label, fold_idx (alpha fixed)
    # Note: only include alpha=0.5 (or average across alphas -- here we take
    # alpha=0.5 because it is the primary endpoint per protocol §4.3).
    sub = obs[obs.alpha == 0.5].copy() if "alpha" in obs.columns else obs
    formula = "observed_delta ~ 1"
    model = smf.mixedlm(formula, sub, groups=sub["fold_idx"])
    # Random effects for network + label via variance_components
    # MixedLM only takes ONE grouping var. Use fold_idx as that; build network
    # and label effects separately.
    fit = model.fit(reml=True)
    print(fit.summary())

    # Variance components
    var_resid = float(fit.scale)
    var_fold = float(fit.cov_re.iloc[0, 0])
    # network + label + network:label via separate simple regressions
    sub2 = obs[obs.alpha == 0.5].copy() if "alpha" in obs.columns else obs
    network_var = sub2.groupby("network").observed_delta.mean().var()
    label_var = sub2.groupby("label").observed_delta.mean().var()
    network_label_var = 0
    # intra/inter by network
    X = sub2.groupby(["network", "label"]).observed_delta.mean()
    inter_var = X.unstack(level=0).fillna(0).values.var()  # rough
    # report
    total = var_resid + var_fold + network_var + label_var
    return pd.DataFrame([{
        "variance_source": "fold (within cell)",
        "variance": float(var_fold),
        "fraction": float(var_fold / total) if total > 0 else np.nan,
    }, {
        "variance_source": "network (between)",
        "variance": float(network_var),
        "fraction": float(network_var / total) if total > 0 else np.nan,
    }, {
        "variance_source": "label (between)",
        "variance": float(label_var),
        "fraction": float(label_var / total) if total > 0 else np.nan,
    }, {
        "variance_source": "residual",
        "variance": float(var_resid),
        "fraction": float(var_resid / total) if total > 0 else np.nan,
    }])


# ---------------------------------------------------------------------------
# 3. CLUSTER BOOTSTRAP CIS (genes as the cluster)
# ---------------------------------------------------------------------------

def _bootstrap_ci(net: str, label: str, alpha: float = 0.5,
                  n_iter: int = 1000, seed: int = 42,
                  unit: str = "repeat") -> dict:
    """Bootstrap CI of the mean delta for one cell.

    FIX (pre-registration amendment 2): the original version did NOT
    filter on alpha, so `pivot_table` averaged the RWR arm over
    alpha in {0.3, 0.5, 0.7} while the degree arm (alpha = NaN) stayed
    single-valued. Every released CI was therefore computed on the mean
    over the three alphas. The alpha filter below keeps the degree rows
    unconditionally (their alpha is NaN) and filters the RWR rows.

    The resampling unit is explicit: unit="repeat" (primary) clusters
    on the repeat, unit="fold" on the individual folds. The gene-cluster
    bootstrap of pre-registration section 9.3 is computed separately by
    scripts/run_gene_cluster_ci.py.
    """
    a4 = pd.read_csv(a4_path or RUN_A4, sep="\t")
    obs = a4[(a4.network == net) & (a4.label_set == label)].copy()
    # FIX: keep the degree rows unconditionally (their alpha is NaN) and
    # filter the RWR rows to the requested alpha.
    obs = obs[(obs.method == "degree") | (obs.alpha == alpha)]
    if obs.empty:
        return {"network": net, "label": label, "alpha": alpha, "unit": unit,
                "ci_low": np.nan, "ci_high": np.nan, "n_iter": 0,
                "n": 0, "mean": np.nan}
    pivot = (obs.pivot_table(index=["repeat", "fold"], values="auroc",
                              columns="method").reset_index())
    pivot["delta"] = pivot["rwr"] - pivot["degree"]
    if unit == "repeat":
        obs_delta = pivot.groupby("repeat")["delta"].mean().to_numpy()
    else:
        obs_delta = pivot["delta"].dropna().to_numpy()
    rng = np.random.default_rng(seed)
    deltas = []
    for _ in range(n_iter):
        idx = rng.integers(0, len(obs_delta), size=len(obs_delta))
        deltas.append(float(np.mean(obs_delta[idx])))
    deltas = np.array(deltas)
    return {
        "network": net, "label": label, "alpha": alpha, "unit": unit,
        "n_iter": n_iter, "n": len(obs_delta),
        "mean": float(np.mean(obs_delta)),
        "ci_low": float(np.quantile(deltas, 0.025)),
        "ci_high": float(np.quantile(deltas, 0.975)),
    }


def _bootstrap_all():
    rows = []
    for net in NETS:
        for lab in LABELS:
            rows.append(_bootstrap_ci(net, lab, 0.5, n_iter=1000, unit="repeat"))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 4. EQUIVALENCE TEST (TOST) for the borderline cell
# ---------------------------------------------------------------------------

def _tost(
    net: str, label: str,
    margin: float = 0.05,
    alpha: float = 0.05,
    rwr_alpha: float = 0.5,
    unit: str = "repeat",
) -> dict:
    """Two one-sided tests for equivalence of zero, within +-margin.

    H0: |delta| >= margin (the effect is at least the margin).
    H1: |delta| <  margin (the effect is indistinguishable from zero).

    FIX (pre-registration amendment 2): the original implementation
    tested `-margin <= mean <= margin` BEFORE the equivalence branch,
    so any cell whose point estimate sat inside the margin was labelled
    BOUND and never reached the equivalence test -- all 24 released
    rows are BOUND for that reason. It also never wrote the
    equivalence interval, so the released table could not be checked
    against pre-registration section 10.

    The TOST rule is: reject H0 when BOTH one-sided p-values are below
    alpha, which is exactly the condition that the (1-2*alpha) CI lies
    entirely inside [-margin, +margin]. The (1-alpha) CI is written out
    as well.

    A cell that is EQUIVALENT to zero is a NULL: there is no evidence
    of an effect in either direction.

    Parameters
    ----------
    alpha     : significance level of the TOST (0.05)
    rwr_alpha : the RWR restart parameter of the arm (0.5, the primary
                endpoint). NOTE this is a different quantity from
                `alpha`; the original implementation reused one name for
                both, which is why it never filtered on it at all.
    """
    a4 = pd.read_csv(a4_path or RUN_A4, sep="\t")
    obs = a4[(a4.network == net) & (a4.label_set == label)].copy()
    # FIX: alpha filter (the degree rows carry alpha = NaN and are kept
    # unconditionally; the RWR rows are filtered to the requested alpha).
    obs = obs[(obs.method == "degree") | (obs.alpha == rwr_alpha)]
    if obs.empty:
        return {"network": net, "label": label, "margin": margin,
                "alpha": alpha, "rwr_alpha": rwr_alpha, "unit": unit,
                "t_low": np.nan, "p_low": np.nan,
                "t_high": np.nan, "p_high": np.nan,
                "verdict": "no data"}
    pivot = (obs.pivot_table(index=["repeat", "fold"], values="auroc",
                              columns="method").reset_index())
    pivot["delta"] = pivot["rwr"] - pivot["degree"]
    if unit == "repeat":
        delta = pivot.groupby("repeat")["delta"].mean().to_numpy()
    else:
        delta = pivot["delta"].dropna().to_numpy()
    n = len(delta)
    mean = float(delta.mean())
    sd = float(delta.std(ddof=1))
    if sd == 0:
        return {"network": net, "label": label, "margin": margin,
                "alpha": alpha, "rwr_alpha": rwr_alpha, "unit": unit,
                "t_low": np.nan, "p_low": np.nan,
                "t_high": np.nan, "p_high": np.nan,
                "verdict": "zero variance"}
    se = sd / np.sqrt(n)
    # Lower test: H0: mean <= -margin; reject when the mean is above it
    t_low = (mean + margin) / se
    p_low = 1 - stats.t.cdf(t_low, df=n - 1)
    # Upper test: H0: mean >= +margin; reject when the mean is below it
    t_high = (mean - margin) / se
    p_high = stats.t.cdf(t_high, df=n - 1)
    reject = (p_low < alpha) and (p_high < alpha)
    t90 = stats.t.ppf(1 - alpha, df=n - 1)
    t95 = stats.t.ppf(1 - alpha / 2, df=n - 1)
    ci90 = (mean - t90 * se, mean + t90 * se)
    ci95 = (mean - t95 * se, mean + t95 * se)
    if reject:
        verdict = (f"EQUIVALENT to zero (mean={mean:.4f}, margin=+-{margin}, "
                   f"90% CI [{ci90[0]:+.4f}, {ci90[1]:+.4f}] inside)")
    else:
        verdict = (f"NOT equivalent to zero (mean={mean:.4f}, margin=+-{margin}, "
                   f"90% CI [{ci90[0]:+.4f}, {ci90[1]:+.4f}])")
    return {"network": net, "label": label, "margin": margin,
            "alpha": alpha, "rwr_alpha": rwr_alpha, "unit": unit,
            "mean": mean, "sd": sd, "n": n,
            "t_low": t_low, "p_low": p_low,
            "t_high": t_high, "p_high": p_high,
            "ci90_low": ci90[0], "ci90_high": ci90[1],
            "ci95_low": ci95[0], "ci95_high": ci95[1],
            "equivalent_to_zero": bool(reject),
            "verdict": verdict}


def _equivalence_all():
    rows = []
    # run TOST on every (network, label) at alpha=0.5
    for net in NETS:
        for lab in LABELS:
            r = _tost(net, lab, margin=0.05, alpha=0.05)
            rows.append(r)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="fine-z",
                    choices=["fine-z", "all"],
                    help="which steps to run. Steps 2-4 (variance "
                         "decomposition, gene-cluster CI, TOST) are "
                         "SUPERSEDED by scripts/run_variance_decomp.py, "
                         "scripts/run_statistics.py and "
                         "scripts/run_ci_variants.py and must not be "
                         "written from this script; 'all' is kept only "
                         "for archaeology and needs --allow-superseded")
    ap.add_argument("--allow-superseded", action="store_true")
    ap.add_argument("--a4", default=None,
                    help="a4_combined.tsv for the observed values "
                         "(default results/runs/a4_combined.tsv; use "
                         "results/runs_v40/a4_combined.tsv for the 40-repeat "
                         "design)")
    ap.add_argument("--a5-dir", default=None,
                    help="directory holding the control-2 per-resample "
                         "traces (default results/runs/a5)")
    args = ap.parse_args()

    a4_path = pathlib.Path(args.a4) if args.a4 else None
    if a4_path is not None and not a4_path.is_absolute():
        a4_path = ROOT / a4_path
    trace_dir = pathlib.Path(args.a5_dir) if args.a5_dir else None
    if trace_dir is not None and not trace_dir.is_absolute():
        trace_dir = ROOT / trace_dir

    print("=" * 70)
    print("§4.6 statistical analysis")
    print("=" * 70)

    print("\n--- 1. fine-mode z-score (per-fold paired) ---")
    fz = _fine_z_all(a4_path=a4_path, trace_dir=trace_dir)
    fz.to_csv(TABLE_DIR / "stats_fine_z.tsv", sep="\t", index=False)
    fz_summary = _summarize_fine_z(fz)
    fz_summary.to_csv(TABLE_DIR / "stats_fine_z_summary.tsv", sep="\t",
                        index=False)
    print(fz_summary.round(4).to_string(index=False))

    if args.only == "fine-z":
        print("\n--- steps 2-4 skipped: they are superseded by "
              "run_variance_decomp.py / run_statistics.py / "
              "run_ci_variants.py (see amendment 2) ---")
        return
    if not args.allow_superseded:
        raise SystemExit("--only all would overwrite the corrected tables; "
                         "pass --allow-superseded only for archaeology")

    print("\n--- 2. mixed-effects variance decomposition ---")
    vd = _variance_decomp()
    vd.to_csv(TABLE_DIR / "stats_variance_decomp.tsv", sep="\t", index=False)
    print(vd.round(6).to_string(index=False))

    print("\n--- 3. cluster bootstrap CIs (genes as cluster) ---")
    bc = _bootstrap_all()
    bc.to_csv(TABLE_DIR / "stats_bootstrap_ci.tsv", sep="\t", index=False)
    print(bc.round(4).to_string(index=False))

    print("\n--- 4. equivalence tests (TOST, margin=0.05) ---")
    eq = _equivalence_all()
    eq.to_csv(TABLE_DIR / "stats_equivalence.tsv", sep="\t", index=False)
    print(eq.round(4).to_string(index=False))


if __name__ == "__main__":
    main()