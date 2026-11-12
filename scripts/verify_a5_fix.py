"""Verify the two P0-4b fixes by making the two code paths meet.

FIX 1 (resumed arms used the pooled fold deltas instead of the
per-resample means). Test: recompute one arm FRESH through
run_a5._run_one_arm() and compare it with the same arm aggregated from
its checkpoint file -- they must agree elementwise.

FIX 2 (the aggregation loop unpacked one level too few). Test: the
release tables must carry, for every arm, (a) the observed delta of
that exact cell from the 40-repeat primary run, and (b) a null mean /
SD computed from that exact arm's 100 per-resample means. The check
recomputes both from the raw files and compares with the published
tables.

Also checks the alignment assertion: every cell's published row must
come from the trace file with the matching name.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import scripts.run_a5 as a5  # noqa: E402

CHECK_ARM = ("reactome", "intogen_temporal_new")   # cheapest arm (~6 CPU-min)
TOL = 1e-12

ok = True


def report(name, passed, detail=""):
    global ok
    ok = ok and passed
    print(f"[{'PASS' if passed else 'FAIL'}] {name}"
          + (f"  {detail}" if detail else ""))


# ---------------------------------------------------------------- FIX 1
print("=" * 72)
print("FIX 1: fresh path vs resume path must give the same null statistic")
print("=" * 72)
net, lab = CHECK_ARM
print(f"recomputing {net} x {lab} fresh (100 resamples x 200 folds) ...")
a5.PROCESSED_SPLIT = ROOT / "data/processed/splits_v40"
fresh_df, fresh_means = a5._run_one_arm(net, lab, 100)

trace = ROOT / "results/runs/a5" / f"{net}__{lab}__seed_resamples.tsv"
disk = pd.read_csv(trace, sep="\t")
resume_means = disk.groupby("resample")["delta_auroc"].mean().to_numpy()

report("same number of resample means",
       len(fresh_means) == len(resume_means),
       f"fresh={len(fresh_means)} resume={len(resume_means)}")
if len(fresh_means) == len(resume_means):
    d = np.max(np.abs(np.asarray(fresh_means) - resume_means))
    report("per-resample means identical", d < TOL, f"max|diff| = {d:.3e}")

# the fold-level rows must match too
f_sorted = fresh_df.sort_values(["resample", "repeat", "fold"]).reset_index(drop=True)
d_sorted = disk.sort_values(["resample", "repeat", "fold"]).reset_index(drop=True)
d_rows = float(np.max(np.abs(f_sorted["delta_auroc"].to_numpy()
                             - d_sorted["delta_auroc"].to_numpy())))
report("fold-level deltas identical", d_rows < TOL, f"max|diff| = {d_rows:.3e}")

# and the pooled statistic (the BUG) must be visibly different, i.e. the
# bug would have been detectable
pooled = disk["delta_auroc"].to_numpy()
report("the old (buggy) statistic differs from the fixed one",
       abs(np.std(pooled, ddof=1) - np.std(resume_means, ddof=1)) > 1e-4,
       f"sd(pooled)={np.std(pooled, ddof=1):.6f} vs "
       f"sd(per-resample means)={np.std(resume_means, ddof=1):.6f}")

# ---------------------------------------------------------------- FIX 2
print()
print("=" * 72)
print("FIX 2: every published row must come from its own arm's files")
print("=" * 72)
p_tab = pd.read_csv(ROOT / "results/tables/a5_p_values.tsv", sep="\t")
n_tab = pd.read_csv(ROOT / "results/tables/a5_null_distribution.tsv", sep="\t")
a4 = pd.read_csv(ROOT / "results/runs_v40/a4_combined.tsv", sep="\t")
sub = a4[(a4["method"] == "degree") | (a4["alpha"] == 0.5)]
sub = sub[sub["method"].isin(["rwr", "degree"])]
piv = sub.pivot_table(index=["network", "label_set", "repeat", "fold"],
                      columns="method", values="auroc").reset_index()
piv = piv.dropna(subset=["rwr", "degree"])
piv["delta"] = piv["rwr"] - piv["degree"]
obs = piv.groupby(["network", "label_set"])["delta"].mean()

bad_rows = []
for _, row in p_tab.iterrows():
    n, l = row["network"], row["label"]
    fp = ROOT / "results/runs/a5" / f"{n}__{l}__seed_resamples.tsv"
    if not fp.exists():
        bad_rows.append((n, l, "missing trace"))
        continue
    d = pd.read_csv(fp, sep="\t")
    m = d.groupby("resample")["delta_auroc"].mean().to_numpy()
    if abs(m.mean() - row["null_mean"]) > TOL:
        bad_rows.append((n, l, f"null_mean {m.mean():.6f} != {row['null_mean']:.6f}"))
    if abs(m.std(ddof=1) - row["null_sd"]) > TOL:
        bad_rows.append((n, l, f"null_sd {m.std(ddof=1):.6f} != {row['null_sd']:.6f}"))
    if abs(obs.loc[(n, l)] - row["observed_delta"]) > 1e-12:
        bad_rows.append((n, l, f"observed {obs.loc[(n,l)]:.6f} != "
                                  f"{row['observed_delta']:.6f}"))
report(f"all {len(p_tab)} published rows match their own arm's files",
       not bad_rows, f"{len(bad_rows)} mismatches")
for b in bad_rows[:8]:
    print("      ", b)

report("the two control-2 tables agree cell for cell",
       len(p_tab) == len(n_tab) == 24
       and set(zip(p_tab.network, p_tab.label)) == set(zip(n_tab.network, n_tab.label)),
       f"{len(p_tab)} p-value rows / {len(n_tab)} null rows")

# z-score consistency
z = (p_tab["observed_delta"] - p_tab["null_mean"]) / p_tab["null_sd"]
report("published z-scores are the published numbers' own z",
       float(np.max(np.abs(z.to_numpy() - p_tab["z_score"].to_numpy()))) < 1e-9,
       f"max|diff| = {float(np.max(np.abs(z.to_numpy() - p_tab['z_score'].to_numpy()))):.2e}")

# the fine-z table must be built from the same traces
fz = pd.read_csv(ROOT / "results/tables/stats_fine_z.tsv", sep="\t")
report("per-fold z table covers 24 cells x 3 alphas x 200 folds",
       len(fz) == 24 * 3 * 200, f"{len(fz)} rows")

print()
print("=" * 72)
print("RESULT:", "ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
print("=" * 72)
sys.exit(0 if ok else 1)
