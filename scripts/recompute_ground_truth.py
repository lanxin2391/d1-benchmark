#!/usr/bin/env python3
"""Recompute every headline number of the paper from the released result tables.

This is the ground-truth check. The manuscript is written by hand in three
Markdown files, and hand-written numbers drift: the 10-repeat design was
replaced by the 40-repeat design, several control tables were regenerated after
bugs were fixed, and the prose kept the old values.

This script does not check the prose (that is `scripts/verify_claims.py`, which
consumes this script's output). It only computes the numbers, from the
released tables, with no hand-entered values. Every number in the paper must
be traceable to a row here.

Usage
-----
    python scripts/recompute_ground_truth.py
    python scripts/recompute_ground_truth.py --out review

Outputs
-------
    review/ground_truth.tsv      one row per (quantity, scope, value, source)
    review/ground_truth.json     the same, keyed
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
TABLES = REPO / "results" / "tables"

ROWS: list[dict] = []


def record(quantity: str, scope: str, value, source: str, note: str = "") -> None:
    ROWS.append({
        "quantity": quantity,
        "scope": scope,
        "value": value,
        "source": source,
        "note": note,
    })


def t(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLES / name, sep="\t")


def r4(x) -> float:
    return float(np.round(float(x), 4))


def r6(x) -> float:
    return float(np.round(float(x), 6))


# ----------------------------------------------------------------- design
def design() -> None:
    combined = pd.read_csv(REPO / "results" / "runs_v40" / "a4_combined.tsv", sep="\t")
    rwr = combined[combined.method == "rwr"]
    record("cv_design", "primary factorial",
           f"{int(rwr.repeat.nunique())} repeats x 5 folds",
           "results/runs_v40/a4_combined.tsv")
    record("cv_design", "folds per arm",
           int(rwr.groupby(["network", "label_set"]).size().iloc[0] / 3),
           "results/runs_v40/a4_combined.tsv",
           "one alpha; the released design has 40 repeats, not the 10 of the "
           "pre-registration")


# ----------------------------------------------------------------- primary endpoint
def primary() -> None:
    a4 = t("a4_table1_delta_auroc_v40.tsv")
    a4 = a4.set_index(["network", "label_set"])
    for (net, lab), r in a4.iterrows():
        record("delta_auroc", f"{net} x {lab}", r6(r["a0.5_mean"]),
               "a4_table1_delta_auroc_v40.tsv")
    vals = a4["a0.5_mean"]
    record("delta_auroc", "panel range (min)", r6(vals.min()),
           "a4_table1_delta_auroc_v40.tsv")
    record("delta_auroc", "panel range (max)", r6(vals.max()),
           "a4_table1_delta_auroc_v40.tsv")
    neg = a4[vals < 0]
    record("delta_auroc", "n negative cells", int(len(neg)),
           "a4_table1_delta_auroc_v40.tsv",
           "; ".join(f"{i[0]} x {i[1]}" for i in neg.index))
    record("delta_auroc", "n positive cells", int((vals > 0).sum()),
           "a4_table1_delta_auroc_v40.tsv")
    for a in ("0.3", "0.7"):
        record("delta_auroc_alpha", f"alpha={a} panel mean",
               r6(a4[f"a{a}_mean"].mean()), "a4_table1_delta_auroc_v40.tsv")
    # per-label means (the numbers the manuscript's label table quotes)
    for lab, g in a4.groupby(level="label_set"):
        record("label_mean", lab, r6(g["a0.5_mean"].mean()),
               "a4_table1_delta_auroc_v40.tsv")
    record("label_mean", "panel median cell", r6(vals.median()),
           "a4_table1_delta_auroc_v40.tsv")
    record("label_mean", "panel mean cell", r6(vals.mean()),
           "a4_table1_delta_auroc_v40.tsv")


# ----------------------------------------------------------------- statistics
def statistics() -> None:
    rs = t("results_stats.tsv")
    record("n_primary_cells", "results_stats.tsv", int(len(rs)),
           "results_stats.tsv")
    record("BH_FDR", "cells rejecting at q<0.05 (fold unit)",
           int(rs["reject_H0_at_q_0.05"].sum()), "results_stats.tsv")
    record("BH_FDR", "cells rejecting at q<0.05 (repeat unit)",
           int(rs["reject_H0_at_q_0.05_repeat"].sum()), "results_stats.tsv")
    fail = rs[~rs["reject_H0_at_q_0.05_repeat"]]
    record("BH_FDR", "failing cell(s)",
           "; ".join(f"{r.network} x {r.label}" for r in fail.itertuples()),
           "results_stats.tsv")

    vd = t("stats_variance_decomp.tsv")
    for r in vd.itertuples():
        record("variance_fraction_between_arm", r.variance_source,
               r6(r.fraction_of_between_arm_spread), "stats_variance_decomp.tsv")
        record("variance_fraction_total", r.variance_source,
               r6(r.fraction_of_total_variance), "stats_variance_decomp.tsv")

    ci = t("stats_bootstrap_ci.tsv").set_index(["network", "label"])
    for unit in ("fold_pct", "repeat_pct", "repeat_t", "repeat_bca", "gene"):
        col = f"excl0_{unit}"
        n = int(ci[col].sum())
        record("H1_CI_excludes_zero", f"{unit} unit", f"{n}/24",
               "stats_bootstrap_ci.tsv",
               "; ".join(f"{i[0]} x {i[1]}" for i in ci[~ci[col].astype(bool)].index))
    rec = ci.loc[("reactome", "intogen_temporal_new")]
    for unit in ("fold_pct", "repeat_pct", "gene"):
        record("negative_cell_CI", f"reactome x intogen_temporal_new ({unit})",
               f"[{r6(rec[f'ci_low_{unit}'])}, {r6(rec[f'ci_high_{unit}'])}]",
               "stats_bootstrap_ci.tsv",
               f"excludes zero: {bool(rec[f'excl0_{unit}'])}")

    eq = t("stats_equivalence.tsv")
    eqv = eq[eq.equivalent_to_zero.astype(bool)]
    record("TOST", "cells equivalent to zero", int(len(eqv)),
           "stats_equivalence.tsv",
           "; ".join(f"{r.network} x {r.label}" for r in eqv.itertuples()))
    record("TOST", "cells not equivalent to zero", int(len(eq) - len(eqv)),
           "stats_equivalence.tsv")

    verdict = t("stats_h1_verdict_v40.tsv")
    for r in verdict.itertuples():
        record("H1_rule", r.interval, f"{r.n_excl0}/{r.n_arms}",
               "stats_h1_verdict_v40.tsv",
               "rule met" if r.rule_met else f"NOT met: {r.failing_arms}")


# ----------------------------------------------------------------- controls
def controls() -> None:
    a4 = t("a4_table1_delta_auroc_v40.tsv")
    a5 = t("a5_p_values.tsv")
    record("control2_degree_matched_seed_null", "z range",
           f"{r4(a5.z_score.min())} to {r4(a5.z_score.max())}", "a5_p_values.tsv")
    record("control2_degree_matched_seed_null", "n cells with z>5",
           int((a5.z_score > 5).sum()), "a5_p_values.tsv")
    record("control2_degree_matched_seed_null", "n null means <= 0",
           int((a5.null_mean <= 0).sum()), "a5_p_values.tsv")
    record("control2_degree_matched_seed_null", "null mean range",
           f"{r4(a5.null_mean.min())} to {r4(a5.null_mean.max())}",
           "a5_p_values.tsv")

    w5 = t("w5_p_values.tsv")
    record("control3_rewiring_null", "z range",
           f"{r4(w5.z_score.min())} to {r4(w5.z_score.max())}", "w5_p_values.tsv")
    record("control3_rewiring_null", "n cells with z>5",
           int((w5.z_score > 5).sum()), "w5_p_values.tsv")
    record("control3_rewiring_null", "n cells with z<=5",
           int((w5.z_score <= 5).sum()), "w5_p_values.tsv")
    per_net = w5.groupby("network").n_rewire.agg(["min", "max"])
    for net, r in per_net.iterrows():
        record("control3_rewiring_null", f"{net}: n rewirings",
               int(r["max"]), "w5_p_values.tsv")
    record("control3_rewiring_null", "total rewired-network runs",
           int(w5.n_rewire.sum()), "w5_p_values.tsv")

    d12 = t("d12_table.tsv")
    record("control4_shared_universe", "shared genes",
           int(d12.shared_genes.iloc[0]), "d12_table.tsv")
    record("control4_shared_universe", "n_test",
           int(d12.n_test.iloc[0]), "d12_table.tsv")
    record("control4_shared_universe", "delta range",
           f"{r4(d12.delta_auroc.min())} to {r4(d12.delta_auroc.max())}",
           "d12_table.tsv")
    record("control4_shared_universe", "n positive cells",
           int((d12.delta_auroc > 0).sum()), "d12_table.tsv")

    tr = t("results_stats_transfer.tsv")
    record("control5_cross_catalogue", "n transfer cells", int(len(tr)),
           "results_stats_transfer.tsv")
    record("control5_cross_catalogue", "n positive cells",
           int((tr.mean_delta_auroc > 0).sum()), "results_stats_transfer.tsv")
    record("control5_cross_catalogue", "n negative cells",
           int((tr.mean_delta_auroc < 0).sum()), "results_stats_transfer.tsv",
           "; ".join(f"{r.network} {r.label_a}->{r.label_b}"
                     for r in tr[tr.mean_delta_auroc < 0].itertuples()))
    record("control5_cross_catalogue", "mean delta",
           r4(tr.mean_delta_auroc.mean()), "results_stats_transfer.tsv")
    within = a4["a0.5_mean"].mean()
    retention = 100.0 * tr.mean_delta_auroc.mean() / within
    record("control5_cross_catalogue", "retention percent",
           round(float(retention), 2), "results_stats_transfer.tsv",
           f"cross-catalogue mean / within-label mean; within-label mean = "
           f"{r4(within)}")

    d09 = t("d09_table.tsv")
    d09 = d09[d09.alpha == 0.5]
    record("control6_temporal_split", "delta range",
           f"{r4(d09.delta_auroc_mean.min())} to "
           f"{r4(d09.delta_auroc_mean.max())}", "d09_table.tsv")
    record("control6_temporal_split", "n positive cells",
           int((d09.delta_auroc_mean > 0).sum()), "d09_table.tsv")
    record("control6_temporal_split", "n effective seed sets (min)",
           int(d09.n_eff_seed_sets.min()), "d09_table.tsv")
    for r in d09.itertuples():
        record("control6_temporal_split", f"{r.network} delta",
               r6(r.delta_auroc_mean), "d09_table.tsv")

    d11 = t("d11_table.tsv")
    record("control_d11_cap", "N per network",
           "; ".join(f"{r.network}={int(r.N)}" for r in d11.itertuples()),
           "d11_table.tsv")

    c7 = t("control7_matched_null.tsv")
    record("control7_publication_count", "n cells", int(len(c7)),
           "control7_matched_null.tsv")
    record("control7_publication_count", "n labels",
           int(c7.label.nunique()), "control7_matched_null.tsv")
    record("control7_publication_count", "mean matched-minus-standard",
           r4(c7.matched_minus_standard.mean()), "control7_matched_null.tsv")

    cx = t("d4_v1_vs_v2_comparison.tsv")
    record("control_complex_aware_cv", "median |relative change|",
           r4(cx.abs_rel_change.median()), "d4_v1_vs_v2_comparison.tsv")
    record("control_complex_aware_cv", "n sign flips",
           int((cx.v1_mean * cx.v2_mean < 0).sum()),
           "d4_v1_vs_v2_comparison.tsv")


# ---------------------------------------------------------------- new controls
def control8() -> None:
    p = TABLES / "control8_seed_degree_summary.tsv"
    if not p.exists():
        return
    s = pd.read_csv(p, sep="\t")
    for stratum in ("low", "mid", "high"):
        record("control8_seed_degree", f"mean {stratum}",
               r4(s[stratum].mean()), "control8_seed_degree_summary.tsv")
        record("control8_seed_degree", f"n positive {stratum}",
               int((s[stratum] > 0).sum()), "control8_seed_degree_summary.tsv")
    record("control8_seed_degree", "mean high-minus-low",
           r4(s["high_minus_low"].mean()), "control8_seed_degree_summary.tsv")


def control9() -> None:
    m = TABLES / "control9_size_model.tsv"
    a = TABLES / "control9_size_adjustment.tsv"
    if not (m.exists() and a.exists()):
        return
    mod = pd.read_csv(m, sep="\t").set_index("term")
    for term in ("log10_edges", "coverage"):
        if term in mod.index:
            record("control9_size", f"{term} coef",
                   r4(mod.loc[term, "coef"]), "control9_size_model.tsv")
            record("control9_size", f"{term} p",
                   float(f"{mod.loc[term, 'p']:.3g}"), "control9_size_model.tsv")
            record("control9_size", f"{term} partial R2",
                   r4(mod.loc[term, "r2_partial"]), "control9_size_model.tsv")
    df = pd.read_csv(a, sep="\t")
    if "spearman_rho_pooled_vs_log10_edges" in df.columns:
        record("control9_size", "pooled Spearman rho",
               r4(df["spearman_rho_pooled_vs_log10_edges"].iloc[0]),
               "control9_size_adjustment.tsv")
        record("control9_size", "permutation p",
               round(float(df["spearman_permutation_p"].iloc[0]), 2),
               "control9_size_adjustment.tsv")
    record("control9_size", "adjusted positive cells",
           int((df["fitted_delta_auroc"] > 0).sum()), "control9_size_adjustment.tsv")
    record("control9_size", "raw positive cells",
           int((df["delta_auroc"] > 0).sum()), "control9_size_adjustment.tsv")
    record("control9_size", "sign preserved",
           int((np.sign(df["delta_auroc"]) == np.sign(df["fitted_delta_auroc"])).sum()),
           "control9_size_adjustment.tsv")
    record("control9_size", "mean abs residual",
           r4(df["residual_delta_auroc"].abs().mean()), "control9_size_adjustment.tsv")


def control10() -> None:
    p = TABLES / "control10_random_input_summary.tsv"
    if not p.exists():
        return
    s = pd.read_csv(p, sep="\t")
    for col, name in [
        ("n_persistent_topk_genes__real", "persistent genes, real input"),
        ("n_persistent_topk_genes__random", "persistent genes, random input"),
        ("mean_degree_percentile_of_topk__real", "degree percentile, real"),
        ("mean_degree_percentile_of_topk__random", "degree percentile, random"),
        ("mean_hub_overlap_topk__real", "hub overlap, real"),
        ("mean_hub_overlap_topk__random", "hub overlap, random"),
        ("input_sensitivity_spearman_real_vs_random", "input sensitivity Spearman"),
        ("input_sensitivity_topk_jaccard_real_vs_random", "input sensitivity Jaccard"),
    ]:
        if col in s.columns:
            record("control10_random_input", name, r4(s[col].mean()),
                   "control10_random_input_summary.tsv")


# ----------------------------------------------------------------- provenance
def provenance() -> None:
    ps = t("provenance_summary.tsv")
    ps = ps[ps.alpha == 0.5]
    for r in ps.itertuples():
        record("provenance_eta2", r.label, r6(r.eta_squared_class),
               "provenance_summary.tsv")
        record("provenance_kruskal_p", r.label, r6(r.p_kruskal),
               "provenance_summary.tsv")


# ----------------------------------------------------------------- GNN arm
def gnn() -> None:
    g = t("gnn_courtesy_arm.tsv")
    a4 = t("a4_table1_delta_auroc_v40.tsv").set_index(["network", "label_set"])
    record("gnn_arm", "n cells", int(len(g)), "gnn_courtesy_arm.tsv")
    record("gnn_arm", "mean GNN delta", r4(g.mean_delta.mean()),
           "gnn_courtesy_arm.tsv")
    record("gnn_arm", "n GNN cells positive", int((g.mean_delta > 0).sum()),
           "gnn_courtesy_arm.tsv")
    rwr_mean = float(np.mean([a4.loc[(r.network, r.label), "a0.5_mean"]
                               for r in g.itertuples()]))
    record("gnn_arm", "mean RWR delta (same 24 cells)", r4(rwr_mean),
           "a4_table1_delta_auroc_v40.tsv")
    for r in g.itertuples():
        record("gnn_arm", f"{r.network} x {r.label}: GNN delta",
               r6(r.mean_delta), "gnn_courtesy_arm.tsv")
        record("gnn_arm", f"{r.network} x {r.label}: RWR delta",
               r6(a4.loc[(r.network, r.label), "a0.5_mean"]),
               "a4_table1_delta_auroc_v40.tsv")


# ----------------------------------------------------------------- networks
def networks() -> None:
    for meta in sorted((REPO / "data" / "processed" / "networks").glob("*.meta.json")):
        m = json.loads(meta.read_text(encoding="utf-8"))
        record("network_nodes", m["net_id"], int(m["n_nodes_before_lcc"]),
               f"data/processed/networks/{m['net_id']}.meta.json",
               "nodes before the LCC restriction")
        record("network_lcc_nodes", m["net_id"], int(m["n_nodes_lcc"]),
               f"data/processed/networks/{m['net_id']}.meta.json")
        record("network_lcc_edges", m["net_id"], int(m["n_edges_lcc"]),
               f"data/processed/networks/{m['net_id']}.meta.json")


# ----------------------------------------------------------------- labels
def labels() -> None:
    lab_dir = REPO / "data" / "processed" / "labels"
    base = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new",
            "clingen_label"]
    sets = {}
    for lab in base:
        p = lab_dir / f"{lab}.tsv"
        if not p.exists():
            continue
        sets[lab] = set(pd.read_csv(p, sep="\t", dtype=str).gene)
        record("label_positives", lab, len(sets[lab]),
               f"data/processed/labels/{lab}.tsv")
    for a in base:
        for b in base:
            if a >= b or a not in sets or b not in sets:
                continue
            record("label_overlap", f"{a} x {b}", len(sets[a] & sets[b]),
                   "data/processed/labels/*.tsv")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(REPO / "review"))
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    design()
    primary()
    statistics()
    controls()
    control8()
    control9()
    control10()
    provenance()
    gnn()
    networks()
    labels()

    df = pd.DataFrame(ROWS)
    df.to_csv(out / "ground_truth.tsv", sep="\t", index=False, lineterminator="\n")
    (out / "ground_truth.json").write_text(
        json.dumps(ROWS, indent=1, default=str), encoding="utf-8")

    print(f"{len(df)} ground-truth quantities written to {out/'ground_truth.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
