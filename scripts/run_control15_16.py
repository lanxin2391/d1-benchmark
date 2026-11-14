"""Control 15 -- why exactly do four arms fail the gene-cluster
interval, and Control 16 -- how many independent label opinions are
there really?

CONTROL 15 (the size-matched diagnostic)

The pre-registered gene-cluster interval excludes zero in 20 of 24 arms,
and the four failures are all `intogen_temporal_new` -- the label with
the smallest positive set (63-150 positives per network against
350-700 for the others). Two explanations are possible: the temporal
label's effect is simply too small, or the gene-cluster interval is
wide because the positive set is small, and every arm with a small
positive set would fail it.

The D-11 capped panel decides between them: there all four labels have
the SAME number of positives inside each arm, so the size is held
constant. If the temporal arms still fail there, the label's effect
size is the reason; if all the arms fail equally, the size is.

CONTROL 16 (the effective number of labels)

Five label sets, but they are not five independent opinions: the
Jaccard matrix says IntOGen and NCG overlap by 0.522 and the two
Open Targets panels by 0.558. Two characterisations:

A. the participation ratio of the eigenvalues of the (centred) overlap
   matrix -- the effective number of independent curation opinions;
B. the principal components of the 6 x 4 matrix of arm means, which
   show which labels travel together; and the R-squared of each label's
   arm means predicted by the other three, a per-label redundancy.

Outputs:
  results/tables/c15_size_matched_gene_ci.tsv
  results/tables/c16_label_redundancy.tsv
  results/tables/c16_label_pca.tsv
  docs/label_size_and_redundancy.md
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

T = ROOT / "results/tables"
SPLITS = ROOT / "data/processed/splits_v40"
NETS = ROOT / "data/processed/networks"
SCORES_CAP = ROOT / "results/scores_v40_cap"
PY = r"E:\anaconda\python.exe"

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
DOC = ROOT / "docs/label_size_and_redundancy.md"


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ------------------------------------------------------------------ 15
def control_15():
    """Re-run the D-11 capped arms with the per-gene dumps, then the
    gene-cluster bootstrap on them."""
    SCORES_CAP.mkdir(parents=True, exist_ok=True)
    runs = ROOT / "results/runs_v40_cap"
    runs.mkdir(parents=True, exist_ok=True)
    # the capped label files are per (network, label)
    for net in NETWORKS:
        for lab in LABELS:
            cap = f"{lab}_cap_{net}"
            sp = SPLITS / f"{cap}__native_{net}.tsv"
            if not sp.exists():
                continue
            r = subprocess.run(
                [PY, "scripts/run_a4.py", "--only-net", net,
                 "--only-label", cap,
                 "--splits-dir", str(SPLITS),
                 "--runs-dir", str(runs),
                 "--scores-dir", str(SCORES_CAP),
                 "--table-suffix", "_cap", "--n-rep", "40",
                 "--alpha", "0.5"],
                capture_output=True, text=True, cwd=str(ROOT))
            if r.returncode != 0:
                print((r.stdout or r.stderr)[-800:])
                raise SystemExit(f"the capped arm {cap} failed")
    log("the capped arms are re-run with the per-gene dumps")

    r = subprocess.run(
        [PY, "scripts/run_gene_cluster_ci.py",
         "--scores-dir", str(SCORES_CAP),
         "--out", str(T / "stats_bootstrap_ci_gene_cap.tsv")],
        capture_output=True, text=True, cwd=str(ROOT))
    log((r.stdout or r.stderr).strip().splitlines()[-1][:90])

    std = pd.read_csv(T / "stats_bootstrap_ci_gene.tsv", sep="\t")
    cap = pd.read_csv(T / "stats_bootstrap_ci_gene_cap.tsv", sep="\t")
    std = std[std.alpha == 0.5]
    cap = cap[cap.alpha == 0.5]
    # the number of positives per arm, from the run tables
    def n_pos_of(runs_dir: pathlib.Path, net: str, lab: str):
        f = runs_dir / f"a4_{net}__{lab}.tsv"
        if not f.exists():
            return np.nan
        d = pd.read_csv(f, sep="\t")
        d = d[(d.method == "rwr") & (d.alpha == 0.5)]
        return int(d.n_test_pos.iloc[0]) if len(d) else np.nan

    rows = []
    for r in cap.itertuples():
        base_lab = r.label.replace("_cap_" + r.network, "")
        s = std[(std.network == r.network) & (std.label == base_lab)]
        rows.append({
            "network": r.network,
            "label": base_lab,
            "n_test_pos_standard": n_pos_of(
                ROOT / "results/runs_v40", r.network, base_lab),
            "gene_ci_excludes_0_standard":
                bool(s.excl0.iloc[0]) if not s.empty else None,
            "gene_ci_width_standard": float(s.width.iloc[0])
                if not s.empty and "width" in s.columns else np.nan,
            "margin_standard": float(s.observed_mean.iloc[0])
                if not s.empty else np.nan,
            "n_test_pos_capped": n_pos_of(
                ROOT / "results/runs_v40_cap", r.network, r.label),
            "gene_ci_excludes_0_capped": bool(r.excl0),
            "gene_ci_width_capped": float(r.width),
            "margin_capped": float(r.observed_mean),
        })
    out = pd.DataFrame(rows)
    out.to_csv(T / "c15_size_matched_gene_ci.tsv", sep="\t", index=False)
    log("wrote the size-matched gene-cluster table")
    return out


# ------------------------------------------------------------------ 16
def control_16():
    """The effective number of independent label opinions."""
    # the five-label matrix (the 30-arm grid) so that the NCG panel is
    # in the redundancy analysis as well. The label id is "ncg_known" in
    # the margin table and "ncg" in the overlap matrix.
    t1 = pd.read_csv(T / "a4_table1_delta_auroc_v2.tsv", sep="\t")
    piv = t1.pivot_table(index="network", columns="label_set",
                         values="a0.5_mean")
    piv = piv.rename(columns={"ncg_known": "ncg"})
    jac = pd.read_csv(T / "c2_label_jaccard_matrix.tsv", sep="\t",
                      index_col=0)
    labels = list(jac.index)

    # A. the participation ratio of the overlap matrix
    J = jac.loc[labels, labels].to_numpy(dtype=float)
    Jc = J - np.eye(len(J))      # the self-overlap is 1 by definition
    w = np.linalg.eigvalsh(Jc)
    w = np.clip(w, 0, None)
    pr = float(w.sum() ** 2 / (w ** 2).sum()) if (w ** 2).sum() > 0 else np.nan
    # the same on the mean-centred overlap
    Jcm = Jc - Jc.mean(axis=0, keepdims=True) - Jc.mean(axis=1, keepdims=True) \
        + Jc.mean()
    wm = np.clip(np.linalg.eigvalsh(Jcm), 0, None)
    prm = float(wm.sum() ** 2 / (wm ** 2).sum()) if (wm ** 2).sum() > 0 \
        else np.nan

    # B. the redundancy of each label: how well the other three predict
    #    its arm means
    rows = []
    for lab in labels:
        if lab not in piv.columns:
            continue
        y = piv[lab].dropna()
        others = [c for c in piv.columns if c != lab and c in labels]
        if lab not in piv.columns:
            continue
        X = piv.loc[y.index, others].to_numpy(dtype=float)
        ok = np.isfinite(X).all(axis=1)
        if ok.sum() < 3:
            r2 = np.nan
        else:
            Xo = np.column_stack([np.ones(ok.sum()), X[ok]])
            beta, *_ = np.linalg.lstsq(Xo, y.to_numpy()[ok], rcond=None)
            resid = y.to_numpy()[ok] - Xo @ beta
            r2 = float(1 - resid.var() / y.to_numpy()[ok].var())
        rows.append({
            "label": lab,
            "mean_jaccard_with_the_others": float(
                np.mean([jac.loc[lab, o] for o in others])),
            "r2_of_its_arm_means_from_the_other_labels": r2,
            "mean_arm_margin": float(piv[lab].mean()),
        })
    red = pd.DataFrame(rows)

    # C. the principal components of the arm means (labels as variables)
    have = [l for l in labels if l in piv.columns]
    X = piv[have].dropna()
    Xc = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(Xc.to_numpy(dtype=float), full_matrices=False)
    var = S ** 2 / (S ** 2).sum()
    load = pd.DataFrame(Vt.T, index=have,
                        columns=[f"PC{i + 1}" for i in range(len(var))])
    pca = pd.DataFrame({
        "component": [f"PC{i + 1}" for i in range(len(var))],
        "variance_share": var,
    })
    load_out = load.reset_index().rename(columns={"index": "label"})

    out = {"participation_ratio_overlap": pr,
           "participation_ratio_centred": prm,
           "n_labels": len(labels),
           "top_eigenvalues": [float(x) for x in w],
           "var_share": [float(x) for x in var],
           "loadings": load_out.to_dict("records"),
           "redundancy": red.to_dict("records")}
    red.to_csv(T / "c16_label_redundancy.tsv", sep="\t", index=False)
    load_out.to_csv(T / "c16_label_pca.tsv", sep="\t", index=False)
    log(f"the effective number of labels: {pr:.2f} (participation ratio), "
        f"{prm:.2f} on the centred overlap")
    return out


def write_doc(c15, c16):
    lines = [
        "# Two questions about the labels: why four arms fail the "
        "gene-cluster interval, and how many label opinions there really "
        "are",
        "",
        "**Date:** 2026-11-12  ",
        "**Scripts:** `scripts/run_control15_16.py`  ",
        "**Tables:** `c15_size_matched_gene_ci.tsv`, "
        "`c16_label_redundancy.tsv`, `c16_label_pca.tsv`",
        "",
        "## Control 15: the size-matched diagnostic",
        "",
        "The pre-registered gene-cluster interval excludes zero in 20 of "
        "24 arms and the four failures are all `intogen_temporal_new`, "
        "the label with the smallest positive set. Two explanations: the "
        "temporal label's effect is too small, or the interval is wide "
        "*because* the positive set is small.",
        "",
        "The D-11 capped panel settles it: there all four labels have "
        "the same number of positives inside each arm.",
        "",
        "| network | label | positives (standard) | gene CI excludes 0 "
        "(standard) | positives (capped) | gene CI excludes 0 (capped) |",
        "|---|---|---|---|---|---|",
    ]
    for r in c15.itertuples():
        lines.append(
            f"| {r.network} | {r.label} | {r.n_test_pos_standard} | "
            f"{'yes' if r.gene_ci_excludes_0_standard else 'no'} | "
            f"{r.n_test_pos_capped} | "
            f"{'yes' if r.gene_ci_excludes_0_capped else 'no'} |")
    n_std = int(sum(bool(x) for x in c15.gene_ci_excludes_0_standard))
    n_cap = int(sum(bool(x) for x in c15.gene_ci_excludes_0_capped))
    lines += [
        "",
        f"**Reading.** Standard design: {n_std}/{len(c15)} arms exclude "
        f"zero. Equal-positive-size panel: {n_cap}/{len(c15)}.",
        "",
        "## Control 16: how many independent label opinions are there?",
        "",
        f"- the participation ratio of the overlap matrix: "
        f"**{c16['participation_ratio_overlap']:.2f}**; on the centred "
        f"overlap: **{c16['participation_ratio_centred']:.2f}** "
        f"(for {c16['n_labels']} labels). The panel carries fewer "
        "independent curation opinions than its label count suggests.",
        "- the eigenvalues of the overlap matrix: "
        + ", ".join(f"{x:.2f}" for x in c16["top_eigenvalues"]),
        "",
        "| label | mean Jaccard with the others | R^2 of its arm means "
        "from the other labels | mean arm margin |",
        "|---|---|---|---|",
    ]
    for r in c16["redundancy"]:
        r2 = r["r2_of_its_arm_means_from_the_other_labels"]
        lines.append(
            f"| {r['label']} | {r['mean_jaccard_with_the_others']:.3f} | "
            f"{'n/a' if r2 != r2 else f'{r2:.2f}'} | "
            f"{r['mean_arm_margin']:+.4f} |")
    lines += [
        "",
        "| label | " + " | ".join(f"PC{i + 1}" for i in range(
            len(c16["var_share"]))) + " |",
        "|---|" + "---|" * len(c16["var_share"]),
    ]
    for row in c16["loadings"]:
        lines.append(f"| {row['label']} | " + " | ".join(
            f"{row[f'PC{i + 1}']:+.2f}" for i in range(len(c16["var_share"])))
            + " |")
    lines += [
        "",
        "## Caveats",
        "",
        "1. The capped panel changes the positive sets, so its margins "
        "are not the primary design's; it is used here only to hold the "
        "size constant.",
        "2. The redundancy measures are descriptive: they say how much "
        "independent information the label panel carries, not how much "
        "the results depend on it.",
        "3. Exploratory family (amendment 4).",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    c15 = control_15()
    c16 = control_16()
    write_doc(c15, c16)
    log("wrote the doc")
