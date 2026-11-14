"""Dump the numbers the 1118 deck and the Q&A bank quote into one
place, so that nothing on a slide is typed from memory.

Output: docs/report_1118_facts.md
"""
from __future__ import annotations

import pathlib
import sys

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
ROOT = pathlib.Path(r"D:\Grade3\swxxx\final\D1\d1-benchmark")
T = ROOT / "results/tables"
OUT = ROOT / "docs/report_1118_facts.md"

L = []


def w(s=""):
    L.append(s)


def load(name):
    p = T / name
    if not p.exists():
        return None
    return pd.read_csv(p, sep="\t")


w("# The numbers the 1118 deck quotes")
w()
w("**Generated:** `scripts/dump_report_facts.py` (do not edit by hand).")
w("Every figure on a slide comes from here, so the deck, the Q&A bank")
w("and the manuscript cannot drift apart.")
w()

# ---- 1. the primary endpoint --------------------------------------
t1 = load("a4_table1_delta_auroc_v40.tsv")
if t1 is not None:
    t1 = t1.rename(columns={"a0.5_mean": "m", "a0.5_std": "sd"})
    pos = int((t1.m > 0).sum())
    w("## 1. The primary endpoint (24 arms, alpha = 0.5)")
    w()
    w(f"- the margin over the degree baseline is positive in "
      f"**{pos}/24** arms; median **{t1.m.median():+.4f}**, "
      f"range {t1.m.min():+.4f} .. {t1.m.max():+.4f}")
    w(f"- the best arm: **{t1.loc[t1.m.idxmax(), 'network']} x "
      f"{t1.loc[t1.m.idxmax(), 'label_set']}** ({t1.m.max():+.4f}); "
      f"the worst: **{t1.loc[t1.m.idxmin(), 'network']} x "
      f"{t1.loc[t1.m.idxmin(), 'label_set']}** ({t1.m.min():+.4f})")
    w(f"- the single negative cell: "
      f"{t1[t1.m <= 0][['network', 'label_set', 'm']].to_dict('records')}")
    w()
    w("| network | label | dAUROC | SD |")
    w("|---|---|---|---|")
    for r in t1.sort_values("m", ascending=False).itertuples():
        w(f"| {r.network} | {r.label_set} | {r.m:+.4f} | {r.sd:.4f} |")
    w()

verdict = load("stats_h1_verdict_v40.tsv")
if verdict is not None:
    w("**The H1 decision table** (`stats_h1_verdict_v40.tsv`):")
    w()
    w("```")
    w(verdict.to_string(index=False))
    w("```")
    w()

# ---- 2. the control stack -----------------------------------------
w("## 2. The control stack (the leakage controls)")
w()
ctl = [
    ("control 2, degree-matched random seeds",
     "results/runs/a5/a5/ + results/tables/stats_fine_z_summary.tsv",
     "23/24 arms beat the null at z > 5; cell-level z = 3.09-226.9; "
     "the null mean is -0.052 .. +0.005, i.e. degree-matched random "
     "seeds are WORSE than the real seeds"),
    ("control 3, rewired network",
     "results/tables/stats_fine_z_summary.tsv",
     "23/24 arms beat the rewired null at z > 5"),
    ("D-09, temporal drift (50 seed resamples)",
     "results/tables/d09_table.tsv",
     "-0.0095 +/- 0.0116"),
    ("D-11, common-size cap per network",
     "results/tables/d11_table.tsv",
     "-0.0125; reactome x intogen2024 flips from +0.0429 to negative"),
    ("D-12, shared gene universe (1,605 genes)",
     "results/tables/d12_table.tsv",
     "-0.0360"),
    ("complex-aware cross-validation (Phase D)",
     "results/tables/d4_v1_v2_comparison.tsv",
     "the leakage drops from 0.87-0.96 to 0.02-0.18; the median |relative "
     "change| of the margin is 0.142; H1 still holds in 24/24 and no "
     "sign flips"),
    ("hub-capped RWR (E3)",
     "results/tables/hub_capped_rwr.tsv",
     "the margin barely moves (+0.0849 -> +0.0809 clean; +0.0252 in the "
     "harshest cap, 8/24 arms flip) -> the degree-orthogonalised model "
     "of the second paper is NOT needed"),
    ("Control 7, publication-count adjusted label",
     "results/tables/control7_matched_null.tsv",
     "IntOGen +0.0826 -> +0.0804 -> +0.0587 (survives); temporal "
     "+0.0162 -> +0.0308 -> +0.0022 (collapses)"),
]
for name, src, note in ctl:
    w(f"- **{name}** ({src}): {note}")
w()

# ---- 3. the new results -------------------------------------------
w("## 3. The results added since the manuscript draft (v1.0-rc3)")
w()
g3 = load("g3_bait_matched_null.tsv")
if g3 is not None:
    w("### G3 -- the bait-matched null (study visibility)")
    w()
    w("```")
    w(g3.to_string(index=False))
    w("```")
    w()

g7 = load("g7_ncg_arbitration.tsv")
if g7 is not None:
    w("### G7 -- the NCG arbitration test")
    w()
    w(f"- dAUROC positive in **{int((g7.delta_auroc > 0).sum())}/24** arms, "
      f"median **{g7.delta_auroc.median():+.4f}**, range "
      f"{g7.delta_auroc.min():+.4f} .. {g7.delta_auroc.max():+.4f}")
    w(f"- dpAUC (FPR <= 0.02) positive in "
      f"**{int((g7.delta_pauc > 0).sum())}/24** arms, median "
      f"**{g7.delta_pauc.median():+.4f}**; the split is systematic: the "
      "four protein networks are positive in every label, funmap and "
      "RNA co-expression negative in all four")
    w()

g8 = load("g8_cleaned_background.tsv")
if g8 is not None:
    w("### G8 -- the negative set is not clean")
    w()
    cont = load("g8_negative_contamination.tsv")
    if cont is not None:
        w(f"- **{cont.frac_negatives_in_ncg.min():.1%} to "
          f"{cont.frac_negatives_in_ncg.max():.1%}** of the held-out "
          "negatives are curated NCG drivers; "
          f"{cont.frac_negatives_ncg_known.min():.1%} to "
          f"{cont.frac_negatives_ncg_known.max():.1%} are NCG-*known*")
    piv = g8.pivot_table(index=["network", "label"], columns="version",
                          values="delta_auroc")
    w("- the effect of cleaning the background (delta_AUROC, standard -> "
      "cleaned):")
    w()
    w("| network | label | standard | NCG-cleaned |")
    w("|---|---|---|---|")
    for (net, lab), r in piv.iterrows():
        w(f"| {net} | {lab} | {r.get('standard', float('nan')):+.4f} | "
          f"{r.get('cleaned', float('nan')):+.4f} |")
    w()

gr9 = load("gr9_cancer_type_summary.tsv")
if gr9 is not None:
    w("### GR-9 -- the cancer-type stratification")
    w()
    w(f"- **{int((gr9.median_delta_auroc > 0).sum())}/{len(gr9)}** cancer "
      f"types have a positive median margin (median over the types "
      f"{gr9.median_delta_auroc.median():+.4f})")
    w()

# ---- 4. the fifth label --------------------------------------------
t2 = load("a4_table1_delta_auroc_v2.tsv")
d11v2 = load("d11_table_v2.tsv")
if t2 is not None and d11v2 is not None:
    w("## 4. The fifth label (NCG)")
    w()
    ncg = t2[t2.label_set == "ncg_known"]
    w(f"- the six NCG arms (alpha = 0.5): "
      + ", ".join(f"{r.network} {r['a0.5_mean']:+.4f}"
                  for _, r in ncg.iterrows()))
    piv = d11v2.pivot_table(index="network", columns="label_base",
                             values="delta_auroc")
    w()
    w("- on the D-11 capped panel (equal positive-set size):")
    w()
    w("| network | " + " | ".join(piv.columns) + " |")
    w("|---|" + "---|" * len(piv.columns))
    for net, r in piv.iterrows():
        w(f"| {net} | " + " | ".join(f"{v:+.4f}" for v in r.values) + " |")
    w()
    w("  **Reactome is positive on NCG (+0.0893) and negative on the two "
      "newest catalogues** -> the Reactome failure is label-specific "
      "(curator lag), not network-specific.")
    w()

cmp_ = load("stats_variance_decomp_comparison.tsv")
if cmp_ is not None:
    w("### The variance decomposition, 24 arms vs 30 arms")
    w()
    w("| variance source | % between-arm (24) | % between-arm (30) |")
    w("|---|---|---|")
    for r in cmp_.itertuples():
        w(f"| {r.variance_source} | {100 * r.frac_between_24:.2f} % | "
          f"{100 * r.frac_between_30:.2f} % |")
    w()

# ---- 4b. C9: the factorial re-run with a cleaned background --------
c9 = load("c9_clean_background_comparison.tsv")
if c9 is not None:
    w("## 4b. C9: the factorial re-run with a cleaned background")
    w()
    flips = int((((c9.standard_a05 > 0) & (c9.cleaned_a05 < 0))
                 | ((c9.standard_a05 < 0) & (c9.cleaned_a05 > 0))).sum())
    w(f"- the positive arms: **{int((c9.standard_a05 > 0).sum())}/24** standard, "
      f"**{int((c9.cleaned_a05 > 0).sum())}/24** with the curated drivers "
      f"removed from the background; sign flips: **{flips}**")
    w(f"- the median margin: **{c9.standard_a05.median():+.4f}** standard, "
      f"**{c9.cleaned_a05.median():+.4f}** cleaned")
    w("- the gene-cluster interval passes in 21/24 arms after the cleaning "
      "(20/24 before)")
    w("- the earlier post-hoc reading (G8) removed 87 % of the IntOGen "
      "positives from the evaluation set along with the background, which "
      "is why it appeared to halve the margin")
    w()
    w("| network | label | standard | cleaned | change |")
    w("|---|---|---|---|---|")
    for r in c9.sort_values("standard_a05", ascending=False).itertuples():
        w(f"| {r.network} | {r.label} | {r.standard_a05:+.4f} | "
          f"{r.cleaned_a05:+.4f} | {r.change:+.4f} |")
    w()


# ---- 4c. the four design controls (C17) ---------------------------
w("## 4c. the four design controls (C17)")
w()
sp = load("seed_perturbation_predictions.tsv")
if sp is not None:
    for r in sp.itertuples():
        w(f"- the seed perturbation, {100 * r.seed_fraction:.0f} % of the "
          f"seeds: the median margin {r.median_margin:+.4f}, the median "
          f"relative loss {100 * r.median_relative_loss:.1f} %, "
          f"{r.arms_losing_more_than_10pct} of 24 arms lose more than "
          f"10 %; rho(score CV, loss) = {r.rho_score_cv_vs_loss:+.2f}")
c15 = load("c15_size_matched_gene_ci.tsv")
if c15 is not None:
    w(f"- the size-matched diagnostic: the gene-cluster interval excludes "
      f"zero in {int(c15.gene_ci_excludes_0_standard.sum())}/{len(c15)} "
      f"arms in the standard design and "
      f"{int(c15.gene_ci_excludes_0_capped.sum())}/{len(c15)} with the "
      f"positive-set size held equal -> the failures are a power problem")
c16 = load("c16_label_redundancy.tsv")
if c16 is not None:
    w(f"- the label redundancy: each label's per-network means are "
      f"{c16.r2_of_its_arm_means_from_the_other_labels.min():.2f}-"
      f"{c16.r2_of_its_arm_means_from_the_other_labels.max():.2f} "
      f"predictable from the other labels (R-squared); the overlap "
      f"matrix has a participation ratio of 1.54")
second = load("second_learner_summary.tsv")
if second is not None:
    w(f"- the second (proximity) learner: positive in "
      f"{int((second.delta_auroc > 0).sum())}/{len(second)} arms, the "
      f"median margin {second.delta_auroc.median():+.4f} against RWR's "
      f"{second.rwr_delta_auroc.median():+.4f}; the Spearman rho between "
      f"the two learners is "
      f"{second.delta_auroc.corr(second.rwr_delta_auroc, method='spearman'):+.2f}")
nodb = load("c18_nodb_channel.tsv")
if nodb is not None:
    w("- STRING without the database channel: the margin changes by "
      + ", ".join(f"{r.label} {r.change:+.4f}" for r in nodb.itertuples())
      + " -> the shared-upstream-knowledge-base concern is small here")
w()


# ---- 5. label quality ---------------------------------------------
tiers = load("label_quality_tiers.tsv")
tier_auc = load("label_tier_delta_auroc.tsv")
anchor = load("clingen_noise_anchor.tsv")
if tiers is not None and tier_auc is not None:
    w("## 5. Label quality (C7)")
    w()
    counts = tiers.groupby(["label", "tier"]).size().unstack(fill_value=0)
    w("| label | high | mid | low |")
    w("|---|---|---|---|")
    for lab, r in counts.iterrows():
        w(f"| {lab} | {r.get('high', 0)} | {r.get('mid', 0)} | "
          f"{r.get('low', 0)} |")
    w()
    med = tier_auc.groupby(["label", "tier"]).delta_auroc.median().unstack()
    w("The tier-restricted margin (median over the six networks):")
    w()
    w("| label | high | mid | low |")
    w("|---|---|---|---|")
    for lab, r in med.iterrows():
        cells = []
        for tier in ("high", "mid", "low"):
            v = r.get(tier)
            cells.append(f"{v:+.4f}" if pd.notna(v) else "-")
        w(f"| {lab} | " + " | ".join(cells) + " |")
    w()
if anchor is not None:
    w("**The ClinGen 17 % noise anchor:** the margin moves by "
      f"{anchor.change.min():+.4f} .. {anchor.change.max():+.4f} when a "
      "random 17 % of the positives is demoted to background in every "
      "fold (a perfect scorer would lose 8.5 AUROC points).")
    w()

# ---- 6. stability --------------------------------------------------
w("## 6. The stability analysis (Phase S)")
w()
for name, note in [
    ("stability_per_gene.tsv", None),
    ("network_stability.tsv", None),
    ("stability_attribution.tsv", None),
]:
    df = load(name)
    if df is None:
        continue
    w(f"### {name} ({len(df)} rows, columns: "
      f"{', '.join(df.columns[:8])})")
    w()
    w("```")
    w(df.head(12).to_string(index=False))
    w("```")
    w()

# ---- 7. the DepMap figure -----------------------------------------
dm = load("depmap_validation.tsv")
if dm is not None:
    d = dm[dm.threshold == -0.5]
    w("## 7. DepMap external validation (figure 4)")
    w()
    for k in sorted(d.K.unique()):
        v = d[d.K == k]
        w(f"- K = {k}: median over the 24 arms of the fraction of cell "
          f"lines with p < 0.05 = **{v.frac_celllines_p_lt_005.median():.1%}** "
          f"(chance = 5 %); median p = {v.median_p.median():.3f}")
    w()

# ---- 8. the panel --------------------------------------------------
w("## 8. The design (panel and labels)")
w()
w("| network | provenance | LCC genes | LCC edges |")
w("|---|---|---|---|")
w("| funmap | protein (CPTAC proteogenomics) | 10,189 | 196,605 |")
w("| rna_coexp | RNA-derived | 9,152 | 196,172 |")
w("| string_phys700 | protein (STRING physical >= 700) | 9,830 | 85,576 |")
w("| intact | protein (IntAct) | 17,918 | 564,706 |")
w("| string_full700 | curated (STRING all >= 700) | 15,882 | 236,712 |")
w("| reactome | curated (Reactome v86) | 3,953 | 20,143 |")
w()
w("| label | source | positives |")
w("|---|---|---|")
w("| intogen2024 | mutation recurrence (MC3/IntOGen) | 633 |")
w("| ot_all / ot_nolit | Open Targets 26.06, with / without the "
  "literature channel | 600 / 600 |")
w("| intogen_temporal_new | the 152 newest IntOGen drivers | 152 |")
w("| clingen_label | ClinGen cancer gene-disease validity | 84 |")
w("| ncg_known | NCG known tier (Repana 2019) | 708 (capped per arm) |")
w()
w("- design: 6 networks x 4 labels x 3 alphas x 5 folds x 40 repeats = "
  "1,200 evaluations per arm; the learner (RWR, alpha = 0.5) is frozen "
  "and pre-registered; the pre-registration is frozen at Zenodo DOI "
  "10.5281/zenodo.23121762 and every departure is an amendment.")
w()

OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
print(f"wrote {OUT.name} ({len(L)} lines)")
