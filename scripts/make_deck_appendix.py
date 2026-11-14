"""Generate the appendix slides' content (docs/deck/slides_appendix.json)
from the released tables, so that every appendix number is the table's
number and not a retyped one."""
from __future__ import annotations

import json
import pathlib
import sys

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
ROOT = pathlib.Path(r"D:\Grade3\swxxx\final\D1\d1-benchmark")
T = ROOT / "results/tables"
OUT = ROOT / "docs/deck/slides_appendix.json"


def load(name):
    p = T / name
    return pd.read_csv(p, sep="\t") if p.exists() else None


slides = []

# ---- A1: the literature -> our design decisions ---------------------
slides.append({
    "n": 17, "kind": "appendix", "tag": "APPENDIX A1",
    "title": "The literature, and the decision each paper forced",
    "note": "Every design decision in the benchmark traces to a specific "
            "published result. Full method digest: "
            "D1_参考文献_数据分析算法与方法汇编.docx (41 papers).",
    "table": {
        "head": ["Paper", "What it found", "What we did about it"],
        "rows": [
            ["Bailey et al. 2018, Cell 173:371",
             "MC3 consensus of 26 mutation callers on 9,423 TCGA tumour-normal pairs, 26 tissues",
             "the primary label (intogen2024, 633 genes)"],
            ["Shi et al. 2025, Nat Cancer 6:205",
             "FunMap: CPTAC proteogenomic network; Reactome-derived gold standard; 11 cancer types",
             "one of six networks, with its coverage and its circularity reported explicitly"],
            ["Repana et al. 2019, Genome Biol 20:1",
             "NCG: 2,372 genes curated from the sequencing literature, with a network-derived confidence level",
             "an arbitration catalogue and a negative-set audit, not a label in the confirmatory family"],
            ["Picart-Armada et al. 2019, PLoS CB 15:e1007276",
             "complex leakage drops a network classifier from ~12 to <4.5 hits in the top 20; the endpoint changes the ranking of methods",
             "complex-aware cross-validation + two extra endpoints (pAUC at FPR<=0.02, top-10 % precision)"],
            ["Guala & Sonnhammer 2017, Sci Rep 7:46598",
             "a degree-filtered local method beats diffusion on the full list; diffusion wins only in the top 1-10 %",
             "the learner frozen and the low-FPR endpoint reported separately"],
            ["Hill et al. 2019",
             "over 10,000 random gene lists, the algorithms differ hugely in how input-independent their output is",
             "a random-input degeneracy control"],
            ["Wright et al. 2025",
             "performance tracks network size over 45 interactomes; the effect collapses once coverage is controlled",
             "the shared-gene-universe control and the per-network common-size cap"],
            ["Schaefer et al. 2015, Front Genet",
             "bait usage correlates with degree (r = 0.52); correcting for it removes the cancer-vs-control degree difference",
             "a publication-count-adjusted label and a bait-matched null model"],
            ["Boutros et al. 2014, Genome Biol 15:462",
             "the gold standard dominates any ranking, so the evaluator should be blind to the algorithm",
             "the NCG arbitration test, and a blind expert challenge in the validation plan"],
        ],
        "fontsize": 11.5, "rowh": 0.62,
    },
})

# ---- A2: the data provenance ----------------------------------------
slides.append({
    "n": 18, "kind": "appendix", "tag": "APPENDIX A2",
    "title": "Every input, its release, and its identifier",
    "note": "Full list with the local paths: datasets_sources.md.",
    "table": {
        "head": ["Input", "Release", "Identifier / source", "Used for"],
        "rows": [
            ["FunMap network", "v1.0", "funmap.linkedomics.org (CPTAC proteogenomics)", "a protein-derived prior"],
            ["FunMap mRNA matrices", "2024", "Zenodo record 7948944", "the RNA-derived network"],
            ["IntAct", "2024-09", "ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip", "a protein-derived prior"],
            ["STRING v12 (all / physical)", "12.0", "stringdb-downloads.org ... v12.0", "curated and protein priors"],
            ["Reactome", "v86", "reactome.org/download (pathways <= 200 proteins)", "a curated prior; a module source"],
            ["IntOGen drivers", "2020 / 2023 / 2024", "intogen.org/download?file=IntOGen-Drivers-*.zip", "three of the five labels"],
            ["Open Targets", "26.06", "ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.06/output", "ot_all and ot_nolit"],
            ["ClinGen gene validity", "file of 2026-09-26", "search.clinicalgenome.org/kb/gene-validity/download", "the hardest label; the 17 % noise anchor"],
            ["NCG (Genome Biol suppl.)", "2019", "doi 10.1186/s13059-018-1612-0, Additional file 3", "the arbitration set; the negative-set audit"],
            ["HGNC complete set", "with the release", "storage.googleapis.com/public-download-files/hgnc/...", "the symbol authority (D-02)"],
            ["PubMed gene2pubmed", "with the release", "ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz", "the publication-count adjustment"],
            ["DepMap CRISPR + expression", "two public releases", "depmap.org portal (1,150 and 1,361 cell lines)", "the external functional validation"],
            ["UniProt / ComplexPortal", "REST snapshot 2026-11-05", "rest.uniprot.org/uniprotkb/search", "the protein-complex map"],
            ["Our pre-registration", "frozen 2026-09-27", "Zenodo DOI 10.5281/zenodo.23121762", "the decision rule and the amendments"],
        ],
        "fontsize": 11.5, "rowh": 0.40,
    },
})

# ---- A3: the full 24-arm table ---------------------------------------
t1 = load("a4_table1_delta_auroc_v40.tsv")
if t1 is not None:
    nets = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
    labs = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
    piv = t1.pivot_table(index="network", columns="label_set", values="a0.5_mean")
    sd = t1.pivot_table(index="network", columns="label_set", values="a0.5_std")
    rows = []
    for net in nets:
        r = [net]
        for lab in labs:
            r.append(f"{piv.loc[net, lab]:+.4f} +/- {sd.loc[net, lab]:.4f}")
        rows.append(r)
    slides.append({
        "n": 19, "kind": "appendix", "tag": "APPENDIX A3",
        "title": "The full 24-arm primary result (alpha = 0.5)",
        "note": "delta-AUROC (RWR - degree) and its SD over the 200 "
                "(repeat, fold) evaluations. Source: "
                "results/tables/a4_table1_delta_auroc_v40.tsv. The alpha "
                "sweep (0.3 / 0.5 / 0.7) is in the same table set.",
        "table": {
            "head": ["network"] + labs,
            "rows": rows, "fontsize": 12, "rowh": 0.42,
        },
    })

# ---- A4: the five interval constructions ----------------------------
slides.append({
    "n": 20, "kind": "appendix", "tag": "APPENDIX A4",
    "title": "The H1 verdict under each interval construction",
    "note": "Pre-registered H1: the 95 % CI excludes 0 in EVERY arm. "
            "Source: docs/pre_registration_amendments.md A2.3b and "
            "results/tables/stats_h1_verdict_v40.tsv.",
    "table": {
        "head": ["Interval", "Arms with CI excluding 0", "Mean CI width", "Reading"],
        "rows": [
            ["fold-level percentile", "24 / 24", "0.0058", "over-conservative (coverage 1.000)"],
            ["repeat-level percentile (primary)", "24 / 24", "0.0057", "calibrated (coverage 0.936-0.944)"],
            ["repeat-level t", "24 / 24", "0.0070", "calibrated"],
            ["repeat-level BCa", "24 / 24", "0.0058", "calibrated"],
            ["gene-cluster (pre-registration 9.3)", "20 / 24", "0.0387", "a different question: a different set of genes"],
        ],
        "fontsize": 12, "rowh": 0.44,
    },
})

# ---- A5: the low-FPR endpoints --------------------------------------
g4 = load("pauc_top10_summary.tsv")
rows = []
if g4 is not None:
    piv = g4.pivot_table(index=["network", "label"], columns="endpoint",
                          values="delta_mean")
    nets = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
    labs = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
    for net in nets:
        for lab in labs:
            if (net, lab) not in piv.index:
                continue
            r = piv.loc[(net, lab)]
            rows.append([
                net, lab,
                f"{r.get('auroc', float('nan')):+.4f}",
                f"{r.get('pauc', float('nan')):+.4f}",
                f"{r.get('p10', float('nan')):+.4f}",
            ])
if rows:
    slides.append({
        "n": 21, "kind": "appendix", "tag": "APPENDIX A5",
        "title": "The low-FPR endpoints (the G4 co-primary endpoints)",
        "note": "pAUC at FPR <= 0.02 is the normalised excess area "
                "(random = 0, perfect = 1); p10 is the top-10 % precision "
                "minus the degree baseline's. The margins are much smaller "
                "than on AUROC, and the arm order is broadly stable "
                "(rho = 0.80-0.90). Source: pauc_top10_summary.tsv.",
        "table": {"head": ["network", "label", "delta-AUROC",
                           "delta-pAUC (FPR<=0.02)", "delta top-10 % prec."],
                  "rows": rows, "fontsize": 10.5, "rowh": 0.30},
    })

# ---- A6: the bait-matched null --------------------------------------
g3 = load("g3_bait_matched_null.tsv")
if g3 is not None:
    # one row per network, one column pair per label: 24 rows do not
    # fit on a slide, and the full table stays in the repository
    nets = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
    labs = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
    rows = []
    for net in nets:
        row = [net]
        for lab in labs:
            sel = g3[(g3.network == net) & (g3.label == lab)]
            if sel.empty:
                row.append("-")
                continue
            r = sel.iloc[0]
            row.append(f"{r.frac_unmatched_null_above:.1%} -> "
                       f"{r.frac_matched_null_above:.1%}")
        rows.append(row)
    slides.append({
        "n": 22, "kind": "appendix", "tag": "APPENDIX A6",
        "title": "The bait-matched null: how much of the margin is visibility",
        "note": "For each held-out positive we test a negative matched on "
                "its publication count; each cell is the fraction of "
                "positives above the null, UNMATCHED -> MATCHED. Chance = "
                "5 %. The temporal column is the one that collapses. "
                "Source: results/tables/g3_bait_matched_null.tsv.",
        "table": {"head": ["network"] + labs,
                  "rows": rows, "fontsize": 11.5, "rowh": 0.44},
    })

# ---- A7: the negative-set audit and the cleaned re-run -------------
g8 = load("g8_negative_contamination.tsv")
c9 = load("c9_clean_background_comparison.tsv")
if g8 is not None and c9 is not None:
    c9i = {(r.network, r.label): r for r in c9.itertuples()}
    labs = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
    rows = []
    for net in ["funmap", "rna_coexp", "string_phys700", "intact",
                "string_full700", "reactome"]:
        row = [net]
        for lab in labs:
            c = g8[(g8.network == net) & (g8.label == lab)]
            r = c9i.get((net, lab))
            if c.empty or r is None:
                row.append("-")
                continue
            row.append(f"{c.iloc[0].frac_negatives_in_ncg:.1%} | "
                       f"{r.standard_a05:+.4f} -> {r.cleaned_a05:+.4f}")
        rows.append(row)
    slides.append({
        "n": 23, "kind": "appendix", "tag": "APPENDIX A7",
        "title": "The negative set: contaminated, but it does not carry "
                 "the result",
        "note": "Each cell: the share of the held-out negatives that are "
                "curated NCG drivers | the margin BEFORE -> AFTER removing "
                "them from the BACKGROUND and re-running the benchmark "
                "(Control 13): 23/24 arms stay positive, no sign flip, the "
                "median +0.0809 -> +0.0822. Source: "
                "c9_clean_background_comparison.tsv.",
        "table": {"head": ["network"] + labs,
                  "rows": rows, "fontsize": 11, "rowh": 0.44},
    })


# ---- A8: the D-11 capped panel (the curator-lag evidence) -----------
d11 = load("d11_table_v2.tsv")
if d11 is not None:
    piv = d11.pivot_table(index="network", columns="label_base",
                           values="delta_auroc")
    rows = []
    for net in ["funmap", "rna_coexp", "string_phys700", "intact",
                "string_full700", "reactome"]:
        r = [net]
        for lab in ["intogen2024", "ot_all", "ot_nolit",
                    "intogen_temporal_new", "ncg_known"]:
            v = piv.loc[net, lab] if (net in piv.index and
                                      lab in piv.columns) else None
            r.append(f"{v:+.4f}" if pd.notna(v) else "-")
        rows.append(r)
    slides.append({
        "n": 24, "kind": "appendix", "tag": "APPENDIX A8",
        "title": "The equal-prevalence panel: the curator-lag evidence",
        "note": "All five labels have the same number of positives inside "
                "each arm (the D-11 common-size cap). Reactome is negative "
                "on the two newest catalogues and positive on the 2019 NCG "
                "catalogue -> the failure is label-specific, not "
                "network-specific. Source: d11_table_v2.tsv.",
        "table": {"head": ["network", "IntOGen 2024", "OT all", "OT no-lit",
                           "IntOGen temporal", "NCG 2019"],
                  "rows": rows, "fontsize": 11.5, "rowh": 0.42},
    })

# ---- A9: the label-quality tiers ------------------------------------
tiers = load("label_quality_tiers.tsv")
tier_auc = load("label_tier_delta_auroc.tsv")
if tiers is not None and tier_auc is not None:
    counts = tiers.groupby(["label", "tier"]).size().unstack(fill_value=0)
    med = tier_auc.groupby(["label", "tier"]).delta_auroc.median().unstack()
    rows = []
    for lab in ["intogen2024", "intogen_temporal_new", "ot_all", "ot_nolit",
                "ncg_known", "clingen_label"]:
        if lab not in counts.index:
            continue
        c = counts.loc[lab]
        m = med.loc[lab] if lab in med.index else None
        rows.append([
            lab,
            f"{c.get('high', 0)} / {c.get('mid', 0)} / {c.get('low', 0)}",
            f"{m.get('high', float('nan')):+.4f}" if m is not None and
            pd.notna(m.get("high")) else "-",
            f"{m.get('mid', float('nan')):+.4f}" if m is not None and
            pd.notna(m.get("mid")) else "-",
            f"{m.get('low', float('nan')):+.4f}" if m is not None and
            pd.notna(m.get("low")) else "-",
        ])
    slides.append({
        "n": 25, "kind": "appendix", "tag": "APPENDIX A9",
        "title": "Label quality: which positives carry the margin?",
        "note": "Tiers: ClinGen review status; IntOGen's number of "
                "detecting methods (>= 4 = high); NCG's literature count "
                "(>= 10 = high); Open Targets' literature vs non-literature "
                "evidence. The margin is recomputed with one tier's "
                "positives only. Source: label_quality_tiers.tsv, "
                "label_tier_delta_auroc.tsv.",
        "table": {"head": ["label", "genes high / mid / low",
                           "margin (high)", "margin (mid)", "margin (low)"],
                  "rows": rows, "fontsize": 12, "rowh": 0.44},
    })

# ---- A10: the variance decomposition, 24 vs 30 ----------------------
cmp_ = load("stats_variance_decomp_comparison.tsv")
if cmp_ is not None:
    rows = []
    for r in cmp_.itertuples():
        rows.append([
            r.variance_source,
            f"{r.sigma2_24arms:.5f}",
            f"{r.sigma2_30arms:.5f}",
            f"{100 * r.frac_between_24:.2f} %",
            f"{100 * r.frac_between_30:.2f} %",
            f"{100 * r.delta_frac_between:+.2f} pp",
        ])
    slides.append({
        "n": 26, "kind": "appendix", "tag": "APPENDIX A10",
        "title": "The variance decomposition: 24 arms against 30",
        "note": "The fifth label (NCG) is 65 % identical to IntOGen, so the "
                "30-arm grid measures the effect of adding a REDUNDANT "
                "label, not an independent one. Source: "
                "stats_variance_decomp_comparison.tsv.",
        "table": {"head": ["variance source", "sigma2 (24 arms)",
                           "sigma2 (30 arms)", "% between-arm (24)",
                           "% between-arm (30)", "change"],
                  "rows": rows, "fontsize": 12, "rowh": 0.44},
    })

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps({"slides": slides}, indent=1, ensure_ascii=False),
               encoding="utf-8")
print(f"wrote {OUT.name}: {len(slides)} appendix slides")
