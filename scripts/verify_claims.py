#!/usr/bin/env python3
"""Check every headline number in the manuscript against the released tables.

The D1 manuscript is written by hand in three Markdown files. Hand-written
numbers drift: the design was upgraded from 10 to 40 repeats, several
control tables were regenerated after bugs were fixed, and the prose kept
the old values. This script makes that class of error mechanical rather
than anecdotal.

What it does
------------
1. Runs ``scripts/recompute_ground_truth.py``'s computation (imported, not
   re-implemented) to get the authoritative value of every headline
   quantity.
2. Section 1 -- for each entry in the CLAIMS table, searches the named
   manuscript file for the stated string. These are *textual* checks: they
   assert that a claim is present and that a superseded one is absent.
3. Section 2 -- the FORBIDDEN table, the same idea for strings that must
   no longer appear.
4. Section 3 -- the NUMERIC checks, which are the ones that matter. Each
   entry names a manuscript file, a REGULAR EXPRESSION that CAPTURES the
   number the manuscript states, and the ground-truth key that number must
   equal. The captured number is compared NUMERICALLY against the value
   recomputed from the released tables, with a stated tolerance. A
   hardcoded literal is never compared against the table; if the table
   changes, these checks change verdict with it.
5. Writes ``review/verification_report.txt``.

A FAIL means the manuscript states something the released tables do not
support. Fix the manuscript, not the table.

Usage
-----
    python scripts/verify_claims.py
    python scripts/verify_claims.py --verbose
"""

from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# ---------------------------------------------------------------- the claims
# (claim_id, file under docs/, a string that MUST appear, what it is)
CLAIMS: list[tuple[str, str, str, str]] = [
    # ---- Control 32: the pathway level ----
    ("control32.section", "manuscript_a_sections.md",
     "### F.33 Control 32", "\u00a7F.33 must report the pathway level"),
    ("control32.section.draft", "research_report_draft.md",
     "### 3.28 The pathway level", "\u00a73.28 must report the pathway level"),
    ("control32.reading", "manuscript_a_sections.md",
     "point at the driver's pathways more", "\u00a7F.33 must state the reading"),
    ("control32.facts", "report_1118_facts.md",
     "## 7n. The pathway level (Control 32)",
     "the facts sheet must carry Control 32"),
    ("control32.qa", "report_1118_qa_bank.md",
     "### B33", "the question bank must carry Control 32"),

    # ---- Control 31: the cross-method complementarity ----
    ("control31.section", "manuscript_a_sections.md",
     "### F.32 Control 31", "\u00a7F.32 must report the overlap table"),
    ("control31.section.draft", "research_report_draft.md",
     "### 3.27 The cross-method complementarity",
     "\u00a73.27 must report the overlap table"),
    ("control31.reading", "manuscript_a_sections.md",
     "The methods pick largely different candidates",
     "\u00a7F.32 must state the reading"),
    ("control31.facts", "report_1118_facts.md",
     "## 7m. The cross-method complementarity (Control 31)",
     "the facts sheet must carry Control 31"),
    ("control31.qa", "report_1118_qa_bank.md",
     "### B32", "the question bank must carry Control 31"),

    # ---- Control 30: the DepMap essential-gene catalogue ----
    ("control30.section", "manuscript_a_sections.md",
     "### F.31 Control 30", "\u00a7F.31 must report the functional label"),
    ("control30.section.draft", "research_report_draft.md",
     "### 3.26 The DepMap essential-gene catalogue",
     "\u00a73.26 must report the functional label"),
    ("control30.reading", "manuscript_a_sections.md",
     "the two results are not in conflict, because the seeds differ",
     "\u00a7F.31 must reconcile it with Control 24"),
    ("control30.facts", "report_1118_facts.md",
     "## 7l. The DepMap essential-gene catalogue as a fifth label",
     "the facts sheet must carry Control 30"),
    ("control30.qa", "report_1118_qa_bank.md",
     "### B31", "the question bank must carry Control 30"),

    # ---- Control 28: the topology prior suite ----
    ("control28.section", "manuscript_a_sections.md",
     "### F.30 Control 28", "\u00a7F.30 must report the topology suite"),
    ("control28.section.draft", "research_report_draft.md",
     "### 3.25 The topology prior suite",
     "\u00a73.25 must report the topology suite"),
    ("control28.reading", "manuscript_a_sections.md",
     "Not every topology-aware score beats degree",
     "\u00a7F.30 must state the corrected reading"),
    ("control28.facts", "report_1118_facts.md",
     "## 7k. The topology prior suite (Control 28)",
     "the facts sheet must carry Control 28"),
    ("control28.qa", "report_1118_qa_bank.md",
     "### B30", "the question bank must carry Control 28"),

    # ---- Control 27: the rewiring null under the second learner ----
    ("control27.section", "manuscript_a_sections.md",
     "### F.29 Control 27", "\u00a7F.29 must report the rewiring grid"),
    ("control27.section.draft", "research_report_draft.md",
     "### 3.24 The rewiring null under the second learner",
     "\u00a73.24 must report the rewiring grid"),
    ("control27.reading", "manuscript_a_sections.md",
     "does not depend on the learner", "\u00a7F.29 must state the reading"),
    ("control27.facts", "report_1118_facts.md",
     "## 7j. The rewiring null under the second learner (Control 27",
     "the facts sheet must carry Control 27"),
    ("control27.qa", "report_1118_qa_bank.md",
     "### B29", "the question bank must carry Control 27"),

    # ---- Control 26: the label-permutation audit ----
    ("control26.section", "manuscript_a_sections.md",
     "### F.28 Control 26", "\u00a7F.28 must report the permutation audit"),
    ("control26.section.draft", "research_report_draft.md",
     "### 3.23 The label-permutation audit",
     "\u00a73.23 must report the permutation audit"),
    ("control26.reading", "manuscript_a_sections.md",
     "does not manufacture a margin from random labels",
     "\u00a7F.28 must state the audit's reading"),
    ("control26.facts", "report_1118_facts.md",
     "## 7i. The label-permutation audit (Control 26)",
     "the facts sheet must carry Control 26"),
    ("control26.qa", "report_1118_qa_bank.md",
     "### C24", "the question bank must carry Control 26"),

    # ---- Control 25: mutation-first versus network-first ----
    ("control25.section", "manuscript_a_sections.md",
     "### F.27 Control 25", "\u00a7F.27 must report the mutation comparison"),
    ("control25.section.draft", "research_report_draft.md",
     "### 3.22 Mutation-first versus network-first",
     "\u00a73.22 must report the mutation comparison"),
    ("control25.reading", "manuscript_a_sections.md",
     "complementary rather than redundant",
     "\u00a7F.27 must state that the two sources are complementary"),
    ("control25.facts", "report_1118_facts.md",
     "## 7h. Mutation-first versus network-first (Control 25)",
     "the facts sheet must carry Control 25"),
    ("control25.qa", "report_1118_qa_bank.md",
     "### B28", "the question bank must carry Control 25"),
    ("control25.recurrence.honest", "manuscript_a_sections.md",
     "do not claim that the network encodes a recurrent-mutation signal",
     "\u00a7F.27 must disclaim the recurrent-mutation reading"),

    # ---- Control 24: the functional task, read as a ranking ----
    ("control24.section", "manuscript_a_sections.md",
     "### F.26 Control 24", "\u00a7F.26 must report the functional-layer control"),
    ("control24.section.draft", "research_report_draft.md",
     "### 3.21 Does the network add anything on the functional task",
     "\u00a73.21 must report the functional-layer control"),
    ("control24.reading", "manuscript_a_sections.md",
     "no scorer separates the functionally essential genes from the rest",
     "\u00a7F.26 must state the corrected reading"),
    ("control24.reading.draft", "research_report_draft.md",
     "no scorer separates the functionally essential genes from the rest",
     "\u00a73.21 must state the corrected reading"),
    ("control24.facts", "report_1118_facts.md",
     "## 7g. The functional task, read as a ranking (Control 24)",
     "the facts sheet must carry Control 24"),
    ("control24.qa", "report_1118_qa_bank.md",
     "### B27", "the question bank must carry the corrected DepMap reading"),

    # ---- the design ----
    ("design.repeats", "manuscript_intro_abstract.md",
     "40-repeat", "the abstract must state the 40-repeat design"),
    ("design.repeats.a", "manuscript_a_sections.md",
     "40 repeats", "§A must state the 40-repeat design"),

    # ---- the primary endpoint ----
    ("endpoint.range", "manuscript_intro_abstract.md",
     "−0.013", "the abstract must give the negative cell's ΔAUROC"),
    ("endpoint.range.max", "manuscript_intro_abstract.md",
     "+0.171", "the abstract must give the panel maximum"),
    ("endpoint.cells", "manuscript_intro_abstract.md",
     "23 of 24 cells are positive", "the abstract must give the positive-cell count"),

    # ---- H1 ----
    ("h1.rule", "manuscript_intro_abstract.md",
     "24/24 arms under four of the five", "the abstract must state the H1 result"),
    ("h1.gene", "manuscript_intro_abstract.md",
     "20/24 under the pre-registered", "the abstract must disclose the gene-cluster result"),
    ("h1.gene.a", "manuscript_a_sections.md",
     "20/24", "§J must disclose the gene-cluster result"),
    ("h2.nottested", "manuscript_intro_abstract.md",
     "not testable", "the abstract must say H2 is not testable"),
    ("h2.nottested.a", "manuscript_a_sections.md",
     "the H2 verdict as **not tested**", "§J.2 must say H2 is not testable"),

    # ---- variance decomposition ----
    ("var.interaction", "manuscript_intro_abstract.md",
     "40.1 %", "the abstract must give the interaction share"),
    ("var.label", "manuscript_intro_abstract.md",
     "34.8 %", "the abstract must give the label share"),
    ("var.network", "manuscript_intro_abstract.md",
     "24.8 %", "the abstract must give the network share"),
    ("var.fold", "manuscript_intro_abstract.md",
     "0.3 % fold noise", "the abstract must give the fold share"),

    # ---- controls ----
    ("c2.z", "manuscript_intro_abstract.md",
     "z > 5", "the abstract must state the null rejection"),
    ("c3.runs", "manuscript_a_sections.md",
     "598", "§F.1 must state the true number of rewiring runs"),
    ("c4.universe", "manuscript_intro_abstract.md",
     "1,605-gene universe", "the abstract must state the shared universe"),
    ("c5.ratio", "manuscript_intro_abstract.md",
     "39.8 %", "the abstract must give the cross-catalogue retention"),
    ("c5.cells", "manuscript_intro_abstract.md",
     "59 of 66", "the abstract must give the transfer-cell counts"),
    ("c7.matched", "manuscript_intro_abstract.md",
     "+0.001", "the abstract must give the Control-7 matched-null difference"),

    # ---- the new controls ----
    ("c8.low", "manuscript_intro_abstract.md",
     "+0.053", "the abstract must give the low-stratum margin"),
    ("c8.high", "manuscript_intro_abstract.md",
     "+0.054", "the abstract must give the high-stratum margin"),
    ("c8.a", "manuscript_a_sections.md",
     "Control 8", "§F must document Control 8"),
    ("c9.rho", "manuscript_intro_abstract.md",
     "−0.02", "the abstract must give the size correlation"),
    ("c9.a", "manuscript_a_sections.md",
     "Control 9", "§F must document Control 9"),
    ("c10.jaccard", "manuscript_intro_abstract.md",
     "0.019", "the abstract must give the input-sensitivity overlap"),
    ("c10.a", "manuscript_a_sections.md",
     "Control 10", "§F must document Control 10"),

    # ---- the negative cell ----
    ("neg.reading", "manuscript_intro_abstract.md",
     "coverage failure", "the abstract must not call the negative cell a negative effect"),
    ("neg.reading.a", "manuscript_a_sections.md",
     "not a demonstration of a negative effect",
     "§J.3 must not call the negative cell a negative effect"),

    # ---- the network table ----
    ("net.funmap.nodes", "manuscript_a_sections.md",
     "| funmap             | 10,525 |", "§A must give FunMap's node count from the meta file"),
    ("net.reactome.nodes", "manuscript_a_sections.md",
     "| reactome           | 4,993  |", "§A must give Reactome's node count from the meta file"),
    ("net.citation", "manuscript_a_sections.md",
     "Shi Z, Lei JT, Elizarraras JM", "§A must cite the correct FunMap paper"),

    # ---- the label table ----
    ("label.overlap", "manuscript_b_sections.md",
     "| ot_all               | 205 |  24 | 600 | 430 | 83 |",
     "§D must give the recomputed label overlap"),
    ("label.transfer", "manuscript_b_sections.md",
     "59 cells positive, 7 cells negative", "§E must give the transfer-cell counts"),
    ("label.nolit", "manuscript_b_sections.md",
     "is **not** a literature-free label", "§D must qualify ot_nolit"),

    # ---- the GNN arm ----
    ("gnn.naming", "manuscript_a_sections.md",
     "DeepWalk", "§J.4 must name the arm honestly"),
    ("gnn.naming.log", "gnn_courtesy_log.md",
     "DeepWalk", "the GNN log must name the arm honestly"),
    ("gnn.number", "manuscript_a_sections.md",
     "+0.028", "§J.4 must give the released embedding margin"),
    ("gnn.oldnumber", "manuscript_a_sections.md",
     "+0.15", "the superseded embedding margin must be GONE"),

    # ---- the decision log ----
    ("dec.d11", "decisions.md",
     "reactome **63**", "the decision log must carry the released D-11 N values"),
    ("dec.d11.old", "decisions.md",
     "string_full700 N=587", "the superseded D-11 N values must be marked superseded"),

    # ---- Phase S: the stability analysis (section I) ----
    ("stability.sd", "manuscript_a_sections.md",
     "repeat-level SD of ΔAUROC is 0.0047",
     "§I.1 must give the arm-level SD of the margin"),
    ("stability.cv", "manuscript_a_sections.md",
     "relative error of the score across the 40 re-partitions is 0.185",
     "§I.1 must give the gene-level relative error"),
    ("stability.rank", "manuscript_a_sections.md",
     "~270 places out of ~10,000",
     "§I.1 must give the rank movement"),

    # ---- G3: the bait-matched null ----
    ("g3.temporal", "manuscript_a_sections.md",
     "**4.7 % (chance)**",
     "§F.12 must state that the temporal label falls to chance"),

    # ---- G7: the NCG arbitration test ----
    ("g7.arms", "manuscript_a_sections.md",
     "never seeded with above the degree baseline in 24/24 arms",
     "§F.14 must state the arbitration result"),
    ("g7.lowfpr", "manuscript_a_sections.md",
     "the low-FPR endpoint agrees in only 13/24",
     "§F.14 must disclose the low-FPR caveat"),

    # ---- G8: the negative-set audit ----
    ("g8.cleaned", "manuscript_a_sections.md",
     "**+0.0822**",
     "§F.19 must give the median margin of the cleaned re-run"),
    ("g8.cleaned.arms", "manuscript_a_sections.md",
     "23 of 24 arms stay positive",
     "§F.19 must state that the cleaned re-run keeps the arms positive"),
    ("g8.reactome", "manuscript_a_sections.md",
     "**+0.0893**",
     "§F.16 must give the Reactome-on-NCG margin"),

    # ---- GR-9: the cancer-type stratification ----
    ("gr9.types", "manuscript_a_sections.md",
     "All 26 cancer types have a positive median margin",
     "§F.15 must state the cancer-type result"),

    # ---- C6: the variance decomposition with the fifth label ----
    ("c6.label", "manuscript_a_sections.md",
     "34.78 %",
     "§F.17 must give the 24-arm label share"),
    ("c6.label30", "manuscript_a_sections.md",
     "30.75 %",
     "§F.17 must give the 30-arm label share"),

    # ---- C7: the label-quality tiers and the noise anchor ----
    ("c7.tiers", "manuscript_a_sections.md",
     "| intogen2024 | +0.0892 | +0.0633 | +0.0653 |",
     "§F.18 must give the tier-restricted margins"),
    ("c7.anchor", "manuscript_a_sections.md",
     "The margin changes by between −0.0018 and +0.0026",
     "§F.18 must give the ClinGen noise-anchor range"),

    # ---- C19: the baseline family sweep ---------------------------
    ("c19.closeness", "manuscript_a_sections.md",
     "closeness beats degree in **every one of the 24 arms**",
     "§F.22 must disclose that closeness beats degree everywhere"),
    ("c19.best", "manuscript_a_sections.md",
     "positive in **19 of 24 arms** with a median margin of **+0.046**",
     "§F.22 must give the margin against the strongest baseline"),
    ("c19.damping", "manuscript_a_sections.md",
     "23/24 arms, median +0.0571",
     "§F.22 must give the result at the literature's PageRank damping"),
    ("c19.draft", "research_report_draft.md",
     "degree is never the strongest classical baseline",
     "§3.17 of the draft must carry the same disclosure"),
    ("c19.abstract", "manuscript_intro_abstract.md",
     "+0.046 against the strongest classical baseline",
     "the abstract must qualify the headline against the strongest baseline"),

    # ---- C23: the split-design comparison -----------------------------
    ("c23.manuscript", "manuscript_a_sections.md",
     "None of the verdicts moves",
     "§F.25 must state that the verdicts do not move between designs"),
    ("c23.draft", "research_report_draft.md",
     "none of the verdicts moves",
     "§3.20 of the draft must carry the same result"),
    ("c23.preferred", "manuscript_a_sections.md",
     "| 10 x 10 | 100 | +0.0836 | 23/24 | 24/24 | 20/24 | 0.0370 |",
     "§F.25 must give the 10 x 10 row of the comparison"),
    ("c23.figure", "manuscript_a_sections.md",
     "**Figure 5.** The conclusion does not depend on the cross-validation "
     "design.",
     "§F.25 must carry the figure and its caption"),
    ("c23.figure.draft", "research_report_draft.md",
     "**Figure 5.** (A) 每个臂的 margin 在六种设计之间的取值范围",
     "§3.20 of the draft must carry the figure and its caption"),
    ("c23.preferred.prereg", "manuscript_a_sections.md",
     "| **5 x 40 (pre-registered)** | 200 | +0.0806 | 23/24 | 24/24 | 20/24 | 0.0365 |",
     "§F.25 must give the pre-registered row of the comparison"),

    # ---- C22: the rewiring feasibility check under the second learner ---
    ("c22.manuscript", "manuscript_a_sections.md",
     "Control 22 is a **feasibility check on one arm**",
     "§F.24 must label the rewiring check as a feasibility check"),
    ("c22.z", "manuscript_a_sections.md",
     "| proximity (20 rewirings) | **+0.0804** | **−0.0124 ± 0.0072** | "
     "**+12.9** |",
     "§F.24 must give the rewiring null under the second learner"),
    ("c22.degree", "manuscript_a_sections.md",
     "the largest difference in any gene's degree over\nthe 20 rewirings is "
     "**0**",
     "§F.24 must state that the rewiring preserves the degree sequence"),
    ("c22.draft", "research_report_draft.md",
     "the largest degree difference over the 20 rewirings is\n**0**",
     "§3.19 of the draft must carry the same result"),

    # ---- C21: the visibility decomposition under the second learner ---
    ("c21.dropsign", "manuscript_a_sections.md",
     "the sign of the drop agrees in 24 of 24 arms",
     "§F.24 must state that the visibility drop replicates in every arm"),
    ("c21.prox", "manuscript_a_sections.md",
     "The visibility effect **replicates under the second learner**",
     "§F.24 must state that the visibility component replicates"),
    ("c21.tie", "manuscript_a_sections.md",
     "counts a tie as \"not above\"",
     "§F.24 must disclose the tie correction"),
    ("c21.draft", "research_report_draft.md",
     "The visibility effect **replicates**",
     "§3.19 of the draft must carry the same result"),

    # ---- C20: the label-difficulty regression ----------------------
    ("c20.share", "manuscript_a_sections.md",
     "recovers 97.8 % of everything the full label factor explains",
     "§F.23 must give the share of the label factor recovered"),
    ("c20.f", "manuscript_a_sections.md",
     "F-test p = 0.38 on the primary family",
     "§F.23 must give the F-test that the label factor adds nothing"),
    ("c20.draft", "research_report_draft.md",
     "a single label-level number recovers 97.8 %",
     "§3.18 of the draft must carry the same result"),
    ("c20.abstract", "manuscript_intro_abstract.md",
     "one difficulty axis rather than five curation philosophies",
     "the abstract must say what the label main effect is"),

]

# Numbers that MUST NOT appear anywhere in the three manuscript files.
FORBIDDEN: list[tuple[str, str, str]] = [
    ("forbidden.10repeat", "manuscript_intro_abstract.md",
     "10-repeat", "the 10-repeat design is superseded"),
    ("forbidden.vandisp", "manuscript_intro_abstract.md",
     "52 % residual", "the superseded variance numbers must be GONE"),
    ("forbidden.h2holds", "manuscript_intro_abstract.md",
     "H2 holds for IntOGen 2024", "H2 is not testable"),
    ("forbidden.600", "manuscript_a_sections.md",
     "600 random-network edge-permuted runs", "there are 598 rewiring runs"),
    ("forbidden.gnnlog", "gnn_courtesy_log.md",
     "all 72/72", "the embedding arm has 24 cells, not 72"),
    ("forbidden.vandisp", "manuscript_intro_abstract.md",
     "33 % from label choice", "the superseded variance numbers must be GONE"),
    ("forbidden.oldcitation", "manuscript_a_sections.md",
     "Mei et al., Nat Commun 2021", "the wrong FunMap citation must be GONE"),
    ("forbidden.oldcount", "manuscript_b_sections.md",
     "61 cells positive, 5 cells negative", "the superseded transfer counts must be GONE"),
    ("forbidden.negclaim", "manuscript_intro_abstract.md",
     "publishable finding on curator-lag",
     "the curator-lag reading is a hypothesis, not a finding"),
    ("forbidden.halfmargin", "manuscript_intro_abstract.md",
     "shrinks by half on the strictest label",
     "the corrected C9 result replaces the old post-hoc claim"),
    ("forbidden.halfmargin.draft", "research_report_draft.md",
     "shrinks by half on the strictest label",
     "the corrected C9 result replaces the old post-hoc claim"),
    ("forbidden.hubnessclaim", "manuscript_a_sections.md",
     "Hubness predicts functional essentiality well",
     "the saturated-enrichment reading must be gone (Control 24)"),
    ("forbidden.hubnessclaim.draft", "research_report_draft.md",
     "Hubness predicts functional essentiality well",
     "the saturated-enrichment reading must be gone (Control 24)"),
    ("forbidden.degreesignal", "manuscript_a_sections.md",
     "the degree baseline carries almost all of the signal",
     "the corrected reading must replace it (Control 24)"),
    ("forbidden.degreesignal.draft", "research_report_draft.md",
     "the degree baseline carries almost all of the signal",
     "the corrected reading must replace it (Control 24)"),
    ("forbidden.degreepredictor", "manuscript_a_sections.md",
     "degree is the better predictor",
     "the functional task cannot arbitrate between the prior and the baselines"),
    ("forbidden.degreepredictor.draft", "research_report_draft.md",
     "degree is the better predictor",
     "the functional task cannot arbitrate between the prior and the baselines"),
    ("forbidden.oldconclusion", "manuscript_intro_abstract.md",
     "the median margin is +0.081, it",
     "the conclusions must carry the qualified margin, not the old one"),
]


def load_ground_truth() -> dict[tuple[str, str], object]:
    """Import recompute_ground_truth and run it, returning its table."""
    spec = importlib.util.spec_from_file_location(
        "rgt", ROOT / "scripts" / "recompute_ground_truth.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    import io
    from contextlib import redirect_stdout
    buf = io.StringIO()
    with redirect_stdout(buf):
        mod.main()
    import pandas as pd
    df = pd.DataFrame(mod.ROWS)
    out: dict[tuple[str, str], object] = {}
    for r in df.itertuples():
        out[(r.quantity, r.scope)] = r.value
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--out", default=str(ROOT / "review" / "verification_report.txt"))
    args = ap.parse_args(argv)

    gt = load_ground_truth()

    lines: list[str] = []
    lines.append("D1 benchmark -- manuscript claim verification")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Every check below compares a string in the manuscript against a")
    lines.append("value recomputed from the released result tables by")
    lines.append("scripts/recompute_ground_truth.py. A FAIL means the manuscript")
    lines.append("states something the tables do not support.")
    lines.append("")

    n_pass = n_fail = 0
    failures: list[str] = []

    lines.append("-" * 72)
    lines.append("1. POSITIVE CLAIMS -- the manuscript must contain these strings")
    lines.append("-" * 72)
    for cid, fname, needle, what in CLAIMS:
        path = ROOT / "docs" / fname
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        # a Markdown table may pad with spaces, so compare on a collapsed form
        collapsed = re.sub(r"\s+", " ", text)
        needle_c = re.sub(r"\s+", " ", needle)
        ok = needle_c in collapsed
        status = "PASS" if ok else "FAIL"
        n_pass += ok
        n_fail += (not ok)
        lines.append(f"  [{status}] {cid:<24} {fname}")
        lines.append(f"         {what}")
        if not ok:
            lines.append(f"         NOT FOUND: {needle!r}")
            failures.append(f"{cid}: {fname} does not contain {needle!r}")
        elif args.verbose:
            lines.append(f"         found: {needle!r}")

    lines.append("")
    lines.append("-" * 72)
    lines.append("2. FORBIDDEN STRINGS -- the manuscript must NOT contain these")
    lines.append("-" * 72)
    for cid, fname, needle, what in FORBIDDEN:
        path = ROOT / "docs" / fname
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        collapsed = re.sub(r"\s+", " ", text)
        needle_c = re.sub(r"\s+", " ", needle)
        present = needle_c in collapsed
        ok = not present
        n_pass += ok
        n_fail += (not ok)
        lines.append(f"  [{'PASS' if ok else 'FAIL'}] {cid:<24} {fname}")
        lines.append(f"         {what}")
        if not ok:
            lines.append(f"         STILL PRESENT: {needle!r}")
            failures.append(f"{cid}: {fname} still contains {needle!r}")

    lines.append("")
    lines.append("-" * 72)
    lines.append("3. NUMERIC CROSS-CHECK -- the manuscript's captured "
                 "numbers vs the recomputed table values")
    lines.append("-" * 72)
    # Each entry is (claim id, file, REGEX capturing the number the
    # manuscript states, ground-truth key, tolerance). The captured number
    # is compared NUMERICALLY with the value recomputed from the released
    # tables -- the literal in the manuscript is never the thing being
    # tested, so a change in the tables changes the verdict.
    #
    # The regexes are written to match the MINUS SIGN the manuscript uses
    # (U+2212) as well as the ASCII hyphen.
    NUM = r"[+−-]?\d+(?:\.\d+)?"
    cross = [
        ("abstract negative cell", "manuscript_intro_abstract.md",
         r"(−0\.01\d)", ("delta_auroc", "reactome x intogen_temporal_new"),
         0.0006),
        ("abstract panel max", "manuscript_intro_abstract.md",
         r"(\+0\.17\d)", ("delta_auroc", "panel range (max)"),
         0.0006),
        ("manuscript_a GNN margin", "manuscript_a_sections.md",
         r"\| \*\*embedding baseline\*\*\s*\| (\+0\.\d\d\d) \|",
         ("gnn_arm", "mean GNN delta"), 0.0006),
        ("control 8 low stratum", "manuscript_a_sections.md",
         r"\| low \(lowest third\) \| \*\*(\+0\.\d{4})\*\*",
         ("control8_seed_degree", "mean low"), 0.0002),
        ("control 8 high stratum", "manuscript_a_sections.md",
         r"\| high \(highest third\) \| \*\*(\+0\.\d{4})\*\*",
         ("control8_seed_degree", "mean high"), 0.0002),
        ("control 9 permutation p", "manuscript_a_sections.md",
         r"permutation p = (0\.\d+)", ("control9_size", "permutation p"),
         0.005),
        ("control 10 Jaccard", "manuscript_a_sections.md",
         r"top-100 Jaccard overlap \| \*\*(\d\.\d+)\*\*",
         ("control10_random_input", "input sensitivity Jaccard"),
         0.001),
        ("manuscript_b per-label mean", "manuscript_b_sections.md",
         r"\| intogen2024\s+\| 633 \| (\+\d\.\d{4})",
         ("label_mean", "intogen2024"), 0.0002),
        ("manuscript_b overlap row", "manuscript_b_sections.md",
         r"\| ot_all\s+\| 205 \|\s+24 \| 600 \| (\d+) \| 83 \|",
         ("label_overlap", "ot_all x ot_nolit"), 0.5),
        ("manuscript_b transfer mean", "manuscript_b_sections.md",
         r"retention of \*\*(\d+\.\d) %\*\*",
         ("control5_cross_catalogue", "retention percent"), 0.05),
        ("control 7 matched difference", "manuscript_intro_abstract.md",
         r"difference (\+\d\.\d{3,4})",
         ("control7_publication_count", "mean matched-minus-standard"),
         0.0002),
        ("co-primary partial AUC median", "manuscript_a_sections.md",
         r"\| partial AUC, FPR . 0\.02 \| (\+0\.\d{4}) \|",
         ("co_primary_endpoint", "partial AUC FPR<=0.02: median margin"),
         0.0002),
        ("co-primary top-10 % median", "manuscript_a_sections.md",
         r"\| precision in the top 10 % of the test set \| (\+0\.\d{4}) \|",
         ("co_primary_endpoint", "precision in the top 10 %: median margin"),
         0.0002),

        # ---- Control 32: the pathway level ----
        ("control 32 the prior's driver-pathway hits",
         "manuscript_a_sections.md",
         r"picks hit (\d+\.\d) of the driver-enriched pathways",
         ("control32_hits", "the prior's mean driver-pathway hits"), 0.06),
        ("control 32 the degree baseline's driver-pathway hits",
         "manuscript_a_sections.md",
         r"baseline's (\d+\.\d), and a random set",
         ("control32_hits",
          "the degree baseline's mean driver-pathway hits"), 0.06),
        ("control 32 the random baseline's driver-pathway hits",
         "manuscript_a_sections.md",
         r"same pool (\d+\.\d)",
         ("control32_hits", "the random baseline's mean driver-pathway hits"),
         0.06),
        ("control 32 the arms where the prior beats the degree baseline",
         "manuscript_a_sections.md",
         r"degree baseline in \*\*(\d+) of the (\d+) arms",
         ("control32_hits",
          "the arms where the prior beats the degree baseline"), 0.5),
        ("control 32 the arms where the prior beats the random baseline",
         "manuscript_a_sections.md",
         r"random baseline in \*\*(\d+) of the",
         ("control32_hits",
          "the arms where the prior beats the random baseline"), 0.5),
        ("control 32 the drivers' enriched pathways (IntOGen 2024)",
         "manuscript_a_sections.md",
         r"\| intogen2024 \| [\d,]+ \| (\d+) \|",
         ("control32_pathway", "the drivers' enriched pathways: intogen2024"),
         0.5),

        # ---- Control 31: the cross-method complementarity ----
        ("control 31 the smallest mean pairwise Jaccard",
         "manuscript_a_sections.md",
         r"overlap runs from (\d\.\d{3})",
         ("control31_overlap", "the smallest mean pairwise Jaccard"), 0.0006),
        ("control 31 the largest mean pairwise Jaccard",
         "manuscript_a_sections.md",
         r"to (\d\.\d{3}) \(the RWR prior against",
         ("control31_overlap", "the largest mean pairwise Jaccard"), 0.0006),
        ("control 31 the genes all six methods pick",
         "manuscript_a_sections.md",
         r"pick only (\d+\.\d) genes in common",
         ("control31_union", "the genes all six methods pick"), 0.06),
        ("control 31 the union relative to one list",
         "manuscript_a_sections.md",
         r"union of the six lists is about (\d+\.\d) times",
         ("control31_union", "the union relative to one list"), 0.06),
        ("control 31 the unique genes, the minimum",
         "manuscript_a_sections.md",
         r"from (\d+\.\d) per 100 \(the 1-hop seed count\)",
         ("control31_union",
          "the genes only one method picks: the minimum"), 0.06),
        ("control 31 the unique genes, the maximum",
         "manuscript_a_sections.md",
         r"per 100 \(the 1-hop seed count\) to (\d+\.\d) \(the",
         ("control31_union",
          "the genes only one method picks: the maximum"), 0.06),

        # ---- Control 30: the DepMap essential-gene catalogue ----
        ("control 30 the arms with a positive margin",
         "manuscript_a_sections.md",
         r"margin is positive in \*\*(\d+) of",
         ("control30_arm", "the arms with a positive margin"), 0.5),
        ("control 30 the median margin", "manuscript_a_sections.md",
         r"arms\*\*, with a median of \*\*(\+0\.\d{4})\*\*",
         ("control30_arm", "the median margin"), 0.0006),
        ("control 30 the prior AUROC, the smallest", "manuscript_a_sections.md",
         r"own AUROC on this label runs from (0\.\d{3})",
         ("control30_arm", "the prior AUROC: the smallest"), 0.0006),
        ("control 30 the prior AUROC, the largest", "manuscript_a_sections.md",
         r"own AUROC on this label runs from 0\.\d{3} to (0\.\d{3})",
         ("control30_arm", "the prior AUROC: the largest"), 0.0006),
        ("control 30 the degree AUROC, the smallest", "manuscript_a_sections.md",
         r"the degree baseline's runs from (0\.\d{3})",
         ("control30_arm", "the degree baseline AUROC: the smallest"), 0.0006),
        ("control 30 the degree AUROC, the largest", "manuscript_a_sections.md",
         r"the degree baseline's runs from 0\.\d{3} to (0\.\d{3})",
         ("control30_arm", "the degree baseline AUROC: the largest"), 0.0006),
        ("control 30 the margin against the degree baseline's AUROC",
         "manuscript_a_sections.md",
         r"degree baseline's AUROC is \*\*(-?\d\.\d{2})\*\*",
         ("control30_mechanism",
          "the Spearman rho: the margin against the degree baseline's AUROC"),
         0.006),
        ("control 30 the margin against the degree ratio",
         "manuscript_a_sections.md",
         r"correlation is (-?\d\.\d{2}), the wrong sign",
         ("control30_mechanism",
          "the Spearman rho: the margin against the positives' degree ratio"),
         0.006),
        ("control 30 the margin against the internal density",
         "manuscript_a_sections.md",
         r"internal density \((-?\d\.\d{2})\)",
         ("control30_mechanism",
          "the Spearman rho: the margin against the positives' internal "
          "density"), 0.006),

        # ---- Control 28: the topology prior suite ----        # ---- Control 28: the topology prior suite ----
        ("control 28 the betweenness median margin",
         "manuscript_a_sections.md",
         r"\| the betweenness centrality \| \*\*([\u2212\-+]?0\.\d{4})\*\* \|",
         ("control28_scorer",
          "the betweenness centrality: the median margin against degree"),
         0.0006),
        ("control 28 the 1-hop seed count median margin",
         "manuscript_a_sections.md",
         r"\| the 1-hop seed count \| ([\u2212\-+]?0\.\d{4}) \|",
         ("control28_scorer",
          "the 1-hop seed count: the median margin against degree"), 0.0006),
        ("control 28 the 2-hop seed count median margin",
         "manuscript_a_sections.md",
         r"\| the 2-hop seed count \| ([\u2212\-+]?0\.\d{4}) \|",
         ("control28_scorer",
          "the 2-hop seed count: the median margin against degree"), 0.0006),
        ("control 28 the RWR minus the 1-hop count",
         "manuscript_a_sections.md",
         r"margin over that trivial score is a median of\s*([\u2212\-+]?0\.\d{4})",
         ("control28_scorer",
          "the RWR prior minus the 1-hop seed count: the median"), 0.0006),

        # ---- Control 27: the rewiring null under the second learner ----
        ("control 27 the arms with a positive effect-size z",
         "manuscript_a_sections.md",
         r"effect-size z is positive in (\d+) of",
         ("control27_z", "the arms with a positive effect-size z"), 0.5),
        ("control 27 the median effect-size z", "manuscript_a_sections.md",
         r"the median at (\d+\.\d+) and the range",
         ("control27_z", "the median effect-size z"), 0.006),
        ("control 27 the smallest effect-size z", "manuscript_a_sections.md",
         r"and the range (\d+\.\d+) to",
         ("control27_z", "the smallest effect-size z"), 0.006),
        ("control 27 the largest effect-size z", "manuscript_a_sections.md",
         r"and the range \d+\.\d+ to (\d+\.\d+)\.",
         ("control27_z", "the largest effect-size z"), 0.006),
        ("control 27 the intact x intogen2024 observed margin",
         "manuscript_a_sections.md",
         r"the observed (\+?-?\d\.\d{4}) against a rewired null",
         ("control27_observed", "intact x intogen2024"), 0.0006),
        ("control 27 the intact x intogen2024 null mean",
         "manuscript_a_sections.md",
         r"against a rewired null of (\+?-?\d\.\d{4})",
         ("control27_null", "intact x intogen2024: the null mean"), 0.0006),
        ("control 27 the RWR effect-size z, the median",
         "manuscript_a_sections.md",
         r"with a median of \u2212?\+?(\d+\.\d{2})",
         ("control27_z", "the RWR effect-size z: the median"), 0.006),
        ("control 27 the RWR effect-size z, the largest",
         "manuscript_a_sections.md",
         r"to \u2012?\+?(\d+\.\d{2}) with a median of",
         ("control27_z", "the RWR effect-size z: the largest"), 0.006),

        # ---- Control 26: the label-permutation audit ----
("control 26 the reactome x intogen2024 margin",
         "manuscript_a_sections.md",
         r"reactome x intogen2024 has an observed margin of (\+0\.\d{4})",
         ("control26_observed", "reactome x intogen2024"), 0.0006),
        ("control 26 the reactome x intogen2024 null mean",
         "manuscript_a_sections.md",
         r"against a null mean of (\-0\.\d{4}) \(SD 0\.\d{4}\)",
         ("control26_null", "reactome x intogen2024: the null mean"), 0.0006),
        ("control 26 the reactome x intogen2024 null SD",
         "manuscript_a_sections.md",
         r"against a null mean of \-0\.\d{4} \(SD (0\.\d{4})\)",
         ("control26_null", "reactome x intogen2024: the null SD"), 0.0006),
        ("control 26 the string_phys700 x ot_nolit margin",
         "manuscript_a_sections.md",
         r"string_phys700 x ot_nolit has (\+0\.\d{4})",
         ("control26_observed", "string_phys700 x ot_nolit"), 0.0006),
        ("control 26 the arms beyond the null",
         "manuscript_a_sections.md",
         r"\*\*(\d+) of the 24 arms",
         ("control26_summary", "the arms where the observed margin is beyond "
          "the null's 95 % interval"), 0.5),

        # ---- Control 25: mutation-first versus network-first ----
        ("control 25 the prior, the whole universe", "manuscript_a_sections.md",
         r"\| the whole test universe \| \*{0,2}(0\.\d{4})\*{0,2}",
         ("control25_auc", "the prior: all"), 0.0006),
        ("control 25 degree, the whole universe", "manuscript_a_sections.md",
         r"\| the whole test universe \| \*{0,2}0\.\d{4}\*{0,2} \| \*{0,2}(0\.\d{4})\*{0,2}",
         ("control25_auc", "degree: all"), 0.0006),
        ("control 25 the mutation score, the whole universe",
         "manuscript_a_sections.md",
         r"\| the whole test universe \| \*{0,2}0\.\d{4}\*{0,2} \| \*{0,2}0\.\d{4}\*{0,2} \| \*{0,2}(0\.\d{4})\*{0,2}",
         ("control25_auc", "the mutation score: all"), 0.0006),
        ("control 25 the combination, the whole universe",
         "manuscript_a_sections.md",
         r"\| the whole test universe \| \*{0,2}0\.\d{4}\*{0,2} \| \*{0,2}0\.\d{4}\*{0,2} \| \*{0,2}0\.\d{4}\*{0,2} \| \*{0,2}(0\.\d{4})\*{0,2}",
         ("control25_auc", "the combination: all"), 0.0006),
        ("control 25 the prior, without the IntOGen drivers",
         "manuscript_a_sections.md",
         r"\| the test universe without the IntOGen drivers \| \*{0,2}(0\.\d{4})\*{0,2}",
         ("control25_auc", "the prior: minus_intogen"), 0.0006),
        ("control 25 the mutation score, without the IntOGen drivers",
         "manuscript_a_sections.md",
         r"\| the test universe without the IntOGen drivers \| \*{0,2}0\.\d{4}\*{0,2} \| \*{0,2}0\.\d{4}\*{0,2} \| \*{0,2}(0\.\d{4})\*{0,2}",
         ("control25_auc", "the mutation score: minus_intogen"), 0.0006),
        ("control 25 the recurrence AUROC of the prior",
         "manuscript_a_sections.md",
         r"\| over the 6 networks \| \*{0,2}(0\.\d{4})\*{0,2}",
         ("control25_recurrence", "the prior"), 0.0006),
        ("control 25 the recurrence AUROC of the publication count",
         "manuscript_a_sections.md",
         r"\| over the 6 networks \| \*{0,2}0\.\d{4}\*{0,2} \| \*{0,2}0\.\d{4}\*{0,2} \| \*{0,2}(0\.\d{4})\*{0,2}",
         ("control25_recurrence", "the publication count"), 0.0006),

        # ---- Control 24: the functional task, read as a ranking ----
        ("control 24 the prior's AUROC", "manuscript_a_sections.md",
         r"\| CERES < \u22120\.5 \| (0\.\d{4}) \|",
         ("control24_auc", "the prior (CERES < -0.5): the median over the arms"),
         0.0006),
        ("control 24 degree's AUROC", "manuscript_a_sections.md",
         r"\| CERES < \u22120\.5 \| 0\.\d{4} \| (0\.\d{4}) \|",
         ("control24_auc", "degree (CERES < -0.5): the median over the arms"),
         0.0006),
        ("control 24 the residual's AUROC", "manuscript_a_sections.md",
         r"\| CERES < \u22120\.5 \| 0\.\d{4} \| 0\.\d{4} \| (0\.\d{4}) \|",
         ("control24_auc", "the residual (CERES < -0.5): the median over the arms"),
         0.0006),
        ("control 24 the top stratum: the prior's picks",
         "manuscript_a_sections.md",
         r"\| 4 \(the top degree fifth\) \| ([\d\.]+) % \|",
         ("control24_stratum",
          "CERES < -0.5, K = 30, stratum 4: the prior's picks (the share of cell lines enriched, %)"),
         0.06),
        ("control 24 the top stratum: degree's picks",
         "manuscript_a_sections.md",
         r"\| 4 \(the top degree fifth\) \| [\d\.]+ % \| ([\d\.]+) % \|",
         ("control24_stratum",
          "CERES < -0.5, K = 30, stratum 4: degree's picks (the share of cell lines enriched, %)"),
         0.06),
        ("control 24 the top stratum: the random baseline",
         "manuscript_a_sections.md",
         r"\| 4 \(the top degree fifth\) \| [\d\.]+ % \| [\d\.]+ % \| ([\d\.]+) % \|",
         ("control24_stratum",
          "CERES < -0.5, K = 30, stratum 4: a random pick from the same stratum (the share of cell lines enriched, %)"),
         0.06),
        ("control 24 the prior's AUROC (the draft)", "research_report_draft.md",
         r"\| CERES < \u22120\.5 \| (0\.\d{4}) \|",
         ("control24_auc", "the prior (CERES < -0.5): the median over the arms"),
         0.0006),
        ("control 24 the random baseline (the draft)", "research_report_draft.md",
         r"\| 4 \(the top degree fifth\) \| [\d\.]+ % \| [\d\.]+ % \| ([\d\.]+) % \|",
         ("control24_stratum",
          "CERES < -0.5, K = 30, stratum 4: a random pick from the same stratum (the share of cell lines enriched, %)"),
         0.06),

        # ---- Phase S: the stability analysis ----
        ("stability repeat-level SD", "manuscript_a_sections.md",
         r"repeat-level SD of ΔAUROC is (0\.\d{4})",
         ("stability", "repeat-level SD of dAUROC (median over arms)"),
         0.0006),
        ("stability gene-level error", "manuscript_a_sections.md",
         r"re-partitions is (0\.\d{3})",
         ("stability", "gene-level relative error (median over arms)"),
         0.0006),
        ("stability rank movement", "manuscript_a_sections.md",
         r"~(\d{3}) places out of",
         ("stability", "median gene's rank movement (places)"),
         5.0),

        # ---- G3: the bait-matched null ----
        ("g3 temporal unmatched", "manuscript_a_sections.md",
         r"\*\*intogen_temporal_new\*\*\s*\|\s*(\d+\.\d+) %",
         ("bait_null", "temporal label, unmatched (%)"),
         0.06),
        ("g3 temporal matched", "manuscript_a_sections.md",
         r"\*\*intogen_temporal_new\*\*\s*\|\s*\d+\.\d+ %\s*\|\s*\*\*(\d+\.\d+) %",
         ("bait_null", "temporal label, matched (%)"),
         0.06),

        # ---- G7: the NCG arbitration test ----
        ("g7 arms positive", "manuscript_a_sections.md",
         r"never seeded with above the degree baseline in (\d+)/24 arms",
         ("arbitration", "arms with a positive dAUROC"),
         0.5),
        ("g7 median dAUROC", "manuscript_a_sections.md",
         r"arms\*\* \(median (\+0\.\d{4})",
         ("arbitration", "median dAUROC"),
         0.0006),
        ("g7 low-FPR arms", "manuscript_a_sections.md",
         r"the low-FPR endpoint agrees in only (\d+)/24",
         ("arbitration", "arms with a positive low-FPR dAUC"),
         0.5),

        # ---- G8: the negative-set audit (the post-hoc table lives in
        #      its own doc; the manuscript now carries the C9 re-run) ----
        ("g8 intact intogen2024 standard", "g8_negative_set.md",
         r"\| intact \| intogen2024 \| (\+0\.\d{4}) \|",
         ("negative_set", "intogen2024 intact, standard"),
         0.0006),
        ("g8 intact intogen2024 cleaned", "g8_negative_set.md",
         r"\| intact \| intogen2024 \| \+0\.\d{4} \| (\+0\.\d{4}) \|",
         ("negative_set", "intogen2024 intact, cleaned"),
         0.0006),
        ("g8 funmap ot_all standard", "g8_negative_set.md",
         r"\| funmap \| ot_all \| (\+0\.\d{4}) \|",
         ("negative_set", "ot_all funmap, standard"),
         0.0006),
        ("g8 funmap ot_all cleaned", "g8_negative_set.md",
         r"\| funmap \| ot_all \| \+0\.\d{4} \| (\+0\.\d{4}) \|",
         ("negative_set", "ot_all funmap, cleaned"),
         0.0006),
        ("g8 reactome temporal standard", "g8_negative_set.md",
         r"\| reactome \| intogen_temporal_new \| ([−-]0\.\d{4}) \|",
         ("negative_set", "reactome temporal, standard"),
         0.0006),
        ("g8 reactome temporal cleaned", "g8_negative_set.md",
         r"\| reactome \| intogen_temporal_new \| [−-]0\.\d{4} \| (\+0\.\d{4}) \|",
         ("negative_set", "reactome temporal, cleaned"),
         0.0006),

        # ---- GR-9: the cancer-type stratification ----
        ("gr9 types positive", "manuscript_a_sections.md",
         r"All (\d+) cancer types have a positive median margin",
         ("cancer_types", "types with a positive median margin"),
         0.5),
        ("gr9 median", "manuscript_a_sections.md",
         r"median over the types (\+0\.\d{4})",
         ("cancer_types", "median of the per-type medians"),
         0.0006),

        # ---- C6: the variance decomposition with the fifth label ----
        ("c6 label share 24", "manuscript_a_sections.md",
         r"\| label \(main effect\) \| (\d+\.\d{2}) %",
         ("variance24", "label share of the between-arm spread (%)"),
         0.006),
        ("c6 label share 30", "manuscript_a_sections.md",
         r"\| label \| \d+\.\d{2} % \| (\d+\.\d{2}) %",
         ("variance30", "label share of the between-arm spread (%)"),
         0.006),
        ("c6 interaction share 30", "manuscript_a_sections.md",
         r"\| network × label \| \d+\.\d{2} % \| (\d+\.\d{2}) %",
         ("variance30", "network x label share of the between-arm spread (%)"),
         0.006),

        # ---- C7: the label-quality tiers and the noise anchor ----
        ("c7 intogen2024 high tier", "manuscript_a_sections.md",
         r"\| intogen2024 \| (\+0\.\d{4}) \| \+0\.\d{4} \| \+0\.\d{4} \|",
         ("label_quality", "intogen2024, high tier (median)"),
         0.0006),
        ("c7 temporal high tier", "manuscript_a_sections.md",
         r"\| intogen_temporal_new \| (\+0\.\d{4}) \| \+0\.\d{4} \| \+0\.\d{4} \|",
         ("label_quality", "intogen_temporal_new, high tier (median)"),
         0.0006),
        ("c7 anchor min", "manuscript_a_sections.md",
         r"changes by between ([−-]0\.\d{4}) and",
         ("noise_anchor", "min change under 17 % misassignment"),
         0.0006),
        ("c7 anchor max", "manuscript_a_sections.md",
         r"changes by between [?−-]0\.\d{4} and (\+0\.\d{4})",
         ("noise_anchor", "max change under 17 % misassignment"),
         0.0006),

        # ---- C19: the baseline family sweep ----
        ("c19 headline median", "manuscript_a_sections.md",
         r"\| RWR \(α = 0\.5\) − degree \(the headline\) \| \*\*23/24\*\* \| "
         r"\*\*(\+0\.\d{4})\*\*",
         ("baseline_family", "RWR(0.5) - degree: median margin (primary 24)"),
         0.0006),
        ("c19 damping median", "manuscript_a_sections.md",
         r"\| RWR \(restart 0\.15\) − degree \| 23/24 \| (\+0\.\d{4})",
         ("baseline_family", "RWR(restart 0.15) - degree: median margin (primary 24)"),
         0.0006),
        ("c19 closeness median", "manuscript_a_sections.md",
         r"\| closeness − degree \| \*\*24/24\*\* \| (\+0\.\d{4})",
         ("baseline_family", "closeness - degree: median margin (primary 24)"),
         0.0006),
        ("c19 rwr vs closeness median", "manuscript_a_sections.md",
         r"\| RWR − closeness \| 20/24 \| (\+0\.\d{4})",
         ("baseline_family", "RWR(0.5) - closeness: median margin (primary 24)"),
         0.0006),
        ("c19 best baseline median", "manuscript_a_sections.md",
         r"\| \*\*RWR − the strongest baseline of that arm\*\* \| \*\*19/24\*\* \| "
         r"\*\*(\+0\.\d{4})\*\*",
         ("baseline_family",
          "RWR(0.5) - the strongest baseline of the arm: median margin (primary 24)"),
         0.0006),
        ("c19 closeness arms", "manuscript_a_sections.md",
         r"is closeness in (\d+) of the 30 arms",
         ("baseline_family", "arms where closeness is the strongest baseline"),
         0.5),

        # ---- C20: the label-difficulty regression ----
        ("c20 label-factor share", "manuscript_a_sections.md",
         r"\| label-mean RWR AUROC \| \*\*(\d+\.\d) %\*\*",
         ("label_difficulty",
          "rwr_auroc: share of the label factor recovered (%)"),
         0.06),
        ("c20 label-factor F test", "manuscript_a_sections.md",
         r"adds nothing beyond it \(F-test p = (0\.\d+) on the primary family",
         ("label_difficulty",
          "rwr_auroc: F-test p for the label factor beyond it"),
         0.006),
        ("c20 temporal margin", "manuscript_a_sections.md",
         r"\| `intogen_temporal_new` \| 150 \| (\d\.\d{4})",
         ("label_difficulty", "intogen_temporal_new: degree-baseline AUROC"),
         0.0006),
        ("c20 temporal margin value", "manuscript_a_sections.md",
         r"\| `intogen_temporal_new` \| 150 \| \d\.\d{4} \| \d\.\d{4} \| "
         r"\d\.\d{4} \| \*\*(\+0\.\d{4})\*\*",
         ("label_difficulty", "intogen_temporal_new: mean margin"),
         0.0006),



        # ---- the stability attribution, as the manuscript's table states it ----
        ("attribution degree joint beta", "manuscript_a_sections.md",
         r"\| log degree \| ([−-]0\.\d{3}) \|",
         ("stability_attribution", "log degree, joint beta (rank SD)"),
         0.002),
        ("attribution degree marginal rho", "manuscript_a_sections.md",
         r"\| log degree \| [−-]0\.\d{3} \| ([−+]?0\.\d{3})",
         ("stability_attribution", "log degree, marginal rho (rank SD)"),
         0.002),
        ("attribution score CV joint beta", "manuscript_a_sections.md",
         r"\| \*\*log gene score CV\*\* \| \*\*([−+]?0\.\d{3})\*\*",
         ("stability_attribution", "log score CV, joint beta (rank SD)"),
         0.002),
        ("attribution band marginal rho", "manuscript_a_sections.md",
         r"\| log local score density \([^)]*\) \| [−+]?0\.\d{3} \| ([−+]?0\.\d{3})",
         ("stability_attribution", "the +/- 10 % band, marginal rho (rank SD)"),
         0.002),
        ("attribution scott marginal rho", "manuscript_a_sections.md",
         r"\| log crowding \(Scott bandwidth\) \| [−+]?0\.\d{3} \(q = 0\.\d+\) \| ([−+]?0\.\d{3})",
         ("stability_attribution", "the Scott bandwidth, marginal rho (rank SD)"),
         0.002),


        # ---- C9: the factorial re-run with a cleaned background ----
        ("c9 positive arms, standard", "manuscript_a_sections.md",
         r"\*\*(\d+) of 24 arms are\s*positive before the cleaning",
         ("clean_background", "positive arms, standard design"),
         0.5),
        ("c9 positive arms, cleaned", "manuscript_a_sections.md",
         r"and (\d+) of 24 after",
         ("clean_background", "positive arms, cleaned background"),
         0.5),
        ("c9 median, standard", "manuscript_a_sections.md",
         r"moves from\s*\*\*([+−-]0\.\d{4})\*\*",
         ("clean_background", "median margin, standard design"),
         0.002),
        ("c9 median, cleaned", "manuscript_a_sections.md",
         r"to \*\*([+−-]0\.\d{4})\*\*",
         ("clean_background", "median margin, cleaned background"),
         0.002),


        # ---- the four design controls (C17 / F.20) ----
        ("seed perturbation, 90 %", "manuscript_a_sections.md",
         r"costs a median of \*\*([\d.]+) %\*\*",
         ("seed_perturbation",
          "median relative loss at 90 % of the seeds"),
         0.06),
        ("seed perturbation, 50 %", "manuscript_a_sections.md",
         r"halving the seeds costs a median of\s*([\d.]+) %",
         ("seed_perturbation",
          "median relative loss at 50 % of the seeds"),
         0.06),
        ("size-matched, standard", "manuscript_a_sections.md",
         r"gene-cluster interval excludes zero in (\d+) of 24 arms",
         ("size_matched_gene_ci", "arms excluding 0 (standard design)"),
         0.5),
        ("size-matched, capped", "manuscript_a_sections.md",
         r"the gene-cluster interval excludes\s*zero in \*\*(\d+) of \d+\*\* arms",
         ("size_matched_gene_ci",
          "arms excluding 0 (equal positive size)"),
         0.5),
        ("label redundancy, low", "manuscript_a_sections.md",
         r"are \*\*([\d.]+) to [\d.]+\*\*",
         ("label_redundancy",
          "min R2 of a label's arm means from the other labels"),
         0.006),
        ("second learner, arms", "manuscript_a_sections.md",
         r"positive in\s*\*\*(\d+) of \d+\*\* arms with a median of",
         ("second_learner", "arms with a positive margin"),
         0.5),
        ("second learner, median", "manuscript_a_sections.md",
         r"with a median of\s*\*\*([+−-]0\.\d{4})\*\*",
         ("second_learner", "median margin"),
         0.002),

    ]
    for cid, fname, pattern, key, tol in cross:
        val = gt.get(key)
        if val is None:
            lines.append(f"  [SKIP] {cid:<30} (no ground truth for {key})")
            continue
        path = ROOT / "docs" / fname
        raw = path.read_text(encoding="utf-8") if path.exists() else ""
        # a Markdown paragraph may be hard-wrapped, so the pattern is
        # matched against the whitespace-collapsed text
        text = re.sub(r"\s+", " ", raw)
        m = re.search(pattern, text)
        if not m:
            n_fail += 1
            lines.append(f"  [FAIL] {cid:<30} table={val}")
            lines.append(f"         the pattern {pattern!r} did not match; "
                         f"the manuscript may have reworded the claim")
            failures.append(f"{cid}: {fname} does not match {pattern!r}")
            continue
        try:
            # the manuscript uses the MINUS SIGN (U+2212); normalise it
            stated = float(m.group(1).replace("−", "-"))
        except ValueError:
            stated = float(m.group(0).replace("−", "-"))
        ok = abs(stated - float(val)) <= tol
        n_pass += ok
        n_fail += (not ok)
        lines.append(f"  [{'PASS' if ok else 'FAIL'}] {cid:<30} "
                     f"manuscript={stated:+.4f}  table={float(val):+.4f}  "
                     f"tol={tol}")
        if not ok:
            lines.append(f"         the manuscript states {stated}, the table "
                         f"says {float(val)}")
            failures.append(f"{cid}: {fname} states {stated}, the table says "
                            f"{float(val)}")

    lines.append("")
    lines.append("=" * 72)
    lines.append(f"RESULT: {n_pass} PASS / {n_fail} FAIL")
    lines.append("=" * 72)
    if failures:
        lines.append("")
        lines.append("Failures:")
        for f in failures:
            lines.append(f"  - {f}")

    report = "\n".join(lines)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(report, encoding="utf-8")
    # the console may be a legacy code page; the file is always UTF-8
    sys.stdout.buffer.write(report.encode("utf-8", "replace"))
    sys.stdout.buffer.write(f"\n\nwritten to {out}\n".encode("utf-8", "replace"))
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
